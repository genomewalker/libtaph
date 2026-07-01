#include "taph/fraction_fitting.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace taph {

BgEstimate pool_interior_bg(const int64_t* t, const int64_t* tc,
                             int n_pos, double fallback_bg, int start_pos)
{
    constexpr int64_t MIN_COV = 50;
    double sum = 0.0, var_sum = 0.0; int cnt = 0;
    for (int p = start_pos; p < n_pos; ++p) {
        if (tc[p] < MIN_COV) continue;
        double r = static_cast<double>(t[p]) / tc[p];
        sum += r;
        var_sum += r * (1.0 - r) / tc[p];
        ++cnt;
    }
    if (cnt < 2) return {fallback_bg, 0.0};
    double m = static_cast<double>(cnt);
    return {sum / m, var_sum / (m * m)};  // var = (1/m^2) * sum Var(rate_i)
}

// Design notes (multi-model review):
// 1. Arithmetic mixing IS the correct generative model for a two-component binomial
//    mixture: E[T/(C+T)] = pi*d_anc*exp(-lam*p) + bg  (linear, not logit).
//    The bulk GLM's logit link is a single-component convenience approximation; near
//    bg=0.5 (Chargaff) its curvature is negligible (<0.3% bias on d_max) because
//    sigmoid is approximately linear at its inflection point. pi = d_bulk/d_anc is
//    correct as long as both use the same arithmetic normalization.
// 2. bg is shared across all fitting positions, introducing a rank-1 covariance
//    diag(Var(rate)) + Var(bg)*J. The diagonal fix below (adding bg.var to each
//    per-position variance) corrects lambda's SE (~7% wider) but underestimates
//    Var(d_max) because d_max's OLS coefficient vector projects heavily onto the
//    constant direction: true SE gain for d_max is ~25-35%, needing GLS with
//    Sherman-Morrison: Sigma^{-1} = diag(1/v) - Var(bg)*diag(1/v)*J*diag(1/v) /
//    (1 + Var(bg)*sum_p(1/v_p)). Tracked as a future improvement.
std::pair<double,double> fit_exp_decay_irls(
    const int64_t* t, const int64_t* tc, int n_pos,
    BgEstimate bg, bool skip_pos0)
{
    constexpr int64_t MIN_FIT_COV = 50;
    constexpr int FIT_START = 1;
    constexpr int FIT_END   = 7;   // fit positions 1..6 (exclusive end)

    auto ols = [&](const double* w) -> std::pair<double,double> {
        double sx=0,sy=0,sxx=0,sxy=0,sw=0;
        for (int p = FIT_START; p < FIT_END && p < n_pos; ++p) {
            if (tc[p] < MIN_FIT_COV) continue;
            double rate = static_cast<double>(t[p]) / tc[p];
            double excess = rate - bg.mean;
            // Var(excess) = Var(rate[p]) + Var(bg); bg.var is the variance of the
            // background mean estimate (propagated from interior positions).
            double se2 = rate * (1.0 - rate) / tc[p] + bg.var;
            if (se2 <= 0.0) continue;
            if (excess <= 2.0 * std::sqrt(se2)) continue;
            double logy = std::log(excess);
            double wi = (w ? w[p - FIT_START] : 1.0);
            sx  += wi * p;
            sy  += wi * logy;
            sxx += wi * p * p;
            sxy += wi * p * logy;
            sw  += wi;
        }
        if (sw <= 0.0) return {0.0, 0.0};
        double denom = sw * sxx - sx * sx;
        if (denom <= 0.0) return {0.0, 0.0};
        double slope     = (sw * sxy - sx * sy) / denom;
        double intercept = (sy - slope * sx) / sw;
        double lam  = -slope;
        double dmax = std::exp(intercept);
        if (lam < 0.0 || dmax <= 0.0 || dmax > 1.0) return {0.0, 0.0};
        return {dmax, lam};
    };

    auto [dmax0, lam0] = ols(nullptr);
    if (dmax0 == 0.0)
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};

    // Second pass: rank-1 GLS for Sigma_y = D + bg.var * u*u^T
    // D_p = var_rate_p / exc_p^2,  u_p = 1/exc_p  =>  D^{-1}u = g_p = exc_p/var_rate_p
    // Sherman-Morrison: Sigma_y^{-1} = D^{-1} - f*(D^{-1}u)(D^{-1}u)^T, f=bg.var/(1+bg.var*k)
    // IMPORTANT: bg.var enters ONLY via the rank-1 term; adding it to var_rate would double-count.
    {
        double A00=0,A01=0,A11=0, B0=0,B1=0;   // X^T D^{-1} X,  X^T D^{-1} y
        double G0=0,G1=0, Gy=0, kappa=0;       // X^T g,  g^T y,  sum(1/var_rate)
        int n_pts = 0;
        for (int p = FIT_START; p < FIT_END && p < n_pos; ++p) {
            if (tc[p] < MIN_FIT_COV) continue;
            double rate    = static_cast<double>(t[p]) / tc[p];
            double excess  = rate - bg.mean;
            double var_rate = rate * (1.0 - rate) / tc[p];  // diagonal only — NOT +bg.var
            double se2_tot  = var_rate + bg.var;
            if (se2_tot <= 0.0 || excess <= 2.0 * std::sqrt(se2_tot)) continue;
            double y       = std::log(excess);
            double exc_fit = dmax0 * std::exp(-lam0 * p);
            if (exc_fit <= 0.0 || var_rate <= 0.0) continue;
            double a  = exc_fit * exc_fit / var_rate;   // D^{-1} weight
            double g  = exc_fit / var_rate;             // D^{-1} u_p
            double pp = static_cast<double>(p);
            A00 += a;    A01 -= a*pp;   A11 += a*pp*pp;
            B0  += a*y;  B1  -= a*pp*y;
            G0  += g;    G1  -= g*pp;   Gy  += g*y;
            kappa += 1.0 / var_rate;
            ++n_pts;
        }
        if (n_pts < 2 || A00 <= 0.0)
            return {std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN()};
        double f   = bg.var / (1.0 + bg.var * kappa);
        double M00 = A00 - f*G0*G0,  M01 = A01 - f*G0*G1,  M11 = A11 - f*G1*G1;
        double C0  = B0  - f*G0*Gy,  C1  = B1  - f*G1*Gy;
        double det = M00*M11 - M01*M01;
        if (det <= 0.0)
            return {std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN()};
        double alpha = (C0*M11 - C1*M01) / det;   // log(d_max)
        double lam1  = (C1*M00 - C0*M01) / det;   // decay rate
        double dmax1 = std::exp(alpha);
        if (lam1 <= 0.0 || dmax1 <= 0.0 || dmax1 > 1.0)
            return {std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN()};
        if (skip_pos0) return {dmax1 * std::exp(-lam1), lam1};
        return {dmax1, lam1};
    }
}

CtSpecificDecayFit fit_ct_specific_decay(
    const std::array<uint64_t, 64>* tri_pos, int n_pos)
{
    const CtSpecificDecayFit FAIL{std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::quiet_NaN(), 0};
    constexpr int MIN_N = 10;  // minimum FROM+TO count to trust a position's rate

    auto rate_at = [&](int p, int from, int to) -> std::pair<double,double> {
        // returns {rate, weight=nf+nt}, weight=-1 if under-covered
        uint64_t nf = 0, nt = 0;
        for (int prev = 0; prev < 4; ++prev)
            for (int next = 0; next < 4; ++next) {
                nf += tri_pos[p][prev * 16 + from * 4 + next];
                nt += tri_pos[p][prev * 16 + to   * 4 + next];
            }
        if (nf + nt < MIN_N) return {0.0, -1.0};
        return {static_cast<double>(nt) / (nf + nt), static_cast<double>(nf + nt)};
    };

    // Weighted log-linear regression: log(g(p)) ~ intercept - lambda*p,
    // over positions where g(p) = CT(p) - 0.5*(CA(p)+CG(p)) is positive.
    // Weight by the harmonic-ish combined coverage of the three channels at p.
    double sw=0, sx=0, sy=0, sxx=0, sxy=0;
    int n_pts = 0;
    for (int p = 1; p < n_pos; ++p) {
        auto [ct, w_ct] = rate_at(p, 1, 3);
        auto [ca, w_ca] = rate_at(p, 1, 0);
        auto [cg, w_cg] = rate_at(p, 1, 2);
        if (w_ct < 0 || w_ca < 0 || w_cg < 0) continue;
        double g = ct - 0.5 * (ca + cg);
        if (g <= 0.005) continue;
        double w = std::min({w_ct, w_ca, w_cg});
        double logy = std::log(g);
        sw += w; sx += w*p; sy += w*logy; sxx += w*p*p; sxy += w*p*logy;
        ++n_pts;
    }
    if (n_pts < 3) return FAIL;
    double denom = sw * sxx - sx * sx;
    if (std::abs(denom) < 1e-12) return FAIL;
    double slope     = (sw * sxy - sx * sy) / denom;
    double intercept = (sy - slope * sx) / sw;
    double lambda = -slope;
    double A = std::exp(intercept);
    if (lambda <= 0.0 || A <= 0.0 || A > 1.0) return FAIL;
    return {A, lambda, n_pts};
}

} // namespace taph
