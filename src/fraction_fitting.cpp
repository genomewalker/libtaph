#include "taph/fraction_fitting.hpp"
#include <cmath>
#include <limits>

namespace taph {

double pool_interior_bg(const int64_t* t, const int64_t* tc,
                        int n_pos, double fallback_bg, int start_pos)
{
    constexpr int64_t MIN_COV = 50;
    double sum = 0.0; int cnt = 0;
    for (int p = start_pos; p < n_pos; ++p) {
        if (tc[p] < MIN_COV) continue;
        sum += static_cast<double>(t[p]) / tc[p];
        ++cnt;
    }
    return cnt >= 2 ? sum / cnt : fallback_bg;
}

std::pair<double,double> fit_exp_decay_irls(
    const int64_t* t, const int64_t* tc, int n_pos,
    double bg, bool skip_pos0)
{
    constexpr int64_t MIN_FIT_COV = 50;
    constexpr int FIT_START = 1;
    constexpr int FIT_END   = 7;   // fit positions 1..6 (exclusive end)

    auto ols = [&](const double* w) -> std::pair<double,double> {
        double sx=0,sy=0,sxx=0,sxy=0,sw=0;
        for (int p = FIT_START; p < FIT_END && p < n_pos; ++p) {
            if (tc[p] < MIN_FIT_COV) continue;
            double rate = static_cast<double>(t[p]) / tc[p];
            double excess = rate - bg;
            double se2 = rate * (1.0 - rate) / tc[p];
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

    double w2[FIT_END - FIT_START] = {};
    for (int p = FIT_START; p < FIT_END && p < n_pos; ++p) {
        if (tc[p] < MIN_FIT_COV) { w2[p - FIT_START] = 0.0; continue; }
        double rate = static_cast<double>(t[p]) / tc[p];
        double exc_fit = dmax0 * std::exp(-lam0 * p);
        double var_rate = rate * (1.0 - rate) / tc[p];
        if (var_rate <= 0.0 || exc_fit <= 0.0) { w2[p - FIT_START] = 0.0; continue; }
        w2[p - FIT_START] = (exc_fit * exc_fit) * tc[p] / var_rate;
    }
    auto [dmax1, lam1] = ols(w2);
    if (dmax1 == 0.0)
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};

    if (skip_pos0) {
        double d1 = dmax1 * std::exp(-lam1);
        return {d1, lam1};
    }
    return {dmax1, lam1};
}

} // namespace taph
