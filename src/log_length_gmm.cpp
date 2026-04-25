// log_length_gmm.cpp — 1D Gaussian mixture on a log-length histogram for
// auto-detecting natural bin edges. Pairs with taph::LengthBinStats.
//
// Pass 1 of a length-stratified workflow builds the histogram; this helper
// fits K∈[1, max] components via EM, selects K by BIC, and returns integer
// length edges at the density minima between consecutive component means.

#include "taph/log_length_gmm.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace taph {

namespace {

constexpr int    EM_MAX_ITER   = 80;
constexpr double EM_LL_EPS     = 1e-5;
constexpr double MIN_SIGMA     = 0.02;   // guard against collapse to a single bin
constexpr int    EDGE_GRID_PTS = 64;

struct GaussComp {
    double w    = 0.0;
    double mu   = 0.0;
    double sig2 = 0.0;
};

double gauss_pdf(double x, double mu, double sig2) {
    double d = x - mu;
    return std::exp(-0.5 * d * d / sig2) / std::sqrt(2.0 * M_PI * sig2);
}

// EM on binned data. `centers[i]` = log-length at bin i center; `counts[i]` = n.
// Returns (components, loglik). n_comp >= 1; centers.size() == counts.size().
std::pair<std::vector<GaussComp>, double> fit_gmm_em(
    const std::vector<double>& centers,
    const std::vector<double>& counts,
    int K,
    double x_lo,
    double x_hi)
{
    std::vector<GaussComp> C(K);
    double span = x_hi - x_lo;
    if (span < MIN_SIGMA) span = MIN_SIGMA;
    double sig0 = (span / std::max(K, 1)) * 0.5;
    double sig2_0 = std::max(sig0 * sig0, MIN_SIGMA * MIN_SIGMA);
    for (int k = 0; k < K; ++k) {
        C[k].w    = 1.0 / K;
        C[k].mu   = x_lo + (span * (k + 0.5)) / K;
        C[k].sig2 = sig2_0;
    }

    const std::size_t nb = centers.size();
    double total = 0.0;
    for (double c : counts) total += c;
    if (total <= 0.0) return {C, -std::numeric_limits<double>::infinity()};

    std::vector<double> resp(nb * K, 0.0);
    double prev_ll = -std::numeric_limits<double>::infinity();
    double ll = prev_ll;

    for (int it = 0; it < EM_MAX_ITER; ++it) {
        ll = 0.0;
        for (std::size_t i = 0; i < nb; ++i) {
            if (counts[i] <= 0.0) continue;
            double px = 0.0;
            for (int k = 0; k < K; ++k) {
                double p = C[k].w * gauss_pdf(centers[i], C[k].mu, C[k].sig2);
                resp[i * K + k] = p;
                px += p;
            }
            if (px <= 0.0) {
                for (int k = 0; k < K; ++k) resp[i * K + k] = 1.0 / K;
                continue;
            }
            for (int k = 0; k < K; ++k) resp[i * K + k] /= px;
            ll += counts[i] * std::log(px);
        }
        for (int k = 0; k < K; ++k) {
            double Nk = 0.0, muk = 0.0;
            for (std::size_t i = 0; i < nb; ++i) {
                double w = counts[i] * resp[i * K + k];
                Nk  += w;
                muk += w * centers[i];
            }
            if (Nk < 1e-9) {
                C[k].w = 1e-9;
                continue;
            }
            muk /= Nk;
            double s2 = 0.0;
            for (std::size_t i = 0; i < nb; ++i) {
                double w = counts[i] * resp[i * K + k];
                double d = centers[i] - muk;
                s2 += w * d * d;
            }
            s2 /= Nk;
            C[k].w    = Nk / total;
            C[k].mu   = muk;
            C[k].sig2 = std::max(s2, MIN_SIGMA * MIN_SIGMA);
        }
        if (std::abs(ll - prev_ll) < EM_LL_EPS * std::max(1.0, std::abs(ll))) break;
        prev_ll = ll;
    }

    std::sort(C.begin(), C.end(),
              [](const GaussComp& a, const GaussComp& b){ return a.mu < b.mu; });
    return {C, ll};
}

double mixture_density(double x, const std::vector<GaussComp>& C) {
    double s = 0.0;
    for (const auto& c : C) s += c.w * gauss_pdf(x, c.mu, c.sig2);
    return s;
}

double find_density_valley(double mu_a, double mu_b, const std::vector<GaussComp>& C) {
    double best_x = 0.5 * (mu_a + mu_b);
    double best_d = std::numeric_limits<double>::infinity();
    for (int i = 1; i < EDGE_GRID_PTS; ++i) {
        double t = static_cast<double>(i) / EDGE_GRID_PTS;
        double x = mu_a + t * (mu_b - mu_a);
        double d = mixture_density(x, C);
        if (d < best_d) { best_d = d; best_x = x; }
    }
    return best_x;
}

}  // namespace

LogLengthGmmResult detect_log_length_gmm_edges(
    const std::vector<uint64_t>& histogram,
    double log_min,
    double log_max,
    int min_length,
    int max_length,
    int max_components)
{
    LogLengthGmmResult out;
    out.n_components = 1;
    out.converged = false;
    out.bic = 0.0;

    const std::size_t nb = histogram.size();
    if (nb < 4 || max_length <= min_length || log_max <= log_min) {
        return out;
    }

    std::vector<double> centers(nb), counts(nb);
    double step = (log_max - log_min) / static_cast<double>(nb);
    double total = 0.0;
    for (std::size_t i = 0; i < nb; ++i) {
        centers[i] = log_min + (i + 0.5) * step;
        counts[i]  = static_cast<double>(histogram[i]);
        total += counts[i];
    }
    if (total < 256.0) return out;

    int Kmax = std::clamp(max_components, 1, 4);
    double best_bic  = std::numeric_limits<double>::infinity();
    int    best_K    = 1;
    std::vector<GaussComp> best_C;
    double best_ll   = -std::numeric_limits<double>::infinity();

    for (int K = 1; K <= Kmax; ++K) {
        auto [C, ll] = fit_gmm_em(centers, counts, K, log_min, log_max);
        double p  = 3.0 * K - 1.0;
        double bic = -2.0 * ll + p * std::log(total);
        if (bic < best_bic) {
            best_bic = bic;
            best_K   = K;
            best_C   = std::move(C);
            best_ll  = ll;
        }
    }

    out.bic = best_bic;
    out.converged = std::isfinite(best_ll);
    out.n_components = best_K;

    if (best_K <= 1 || best_C.empty()) {
        return out;
    }

    // Require sufficient separation: means differ by >= 0.25 in log space
    // (≈ 30% length change) AND each component weight >= 3%.
    std::vector<GaussComp> kept;
    for (const auto& c : best_C) {
        if (c.w >= 0.03) kept.push_back(c);
    }
    if (kept.size() < 2) {
        out.n_components = 1;
        return out;
    }
    std::vector<GaussComp> merged;
    merged.push_back(kept.front());
    for (std::size_t i = 1; i < kept.size(); ++i) {
        if (kept[i].mu - merged.back().mu < 0.25) {
            double w = merged.back().w + kept[i].w;
            merged.back().mu   = (merged.back().w * merged.back().mu
                                  + kept[i].w * kept[i].mu) / w;
            merged.back().sig2 = (merged.back().w * merged.back().sig2
                                  + kept[i].w * kept[i].sig2) / w;
            merged.back().w    = w;
        } else {
            merged.push_back(kept[i]);
        }
    }
    if (merged.size() < 2) {
        out.n_components = 1;
        return out;
    }
    out.n_components = static_cast<int>(merged.size());

    std::vector<int> edges;
    for (std::size_t i = 0; i + 1 < merged.size(); ++i) {
        double x_star = find_density_valley(merged[i].mu, merged[i + 1].mu, merged);
        int L = static_cast<int>(std::round(std::exp(x_star)));
        L = std::clamp(L, min_length + 1, max_length - 1);
        if (edges.empty() || L > edges.back()) edges.push_back(L);
    }
    // After clamp/dedup, the realized component count is edges+1. Without this
    // sync, n_components could disagree with the number of usable bins
    // downstream (silent metadata drift).
    out.edges = std::move(edges);
    out.n_components = out.edges.empty()
                       ? 1
                       : static_cast<int>(out.edges.size()) + 1;
    return out;
}

std::vector<int> detect_quantile_length_edges(
    const std::vector<uint64_t>& histogram,
    double log_min,
    double log_max,
    int min_length,
    int max_length,
    int n_bins)
{
    std::vector<int> edges;
    if (n_bins <= 1) return edges;
    const std::size_t nb = histogram.size();
    if (nb == 0 || max_length <= min_length || log_max <= log_min) return edges;

    double total = 0.0;
    for (uint64_t c : histogram) total += static_cast<double>(c);
    if (total <= 0.0) return edges;

    double step = (log_max - log_min) / static_cast<double>(nb);
    double cum = 0.0;
    int next_edge_k = 1;
    for (std::size_t i = 0; i < nb; ++i) {
        cum += static_cast<double>(histogram[i]);
        while (next_edge_k < n_bins && cum / total >= static_cast<double>(next_edge_k) / n_bins) {
            double x = log_min + (i + 0.5) * step;
            int L = static_cast<int>(std::round(std::exp(x)));
            L = std::clamp(L, min_length + 1, max_length - 1);
            if (edges.empty() || L > edges.back()) edges.push_back(L);
            ++next_edge_k;
        }
    }
    return edges;
}

}  // namespace taph
