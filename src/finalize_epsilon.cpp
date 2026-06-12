// Reference-free 8-oxoG floor ε₀ from interior Chargaff contrast T/(T+G) − A/(A+C).
//
// Estimator: per LenBinStats bin b, using interior counts only:
//   ε_l(b) = T/(T+G) − A/(A+C)
// where T = t_interior[b], G = g_interior[b], A = a_interior[b], C = c_interior[b].
//
// The Chargaff control A/(A+C) cancels genomic base composition to first order (PR2: within-strand
// A≈T, C≈G in undamaged DNA). Undamaged → ε_l ≈ 0. 8-oxoG raises T/(T+G) while leaving A/(A+C)
// flat → ε_l > 0. C→T deamination is excluded by geometry: interior positions (≥ L/3) have
// negligible overhang deamination; f0 partially cancels in the difference (T up, C down moves both
// terms the same direction).
//
// Library-level ε₀ = weighted mean over bins with min_counts ≥ 200, weight = min(T+G, A+C).
// Flatness diagnostic: ε_term = mean of T/(T+G)−A/(A+C) at positions {0,1,2} pooled across bins.
// Authentic oxidation: ε_term ≈ ε₀ (position-flat). C→T leakage: ε_term > ε₀ (end-elevated).
//
// phi_share = ε₀/(ε₀+f0) is set by damage_estimation_finalize after finalize_tau completes.

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "taph/sample_damage_profile.hpp"
#include "taph/read_ancient_llr.hpp"

namespace taph {

namespace {
constexpr uint64_t MIN_INTERIOR_COUNTS = 200;
constexpr uint64_t MIN_TERM_COUNTS     = 50;
}

EpsilonEstimate finalize_epsilon(const SampleDamageProfile& profile) {
    EpsilonEstimate out;

    // Per-bin ε_l from interior counts
    double sum_w = 0.0, sum_we = 0.0, sum_we2 = 0.0;
    int n_bins = 0;

    for (const auto& b : profile.len_bins) {
        const uint64_t denom_tg = b.t_interior + b.g_interior;
        const uint64_t denom_ac = b.a_interior + b.c_interior;
        if (denom_tg < MIN_INTERIOR_COUNTS || denom_ac < MIN_INTERIOR_COUNTS) continue;

        const double tg_rate = static_cast<double>(b.t_interior) / static_cast<double>(denom_tg);
        const double ac_rate = static_cast<double>(b.a_interior) / static_cast<double>(denom_ac);
        const double eps_l   = tg_rate - ac_rate;

        // Weight by effective sample size: min of the two denominators
        const double w = static_cast<double>(std::min(denom_tg, denom_ac));

        // Running sufficient stats for weighted mean + variance
        sum_w   += w;
        sum_we  += w * eps_l;
        sum_we2 += w * eps_l * eps_l;
        ++n_bins;
    }

    if (n_bins == 0 || sum_w <= 0.0) return out;

    const double eps0    = sum_we / sum_w;
    // Weighted sample variance (population-weighted, not Bessel-corrected — large n)
    const double var_eps = std::max(0.0, sum_we2 / sum_w - eps0 * eps0);
    const double se_eps  = std::sqrt(var_eps / static_cast<double>(n_bins));

    out.epsilon_floor = eps0;
    out.epsilon_lo    = eps0 - 1.96 * se_eps;
    out.epsilon_hi    = eps0 + 1.96 * se_eps;
    out.n_bins        = n_bins;

    // Flatness diagnostic: ε at positions 0–2 pooled across all bins
    // Use 5' t_counts[p], g_counts[p], a_counts[p], c_counts[p], p in {0,1,2}
    uint64_t term_T = 0, term_G = 0, term_A = 0, term_C = 0;
    for (const auto& b : profile.len_bins) {
        for (int p = 0; p < 3; ++p) {
            if (b.t_counts[p] + b.g_counts[p] < MIN_TERM_COUNTS) continue;
            term_T += b.t_counts[p];
            term_G += b.g_counts[p];
            term_A += b.a_counts[p];
            term_C += b.c_counts[p];
        }
    }
    if (term_T + term_G > 0 && term_A + term_C > 0) {
        const double tg = static_cast<double>(term_T) / static_cast<double>(term_T + term_G);
        const double ac = static_cast<double>(term_A) / static_cast<double>(term_A + term_C);
        out.epsilon_term = tg - ac;
    }

    out.fitted = true;
    return out;
}

} // namespace taph
