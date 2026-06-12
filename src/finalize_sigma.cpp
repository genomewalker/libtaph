// Reference-free GC→AT depletion pressure σ₀ = T/(T+G) + A/(A+C) − 1.
//
// σ is the SYMMETRIC channel of the base-count quadruplet: it captures both 8-oxoG (G→T, C→A
// equally) and interior C→T deamination (f0 floor), along with intrinsic GC composition. The
// ASYMMETRIC channel ε = T/(T+G) − A/(A+C) cancels symmetric oxidation by construction and is
// computed separately in finalize_epsilon as a leakage control.
//
// σ₀ is pooled from interior counts across all length bins (MLE: sum numerators/denominators,
// not mean-of-ratios). It is composition-confounded: damage-free GC-rich libraries yield σ₀ < 0,
// AT-rich yield σ₀ > 0. ALWAYS emit alongside gc_interior and residualize across samples before
// interpreting as an oxidation proxy.
//
// length_slope: WLS of per-bin σ_l vs bin midpoint. Within a sample all fragments share burial
// time → per-base oxidation rate is length-independent → slope ≈ 0 under null. Nonzero slope
// indicates length-stratified composition or deamination-length coupling (confound flag, not
// oxidation axis).
//
// sigma_term: σ at pooled positions 0–2. σ_term ≫ σ₀ = terminal C→T contamination driving σ up
// at the ends (same confound that makes ε_term strongly negative).

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "taph/sample_damage_profile.hpp"
#include "taph/read_ancient_llr.hpp"

namespace taph {

namespace {
constexpr uint64_t MIN_INTERIOR = 200;
constexpr uint64_t MIN_TERM     = 50;
constexpr int      MIN_BINS_SLOPE = 4;
}

OxidationEstimate finalize_sigma(const SampleDamageProfile& profile) {
    OxidationEstimate out;

    // Pooled interior counts across all contributing bins
    uint64_t pT = 0, pG = 0, pA = 0, pC = 0;

    // WLS sufficient stats for length slope: σ_l vs bin midpoint
    double Sw = 0.0, Swx = 0.0, Swxx = 0.0, Swy = 0.0, Swxy = 0.0;
    int n_bins = 0;
    int last_bin_i = -1;  // index of highest contributing bin (approximates long-fragment composition)

    for (int i = 0; i < SampleDamageProfile::N_LEN_FINE; ++i) {
        const auto& b = profile.len_bins[i];
        const uint64_t dTG = b.t_interior + b.g_interior;
        const uint64_t dAC = b.a_interior + b.c_interior;
        if (dTG < MIN_INTERIOR || dAC < MIN_INTERIOR) continue;

        pT += b.t_interior;  pG += b.g_interior;
        pA += b.a_interior;  pC += b.c_interior;
        last_bin_i = i;
        ++n_bins;

        const double mid = static_cast<double>(SampleDamageProfile::LEN_FINE_MIN)
                         + static_cast<double>(SampleDamageProfile::LEN_FINE_W) * i
                         + 0.5 * static_cast<double>(SampleDamageProfile::LEN_FINE_W);
        const double s_l = static_cast<double>(b.t_interior) / static_cast<double>(dTG)
                         + static_cast<double>(b.a_interior) / static_cast<double>(dAC)
                         - 1.0;
        const double w   = static_cast<double>(std::min(dTG, dAC));

        Sw   += w;
        Swx  += w * mid;
        Swxx += w * mid * mid;
        Swy  += w * s_l;
        Swxy += w * mid * s_l;
    }

    if (n_bins == 0 || (pT + pG) == 0 || (pA + pC) == 0) return out;

    const double dTG_pool = static_cast<double>(pT + pG);
    const double dAC_pool = static_cast<double>(pA + pC);
    const double pTG      = static_cast<double>(pT) / dTG_pool;
    const double pAC      = static_cast<double>(pA) / dAC_pool;

    out.sigma0    = pTG + pAC - 1.0;
    // Delta-method SE: Var(T/(T+G)) ≈ p(1-p)/(T+G), same for A/(A+C); independent.
    out.sigma0_se = std::sqrt(pTG * (1.0 - pTG) / dTG_pool
                            + pAC * (1.0 - pAC) / dAC_pool);

    const uint64_t total = pT + pG + pA + pC;
    out.gc_interior = static_cast<double>(pG + pC) / static_cast<double>(total);
    out.n_counts    = total;
    out.n_bins      = n_bins;

    // sigma_long: σ of the longest contributing bin (approximates intrinsic composition).
    // delta_sigma = sigma0 - sigma_long: partial composition correction; >0 = excess GC→AT
    // beyond the longest-fragment baseline.
    if (last_bin_i >= 0) {
        const auto& lb = profile.len_bins[last_bin_i];
        const double dTGl = static_cast<double>(lb.t_interior + lb.g_interior);
        const double dACl = static_cast<double>(lb.a_interior + lb.c_interior);
        out.sigma_long  = static_cast<double>(lb.t_interior) / dTGl
                        + static_cast<double>(lb.a_interior) / dACl - 1.0;
        out.delta_sigma = out.sigma0 - out.sigma_long;
    }

    // Length slope (confound diagnostic)
    if (n_bins >= MIN_BINS_SLOPE) {
        const double det = Sw * Swxx - Swx * Swx;
        if (det > 0.0) {
            const double slope     = (Sw * Swxy - Swx * Swy) / det;
            const double intercept = (Swy - slope * Swx) / Sw;

            // Residual sum of squares for SE estimate
            double rss = 0.0;
            for (int i = 0; i < SampleDamageProfile::N_LEN_FINE; ++i) {
                const auto& b = profile.len_bins[i];
                const uint64_t dTG = b.t_interior + b.g_interior;
                const uint64_t dAC = b.a_interior + b.c_interior;
                if (dTG < MIN_INTERIOR || dAC < MIN_INTERIOR) continue;
                const double mid = static_cast<double>(SampleDamageProfile::LEN_FINE_MIN)
                                 + static_cast<double>(SampleDamageProfile::LEN_FINE_W) * i
                                 + 0.5 * static_cast<double>(SampleDamageProfile::LEN_FINE_W);
                const double s_l = static_cast<double>(b.t_interior) / static_cast<double>(dTG)
                                 + static_cast<double>(b.a_interior) / static_cast<double>(dAC)
                                 - 1.0;
                const double w   = static_cast<double>(std::min(dTG, dAC));
                const double r   = s_l - (intercept + slope * mid);
                rss += w * r * r;
            }
            const double sxx_w = Swxx - Swx * Swx / Sw;
            if (n_bins > 2 && sxx_w > 0.0) {
                out.length_slope    = slope;
                out.length_slope_se = std::sqrt((rss / static_cast<double>(n_bins - 2)) / sxx_w);
            }
        }
    }

    // Terminal σ at positions 0–2 (leakage QC)
    uint64_t tT = 0, tG = 0, tA = 0, tC = 0;
    for (const auto& b : profile.len_bins) {
        for (int p = 0; p < 3; ++p) {
            if (b.t_counts[p] + b.g_counts[p] < MIN_TERM) continue;
            tT += b.t_counts[p];
            tG += b.g_counts[p];
            tA += b.a_counts[p];
            tC += b.c_counts[p];
        }
    }
    if (tT + tG > 0 && tA + tC > 0) {
        out.sigma_term = static_cast<double>(tT) / static_cast<double>(tT + tG)
                       + static_cast<double>(tA) / static_cast<double>(tA + tC)
                       - 1.0;
    }

    out.fitted = true;
    return out;
}

} // namespace taph
