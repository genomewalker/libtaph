#pragma once
// Private helpers shared across damage_estimation_*.cpp TUs.
// Not part of the public taph API.
#include "taph/sample_damage_profile.hpp"
#include "taph/channel_count_table.hpp"
#include "taph/channel_registry.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace taph {
static inline double binomial_ll(double k, double n, double p) {
    if (n < 1 || p <= 0 || p >= 1) return 0.0;
    // Log-likelihood: k*log(p) + (n-k)*log(1-p) + constant (ignored for LLR)
    return k * std::log(p) + (n - k) * std::log(1.0 - p);
}

// Unclamped two-sample pooled-proportion z (terminal vs interior). The clamped binom_z used in the
// channel blocks is clamp(binom_z_raw, +/-kZCap); exposing the raw value lets the Layer-0 count
// table record pre_clamp_z so the golden gate sees z regressions the clamp would otherwise mask.
static inline double binom_z_raw(double k_t, double n_t, double k_i, double n_i) {
    if (n_t < 10.0 || n_i < 10.0) return std::numeric_limits<double>::quiet_NaN();
    double p_pool = (k_t + k_i) / (n_t + n_i);
    double var = p_pool * (1.0 - p_pool) * (1.0 / n_t + 1.0 / n_i);
    if (var < 1e-12) return std::numeric_limits<double>::quiet_NaN();
    return ((k_t / n_t) - (k_i / n_i)) / std::sqrt(var);
}

// The ONE clamped pooled-proportion z used by every stop-channel path (the primary pass AND the
// adapter-prefix-exclusion recompute): clamp(binom_z_raw, +/-kZCap) with NaN preserved.
static inline float binom_z_clamped(double k_t, double n_t, double k_i, double n_i) {
    double z = binom_z_raw(k_t, n_t, k_i, n_i);
    if (!std::isfinite(z)) return std::numeric_limits<float>::quiet_NaN();
    const double cap = static_cast<double>(SampleDamageProfile::kZCap);
    return static_cast<float>(std::clamp(z, -cap, cap));
}

// Layer-0 producer (file-scope so BOTH finalize_sample_profile and the adapter-prefix-exclusion
// recompute build IDENTICAL count tables from one code path): freezes the shadow-free 2x2 plus the
// exact numerator/denominator that fed the pooled z, the pre-clamp z, and the cap decision. The
// golden gate diffs these, so every emitted stop-channel z/rate must trace to a row produced here.
static inline StopChannelCountTable make_stop_count_table(const ChannelSpec& spec,
        double tp, double ts, double ip, double is, double tsh, double ish) {
    StopChannelCountTable ct;
    ct.channel_id = spec.channel_id;  ct.channel_type = spec.channel_type;
    ct.term_pre = tp;  ct.term_stop = ts;  ct.int_pre = ip;  ct.int_stop = is;
    ct.has_shadow = spec.has_deam_shadow;  ct.term_shadow = tsh;  ct.int_shadow = ish;
    ct.shadow_in_z = spec.shadow_in_z;  ct.shadow_in_rate = spec.shadow_in_rate;
    const bool shadow = spec.shadow_in_z;
    ct.z_num_term = ts;  ct.z_den_term = tp + ts + (shadow ? tsh : 0.0);
    ct.z_num_int  = is;  ct.z_den_int  = ip + is + (shadow ? ish : 0.0);
    ct.raw_rate_term = (tp + ts > 0.0) ? ts / (tp + ts) : 0.0;
    ct.raw_rate_int  = (ip + is > 0.0) ? is / (ip + is) : 0.0;
    ct.pre_clamp_z = binom_z_raw(ct.z_num_term, ct.z_den_term, ct.z_num_int, ct.z_den_int);
    // Cap policy is registry-driven: CLAMP_ZCAP (every current channel) == the prior hardcoded kZCap,
    // so this is byte-identical today, but a NONE channel would now be uncapped instead of silently 12.
    ct.z_cap = (spec.cap == CapPolicy::CLAMP_ZCAP)
               ? static_cast<double>(SampleDamageProfile::kZCap)
               : std::numeric_limits<double>::infinity();
    ct.cap_applied = std::isfinite(ct.pre_clamp_z) && std::abs(ct.pre_clamp_z) > ct.z_cap;
    ct.post_clamp_z = std::isfinite(ct.pre_clamp_z)
                      ? std::clamp(ct.pre_clamp_z, -ct.z_cap, ct.z_cap)
                      : ct.pre_clamp_z;
    ct.has_strata = spec.has_mh_stratification;
    return ct;
}

// Single-stratum 2x2 odds ratio (terminal stop vs interior stop) with a Haldane-Anscombe +0.5
// continuity correction for sparse cells. The four cells come straight from the Layer-0 count
// table, so the OR is reproducible from the gated counts. Unlike the pooled-Bernoulli z (which is
// inflated by correlated reads), an odds ratio is a descriptive effect size, so it is the primary
// statistic for the single-stratum channels (G/H); F uses the Mantel-Haenszel common OR instead.
static double stop_channel_or_haldane(const StopChannelCountTable& ct) {
    double a = ct.z_num_term;                    // terminal stop
    double b = ct.z_den_term - ct.z_num_term;    // terminal non-stop
    double c = ct.z_num_int;                     // interior stop
    double d = ct.z_den_int  - ct.z_num_int;     // interior non-stop
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5));
}

// LLR of exponential decay vs constant model over positions 1-9 (excludes pos 0 artifacts).
} // namespace taph
