// QC flags, preservation summary, context profile, oxoG estimates.
#include "taph/library_interpretation.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>

namespace taph {
// ── Library QC flags ──────────────────────────────────────────────────────────

LibraryQcFlags compute_library_qc_flags(
        const SampleDamageProfile& dp,
        bool is_ss,
        bool flag_hex_artifact,
        double jsd,
        double h_term,
        double short_read_frac) {
    LibraryQcFlags r;

    double tot_term = static_cast<double>(dp.n_hexamers_5prime);

    r.adapter_remnant_5prime   = dp.fit_offset_5prime > 1 || dp.position_0_artifact_5prime;
    r.adapter_remnant_3prime   = dp.fit_offset_3prime > 1 || dp.position_0_artifact_3prime;
    r.hexamer_composition_bias = tot_term > 0 && h_term < 4.5;
    r.hexamer_terminal_shift   = tot_term > 0 && dp.n_hexamers_interior > 0 && jsd > 0.05;
    r.short_read_spike         = short_read_frac > 0.20;
    r.depurination             = dp.depurination_detected;
    r.ds_3prime_signal_absent  = !is_ss
                                  && dp.d_max_5prime > 0.05f
                                  && dp.d_max_3prime < 0.03f;
    r.hexamer_artifact_bias    = flag_hex_artifact;

    // Inward-displaced G→A: pos0 depleted while pos1-3 elevated (trimmer artifact)
    if (!is_ss && dp.d_max_5prime > 0.05f) {
        constexpr double MIN_GA = 100.0;
        double bg_n = 0.0, bg_d = 0.0;
        for (int p = 8; p < 15; ++p) {
            bg_n += dp.a_freq_3prime[p];
            bg_d += dp.ag_total_3prime[p];
        }
        double bg = (bg_d >= 500.0) ? bg_n / bg_d : 0.47;
        double ga0 = (dp.ag_total_3prime[0] >= MIN_GA)
            ? dp.a_freq_3prime[0] / dp.ag_total_3prime[0] : -1.0;
        double ga_peak = 0.0;
        for (int p = 1; p <= 3; ++p)
            if (dp.ag_total_3prime[p] >= MIN_GA)
                ga_peak = std::max(ga_peak,
                                   dp.a_freq_3prime[p] / dp.ag_total_3prime[p]);
        if (ga0 >= 0.0)
            r.ga3_inward_displaced = ((ga0 - bg) < -0.02) && ((ga_peak - bg) > 0.02);
    }
    return r;
}

// ── Preservation / authenticity ───────────────────────────────────────────────

PreservationSummary compute_preservation_summary(
        const SampleDamageProfile& dp,
        bool   is_ss,
        bool   adapter_clipped,
        bool   flag_hex_artifact,
        double cpg_score_z,
        double oxog_score_z,
        double oxog_context_cosine,
        double hex_shift_p) {

    (void)oxog_score_z;
    (void)oxog_context_cosine;

    PreservationSummary r;

    auto clamp01 = [](double x) -> double { return std::clamp(x, 0.0, 1.0); };
    auto sat = [](double x, double half_sat) -> double {
        x = std::max(0.0, x);
        return x / (x + half_sat);
    };
    auto z_evidence = [](double z, double scale = 8.0) -> double {
        return 1.0 - std::exp(-std::max(0.0, z) / scale);
    };
    auto p_evidence = [](double p) -> double {
        if (std::isnan(p)) return 0.0;
        p = std::clamp(p, 1e-300, 1.0);
        double x = -std::log10(p);
        return x / (x + 3.0);
    };

    // d5_hexamer_corrected: when adapter artifact suppresses pos0 C→T,
    // extrapolate the decay curve back: d5_corr = d5 * exp(λ × fit_offset).
    double d5_raw  = static_cast<double>(dp.d_max_5prime);
    double d5_corr = d5_raw;
    if (!is_ss && !adapter_clipped && flag_hex_artifact
            && dp.position_0_artifact_5prime
            && dp.lambda_5prime > 0.0f && dp.fit_offset_5prime >= 1) {
        d5_corr = std::min(d5_raw * std::exp(static_cast<double>(dp.lambda_5prime)
                                              * dp.fit_offset_5prime), 0.50);
    }
    r.d5_raw               = d5_raw;
    r.d5_hexamer_corrected = d5_corr;
    r.d5_was_corrected     = (d5_corr != d5_raw);

    double a5 = sat(d5_corr, 0.05);
    double a3 = sat(static_cast<double>(dp.d_max_3prime), 0.05);
    double w5 = is_ss ? 1.0 : 2.0, w3 = is_ss ? 1.0 : 0.5;

    // Terminal-rate authenticity (mixture-free): the reference-free mixture pi floored
    // at H0, so its 0.80·pi branch reported ~0.87-authentic on true nulls. d_max is
    // already w_length-gated by finalize_dmax (SOLUTION_pi_delta_dmax.md §6).
    double authenticity_eff = clamp01((w5*a5 + w3*a3) / (w5 + w3));
    r.authenticity_eff = authenticity_eff;

    // authenticity_evidence
    double cpg_cov = static_cast<double>(dp.effcov_ct5_cpg_like_terminal) +
                     static_cast<double>(dp.effcov_ct5_noncpg_like_terminal);
    double wcpg = cpg_cov / (cpg_cov + 2000.0);

    double z5e = z_evidence(static_cast<double>(dp.terminal_z_5prime));
    double z3e = is_ss ? 0.0 : z_evidence(static_cast<double>(dp.terminal_z_3prime));
    // bulk authenticity evidence: validated terminal magnitude × length-coupling
    // (w_length), replacing the reference-free-non-identifiable mixture term (which
    // floored at H0). w_length maps 0.5→0 (null cluster) … 0.9→1 (SOLUTION §6).
    double w_coupling = std::clamp(
        (static_cast<double>(dp.bulk_damage.w_length) - 0.5) / 0.4, 0.0, 1.0);
    double bulk_e = sat(static_cast<double>(dp.d_max_combined), 0.05) * w_coupling;
    double cpg_e = wcpg * z_evidence(std::max(0.0, cpg_score_z), 5.0);
    double wbulk_e = dp.bulk_attempted ? 1.0 : 0.0;
    double auth_n_terms = 1.0 + (is_ss ? 0.0 : 1.0) + wbulk_e + wcpg;
    r.authenticity_evidence = (auth_n_terms > 1e-9)
        ? clamp01((z5e + z3e + bulk_e + cpg_e) / auth_n_terms)
        : 0.0;

    // Oxidation-like signal: reference-free, length/GC-stratified composition
    // contrast between deamination-weighted ancient-like and background-like reads.
    {
        const double signal = static_cast<double>(dp.oxidation_like_signal);
        const double control = static_cast<double>(dp.oxidation_like_control);
        const double adjusted = static_cast<double>(dp.oxidation_like_adjusted);
        const double excess = static_cast<double>(dp.oxidation_like_excess);
        const double signal_se = std::max(0.0, static_cast<double>(dp.oxidation_like_signal_se));
        const double control_se = std::max(0.0, static_cast<double>(dp.oxidation_like_control_se));
        const double adjusted_se = std::max(0.0, static_cast<double>(dp.oxidation_like_se));

        r.oxidation.raw_rate = {
            signal,
            signal - 1.96 * signal_se,
            signal + 1.96 * signal_se
        };
        r.oxidation.control_rate = {
            control,
            std::max(0.0, control - 1.96 * control_se),
            control + 1.96 * control_se
        };
        r.oxidation.excess_rate = {
            excess,
            std::max(0.0, adjusted - 1.96 * adjusted_se),
            std::max(0.0, adjusted + 1.96 * adjusted_se)
        };
        // Signed contrast (unclamped) so anti-oxidative vs genuine-zero vs
        // not-computed are distinguishable; excess_rate (clamped) is {0,0,0} for
        // all anti-oxidative libraries, an ambiguous sentinel on its own.
        r.oxidation.adjusted_rate = {
            adjusted,
            adjusted - 1.96 * adjusted_se,
            adjusted + 1.96 * adjusted_se
        };
        r.oxidation.excess_rate_computed = (adjusted > 0.0);
        r.oxidation.z_score = static_cast<double>(dp.oxidation_like_z);
        r.oxidation.bins_used = static_cast<double>(dp.oxidation_like_bins_used);
        r.oxidation.effective_bins = static_cast<double>(dp.oxidation_like_effective_bins);
        r.oxidation.heterogeneity = static_cast<double>(dp.oxidation_like_heterogeneity);
        r.oxidation.reliability_score = static_cast<double>(dp.oxidation_like_reliability);

        // BEHAVIORAL CHANGE (C5): "pass" previously certified only bin count and
        // artifact_suspect, so a strongly anti-oxidative library (z<0, excess=0)
        // earned "pass" — misread as "estimate is trustworthy / oxidation found".
        // Split the directional verdict out: when the channel is methodologically
        // operational but the contrast is anti-oxidative (z<=0), report "negative"
        // rather than "pass" (and rather than "fail", which would hide that the
        // channel was fully evaluated).
        const bool quality_ok = dp.oxidation_like_bins_used >= 6 &&
                                dp.oxidation_like_effective_bins >= 4.0f &&
                                !dp.oxidation_like_artifact_suspect;
        if (quality_ok && dp.oxidation_like_z > 0.0f)
            r.oxidation.reliability = "pass";
        else if (quality_ok)
            r.oxidation.reliability = "negative";
        else if (dp.oxidation_like_bins_used >= 2 && dp.oxidation_like_effective_bins >= 1.5f)
            r.oxidation.reliability = "warning";
        else
            r.oxidation.reliability = "fail";

        r.oxidation_eff = r.oxidation.excess_rate.estimate;
        r.oxidation_evidence = r.oxidation.reliability_score;
    }

    // qc_risk_eff / qc_evidence
    double d5_eff = static_cast<double>(dp.d_max_5prime);
    double d3_eff = static_cast<double>(dp.d_max_3prime);
    bool adapter_5 = (dp.fit_offset_5prime > 1) ||
                     (dp.position_0_artifact_5prime && d5_eff < 0.05 && !is_ss);
    bool adapter_3 = (dp.fit_offset_3prime > 1) ||
                     (dp.position_0_artifact_3prime && d3_eff < 0.05);
    double sig_eff     = is_ss ? d3_eff : d5_eff;
    double hex_dam_disc = sat(sig_eff, 0.05);
    double qhex  = sat(std::abs(static_cast<double>(dp.hexamer_excess_tc)), 0.05)
                   * (1.0 - hex_dam_disc);
    double qart  = dp.damage_artifact ? 1.0 : 0.0;
    double qadapt = (adapter_5 || adapter_3) ? 1.0 : 0.0;
    r.qc_risk_eff = clamp01((2.0*qhex + qart + qadapt) / 4.0);

    double hex_ev = qhex * p_evidence(hex_shift_p);
    r.qc_evidence = clamp01(std::max({hex_ev, qart, qadapt * 0.75}));

    // label
    if      (authenticity_eff < 0.10) r.label = "modern-like";
    else if (authenticity_eff < 0.30) r.label = "weak";
    else if (authenticity_eff < 0.55) r.label = "moderate";
    else                               r.label = "ancient";

    return r;
}

// ── Damage-context profile ────────────────────────────────────────────────────

namespace {

constexpr float kDeamNorm      = 0.10f;  // d_max at which deamination score ≈ 0.63
constexpr float kDipyrNorm     = 0.05f;  // dipyr contrast at which score saturates
                                          // (dp.dipyr_contrast uses 0.5 factors)
constexpr float kOxidativeNorm = 0.05f;  // oxidative signal at which score saturates
constexpr float kFragNorm      = 0.15f;  // purine enrichment at which score saturates
constexpr uint64_t kMinReads   = 1000;   // SampleDamageProfile::is_valid() threshold

inline float clamp01f(float x) {
    if (std::isnan(x)) return x;
    return x < 0.0f ? 0.0f : (x > 1.0f ? 1.0f : x);
}
inline float sigmoidf(float x) { return 1.0f / (1.0f + std::exp(-x)); }

} // namespace

const char* to_string(DamageContextProfile::DominantProcess p) {
    using D = DamageContextProfile::DominantProcess;
    switch (p) {
        case D::None:                   return "none";
        case D::LowDamage:              return "low_damage";
        case D::CytosineDeamination:    return "cytosine_deamination";
        case D::CpgEnrichedDeamination: return "cpg_enriched_deamination";
        case D::DipyrimidineBiased:     return "dipyrimidine_biased";
        case D::OxidativeLike:          return "oxidative_like";
        case D::FragmentationBias:      return "fragmentation_bias";
        case D::LibraryArtifactLikely:  return "library_artifact_likely";
    }
    return "none";
}

DamageContextProfile compute_damage_context_profile(
        const SampleDamageProfile& dp,
        double cpg_z,
        double hex_shift_z,
        bool   adapter_clipped,
        bool   adapter3_clipped,
        bool   flag_hex_artifact) {

    DamageContextProfile r;

    // Evidence block: raw underlying numbers for downstream auditing.
    r.evidence.d_max_5                    = dp.d_max_5prime;
    r.evidence.d_max_3                    = dp.d_max_3prime;
    r.evidence.lambda_5                   = dp.lambda_5prime;
    r.evidence.lambda_3                   = dp.lambda_3prime;
    r.evidence.log2_cpg_ratio             = dp.log2_cpg_ratio;
    r.evidence.cpg_z                      = static_cast<float>(cpg_z);
    r.evidence.ox_gt_asymmetry            = dp.ox_gt_asymmetry;
    r.evidence.purine_enrichment_5prime   = dp.purine_enrichment_5prime;
    r.evidence.hex_shift_z                = static_cast<float>(hex_shift_z);
    r.evidence.adapter_clipped            = adapter_clipped;
    r.evidence.adapter3_clipped           = adapter3_clipped;
    r.evidence.flag_hex_artifact          = flag_hex_artifact;
    r.evidence.position_0_artifact_5prime = dp.position_0_artifact_5prime;
    r.evidence.position_0_artifact_3prime = dp.position_0_artifact_3prime;
    r.evidence.fit_offset_5prime          = dp.fit_offset_5prime;
    r.evidence.fit_offset_3prime          = dp.fit_offset_3prime;
    r.evidence.n_reads                    = dp.n_reads;

    const auto NaN = std::numeric_limits<float>::quiet_NaN();

    // 8-oxoG 16-context panel summary (mean/max of the per-context signal).
    double s_sum = 0.0; float s_max = 0.0f; int n_ctx = 0;
    for (int i = 0; i < SampleDamageProfile::N_OXOG16; ++i) {
        float v = dp.s_oxog_16ctx[i];
        if (!std::isfinite(v)) continue;
        s_sum += v;
        if (v > s_max) s_max = v;
        ++n_ctx;
    }
    // NaN (not 0.0) when no context is evaluable, so "not computed" stays
    // distinguishable from a genuine zero signal. The isnan guard at the
    // oxidative_context_score step then correctly propagates NaN downstream.
    r.evidence.s_oxog_mean = n_ctx > 0 ? static_cast<float>(s_sum / n_ctx) : NaN;
    r.evidence.s_oxog_max  = n_ctx > 0 ? s_max : NaN;

    // Use the already-finalized contrast from SampleDamageProfile:
    //   dipyr_contrast = 0.5*(CC + TC) - 0.5*(AC + GC).
    r.evidence.dipyr_contrast = dp.dipyr_contrast;

    const bool any_art_flag = flag_hex_artifact || adapter_clipped ||
                              adapter3_clipped ||
                              dp.position_0_artifact_5prime ||
                              dp.position_0_artifact_3prime ||
                              dp.fit_offset_5prime > 1 ||
                              dp.fit_offset_3prime > 1;
    auto finite_max2 = [](float a, float b) {
        bool fa = std::isfinite(a), fb = std::isfinite(b);
        if (fa && fb) return std::max(a, b);
        if (fa)       return a;
        if (fb)       return b;
        return std::numeric_limits<float>::quiet_NaN();
    };

    // Six scores. Computed unconditionally so that low-coverage samples still
    // surface any evaluable raw signal; dominant_process is set to None below.

    // terminal deamination: 1 - exp(-max(d5, d3, d_combined) / kDeamNorm).
    // BEHAVIORAL CHANGE (C5): include the canonically rescue-corrected
    // d_max_combined so structurally-rescued libraries (e.g. marine DS) are not
    // under-scored; combined can only raise TDS above the raw terminal max,
    // never lower it (which would mis-flag confirmed-damage libraries). NaN-safe.
    float dmax_term = finite_max2(finite_max2(dp.d_max_5prime, dp.d_max_3prime),
                                  dp.d_max_combined);
    r.terminal_deamination_score = std::isnan(dmax_term)
        ? NaN
        : clamp01f(1.0f - std::exp(-dmax_term / kDeamNorm));

    // CpG context: sigmoid on the CpG/non-CpG z-score from compute_cpg_score.
    r.cpg_context_score = (!std::isfinite(dp.log2_cpg_ratio) || !std::isfinite(cpg_z))
        ? NaN
        : sigmoidf(static_cast<float>(cpg_z));

    // Dipyrimidine context: normalized upstream-context excess.
    r.dipyrimidine_context_score = !std::isfinite(r.evidence.dipyr_contrast)
        ? NaN
        : clamp01f(r.evidence.dipyr_contrast / kDipyrNorm);

    // Oxidative context: max of strand-asymmetric G->T signal and mean NGN panel.
    float ox_signed = finite_max2(std::fabs(dp.ox_gt_asymmetry),
                                  r.evidence.s_oxog_mean);
    r.oxidative_context_score = std::isnan(ox_signed)
        ? NaN
        : clamp01f(ox_signed / kOxidativeNorm);

    // Fragmentation context: purine enrichment at read starts vs interior.
    r.fragmentation_context_score = !std::isfinite(dp.purine_enrichment_5prime)
        ? NaN
        : clamp01f(dp.purine_enrichment_5prime / kFragNorm);

    // Library artifact: continuous sigmoid on the composition shift z. The
    // previous `any_art_flag ? 1.0f` floor saturated the score at 1.0 for nearly
    // every real ancient library (position-0 / adapter flags are near-universal),
    // making it useless as a [0,1] discriminator. The boolean flag presence is
    // surfaced separately in evidence (adapter_clipped, flag_hex_artifact,
    // position_0_artifact_*, fit_offset_*), so callers retain that signal; the
    // score now reflects artifact *magnitude*. If hex_shift_z is NaN and no flags
    // are set the signal is unknown — keep NaN rather than reading as 0.
    float art_sig  = std::isfinite(hex_shift_z)
        ? sigmoidf(static_cast<float>(hex_shift_z) - 4.0f)
        : NaN;
    if (std::isnan(art_sig)) {
        r.library_artifact_score = any_art_flag ? 0.0f : NaN;
    } else {
        r.library_artifact_score = clamp01f(art_sig);
    }

    // Deterministic dominant-process rule. Unevaluable signals (NaN) are not
    // allowed to trigger a branch; a NaN terminal-deamination score in
    // particular must not collapse to low_damage.
    using D = DamageContextProfile::DominantProcess;
    r.damage_status_present =
        (dp.damage_status == SampleDamageProfile::DamageStatus::PRESENT);
    float td  = r.terminal_deamination_score;
    float cpg = r.cpg_context_score;
    float dip = r.dipyrimidine_context_score;
    float ox  = r.oxidative_context_score;
    float fr  = r.fragmentation_context_score;
    auto gt = [](float x, float t){ return std::isfinite(x) && x >  t; };
    auto lt = [](float x, float t){ return std::isfinite(x) && x <  t; };

    if (dp.n_reads < kMinReads) {
        r.dominant_process = D::None;
        r.interpretation   = "insufficient coverage for damage-context scoring";
    } else if (std::isnan(td)) {
        r.dominant_process = D::None;
        r.interpretation   = "terminal deamination signal not evaluable";
    } else if (any_art_flag && lt(td, 0.5f) && !r.damage_status_present) {
        // Rule fires when a boolean flag is set (hex-artifact detector, adapter
        // stub, position-0 artifact, or fit-offset adapter remnant) AND the
        // terminal deamination signal is not dominating. The score itself
        // (library_artifact_score) saturates at 1.0 whenever a flag is present,
        // so `gt(art, 0.7)` would be tautological; the meaningful second gate
        // is "damage is not the primary driver", i.e. td < 0.5.
        // BEHAVIORAL CHANGE (C5): also require !damage_status_present so a
        // CI-confirmed-damage library (e.g. TDS just below 0.5 with adapter
        // flags) is not stamped library_artifact_likely while damage_status=
        // present — that was a direct internal contradiction. The damage label
        // then wins; artifact booleans remain visible in evidence.
        r.dominant_process = D::LibraryArtifactLikely;
        r.interpretation = "composition or adapter-stub evidence dominates over damage signal";
    } else if (gt(fr, 0.5f) && lt(td, 0.3f)) {
        r.dominant_process = D::FragmentationBias;
        r.interpretation = "purine enrichment at fragment starts without matching terminal deamination";
    } else if (lt(td, 0.10f)) {
        r.dominant_process = D::LowDamage;
        r.interpretation = "terminal deamination signal below detection threshold";
    } else if (gt(cpg, 0.7f) && gt(td, 0.3f)
               && gt(r.evidence.log2_cpg_ratio, 0.15f)) {
        // CpG z-scores can exceed the sigmoid threshold at modest effect sizes
        // on large samples; require a minimum log2 ratio so the label reflects
        // a real CpG-vs-non-CpG enrichment, not just statistical significance.
        r.dominant_process = D::CpgEnrichedDeamination;
        r.interpretation = "terminal deamination with elevated CpG-context contribution "
                           "consistent with methylated-cytosine deamination";
    } else if (gt(ox, 0.5f) && lt(td, 0.5f)) {
        r.dominant_process = D::OxidativeLike;
        r.interpretation = "strand-asymmetric G->T / C->A excess consistent with oxidative damage";
    } else if (gt(dip, 0.4f)) {
        r.dominant_process = D::DipyrimidineBiased;
        r.interpretation = "dipyrimidine upstream-context excess in terminal C->T rates";
    } else {
        r.dominant_process = D::CytosineDeamination;
        r.interpretation = "terminal C->T / G->A enrichment consistent with post-mortem cytosine deamination";
    }
    r.dominant_process_str = to_string(r.dominant_process);
    return r;
}


OxoGEstimate compute_oxog_estimate(const SampleDamageProfile& dp, bool is_ss) {
    OxoGEstimate r;

    r.ox_theta          = static_cast<double>(dp.g_bg_fitted);
    r.ox_theta_ci_lo    = static_cast<double>(dp.g_bg_fitted_ci_lo);
    r.ox_theta_ci_hi    = static_cast<double>(dp.g_bg_fitted_ci_hi);
    r.ox_theta_interior = static_cast<double>(dp.g_bg_interior_mean);
    r.fit_degenerate    = dp.g_fit_degenerate;

    r.ox_like_excess = static_cast<double>(dp.oxidation_like_excess);
    r.ox_like_z      = static_cast<double>(dp.oxidation_like_z);
    const double ox_se = std::max(0.0, static_cast<double>(dp.oxidation_like_se));
    const double ex    = std::max(0.0, static_cast<double>(dp.oxidation_like_adjusted));
    r.ox_like_ci_lo  = std::max(0.0, ex - 1.96 * ox_se);
    r.ox_like_ci_hi  = ex + 1.96 * ox_se;

    r.ox_uniformity_ratio = static_cast<double>(dp.ox_uniformity_ratio);

    if (!is_ss) {
        r.control_mode = "chargaff";
    } else if (dp.oxidation_like_bins_used >= 4 && dp.oxidation_like_effective_bins >= 3.0f) {
        r.control_mode = "trinuc";
    } else {
        r.control_mode = "auto";
    }

    const double asym = std::abs(static_cast<double>(dp.ox_gt_asymmetry));
    r.gc_skew_warning = (asym > 0.05 && dp.ox_gt_baseline < 0.01);

    return r;
}


// WLS regression helper: solve 5x5 system A*x = b via Gaussian elimination.
// Returns false if matrix is singular.
static bool solve5x5(double A[5][5], double b[5], double x[5]) {
    double M[5][6];
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) M[i][j] = A[i][j];
        M[i][5] = b[i];
    }
    for (int col = 0; col < 5; ++col) {
        int piv = col;
        for (int row = col + 1; row < 5; ++row)
            if (std::abs(M[row][col]) > std::abs(M[piv][col])) piv = row;
        if (std::abs(M[piv][col]) < 1e-14) return false;
        for (int j = 0; j <= 5; ++j) std::swap(M[col][j], M[piv][j]);
        for (int row = col + 1; row < 5; ++row) {
            double f = M[row][col] / M[col][col];
            for (int j = col; j <= 5; ++j) M[row][j] -= f * M[col][j];
        }
    }
    for (int row = 4; row >= 0; --row) {
        x[row] = M[row][5];
        for (int col = row + 1; col < 5; ++col) x[row] -= M[row][col] * x[col];
        x[row] /= M[row][row];
    }
    return true;
}

OxoTwoMarkerResult compute_oxo_two_marker(const SampleDamageProfile& dp, bool is_ss) {
    using Bins = SampleDamageProfile::OxoTwoMarkerBins;
    OxoTwoMarkerResult r;
    r.consistency_basis = is_ss ? OxoConsistencyBasis::SS_END_SYMMETRY
                                : OxoConsistencyBasis::DS_STRAND_SYMMETRY;
    const auto& panel = is_ss ? dp.oxo_two_marker_ss : dp.oxo_two_marker;

    // Bin midpoints for predictors
    // s1/s2: bins 0-3 → normalised rate 0/3, 1/3, 2/3, 3/3
    // GC bins: 0-3 → mid 12.5%, 37.5%, 62.5%, 87.5%
    // len bins: 0-3 → mid 37.5, 57.5, 90, 130
    static const double s_mid[4]   = {0.0, 1.0/3, 2.0/3, 1.0};
    static const double gc_mid[4]  = {0.125, 0.375, 0.625, 0.875};
    static const double len_mid[4] = {37.5,  57.5,  90.0, 130.0};

    // Accumulate WLS sufficient statistics
    double XtWX[5][5] = {};
    double XtWy[5]    = {};
    int    n_cells    = 0;

    for (int s1 = 0; s1 < Bins::N_S; ++s1) {
    for (int s2 = 0; s2 < Bins::N_S; ++s2) {
    for (int gc = 0; gc < Bins::N_GC; ++gc) {
    for (int l  = 0; l  < Bins::N_L;  ++l) {
        const auto& c = panel.cells[Bins::idx(s1, s2, gc, l)];
        if (c.sum_nGT < 200 || c.sum_nAC < 200) continue;

        const double D = static_cast<double>(c.sum_T)   / c.sum_nGT
                       - static_cast<double>(c.sum_A)   / c.sum_nAC;
        // Harmonic-mean weight (effective sample size)
        const double w = 2.0 * c.sum_nGT * c.sum_nAC
                       / static_cast<double>(c.sum_nGT + c.sum_nAC);

        const double x[5] = {1.0, s_mid[s1], s_mid[s2], gc_mid[gc], len_mid[l]};
        for (int i = 0; i < 5; ++i) {
            XtWy[i] += w * x[i] * D;
            for (int j = 0; j < 5; ++j)
                XtWX[i][j] += w * x[i] * x[j];
        }
        ++n_cells;
    }}}}

    r.n_cells_used = n_cells;
    if (n_cells < 6) return r;  // under-determined

    double beta[5] = {};
    if (!solve5x5(XtWX, XtWy, beta)) return r;

    r.alpha  = beta[0];
    r.beta1  = beta[1];
    r.beta2  = beta[2];

    // Residual variance
    double sse = 0.0, sw = 0.0;
    for (int s1 = 0; s1 < Bins::N_S; ++s1) {
    for (int s2 = 0; s2 < Bins::N_S; ++s2) {
    for (int gc = 0; gc < Bins::N_GC; ++gc) {
    for (int l  = 0; l  < Bins::N_L;  ++l) {
        const auto& c = panel.cells[Bins::idx(s1, s2, gc, l)];
        if (c.sum_nGT < 200 || c.sum_nAC < 200) continue;
        const double D = static_cast<double>(c.sum_T) / c.sum_nGT
                       - static_cast<double>(c.sum_A) / c.sum_nAC;
        const double w = 2.0 * c.sum_nGT * c.sum_nAC
                       / static_cast<double>(c.sum_nGT + c.sum_nAC);
        const double x[5] = {1.0, s_mid[s1], s_mid[s2], gc_mid[gc], len_mid[l]};
        double pred = 0.0;
        for (int i = 0; i < 5; ++i) pred += beta[i] * x[i];
        sse += w * (D - pred) * (D - pred);
        sw  += w;
    }}}}
    r.sigma2 = (n_cells > 5) ? sse / (n_cells - 5) : 0.0;

    // SE of beta1 and beta2: solve XtWX * v = e_j; SE²(j) = sigma2 * v[j]
    for (int j = 1; j <= 2; ++j) {
        double ej[5] = {}, vj[5] = {};
        ej[j] = 1.0;
        if (solve5x5(XtWX, ej, vj)) {
            double se = std::sqrt(std::max(0.0, r.sigma2 * vj[j]));
            if (j == 1) { r.beta1_se = se; r.beta1_z = (se > 0) ? r.beta1 / se : 0.0; }
            else        { r.beta2_se = se; r.beta2_z = (se > 0) ? r.beta2 / se : 0.0; }
        }
    }

    r.delta_beta = r.beta1 - r.beta2;
    // B1: propagated SE of delta_beta (independent-fit approximation). Same quantity the DS
    // consistency test below already forms inline; surfaced as a field for the downstream gate.
    double se_delta = std::sqrt(r.beta1_se * r.beta1_se + r.beta2_se * r.beta2_se);
    r.delta_beta_se = se_delta;
    // Library-appropriate consistency test:
    //   DS: magnitude equality ||beta1|-|beta2|| < 2*se_delta.
    //       beta2 is structurally negative (G→A couples with opposite sign), so the
    //       test compares absolute strengths. Catches strand-asymmetric artefacts.
    //   SS: both markers independently significant at z > 3.0 AND same sign (both > 0).
    //       5' C→T and 3' C→T have systematically different magnitudes (~0.72 ratio)
    //       so magnitude equality is the wrong criterion; independent significance of
    //       both deamination channels is the physically meaningful check.
    if (!is_ss) {
        r.markers_consistent =
            (se_delta > 0 &&
             std::abs(std::abs(r.beta1) - std::abs(r.beta2)) < 2.0 * se_delta);
    } else {
        constexpr double kSSZThreshold = 3.0;
        r.markers_consistent = (r.beta1_z > kSSZThreshold &&
                                 r.beta2_z > kSSZThreshold &&
                                 r.beta1 > 0.0 && r.beta2 > 0.0);
    }
    r.valid = true;
    return r;
}

} // namespace taph
