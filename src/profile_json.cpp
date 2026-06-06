#include "taph/profile_json.hpp"
#include "taph/library_interpretation.hpp"
#include "taph/length_bin_damage_profile.hpp"
#include "taph/channel_registry.hpp"
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>

namespace taph {

static constexpr double kMinCov          = 100.0;
static constexpr float  kOxChannelZDetect = 3.0f;
// C2: emitted-statistic clamp. z/G-test stats built from per-read/per-site counts as
// independent Bernoulli trials scale ~sqrt(N) (correlated reads), reaching hundreds–thousands.
// Detection gates use z>3, well inside ±12, so clamping the EMITTED value is behaviour-neutral.
static inline double clamp_z(double z) {
    if (!std::isfinite(z)) return z;       // NaN/Inf handled by nan_or at emission
    const double cap = static_cast<double>(SampleDamageProfile::kZCap);  // single source of truth
    return z < -cap ? -cap : (z > cap ? cap : z);
}

static const char* libtype_cstr(SampleDamageProfile::LibraryType t) {
    switch (t) {
        case SampleDamageProfile::LibraryType::DOUBLE_STRANDED: return "DOUBLE_STRANDED";
        case SampleDamageProfile::LibraryType::SINGLE_STRANDED: return "SINGLE_STRANDED";
        default: return "UNKNOWN";
    }
}

static const char* libtype_human(SampleDamageProfile::LibraryType t) {
    switch (t) {
        case SampleDamageProfile::LibraryType::DOUBLE_STRANDED: return "double-stranded";
        case SampleDamageProfile::LibraryType::SINGLE_STRANDED: return "single-stranded";
        default: return "unknown";
    }
}

template<typename T>
static std::string nan_or(T v) {
    // C1: null any non-finite sentinel — NaN AND ±Inf. std::isnan(Inf)==false would
    // otherwise emit a bare `inf`/`-inf` token (invalid JSON) for any of the 11 call sites.
    return std::isfinite(static_cast<double>(v)) ? std::to_string(v) : "null";
}

static std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s) {
        if (c == '"')  out += "\\\"";
        else if (c == '\\') out += "\\\\";
        else out += c;
    }
    return out;
}

void profile_to_json(const SampleDamageProfile& dp,
                     std::ostream& j,
                     const ProfileJsonInput& in) {
    using SP = SampleDamageProfile;
    constexpr int N_POS = SP::N_POS;

    const bool is_ss = (dp.library_type == SP::LibraryType::SINGLE_STRANDED);

    // ── Pre-compute CircLigase selection-bias diagnostics ─────────────────────
    // These must be available before the deamination block (line ~148) and the
    // diagnostic_groups block (~line 1076), so they are computed once here.
    const HexEndAsymmetry& hea_early = in.hex_end_asymmetry;
    const bool d5_suppressed_early    = dp.d5_profile_flat && (dp.d_max_3prime >= 0.05f);
    const bool d3_selection_biased_early = is_ss
                                        && (hea_early.rc_overlap_topk == 0)
                                        && (dp.d_max_3prime > dp.d_max_5prime * 1.5f)
                                        && (dp.d_max_5prime < 0.02f)
                                        && !d5_suppressed_early;
    float d_max_combined_out = dp.d_max_combined;
    float d_metamatch_out     = dp.d_metamatch;
    const char* source_str_out = dp.d_max_source_str();
    if (d3_selection_biased_early) {
        d_max_combined_out = dp.d_max_5prime;
        source_str_out = "5prime_conservative_ss";
        // C5: dp.d_metamatch was computed in damage_estimation.cpp from the pre-override
        // d_max_combined; mirror the CircLigase override here so d_max_combined and d_metamatch
        // stay internally consistent in the JSON. [behavioral change to emitted d_metamatch]
        d_metamatch_out = dp.d_max_5prime;
    }

    // ── Pre-compute all derived scores ────────────────────────────────────────
    auto cpg     = compute_cpg_score(dp);
    auto oxog_is = compute_oxog_interior_score(dp);
    auto hs      = compute_hex_stats(dp);
    auto ds      = compute_depur_score(dp, is_ss);
    auto otr     = compute_oxog_trinuc(dp);
    auto oxe     = compute_oxog_estimate(dp, is_ss);
    auto otm     = compute_oxo_two_marker(dp, is_ss);
    auto pres    = compute_preservation_summary(dp, is_ss,
                       in.adapter_clipped, in.flag_hex_artifact,
                       cpg.z, oxog_is.z, otr.cosine, hs.shift_p);
    auto dcp     = compute_damage_context_profile(dp,
                       cpg.z, hs.shift_z,
                       in.adapter_clipped, in.adapter3_clipped,
                       in.flag_hex_artifact);


    // Artifact reasons
    std::vector<const char*> artifact_reasons;
    if (dp.position_0_artifact_5prime || dp.fit_offset_5prime > 1)
        artifact_reasons.push_back("adapter_remnant_5prime");
    if (dp.position_0_artifact_3prime || dp.fit_offset_3prime > 1)
        artifact_reasons.push_back("adapter_remnant_3prime");
    if (in.flag_hex_artifact)
        artifact_reasons.push_back("hexamer_artifact_bias");
    if (dp.position_0_artifact_5prime)
        artifact_reasons.push_back("position0_artifact_5prime");
    if (dp.position_0_artifact_3prime)
        artifact_reasons.push_back("position0_artifact_3prime");
    if (dp.inverted_pattern_5prime)
        artifact_reasons.push_back("inverted_pattern_5prime");
    if (dp.inverted_pattern_3prime)
        artifact_reasons.push_back("inverted_pattern_3prime");

    j << std::fixed;

    // ── Top-level ─────────────────────────────────────────────────────────────
    j << "{\n";
    j << "  \"input\": \"" << in.sample_name << "\",\n";
    j << "  \"n_reads\": " << in.n_reads << ",\n";
    j << "  \"library_type\": \"" << dp.library_type_str() << "\",\n";
    j << "  \"library_bic_call\": \"" << libtype_cstr(dp.library_bic_call) << "\",\n";  // Wave-2: frozen pure-BIC call
    j << "  \"library_call_overridden\": "
      << ((dp.library_bic_call != SampleDamageProfile::LibraryType::UNKNOWN
           && dp.library_bic_call != dp.library_type) ? "true" : "false") << ",\n";
    // Wave-2b: the override ledger -- one {rule, from, to} record per post-hoc rescue/veto
    // that moved library_type away from library_bic_call (empty [] when BIC stands).
    j << "  \"library_overrides\": [";
    for (size_t i = 0; i < dp.library_overrides.size(); ++i) {
        const auto& ov = dp.library_overrides[i];
        if (i) j << ", ";
        j << "{\"rule\": \"" << ov.rule_id
          << "\", \"from\": \"" << libtype_human(ov.from)
          << "\", \"to\": \"" << libtype_human(ov.to) << "\"}";
    }
    j << "],\n";
    j << "  \"library_type_auto\": " << (dp.library_type_auto_detected ? "true" : "false") << ",\n";
    j << "  \"library_type_rescued\": " << (dp.library_type_rescued ? "true" : "false") << ",\n";
    j << "  \"library_type_evaluable\": " << (dp.library_type_evaluable ? "true" : "false") << ",\n";
    j << "  \"library_p_ds\": " << std::setprecision(6) << dp.library_p_ds << ",\n";
    j << "  \"library_p_ss\": " << dp.library_p_ss << ",\n";
    j << "  \"library_p_bias\": " << dp.library_p_bias << ",\n";
    j << "  \"library_p_winner\": " << dp.library_p_winner << ",\n";
    j << "  \"library_p_ds_final\": " << dp.library_p_ds_final << ",\n";
    j << "  \"library_p_ss_final\": " << dp.library_p_ss_final << ",\n";
    j << "  \"library_p_bias_final\": " << dp.library_p_bias_final << ",\n";
    j << "  \"library_p_winner_final\": " << dp.library_p_winner_final << ",\n";
    j << "  \"library_auto_type\": \"" << libtype_cstr(dp.library_auto_type) << "\",\n";
    j << "  \"library_auto_evaluable\": " << (dp.library_auto_evaluable ? "true" : "false") << ",\n";
    j << "  \"library_forced_type\": \"" << libtype_cstr(dp.library_forced_type) << "\",\n";
    j << "  \"library_bic_winner_model\": \"" << dp.library_bic_winner_model << "\",\n";
    j << "  \"library_bic_winner_model_class\": \"" << dp.library_bic_winner_model_class << "\",\n";
    j << "  \"library_p_is_pre_override\": " << (dp.library_p_is_pre_override ? "true" : "false") << ",\n";
    j << "  \"library_bic_second_model\": \"" << dp.library_bic_second_model << "\",\n";
    j << "  \"library_bic_margin\": " << std::setprecision(2) << dp.library_bic_margin << ",\n";
    j << "  \"library_p_ds_class_min\": " << std::setprecision(6) << dp.library_p_ds_class_min << ",\n";
    j << "  \"library_p_ss_class_min\": " << dp.library_p_ss_class_min << ",\n";
    j << "  \"library_p_bias_class_min\": " << dp.library_p_bias_class_min << ",\n";
    j << "  \"library_artifact_contaminated\": " << (artifact_reasons.empty() ? "false" : "true") << ",\n";
    j << "  \"library_artifact_reasons\": [";
    for (size_t i = 0; i < artifact_reasons.size(); ++i) {
        if (i) j << ",";
        j << "\"" << artifact_reasons[i] << "\"";
    }
    j << "],\n";
    j << "  \"library_M_SS_orig_active\": " << (dp.library_M_SS_orig_active ? "true" : "false") << ",\n";
    j << "  \"library_M_SS_asym_active\": " << (dp.library_M_SS_asym_active ? "true" : "false") << ",\n";
    j << "  \"library_spike_is_ss\": " << (dp.library_spike_is_ss ? "true" : "false") << ",\n";
    j << "  \"library_spike_gate_ga0_amp\": " << std::setprecision(6) << dp.library_spike_gate_ga0_amp << ",\n";
    j << "  \"library_spike_gate_structural_bilateral\": "
      << (dp.library_spike_gate_structural_bilateral ? "true" : "false") << ",\n";
    j << "  \"library_joint_lambda_restricted\": "
      << (dp.library_joint_lambda_restricted ? "true" : "false") << ",\n";
    {
        const char* m = "neutral";
        switch (dp.library_type) {
            case SP::LibraryType::DOUBLE_STRANDED: m = "ds_ga"; break;
            case SP::LibraryType::SINGLE_STRANDED: m = "ss_ct"; break;
            default: break;
        }
        j << "  \"damage_3prime_mode\": \"" << m << "\",\n";
    }
    const char* ds_str =
        (dp.damage_status == SP::DamageStatus::PRESENT) ? "present" :
        (dp.damage_status == SP::DamageStatus::WEAK)    ? "weak"    : "absent";
    j << "  \"damage_status\": \"" << ds_str << "\",\n";

    // ── Deamination ───────────────────────────────────────────────────────────
    j << "  \"deamination\": {\n";
    j << "    \"d_max_5prime\": " << std::setprecision(6) << dp.d_max_5prime << ",\n";
    j << "    \"d_max_3prime\": " << dp.d_max_3prime << ",\n";
    j << "    \"d_max_combined\": " << d_max_combined_out << ",\n";
    j << "    \"d_metamatch\": " << d_metamatch_out << ",\n";
    j << "    \"source\": \"" << source_str_out << "\",\n";
    j << "    \"lambda_5prime\": " << (dp.lambda_5prime_fitted ? std::to_string(dp.lambda_5prime) : "null") << ",\n";
    j << "    \"lambda_3prime\": " << (dp.lambda_3prime_fitted ? std::to_string(dp.lambda_3prime) : "null") << ",\n";
    j << "    \"bg_5prime\": " << dp.fit_baseline_5prime << ",\n";
    j << "    \"bg_3prime\": " << dp.fit_baseline_3prime << ",\n";
    j << "    \"validated\": " << (dp.damage_validated ? "true" : "false") << ",\n";
    j << "    \"artifact\": " << (dp.damage_artifact ? "true" : "false") << ",\n";
    j << "    \"joint\": {\n";
    // C2: ΔBIC and Wald z both scale ~O(N) on correlated reads — exploratory, not a
    // calibrated test. Clamp the EMITTED magnitude (ΔBIC→[-200,200], z→±kZCap) and flag
    // saturation; p_damage uses a stable internal logistic on clamped ΔBIC/2, unaffected.
    {
        const double dbic = static_cast<double>(dp.joint_delta_bic);
        const bool dbic_sat = std::isfinite(dbic) && std::abs(dbic) > 200.0;
        const double dbic_emit = !std::isfinite(dbic) ? dbic
                               : (dbic < -200.0 ? -200.0 : (dbic > 200.0 ? 200.0 : dbic));
        const bool zdel_cap = dp.joint_z_delta_capped;
        j << "      \"delta_bic\": " << nan_or(dbic_emit) << ",\n";
        j << "      \"delta_bic_saturated\": " << (dbic_sat ? "true" : "false") << ",\n";
        j << "      \"z_delta\": " << nan_or(clamp_z(dp.joint_z_delta)) << ",\n";
        j << "      \"z_delta_capped\": " << (zdel_cap ? "true" : "false") << ",\n";
    }
    // C1: p_damage and n_informative_positions are reset-value zeros when valid=false;
    // emit null so "not computed" is not read as a genuine probability/count estimate.
    j << "      \"p_damage\": " << (dp.joint_model_valid ? std::to_string(dp.joint_p_damage) : "null") << ",\n";
    j << "      \"n_informative_positions\": " << (dp.joint_model_valid ? std::to_string(dp.joint_n_informative) : "null") << ",\n";
    j << "      \"valid\": " << (dp.joint_model_valid ? "true" : "false") << "\n";
    j << "    },\n";
    const bool mix_gc_identified = dp.mixture_identifiable && dp.mixture_converged;
    j << "    \"mixture_gc\": {\n";
    j << "      \"status\": \"" << (mix_gc_identified ? "identified" : "undetermined") << "\",\n";
    if (mix_gc_identified) {
        j << "      \"d_ancient\": " << dp.mixture_d_ancient << ",\n";
        j << "      \"pi_ancient\": " << dp.mixture_pi_ancient << ",\n";
        j << "      \"d_population_highgc\": " << dp.mixture_d_population_highgc << ",\n";
    } else {
        j << "      \"d_ancient\": null,\n";
        j << "      \"pi_ancient\": null,\n";
        j << "      \"d_population_highgc\": null,\n";
    }
    j << "      \"identifiable\": " << (dp.mixture_identifiable ? "true" : "false") << ",\n";
    j << "      \"converged\": " << (dp.mixture_converged ? "true" : "false") << ",\n";
    j << "      \"n_components\": " << dp.mixture_n_components << "\n";
    j << "    },\n";
    j << "    \"cpg_like\": {\n";
    j << "      \"dmax_ct5_cpg\": "           << nan_or(dp.dmax_ct5_cpg_like)    << ",\n";
    j << "      \"dmax_ct5_noncpg\": "        << nan_or(dp.dmax_ct5_noncpg_like) << ",\n";
    j << "      \"cpg_ratio\": "              << nan_or(dp.cpg_ratio)            << ",\n";
    j << "      \"log2_cpg_ratio\": "         << nan_or(dp.log2_cpg_ratio)       << ",\n";
    j << "      \"cpg_ratio_backwards\": "    << (dp.cpg_ratio_backwards ? "true" : "false") << ",\n";  // P5 QC
    j << "      \"baseline_cpg\": "           << nan_or(dp.fit_baseline_ct5_cpg_like)    << ",\n";
    j << "      \"baseline_noncpg\": "        << nan_or(dp.fit_baseline_ct5_noncpg_like) << ",\n";
    j << "      \"cov_terminal_cpg\": "       << std::setprecision(0) << dp.cov_ct5_cpg_like_terminal    << ",\n";
    j << "      \"cov_terminal_noncpg\": "    << dp.cov_ct5_noncpg_like_terminal << ",\n";
    j << "      \"effcov_terminal_cpg\": "    << dp.effcov_ct5_cpg_like_terminal    << ",\n";
    j << "      \"effcov_terminal_noncpg\": " << dp.effcov_ct5_noncpg_like_terminal << ",\n";
    // C1: when log2_cpg_ratio is NaN (degenerate noncpg fit) compute_cpg_score returns its
    // zero-init {z=0,p=1} default — emit null so "not computed" is not read as z=0/p=1.
    j << "      \"cpg_score_z\": " << std::setprecision(6) << nan_or(clamp_z(cpg.z)) << ",\n";  // P4: capped to +/-kZCap (uniform with the channel z's; detection saturates well inside)
    j << "      \"cpg_score_p\": " << nan_or(cpg.p) << "\n";
    j << "    },\n";
    j << "    \"context_deamination\": {\n";
    j << "      \"dmax_AC\": " << nan_or(dp.dmax_ct5_by_upstream[SP::CTX_AC]) << ",\n";
    j << "      \"dmax_CC\": " << nan_or(dp.dmax_ct5_by_upstream[SP::CTX_CC]) << ",\n";
    j << "      \"dmax_GC\": " << nan_or(dp.dmax_ct5_by_upstream[SP::CTX_GC]) << ",\n";
    j << "      \"dmax_TC\": " << nan_or(dp.dmax_ct5_by_upstream[SP::CTX_TC]) << ",\n";
    j << "      \"dipyr_contrast\": " << nan_or(dp.dipyr_contrast) << ",\n";

    // C1: emit null when context test was not run (valid_ctx_count<4 or mean_d<=0.001)
    j << "      \"heterogeneity_chi2\": " << (dp.context_heterogeneity_computed ? nan_or(dp.context_heterogeneity_chi2) : std::string("null")) << ",\n";
    j << "      \"heterogeneity_p\": " << (dp.context_heterogeneity_computed ? nan_or(dp.context_heterogeneity_p) : std::string("null")) << ",\n";
    j << "      \"heterogeneity_detected\": " << (dp.context_heterogeneity_computed ? (dp.context_heterogeneity_detected ? "true" : "false") : "null") << "\n";
    j << "    },\n";
    j << "    \"per_pos_5prime_ct\": [";
    for (int p = 0; p < N_POS; ++p) {
        // C1: emit null (not -1.0) when below kMinCov — a probability can't be negative.
        if (dp.tc_total_5prime[p] >= kMinCov) j << std::setprecision(6) << dp.t_freq_5prime[p];
        else                                  j << "null";
        if (p < N_POS - 1) j << ",";
    }
    j << "],\n";
    j << "    \"per_pos_3prime\": [";
    for (int p = 0; p < N_POS; ++p) {
        // C1: emit null when coverage is below kMinCov instead of the -1.0 sentinel.
        bool ok = is_ss ? (dp.tc_total_3prime[p] >= kMinCov)
                        : (dp.ag_total_3prime[p] >= kMinCov);
        if (!ok) j << "null";
        else if (is_ss) j << std::setprecision(6) << (dp.t_freq_3prime[p] / dp.tc_total_3prime[p]);
        else            j << std::setprecision(6) << dp.a_freq_3prime[p];
        if (p < N_POS - 1) j << ",";
    }
    j << (in.lsd && !in.lsd->bins.empty() ? "],\n" : "]\n");

    // by_length
    if (in.lsd && !in.lsd->bins.empty()) {
        const auto& lsd = *in.lsd;
        j << "    \"by_length_method\": \"" << lsd.method << "\",\n";
        j << "    \"by_length\": [";
        for (size_t b = 0; b < lsd.bins.size(); ++b) {
            const auto& lb = lsd.bins[b];
            if (b > 0) j << ",";
            j << "\n      {";
            j << "\"length_lo\":" << lb.length_lo
              << ",\"length_hi\":" << lb.length_hi
              << ",\"n_reads\":" << lb.n_reads
              << ",\"d_max_5prime\":" << std::setprecision(6) << lb.d_max_5prime
              << ",\"d_max_3prime\":" << lb.d_max_3prime
              << ",\"lambda_5prime\":" << lb.lambda_5prime
              << ",\"lambda_3prime\":" << lb.lambda_3prime
              << ",\"bg_5prime\":" << lb.bg_5prime
              << ",\"bg_3prime\":" << lb.bg_3prime

              << ",\"validated\":" << (lb.validated ? "true" : "false")
              << ",\"ss_mode\":" << (lb.ss_mode ? "true" : "false")
              << ",\"source\":\"" << lb.source << "\"";
            j << ",\"per_pos_5prime_ct\":[";
            for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                j << std::setprecision(6) << lb.per_pos_5prime_ct[p];
                if (p + 1 < LengthBinDamageProfile::N_POS) j << ",";
            }
            j << "],\"per_pos_3prime\":[";
            for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                j << std::setprecision(6) << lb.per_pos_3prime[p];
                if (p + 1 < LengthBinDamageProfile::N_POS) j << ",";
            }
            j << "]";
            const bool mix_lb_identified = lb.mixture_identifiable && lb.mixture_converged;
            j << ",\"mixture\":{\"status\":\"" << (mix_lb_identified ? "identified" : "undetermined") << "\"";
            if (mix_lb_identified) {
                j << ",\"d_ancient\":" << std::setprecision(6) << lb.mixture_d_damaged
                  << ",\"d_population_highgc\":" << lb.mixture_d_population_highgc
                  << ",\"d_population\":" << lb.mixture_d_population
                  << ",\"pi_ancient\":" << lb.mixture_pi_damaged;
            } else {
                j << ",\"d_ancient\":null,\"d_population_highgc\":null,\"d_population\":null,\"pi_ancient\":null";
            }
            j << ",\"n_components\":" << lb.mixture_n_components
              << ",\"converged\":" << (lb.mixture_converged ? "true" : "false")
              << ",\"identifiable\":" << (lb.mixture_identifiable ? "true" : "false")
              << "}";
            j << ",\"gc_d_max\":[";
            for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
                // C1: -1.0 is the insufficient-coverage sentinel; d_max is in [0,1] → null.
                if (lb.gc_d_max[g] < 0.0) j << "null";
                else j << std::setprecision(6) << lb.gc_d_max[g];
                if (g + 1 < LengthBinDamageProfile::N_GC_BINS) j << ",";
            }
            j << "],\"gc_n_reads\":[";
            for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
                j << lb.gc_n_reads[g];
                if (g + 1 < LengthBinDamageProfile::N_GC_BINS) j << ",";
            }
            j << "],\"gc_p_damaged\":[";
            for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
                // C1: -1.0 is the insufficient-coverage sentinel (gb.valid==false). A
                // probability in [0,1] is never legitimately negative → emit null.
                if (lb.gc_p_damaged[g] < 0.0f) j << "null";
                else j << std::setprecision(6) << lb.gc_p_damaged[g];
                if (g + 1 < LengthBinDamageProfile::N_GC_BINS) j << ",";
            }
            j << "],\"gc_per_pos_ct\":[";
            for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
                j << "[";
                for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                    double v = lb.gc_per_pos_ct[g][p];
                    if (v < 0.0) j << "null"; else j << std::setprecision(6) << v;
                    if (p + 1 < LengthBinDamageProfile::N_POS) j << ",";
                }
                j << "]";
                if (g + 1 < LengthBinDamageProfile::N_GC_BINS) j << ",";
            }
            j << "]";
            auto dnan = [](double v) -> std::string {
                // C1: null NaN AND ±Inf — a bare inf/-inf token is invalid JSON (isnan(Inf)==false).
                return std::isfinite(v) ? std::to_string(v) : "null";
            };
            auto write_dnan_arr = [&](const std::array<double, LengthBinDamageProfile::N_POS>& arr) {
                for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                    j << dnan(arr[p]);
                    if (p + 1 < LengthBinDamageProfile::N_POS) j << ",";
                }
            };
            j << ",\"llr\":{\"n_damaged\":" << lb.n_damaged
              << ",\"n_undamaged\":" << lb.n_undamaged
              << ",\"d_max_5_damaged\":"          << dnan(lb.d_max_5_damaged)
              << ",\"d_max_3_damaged\":"          << dnan(lb.d_max_3_damaged)
              << ",\"d_max_5_cpg_damaged\":"      << dnan(lb.d_max_5_cpg_damaged)
              << ",\"d_max_5_noncpg_damaged\":"   << dnan(lb.d_max_5_noncpg_damaged)
              << ",\"s_gt_5_damaged_vs_undamaged\":" << dnan(lb.s_gt_5_damaged_vs_undamaged)
              << ",\"g_to_t_5_damaged\":"         << dnan(lb.g_to_t_5_damaged)
              << ",\"pG_terminal_5_damaged\":"    << dnan(lb.pG_terminal_5_damaged)
              << ",\"pG_interior_5_damaged\":"    << dnan(lb.pG_interior_5_damaged)
              << ",\"per_pos_5prime_ct_damaged\":[";
            write_dnan_arr(lb.per_pos_5prime_ct_damaged);
            j << "],\"per_pos_3prime_damaged\":[";
            write_dnan_arr(lb.per_pos_3prime_damaged);
            j << "],\"per_pos_5prime_ct_cpg_damaged\":[";
            write_dnan_arr(lb.per_pos_5prime_ct_cpg_damaged);
            j << "],\"per_pos_5prime_ct_noncpg_damaged\":[";
            write_dnan_arr(lb.per_pos_5prime_ct_noncpg_damaged);
            j << "],\"per_pos_5prime_gt_damaged\":[";
            write_dnan_arr(lb.per_pos_5prime_gt_damaged);
            j << "],\"per_pos_5prime_gt_undamaged\":[";
            write_dnan_arr(lb.per_pos_5prime_gt_undamaged);
            j << "]}";
            j << ",\"trinuc\":{\"tri_5prime_terminal\":[";
            for (int i = 0; i < 64; ++i) { j << lb.tri_5prime_terminal[i]; if (i < 63) j << ","; }
            j << "],\"tri_5prime_interior\":[";
            for (int i = 0; i < 64; ++i) { j << lb.tri_5prime_interior[i]; if (i < 63) j << ","; }
            j << "],\"tri_3prime_terminal\":[";
            for (int i = 0; i < 64; ++i) { j << lb.tri_3prime_terminal[i]; if (i < 63) j << ","; }
            j << "],\"tri_3prime_interior\":[";
            for (int i = 0; i < 64; ++i) { j << lb.tri_3prime_interior[i]; if (i < 63) j << ","; }
            j << "]}";
            j << "}";
        }
        if (!lsd.bins.empty()) j << "\n    ";
        j << "],\n";
        j << "    \"by_length_joint\": {\n";
        j << "      \"d_ancient\": "   << std::setprecision(6) << lsd.d_joint_ancient << ",\n";
        j << "      \"pi_ancient\": "  << lsd.pi_joint_ancient << ",\n";
        j << "      \"d_population\": " << lsd.d_joint_population << ",\n";
        j << "      \"converged\": "   << (lsd.joint_converged ? "true" : "false") << ",\n";
        j << "      \"separated\": "   << (lsd.joint_separated ? "true" : "false") << ",\n";
        j << "      \"cell_w_ancient\": [";
        for (size_t b = 0; b < lsd.cell_w_ancient.size(); ++b) {
            if (b > 0) j << ",";
            j << "[";
            const auto& row = lsd.cell_w_ancient[b];
            for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
                // C1: -1.0 is the insufficient-coverage sentinel; w_ancient in [0,1] → null.
                if (row[g] < 0.0) j << "null";
                else j << std::setprecision(6) << row[g];
                if (g + 1 < LengthBinDamageProfile::N_GC_BINS) j << ",";
            }
            j << "]";
        }
        j << "]\n    }\n";
    }

    j << "  },\n";  // end deamination

    // ── Complement asymmetry ──────────────────────────────────────────────────
    j << "  \"complement_asymmetry\": {\n";
    j << "    \"D\": " << std::setprecision(6) << dp.ox_gt_asymmetry << ",\n";
    j << "    \"tg_interior\": " << dp.ox_gt_baseline << ",\n";
    j << "    \"ac_interior\": " << dp.ox_ca_baseline << ",\n";
    j << "    \"tg_terminal\": " << dp.ox_gt_rate_terminal << ",\n";
    j << "    \"ac_terminal\": " << dp.ox_ca_rate_terminal << ",\n";
    j << "    \"gt_bg_fitted\": " << dp.g_bg_fitted << ",\n";
    j << "    \"gt_term_fitted\": " << dp.g_term_fitted << ",\n";
    j << "    \"gt_decay_fitted\": " << dp.g_decay_fitted << ",\n";
    j << "    \"gt_bg_ci_lo\": " << std::setprecision(6) << dp.g_bg_fitted_ci_lo << ",\n";
    j << "    \"gt_bg_ci_hi\": " << dp.g_bg_fitted_ci_hi << ",\n";
    j << "    \"gt_fit_degenerate\": " << (dp.g_fit_degenerate ? "true" : "false") << ",\n";
    j << "    \"gt_bg_fitted_unclamped\": " << dp.g_bg_fitted_unclamped << ",\n";
    j << "    \"gt_bg_at_upper_boundary\": " << (dp.gt_bg_at_upper_boundary ? "true" : "false") << ",\n";
    j << "    \"gt_decay_at_upper_boundary\": " << (dp.gt_decay_at_upper_boundary ? "true" : "false") << ",\n";
    j << "    \"gt_term_zero_clamped\": " << (dp.gt_term_zero_clamped ? "true" : "false") << ",\n";
    j << "    \"ox_theta_at_clamp\": " << (dp.ox_theta_at_clamp ? "true" : "false") << ",\n";
    j << "    \"gt_bg_interior_mean\": " << dp.g_bg_interior_mean << ",\n";
    j << "    \"s_gt\": " << dp.s_gt << ",\n";
    // D/s_gt are applicable only to SS (Chargaff cancellation makes them ≈0 for DS).
    j << "    \"d_applicable\": " << (is_ss ? "true" : "false") << ",\n";
    j << "    \"per_pos_5prime_gt\": [";
    for (int p = 0; p < N_POS; ++p) {
        double denom = dp.t_from_g_5prime[p] + dp.g_count_5prime[p];
        // C1: emit null (not -1.0) for sub-kMinCov positions — out-of-range for a probability.
        if (denom >= kMinCov) j << std::setprecision(6) << (dp.t_from_g_5prime[p] / denom);
        else                  j << "null";
        if (p < N_POS - 1) j << ",";
    }
    j << "],\n";
    // C1+C2: FGH binom z-scores. Emit null when the channel was not computed (validity flag false:
    // the 0.0 default is indistinguishable from a genuine z=0), otherwise clamp to ±kZCap — these
    // z's scale ~sqrt(N) on correlated reads (magnitudes of hundreds–thousands) and are exploratory,
    // not calibrated p-values; the rate-difference / common_or fields carry the effect size.
    auto emit_z = [&](double z, bool computed) {
        if (computed && std::isfinite(z)) j << std::setprecision(6) << clamp_z(z);
        else                              j << "null";
    };
    j << "    \"channel_c_valid\": " << (dp.channel_c_valid ? "true" : "false") << ",\n";
    j << "    \"channel_c_detected\": " << (dp.ox_damage_detected ? "true" : "false") << ",\n";
    j << "    \"ox_is_artifact\": " << (dp.ox_is_artifact ? "true" : "false") << ",\n";
    j << "    \"ox_d_max\": " << std::setprecision(6) << dp.ox_d_max << ",\n";
    j << "    \"ox_stop_rate_terminal\": " << dp.ox_stop_rate_terminal << ",\n";
    j << "    \"ox_stop_rate_interior\": " << dp.ox_stop_rate_interior << ",\n";
    j << "    \"ox_stop_rate_baseline\": " << dp.ox_stop_conversion_rate_baseline << ",\n";
    j << "    \"ox_uniformity_ratio\": " << dp.ox_uniformity_ratio << ",\n";
    j << "    \"channel_c3_valid\": " << (dp.channel_c3_valid ? "true" : "false") << ",\n";
    j << "    \"ox_stop_baseline_3prime\": " << std::setprecision(6) << dp.ox_stop_baseline_3prime << ",\n";
    j << "    \"ox_stop_rate_terminal_3prime\": " << dp.ox_stop_rate_terminal_3prime << ",\n";
    j << "    \"ox_stop_rate_interior_3prime\": " << dp.ox_stop_rate_interior_3prime << ",\n";
    j << "    \"ox_uniformity_ratio_3prime\": " << dp.ox_uniformity_ratio_3prime << ",\n";
    j << "    \"channel_f_valid\": " << (dp.channel_f_valid ? "true" : "false") << ",\n";
    j << "    \"ca_stop_rate_baseline\": " << std::setprecision(6) << dp.ca_stop_rate_baseline << ",\n";
    j << "    \"ca_stop_rate_terminal\": " << dp.ca_stop_rate_terminal << ",\n";
    j << "    \"ca_uniformity_ratio\": " << dp.ca_uniformity_ratio << ",\n";
    j << "    \"channel_f_z\": ";    emit_z(dp.channel_f_z, dp.channel_f_valid);    j << ",\n";
    j << "    \"channel_f_mh_z\": "; emit_z(dp.channel_f_mh_z, dp.channel_f_valid); j << ",\n";
    // Transparency: the pooled channel_f_z pools across context strata and keeps the deamination
    // shadow in its denominator, so it can sign-reverse vs the stratified MH z + odds ratio
    // (Simpson's paradox). Detection sign-gates on agreement; this flag warns consumers not to read
    // the raw pooled z as the effect direction when it disagrees with the MH-stratified result.
    j << "    \"channel_f_pooled_mh_disagree\": "
      << ((std::isfinite(dp.channel_f_z) && std::isfinite(dp.channel_f_mh_z)
           && ((dp.channel_f_z >= 0.0f) != (dp.channel_f_mh_z >= 0.0f))) ? "true" : "false") << ",\n";
    j << "    \"channel_f_common_or\": " << dp.channel_f_common_or << ",\n";
    j << "    \"channel_f_applicable\": " << (!is_ss ? "true" : "false") << ",\n";
    j << "    \"ca_stop_rate_interior\": " << dp.ca_stop_rate_interior << ",\n";
    j << "    \"channel_f3_valid\": " << (dp.channel_f3_valid ? "true" : "false") << ",\n";
    j << "    \"ca_stop_rate_terminal_3prime\": " << dp.ca_stop_rate_terminal_3prime << ",\n";
    j << "    \"ca_stop_rate_baseline_3prime\": " << dp.ca_stop_rate_baseline_3prime << ",\n";
    j << "    \"ca_stop_rate_interior_3prime\": " << dp.ca_stop_rate_interior_3prime << ",\n";
    j << "    \"ca_uniformity_ratio_3prime\": " << dp.ca_uniformity_ratio_3prime << ",\n";
    j << "    \"channel_g_valid\": " << (dp.channel_g_valid ? "true" : "false") << ",\n";
    j << "    \"cg_stop_rate_baseline\": " << std::setprecision(6) << dp.cg_stop_rate_baseline << ",\n";
    j << "    \"cg_stop_rate_terminal\": " << dp.cg_stop_rate_terminal << ",\n";
    j << "    \"cg_uniformity_ratio\": " << dp.cg_uniformity_ratio << ",\n";
    j << "    \"channel_g_z\": "; emit_z(dp.channel_g_z, dp.channel_g_valid); j << ",\n";
    j << "    \"channel_g_or\": " << nan_or(dp.channel_g_or) << ",\n";  // P4: 2x2 Haldane-Anscombe OR
    j << "    \"channel_g_applicable\": " << (!is_ss ? "true" : "false") << ",\n";
    j << "    \"cg_stop_rate_interior\": " << dp.cg_stop_rate_interior << ",\n";
    j << "    \"channel_g3_valid\": " << (dp.channel_g3_valid ? "true" : "false") << ",\n";
    j << "    \"cg_stop_rate_terminal_3prime\": " << dp.cg_stop_rate_terminal_3prime << ",\n";
    j << "    \"cg_stop_rate_baseline_3prime\": " << dp.cg_stop_rate_baseline_3prime << ",\n";
    j << "    \"cg_stop_rate_interior_3prime\": " << dp.cg_stop_rate_interior_3prime << ",\n";
    j << "    \"cg_uniformity_ratio_3prime\": " << dp.cg_uniformity_ratio_3prime << ",\n";
    j << "    \"channel_h_valid\": " << (dp.channel_h_valid ? "true" : "false") << ",\n";
    j << "    \"at_stop_rate_baseline\": " << std::setprecision(6) << dp.at_stop_rate_baseline << ",\n";
    j << "    \"at_stop_rate_terminal\": " << dp.at_stop_rate_terminal << ",\n";
    j << "    \"at_uniformity_ratio\": " << dp.at_uniformity_ratio << ",\n";
    j << "    \"channel_h_z\": ";        emit_z(dp.channel_h_z, dp.channel_h_valid);        j << ",\n";
    j << "    \"channel_h_or\": " << nan_or(dp.channel_h_or) << ",\n";  // P4: 2x2 Haldane-Anscombe OR
    j << "    \"channel_h_z_p2plus\": "; emit_z(dp.channel_h_z_p2plus, dp.channel_h_valid); j << ",\n";
    j << "    \"at_stop_rate_interior\": " << dp.at_stop_rate_interior << ",\n";
    j << "    \"channel_h3_valid\": " << (dp.channel_h3_valid ? "true" : "false") << ",\n";
    j << "    \"at_stop_rate_terminal_3prime\": " << dp.at_stop_rate_terminal_3prime << ",\n";
    j << "    \"at_stop_rate_baseline_3prime\": " << dp.at_stop_rate_baseline_3prime << ",\n";
    j << "    \"at_stop_rate_interior_3prime\": " << dp.at_stop_rate_interior_3prime << ",\n";
    j << "    \"at_uniformity_ratio_3prime\": " << dp.at_uniformity_ratio_3prime << ",\n";
    j << "    \"ox_gt_rate_interior\": " << dp.ox_gt_rate_interior << ",\n";
    j << "    \"ox_gt_uniformity\": " << dp.ox_gt_uniformity << ",\n";
    j << "    \"ox_ca_rate_interior\": " << dp.ox_ca_rate_interior << ",\n";
    j << "    \"ox_ca_uniformity\": " << dp.ox_ca_uniformity << ",\n";
    // C1: gate the oxoG-score fields on has_oxog_score — emit null when uncomputed so a
    // genuine 0.0 (no 8-oxoG signal) is distinguishable from "not computed". se_s_oxog=0.0
    // as a fallback was especially misleading (signals zero-width CI). Surface the flag too.
    auto emit_oxog = [&](double v) {
        if (in.has_oxog_score && std::isfinite(v)) j << std::setprecision(6) << v; else j << "null";
    };
    j << "    \"has_oxog_score\": " << (in.has_oxog_score ? "true" : "false") << ",\n";
    j << "    \"s_oxog\": ";        emit_oxog(in.s_oxog);    j << ",\n";
    j << "    \"se_s_oxog\": ";     emit_oxog(in.se_s_oxog); j << ",\n";
    j << "    \"ox_d_oriented\": "; emit_oxog(in.d_oriented); j << ",\n";
    j << "    \"s_oxog_16ctx\": [";
    for (int i = 0; i < 16; ++i) {
        float v = dp.s_oxog_16ctx[i];
        if (std::isnan(v)) j << "null"; else j << std::setprecision(6) << v;
        if (i < 15) j << ",";
    }
    j << "],\n";
    j << "    \"cov_oxog_16ctx\": [";
    for (int i = 0; i < 16; ++i) {
        j << std::setprecision(0) << dp.cov_oxog_16ctx[i];
        if (i < 15) j << ",";
    }
    j << "],\n";
    {
        // C2 weight-fix: s_i is a Binomial rate-difference with Var(s_i)≈1/cov_i, so the correct
        // inverse-variance weight is w_i=cov_i and the IVW z = Σ(s_i·cov_i)/√(Σcov_i) (was √cov_i,
        // which is not a z-score under any standard definition). The IVW z then scales ~√N on
        // correlated reads, so it is also clamped to ±kZCap. C1: when no context is evaluable
        // (all cov<500 ⇒ ctx_w2==0) emit null, not a bare 0.0 that reads as a real zero aggregate.
        double ctx_wsum = 0.0, ctx_w2 = 0.0;
        int ctx_pos = 0, ctx_n_valid = 0;
        for (int i = 0; i < 16; ++i) {
            float v = dp.s_oxog_16ctx[i];
            double cov = static_cast<double>(dp.cov_oxog_16ctx[i]);
            if (std::isnan(v) || cov <= 0) continue;
            ctx_wsum += v * cov;   // w_i = cov_i (inverse-variance)
            ctx_w2   += cov;
            ++ctx_n_valid;
            if (v > 0) ++ctx_pos;
        }
        j << "    \"oxog_ctx_n_positive\": " << ctx_pos << ",\n";
        j << "    \"oxog_ctx_n_valid\": " << ctx_n_valid << ",\n";
        // exploratory; clamped, not a calibrated p-value (correlated reads)
        if (ctx_w2 > 0) j << "    \"oxog_ctx_z\": " << std::setprecision(4)
                          << clamp_z(ctx_wsum / std::sqrt(ctx_w2)) << ",\n";
        else            j << "    \"oxog_ctx_z\": null,\n";
        // C1/C3: gate on has_oxog_score (not computed) and is_ss (RC parity assumption fails
        // for single-strand reads, producing systematic false-positive scores).
        j << "    \"oxog_score_rc_parity_ok\": " << (!is_ss ? "true" : "false") << ",\n";
        j << "    \"oxog_score_z\": "; emit_oxog(!is_ss ? clamp_z(oxog_is.z) : std::numeric_limits<double>::quiet_NaN()); j << ",\n";
        j << "    \"oxog_score_p\": "; emit_oxog(!is_ss ? oxog_is.p : std::numeric_limits<double>::quiet_NaN()); j << ",\n";
    }
    if (std::isnan(otr.cosine))
        j << "    \"oxog_context_cosine\": null,\n";
    else
        j << "    \"oxog_context_cosine\": " << std::setprecision(6) << otr.cosine << ",\n";
    j << "    \"oxog_trinuc_n_context\": " << otr.n_ctx << ",\n";
    if (std::isnan(otr.gt_asymmetry))
        j << "    \"oxog_gt_asymmetry\": null,\n    \"oxog_gt_rate\": null,\n";
    else
        j << "    \"oxog_gt_asymmetry\": " << std::setprecision(6) << otr.gt_asymmetry << ",\n"
          << "    \"oxog_gt_rate\": " << std::setprecision(6) << otr.gt_rate << ",\n";

    // complement_asymmetry by_length (oxog channels per bin)
    if (in.lsd && !in.lsd->bins.empty()) {
        const auto& lsd = *in.lsd;
        j << "    \"by_length\": [";
        for (size_t b = 0; b < lsd.bins.size(); ++b) {
            const auto& lb = lsd.bins[b];
            if (b > 0) j << ",";
            j << std::setprecision(6);
            j << "{\"length_lo\":" << lb.length_lo << ",\"length_hi\":" << lb.length_hi
              << ",\"n_reads\":" << lb.n_reads
              << ",\"channel_c\":{\"baseline_rate\":" << lb.ox_stop_rate_baseline
              << ",\"terminal_rate\":" << lb.ox_stop_rate_terminal
              << ",\"uniformity_ratio\":" << lb.ox_uniformity_ratio
              << ",\"valid\":" << (lb.channel_c_valid ? "true" : "false") << "}"
              << ",\"channel_f\":{\"baseline_rate\":" << lb.ca_stop_rate_baseline
              << ",\"terminal_rate\":" << lb.ca_stop_rate_terminal
              << ",\"uniformity_ratio\":" << lb.ca_uniformity_ratio
              << ",\"z_score\":" << lb.channel_f_z
              << ",\"valid\":" << (lb.channel_f_valid ? "true" : "false")
              << ",\"detected\":" << (lb.channel_f_valid && lb.channel_f_z > kOxChannelZDetect ? "true" : "false") << "}"
              << ",\"channel_g\":{\"baseline_rate\":" << lb.cg_stop_rate_baseline
              << ",\"terminal_rate\":" << lb.cg_stop_rate_terminal
              << ",\"uniformity_ratio\":" << lb.cg_uniformity_ratio
              << ",\"z_score\":" << lb.channel_g_z
              << ",\"valid\":" << (lb.channel_g_valid ? "true" : "false")
              << ",\"detected\":" << (lb.channel_g_valid && lb.channel_g_z > kOxChannelZDetect ? "true" : "false") << "}"
              << ",\"channel_h\":{\"baseline_rate\":" << lb.at_stop_rate_baseline
              << ",\"terminal_rate\":" << lb.at_stop_rate_terminal
              << ",\"uniformity_ratio\":" << lb.at_uniformity_ratio
              << ",\"z_score\":" << lb.channel_h_z
              << ",\"z_score_p2plus\":" << lb.channel_h_z_p2plus
              << ",\"valid\":" << (lb.channel_h_valid ? "true" : "false")
              << ",\"detected\":" << (lb.channel_h_valid && lb.channel_h_z > kOxChannelZDetect && lb.channel_h_z_p2plus > kOxChannelZDetect ? "true" : "false") << "}"  // C5: OR→AND, see top-level gate [behavioral change]
              << "}";
        }
        j << "],\n";
    }
    j << "    \"fgh_adapter_prefixes_excluded\": " << dp.fgh_adapter_prefixes_excluded << "\n";
    j << "  },\n";

    // ── OxoG unified estimate ─────────────────────────────────────────────────
    j << "  \"oxog_estimate\": {\n";
    j << "    \"oxo_schema\": " << OxoGEstimate::oxo_schema << ",\n";
    j << "    \"ox_theta\": " << std::setprecision(6) << oxe.ox_theta << ",\n";
    j << "    \"ox_theta_ci_lo\": " << oxe.ox_theta_ci_lo << ",\n";
    j << "    \"ox_theta_ci_hi\": " << oxe.ox_theta_ci_hi << ",\n";
    j << "    \"ox_theta_interior\": " << oxe.ox_theta_interior << ",\n";
    j << "    \"fit_degenerate\": " << (oxe.fit_degenerate ? "true" : "false") << ",\n";
    j << "    \"ox_like_excess\": " << oxe.ox_like_excess << ",\n";
    j << "    \"ox_like_z\": " << oxe.ox_like_z << ",\n";
    j << "    \"ox_like_ci_lo\": " << oxe.ox_like_ci_lo << ",\n";
    j << "    \"ox_like_ci_hi\": " << oxe.ox_like_ci_hi << ",\n";
    j << "    \"ox_uniformity_ratio\": " << oxe.ox_uniformity_ratio << ",\n";
    j << "    \"control_mode\": \"" << oxe.control_mode << "\",\n";
    j << "    \"gc_skew_warning\": " << (oxe.gc_skew_warning ? "true" : "false") << "\n";
    j << "  },\n";

    // ── OxoG two-marker ───────────────────────────────────────────────────────
    j << "  \"oxo_two_marker\": {\n";
    // Wave-1 item 6: oxo_two_marker is the PRIMARY reference-free oxidation readout (interior
    // Chargaff D regressed on ssDNA-overhang markers). role/interpretation are non-numeric labels.
    j << "    \"role\": \"primary_oxidation_readout\",\n";
    j << "    \"interpretation\": \"Interior-Chargaff D regressed on ssDNA-overhang markers; recovers modest interior oxidation (G to T strong, C to A weak) that the terminal-scope channels C/D/F miss. C/D/F are scope-limited detectors whose nulls are informative, not the primary oxidation verdict.\",\n";
    j << "    \"valid\": " << (otm.valid ? "true" : "false") << ",\n";
    j << "    \"n_cells_used\": " << otm.n_cells_used << ",\n";
    j << "    \"beta1\": " << std::setprecision(6) << otm.beta1 << ",\n";
    j << "    \"beta1_se\": " << otm.beta1_se << ",\n";
    // exploratory; clamped, not a calibrated p-value (correlated reads)
    j << "    \"beta1_z\": " << clamp_z(otm.beta1_z) << ",\n";
    j << "    \"beta2\": " << otm.beta2 << ",\n";
    j << "    \"beta2_se\": " << otm.beta2_se << ",\n";
    // exploratory; clamped, not a calibrated p-value (correlated reads)
    j << "    \"beta2_z\": " << clamp_z(otm.beta2_z) << ",\n";
    j << "    \"alpha\": " << otm.alpha << ",\n";
    j << "    \"sigma2\": " << otm.sigma2 << ",\n";   // residual variance — makes the SEs auditable
    j << "    \"delta_beta\": " << otm.delta_beta << ",\n";
    j << "    \"markers_consistent\": " << (otm.markers_consistent ? "true" : "false") << "\n";
    j << "  },\n";

    // ── Interior CT cluster ───────────────────────────────────────────────────
    j << "  \"interior_ct_cluster\": {\n";
    j << "    \"short_asym_log2oe\": " << std::setprecision(6) << dp.interior_ct_cluster_short_asym_log2oe << ",\n";
    j << "    \"short_log2oe\": " << dp.interior_ct_cluster_short_log2oe << ",\n";
    j << "    \"short_z\": " << dp.interior_ct_cluster_short_z << ",\n";
    j << "    \"short_obs\": " << dp.interior_ct_cluster_short_obs << ",\n";
    j << "    \"short_exp\": " << dp.interior_ct_cluster_short_exp << ",\n";
    j << "    \"reads_used\": " << dp.interior_ct_cluster_reads_used << ",\n";
    j << "    \"sep_log2oe\": [";
    for (int i = 0; i < 10; ++i) {
        j << std::setprecision(6) << dp.interior_ct_cluster_sep_log2oe[i];
        if (i < 9) j << ",";
    }
    j << "],\n";
    {
        // C2: signed-root G-test on aggregate (correlated) pair counts O/E — scales ~sqrt(pairs),
        // not a calibrated p-value. Clamp the emitted z to ±kZCap and floor p at DBL_MIN so it
        // never underflows to a literal 0.0 (which happened for every >50M-read library, damaged
        // or modern alike). short_log2oe (emitted above) is the primary interpretable effect size.
        double O = dp.interior_ct_cluster_short_obs;
        double E = dp.interior_ct_cluster_short_exp;
        double sz = 0.0, sp = 1.0;
        if (O > 0 && E > 1e-9) {
            double sign_val = (O >= E) ? 1.0 : -1.0;
            sz = sign_val * std::sqrt(2.0 * (O * std::log(O / E) - (O - E)));
            sp = 0.5 * std::erfc(sz / std::sqrt(2.0));
            // floor at DBL_MIN: |sz|>~38 underflows erfc to exactly 0.0 for any large library
            if (sp < std::numeric_limits<double>::min()) sp = std::numeric_limits<double>::min();
        }
        // exploratory; clamped, not a calibrated p-value (correlated reads)
        // C5 label fix: this is the signed-root G-test (likelihood-ratio statistic), NOT the
        // within-read Bernoulli z (short_z). Same value was published under two `*_z` names with
        // wildly different magnitudes (e.g. short_z=1.34 vs short_score_z=248.7). Emit honest
        // short_g_stat/short_lr_p; keep short_score_z/short_score_p as deprecated aliases.
        const double sz_emit = clamp_z(sz);
        j << "    \"short_g_stat\": " << std::setprecision(6) << sz_emit << ",\n";
        j << "    \"short_lr_p\": " << sp << ",\n";
        j << "    \"short_score_z\": " << std::setprecision(6) << sz_emit << ",\n";
        j << "    \"short_score_p\": " << sp << "\n";
    }
    j << "  },\n";

    // ── Depurination ──────────────────────────────────────────────────────────
    j << "  \"depurination\": {\n";
    j << "    \"detected\": " << (dp.depurination_detected ? "true" : "false") << ",\n";
    j << "    \"enrichment_5prime\": " << std::setprecision(6) << dp.purine_enrichment_5prime << ",\n";
    j << "    \"enrichment_3prime\": " << dp.purine_enrichment_3prime << ",\n";
    j << "    \"rate_interior\": " << std::setprecision(6) << dp.purine_rate_interior << ",\n";
    j << "    \"depur_ctrl_shift_5prime\": " << ds.shift5 << ",\n";
    j << "    \"depur_ctrl_shift_3prime\": " << ds.shift3 << ",\n";
    j << "    \"depur_score_z_5prime\": " << ds.z5 << ",\n";
    j << "    \"depur_score_z_3prime\": "
      << (std::isnan(ds.z3) ? "null" : std::to_string(ds.z3)) << ",\n";
    j << "    \"depur_score_z\": " << ds.z << ",\n";
    j << "    \"depur_score_p\": " << ds.p << "\n";
    j << "  },\n";

    // ── Trinucleotide spectrum ─────────────────────────────────────────────────
    {
        auto emit_arr = [&](const char* name, const std::array<uint64_t, 64>& v, bool trailing) {
            j << "    \"" << name << "\": [";
            for (int i = 0; i < 64; ++i) { j << v[i]; if (i < 63) j << ","; }
            j << "]" << (trailing ? "," : "") << "\n";
        };
        j << "  \"trinuc_spectrum\": {\n";
        emit_arr("tri_5prime_terminal", dp.tri_5prime_terminal, true);
        emit_arr("tri_5prime_interior", dp.tri_5prime_interior, true);
        emit_arr("tri_3prime_terminal", dp.tri_3prime_terminal, true);
        emit_arr("tri_3prime_interior", dp.tri_3prime_interior, false);
        j << "  },\n";
    }

    // ── Per-position trinucleotide counts ─────────────────────────────────────
    // [pos][64 contexts], pos 1..N_POS_TRI-1 (pos 0 skipped — no left flank).
    // Downstream: normalise T/(T+C) per XCY context per position for reference-free
    // positional damage spectra analogous to bam2sbs sbs3d output.
    {
        j << "  \"trinuc_pos_spectrum\": {\n";
        for (const char* end : {"5prime", "3prime"}) {
            const auto& arr = (end[0] == '5') ? dp.tri_5prime_pos : dp.tri_3prime_pos;
            j << "    \"" << end << "\": [";
            for (int p = 0; p < SampleDamageProfile::N_POS_TRI; ++p) {
                j << "[";
                for (int i = 0; i < SampleDamageProfile::N_TRINUC; ++i) {
                    j << arr[p][i];
                    if (i < SampleDamageProfile::N_TRINUC - 1) j << ",";
                }
                j << "]";
                if (p < SampleDamageProfile::N_POS_TRI - 1) j << ",";
            }
            j << "]";
            if (end[0] == '5') j << ",\n"; else j << "\n";
        }
        j << "  },\n";
    }

    // ── Per-position substitution rates (all 12 types) ────────────────────────
    // rate(X→Y, pos) = alt/(ref+alt) collapsed over flanking context.
    // Allows detection of damage types beyond C→T/G→A: AP-site A-rule (G→T
    // elevated at internal positions), CPD-like dipyrimidine patterns, etc.
    // pos 0 = read terminus; rate = -1 when ref+alt counts are zero.
    {
        static constexpr struct { int from; int to; const char* name; } SUBS[12] = {
            {1,3,"CT"},{1,0,"CA"},{1,2,"CG"},
            {2,0,"GA"},{2,3,"GT"},{2,1,"GC"},
            {3,1,"TC"},{3,0,"TA"},{3,2,"TG"},
            {0,2,"AG"},{0,3,"AT"},{0,1,"AC"},
        };

        auto emit_end = [&](const char* end_key,
                            const std::array<std::array<uint64_t,64>,15>& pos_arr) {
            j << "    \"" << end_key << "\": {\n";
            for (int s = 0; s < 12; ++s) {
                int F = SUBS[s].from, T = SUBS[s].to;
                j << "      \"" << SUBS[s].name << "\": [";
                for (int p = 0; p < 15; ++p) {
                    uint64_t nf = 0, nt = 0;
                    for (int prev = 0; prev < 4; ++prev)
                        for (int next = 0; next < 4; ++next) {
                            nf += pos_arr[p][prev*16 + F*4 + next];
                            nt += pos_arr[p][prev*16 + T*4 + next];
                        }
                    // C1: pos 0 has no left-flank context (nf+nt==0) — emit null, not -1.0,
                    // so the structurally-missing terminus element is unambiguous downstream.
                    if (nf + nt > 0) j << std::fixed << std::setprecision(6) << ((double)nt / (nf + nt));
                    else             j << "null";
                    if (p < 14) j << ",";
                }
                j << "]" << (s < 11 ? "," : "") << "\n";
            }
            j << "    }";
        };

        j << "  \"substitution_pos_profiles\": {\n";
        j << "    \"n_positions\": " << SampleDamageProfile::N_POS_TRI << ",\n";
        emit_end("5prime", dp.tri_5prime_pos);
        j << ",\n";
        emit_end("3prime", dp.tri_3prime_pos);
        j << "\n  },\n";
    }

    // ── Per-channel damage statistics ─────────────────────────────────────────
    // Three improvements over fixed-window terminal-vs-interior:
    //
    // 1. Adaptive interior baseline: plateau onset detected from C→T profile
    //    per library, not assumed at fixed positions 10-14. Fast-decay libraries
    //    plateau at pos 3; slow-decay or low-damage libraries may not plateau by
    //    14. Uses tail positions 9-14 as anchor; finds first p≥3 where 3
    //    consecutive positions are within 3σ of the tail mean.
    //
    // 2. Per-channel decay lambda: weighted log-linear fit of excess vs position
    //    for each substitution type independently. lambda_ratio_vs_ct < 1 means
    //    flatter than deamination (AP-site A-rule, oxidative internal damage);
    //    > 1 means steeper (very terminal-concentrated). Divergent shape is the
    //    primary evidence for a damage type with different geometry.
    //
    // 3. Coupled donor/acceptor balance: for X→Y signal, checks that FROM is
    //    depleted AND TO is enriched at terminal vs interior in proportion. A
    //    real substitution-like signal has balance ≈ 1.0; endpoint composition
    //    selection produces imbalance (more enrichment than depletion, or vice
    //    versa). balance = -1 when either side is non-positive.
    {
        static constexpr struct { int from; int to; const char* name; } SUBS12[12] = {
            {1,3,"CT"},{1,0,"CA"},{1,2,"CG"},
            {2,0,"GA"},{2,3,"GT"},{2,1,"GC"},
            {3,1,"TC"},{3,0,"TA"},{3,2,"TG"},
            {0,2,"AG"},{0,3,"AT"},{0,1,"AC"},
        };

        // Marginal per-position rate for FROM→TO collapsed over flanking context
        auto pos_rates = [](const std::array<std::array<uint64_t,64>,15>& arr,
                            int F, int T, std::array<double,15>& out) {
            for (int p = 0; p < 15; ++p) {
                uint64_t nf = 0, nt = 0;
                for (int pv = 0; pv < 4; ++pv)
                    for (int nx = 0; nx < 4; ++nx) {
                        nf += arr[p][pv*16 + F*4 + nx];
                        nt += arr[p][pv*16 + T*4 + nx];
                    }
                out[p] = (nf + nt > 0) ? (double)nt / (nf + nt) : -1.0;
            }
        };

        // Adaptive plateau start: first p≥3 where 3 consecutive positions
        // are within max(0.005, 3σ) of the positions-9-14 tail mean.
        auto find_plateau = [](const std::array<double,15>& rates) -> int {
            double tail = 0.0, tail2 = 0.0; int n = 0;
            for (int p = 9; p < 15; ++p)
                if (rates[p] >= 0) { tail += rates[p]; tail2 += rates[p]*rates[p]; n++; }
            if (n < 2) return 5;
            tail /= n;
            double var = std::max(0.0, tail2/n - tail*tail);
            double tol = std::max(0.005, 3.0 * std::sqrt(var));
            for (int p = 3; p < 12; ++p) {
                if (rates[p] < 0) continue;
                bool ok = true;
                for (int q = p; q < std::min(p+3, 15); ++q)
                    if (rates[q] >= 0 && std::abs(rates[q] - tail) > tol) { ok=false; break; }
                if (ok) return p;
            }
            return 5;
        };

        auto interior_mean = [](const std::array<double,15>& r, int plat) -> double {
            double s = 0.0; int n = 0;
            for (int p = plat; p < 15; ++p) if (r[p] >= 0) { s += r[p]; n++; }
            return n > 0 ? s/n : std::numeric_limits<double>::quiet_NaN();
        };

        // Weighted log-linear OLS: log(excess) = log(A) - lambda*p, w = excess.
        auto fit_lambda = [](const std::array<double,15>& rates, int plat, double bg) -> double {
            double sw=0,swx=0,swy=0,swxx=0,swxy=0;
            for (int p = 1; p < plat; ++p) {
                if (rates[p] < 0) continue;
                double ex = rates[p] - bg;
                if (ex < 1e-4) continue;
                double lx = std::log(ex);
                sw+=ex; swx+=ex*p; swy+=ex*lx; swxx+=ex*p*p; swxy+=ex*p*lx;
            }
            if (sw < 1e-12) return std::numeric_limits<double>::quiet_NaN();
            double det = sw*swxx - swx*swx;
            if (std::abs(det) < 1e-18) return std::numeric_limits<double>::quiet_NaN();
            return -(sw*swxy - swx*swy) / det;
        };

        // Coupled balance: donor_depletion(FROM) / acceptor_enrichment(TO)
        // over terminal zone [0, min(3,plat-1)] vs interior zone [plat,15).
        auto coupled_bal = [](const std::array<std::array<uint64_t,64>,15>& arr,
                              int F, int T, int plat) -> double {
            uint64_t tc[4]{}, ic[4]{};
            for (int p = 0; p <= std::min(3, plat-1); ++p)
                for (int b = 0; b < 4; ++b)
                    for (int pv = 0; pv < 4; ++pv)
                        for (int nx = 0; nx < 4; ++nx)
                            tc[b] += arr[p][pv*16 + b*4 + nx];
            for (int p = plat; p < 15; ++p)
                for (int b = 0; b < 4; ++b)
                    for (int pv = 0; pv < 4; ++pv)
                        for (int nx = 0; nx < 4; ++nx)
                            ic[b] += arr[p][pv*16 + b*4 + nx];
            uint64_t tt=tc[0]+tc[1]+tc[2]+tc[3], it=ic[0]+ic[1]+ic[2]+ic[3];
            if (!tt || !it) return std::numeric_limits<double>::quiet_NaN();
            double depl   = (double)ic[F]/it - (double)tc[F]/tt;
            double enrich = (double)tc[T]/tt  - (double)ic[T]/it;
            if (depl <= 0.0 || enrich <= 0.0) return std::numeric_limits<double>::quiet_NaN();
            return depl / enrich;
        };

        auto emit_dcs_end = [&](const char* end_key,
                                const std::array<std::array<uint64_t,64>,15>& arr) {
            std::array<double,15> ct_r;
            pos_rates(arr, 1, 3, ct_r);
            int   plat   = find_plateau(ct_r);
            double ct_bg  = interior_mean(ct_r, plat);
            double ct_lam = std::isfinite(ct_bg) ? fit_lambda(ct_r, plat, ct_bg) : std::numeric_limits<double>::quiet_NaN();

            j << "    \"" << end_key << "\": {\n";
            j << "      \"interior_start_pos\": " << plat << ",\n";
            j << "      \"ct_interior_rate\": "   << std::fixed << std::setprecision(6) << ct_bg  << ",\n";
            j << "      \"ct_decay_lambda\": "    << nan_or(ct_lam) << ",\n";
            j << "      \"channels\": {";
            bool first = true;
            for (int s = 0; s < 12; ++s) {
                int F = SUBS12[s].from, T = SUBS12[s].to;
                std::array<double,15> r;
                pos_rates(arr, F, T, r);
                double bg  = interior_mean(r, plat);
                double tex = (r[1] >= 0 && std::isfinite(bg)) ? r[1] - bg : std::numeric_limits<double>::quiet_NaN();
                double lam = std::isfinite(bg) ? fit_lambda(r, plat, bg) : std::numeric_limits<double>::quiet_NaN();
                double lr  = (std::isfinite(lam) && lam > 0 && std::isfinite(ct_lam) && ct_lam > 0) ? lam / ct_lam : std::numeric_limits<double>::quiet_NaN();
                double bal = coupled_bal(arr, F, T, plat);
                if (!first) j << ",";
                first = false;
                j << "\n        \"" << SUBS12[s].name << "\":"
                  << "{\"interior_rate\":"      << std::fixed << std::setprecision(6) << nan_or(bg)
                  << ",\"terminal_excess\":"    << nan_or(tex)
                  << ",\"decay_lambda\":"       << nan_or(lam)
                  << ",\"lambda_ratio_vs_ct\":" << nan_or(lr)
                  << ",\"coupled_balance\":"    << nan_or(bal) << "}";
            }
            j << "\n      }\n    }";
        };

        j << "  \"damage_channel_stats\": {\n";
        emit_dcs_end("5prime", dp.tri_5prime_pos);
        j << ",\n";
        emit_dcs_end("3prime", dp.tri_3prime_pos);
        j << "\n  },\n";
    }

    // ── Context-specific damage rates ─────────────────────────────────────────
    // For each of the 16 XCY (5' C→T) and 16 XGY (3' G→A) contexts:
    // terminal_rate = observed rate at read positions 1-4 from end
    // interior_rate = baseline at positions 10-14
    // excess = terminal_rate - interior_rate  (-1 = insufficient data)
    {
        static constexpr char B[4] = {'A','C','G','T'};
        // encoding: prev*16 + mid*4 + next  (A=0,C=1,G=2,T=3)
        j << "  \"trinuc_damage_rates\": {\n";

        auto emit_ctx_rates = [&](const char* key,
                                  const std::array<uint64_t,64>& term,
                                  const std::array<uint64_t,64>& intr,
                                  int mid_ref, int mid_dam) {
            j << "    \"" << key << "\": {";
            bool first = true;
            for (int p = 0; p < 4; ++p) {
                for (int n = 0; n < 4; ++n) {
                    int ir = p*16 + mid_ref*4 + n;
                    int id = p*16 + mid_dam*4 + n;
                    uint64_t tn = term[ir] + term[id];
                    uint64_t xn = intr[ir] + intr[id];
                    double tr = tn > 0 ? (double)term[id] / tn : -1.0;
                    double ir2= xn > 0 ? (double)intr[id] / xn : -1.0;
                    double ex = (tr >= 0 && ir2 >= 0) ? tr - ir2 : -1.0;
                    char ctx[4] = {B[p], B[mid_ref], B[n], 0};
                    if (!first) j << ",";
                    first = false;
                    j << "\"" << ctx << "\":{"
                      << "\"terminal_rate\":" << std::fixed << std::setprecision(6) << tr
                      << ",\"interior_rate\":" << ir2
                      << ",\"excess\":" << ex
                      << ",\"terminal_n\":" << tn
                      << ",\"interior_n\":" << xn << "}";
                }
            }
            j << "}";
        };

        emit_ctx_rates("ct_5prime",
                       dp.tri_5prime_terminal, dp.tri_5prime_interior,
                       1 /*C*/, 3 /*T*/);
        j << ",\n";
        emit_ctx_rates("ga_3prime",
                       dp.tri_3prime_terminal, dp.tri_3prime_interior,
                       2 /*G*/, 0 /*A*/);
        j << ",\n";
        emit_ctx_rates("ct_3prime",
                       dp.tri_3prime_terminal, dp.tri_3prime_interior,
                       1 /*C*/, 3 /*T*/);
        j << "\n  },\n";
    }

    // ── Context-modulated deamination profile ─────────────────────────────────
    // Reference-free per-NXN terminal deamination excess (terminal_rate −
    // interior_rate). Interior positions are the self-normalising background.
    // comparison_vector is positive-excess L1-normalised: encodes context shape
    // independent of total damage amplitude, for cross-library comparison.
    // SS libraries: 3′ C→T arm (dominant end in single-stranded prep).
    // DS libraries: mean of 5′ C→T and 3′ G→A RC-mapped to C-equivalent.
    {
        static constexpr char B4[4]  = {'A','C','G','T'};
        static constexpr int  RC4[4] = {3,2,1,0};  // A↔T, C↔G

        auto ctx_ex = [&](const std::array<uint64_t,64>& term,
                          const std::array<uint64_t,64>& intr,
                          int mid_ref, int mid_dam, int p, int n) -> double {
            int ir = p*16 + mid_ref*4 + n;
            int id = p*16 + mid_dam*4 + n;
            uint64_t tn = term[ir] + term[id];
            uint64_t xn = intr[ir] + intr[id];
            if (!tn || !xn) return 0.0;
            return (double)term[id]/tn - (double)intr[id]/xn;
        };

        // 16 canonical C-centred contexts: i = prev*4 + next  (p=i/4, n=i%4)
        // ACA ACC ACG ACT | CCA CCC CCG CCT | GCA GCC GCG GCT | TCA TCC TCG TCT
        double ex5[16], ex3c[16], ex3g_rc[16];
        for (int i = 0; i < 16; ++i) {
            int p = i/4, n = i%4;
            ex5[i]     = ctx_ex(dp.tri_5prime_terminal, dp.tri_5prime_interior, 1,3, p,n);
            ex3c[i]    = ctx_ex(dp.tri_3prime_terminal, dp.tri_3prime_interior, 1,3, p,n);
            // DS arm: 3′ G→A context RC(n)GRC(p) is the complement of canonical NCN (p,n)
            ex3g_rc[i] = ctx_ex(dp.tri_3prime_terminal, dp.tri_3prime_interior,
                                2,0, RC4[n],RC4[p]);
        }

        double cmp[16];
        for (int i = 0; i < 16; ++i)
            cmp[i] = is_ss ? ex3c[i] : 0.5*(ex5[i] + ex3g_rc[i]);

        // L1-positive normalisation (comparison_vector): shape of positive excess only.
        // Useful for DS where most channels are positive.
        double pos_sum = 0.0;
        for (int i = 0; i < 16; ++i) if (cmp[i] > 0.0) pos_sum += cmp[i];
        double norm[16];
        for (int i = 0; i < 16; ++i)
            norm[i] = (pos_sum > 0.0) ? std::max(0.0, cmp[i]) / pos_sum : 0.0;

        // L2-signed normalisation (signed_comparison_vector): preserves negative
        // channels (non-CpG in SS). Correct for ordination mixing SS and DS.
        double l2 = 0.0;
        for (int i = 0; i < 16; ++i) l2 += cmp[i]*cmp[i];
        l2 = std::sqrt(l2);
        double snorm[16];
        for (int i = 0; i < 16; ++i)
            snorm[i] = (l2 > 1e-9) ? cmp[i] / l2 : 0.0;

        auto emit16 = [&](const char* name, const double* v, bool comma) {
            j << "    \"" << name << "\": [";
            for (int i = 0; i < 16; ++i) {
                j << std::fixed << std::setprecision(6) << v[i];
                if (i < 15) j << ",";
            }
            j << "]" << (comma ? "," : "") << "\n";
        };

        j << "  \"deam_context_spectrum\": {\n";
        j << "    \"channels\": [";
        for (int i = 0; i < 16; ++i) {
            j << "\"" << B4[i/4] << "C" << B4[i%4] << "\"";
            if (i < 15) j << ",";
        }
        j << "],\n";
        emit16("ct_5prime_excess",          ex5,     true);
        emit16("ct_3prime_excess",          ex3c,    true);
        emit16("ga_3prime_rc_excess",       ex3g_rc, true);
        emit16("comparison_vector",         norm,    true);
        emit16("signed_comparison_vector",  snorm,   true);
        j << "    \"primary_arm\": \""
          << (is_ss ? "ct_3prime" : "ct_5prime_ga_3prime_mean") << "\",\n";
        j << "    \"sum_positive_excess\": "
          << std::fixed << std::setprecision(6) << pos_sum << "\n";
        j << "  },\n";
    }

    // ── Ancient-fraction d_max ─────────────────────────────────────────────────
    // d_max estimated only from reads classified as ancient by the per-read LLR
    // scorer (fused into the oxoG pass). Comparable to metaDMG "damaged fraction".
    j << "  \"ancient_fraction\": {\n";
    j << "    \"valid\": " << (dp.damaged_fraction_valid ? "true" : "false") << ",\n";
    // C1: d_max_5prime/3prime use the mixture identity d_anc = d_bulk / π, which amplifies
    // noise as π → 0. Below π=0.1 the identity estimate is unreliable → null + flag; the
    // OLS-fit fields stay the trustworthy source. *_fit/lambda are nulled when the fit failed.
    const bool identity_amplified =
        dp.damaged_fraction_pi > 0.0f && dp.damaged_fraction_pi < 0.1f;
    j << "    \"identity_amplified\": " << (identity_amplified ? "true" : "false") << ",\n";
    j << std::setprecision(6);
    j << "    \"ancient\": {\n";
    j << "      \"d_max_5prime\": " << (identity_amplified ? std::string("null") : nan_or(dp.damaged_fraction_d5)) << ",\n";
    j << "      \"d_max_3prime\": " << (identity_amplified ? std::string("null") : nan_or(dp.damaged_fraction_d3)) << ",\n";
    j << "      \"d_max_5prime_fit\": " << nan_or(dp.damaged_fraction_d5_fit) << ",\n";
    j << "      \"lambda_5prime\": " << nan_or(dp.damaged_fraction_lambda5) << ",\n";
    j << "      \"d_max_3prime_fit\": " << nan_or(dp.damaged_fraction_d3_fit) << ",\n";
    j << "      \"lambda_3prime\": " << nan_or(dp.damaged_fraction_lambda3) << ",\n";
    j << "      \"fraction\": " << dp.damaged_fraction_pi << ",\n";
    j << "      \"n_reads\": " << dp.damaged_fraction_n << ",\n";
    j << "      \"rate_5prime\": [";
    for (int p = 0; p < 15; ++p) { if (p) j << ","; j << dp.damaged_fraction_rate5[p]; }
    j << "],\n";
    j << "      \"rate_3prime\": [";
    for (int p = 0; p < 15; ++p) { if (p) j << ","; j << dp.damaged_fraction_rate3[p]; }
    j << "]\n";
    j << "    },\n";
    j << "    \"modern\": {\n";
    j << "      \"d_max_5prime\": " << (dp.modern_fraction_d5_computed ? std::to_string(dp.modern_fraction_d5) : "null") << ",\n";
    j << "      \"d_max_3prime\": " << (dp.modern_fraction_d3_computed ? std::to_string(dp.modern_fraction_d3) : "null") << ",\n";
    j << "      \"d_max_5prime_fit\": " << nan_or(dp.modern_fraction_d5_fit) << ",\n";
    j << "      \"lambda_5prime\": " << nan_or(dp.modern_fraction_lambda5) << ",\n";
    j << "      \"d_max_3prime_fit\": " << nan_or(dp.modern_fraction_d3_fit) << ",\n";
    j << "      \"lambda_3prime\": " << nan_or(dp.modern_fraction_lambda3) << ",\n";
    j << "      \"rate_5prime\": [";
    for (int p = 0; p < 15; ++p) { if (p) j << ","; j << dp.modern_fraction_rate5[p]; }
    j << "],\n";
    j << "      \"rate_3prime\": [";
    for (int p = 0; p < 15; ++p) { if (p) j << ","; j << dp.modern_fraction_rate3[p]; }
    j << "],\n";
    j << "      \"leakage_5prime\": " << (dp.modern_fraction_leakage_5prime ? "true" : "false") << ",\n";
    j << "      \"leakage_3prime\": " << (dp.modern_fraction_leakage_3prime ? "true" : "false") << "\n";
    j << "    }\n";
    j << "  },\n";

    // ── Preservation ──────────────────────────────────────────────────────────
    j << "  \"preservation\": {\n";
    j << "    \"score\": " << std::setprecision(6) << dp.preservation_score << ",\n";
    j << "    \"evidence\": " << dp.preservation_evidence << ",\n";
    j << "    \"reliability\": " << dp.preservation_reliability << ",\n";
    j << "    \"f5\": " << dp.preservation_f5 << ",\n";
    j << "    \"f3\": " << dp.preservation_f3 << ",\n";
    j << "    \"f_coh\": " << dp.preservation_f_coh << ",\n";
    j << "    \"f_cpg\": " << dp.preservation_f_cpg << ",\n";
    j << "    \"authenticity_eff\": " << std::setprecision(6) << pres.authenticity_eff << ",\n";
    j << "    \"authenticity_evidence\": " << pres.authenticity_evidence << ",\n";
    if (pres.d5_was_corrected)
        j << "    \"d5_hexamer_corrected\": " << pres.d5_hexamer_corrected << ",\n";
    j << "    \"oxidation_eff\": " << pres.oxidation_eff << ",\n";
    j << "    \"oxidation_evidence\": " << pres.oxidation_evidence << ",\n";
    j << "    \"oxidation\": {\n";
    j << "      \"raw_rate\": {\"estimate\": " << pres.oxidation.raw_rate.estimate
      << ", \"ci95_low\": " << pres.oxidation.raw_rate.ci95_low
      << ", \"ci95_high\": " << pres.oxidation.raw_rate.ci95_high << "},\n";
    j << "      \"control_rate\": {\"estimate\": " << pres.oxidation.control_rate.estimate
      << ", \"ci95_low\": " << pres.oxidation.control_rate.ci95_low
      << ", \"ci95_high\": " << pres.oxidation.control_rate.ci95_high << "},\n";
    j << "      \"excess_rate\": {\"estimate\": " << pres.oxidation.excess_rate.estimate
      << ", \"ci95_low\": " << pres.oxidation.excess_rate.ci95_low
      << ", \"ci95_high\": " << pres.oxidation.excess_rate.ci95_high << "},\n";
    j << "      \"z_score\": " << pres.oxidation.z_score << ",\n";
    j << "      \"bins_used\": " << pres.oxidation.bins_used << ",\n";
    j << "      \"effective_bins\": " << pres.oxidation.effective_bins << ",\n";
    j << "      \"heterogeneity\": " << pres.oxidation.heterogeneity << ",\n";
    j << "      \"reliability_score\": " << pres.oxidation.reliability_score << ",\n";
    j << "      \"reliability\": \"" << pres.oxidation.reliability << "\"\n";
    j << "    },\n";
    j << "    \"qc_risk_eff\": " << pres.qc_risk_eff << ",\n";
    j << "    \"qc_evidence\": " << pres.qc_evidence << "\n";
    j << "  },\n";

    // ── Library QC ────────────────────────────────────────────────────────────
    {
        j << "  \"library_qc\": {\n";
        if (in.adapter_clipped && !in.adapter_stubs_5prime.empty()) {
            j << "    \"adapter_stubs_clipped\": [";
            for (size_t i = 0; i < in.adapter_stubs_5prime.size(); ++i) {
                if (i) j << ",";
                j << "\"" << in.adapter_stubs_5prime[i] << "\"";
            }
            j << "],\n";
        }
        if (in.adapter3_clipped && !in.adapter_stubs_3prime.empty()) {
            j << "    \"adapter_stubs_clipped_3prime\": [";
            for (size_t i = 0; i < in.adapter_stubs_3prime.size(); ++i) {
                if (i) j << ",";
                j << "\"" << in.adapter_stubs_3prime[i] << "\"";
            }
            j << "],\n";
        }
        if (in.adapter_stub_reads_checked > 0) {
            j << "    \"adapter_stub5_read_fraction\": " << std::fixed << std::setprecision(6) << in.adapter_stub5_read_fraction << ",\n";
            j << "    \"adapter_stub3_read_fraction\": " << in.adapter_stub3_read_fraction << ",\n";
            j << "    \"adapter_stub_reads_checked\": " << in.adapter_stub_reads_checked << ",\n";
        }
        j << "    \"adapter_offset_5prime\": " << dp.fit_offset_5prime << ",\n";
        j << "    \"adapter_offset_3prime\": " << dp.fit_offset_3prime << ",\n";
        j << "    \"position0_artifact_5prime\": " << (dp.position_0_artifact_5prime ? "true" : "false") << ",\n";
        j << "    \"position0_artifact_3prime\": " << (dp.position_0_artifact_3prime ? "true" : "false") << ",\n";
        j << "    \"hexamer_entropy_5prime\": " << std::setprecision(4) << hs.entropy_terminal << ",\n";
        j << "    \"hexamer_entropy_interior\": " << hs.entropy_interior << ",\n";
        j << "    \"hexamer_terminal_interior_jsd\": " << hs.jsd << ",\n";
        j << "    \"hex_shift_g\": " << std::setprecision(4) << hs.shift_g << ",\n";
        j << "    \"hex_shift_z\": " << hs.shift_z << ",\n";
        j << "    \"hex_shift_p\": " << std::setprecision(6) << hs.shift_p << ",\n";
        j << "    \"hexamer_excess_tc\": " << std::setprecision(6) << dp.hexamer_excess_tc << ",\n";
        j << "    \"top_hexamers_5prime\": [";
        {
            int n_out = 0;
            for (const auto& hr : in.top_hex_enriched) {
                if (n_out >= 5) break;
                auto seq = decode_hex(hr.idx);
                if (n_out) j << ",";
                j << "{\"seq\":\"" << seq.data() << "\","
                  << "\"log2fc\":" << std::setprecision(3) << hr.log2fc << ","
                  << "\"damage_consistent\":" << (hr.damage_consistent ? "true" : "false") << "}";
                ++n_out;
            }
        }
        j << "],\n";
        j << "    \"top_hexamers_3prime\": [";
        {
            int n_out = 0;
            for (const auto& hr : in.top_hex_enriched_3prime) {
                if (n_out >= 5) break;
                auto seq = decode_hex(hr.idx);
                if (n_out) j << ",";
                j << "{\"seq\":\"" << seq.data() << "\","
                  << "\"log2fc\":" << std::setprecision(3) << hr.log2fc << ","
                  << "\"damage_consistent\":" << (hr.damage_consistent ? "true" : "false") << "}";
                ++n_out;
            }
        }
        j << "],\n";
        // Adapter prefix identification
        if (in.adapter_clipped && !in.adapter_stubs_5prime.empty()) {
            static const std::pair<const char*, const char*> kAdapters[] = {
                {"ACACTC", "TruSeq/P5"},
                {"AATGAT", "TruSeq/Universal"},
                {"GATCGG", "TruSeq/i7"},
                {"CTGTCT", "Nextera/Tn5"},
                {"AGATCG", "TruSeq/R1"},
                {"TGGAAT", "TruSeq/R2"},
                {"GCGAAT", "TruSeq/R2alt"},
            };
            const std::string& top_seq = in.adapter_stubs_5prime[0];
            const char* name = "unknown";
            for (const auto& kv : kAdapters)
                if (top_seq == kv.first) { name = kv.second; break; }
            j << "    \"adapter_prefix_identified\": {\"seq\":\"" << top_seq
              << "\",\"name\":\"" << name << "\"},\n";
        }
        j << "    \"depurination_detected\": " << (dp.depurination_detected ? "true" : "false") << ",\n";
        j << "    \"short_read_frac\": " << std::setprecision(4)
          << (in.short_read_frac < 0 ? 0.0 : in.short_read_frac) << ",\n";
        // ── diagnostic_groups ─────────────────────────────────────────────────
        // Convenience fractions derived from top enriched hexamer lists
        const auto& t5 = in.top_hex_enriched;
        const auto& t3 = in.top_hex_enriched_3prime;
        auto dmg_frac = [](const std::vector<HexEnrichment>& v, int k) -> double {
            if (v.empty() || k <= 0) return std::numeric_limits<double>::quiet_NaN();
            int n = std::min(k, (int)v.size()), m = 0;
            for (int i = 0; i < n; ++i) if (v[i].damage_consistent) ++m;
            return static_cast<double>(m) / n;
        };
        double dmg_frac_5 = dmg_frac(t5, 5);
        double dmg_frac_3 = dmg_frac(t3, 5);

        const HexEndAsymmetry& hea = in.hex_end_asymmetry;

        // Output certification for d_max_3:
        // confounded when ends are asymmetric (rc_overlap_topk==0 and hea.rc_excess_jsd high)
        // AND the 3' enriched hexamers are mostly non-damage-consistent.
        bool d3_confounded = (hea.rc_overlap_topk == 0)
                          && (!std::isnan(dmg_frac_3) && dmg_frac_3 < 0.5)
                          && (dp.fit_offset_3prime >= 1);
        bool d3_corrected  = !d3_confounded && dp.fit_offset_3prime >= 1;

        // Aliases so code below reads clearly (pre-computed at function top).
        const bool d5_suppressed        = d5_suppressed_early;
        const bool d3_selection_biased  = d3_selection_biased_early;
        // Both ends unusable: d5 blunted AND d3 also flat (no real decay anywhere).
        const bool no_reliable_estimate = is_ss
                                       && (hea.rc_overlap_topk == 0)
                                       && d5_suppressed
                                       && dp.d3_profile_flat;
        const bool ss_extreme_asym = is_ss && (hea.rc_overlap_topk == 0)
                                  && ((dp.d_max_5prime < 0.01f && dp.d_max_3prime > 0.05f)
                                   || (dp.d_max_3prime < 0.01f && dp.d_max_5prime > 0.05f));

        j << "    \"diagnostic_groups\": {\n";

        // adapter_position_effects
        j << "      \"adapter_position_effects\": {\n";
        j << "        \"adapter_offset_5prime\": " << dp.fit_offset_5prime << ",\n";
        j << "        \"adapter_offset_3prime\": " << dp.fit_offset_3prime << ",\n";
        j << "        \"position0_artifact_5prime\": " << (dp.position_0_artifact_5prime ? "true" : "false") << ",\n";
        j << "        \"position0_artifact_3prime\": " << (dp.position_0_artifact_3prime ? "true" : "false") << ",\n";
        if (in.adapter_stub_reads_checked > 0) {
            j << "        \"prefix_read_fraction_5prime\": " << std::fixed << std::setprecision(6) << in.adapter_stub5_read_fraction << ",\n";
            j << "        \"prefix_read_fraction_3prime\": " << in.adapter_stub3_read_fraction << ",\n";
        } else {
            j << "        \"prefix_read_fraction_5prime\": null,\n";
            j << "        \"prefix_read_fraction_3prime\": null,\n";
        }
        j << "        \"corrected_outputs\": [";
        {
            bool first = true;
            auto ca = [&](const char* s) { if (!first) j << ","; j << "\"" << s << "\""; first = false; };
            if (dp.fit_offset_5prime > 1) ca("d_max_5");
            if (dp.fit_offset_3prime > 1 && !d3_confounded) ca("d_max_3");
        }
        j << "],\n";
        j << "        \"residual_outputs\": [\"position0_base_composition\"],\n";
        j << "        \"confounded_outputs\": [" << (d3_confounded ? "\"d_max_3\"" : "") << "],\n";
        j << "        \"selection_biased_outputs\": [" << (d3_selection_biased ? "\"d_max_3\"" : "") << "],\n";
        const bool ds_bilateral_inversion_unreliable =
            !is_ss
            && dp.inverted_pattern_5prime && dp.inverted_pattern_3prime
            && !dp.position_0_artifact_5prime && !dp.position_0_artifact_3prime
            && dp.damage_validated;
        j << "        \"suspect_outputs\": ["
          << ((ss_extreme_asym || ds_bilateral_inversion_unreliable) ? "\"d_max_combined\"" : "")
          << "]\n";
        j << "      },\n";

        // terminal_hexamer_bias
        j << "      \"terminal_hexamer_bias\": {\n";
        j << "        \"hexamer_terminal_interior_jsd\": " << std::setprecision(6) << hs.jsd << ",\n";
        j << "        \"hexamer_excess_tc\": " << dp.hexamer_excess_tc << ",\n";
        j << "        \"top_damage_consistent_fraction_5prime\": ";
        if (std::isnan(dmg_frac_5)) j << "null"; else j << std::setprecision(4) << dmg_frac_5;
        j << ",\n";
        j << "        \"top_damage_consistent_fraction_3prime\": ";
        if (std::isnan(dmg_frac_3)) j << "null"; else j << std::setprecision(4) << dmg_frac_3;
        j << ",\n";
        j << "        \"residual_outputs\": [\"terminal_NXN_context_spectrum\",\"top_hexamers\"]\n";
        j << "      },\n";

        // end_hexamer_asymmetry
        j << "      \"end_hexamer_asymmetry\": {\n";
        j << "        \"hexamer_end_rc_excess_jsd_status\": \"" << hea.status << "\",\n";
        j << "        \"hexamer_end_rc_excess_jsd\": ";
        if (std::isnan(hea.rc_excess_jsd)) j << "null"; else j << std::setprecision(6) << hea.rc_excess_jsd;
        j << ",\n";
        j << "        \"hexamer_end_fwd_excess_jsd\": ";
        if (std::isnan(hea.fwd_excess_jsd)) j << "null"; else j << std::setprecision(6) << hea.fwd_excess_jsd;
        j << ",\n";
        j << "        \"rc_overlap_topk\": " << hea.rc_overlap_topk << ",\n";
        j << "        \"topk\": " << hea.topk << ",\n";
        j << "        \"terminal_excess_mass_5prime\": " << std::setprecision(6) << hea.excess_mass_5prime << ",\n";
        j << "        \"terminal_excess_mass_3prime\": " << hea.excess_mass_3prime << ",\n";
        j << "        \"terminal_hexamer_n_5prime\": " << hea.n_5prime << ",\n";
        j << "        \"terminal_hexamer_n_3prime\": " << hea.n_3prime << "\n";
        j << "      },\n";

        // end_damage_profile: per-position decay flatness diagnostics.
        // d5_profile_flat / d3_profile_flat indicate the end signal is indistinguishable
        // from noise — likely library construction erasure, not genuinely zero damage.
        j << "      \"end_damage_profile\": {\n";
        j << "        \"d5_profile_flat\": " << (dp.d5_profile_flat ? "true" : "false") << ",\n";
        j << "        \"d3_profile_flat\": " << (dp.d3_profile_flat ? "true" : "false") << ",\n";
        j << std::setprecision(6);
        j << "        \"d5_max_rate_pos0_4\": " << dp.d5_max_rate_pos0_4 << ",\n";
        j << "        \"d3_max_rate_pos0_4\": " << dp.d3_max_rate_pos0_4 << ",\n";
        j << "        \"d5_var_pos0_9\": " << dp.d5_profile_var_pos0_9 << ",\n";
        j << "        \"d3_var_pos0_9\": " << dp.d3_profile_var_pos0_9 << ",\n";
        j << "        \"d5_amp_over_bg\": " << (dp.fit_baseline_5prime > 1e-4f
              ? dp.fit_amplitude_5prime / dp.fit_baseline_5prime : 0.0f) << ",\n";
        j << "        \"d3_amp_over_bg\": " << (dp.fit_baseline_3prime > 1e-4f
              ? dp.fit_amplitude_3prime / dp.fit_baseline_3prime : 0.0f) << ",\n";
        j << "        \"d5_blunting_suspected\": " << (dp.d5_blunting_suspected ? "true" : "false") << "\n";
        j << "      }" << (is_ss ? "," : "") << "\n";

        // ss_end_asymmetry: CircLigase selection bias and extreme 5'/3' asymmetry.
        // Only emitted for SS libraries (CircLigase is the relevant ligation chemistry).
        if (is_ss) {
            const char* recommended =
                no_reliable_estimate ? "none" :
                d5_suppressed        ? "d_max_3prime" :
                d3_selection_biased  ? "d_max_5prime" :
                (dp.d_max_5prime < 0.01f && dp.d_max_3prime > 0.05f) ? "d_max_3prime" :
                (dp.d_max_3prime < 0.01f && dp.d_max_5prime > 0.05f) ? "d_max_5prime" :
                "none";

            const char* note =
                no_reliable_estimate ? "d5 suppressed by 5' blunting and d3 selection-biased; no reliable population estimate" :
                d5_suppressed        ? "d5 flat — 5' overhang likely blunted before adapter ligation; d3 is best available" :
                d3_selection_biased  ? "d3 inflated by CircLigase selection of deaminated 3' termini; d5 is conservative estimate" :
                ss_extreme_asym      ? "extreme 5'/3' asymmetry in SS library; interpret both ends with caution" :
                "";

            j << "      \"ss_end_asymmetry\": {\n";
            j << "        \"d5_suppressed\": " << (d5_suppressed ? "true" : "false") << ",\n";
            j << "        \"d3_selection_biased\": " << (d3_selection_biased ? "true" : "false") << ",\n";
            j << "        \"no_reliable_estimate\": " << (no_reliable_estimate ? "true" : "false") << ",\n";
            j << "        \"recommended_estimate\": \"" << recommended << "\",\n";
            j << "        \"note\": \"" << note << "\"\n";
            j << "      }\n";
        }

        j << "    },\n";

        // output_effects summary (d_max_3 certification)
        j << "    \"output_effects\": {\n";
        j << "      \"d_max_3\": {\n";
        {
            const char* d3_status = d3_confounded       ? "confounded" :
                                    d3_selection_biased  ? "selection_biased" :
                                    d3_corrected         ? "corrected" : "residual";
            j << "        \"status\": \"" << d3_status << "\",\n";
            j << "        \"evidence\": [";
            if (d3_confounded) {
                j << "\"top_damage_consistent_fraction_3prime\",\"rc_overlap_topk\"";
            } else if (d3_selection_biased) {
                j << "\"rc_overlap_topk\",\"d_max_3prime_over_d_max_5prime\"";
            } else if (d3_corrected) {
                j << "\"adapter_offset_3prime\"";
            }
            j << "]\n";
        }
        j << "      }\n";
        j << "    },\n";

        j << "    \"flags\": []\n";
        j << "  },\n";
    }

    // ── Damage context ────────────────────────────────────────────────────────
    {
        j << "  \"damage_context\": {\n";
        j << "    \"method\": \"" << dcp.method << "\",\n";
        j << "    \"reference_required\": " << (dcp.reference_required ? "true" : "false") << ",\n";
        j << "    \"alignment_required\": " << (dcp.alignment_required ? "true" : "false") << ",\n";
        j << "    \"dominant_process\": \"" << dcp.dominant_process_str << "\",\n";
        j << "    \"interpretation\": \"" << json_escape(dcp.interpretation) << "\",\n";
        auto emit_score = [&](const char* name, float v, bool trailing) {
            j << "    \"" << name << "\": ";
            if (std::isnan(v)) j << "null"; else j << std::setprecision(6) << v;
            j << (trailing ? ",\n" : "\n");
        };
        emit_score("terminal_deamination_score",  dcp.terminal_deamination_score,  true);
        emit_score("cpg_context_score",           dcp.cpg_context_score,           true);
        emit_score("dipyrimidine_context_score",  dcp.dipyrimidine_context_score,  true);
        emit_score("oxidative_context_score",     dcp.oxidative_context_score,     true);
        emit_score("fragmentation_context_score", dcp.fragmentation_context_score, true);
        emit_score("library_artifact_score",      dcp.library_artifact_score,      true);
        j << "    \"evidence\": {\n";
        auto emit_ef = [&](const char* name, float v, bool trailing) {
            j << "      \"" << name << "\": ";
            if (std::isnan(v) || std::isinf(v)) j << "null";
            else j << std::setprecision(6) << v;
            j << (trailing ? ",\n" : "\n");
        };
        emit_ef("d_max_5",               dcp.evidence.d_max_5,               true);
        emit_ef("d_max_3",               dcp.evidence.d_max_3,               true);
        emit_ef("lambda_5",              dcp.evidence.lambda_5,              true);
        emit_ef("lambda_3",              dcp.evidence.lambda_3,              true);
        emit_ef("log2_cpg_ratio",        dcp.evidence.log2_cpg_ratio,        true);
        emit_ef("cpg_z",                 dcp.evidence.cpg_z,                 true);
        emit_ef("dipyr_contrast",        dcp.evidence.dipyr_contrast,        true);
        emit_ef("ox_gt_asymmetry",       dcp.evidence.ox_gt_asymmetry,       true);
        emit_ef("s_oxog_mean",           dcp.evidence.s_oxog_mean,           true);
        emit_ef("s_oxog_max",            dcp.evidence.s_oxog_max,            true);
        emit_ef("purine_enrichment_5prime", dcp.evidence.purine_enrichment_5prime, true);
        emit_ef("hex_shift_z",           dcp.evidence.hex_shift_z,           true);
        j << "      \"adapter_clipped\": " << (dcp.evidence.adapter_clipped ? "true" : "false") << ",\n";
        j << "      \"adapter3_clipped\": " << (dcp.evidence.adapter3_clipped ? "true" : "false") << ",\n";
        j << "      \"flag_hex_artifact\": " << (dcp.evidence.flag_hex_artifact ? "true" : "false") << ",\n";
        j << "      \"position_0_artifact_5prime\": " << (dcp.evidence.position_0_artifact_5prime ? "true" : "false") << ",\n";
        j << "      \"position_0_artifact_3prime\": " << (dcp.evidence.position_0_artifact_3prime ? "true" : "false") << ",\n";
        j << "      \"fit_offset_5prime\": " << dcp.evidence.fit_offset_5prime << ",\n";
        j << "      \"fit_offset_3prime\": " << dcp.evidence.fit_offset_3prime << ",\n";
        j << "      \"n_reads\": " << dcp.evidence.n_reads << "\n";
        j << "    }\n";
        j << "  },\n";
    }

    // ── BIC tournament metadata ───────────────────────────────────────────────
    j << "  \"bic\": {\n";
    j << "    \"bias\": " << std::setprecision(2) << dp.library_bic_bias << ",\n";
    j << "    \"ds\": " << dp.library_bic_ds << ",\n";
    j << "    \"ss\": " << dp.library_bic_ss << ",\n";
    j << "    \"ct5_amp\": " << std::setprecision(6) << dp.libtype_amp_ct5 << ",\n";
    j << "    \"ga3_amp\": " << dp.libtype_amp_ga3 << ",\n";
    j << "    \"ga0_amp\": " << dp.libtype_amp_ga0 << ",\n";
    j << "    \"ct3_amp\": " << dp.libtype_amp_ct3 << ",\n";
    j << "    \"p_ds\": " << dp.library_p_ds << ",\n";
    j << "    \"p_ss\": " << dp.library_p_ss << ",\n";
    j << "    \"p_bias\": " << dp.library_p_bias << ",\n";
    j << "    \"p_winner\": " << dp.library_p_winner << ",\n";
    j << "    \"p_ds_final\": " << dp.library_p_ds_final << ",\n";
    j << "    \"p_ss_final\": " << dp.library_p_ss_final << ",\n";
    j << "    \"p_bias_final\": " << dp.library_p_bias_final << ",\n";
    j << "    \"p_winner_final\": " << dp.library_p_winner_final << ",\n";
    j << "    \"evaluable\": " << (dp.library_type_evaluable ? "true" : "false") << ",\n";
    // C1: a failed-validity submodel stores kInvalidBIC (1e300); emitted raw it is a valid JSON
    // number (10^299) that external softmax/plot/threshold tools cannot tell from a real BIC.
    // Emit null for any non-finite or sentinel-magnitude (>1e200) value.
    auto jbic = [&](double v) {
        if (std::isfinite(v) && std::abs(v) < 1e200) j << std::setprecision(2) << v;
        else                                         j << "null";
    };
    j << "    \"library_submodel_bic\": {\n";
    j << "      \"M_bias\": ";        jbic(dp.library_bic_M_bias);        j << ",\n";
    j << "      \"M_DS_symm\": ";     jbic(dp.library_bic_M_DS_symm);     j << ",\n";
    j << "      \"M_DS_spike\": ";    jbic(dp.library_bic_M_DS_spike);    j << ",\n";
    j << "      \"M_DS_symm_art\": "; jbic(dp.library_bic_M_DS_symm_art); j << ",\n";
    j << "      \"M_SS_comp\": ";     jbic(dp.library_bic_M_SS_comp);     j << ",\n";
    j << "      \"M_SS_orig\": ";     jbic(dp.library_bic_M_SS_orig);     j << ",\n";
    j << "      \"M_SS_asym\": ";     jbic(dp.library_bic_M_SS_asym);     j << ",\n";
    j << "      \"M_SS_full\": ";     jbic(dp.library_bic_M_SS_full);     j << ",\n";
    j << "      \"M_DS_asym_art\": "; jbic(dp.library_bic_M_DS_asym_art); j << "\n";
    j << "    },\n";
    j << "    \"library_bic_winner_model\": \"" << dp.library_bic_winner_model << "\",\n";
    j << "    \"library_bic_second_model\": \"" << dp.library_bic_second_model << "\",\n";
    j << "    \"library_bic_margin\": " << dp.library_bic_margin << ",\n";
    j << "    \"protocol_tag_5prime\": \""   << dp.protocol_tag_5prime   << "\",\n";
    j << "    \"protocol_tag_protocol\": \"" << dp.protocol_tag_protocol << "\",\n";
    j << "    \"protocol_tag_class\": \""    << libtype_human(dp.protocol_tag_class) << "\",\n";
    j << "    \"protocol_tag_log2fc\": "     << std::setprecision(4) << dp.protocol_tag_log2fc << ",\n";
    j << "    \"protocol_tag_log_lr\": "     << dp.protocol_tag_log_lr << ",\n";
    j << "    \"protocol_tag_applied\": "    << (dp.protocol_tag_applied ? "true" : "false") << ",\n";
    j << "    \"briggs_pos0_masked_5prime\": " << (dp.briggs_pos0_masked_5prime ? "true" : "false") << ",\n";
    j << "    \"briggs_pos0_masked_3prime\": " << (dp.briggs_pos0_masked_3prime ? "true" : "false") << ",\n";
    // Wave-3: ss 5' single-strand overhang kernel status (bulk damage model). modeled=true when the
    // terminus was included as overhang (r(0)=1); degenerate=true when an ss library fell back to the
    // ds exp kernel because the overhang was not identifiable (p0 artifact). Both false for ds.
    j << "    \"ss_overhang_modeled\": " << (dp.ss_overhang_modeled ? "true" : "false") << ",\n";
    j << "    \"ss_overhang_degenerate\": " << (dp.ss_overhang_degenerate ? "true" : "false") << ",\n";
    j << "    \"bg_5prime_anchored\": " << std::setprecision(6) << dp.bg_5prime_anchored << ",\n";
    j << "    \"bg_3prime_anchored\": " << dp.bg_3prime_anchored << ",\n";
    j << "    \"bg_n_positions_5prime\": " << dp.bg_n_positions_5prime << ",\n";
    j << "    \"bg_n_positions_3prime\": " << dp.bg_n_positions_3prime << ",\n";
    j << "    \"bg_denominator_5prime\": " << static_cast<long long>(dp.bg_denominator_5prime) << ",\n";
    j << "    \"bg_denominator_3prime\": " << static_cast<long long>(dp.bg_denominator_3prime) << ",\n";
    j << "    \"damage_5prime_area_excess\": " << std::setprecision(6) << dp.damage_5prime_area_excess << ",\n";
    j << "    \"damage_3prime_area_excess\": " << dp.damage_3prime_area_excess << ",\n";
    j << "    \"damage_5prime_lr\": " << std::setprecision(4) << dp.damage_5prime_lr << ",\n";
    j << "    \"damage_3prime_lr\": " << dp.damage_3prime_lr << ",\n";
    j << "    \"input_mode\": \""
      << (dp.input_mode == SP::InputMode::PAIRED ? "paired" : "single") << "\",\n";
    j << "    \"pe_short_insert_skipped\": " << static_cast<long long>(dp.pe_short_insert_skipped) << ",\n";
    j << "    \"library_bic_M_DS_asym_art\": " << dp.library_bic_M_DS_asym_art << ",\n";
    j << "    \"library_bic_M_DS_symm_art_no_offset\": " << dp.library_bic_M_DS_symm_art_no_offset << ",\n";
    j << "    \"library_bic_raw_winner_model\": \"" << dp.library_bic_raw_winner_model << "\",\n";
    j << "    \"library_bic_raw_winner_class\": \"" << dp.library_bic_raw_winner_class << "\",\n";
    j << "    \"library_bic_raw_second_model\": \"" << dp.library_bic_raw_second_model << "\",\n";
    j << "    \"library_bic_raw_margin\": " << dp.library_bic_raw_margin << ",\n";
    j << "    \"library_bic_raw_winner_in_cascade\": "
      << (dp.library_bic_raw_winner_in_cascade ? "true" : "false") << ",\n";
    j << "    \"library_bic_excl_structural_bilateral\": "
      << (dp.library_bic_excl_structural_bilateral ? "true" : "false") << ",\n";
    j << "    \"library_bic_excl_spike_is_ss\": "
      << (dp.library_bic_excl_spike_is_ss ? "true" : "false") << ",\n";
    j << "    \"library_bic_excl_ct3_zero\": "
      << (dp.library_bic_excl_ct3_zero ? "true" : "false") << ",\n";
    j << "    \"library_bic_excl_M_SS_full_hardcoded\": "
      << (dp.library_bic_excl_M_SS_full_hardcoded ? "true" : "false") << ",\n";
    j << "    \"library_bic_excl_in_cascade\": "
      << (dp.library_bic_excl_in_cascade ? "true" : "false") << ",\n";
    j << "    \"final_library_bic_winner_model\": \"" << dp.final_library_bic_winner_model << "\",\n";
    j << "    \"final_library_bic_override_reason\": \"" << dp.final_library_bic_override_reason << "\",\n";
    j << "    \"library_gate_artifact_5\": " << (dp.library_gate_artifact_5 ? "true" : "false") << ",\n";
    j << "    \"library_gate_artifact_3\": " << (dp.library_gate_artifact_3 ? "true" : "false") << ",\n";
    j << "    \"library_gate_position_0_artifact_5prime\": "
      << (dp.library_gate_position_0_artifact_5prime ? "true" : "false") << ",\n";
    j << "    \"library_gate_position_0_artifact_3prime\": "
      << (dp.library_gate_position_0_artifact_3prime ? "true" : "false") << ",\n";
    j << "    \"library_gate_inverted_pattern_5prime\": "
      << (dp.library_gate_inverted_pattern_5prime ? "true" : "false") << ",\n";
    j << "    \"library_gate_inverted_pattern_3prime\": "
      << (dp.library_gate_inverted_pattern_3prime ? "true" : "false") << ",\n";
    j << "    \"library_gate_max_damage_5prime\": " << std::setprecision(6)
      << dp.library_gate_max_damage_5prime << ",\n";
    j << "    \"library_gate_ss_orientation_evidence\": "
      << (dp.library_gate_ss_orientation_evidence ? "true" : "false") << ",\n";
    j << "    \"library_gate_ga_spike_dominant\": "
      << (dp.library_gate_ga_spike_dominant ? "true" : "false") << ",\n";
    j << "    \"library_gate_ds_spike_won\": "
      << (dp.library_gate_ds_spike_won ? "true" : "false") << ",\n";
    j << "    \"library_gate_ga0_dominates_ct5\": "
      << (dp.library_gate_ga0_dominates_ct5 ? "true" : "false") << ",\n";
    j << "    \"library_gate_structural_bilateral\": "
      << (dp.library_gate_structural_bilateral ? "true" : "false") << ",\n";
    j << "    \"library_joint_lambda_restricted\": "
      << (dp.library_joint_lambda_restricted ? "true" : "false") << ",\n";
    j << "    \"fit_offset_5prime\": " << dp.fit_offset_5prime << ",\n";
    j << "    \"fit_offset_3prime\": " << dp.fit_offset_3prime << ",\n";
    j << "    \"ct3_offset\": " << dp.library_ct3_offset << ",\n";
    j << "    \"ds_symm_offset\": " << dp.library_ds_symm_offset << ",\n";
    j << "    \"ds_symm_lambda_used\": " << std::setprecision(6) << dp.library_ds_symm_lambda_used << ",\n";
    j << "    \"ds_symm_amp\": " << dp.library_ds_symm_amp << ",\n";
    j << "    \"ds_symm_ct5_resid\": " << dp.library_ds_symm_ct5_resid << ",\n";
    j << "    \"ds_symm_ga3_resid\": " << dp.library_ds_symm_ga3_resid << "\n";
    j << "  },\n";

    // ── Bulk damage law (Phase 1): threshold-free δ(L) ─────────────────────────
    // Count-level binomial GLM over data-driven length bins. Emits continuous,
    // uncertainty-bearing quantities only (δ̂, K_eff, R_damage, Δℓ0/S0 and the
    // FULL profile-likelihood curve) — no cutoffs, no verdict.
    {
        const auto& bd = dp.bulk_damage;
        // Finite-safe number emitter: a non-finite fit value (NaN/inf) prints as a bare
        // `nan`/`inf` token under std::fixed and would corrupt the whole document — emit `null`
        // instead (the validator treats null as non-finite and flags it loudly, the correct signal).
        auto jn = [&](double v, int prec = 6) {
            if (std::isfinite(v)) j << std::setprecision(prec) << v;
            else                  j << "null";
        };
        j << "  \"bulk_damage\": {\n";
        j << "    \"attempted\": " << (dp.bulk_attempted ? "true" : "false") << ",\n";
        j << "    \"valid\": "      << (bd.valid ? "true" : "false") << ",\n";
        j << "    \"converged\": "  << (bd.converged ? "true" : "false") << ",\n";
        j << "    \"n_sweeps\": "   << bd.n_sweeps << ",\n";
        j << "    \"log_lik\": ";        jn(bd.log_lik);              j << ",\n";
        j << "    \"lambda\": ";         jn(bd.lambda);               j << ",\n";
        j << "    \"lambda_at_boundary\": " << (bd.lambda_at_boundary ? "true" : "false") << ",\n";
        j << "    \"headline_delta\": "; jn(dp.bulk_headline_delta);  j << ",\n";
        // threshold-free length-coupling weight: w_length∈[0,1] = P(terminal-damage mass falls with read
        // length); slope_m = OLS slope of m(L)=δ_L·L (≤0 authentic terminal decay, ≫0 pervasive artifact).
        // δ_auth = δ·w_length.
        j << "    \"slope_m\": ";        jn(bd.slope_m);              j << ",\n";
        j << "    \"w_length\": ";       jn(bd.w_length);             j << ",\n";
        j << "    \"headline_delta_auth\": "; jn(dp.bulk_headline_delta * bd.w_length); j << ",\n";
        // length-coupling AGE/AUTHENTICITY indicator (NOT a damage estimate): sign of OLS slope of the
        // recovered δ_l vs read length. −1 = δ falls with length (fragmented ancient DNA, the classic
        // aDNA signature); 0 = flat (young or contamination-dominated); +1 = rises (artifact-suspect).
        // Read ALONGSIDE: deamination.d_max_5/3prime (ABSOLUTE terminal C→T — catches flat-young damage
        // like 1610 CE samples that the length-coupled bulk δ conservatively reports as 0), and the
        // library_artifact_* flags. The bulk δ is the length-coupled, artifact-screened view; d_max is
        // the absolute view; together they are complementary (each catches libraries the other misses).
        j << "    \"length_coupling\": ";       jn(static_cast<double>(bd.length_coupling)); j << ",\n";
        j << "    \"length_coupling_slope\": "; jn(bd.length_coupling_slope);                j << ",\n";
        // mixture P2: reference-free per-ancient-read deamination intensity (analog of metaDMG A_b);
        // −1 ⇒ undetermined (no usable co-occurrence signal or artifact-gated). d_max_se = its SE.
        // C1: gate on the producer's validity flags (bulk_damage_model.hpp). When the mixture
        // split is not identified d_max_ancient=-1.0 and d_max_se=0.0 are NOT-COMPUTED sentinels
        // (se=0.0 deceptively reads as a zero-width CI) — emit null. d_max_raw>1 is non-physical
        // (railed); emit null and surface the railed flag so the state is machine-readable.
        if (bd.d_max_ancient_valid) { j << "    \"d_max_ancient\": "; jn(bd.d_max_ancient); j << ",\n";
                                      j << "    \"d_max_se\": ";      jn(bd.d_max_se);      j << ",\n"; }
        else                        { j << "    \"d_max_ancient\": null,\n";
                                      j << "    \"d_max_se\": null,\n"; }
        j << "    \"d_max_ancient_valid\": " << (bd.d_max_ancient_valid ? "true" : "false") << ",\n";
        if (bd.d_max_raw_railed) j << "    \"d_max_raw\": null,\n";
        else                   { j << "    \"d_max_raw\": "; jn(bd.d_max_raw); j << ",\n"; }
        j << "    \"d_max_raw_railed\": " << (bd.d_max_raw_railed ? "true" : "false") << ",\n";
        j << "    \"kappa\": [";
        for (int r = 0; r < 2; ++r) {
            j << "[";
            for (int c = 0; c < 2; ++c) {
                jn(bd.kappa[r][c]);
                if (c < 1) j << ", ";
            }
            j << "]";
            if (r < 1) j << ", ";
        }
        j << "],\n";
        j << "    \"bins\": [";
        for (size_t b = 0; b < bd.bins.size(); ++b) {
            const auto& bb = bd.bins[b];
            if (b > 0) j << ",";
            j << "\n      {";
            j << "\"length_lo\": "    << bb.length_lo
              << ", \"length_hi\": "  << bb.length_hi
              << ", \"median_len\": "; jn(bb.median_len, 4);
            j << ", \"n_reads\": "    << bb.n_reads
              << ", \"delta\": ";     jn(bb.delta);
            j << ", \"delta_auth\": "; jn(bb.delta * bd.w_length);
            j << ", \"identified\": " << (bb.identified ? "true" : "false")
              << ", \"borrowed\": "   << (bb.borrowed ? "true" : "false");
            j << ", \"k_eff\": [";    jn(bb.k_eff[0], 4);  j << ", ";  jn(bb.k_eff[1], 4);  j << "]";
            j << ", \"r_damage\": ["; jn(bb.r_damage[0]);  j << ", ";  jn(bb.r_damage[1]);  j << "]";
            j << ", \"delta_ell0\": "; jn(bb.delta_ell0);
            j << ", \"s0\": ";         jn(bb.s0);
            j << ", \"profile_delta\": [";
            for (size_t p = 0; p < bb.profile_delta.size(); ++p) {
                jn(bb.profile_delta[p]);
                if (p + 1 < bb.profile_delta.size()) j << ", ";
            }
            j << "], \"profile_loglik\": [";
            for (size_t p = 0; p < bb.profile_loglik.size(); ++p) {
                jn(bb.profile_loglik[p]);
                if (p + 1 < bb.profile_loglik.size()) j << ", ";
            }
            // mixture P2: eligibility-conditioned ancient/modern split. pi_ancient = δ_l/d_max (ancient
            // read fraction); [pi_lo,pi_hi] = 95% interval (−1 ⇒ undetermined). joint_cov is the raw
            // (GC-confounded) co-occurrence diagnostic; the verdict uses the conditioned d_max/π.
            j << "], \"joint_n\": " << static_cast<long long>(bb.joint_n);
            j << ", \"joint_mean\": ["; jn(bb.joint_mean[0], 4); j << ", "; jn(bb.joint_mean[1], 4); j << "]";
            j << ", \"joint_cov\": "; jn(bb.joint_cov, 6);
            // C1: gate on the producer's per-bin pi_identified flag. The -1.0 defaults are
            // out-of-range for a probability and would read as a negative fraction; emit null.
            if (bb.pi_identified) {
                j << ", \"pi_ancient\": "; jn(bb.pi_ancient, 4);
                j << ", \"pi_lo\": ";      jn(bb.pi_lo, 4);
                j << ", \"pi_hi\": ";      jn(bb.pi_hi, 4);
            } else {
                j << ", \"pi_ancient\": null, \"pi_lo\": null, \"pi_hi\": null";
            }
            j << ", \"pi_identified\": " << (bb.pi_identified ? "true" : "false");
            j << "}";
        }
        if (!bd.bins.empty()) j << "\n    ";
        j << "]\n";
        j << "  },\n";
    }

    // ── Damage types ──────────────────────────────────────────────────────────
    {
        auto jbool = [](bool v) -> const char* { return v ? "true" : "false"; };
        auto jfloat = [&](double v) { j << std::setprecision(6) << v; };

        double cf_ratio = (dp.channel_c_valid && dp.channel_f_valid
                           && dp.ca_stop_rate_baseline > 1e-6f)
            ? static_cast<double>(dp.ox_stop_conversion_rate_baseline)
              / static_cast<double>(dp.ca_stop_rate_baseline)
            : -1.0;

        // C5: require channel_f_z (pooled) and channel_f_mh_z (MH stratified) to agree in sign
        // before firing detection. The pooled z is confounded by terminal context-composition
        // differences (Simpson's paradox); MH removes that bias. Mirroring the channel_h OR→AND fix.
        bool ch_f_z_consistent = std::isfinite(dp.channel_f_z) && std::isfinite(dp.channel_f_mh_z)
                                 && ((dp.channel_f_z >= 0.0f) == (dp.channel_f_mh_z >= 0.0f));
        bool ch_f_detected = dp.channel_f_valid && dp.channel_f_z > kOxChannelZDetect
                             && ch_f_z_consistent;
        bool ch_g_detected = dp.channel_g_valid && dp.channel_g_z > kOxChannelZDetect;
        // C5: OR→AND. channel_h_z_p2plus (positions 2-4 only) exists to exclude the p0/p1
        // artifact; the OR gate let an artifact-driven positive channel_h_z fire detection while
        // p2plus was strongly negative (kapk, rocs). Require BOTH windows to agree. h_z_consistent
        // surfaces sign agreement so callers see the (now-required) consistency. [behavioral change]
        bool ch_h_detected = dp.channel_h_valid &&
            dp.channel_h_z > kOxChannelZDetect && dp.channel_h_z_p2plus > kOxChannelZDetect;
        bool ch_h_z_consistent =
            (dp.channel_h_z >= 0.0f) == (dp.channel_h_z_p2plus >= 0.0f);
        bool ch_d_detected = std::abs(dp.ox_gt_asymmetry) > 0.01f;

        j << "  \"damage_types\": [\n";

        // Channel A
        j << "    {\n";
        j << "      \"channel\": \"A\",\n";
        j << "      \"name\": \"cytosine_deamination\",\n";
        j << "      \"description\": \"C to T post-mortem deamination at terminal positions\",\n";
        j << "      \"mechanism\": \"hydrolytic_deamination\",\n";
        // C5: `detected` was a tautological alias for damage_validated (the joint-model verdict),
        // which returns False on heavily-damaged libraries (d_max>10%) due to the position-0
        // artifact interaction — contradicting d_max_5prime/3prime printed right beside it. Define
        // `detected` as a direct effect-size gate (consistent with damage_status), and keep the
        // strict joint-model verdict as joint_model_validated. [behavioral change to `detected`]
        bool ch_a_detected = (dp.d_max_5prime > 0.02f || dp.d_max_3prime > 0.02f);
        j << "      \"detected\": " << jbool(ch_a_detected) << ",\n";
        j << "      \"d_max_5prime\": "; jfloat(dp.d_max_5prime); j << ",\n";
        j << "      \"d_max_3prime\": "; jfloat(dp.d_max_3prime); j << ",\n";
        j << "      \"joint_model_validated\": " << jbool(dp.damage_validated) << ",\n";
        j << "      \"validated\": " << jbool(dp.damage_validated) << "\n";
        j << "    },\n";

        // Channel B
        j << "    {\n";
        j << "      \"channel\": \"B\",\n";
        j << "      \"name\": \"stop_codon_conversion\",\n";
        j << "      \"description\": \"CAA/CAG/CGA to TAA/TAG/TGA via C to T; reference-free damage validator\",\n";
        j << "      \"mechanism\": \"hydrolytic_deamination\",\n";
        // C5: `detected` was `channel_b_valid && !channel_b_inverted`, but channel_b_inverted
        // only becomes true inside the WLS block (terminal total_exposure>1000). When valid is
        // set on interior coverage but the WLS was skipped, detected fired with no fit done.
        // channel_b_quantifiable is set true ONLY when the WLS ran and slope c>0. [behavioral change]
        j << "      \"detected\": " << jbool(dp.channel_b_quantifiable) << ",\n";
        j << "      \"valid\": " << jbool(dp.channel_b_valid) << ",\n";
        j << "      \"lrt_5prime\": "; jfloat(dp.stop_decay_llr_5prime); j << ",\n";
        j << "      \"lrt_5prime_capped\": " << (std::abs(dp.stop_decay_llr_5prime) >= SampleDamageProfile::kLlrCap ? "true" : "false") << ",\n";
        j << "      \"d_max_estimate\": "; jfloat(dp.d_max_from_channel_b); j << ",\n";
        j << "      \"inverted\": " << jbool(dp.channel_b_inverted) << ",\n";
        j << "      \"per_type\": {\n";
        j << "        \"lrt_taa\": "; jfloat(dp.stop_decay_llr_taa_5prime); j << ",\n";
        j << "        \"lrt_taa_capped\": " << (std::abs(dp.stop_decay_llr_taa_5prime) >= SampleDamageProfile::kLlrCap ? "true" : "false") << ",\n";
        j << "        \"lrt_tag\": "; jfloat(dp.stop_decay_llr_tag_5prime); j << ",\n";
        j << "        \"lrt_tag_capped\": " << (std::abs(dp.stop_decay_llr_tag_5prime) >= SampleDamageProfile::kLlrCap ? "true" : "false") << ",\n";
        j << "        \"lrt_tga\": "; jfloat(dp.stop_decay_llr_tga_5prime); j << ",\n";
        j << "        \"lrt_tga_capped\": " << (std::abs(dp.stop_decay_llr_tga_5prime) >= SampleDamageProfile::kLlrCap ? "true" : "false") << ",\n";
        j << "        \"n_caa\": " << static_cast<long long>(dp.n_convertible_caa) << ",\n";
        j << "        \"n_cag\": " << static_cast<long long>(dp.n_convertible_cag) << ",\n";
        j << "        \"n_cga\": " << static_cast<long long>(dp.n_convertible_cga) << ",\n";
        j << "        \"valid_taa\": " << jbool(dp.channel_b_valid_taa) << ",\n";
        j << "        \"valid_tag\": " << jbool(dp.channel_b_valid_tag) << ",\n";
        j << "        \"valid_tga\": " << jbool(dp.channel_b_valid_tga) << "\n";
        j << "      }\n";
        j << "    },\n";

        // Channel C
        j << "    {\n";
        j << "      \"channel\": \"C\",\n";
        j << "      \"name\": \"8_oxog_top_strand\",\n";
        j << "      \"description\": \"G to T oxidative stop codons (GAA/GAG/GGA); top-strand 8-oxoguanine. Scope: top-strand G to T stop codons only; an informative null does not exclude interior oxidation (see oxo_two_marker), and the channel is diluted by modern-read admixture.\",\n";
        j << "      \"mechanism\": \"oxidative_guanine_8_oxog\",\n";
        j << "      \"detected\": " << jbool(dp.ox_damage_detected) << ",\n";
        j << "      \"valid\": " << jbool(dp.channel_c_valid) << ",\n";
        j << "      \"baseline_rate\": "; jfloat(dp.ox_stop_conversion_rate_baseline); j << ",\n";
        j << "      \"terminal_rate\": "; jfloat(dp.ox_stop_rate_terminal); j << ",\n";
        j << "      \"uniformity_ratio\": "; jfloat(dp.ox_uniformity_ratio); j << "\n";
        j << "    },\n";

        // Channel D
        j << "    {\n";
        j << "      \"channel\": \"D\",\n";
        j << "      \"name\": \"chargaff_gt_asymmetry\",\n";
        j << "      \"description\": \"Strand-asymmetric G to T / C to A excess; reference-free 8-oxoG cross-check. SS-only: for DS libraries D cancels by Chargaff second-parity (an expected null, not absence of oxidation). The marker-regressed primary form is oxo_two_marker.\",\n";
        j << "      \"mechanism\": \"oxidative_guanine_8_oxog\",\n";
        // C4-style applicability: D and s_gt cancel for DS libraries by Chargaff's second parity
        // rule (T/(T+G)≈A/(A+C) on both strands), so D is applicable only to SS. Parallel to the
        // F/G `applicable` field, but inverted: applicable = is_ss (true SS, false DS).
        j << "      \"applicable\": " << (is_ss ? "true" : "false") << ",\n";
        j << "      \"detected\": " << jbool(ch_d_detected) << ",\n";
        j << "      \"D\": "; jfloat(dp.ox_gt_asymmetry); j << ",\n";
        j << "      \"s_gt\": "; jfloat(dp.s_gt); j << "\n";
        j << "    },\n";

        // Channel E
        j << "    {\n";
        j << "      \"channel\": \"E\",\n";
        j << "      \"name\": \"purine_enrichment\",\n";
        j << "      \"description\": \"A/G enrichment at 5' read starts; fragmentation-bias / depurination\",\n";
        j << "      \"mechanism\": \"depurination_fragmentation\",\n";
        j << "      \"detected\": " << jbool(dp.depurination_detected) << ",\n";
        j << "      \"enrichment_5prime\": "; jfloat(dp.purine_enrichment_5prime); j << ",\n";
        j << "      \"enrichment_3prime\": "; jfloat(dp.purine_enrichment_3prime); j << "\n";
        j << "    },\n";

        // Layer-3 emitter: F/G/H legend metadata (name/description/mechanism/observable/status/
        // lesion) is sourced from the ChannelSpec registry, so the registry is the single authority
        // for channel semantics; the per-channel numeric fields are emitted inline below.
        auto emit_channel_meta = [&](const ChannelSpec& s) {
            j << "      \"channel\": \"" << s.channel_type << "\",\n";
            j << "      \"name\": \"" << s.json_name << "\",\n";
            j << "      \"description\": \"" << s.json_description << "\",\n";
            j << "      \"mechanism\": \"" << s.json_mechanism << "\",\n";
            j << "      \"observable\": \"" << s.observable_name << "\",\n";
            j << "      \"mechanism_status\": \""
              << (s.mechanism_status == MechanismStatus::ESTABLISHED ? "established" : "empirical") << "\",\n";
            j << "      \"inferred_lesion\": ";
            if (s.inferred_lesion) j << "\"" << s.inferred_lesion << "\""; else j << "null";
            j << ",\n";
            j << "      \"strand\": \""
              << (s.strand == StrandFrame::TOP_5P ? "top_5prime" :
                  s.strand == StrandFrame::TOP_3P ? "top_3prime" : "complement") << "\",\n";
            // Biochem A: the deamination-shadow confound travels IN-BAND with the channel's claim,
            // not only in the count table. True for F (an established bottom-strand 8-oxoG lesion whose
            // z denominator folds in the C->T deamination shadow); false for the empirical channels.
            j << "      \"z_deamination_shadow_in_denominator\": " << (s.shadow_in_z ? "true" : "false") << ",\n";
            j << "      \"inference_class\": \""
              << (s.inference_class == InferenceClass::STOP_CHANNEL ? "stop_channel" : "rate_only") << "\",\n";
        };

        // Channel F
        j << "    {\n";
        emit_channel_meta(stop_channel_spec('F'));
        j << "      \"detected\": " << jbool(ch_f_detected) << ",\n";
        j << "      \"applicable\": " << ((!is_ss || stop_channel_spec('F').applicable_in_ss) ? "true" : "false") << ",\n";
        j << "      \"valid\": " << jbool(dp.channel_f_valid) << ",\n";
        j << "      \"baseline_rate\": "; jfloat(dp.ca_stop_rate_baseline); j << ",\n";
        j << "      \"terminal_rate\": "; jfloat(dp.ca_stop_rate_terminal); j << ",\n";
        j << "      \"uniformity_ratio\": "; jfloat(dp.ca_uniformity_ratio); j << ",\n";
        j << "      \"z_score\": "; emit_z(dp.channel_f_z, dp.channel_f_valid); j << ",\n";
        j << "      \"mh_z\": "; emit_z(dp.channel_f_mh_z, dp.channel_f_valid); j << ",\n";
        j << "      \"z_consistent\": " << jbool(ch_f_z_consistent) << ",\n";
        j << "      \"common_or\": "; jfloat(dp.channel_f_common_or); j << ",\n";
        j << "      \"cf_ratio\": "; jfloat(cf_ratio); j << "\n";
        j << "    },\n";

        // Channel G
        j << "    {\n";
        emit_channel_meta(stop_channel_spec('G'));
        j << "      \"detected\": " << jbool(ch_g_detected) << ",\n";
        j << "      \"applicable\": " << ((!is_ss || stop_channel_spec('G').applicable_in_ss) ? "true" : "false") << ",\n";
        j << "      \"valid\": " << jbool(dp.channel_g_valid) << ",\n";
        j << "      \"baseline_rate\": "; jfloat(dp.cg_stop_rate_baseline); j << ",\n";
        j << "      \"terminal_rate\": "; jfloat(dp.cg_stop_rate_terminal); j << ",\n";
        j << "      \"uniformity_ratio\": "; jfloat(dp.cg_uniformity_ratio); j << ",\n";
        j << "      \"z_score\": "; emit_z(dp.channel_g_z, dp.channel_g_valid); j << ",\n";
        j << "      \"odds_ratio\": " << nan_or(dp.channel_g_or) << "\n";  // P4: 2x2 Haldane-Anscombe OR (primary effect size)
        j << "    },\n";

        // Channel H
        j << "    {\n";
        emit_channel_meta(stop_channel_spec('H'));
        j << "      \"detected\": " << jbool(ch_h_detected) << ",\n";
        j << "      \"valid\": " << jbool(dp.channel_h_valid) << ",\n";
        j << "      \"baseline_rate\": "; jfloat(dp.at_stop_rate_baseline); j << ",\n";
        j << "      \"terminal_rate\": "; jfloat(dp.at_stop_rate_terminal); j << ",\n";
        j << "      \"uniformity_ratio\": "; jfloat(dp.at_uniformity_ratio); j << ",\n";
        j << "      \"z_score\": "; emit_z(dp.channel_h_z, dp.channel_h_valid); j << ",\n";
        j << "      \"z_score_p2plus\": "; emit_z(dp.channel_h_z_p2plus, dp.channel_h_valid); j << ",\n";
        j << "      \"z_consistent\": " << jbool(ch_h_z_consistent) << ",\n";
        j << "      \"odds_ratio\": " << nan_or(dp.channel_h_or) << "\n";  // P4: 2x2 Haldane-Anscombe OR (primary effect size)
        j << "    },\n";

        // Channels F3/G3/H3 — 3' rate-only (registry inference_class=rate_only): terminal/baseline
        // rates at the 3' end, no stop-enrichment z/OR. Declared so "no z here" is a stated contract.
        j << "    {\n";
        emit_channel_meta(kStopChannels3p[0]);
        j << "      \"valid\": " << jbool(dp.channel_f3_valid) << ",\n";
        j << "      \"baseline_rate\": "; jfloat(dp.ca_stop_rate_baseline_3prime); j << ",\n";
        j << "      \"terminal_rate\": "; jfloat(dp.ca_stop_rate_terminal_3prime); j << ",\n";
        j << "      \"uniformity_ratio\": "; jfloat(dp.ca_uniformity_ratio_3prime); j << "\n";
        j << "    },\n";
        j << "    {\n";
        emit_channel_meta(kStopChannels3p[1]);
        j << "      \"valid\": " << jbool(dp.channel_g3_valid) << ",\n";
        j << "      \"baseline_rate\": "; jfloat(dp.cg_stop_rate_baseline_3prime); j << ",\n";
        j << "      \"terminal_rate\": "; jfloat(dp.cg_stop_rate_terminal_3prime); j << ",\n";
        j << "      \"uniformity_ratio\": "; jfloat(dp.cg_uniformity_ratio_3prime); j << "\n";
        j << "    },\n";
        j << "    {\n";
        emit_channel_meta(kStopChannels3p[2]);
        j << "      \"valid\": " << jbool(dp.channel_h3_valid) << ",\n";
        j << "      \"baseline_rate\": "; jfloat(dp.at_stop_rate_baseline_3prime); j << ",\n";
        j << "      \"terminal_rate\": "; jfloat(dp.at_stop_rate_terminal_3prime); j << ",\n";
        j << "      \"uniformity_ratio\": "; jfloat(dp.at_uniformity_ratio_3prime); j << "\n";
        j << "    }\n";

        j << "  ]\n";
    }
    j << "}\n";
}

void count_tables_to_jsonl(const SampleDamageProfile& dp, std::ostream& out) {
    out << std::setprecision(17);
    auto e = [&](double v) { if (std::isfinite(v)) out << v; else out << "null"; };
    auto b = [&](bool v) { out << (v ? "true" : "false"); };
    for (const auto& ct : dp.count_tables) {
        out << "{\"channel_id\":\"" << ct.channel_id << "\""
            << ",\"channel_type\":\"" << ct.channel_type << "\"";
        out << ",\"term_pre\":";       e(ct.term_pre);
        out << ",\"term_stop\":";      e(ct.term_stop);
        out << ",\"int_pre\":";        e(ct.int_pre);
        out << ",\"int_stop\":";       e(ct.int_stop);
        out << ",\"has_shadow\":";     b(ct.has_shadow);
        out << ",\"term_shadow\":";    e(ct.term_shadow);
        out << ",\"int_shadow\":";     e(ct.int_shadow);
        out << ",\"shadow_in_z\":";    b(ct.shadow_in_z);
        out << ",\"shadow_in_rate\":"; b(ct.shadow_in_rate);
        out << ",\"z_num_term\":";     e(ct.z_num_term);
        out << ",\"z_den_term\":";     e(ct.z_den_term);
        out << ",\"z_num_int\":";      e(ct.z_num_int);
        out << ",\"z_den_int\":";      e(ct.z_den_int);
        out << ",\"raw_rate_term\":";  e(ct.raw_rate_term);
        out << ",\"raw_rate_int\":";   e(ct.raw_rate_int);
        out << ",\"pre_clamp_z\":";    e(ct.pre_clamp_z);
        out << ",\"z_cap\":";          e(ct.z_cap);
        out << ",\"cap_applied\":";    b(ct.cap_applied);
        out << ",\"post_clamp_z\":";   e(ct.post_clamp_z);
        if (ct.has_strata) {
            out << ",\"strata\":[";
            for (size_t k = 0; k < ct.strata.size(); ++k) {
                const auto& s = ct.strata[k];
                if (k) out << ",";
                out << "{\"stratum_id\":" << s.stratum_id;
                out << ",\"a_term_stop\":";  e(s.a_term_stop);
                out << ",\"b_term_other\":"; e(s.b_term_other);
                out << ",\"c_int_stop\":";   e(s.c_int_stop);
                out << ",\"d_int_other\":";  e(s.d_int_other);
                out << "}";
            }
            out << "]";
        }
        out << "}\n";
    }
}

} // namespace taph
