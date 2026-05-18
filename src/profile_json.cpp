#include "taph/profile_json.hpp"
#include "taph/library_interpretation.hpp"
#include "taph/length_bin_damage_profile.hpp"
#include <cmath>
#include <iomanip>
#include <ostream>
#include <string>

namespace taph {

static constexpr double kMinCov          = 100.0;
static constexpr float  kOxChannelZDetect = 3.0f;

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
    return std::isnan(static_cast<double>(v)) ? "null" : std::to_string(v);
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
    auto qcf     = compute_library_qc_flags(dp, is_ss,
                       in.flag_hex_artifact,
                       hs.jsd, hs.entropy_terminal,
                       in.short_read_frac);

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
    j << "    \"d_max_combined\": " << dp.d_max_combined << ",\n";
    j << "    \"d_metamatch\": " << dp.d_metamatch << ",\n";
    j << "    \"source\": \"" << dp.d_max_source_str() << "\",\n";
    j << "    \"lambda_5prime\": " << dp.lambda_5prime << ",\n";
    j << "    \"lambda_3prime\": " << dp.lambda_3prime << ",\n";
    j << "    \"bg_5prime\": " << dp.fit_baseline_5prime << ",\n";
    j << "    \"bg_3prime\": " << dp.fit_baseline_3prime << ",\n";
    j << "    \"validated\": " << (dp.damage_validated ? "true" : "false") << ",\n";
    j << "    \"artifact\": " << (dp.damage_artifact ? "true" : "false") << ",\n";
    j << "    \"mixture_gc\": {\n";
    j << "      \"d_ancient\": " << dp.mixture_d_ancient << ",\n";
    j << "      \"pi_ancient\": " << dp.mixture_pi_ancient << ",\n";
    j << "      \"d_reference\": " << dp.mixture_d_reference << ",\n";
    j << "      \"identifiable\": " << (dp.mixture_identifiable ? "true" : "false") << ",\n";
    j << "      \"converged\": " << (dp.mixture_converged ? "true" : "false") << ",\n";
    j << "      \"n_components\": " << dp.mixture_n_components << "\n";
    j << "    },\n";
    j << "    \"cpg_like\": {\n";
    j << "      \"dmax_ct5_cpg\": "           << nan_or(dp.dmax_ct5_cpg_like)    << ",\n";
    j << "      \"dmax_ct5_noncpg\": "        << nan_or(dp.dmax_ct5_noncpg_like) << ",\n";
    j << "      \"cpg_ratio\": "              << nan_or(dp.cpg_ratio)            << ",\n";
    j << "      \"log2_cpg_ratio\": "         << nan_or(dp.log2_cpg_ratio)       << ",\n";
    j << "      \"baseline_cpg\": "           << nan_or(dp.fit_baseline_ct5_cpg_like)    << ",\n";
    j << "      \"baseline_noncpg\": "        << nan_or(dp.fit_baseline_ct5_noncpg_like) << ",\n";
    j << "      \"cov_terminal_cpg\": "       << std::setprecision(0) << dp.cov_ct5_cpg_like_terminal    << ",\n";
    j << "      \"cov_terminal_noncpg\": "    << dp.cov_ct5_noncpg_like_terminal << ",\n";
    j << "      \"effcov_terminal_cpg\": "    << dp.effcov_ct5_cpg_like_terminal    << ",\n";
    j << "      \"effcov_terminal_noncpg\": " << dp.effcov_ct5_noncpg_like_terminal << ",\n";
    j << "      \"cpg_score_z\": " << std::setprecision(6) << cpg.z << ",\n";
    j << "      \"cpg_score_p\": " << cpg.p << "\n";
    j << "    },\n";
    j << "    \"context_deamination\": {\n";
    j << "      \"dmax_AC\": " << nan_or(dp.dmax_ct5_by_upstream[SP::CTX_AC]) << ",\n";
    j << "      \"dmax_CC\": " << nan_or(dp.dmax_ct5_by_upstream[SP::CTX_CC]) << ",\n";
    j << "      \"dmax_GC\": " << nan_or(dp.dmax_ct5_by_upstream[SP::CTX_GC]) << ",\n";
    j << "      \"dmax_TC\": " << nan_or(dp.dmax_ct5_by_upstream[SP::CTX_TC]) << ",\n";
    j << "      \"dipyr_contrast\": " << nan_or(dp.dipyr_contrast) << ",\n";
    j << "      \"cpg_contrast\": " << nan_or(dp.cpg_contrast) << ",\n";
    j << "      \"heterogeneity_chi2\": " << std::setprecision(2) << dp.context_heterogeneity_chi2 << ",\n";
    j << "      \"heterogeneity_p\": " << std::setprecision(4) << dp.context_heterogeneity_p << ",\n";
    j << "      \"heterogeneity_detected\": " << (dp.context_heterogeneity_detected ? "true" : "false") << "\n";
    j << "    },\n";
    j << "    \"per_pos_5prime_ct\": [";
    for (int p = 0; p < N_POS; ++p) {
        double v = (dp.tc_total_5prime[p] >= kMinCov) ? dp.t_freq_5prime[p] : -1.0;
        j << std::setprecision(6) << v;
        if (p < N_POS - 1) j << ",";
    }
    j << "],\n";
    j << "    \"per_pos_3prime\": [";
    for (int p = 0; p < N_POS; ++p) {
        double v;
        if (is_ss)
            v = (dp.tc_total_3prime[p] >= kMinCov)
                    ? dp.t_freq_3prime[p] / dp.tc_total_3prime[p] : -1.0;
        else
            v = (dp.ag_total_3prime[p] >= kMinCov) ? dp.a_freq_3prime[p] : -1.0;
        j << std::setprecision(6) << v;
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
              << ",\"cpg_contrast\":"
              << (std::isnan(lb.cpg_contrast) ? std::string("null") : std::to_string(lb.cpg_contrast))
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
            j << ",\"mixture\":{\"d_ancient\":" << std::setprecision(6) << lb.mixture_d_damaged
              << ",\"d_reference\":" << lb.mixture_d_reference
              << ",\"d_population\":" << lb.mixture_d_population
              << ",\"pi_ancient\":" << lb.mixture_pi_damaged
              << ",\"n_components\":" << lb.mixture_n_components
              << ",\"converged\":" << (lb.mixture_converged ? "true" : "false")
              << ",\"identifiable\":" << (lb.mixture_identifiable ? "true" : "false")
              << "}";
            j << ",\"gc_d_max\":[";
            for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
                j << std::setprecision(6) << lb.gc_d_max[g];
                if (g + 1 < LengthBinDamageProfile::N_GC_BINS) j << ",";
            }
            j << "],\"gc_n_reads\":[";
            for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
                j << lb.gc_n_reads[g];
                if (g + 1 < LengthBinDamageProfile::N_GC_BINS) j << ",";
            }
            j << "],\"gc_p_damaged\":[";
            for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
                j << std::setprecision(6) << lb.gc_p_damaged[g];
                if (g + 1 < LengthBinDamageProfile::N_GC_BINS) j << ",";
            }
            j << "]";
            auto dnan = [](double v) -> std::string {
                return std::isnan(v) ? "null" : std::to_string(v);
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
                j << std::setprecision(6) << row[g];
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
    j << "    \"gt_bg_interior_mean\": " << dp.g_bg_interior_mean << ",\n";
    j << "    \"s_gt\": " << dp.s_gt << ",\n";
    j << "    \"per_pos_5prime_gt\": [";
    for (int p = 0; p < N_POS; ++p) {
        double denom = dp.t_from_g_5prime[p] + dp.g_count_5prime[p];
        double v = (denom >= kMinCov) ? dp.t_from_g_5prime[p] / denom : -1.0;
        j << std::setprecision(6) << v;
        if (p < N_POS - 1) j << ",";
    }
    j << "],\n";
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
    j << "    \"channel_f_z\": " << dp.channel_f_z << ",\n";
    j << "    \"channel_f_mh_z\": " << dp.channel_f_mh_z << ",\n";
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
    j << "    \"channel_g_z\": " << dp.channel_g_z << ",\n";
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
    j << "    \"channel_h_z\": " << dp.channel_h_z << ",\n";
    j << "    \"channel_h_z_p2plus\": " << dp.channel_h_z_p2plus << ",\n";
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
    j << "    \"s_oxog\": " << std::setprecision(6)
      << (in.has_oxog_score ? in.s_oxog : 0.0) << ",\n";
    j << "    \"se_s_oxog\": " << (in.has_oxog_score ? in.se_s_oxog : 0.0) << ",\n";
    j << "    \"ox_d_oriented\": " << std::setprecision(6) << in.d_oriented << ",\n";
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
        double ctx_wsum = 0.0, ctx_w2 = 0.0;
        int ctx_pos = 0;
        for (int i = 0; i < 16; ++i) {
            float v = dp.s_oxog_16ctx[i];
            double cov = static_cast<double>(dp.cov_oxog_16ctx[i]);
            if (std::isnan(v) || cov <= 0) continue;
            double w = std::sqrt(cov);
            ctx_wsum += v * w;
            ctx_w2   += cov;
            if (v > 0) ++ctx_pos;
        }
        double ctx_z = (ctx_w2 > 0) ? ctx_wsum / std::sqrt(ctx_w2) : 0.0;
        j << "    \"oxog_ctx_n_positive\": " << ctx_pos << ",\n";
        j << "    \"oxog_ctx_z\": " << std::setprecision(4) << ctx_z << ",\n";
        j << "    \"oxog_score_z\": " << std::setprecision(6) << oxog_is.z << ",\n";
        j << "    \"oxog_score_p\": " << std::setprecision(6) << oxog_is.p << ",\n";
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
              << ",\"detected\":" << (lb.channel_h_valid && (lb.channel_h_z > kOxChannelZDetect || lb.channel_h_z_p2plus > kOxChannelZDetect) ? "true" : "false") << "}"
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
    j << "    \"valid\": " << (otm.valid ? "true" : "false") << ",\n";
    j << "    \"n_cells_used\": " << otm.n_cells_used << ",\n";
    j << "    \"beta1\": " << std::setprecision(6) << otm.beta1 << ",\n";
    j << "    \"beta1_se\": " << otm.beta1_se << ",\n";
    j << "    \"beta1_z\": " << otm.beta1_z << ",\n";
    j << "    \"beta2\": " << otm.beta2 << ",\n";
    j << "    \"beta2_se\": " << otm.beta2_se << ",\n";
    j << "    \"beta2_z\": " << otm.beta2_z << ",\n";
    j << "    \"alpha\": " << otm.alpha << ",\n";
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
        double O = dp.interior_ct_cluster_short_obs;
        double E = dp.interior_ct_cluster_short_exp;
        double sz = 0.0, sp = 1.0;
        if (O > 0 && E > 1e-9) {
            double sign_val = (O >= E) ? 1.0 : -1.0;
            sz = sign_val * std::sqrt(2.0 * (O * std::log(O / E) - (O - E)));
            sp = 0.5 * std::erfc(sz / std::sqrt(2.0));
        }
        j << "    \"short_score_z\": " << std::setprecision(6) << sz << ",\n";
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
        j << "    \"flags\": [";
        bool first_flag = true;
        auto emit_flag = [&](const char* name) {
            if (!first_flag) j << ",";
            j << "\"" << name << "\"";
            first_flag = false;
        };
        if (qcf.adapter_remnant_5prime)   emit_flag("adapter_remnant_5prime");
        if (qcf.adapter_remnant_3prime)   emit_flag("adapter_remnant_3prime");
        if (qcf.hexamer_composition_bias) emit_flag("hexamer_composition_bias");
        if (qcf.hexamer_terminal_shift)   emit_flag("hexamer_terminal_shift");
        if (qcf.short_read_spike)         emit_flag("short_read_spike");
        if (qcf.depurination)             emit_flag("depurination");
        if (qcf.ds_3prime_signal_absent)  emit_flag("ds_3prime_signal_absent");
        if (qcf.ga3_inward_displaced)     emit_flag("ga3_inward_displaced");
        if (qcf.hexamer_artifact_bias)    emit_flag("hexamer_artifact_bias");
        j << "]\n";
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
    j << "    \"library_submodel_bic\": {\n";
    j << "      \"M_bias\": "        << std::setprecision(2) << dp.library_bic_M_bias        << ",\n";
    j << "      \"M_DS_symm\": "     << dp.library_bic_M_DS_symm     << ",\n";
    j << "      \"M_DS_spike\": "    << dp.library_bic_M_DS_spike    << ",\n";
    j << "      \"M_DS_symm_art\": " << dp.library_bic_M_DS_symm_art << ",\n";
    j << "      \"M_SS_comp\": "     << dp.library_bic_M_SS_comp     << ",\n";
    j << "      \"M_SS_orig\": "     << dp.library_bic_M_SS_orig     << ",\n";
    j << "      \"M_SS_asym\": "     << dp.library_bic_M_SS_asym     << ",\n";
    j << "      \"M_SS_full\": "     << dp.library_bic_M_SS_full     << ",\n";
    j << "      \"M_DS_asym_art\": " << dp.library_bic_M_DS_asym_art << "\n";
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

    // ── Damage types ──────────────────────────────────────────────────────────
    {
        auto jbool = [](bool v) -> const char* { return v ? "true" : "false"; };
        auto jfloat = [&](double v) { j << std::setprecision(6) << v; };

        double cf_ratio = (dp.channel_c_valid && dp.channel_f_valid
                           && dp.ca_stop_rate_baseline > 1e-6f)
            ? static_cast<double>(dp.ox_stop_conversion_rate_baseline)
              / static_cast<double>(dp.ca_stop_rate_baseline)
            : -1.0;

        bool ch_f_detected = dp.channel_f_valid && dp.channel_f_z > kOxChannelZDetect;
        bool ch_g_detected = dp.channel_g_valid && dp.channel_g_z > kOxChannelZDetect;
        bool ch_h_detected = dp.channel_h_valid &&
            (dp.channel_h_z > kOxChannelZDetect || dp.channel_h_z_p2plus > kOxChannelZDetect);
        bool ch_d_detected = std::abs(dp.ox_gt_asymmetry) > 0.01f;

        j << "  \"damage_types\": [\n";

        // Channel A
        j << "    {\n";
        j << "      \"channel\": \"A\",\n";
        j << "      \"name\": \"cytosine_deamination\",\n";
        j << "      \"description\": \"C to T post-mortem deamination at terminal positions\",\n";
        j << "      \"mechanism\": \"hydrolytic_deamination\",\n";
        j << "      \"detected\": " << jbool(dp.damage_validated) << ",\n";
        j << "      \"d_max_5prime\": "; jfloat(dp.d_max_5prime); j << ",\n";
        j << "      \"d_max_3prime\": "; jfloat(dp.d_max_3prime); j << ",\n";
        j << "      \"validated\": " << jbool(dp.damage_validated) << "\n";
        j << "    },\n";

        // Channel B
        j << "    {\n";
        j << "      \"channel\": \"B\",\n";
        j << "      \"name\": \"stop_codon_conversion\",\n";
        j << "      \"description\": \"CAA/CAG/CGA to TAA/TAG/TGA via C to T; reference-free damage validator\",\n";
        j << "      \"mechanism\": \"hydrolytic_deamination\",\n";
        j << "      \"detected\": " << jbool(dp.channel_b_valid && !dp.channel_b_inverted) << ",\n";
        j << "      \"valid\": " << jbool(dp.channel_b_valid) << ",\n";
        j << "      \"lrt_5prime\": "; jfloat(dp.stop_decay_llr_5prime); j << ",\n";
        j << "      \"d_max_estimate\": "; jfloat(dp.d_max_from_channel_b); j << ",\n";
        j << "      \"inverted\": " << jbool(dp.channel_b_inverted) << ",\n";
        j << "      \"per_type\": {\n";
        j << "        \"lrt_taa\": "; jfloat(dp.stop_decay_llr_taa_5prime); j << ",\n";
        j << "        \"lrt_tag\": "; jfloat(dp.stop_decay_llr_tag_5prime); j << ",\n";
        j << "        \"lrt_tga\": "; jfloat(dp.stop_decay_llr_tga_5prime); j << ",\n";
        j << "        \"n_caa\": " << static_cast<long long>(dp.n_convertible_caa) << ",\n";
        j << "        \"n_cag\": " << static_cast<long long>(dp.n_convertible_cag) << ",\n";
        j << "        \"n_cga\": " << static_cast<long long>(dp.n_convertible_cga) << ",\n";
        j << "        \"valid_tga\": " << jbool(dp.channel_b_valid_tga) << "\n";
        j << "      }\n";
        j << "    },\n";

        // Channel C
        j << "    {\n";
        j << "      \"channel\": \"C\",\n";
        j << "      \"name\": \"8_oxog_top_strand\",\n";
        j << "      \"description\": \"G to T oxidative stop codons (GAA/GAG/GGA); top-strand 8-oxoguanine\",\n";
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
        j << "      \"description\": \"Strand-asymmetric G to T / C to A excess; reference-free 8-oxoG cross-check\",\n";
        j << "      \"mechanism\": \"oxidative_guanine_8_oxog\",\n";
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

        // Channel F
        j << "    {\n";
        j << "      \"channel\": \"F\",\n";
        j << "      \"name\": \"8_oxog_complement\",\n";
        j << "      \"description\": \"C to A oxidative stop codons (TCA/TCG/TAC/TGC); bottom-strand 8-oxoguanine\",\n";
        j << "      \"mechanism\": \"oxidative_guanine_8_oxog\",\n";
        j << "      \"detected\": " << jbool(ch_f_detected) << ",\n";
        j << "      \"applicable\": " << (!is_ss ? "true" : "false") << ",\n";
        j << "      \"valid\": " << jbool(dp.channel_f_valid) << ",\n";
        j << "      \"baseline_rate\": "; jfloat(dp.ca_stop_rate_baseline); j << ",\n";
        j << "      \"terminal_rate\": "; jfloat(dp.ca_stop_rate_terminal); j << ",\n";
        j << "      \"uniformity_ratio\": "; jfloat(dp.ca_uniformity_ratio); j << ",\n";
        j << "      \"z_score\": "; jfloat(dp.channel_f_z); j << ",\n";
        j << "      \"mh_z\": "; jfloat(dp.channel_f_mh_z); j << ",\n";
        j << "      \"common_or\": "; jfloat(dp.channel_f_common_or); j << ",\n";
        j << "      \"cf_ratio\": "; jfloat(cf_ratio); j << "\n";
        j << "    },\n";

        // Channel G
        j << "    {\n";
        j << "      \"channel\": \"G\",\n";
        j << "      \"name\": \"hydantoin_oxidation\",\n";
        j << "      \"description\": \"C to G stop codons (TCA/TAC to TGA/TAG); spiroiminodihydantoin / guanidinohydantoin\",\n";
        j << "      \"mechanism\": \"oxidative_guanine_hydantoin\",\n";
        j << "      \"detected\": " << jbool(ch_g_detected) << ",\n";
        j << "      \"applicable\": " << (!is_ss ? "true" : "false") << ",\n";
        j << "      \"valid\": " << jbool(dp.channel_g_valid) << ",\n";
        j << "      \"baseline_rate\": "; jfloat(dp.cg_stop_rate_baseline); j << ",\n";
        j << "      \"terminal_rate\": "; jfloat(dp.cg_stop_rate_terminal); j << ",\n";
        j << "      \"uniformity_ratio\": "; jfloat(dp.cg_uniformity_ratio); j << ",\n";
        j << "      \"z_score\": "; jfloat(dp.channel_g_z); j << "\n";
        j << "    },\n";

        // Channel H
        j << "    {\n";
        j << "      \"channel\": \"H\",\n";
        j << "      \"name\": \"adenine_oxidation\",\n";
        j << "      \"description\": \"A to T stop codons (AAA/AAG/AGA to TAA/TAG/TGA); adenine oxidation or trans-lesion\",\n";
        j << "      \"mechanism\": \"oxidative_adenine\",\n";
        j << "      \"detected\": " << jbool(ch_h_detected) << ",\n";
        j << "      \"valid\": " << jbool(dp.channel_h_valid) << ",\n";
        j << "      \"baseline_rate\": "; jfloat(dp.at_stop_rate_baseline); j << ",\n";
        j << "      \"terminal_rate\": "; jfloat(dp.at_stop_rate_terminal); j << ",\n";
        j << "      \"uniformity_ratio\": "; jfloat(dp.at_uniformity_ratio); j << ",\n";
        j << "      \"z_score\": "; jfloat(dp.channel_h_z); j << ",\n";
        j << "      \"z_score_p2plus\": "; jfloat(dp.channel_h_z_p2plus); j << "\n";
        j << "    }\n";

        j << "  ]\n";
    }
    j << "}\n";
}

} // namespace taph
