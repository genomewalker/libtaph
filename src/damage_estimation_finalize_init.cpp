#include "damage_estimation_finalize_helpers.hpp"
#include "damage_estimation_finalize_ctx.hpp"
namespace taph {

void finalize_init(SampleDamageProfile& profile, FinalCtx& ctx) {
    if (profile.n_reads == 0) return;
    // Lifecycle guard (see SampleDamageProfile::finalized): finalize mutates
    // raw counts into rates in place — calling it twice would normalize the
    // already-normalized arrays again. Throw rather than silently corrupt.
    if (profile.finalized) {
        throw std::logic_error(
            "finalize_sample_profile: profile already finalized; "
            "reset_sample_profile() must be called before re-running update/finalize.");
    }

    // Empirical GG-breakpoint contrast (GG-vs-A terminal double-difference; see header doc). The
    // trinucleotide spectra are raw counts and are not modified by finalize, so compute here.
    profile.oxidative_scission_delta_5prime =
        compute_oxscission_delta(profile.tri_5prime_terminal, profile.tri_5prime_interior);
    profile.oxidative_scission_delta_3prime =
        compute_oxscission_delta(profile.tri_3prime_terminal, profile.tri_3prime_interior);
    // PRIMARY = 5' only: the 3' terminus G is depleted by G->A deamination, contaminating the 3'
    // GG-G count, so the 3' channel mixes the endpoint-context contrast with deamination/preservation.
    // This is an empirical proxy, not direct oxidative lesion identification.
    profile.oxidative_scission_delta = profile.oxidative_scission_delta_5prime;

    // Capture raw counts before normalization for statistical tests
    double base_tc_total = profile.baseline_t_freq + profile.baseline_c_freq;
    double base_ag_total = profile.baseline_a_freq + profile.baseline_g_freq;

    auto compute_terminal_stats = [](double term_t, double term_c, double base_t, double base_c) {
        double term_total = term_t + term_c;
        double base_total = base_t + base_c;
        if (term_total < 1.0 || base_total < 1.0) {
            return std::pair<float, float>{0.0f, 0.0f};
        }
        double p_term = term_t / term_total;
        double p_base = base_t / base_total;
        double p_pool = (term_t + base_t) / (term_total + base_total);
        double se = std::sqrt(std::max(1e-9, p_pool * (1.0 - p_pool) * (1.0 / term_total + 1.0 / base_total)));
        double z = se > 0.0 ? (p_term - p_base) / se : 0.0;
        return std::pair<float, float>{static_cast<float>(p_term - p_base), static_cast<float>(z)};
    };

    // Terminal enrichment/depletion statistics
    // Check for position-0 artifact: if pos0 shows depletion but pos1 shows enrichment,
    // use pos1 instead (common adapter ligation bias pattern)
    auto stats_5_pos0 = compute_terminal_stats(profile.t_freq_5prime[0], profile.c_freq_5prime[0],
                                               profile.baseline_t_freq, profile.baseline_c_freq);
    auto stats_5_pos1 = compute_terminal_stats(profile.t_freq_5prime[1], profile.c_freq_5prime[1],
                                               profile.baseline_t_freq, profile.baseline_c_freq);

    // Detect position-0 artifact: two criteria
    // 1. Classic: pos0 depleted but pos1 enriched (relative to baseline)
    // 2. Jump pattern: Large pos0-to-pos1 difference (>3%) with pos1 enriched (>2%)
    //    This catches cases where baseline elevation masks the pos0 depletion
    bool pos0_artifact_5 = false;
    if (stats_5_pos0.first < -0.005f && stats_5_pos1.first > 0.01f) {
        // Classic pattern: pos0 below baseline, pos1 above
        pos0_artifact_5 = true;
    } else if (stats_5_pos1.first > 0.02f && (stats_5_pos1.first - stats_5_pos0.first) > 0.03f) {
        // Large pos0-to-pos1 jump pattern: pos1 clearly elevated but pos0 isn't
        // This indicates adapter artifact masking damage at pos0
        pos0_artifact_5 = true;
    }

    if (pos0_artifact_5) {
        profile.terminal_shift_5prime = stats_5_pos1.first;
        profile.terminal_z_5prime = stats_5_pos1.second;
        profile.position_0_artifact_5prime = true;
    } else {
        profile.terminal_shift_5prime = stats_5_pos0.first;
        profile.terminal_z_5prime = stats_5_pos0.second;
        profile.position_0_artifact_5prime = false;
    }

    // Always use pooled positions 1-4 for the 3' terminal G→A signal.
    // Position 0 is excluded: SS library prep introduces an adapter ligation artifact
    // at 3' pos-0 (elevated A/(A+G)) that is unrelated to aDNA damage. At high coverage
    // this artifact drives terminal_z_3prime into the hundreds and falsely blocks SS
    // auto-detection. Report pos-0 separately as position_0_artifact_3prime.
    double sum_a_3_1_4 = 0.0, sum_g_3_1_4 = 0.0;
    for (int p = 1; p <= 4; ++p) {
        sum_a_3_1_4 += profile.a_freq_3prime[p];
        sum_g_3_1_4 += profile.g_freq_3prime[p];
    }
    auto stats_3_pos1_4 = compute_terminal_stats(sum_a_3_1_4, sum_g_3_1_4,
                                                 4.0 * profile.baseline_a_freq,
                                                 4.0 * profile.baseline_g_freq);
    profile.terminal_shift_3prime = stats_3_pos1_4.first;
    profile.terminal_z_3prime     = stats_3_pos1_4.second;

    // Flag pos-0 as an artifact when it is:
    //   (a) substantially elevated above pos1-4 mean — SS ligation artifact spike, or
    //   (b) substantially depleted below pos1-4 mean — DS adapter blunting gap
    //       (mirrors position_0_artifact_5prime depleted-pos0 / enriched-pos1 logic)
    {
        double tot0   = profile.a_freq_3prime[0] + profile.g_freq_3prime[0];
        double tot1_4 = sum_a_3_1_4 + sum_g_3_1_4;
        double ga0    = tot0   > 0.0 ? profile.a_freq_3prime[0] / tot0   : 0.0;
        double ga1_4  = tot1_4 > 0.0 ? sum_a_3_1_4              / tot1_4 : 0.0;
        // Use pos1 alone (mirrors stats_5_pos1) — pos1-4 average diluted by interior positions
        auto stats_3_pos0 = compute_terminal_stats(profile.a_freq_3prime[0],
                                                   profile.g_freq_3prime[0],
                                                   profile.baseline_a_freq,
                                                   profile.baseline_g_freq);
        auto stats_3_pos1 = compute_terminal_stats(profile.a_freq_3prime[1],
                                                   profile.g_freq_3prime[1],
                                                   profile.baseline_a_freq,
                                                   profile.baseline_g_freq);
        bool pos0_spike_3 = (ga0 - ga1_4) > 0.05;
        bool pos0_gap_3   = (stats_3_pos0.first < -0.005f && stats_3_pos1.first > 0.01f) ||
                            (stats_3_pos1.first > 0.02f &&
                             (stats_3_pos1.first - stats_3_pos0.first) > 0.03f);
        profile.position_0_artifact_3prime = pos0_spike_3 || pos0_gap_3;
    }

    // Negative control statistics (should NOT show enrichment if damage is real)
    // 5' control: A/(A+G) at 5' end vs interior
    auto ctrl_5 = compute_terminal_stats(profile.a_freq_5prime[0], profile.g_freq_5prime[0],
                                         profile.baseline_a_freq, profile.baseline_g_freq);
    profile.ctrl_shift_5prime = ctrl_5.first;
    profile.ctrl_z_5prime = ctrl_5.second;

    // 3' control: T/(T+C) at 3' end vs interior
    auto ctrl_3 = compute_terminal_stats(profile.t_freq_3prime[0], profile.c_freq_3prime[0],
                                         profile.baseline_t_freq, profile.baseline_c_freq);
    profile.ctrl_shift_3prime = ctrl_3.first;
    profile.ctrl_z_3prime = ctrl_3.second;

    // Channel divergence: difference between damage shift and control shift
    // High divergence = real damage (only damage channel elevated)
    // Low divergence = composition bias (both channels elevated together)
    profile.channel_divergence_5prime = std::abs(profile.terminal_shift_5prime - profile.ctrl_shift_5prime);
    profile.channel_divergence_3prime = std::abs(profile.terminal_shift_3prime - profile.ctrl_shift_3prime);

    // Channel E: terminal purine enrichment for AP-site/depurination proxy.
    // This uses A+G over all A/C/G/T bases. It must be computed before the
    // per-position arrays are normalized below.
    struct PurineStats {
        float terminal_rate = std::numeric_limits<float>::quiet_NaN();
        float baseline_rate = std::numeric_limits<float>::quiet_NaN();
        float enrichment = std::numeric_limits<float>::quiet_NaN();
        float z = std::numeric_limits<float>::quiet_NaN();
        bool valid = false;
    };
    auto compute_purine_stats = [](double term_a, double term_g, double term_c, double term_t,
                                   double base_a, double base_g, double base_c, double base_t) {
        PurineStats r;
        const double term_pur = term_a + term_g;
        const double base_pur = base_a + base_g;
        const double term_total = term_pur + term_c + term_t;
        const double base_total = base_pur + base_c + base_t;
        if (term_total < 100.0 || base_total < 100.0) return r;
        const double p_term = term_pur / term_total;
        const double p_base = base_pur / base_total;
        if (p_base <= 0.01 || p_base >= 0.99) return r;
        const double p_pool = (term_pur + base_pur) / (term_total + base_total);
        const double se = std::sqrt(std::max(
            1e-12,
            p_pool * (1.0 - p_pool) * (1.0 / term_total + 1.0 / base_total)));
        r.terminal_rate = static_cast<float>(p_term);
        r.baseline_rate = static_cast<float>(p_base);
        r.enrichment = static_cast<float>(p_term - p_base);
        r.z = static_cast<float>((p_term - p_base) / se);
        r.valid = true;
        return r;
    };
    const double base_a = profile.baseline_a_freq;
    const double base_g = profile.baseline_g_freq;
    const double base_c = profile.baseline_c_freq;
    const double base_t = profile.baseline_t_freq;
    const PurineStats pur5 = compute_purine_stats(
        profile.a_freq_5prime[0], profile.g_freq_5prime[0],
        profile.c_freq_5prime[0], profile.t_freq_5prime[0],
        base_a, base_g, base_c, base_t);
    const PurineStats pur3 = compute_purine_stats(
        profile.a_freq_3prime[0], profile.g_freq_3prime[0],
        profile.c_freq_3prime[0], profile.t_freq_3prime[0],
        base_a, base_g, base_c, base_t);
    profile.channel_e_valid = pur5.valid;
    profile.purine_rate_terminal_5prime = pur5.terminal_rate;
    profile.purine_rate_terminal_3prime = pur3.terminal_rate;
    profile.purine_rate_interior = pur5.baseline_rate;
    profile.purine_enrichment_5prime = pur5.enrichment;
    profile.purine_enrichment_3prime = pur3.enrichment;
    profile.purine_z_5prime = pur5.z;
    profile.purine_z_3prime = pur3.z;
    profile.depurination_detected =
        pur5.valid && pur5.enrichment > 0.02f && pur5.z > 3.0f &&
        !profile.position_0_artifact_5prime;

    // Normalize baseline frequencies (using double to preserve precision)
    double mid_total = profile.baseline_t_freq + profile.baseline_c_freq +
                       profile.baseline_a_freq + profile.baseline_g_freq;
    if (mid_total > 0) {
        profile.baseline_t_freq /= mid_total;
        profile.baseline_c_freq /= mid_total;
        profile.baseline_a_freq /= mid_total;
        profile.baseline_g_freq /= mid_total;
    }

    // Compute baseline T/C and A/G ratios from the middle of reads (for inverted pattern detection)
    ctx.baseline_tc = profile.baseline_t_freq /
                        (profile.baseline_t_freq + profile.baseline_c_freq + 0.001);
    ctx.baseline_ag = profile.baseline_a_freq /
                        (profile.baseline_a_freq + profile.baseline_g_freq + 0.001);

    {
        JointDamageSuffStats jstats;

        // BEHAVIORAL CHANGE (C5): mask 5' position 0 from the joint model when an
        // adapter-remnant composition artifact is flagged there. The artifact
        // depletes T at pos0 / enriches pos1, producing a non-monotone profile
        // that defeats the exponential-decay model (M1) and zeroes p_damage even
        // on heavily damaged libraries. Mirrors the p0_tc_5 mask already applied
        // to channels F/G/H. Zeroing the suff-stat counts excludes pos0 from both
        // the M0 and M1 likelihoods without shifting any other position's data.
        const int p0_jt = profile.position_0_artifact_5prime ? 1 : 0;

        // Channel A: T/(T+C) at 5' end (raw counts still in t_freq_5prime, c_freq_5prime)
        for (int p = 0; p < 15; ++p) {
            if (p < p0_jt) { jstats.k_tc[p] = 0; jstats.n_tc[p] = 0; continue; }
            jstats.k_tc[p] = static_cast<uint64_t>(profile.t_freq_5prime[p]);
            jstats.n_tc[p] = static_cast<uint64_t>(profile.t_freq_5prime[p] + profile.c_freq_5prime[p]);
        }

        // Control: A/(A+G) at 5' end
        for (int p = 0; p < 15; ++p) {
            if (p < p0_jt) { jstats.k_ag[p] = 0; jstats.n_ag[p] = 0; continue; }
            jstats.k_ag[p] = static_cast<uint64_t>(profile.a_freq_5prime[p]);
            jstats.n_ag[p] = static_cast<uint64_t>(profile.a_freq_5prime[p] + profile.g_freq_5prime[p]);
        }

        // Channel B: stop codon conversions (TAA+TAG+TGA) / (CAA+CAG+CGA + TAA+TAG+TGA)
        for (int p = 0; p < 15; ++p) {
            if (p < p0_jt) { jstats.k_stop[p] = 0; jstats.n_stop[p] = 0; continue; }
            double stops = profile.convertible_taa_5prime[p] +
                          profile.convertible_tag_5prime[p] +
                          profile.convertible_tga_5prime[p];
            double pre = profile.convertible_caa_5prime[p] +
                        profile.convertible_cag_5prime[p] +
                        profile.convertible_cga_5prime[p];
            jstats.k_stop[p] = static_cast<uint64_t>(stops);
            jstats.n_stop[p] = static_cast<uint64_t>(pre + stops);
        }

        // Interior baselines (from middle of reads)
        jstats.k_tc_interior = static_cast<uint64_t>(base_tc_total * ctx.baseline_tc);
        jstats.n_tc_interior = static_cast<uint64_t>(base_tc_total);
        jstats.k_ag_interior = static_cast<uint64_t>(base_ag_total * ctx.baseline_ag);
        jstats.n_ag_interior = static_cast<uint64_t>(base_ag_total);

        // Channel B interior
        double stop_int = profile.convertible_taa_interior +
                         profile.convertible_tag_interior +
                         profile.convertible_tga_interior;
        double pre_int = profile.convertible_caa_interior +
                        profile.convertible_cag_interior +
                        profile.convertible_cga_interior;
        jstats.k_stop_interior = static_cast<uint64_t>(stop_int);
        jstats.n_stop_interior = static_cast<uint64_t>(pre_int + stop_int);

        // Fit the joint model
        JointDamageResult jresult = JointDamageModel::fit(jstats);

        // Store results in profile
        profile.joint_delta_max = jresult.delta_max;
        profile.joint_lambda = jresult.lambda;
        profile.joint_a_max = jresult.a_max;
        profile.joint_log_lik_m1 = jresult.log_lik_m1;
        profile.joint_log_lik_m0 = jresult.log_lik_m0;
        profile.joint_delta_bic = jresult.delta_bic;
        profile.joint_bayes_factor = jresult.bayes_factor;
        profile.joint_p_damage = jresult.p_damage;
        profile.joint_z_delta = jresult.z_delta;
        profile.joint_z_delta_capped = jresult.z_delta_capped;
        profile.joint_n_informative = jresult.n_informative;
        profile.joint_model_valid = jresult.valid;

        // Use joint model for damage decision
        if (jresult.valid) {
            // Check for strong evidence of damage
            // P(damage) > 0.95 OR ΔBIC > 10 (handles overflow when BF=inf)
            bool strong_damage_evidence = jresult.p_damage > 0.95f ||
                                          jresult.delta_bic > 10.0f ||
                                          std::isinf(jresult.bayes_factor);

            if (strong_damage_evidence) {
                profile.damage_validated = true;
                profile.damage_artifact = false;
            } else if (jresult.p_damage < 0.05f && jresult.a_max > 0.02f) {
                // Strong evidence for artifact (composition bias detected)
                profile.damage_validated = false;
                profile.damage_artifact = true;
            } else if (jresult.delta_bic < -10.0f) {
                // Strong evidence against damage (M0 fits much better)
                profile.damage_validated = false;
                profile.damage_artifact = (jresult.a_max > 0.02f);
            } else {
                // Uncertain - could go either way
                profile.damage_validated = false;
                profile.damage_artifact = false;
            }
        }
    }

    // Snapshot raw ACGT terminal counts BEFORE the in-place normalization below
    // (the artifact_g_excess ->G-overcall flag needs G over the full ACGT total).
    for (int i = 0; i < 15; ++i) {
        profile.a_freq_5prime_raw[i] = profile.a_freq_5prime[i];
        profile.g_freq_5prime_raw[i] = profile.g_freq_5prime[i];
        profile.t_freq_5prime_raw[i] = profile.t_freq_5prime[i];
        profile.c_freq_5prime_raw[i] = profile.c_freq_5prime[i];
        profile.a_freq_3prime_raw[i] = profile.a_freq_3prime[i];
        profile.g_freq_3prime_raw[i] = profile.g_freq_3prime[i];
        profile.t_freq_3prime_raw[i] = profile.t_freq_3prime[i];
        profile.c_freq_3prime_raw[i] = profile.c_freq_3prime[i];
    }

    // Step 1: Normalize position-specific frequencies (but don't compute damage rates yet)
    for (int i = 0; i < 15; ++i) {
        double tc_total = profile.t_freq_5prime[i] + profile.c_freq_5prime[i];
        if (tc_total > 0) {
            profile.t_freq_5prime[i] = profile.t_freq_5prime[i] / tc_total;
            profile.c_freq_5prime[i] = 1.0 - profile.t_freq_5prime[i];
        }
        double ag_total = profile.a_freq_3prime[i] + profile.g_freq_3prime[i];
        if (ag_total > 0) {
            profile.a_freq_3prime[i] = profile.a_freq_3prime[i] / ag_total;
            profile.g_freq_3prime[i] = 1.0 - profile.a_freq_3prime[i];
        }
    }


}

} // namespace taph
