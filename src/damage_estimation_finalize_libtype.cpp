#include "damage_estimation_finalize_helpers.hpp"
#include "damage_estimation_finalize_ctx.hpp"
#include <cmath>
namespace taph {

// Hexamer-stratified terminal normalisation.
// For positions j=0..5, the observed T/(T+C) IS the 5' hexamer-composition-weighted
// T-fraction — adapter protocol biases (SapI, ligation) are fully encoded in the
// hexamer distribution.  The interior hexamer distribution carries no such bias.
// This replaces each t_freq_5prime[j] (j=0..5) with the INTERIOR-hexamer-derived
// expected T/(T+C) at position j, eliminating the adapter bias.  tc_total is
// rescaled to the interior-weighted TC count (scaled to the 5' read total so that
// downstream fits see the correct effective denominator).
// Returns true when the correction meaningfully shifted the profile.
static bool apply_adapter_deconvolution(SampleDamageProfile& profile) {
    if (profile.n_hexamers_5prime < 10000) return false;
    if (profile.n_hexamers_interior < 10000) return false;

    const double n5  = static_cast<double>(profile.n_hexamers_5prime);
    const double ni  = static_cast<double>(profile.n_hexamers_interior);
    const double scale = n5 / ni;  // scale interior counts to 5' read total
    static constexpr char kBases[4] = {'A','C','G','T'};

    // Compute interior-hexamer-based expected T and TC at each terminal position j=0..5.
    // Decode each interior hexamer code and check base at position j.
    double int_T[6]  = {};   // Σ hexamer_count_interior[h] × I(h[j]='T')
    double int_TC[6] = {};   // Σ hexamer_count_interior[h] × I(h[j] in {T,C})
    // Mirror for 3' end (pos j from 3' = base at index 5-j in hexamer)
    double int_A3[6]  = {};  // Σ hexamer_count_interior[h] × I(h[5-j]='A')
    double int_AG3[6] = {};  // Σ hexamer_count_interior[h] × I(h[5-j] in {A,G})

    for (uint32_t code = 0; code < 4096u; ++code) {
        const double cnt = profile.hexamer_count_interior[code];
        if (cnt < 1.0) continue;
        uint32_t tmp = code;
        char bases[6];
        for (int b = 5; b >= 0; --b) { bases[b] = kBases[tmp & 3u]; tmp >>= 2; }
        for (int j = 0; j < 6; ++j) {
            int_T[j]  += cnt * (bases[j] == 'T' ? 1.0 : 0.0);
            int_TC[j] += cnt * (bases[j] == 'T' || bases[j] == 'C' ? 1.0 : 0.0);
            const char b3 = bases[5 - j];
            int_A3[j]  += cnt * (b3 == 'A' ? 1.0 : 0.0);
            int_AG3[j] += cnt * (b3 == 'A' || b3 == 'G' ? 1.0 : 0.0);
        }
    }

    // Measure total adapter bias as the total absolute deviation of 5' terminal TC
    // composition from interior across positions 0-5.
    double n_adapter5 = 0.0;
    for (int j = 0; j < 6; ++j) {
        if (int_TC[j] < 1.0) continue;
        const double int_expected = int_T[j] / int_TC[j];
        // Deviation accounts for both the sign and the TC denominator magnitude
        n_adapter5 += std::abs(profile.t_freq_5prime[j] - int_expected)
                       * profile.tc_total_5prime[j];
    }

    profile.adapter_deconv_n_stub5 = n_adapter5;
    profile.adapter_deconv_p_a = static_cast<float>(n_adapter5 / (n5 + 1.0));

    // Only apply if there is a non-trivial bias (>0.5% of reads affected)
    if (profile.adapter_deconv_p_a < 5e-3f) return false;

    // Replace 5' terminal profile positions j=0..5 with interior-derived expected values.
    bool any_corrected = false;
    for (int j = 0; j < 6; ++j) {
        if (int_TC[j] < 1.0) continue;
        const double tc_new = int_TC[j] * scale;
        const double t_new  = int_T[j]  * scale;
        if (tc_new < 1.0) continue;
        profile.t_freq_5prime[j]   = t_new / tc_new;
        profile.tc_total_5prime[j] = tc_new;
        any_corrected = true;
    }
    if (!any_corrected) return false;

    // Mirror correction on 3' profile: replace a_freq_3prime/ag_total for positions j=0..5.
    if (profile.n_hexamers_3prime >= 10000) {
        const double n3 = static_cast<double>(profile.n_hexamers_3prime);
        const double scale3 = n3 / ni;
        for (int j = 0; j < 6; ++j) {
            if (int_AG3[j] < 1.0) continue;
            const double ag_new = int_AG3[j] * scale3;
            const double a_new  = int_A3[j]  * scale3;
            if (ag_new < 1.0) continue;
            profile.a_freq_3prime[j]   = a_new / ag_new;
            profile.ag_total_3prime[j] = ag_new;
        }
    }

    // Refresh damage_rate_5prime[0..5] from corrected ratios so finalize_dmax sees
    // corrected d_max. Same formula as finalize_context.
    const float bg5  = profile.fit_baseline_5prime;
    const float bg5c = 1.0f - bg5;
    if (bg5c > 0.1f) {
        for (int j = 0; j < 6; ++j) {
            const float raw = static_cast<float>(profile.t_freq_5prime[j]) - bg5;
            profile.damage_rate_5prime[j] = std::max(0.0f, raw / bg5c);
        }
    }
    float peak = 0.0f;
    for (int j = 0; j < 6; ++j) peak = std::max(peak, profile.damage_rate_5prime[j]);
    profile.d_max_5prime_deconv = peak;

    profile.adapter_deconv_applied = true;
    return true;
}

void finalize_libtype(SampleDamageProfile& profile, const FinalCtx& ctx) {
    // Library type detection.
    // DS: C→T at 5' + G→A at 3' (terminal_shift_3prime elevated).
    // SS: C→T at both ends — 3' shows elevated T/(T+C) (ctrl_shift_3prime) with no G→A.
    //
    // Library type detection via 4-model BIC comparison on the 3' end (positions 1-10).
    //
    // Four competing models:
    //   M_bias: no 3' decay in either channel (composition/ligation bias only)
    //   M_DS:   G→A decay at 3' only  (classic double-stranded aDNA)
    //   M_SS:   C→T decay at 3' only  (single-stranded: damage at both ends)
    //   M_mix:  both channels show decay (ambiguous / mixed library)
    //
    // Lambda is fixed to the fitted 5' C→T decay, so coverage inflating z-scores
    // does not affect the classification — only whether the decay shape fits.
    // Position 0 is excluded entirely (known adapter ligation artifact at 3' in SS prep).
    // P0-1: BIC tournament runs UNCONDITIONALLY. The forced-library override is
    // applied AFTER, by overwriting profile.library_type only — diagnostic
    // BIC fields, posteriors, and submodel scores remain populated so KapK-forced
    // runs produce identical numbers to KapK-auto runs.
    {
        float lambda_lib = std::clamp(ctx.fit_lambda_5p, 0.05f, 0.50f);

        // Four channels, all with lambda fixed to the fitted 5' C→T decay rate.
        //
        // ct5: T/(T+C) at 5' pos 1-10  — 5' C→T damage (DS reads and SS-original reads).
        // ga3: A/(A+G) at 3' pos 1-10  — smooth 3' G→A decay (DS reads and SS-complement reads).
        // ga0: A/(A+G) at 3' pos 0     — pos-0 G→A spike (SS-complement ligation artifact).
        // ct3: T/(T+C) at 3' pos 1-10  — 3' C→T (SS-original reads only; absent in DS).
        //
        // Seven competing joint models:
        //
        //   M_bias         = ct5_null + ga3_null  + ga0_null + ct3_null  (no damage anywhere)
        //   M_DS_symm      = ds_symm  + ga0_null  + ct3_null  (DS: ct5+ga3 joint fit, symmetric decay)
        //   M_DS_spike     = ct5_null + ga3_null  + ga0_alt  + ct3_null  (DS: only pos-0 end-repair spike)
        //   M_DS_symm_art  = ds_symm  + ga0_alt   + ct3_null  (DS: symmetric decay + pos-0 artifact)
        //   M_SS_comp      = ct5_null + ga3_alt   + ga0_alt  + ct3_null  (SS complement-orientation reads)
        //   M_SS_orig      = ct5_alt  + ga3_null  + ga0_null + ct3_alt   (SS original-orientation reads)
        //   M_SS_full      = ct5_alt  + ga3_alt   + ga0_alt  + ct3_alt   (SS both orientations)
        //
        // M_DS_symm, M_DS_spike, M_DS_symm_art classify as DS; M_SS_* as SS.
        //
        // M_DS_symm: ct5 and ga3 are jointly fitted with a SINGLE shared amplitude. For genuine DS,
        // ct5 ≈ ga3 and A_joint is near-optimal for both → low BIC. For SS_comp (ct5 << ga3), A_joint
        // is too high for ct5 (over-predicts → LL_ct5 < LL_null) and too low for ga3 → combined BIC
        // exceeds M_SS_comp → M_DS_symm rejected. This replaces the earlier threshold-based approach.
        //
        // M_DS_spike covers DS libraries with composition bias at 5' where positions 1-10 show no
        // usable ct5/ga3 signal (amplitudes clamped to 0), yet a genuine chemical pos-0 signal exists
        // from end-repair. Without M_DS_spike, M_SS_comp would win for such samples.
        //
        // Joint BIC = sum of per-channel BICs (bic_alt for alt channels, bic_null for null).
        // M_DS_symm uses 1 free parameter for both ct5+ga3; all other models use per-channel params.
        // Winner = lowest joint BIC. No post-hoc thresholds.
        // Offset search: only applied when an adapter artifact or terminal inversion was detected
        // (position_0_artifact or inverted_pattern). For normal samples always use start_pos=1
        // to avoid perturbing the BIC scores used in library-type classification.
        const bool artifact_5 = profile.position_0_artifact_5prime || profile.inverted_pattern_5prime;
        const bool artifact_3 = profile.position_0_artifact_3prime || profile.inverted_pattern_3prime;

        // GA-residual junction mask: when the 5' terminal is inverted, compute per-position
        // (CT_shift - GA_shift) residual. The GA channel cannot carry 5' deamination, so a
        // negative CT residual after GA subtraction indicates genuine T depletion at the junction
        // (cut-site GC enrichment / adapter remnant), not deamination signal. Mask the leading
        // prefix of positions where CT is still depressed OR the residual is still negative,
        // then use the mask boundary as the minimum start for all BIC offset searches.
        // Both bic_null and bic_alt in fit_decay_fixed_lambda are computed over the same
        // [min_sp, end_pos] window, so the common-support constraint is automatically satisfied.
        // 3-component adapter deconvolution: if the 5' terminal is inverted and we have
        // exact CTCTTC hexamer matches, correct the aggregate profile before masking.
        // This recovers positions 0-4 where the exponential decay carries >99% of the signal.
        if (profile.inverted_pattern_5prime)
            apply_adapter_deconvolution(profile);

        // After correction, re-evaluate whether inversion is still present.
        // The correction modifies t_freq_5prime/tc_total_5prime in-place; re-check
        // pos0 vs interior to see if the mask is still needed.
        if (profile.adapter_deconv_applied) {
            float ct_bg = 0.0f; int n_bg = 0;
            for (int p = 5; p < 15; ++p) {
                if (profile.tc_total_5prime[p] > 100.0) {
                    ct_bg += static_cast<float>(profile.t_freq_5prime[p]); ++n_bg;
                }
            }
            if (n_bg > 0) ct_bg /= n_bg;
            bool still_inverted = false;
            for (int p = 0; p < 3; ++p) {
                if (profile.tc_total_5prime[p] > 100.0 &&
                    profile.t_freq_5prime[p] < ct_bg - 0.02)
                    still_inverted = true;
            }
            if (!still_inverted) profile.inverted_pattern_5prime = false;
        }

        int ga_mask_5 = 0;
        if (profile.inverted_pattern_5prime) {
            // t_freq_5prime is already a fraction (converted in-place by finalize_init).
            // a_freq_5prime / g_freq_5prime are still raw counts — use as numerator/denominator.
            float ct_bg_m = 0.0f, ga_bg_m = 0.0f; int n_int_m = 0;
            for (int p = 5; p < 15; ++p) {
                double tc = profile.tc_total_5prime[p];
                double ag = profile.a_freq_5prime[p] + profile.g_freq_5prime[p];
                if (tc > 0 && ag > 0) {
                    ct_bg_m += static_cast<float>(profile.t_freq_5prime[p]);  // already T/(T+C)
                    ga_bg_m += static_cast<float>(profile.a_freq_5prime[p] / ag);
                    ++n_int_m;
                }
            }
            if (n_int_m > 0) { ct_bg_m /= n_int_m; ga_bg_m /= n_int_m; }
            for (int p = 0; p < 5; ++p) {
                double tc = profile.tc_total_5prime[p];
                double ag = profile.a_freq_5prime[p] + profile.g_freq_5prime[p];
                if (tc < 100.0 || ag < 100.0) { ga_mask_5 = p + 1; continue; }
                float ct_p  = static_cast<float>(profile.t_freq_5prime[p]);  // already T/(T+C)
                float ga_p  = static_cast<float>(profile.a_freq_5prime[p] / ag);
                float residual = (ct_p - ct_bg_m) - (ga_p - ga_bg_m);
                if ((ct_p - ct_bg_m) < -0.02f || residual < -0.01f)
                    ga_mask_5 = p + 1;
                else
                    break;
            }
            profile.junction_mask_n_5prime = ga_mask_5;
        }
        const int min_sp_5 = std::max(1, ga_mask_5);

        ChannelDecayFit ct5;  int ct5_offset = 1;
        ChannelDecayFit ga3;  int ga3_offset = 1;
        ChannelDecayFit ct3;  int ct3_offset = 1;
        ChannelDecayFit ds_symm;  int ds_symm_offset = 1;
        if (artifact_5) {
            std::tie(ct5, ct5_offset) = fit_decay_best_offset(
                profile.t_freq_5prime, profile.tc_total_5prime,
                static_cast<float>(ctx.baseline_tc), lambda_lib, 10, 3, min_sp_5);
        } else {
            ct5 = fit_decay_fixed_lambda(profile.t_freq_5prime, profile.tc_total_5prime,
                static_cast<float>(ctx.baseline_tc), lambda_lib, 1, 10);
        }
        if (artifact_3) {
            std::tie(ga3, ga3_offset) = fit_decay_best_offset(
                profile.a_freq_3prime, profile.ag_total_3prime,
                static_cast<float>(ctx.baseline_ag), lambda_lib, 10);
            // ct3 does NOT use offset search even when artifact_3 is true: the 3' adapter
            // artifact suppresses GA3 (complement-strand G→A) but not CT3 (original-strand
            // C→T). Offset search for ct3 finds spurious C→T signals in adapter-affected DS
            // libraries and drives false SS-original classifications.
            ct3 = fit_decay_fixed_lambda(ctx.ctrl_freq_3p, ctx.ctrl_total_3p,
                static_cast<float>(ctx.baseline_tc), lambda_lib, 1, 10);
        } else {
            ga3 = fit_decay_fixed_lambda(profile.a_freq_3prime, profile.ag_total_3prime,
                static_cast<float>(ctx.baseline_ag), lambda_lib, 1, 10);
            ct3 = fit_decay_fixed_lambda(ctx.ctrl_freq_3p, ctx.ctrl_total_3p,
                static_cast<float>(ctx.baseline_tc), lambda_lib, 1, 10);
        }
        // ga0: single-position fit at 3' pos-0; no offset search (it IS the pos-0 spike signal)
        ChannelDecayFit ga0 = fit_decay_fixed_lambda(
            profile.a_freq_3prime, profile.ag_total_3prime,
            static_cast<float>(ctx.baseline_ag), lambda_lib, 0, 0, 1);
        // Symmetric DS model: offset search only when either end has artifact.
        // Restrict high-lambda candidates when ga0 dominates ct5 (SS-complement indicator):
        // the high-lambda CT5 fit should not inflate the DS joint BIC when the sample already
        // shows a large GA0 spike (spike_is_ss-like pattern). ga0 is fit before this point.
        const bool restrict_joint_lambda = (ga0.amplitude >= 0.10f && ga0.delta_bic > ct5.delta_bic);
        // P2: spike-gate diagnostic — record the joint-lambda gating decision
        profile.library_joint_lambda_restricted = restrict_joint_lambda;
        if (artifact_5 || artifact_3) {
            std::tie(ds_symm, ds_symm_offset) = fit_decay_joint_best_offset(
                profile.t_freq_5prime, profile.tc_total_5prime, static_cast<float>(ctx.baseline_tc),
                profile.a_freq_3prime, profile.ag_total_3prime, static_cast<float>(ctx.baseline_ag),
                lambda_lib, 10, 3, restrict_joint_lambda, min_sp_5);
        } else {
            ds_symm = fit_decay_fixed_lambda_joint(
                profile.t_freq_5prime, profile.tc_total_5prime, static_cast<float>(ctx.baseline_tc),
                profile.a_freq_3prime, profile.ag_total_3prime, static_cast<float>(ctx.baseline_ag),
                lambda_lib, 1, 10);
        }

        // Store detected offsets; use the single-channel offsets as the canonical per-end values
        // (ds_symm_offset is derived from the joint fit and may differ from individual channels)
        profile.fit_offset_5prime = ct5_offset;
        profile.fit_offset_3prime = ga3_offset;

        // When adapter artifact or terminal inversion is detected, scan positions 1-5 for the peak
        // damage rate. Using damage_rate[fit_offset] is unreliable because the BIC-best offset may
        // point to a position still below baseline (e.g. BPN103cm: fit_offset=3 but pos3 < bg,
        // while actual biological damage is at pos1). Scan-for-peak handles all shift lengths.
        if (profile.position_0_artifact_5prime || profile.inverted_pattern_5prime) {
            float peak = 0.0f;
            for (int p = 1; p <= 5 && p < 15; ++p) {
                if (profile.damage_rate_5prime[p] > peak) peak = profile.damage_rate_5prime[p];
            }
            if (peak > 0.0f) profile.max_damage_5prime = peak;
        }
        if (profile.position_0_artifact_3prime || profile.inverted_pattern_3prime) {
            float peak = 0.0f;
            for (int p = 1; p <= 5 && p < 15; ++p) {
                if (profile.damage_rate_3prime[p] > peak) peak = profile.damage_rate_3prime[p];
            }
            if (peak > 0.0f) profile.max_damage_3prime = peak;
        }
        (void)ct3_offset; (void)ds_symm_offset;

        profile.libtype_amp_ct5  = ct5.amplitude;
        profile.libtype_amp_ga3  = ga3.amplitude;
        profile.libtype_amp_ga0  = ga0.amplitude;
        profile.libtype_amp_ct3  = ct3.amplitude;
        profile.libtype_dbic_ct5 = ct5.delta_bic;
        profile.libtype_dbic_ga3 = ga3.delta_bic;
        profile.libtype_dbic_ga0 = ga0.delta_bic;
        profile.libtype_dbic_ct3 = ct3.delta_bic;

        // Wave-2b: the SINGLE sink for every post-hoc change to library_type. The BIC
        // tournament (below) is the sole producer; thereafter no rule writes library_type
        // directly -- each rescue/veto calls apply_override, which records {rule_id, from, to}
        // so the divergence from library_bic_call is a witnessed ledger. Declared at function
        // scope so it also covers the misspec UNKNOWN->DS rescue outside the BIC-evaluable block.
        auto apply_override = [&profile](const char* rule_id,
                                         SampleDamageProfile::LibraryType to,
                                         bool sets_rescued) {
            profile.library_overrides.push_back({rule_id, profile.library_type, to});
            profile.library_type = to;
            if (sets_rescued) profile.library_type_rescued = true;
        };
        if (ct5.valid && ga3.valid && ga0.valid && ct3.valid && ds_symm.valid) {
            // Hard biological gates (no tunable parameters):
            //   M_DS_symm:     DS damage is symmetric — requires ga3 to have real
            //                  signal. Use the standard Bayes-factor "no evidence"
            //                  threshold: ga3.delta_bic > log(2) ≈ 0.693
            //                  (Kass & Raftery 1995). If ga3 is indistinguishable
            //                  from null, there is nothing for ct5 to be symmetric
            //                  with — the joint ds_symm fit is fitting noise.
            //   M_DS_symm_art: 3' GA spike is the complementary-strand reflection
            //                  of 5' CT damage. A reflection cannot exceed its
            //                  source — require ga0 ≤ ct5.
            // Exception: when position_0_artifact_3prime is flagged, the 3' pos-0
            // base-composition artifact legitimately inflates ga0 above ct5.
            // In that case allow M_DS_symm_art to enter the cascade so it can
            // compete against M_SS_comp (which otherwise wins unopposed).
            // When violated (and no artifact flag), the model is invalid for this
            // sample and excluded from the tournament (BIC = +inf).
            const bool ds_symm_valid     = (ga3.delta_bic > std::log(2.0));
            const bool ds_symm_art_valid = (ga0.amplitude <= ct5.amplitude)
                                        || profile.position_0_artifact_3prime;
            constexpr double kInvalidBIC = 1e300;

            const double bic_M_bias        = ct5.bic_null   + ga3.bic_null + ga0.bic_null + ct3.bic_null;
            const double bic_M_DS_symm     = ds_symm_valid     ? (ds_symm.bic_alt + ga0.bic_null + ct3.bic_null) : kInvalidBIC;
            const double bic_M_DS_spike    = ct5.bic_null   + ga3.bic_null + ga0.bic_alt  + ct3.bic_null;
            const double bic_M_DS_symm_art = ds_symm_art_valid ? (ds_symm.bic_alt + ga0.bic_alt  + ct3.bic_null) : kInvalidBIC;
            const double bic_M_SS_comp     = ct5.bic_null   + ga3.bic_alt  + ga0.bic_alt  + ct3.bic_null;
            const double bic_M_SS_orig     = ct5.bic_alt    + ga3.bic_null + ga0.bic_null + ct3.bic_alt;
            const double bic_M_SS_full     = ct5.bic_alt    + ga3.bic_alt  + ga0.bic_alt  + ct3.bic_alt;
            // S1 telemetry only: asymmetric DS counterfactual (independent ct5/ga3 amps,
            // GA0 artifact, no CT3). Never enters cascade/softmax.
            const double bic_M_DS_asym_art = ct5.bic_alt    + ga3.bic_alt  + ga0.bic_alt  + ct3.bic_null;
            // S1 invariant probe: M_DS_symm_art rebuilt from a forced no-offset joint fit
            // (start_pos=1, no joint best-offset search) regardless of artifact_5/artifact_3.
            // Catches future regressions where joint best-offset search inflates ds_symm.
            ChannelDecayFit ds_symm_no_off = fit_decay_fixed_lambda_joint(
                profile.t_freq_5prime, profile.tc_total_5prime, static_cast<float>(ctx.baseline_tc),
                profile.a_freq_3prime, profile.ag_total_3prime, static_cast<float>(ctx.baseline_ag),
                lambda_lib, 1, 10);
            const double bic_M_DS_symm_art_no_offset = (ds_symm_no_off.valid
                ? ds_symm_no_off.bic_alt + ga0.bic_alt + ct3.bic_null
                : 0.0);
            // SS asymmetric: original-orientation CT5 + complement-orientation GA0 spike; no GA3 smooth decay.
            // Wins over M_DS_symm_art when ga3.delta_bic < log(2), i.e. ga3 has no real signal.
            // Only considered when spike_is_ss=true so it cannot compete with M_DS_spike in the DS-only path.
            const double bic_M_SS_asym     = ct5.bic_alt    + ga3.bic_null + ga0.bic_alt  + ct3.bic_null;

            // ga0 amplitude distinguishes DS end-repair artifact (<0.10) from SS ligation spike (>=0.10).
            // Exception: if the channel-B structural analysis (computed before BIC, line ~1200)
            // confirmed bilateral symmetric damage with d_max_from_channel_b > 0.10, then the
            // 3' pos-0 GA excess is the mirror of the bilateral 5' C→T on the bottom strand —
            // not an SS ligation artifact.  Without this guard, highly damaged DS libraries
            // (d_max > 0.20, steep lambda) produce large GA0_ΔBIC from the 3' bilateral
            // reflection and are wrongly classified as SS.
            const bool structural_bilateral = profile.channel_b_quantifiable
                                           && (profile.d_max_from_channel_b > 0.10f);
            const bool spike_is_ss = (ga0.amplitude >= 0.10f) && !structural_bilateral;
            // P2: spike-gate diagnostics — record the gating inputs and decision
            profile.library_spike_is_ss                    = spike_is_ss;
            profile.library_spike_gate_ga0_amp             = ga0.amplitude;
            profile.library_spike_gate_structural_bilateral = structural_bilateral;
            // P1-B: store all 7 submodel BICs and active-flags so callers can audit the tournament
            profile.library_bic_M_bias        = bic_M_bias;
            profile.library_bic_M_DS_symm     = bic_M_DS_symm;
            profile.library_bic_M_DS_spike    = bic_M_DS_spike;
            profile.library_bic_M_DS_symm_art = bic_M_DS_symm_art;
            profile.library_bic_M_SS_comp     = bic_M_SS_comp;
            profile.library_bic_M_SS_orig     = bic_M_SS_orig;
            profile.library_bic_M_SS_asym     = bic_M_SS_asym;
            profile.library_bic_M_SS_full     = bic_M_SS_full;
            profile.library_M_SS_orig_active  = (ct3.delta_bic > 0.0);
            profile.library_M_SS_asym_active  = spike_is_ss;
            // S1 telemetry: counterfactual + invariant probe BICs (never enter cascade)
            profile.library_bic_M_DS_asym_art          = bic_M_DS_asym_art;
            profile.library_bic_M_DS_symm_art_no_offset = bic_M_DS_symm_art_no_offset;
            // S1 telemetry: gate inputs / per-channel offsets / ds_symm joint diagnostics
            profile.library_gate_artifact_5                   = artifact_5;
            profile.library_gate_artifact_3                   = artifact_3;
            profile.library_gate_position_0_artifact_5prime   = profile.position_0_artifact_5prime;
            profile.library_gate_position_0_artifact_3prime   = profile.position_0_artifact_3prime;
            profile.library_gate_inverted_pattern_5prime      = profile.inverted_pattern_5prime;
            profile.library_gate_inverted_pattern_3prime      = profile.inverted_pattern_3prime;
            profile.library_gate_max_damage_5prime            = profile.max_damage_5prime;
            profile.library_gate_structural_bilateral         = structural_bilateral;
            profile.library_gate_ga0_dominates_ct5            =
                (ga0.amplitude >= 0.10f && ga0.delta_bic > ct5.delta_bic);
            profile.library_ct3_offset      = ct3_offset;
            profile.library_ds_symm_offset  = ds_symm_offset;
            profile.library_ds_symm_lambda_used = ds_symm.lambda;
            profile.library_ds_symm_amp     = ds_symm.amplitude;
            profile.library_ds_symm_ct5_resid = ct5.amplitude - ds_symm.amplitude;
            profile.library_ds_symm_ga3_resid = ga3.amplitude - ds_symm.amplitude;
            const double best_ds = std::min({bic_M_DS_symm,
                                             bic_M_DS_symm_art,
                                             spike_is_ss ? 1e300 : bic_M_DS_spike});
            // M_SS_full excluded: 4-param model unfairly defeats M_DS_symm_art for asymmetric DS.
            // M_SS_asym only enters the SS set when spike_is_ss, preventing competition with M_DS_spike
            // in non-spike samples (which would misclassify one-sided DS libraries).
            const double best_ss = std::min({bic_M_SS_comp,
                                             ct3.delta_bic > 0.0 ? bic_M_SS_orig : 1e300,
                                             spike_is_ss ? bic_M_DS_spike : 1e300,
                                             spike_is_ss ? bic_M_SS_asym  : 1e300});
            profile.library_bic_bias = bic_M_bias;
            profile.library_bic_ds   = best_ds;
            profile.library_bic_ss   = best_ss;
            profile.library_bic_mix  = bic_M_SS_full;

            // Softmax posterior probabilities over all 7 candidate models:
            //   P(M_i | data) ∝ exp(-BIC_i/2)
            // Subtracting min(BIC) before exp() keeps everything finite under
            // the typical BIC range we see (Δ ~ 1e3–1e6 between models). Sum
            // by class — DS includes M_DS_symm/spike/symm_art; SS includes
            // M_SS_comp/orig/full/asym (the latter two only when active).
            // M_DS_spike contributes to SS instead when spike_is_ss=true,
            // matching the cascade's interpretation of bilateral pos-0 GA.
            double best = bic_M_bias;
            // Wave-2b: the BIC tournament writes a LOCAL bic_type; profile.library_type is
            // assigned exactly ONCE below (and frozen as library_bic_call), so the producer is
            // a single writer. Every post-hoc override then flows through apply_override.
            SampleDamageProfile::LibraryType bic_type =
                SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            bool ds_spike_won = false;  // tracks whether M_DS_spike is the current winning model

            // Bilateral gate: M_DS_symm requires both CT5 and GA3 channels to carry
            // positive BIC evidence. A one-sided CT5-only pattern is better captured by
            // M_DS_symm_art (which explicitly models the missing 3' smooth decay).
            if (bic_M_DS_symm < best && ct5.delta_bic > 0.0 && ga3.delta_bic > 0.0) {
                best = bic_M_DS_symm;
                bic_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            }
            if (!spike_is_ss && bic_M_DS_spike < best) {
                best = bic_M_DS_spike;
                bic_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
                ds_spike_won = true;
            }
            if (bic_M_DS_symm_art < best) {
                best = bic_M_DS_symm_art;
                bic_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
                ds_spike_won = false;
            }
            if (bic_M_SS_comp < best) {
                best = bic_M_SS_comp;
                bic_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
                ds_spike_won = false;
            }
            // M_SS_orig requires ct3 signal: SS original-orientation reads produce CT3 whenever
            // they produce CT5. Without CT3, a one-sided DS pattern is more likely.
            if (bic_M_SS_orig < best && ct3.delta_bic > 0.0) {
                best = bic_M_SS_orig;
                bic_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
                ds_spike_won = false;
            }
            if (spike_is_ss && bic_M_DS_spike < best) {
                best = bic_M_DS_spike;
                bic_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
                ds_spike_won = false;
            }
            // M_SS_asym: SS with CT5 from original-orientation + GA0 spike from complement-orientation,
            // but no detectable GA3 smooth decay (ga3.delta_bic ≈ 0). Analytically beats M_DS_symm_art
            // iff ga3.delta_bic < log(2) ≈ 0.693, which is the gap the joint-fit BIC penalty creates.
            if (spike_is_ss && bic_M_SS_asym < best) {
                best = bic_M_SS_asym;
                bic_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            }
            // Wave-2: the single producer write + freeze the pure-BIC argmin verdict BEFORE any
            // rescue mutates the call, so the emitted library_type's divergence from BIC is
            // witnessed (the rescues below are apply_override records, not silent edits).
            profile.library_type     = bic_type;
            profile.library_bic_call = bic_type;
            // Post-hoc symmetry check: DS_symm constrains ct5_amp ≈ ga3_amp.
            // If DS wins but CT5 ΔBIC / GA3 ΔBIC < 0.50, the winning model's own symmetry
            // assumption is violated → reclassify as SS. Guard ga3.delta_bic > 3e4 to avoid
            // misfiring on low-damage DS libraries where small asymmetry is noise.
            // Symmetry veto: DS won BIC, but GA3 >> CT5 asymmetry suggests SS.
            // Only apply when there is positive biological evidence for SS orientation:
            //   ga0.delta_bic > 1e4 — meaningful GA0 ligation spike (amplitude-independent,
            //                         avoids false negatives from samples near the 0.10 amplitude
            //                         threshold that spike_is_ss uses)
            //   inverted_5prime     — 5' CT profile depressed below interior (SS adapter suppression)
            // Without either, the asymmetry likely reflects 3'-adapter inflation of GA3_ΔBIC in a DS
            // library (e.g., high-λ fit capturing sharp pos2+ GA3 decay after 3' adapter suppression).
            const bool ss_orientation_evidence = (ga0.delta_bic > 1e4)
                                              || profile.inverted_pattern_5prime;
            if (profile.library_type == SampleDamageProfile::LibraryType::DOUBLE_STRANDED &&
                ss_orientation_evidence &&
                ga3.delta_bic > 3e4 &&
                ct5.delta_bic / ga3.delta_bic < 0.50) {
                apply_override("symmetry_veto_ds_to_ss",
                               SampleDamageProfile::LibraryType::SINGLE_STRANDED, false);
                ds_spike_won = false;
            }
            // M_DS_spike rescue: a GA0 pos-0 spike with no CT5 and no GA3 smooth decay could be
            // M_DS_spike (DS end-repair, bilateral: both 5' pos0 and 3' pos0 elevated) or
            // complement-orientation SS (3' GA0 spike only, no 5' pos-0 counterpart).
            // Discriminate with max_damage_5prime (background-corrected excess CT at 5' pos-0).
            // Requires ga0.amplitude > 0.02 to exclude near-zero noise spikes (e.g. tiny end-repair
            // artifacts with ga0_amp ≈ 0.001 that also have negligible d_max_5 ≈ 0.005).
            // DS bilateral artifacts have d_max_5 ≈ ga0_amp >> 0.005; SS complement-only has d_max_5 ≈ 0.
            // Restricted to ds_spike_won: M_DS_symm_art can win via joint fit even with marginal
            // ct5/ga3 delta_bic ≤ 0; rescue only fires when M_DS_spike was the actual winner.
            if (profile.library_type == SampleDamageProfile::LibraryType::DOUBLE_STRANDED &&
                ds_spike_won &&
                !spike_is_ss &&
                ct5.delta_bic <= 0.0 &&
                ga3.delta_bic <= 0.0 &&
                ga0.delta_bic > 0.0 &&
                ga0.amplitude > 0.02f &&
                profile.max_damage_5prime <= 0.005f) {
                apply_override("ds_spike_complement_to_ss",
                               SampleDamageProfile::LibraryType::SINGLE_STRANDED, false);
            }
            // GA0 bilateral rescue (spike_is_ss path): when ga0.amplitude >= 0.10 and DS wins,
            // discriminate DS end-repair (bilateral: both 5' pos-0 CT and 3' pos-0 GA elevated)
            // from SS complement-orientation reads (3' GA0 spike only, no 5' CT0 counterpart).
            // Validated on 24 DS controls with ga0_amp >= 0.10: all have d5 >= 0.11.
            // All SS mixed failures with ga0_amp >= 0.10 have d5 = 0. Gap is >20x the threshold.
            if (profile.library_type == SampleDamageProfile::LibraryType::DOUBLE_STRANDED &&
                spike_is_ss &&
                ga0.delta_bic > 0.0 &&
                profile.max_damage_5prime <= 0.005f) {
                apply_override("ga0_bilateral_to_ss",
                               SampleDamageProfile::LibraryType::SINGLE_STRANDED, false);
            }
            // Channel-B structural DS rescue: if channel-B (stop codon conversion, immune to
            // composition bias) confirmed bilateral symmetric damage > 0.10, AND the 3' ligation
            // spike (GA0) dominates the smooth GA3 decay by ≥5×, override SS→DS.
            // Handles highly damaged DS libraries where steep lambda suppresses CT5 in the BIC
            // fit, causing M_SS_comp to win by capturing GA3+GA0 independently. Channel B is
            // only quantifiable when real C→T bilateral damage exists at 5'; SS-complement
            // libraries have ch_b_quant=False and are not affected.
            // The GA0 dominance guard (ga3 * 5 < ga0) excludes SS-orig libraries: they have
            // real bilateral-looking channels (ch_b_quant=True, large GA3) but their GA3 decay
            // is comparable to GA0 rather than being dwarfed by it.
            const bool ga_spike_dominant = (ga3.delta_bic * 5.0 < ga0.delta_bic);
            if (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED
                && structural_bilateral && ga_spike_dominant
                && ct3.delta_bic <= 0.0) {
                apply_override("channel_b_structural_ds",
                               SampleDamageProfile::LibraryType::DOUBLE_STRANDED, true);
            }
            // GA0-spike DS-symm veto: when the 3' ga0 ligation spike unambiguously
            // dominates (ga0 >= 0.10, ga0_dominates_ct5, ga0 > both ct5 and ga3),
            // an apparent CT5/GA3 symmetry absorbed by M_DS_symm is artifact, not
            // real DS damage -- true DS has ct5~=ga3 >> ga0. Restricted to ds_symm
            // winners; independent of max_damage_5prime (the artifact inflates it).
            // Exception: when position_0_artifact_3prime is flagged the inflated ga0
            // is a known composition artifact; M_DS_symm_art already models it via
            // ga0.bic_alt, so the veto must not fire and undo that correct win.
            if (profile.library_type == SampleDamageProfile::LibraryType::DOUBLE_STRANDED &&
                !profile.position_0_artifact_3prime &&
                profile.library_gate_ga0_dominates_ct5 &&
                ga0.amplitude >= 0.10f &&
                ga0.amplitude > std::max(ct5.amplitude, ga3.amplitude) &&
                (best == bic_M_DS_symm || best == bic_M_DS_symm_art)) {
                apply_override("ga0_spike_dssymm_veto",
                               SampleDamageProfile::LibraryType::SINGLE_STRANDED, false);
            }
            // Channel-B DS rescue from M_SS_comp: structural_bilateral confirms
            // bilateral symmetric damage (channel B is artifact-immune); ct5≈ga3
            // (within 30%, both ≥ 0.03) confirms the symmetric DS pair. M_SS_comp
            // wins by capturing ct5+ga3+ga0 as 3 independent bumps; M_DS_symm has
            // no ga0 term and loses on raw BIC despite being the correct model
            // (the ga0 here is end-repair / ligation residue co-occurring with
            // real DS damage, not the SS complement-orientation signature).
            if (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED &&
                structural_bilateral &&
                best == bic_M_SS_comp &&
                ct3.delta_bic <= 0.0 &&
                ct5.amplitude >= 0.03f &&
                ga3.amplitude >= 0.03f &&
                std::abs(ct5.amplitude - ga3.amplitude) <=
                    0.30f * std::max(ct5.amplitude, ga3.amplitude)) {
                apply_override("channel_b_ds_from_sscomp",
                               SampleDamageProfile::LibraryType::DOUBLE_STRANDED, true);
            }
            // Low-amp symmetric DS rescue: zero 3' damage collapses the per-end
            // DS amplitudes so M_SS_full / M_SS_orig wins raw BIC, but the joint
            // M_DS_symm fit is still valid (ct5 ~= ga3, residuals ~= 0). ga0 < 0.005
            // and !spike_is_ss exclude SS-complement orientation.
            if (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED &&
                !spike_is_ss &&
                ds_symm.amplitude >= 0.005f &&
                ct5.amplitude >= 0.005f && ga3.amplitude >= 0.005f &&
                ga0.amplitude < 0.005f &&
                std::abs(ct5.amplitude - ga3.amplitude) <=
                    0.50f * std::max(ct5.amplitude, ga3.amplitude) &&
                profile.max_damage_3prime < 0.005f) {
                apply_override("low_amp_symmetric_ds",
                               SampleDamageProfile::LibraryType::DOUBLE_STRANDED, true);
            }
            // S1 telemetry: capture cascade-derived gate outputs.
            profile.library_gate_ss_orientation_evidence = ss_orientation_evidence;
            profile.library_gate_ga_spike_dominant       = ga_spike_dominant;
            profile.library_gate_ds_spike_won            = ds_spike_won;
            // Uninformative: if nothing beat M_bias (best unchanged), no damage channel
            // provided evidence for either DS or SS. Use exact equality — best is only
            // updated via assignment from another BIC value, so equality is safe here.
            // Does NOT affect low-damage DS libraries where M_DS_symm still beats M_bias.
            if (best == bic_M_bias) {
                apply_override("uninformative_bias_to_unknown",
                               SampleDamageProfile::LibraryType::UNKNOWN, false);
            }
            // Protocol-tag rescue: top 5' hexamer is a chemistry fingerprint
            // (deterministic per library prep, independent of damage shape).
            // Runs AFTER the bias gate so chemistry evidence overrides UNKNOWN
            // for low-damage libraries where BIC has nothing to fit. Only fires
            // when log2fc >= 3.0 (8x enrichment) and the hex matches the curated
            // table, so noise hexamers cannot trip it.
            {
                auto enriched = compute_hex_enriched_5prime(profile, 3.0f);
                if (!enriched.empty()) {
                    auto seq = decode_hex(enriched[0].idx);
                    const ProtocolTag* tag = lookup_protocol_tag(seq.data());
                    if (tag) {
                        profile.protocol_tag_5prime   = std::string(seq.data(), 6);
                        profile.protocol_tag_protocol = tag->protocol;
                        profile.protocol_tag_class    = tag->klass;
                        profile.protocol_tag_log2fc   = static_cast<float>(enriched[0].log2fc);
                        profile.protocol_tag_log_lr   = tag->log_lr;
                        // Do not apply a tentative protocol tag (log_lr < 4.0) to
                        // override a structural_bilateral DS call -- channel B bilateral
                        // confirmation is stronger evidence than a tentative hexamer match.
                        const bool tag_overrides = (profile.library_type != tag->klass)
                            && !(structural_bilateral
                                 && tag->log_lr < 4.0f
                                 && tag->klass == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
                        if (tag_overrides) {
                            apply_override("protocol_tag_5prime", tag->klass, true);
                            profile.protocol_tag_applied = true;
                            profile.library_artifact_reasons.push_back("protocol_tag_5prime");
                        }
                    }
                }
            }
            // Posterior class probabilities P(class | data) ∝ exp(-BIC/2),
            // computed AFTER the cascade so M_DS_spike (and M_SS_asym) are
            // routed to whichever class the cascade's domain rescues actually
            // assigned. Keeps p_ds/p_ss consistent with library_type.
            {
                // Hypothesis set must match the cascade's competition exactly:
                //   M_bias, M_DS_symm, M_DS_spike, M_DS_symm_art, M_SS_comp,
                //   M_SS_orig (only when ct3.delta_bic > 0),
                //   M_SS_asym (only when spike_is_ss).
                // M_SS_full is intentionally NOT in the cascade (its 4-param fit
                // unfairly defeats M_DS_symm_art on asymmetric DS) so it is
                // excluded from the softmax too — otherwise its mass would
                // appear in p_ss without ever being a valid winner.
                const bool inc_orig = (ct3.delta_bic > 0.0);
                const bool inc_asym = spike_is_ss;
                const double BIG = 1e300;
                const double bics[7] = {
                    bic_M_bias,
                    bic_M_DS_symm,
                    bic_M_DS_spike,
                    bic_M_DS_symm_art,
                    bic_M_SS_comp,
                    inc_orig ? bic_M_SS_orig : BIG,
                    inc_asym ? bic_M_SS_asym : BIG
                };
                double bmin = bics[0];
                for (int i = 1; i < 7; ++i) if (bics[i] < bmin) bmin = bics[i];
                double w[7];
                double Z = 0.0;
                for (int i = 0; i < 7; ++i) {
                    double half = (bmin - bics[i]) * 0.5;
                    if (half < -80.0) half = -80.0;
                    w[i] = std::exp(half);
                    Z += w[i];
                }
                if (Z > 0.0) {
                    for (int i = 0; i < 7; ++i) w[i] /= Z;
                    const double p_bias = w[0];
                    double p_ds = w[1] + w[3];                 // M_DS_symm + M_DS_symm_art
                    double p_ss = w[4] + w[5] + w[6];          // M_SS_comp + M_SS_orig (gated) + M_SS_asym (gated)
                    // Route ambiguous M_DS_spike (w[2]) to the class the cascade
                    // (with all rescues applied) ended up choosing. UNKNOWN/bias
                    // falls back to the spike_is_ss heuristic.
                    const auto LT = profile.library_type;
                    if (LT == SampleDamageProfile::LibraryType::SINGLE_STRANDED) {
                        p_ss += w[2];
                    } else if (LT == SampleDamageProfile::LibraryType::DOUBLE_STRANDED) {
                        p_ds += w[2];
                    } else {
                        if (spike_is_ss) p_ss += w[2]; else p_ds += w[2];
                    }
                    double sum = p_bias + p_ds + p_ss;
                    if (sum > 0.0) {
                        profile.library_p_bias   = static_cast<float>(p_bias / sum);
                        profile.library_p_ds     = static_cast<float>(p_ds   / sum);
                        profile.library_p_ss     = static_cast<float>(p_ss   / sum);
                        profile.library_p_winner = std::max({profile.library_p_bias,
                                                             profile.library_p_ds,
                                                             profile.library_p_ss});
                        profile.library_type_evaluable = true;
                    }
                }
            }
            // P1-A: winner / second-best model + margin and class-min softmax.
            // The candidate set MUST match the cascade exactly (same gating as
            // the existing posterior block above): always-on M_bias/M_DS_symm/
            // M_DS_spike/M_DS_symm_art/M_SS_comp; M_SS_orig only when ct3
            // shows positive evidence; M_SS_asym only when spike_is_ss.
            {
                struct CandM { const char* name; double bic; bool active; int klass; };
                // klass: 0=bias, 1=DS, 2=SS. M_DS_spike routes per cascade outcome.
                int spike_klass = 1;  // default DS
                if (spike_is_ss) {
                    spike_klass = 2;
                } else if (profile.library_type ==
                           SampleDamageProfile::LibraryType::SINGLE_STRANDED) {
                    spike_klass = 2;
                }
                CandM cands[7] = {
                    {"M_bias",        bic_M_bias,        true,                       0},
                    {"M_DS_symm",     bic_M_DS_symm,     true,                       1},
                    {"M_DS_spike",    bic_M_DS_spike,    true,                       spike_klass},
                    {"M_DS_symm_art", bic_M_DS_symm_art, true,                       1},
                    {"M_SS_comp",     bic_M_SS_comp,     true,                       2},
                    {"M_SS_orig",     bic_M_SS_orig,     ct3.delta_bic > 0.0,        2},
                    {"M_SS_asym",     bic_M_SS_asym,     spike_is_ss,                2},
                };
                int win_i = -1, sec_i = -1;
                for (int i = 0; i < 7; ++i) {
                    if (!cands[i].active) continue;
                    if (win_i < 0 || cands[i].bic < cands[win_i].bic) {
                        sec_i = win_i; win_i = i;
                    } else if (sec_i < 0 || cands[i].bic < cands[sec_i].bic) {
                        sec_i = i;
                    }
                }
                if (win_i >= 0) {
                    profile.library_bic_winner_model = cands[win_i].name;
                    if (sec_i >= 0) {
                        profile.library_bic_second_model = cands[sec_i].name;
                        profile.library_bic_margin = cands[sec_i].bic - cands[win_i].bic;
                        // MED: raw margin scales with read depth (ΔBIC ≈ n·2·ΔLLR/read)
                        // and is not comparable across libraries. Emit a coverage-
                        // normalised companion (bg_denominator_5prime is the winning
                        // model's 5' trial count, set earlier in finalize).
                        profile.library_bic_margin_per_read =
                            (profile.bg_denominator_5prime > 0.0)
                                ? profile.library_bic_margin / profile.bg_denominator_5prime
                                : 0.0;
                    }
                }
                // Class-min softmax: best-BIC representative per class only.
                // 3-way over { best DS, best SS, M_bias }, subtract-min trick.
                double best_per_class[3] = {1e300, 1e300, 1e300};
                for (int i = 0; i < 7; ++i) {
                    if (!cands[i].active) continue;
                    int k = cands[i].klass;
                    if (cands[i].bic < best_per_class[k]) best_per_class[k] = cands[i].bic;
                }
                double bmin_c = std::min({best_per_class[0], best_per_class[1], best_per_class[2]});
                double w_c[3] = {0.0, 0.0, 0.0};
                double Z_c = 0.0;
                for (int k = 0; k < 3; ++k) {
                    if (best_per_class[k] >= 1e299) continue;
                    double half = (bmin_c - best_per_class[k]) * 0.5;
                    if (half < -80.0) half = -80.0;
                    w_c[k] = std::exp(half);
                    Z_c += w_c[k];
                }
                if (Z_c > 0.0) {
                    profile.library_p_bias_class_min = static_cast<float>(w_c[0] / Z_c);
                    profile.library_p_ds_class_min   = static_cast<float>(w_c[1] / Z_c);
                    profile.library_p_ss_class_min   = static_cast<float>(w_c[2] / Z_c);
                }
            }
            // S1 telemetry: raw 9-candidate ranking — ignores cascade gating, exposes
            // the absolute-best fit so callers can audit cascade exclusions.
            // 8 cascade contenders + M_DS_asym_art (telemetry-only counterfactual).
            {
                struct Raw { const char* name; double bic; const char* klass; bool in_cascade; };
                const Raw raw[9] = {
                    {"M_bias",        bic_M_bias,        "bias", true},
                    {"M_DS_symm",     bic_M_DS_symm,     "ds",   true},
                    {"M_DS_spike",    bic_M_DS_spike,    "ds",   true},
                    {"M_DS_symm_art", bic_M_DS_symm_art, "ds",   true},
                    {"M_SS_comp",     bic_M_SS_comp,     "ss",   true},
                    {"M_SS_orig",     bic_M_SS_orig,     "ss",   true},
                    {"M_SS_asym",     bic_M_SS_asym,     "ss",   true},
                    {"M_SS_full",     bic_M_SS_full,     "ss",   false},  // structurally excluded from cascade
                    {"M_DS_asym_art", bic_M_DS_asym_art, "ds",   false},  // telemetry-only counterfactual
                };
                int rwin = 0;
                for (int i = 1; i < 9; ++i) if (raw[i].bic < raw[rwin].bic) rwin = i;
                int rsec = -1;
                for (int i = 0; i < 9; ++i) {
                    if (i == rwin) continue;
                    if (rsec < 0 || raw[i].bic < raw[rsec].bic) rsec = i;
                }
                profile.library_bic_raw_winner_model = raw[rwin].name;
                profile.library_bic_raw_winner_class = raw[rwin].klass;
                profile.library_bic_raw_winner_in_cascade = raw[rwin].in_cascade;
                if (rsec >= 0) {
                    profile.library_bic_raw_second_model = raw[rsec].name;
                    profile.library_bic_raw_margin       = raw[rsec].bic - raw[rwin].bic;
                }
                // Cascade-exclusion booleans: which gating reasons could have
                // suppressed the raw winner from cascade competition. Multi-valued.
                profile.library_bic_excl_in_cascade            = raw[rwin].in_cascade;
                profile.library_bic_excl_M_SS_full_hardcoded   = (std::string(raw[rwin].name) == "M_SS_full");
                profile.library_bic_excl_ct3_zero              = (std::string(raw[rwin].name) == "M_SS_orig"
                                                                  && ct3.delta_bic <= 0.0);
                profile.library_bic_excl_spike_is_ss           = ((std::string(raw[rwin].name) == "M_SS_asym" && !spike_is_ss)
                                                                  || (std::string(raw[rwin].name) == "M_DS_spike" && spike_is_ss));
                profile.library_bic_excl_structural_bilateral  = (std::string(raw[rwin].name) == "M_DS_asym_art"
                                                                  && !structural_bilateral);
            }
            // Low-confidence override: if the posterior winner sits below the
            // confidence threshold AND no domain rescue rule fired, fall back
            // to UNKNOWN so downstream consumers (fqdup) use a neutral damage
            // prior instead of committing to a marginal call.
            if (!profile.library_type_rescued &&
                profile.library_p_winner > 0.0f &&
                profile.library_p_winner <
                    SampleDamageProfile::kLibraryTypeConfidenceThreshold) {
                apply_override("low_confidence_to_unknown",
                               SampleDamageProfile::LibraryType::UNKNOWN, false);
            }
        } else {
            profile.library_bic_bias = 0.0;
            profile.library_bic_ds   = 0.0;
            profile.library_bic_ss   = 0.0;
            profile.library_bic_mix  = 0.0;
            profile.library_p_bias   = 0.0f;
            profile.library_p_ds     = 0.0f;
            profile.library_p_ss     = 0.0f;
            profile.library_p_winner = 0.0f;
            profile.library_type_evaluable = false;
            profile.library_type = SampleDamageProfile::LibraryType::UNKNOWN;
        }

        // Rescue rule: when BIC classification returns UNKNOWN due to known model
        // misspecification (adapter artifact / composition bias) but empirical CT5 signal
        // clearly shows ancient damage and SS-specific channels (GA0, CT3) are flat,
        // classify as DS. Uses raw per-position rates (not BIC-fitted amplitudes) to
        // avoid circularity. Effect-size thresholds (not z-scores) for depth robustness.
        if (profile.library_type == SampleDamageProfile::LibraryType::UNKNOWN) {
            const bool misspec = profile.composition_bias_5prime
                              || profile.position_0_artifact_5prime
                              || (profile.fit_offset_5prime > 1 && profile.position_0_artifact_5prime);
            // BUG FIX: rescue ran before d_max_5prime was assigned (that
            // happens further below at the raw_d_max_5prime block). Use
            // max_damage_5prime which is set up-stream at line ~1994.
            if (misspec && profile.max_damage_5prime >= 0.03f) {
                // CT5 empirical excess: max T/(T+C) at pos 1-4 over baseline
                float ct5_exc = 0.0f;
                float n_ct5   = 1.0f;
                for (int p = 1; p <= 4; ++p) {
                    double ntc = profile.tc_total_5prime[p];
                    if (ntc < 100.0) continue;
                    float exc = static_cast<float>(profile.t_freq_5prime[p])  // already T/(T+C) fraction
                              - static_cast<float>(ctx.baseline_tc);
                    if (exc > ct5_exc) { ct5_exc = exc; n_ct5 = static_cast<float>(ntc); }
                }
                // GA0 empirical excess at 3' pos-0 (SS complement-orientation indicator)
                float ga0_exc = 0.0f;
                if (profile.ag_total_3prime[0] >= 100.0) {
                    float ga0_rate = static_cast<float>(
                        profile.a_freq_3prime[0] / profile.ag_total_3prime[0]);
                    ga0_exc = std::max(0.0f, ga0_rate - static_cast<float>(ctx.baseline_ag));
                }
                // CT3 empirical excess at 3' pos 1-4 (SS original-orientation indicator)
                float ct3_exc = 0.0f;
                for (int p = 1; p <= 4; ++p) {
                    double ntc3 = profile.tc_total_3prime[p];
                    if (ntc3 < 100.0) continue;
                    float exc = static_cast<float>(ctx.ctrl_freq_3p[p])
                              - static_cast<float>(ctx.baseline_tc);
                    if (exc > ct3_exc) ct3_exc = exc;
                }
                // Overdispersed lower-95 CI for ct5_exc: SE × 2× inflation factor.
                // Avoids z-score dependence (which is depth-driven at 10M+ reads).
                float p_hat    = ct5_exc + static_cast<float>(ctx.baseline_tc);
                float se_od    = std::sqrt(p_hat * (1.0f - p_hat) / n_ct5) * 2.0f;
                float lower95  = ct5_exc - 1.96f * se_od;

                if (ct5_exc  >= 0.025f
                    && lower95  >= 0.01f
                    && ga0_exc  <= 0.01f
                    && ct3_exc  <= 0.01f
                    && std::max(ga0_exc, ct3_exc) <= 0.4f * ct5_exc) {
                    apply_override("misspec_unknown_to_ds",
                                   SampleDamageProfile::LibraryType::DOUBLE_STRANDED, true);
                }
            }
        }

        profile.library_type_auto_detected = true;

        // S1 telemetry: final winner after all post-hoc rescues / vetoes / UNKNOWN
        // overrides. final_library_bic_winner_model echoes the cascade-tournament
        // winner; override_reason is non-empty when post-hoc logic moved
        // library_type to a class that disagrees with the winner model's class.
        {
            const std::string& cw = profile.library_bic_winner_model;
            std::string cw_class;
            if (cw == "M_bias")                                    cw_class = "bias";
            else if (cw == "M_DS_symm" || cw == "M_DS_symm_art")   cw_class = "ds";
            else if (cw == "M_DS_spike")                           cw_class = profile.library_spike_is_ss ? "ss" : "ds";
            else                                                   cw_class = "ss";  // M_SS_*
            std::string lt_class;
            switch (profile.library_type) {
                case SampleDamageProfile::LibraryType::DOUBLE_STRANDED: lt_class = "ds"; break;
                case SampleDamageProfile::LibraryType::SINGLE_STRANDED: lt_class = "ss"; break;
                default:                                                lt_class = "unknown"; break;
            }
            profile.final_library_bic_winner_model = cw;
            // C5: surface the winner-model class at the top level so a
            // library_bic_winner_model vs library_type contradiction (after a veto)
            // is machine-visible without parsing the nested bic section.
            profile.library_bic_winner_model_class = cw_class;
            if (cw_class == lt_class) {
                profile.final_library_bic_override_reason.clear();
            } else if (lt_class == "unknown") {
                profile.final_library_bic_override_reason =
                    profile.library_type_rescued
                        ? "post_hoc_rescue_to_unknown"
                        : "low_confidence_override_to_unknown";
            } else if (profile.library_type_rescued) {
                profile.final_library_bic_override_reason =
                    "post_hoc_rescue_" + cw_class + "_to_" + lt_class;
            } else {
                profile.final_library_bic_override_reason =
                    "post_hoc_veto_" + cw_class + "_to_" + lt_class;
            }
            // C5: top-level mirror of the override reason so the veto signal is
            // co-located with library_bic_winner_model / library_type.
            profile.library_bic_override_reason = profile.final_library_bic_override_reason;
            // C5: raw library_p_ds/p_ss/p_bias are pre-override posteriors; flag when an
            // override fired so consumers read library_p_*_final, not the raw probs, as
            // the type confidence.
            profile.library_p_is_pre_override =
                !profile.final_library_bic_override_reason.empty();
            // F3: post-veto final probabilities. One-hot when override fired
            // (veto/rescue/UNKNOWN), mirror raw probs when class survived.
            if (profile.final_library_bic_override_reason.empty()) {
                profile.library_p_ds_final     = profile.library_p_ds;
                profile.library_p_ss_final     = profile.library_p_ss;
                profile.library_p_bias_final   = profile.library_p_bias;
                profile.library_p_winner_final = profile.library_p_winner;
            } else {
                profile.library_p_ds_final     = (lt_class == "ds")   ? 1.0f : 0.0f;
                profile.library_p_ss_final     = (lt_class == "ss")   ? 1.0f : 0.0f;
                profile.library_p_bias_final   = (lt_class == "bias") ? 1.0f : 0.0f;
                profile.library_p_winner_final = (lt_class == "unknown") ? 0.0f : 1.0f;
            }
        }

        // P0-1: capture the BIC tournament's verdict regardless of any forced override.
        profile.library_auto_type      = profile.library_type;
        profile.library_auto_evaluable = profile.library_type_evaluable;

        // P0-1: forced-library override — applied AFTER the tournament so all
        // diagnostic BICs / posteriors remain populated. Only library_type and
        // the auto-detect flag change; library_auto_type preserves the auto call.
        profile.library_forced_type = profile.forced_library_type;
        if (profile.forced_library_type != SampleDamageProfile::LibraryType::UNKNOWN) {
            profile.library_type               = profile.forced_library_type;
            profile.library_type_auto_detected = false;
        }
    }

    float damage_signal = (profile.max_damage_5prime + profile.max_damage_3prime) / 2.0f;

    float hexamer_boost = 0.0f;
    if (profile.hexamer_damage_llr > 0.02f && !profile.terminal_inversion) {
        // Normal sample with clear hexamer signal
        hexamer_boost = profile.hexamer_damage_llr * 8.0f;
    } else if (profile.terminal_inversion) {
        // Inverted sample: z-score asymmetry distinguishes real damage from composition bias.
        // Real damage: 3' G→A relatively stronger (z_ratio < 1.2).
        // Composition bias: 5' T dominates (z_ratio > 1.5).
        float z5_abs = std::abs(profile.terminal_z_5prime);
        float z3_abs = std::abs(profile.terminal_z_3prime);
        float z_ratio = (z3_abs > 0) ? z5_abs / z3_abs : 10.0f;

        // Use absolute hexamer LLR magnitude for inverted samples
        float abs_llr = std::abs(profile.hexamer_damage_llr);

        if (z_ratio < 1.2f && abs_llr > 0.02f) {
            hexamer_boost = abs_llr * 8.0f;
        } else if (z_ratio > 1.5f) {
            // 5' signal dominates → likely composition bias
        } else {
            hexamer_boost = abs_llr * 4.0f;
        }

        if (hexamer_boost > 0.01f && damage_signal < 0.01f) {
            damage_signal = hexamer_boost;
            profile.max_damage_5prime = hexamer_boost;
            profile.max_damage_3prime = hexamer_boost;
        }
    }

    float cpg_boost = 0.0f;
    float wobble_boost = 0.0f;

    if (damage_signal > 0.01f) {
        if (profile.cpg_ct_fraction > profile.non_cpg_ct_fraction + 0.05f) {
            cpg_boost = 0.03f;
        }

        float wobble_enrichment = profile.codon_pos_t_fraction_5prime[2] -
                                  (profile.codon_pos_t_fraction_5prime[0] + profile.codon_pos_t_fraction_5prime[1]) / 2.0f;
        if (wobble_enrichment > 0.05f) {
            wobble_boost = 0.02f;
        }
    }

    float total_signal = damage_signal + cpg_boost + wobble_boost + hexamer_boost;

    if (total_signal > 0.12f) {
        profile.sample_damage_prob = 0.95f;
    } else if (total_signal > 0.06f) {
        profile.sample_damage_prob = 0.80f;
    } else if (total_signal > 0.03f) {
        profile.sample_damage_prob = 0.50f;
    } else {
        profile.sample_damage_prob = 0.20f;
    }


}

} // namespace taph
