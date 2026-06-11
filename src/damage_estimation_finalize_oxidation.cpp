#include "damage_estimation_finalize_helpers.hpp"
#include "damage_estimation_finalize_ctx.hpp"
namespace taph {

void finalize_oxidation(SampleDamageProfile& profile, const FinalCtx& ctx) {
    // Channel C: oxidative stop codon analysis.
    // Unlike deamination, oxidation is uniform across reads (terminal ≈ interior rate).
    {
        double ox_pre_interior = profile.convertible_gag_interior +
                                 profile.convertible_gaa_interior +
                                 profile.convertible_gga_interior;
        double ox_stop_interior = profile.convertible_tag_ox_interior +
                                  profile.convertible_taa_ox_interior +
                                  profile.convertible_tga_ox_interior;
        double ox_total_interior = ox_pre_interior + ox_stop_interior;

        if (ox_total_interior > 100) {
            profile.ox_stop_conversion_rate_baseline = static_cast<float>(
                ox_stop_interior / ox_total_interior);
            profile.channel_c_valid = true;
        }

        double ox_pre_terminal = 0.0, ox_stop_terminal = 0.0;
        double ox_pre_mid = 0.0, ox_stop_mid = 0.0;

        for (int p = 0; p < 5; ++p) {
            ox_pre_terminal += profile.convertible_gag_5prime[p] +
                              profile.convertible_gaa_5prime[p] +
                              profile.convertible_gga_5prime[p];
            ox_stop_terminal += profile.convertible_tag_ox_5prime[p] +
                               profile.convertible_taa_ox_5prime[p] +
                               profile.convertible_tga_ox_5prime[p];
        }
        for (int p = 5; p < 15; ++p) {
            ox_pre_mid += profile.convertible_gag_5prime[p] +
                         profile.convertible_gaa_5prime[p] +
                         profile.convertible_gga_5prime[p];
            ox_stop_mid += profile.convertible_tag_ox_5prime[p] +
                          profile.convertible_taa_ox_5prime[p] +
                          profile.convertible_tga_ox_5prime[p];
        }

        double ox_total_terminal = ox_pre_terminal + ox_stop_terminal;
        double ox_total_mid = ox_pre_mid + ox_stop_mid;

        profile.ox_uniformity_ratio = 1.0f;
        if (ox_total_terminal > 50 && ox_total_mid > 50) {
            // C1: positional terminal/interior rates actually measured here.
            profile.ox_stop_rate_positional_computed = true;
            profile.ox_stop_rate_terminal = static_cast<float>(ox_stop_terminal / ox_total_terminal);
            profile.ox_stop_rate_interior = static_cast<float>(ox_stop_mid / ox_total_mid);

            // Uniformity ratio: terminal/interior (≈1 for uniform oxidation)
            // F5: raise gate 0.001f -> 0.01f; add log-ratio field.
            if (profile.ox_stop_rate_interior > 0.01f) {
                profile.ox_uniformity_ratio = profile.ox_stop_rate_terminal / profile.ox_stop_rate_interior;
                profile.log_ox_uniformity_ratio =
                    std::log((profile.ox_stop_rate_terminal + 1e-9f) / (profile.ox_stop_rate_interior + 1e-9f));
                // C1: 1.0f is also a genuine uniform value — flag the real measurement
                // so the emitter can distinguish it from the not-computed default.
                profile.ox_uniformity_ratio_computed = true;
            }
        }

        float ox_stop_excess = profile.ox_stop_rate_terminal - profile.ox_stop_conversion_rate_baseline;

        // Use stricter threshold if deamination is present (correlated damage)
        // Deamination enriches ALL terminal damage, including G→T via adjacent effects
        float threshold = 0.02f;  // 2% excess required (was 0.5%)
        // C5: d_max_combined is computed ~2000 lines later (still 0.0 here), so the
        // original guard was dead code. joint_delta_max (the joint-model MLE damage
        // rate, set above at the joint fit) is the live pre-computed deamination
        // level available at this point. Behavioral: the stricter 5% threshold now
        // actually engages when joint deamination >5% (previously only via damage_validated).
        if (profile.joint_delta_max > 0.05f || profile.damage_validated) {
            threshold = 0.05f;  // 5% excess if deamination present
        }

        bool elevated = ox_stop_excess > threshold;
        // Use log-ratio bands: scale-invariant, symmetric around 0, gate-independent.
        // log(0.85)≈-0.163  log(1.15)≈0.140  log(0.7)≈-0.357  log(1.5)≈0.405
        // _computed flag required so the 0.0f default (= log(1.0)) never spuriously fires.
        bool uniform = profile.ox_uniformity_ratio_computed
                       && profile.log_ox_uniformity_ratio > std::log(0.85f)
                       && profile.log_ox_uniformity_ratio < std::log(1.15f);

        if (profile.channel_c_valid && elevated && uniform) {
            // C5: codon-block verdict is coupled to ox_d_max. Record it in the
            // dedicated *_codon field; ox_damage_detected is the legacy field that
            // the GT model-fit block (below, ~line 2615) authoritatively overwrites.
            profile.ox_damage_detected_codon = true;
            profile.ox_damage_detected = true;
            profile.ox_d_max = std::max(0.0f, ox_stop_excess);  // fraction [0,1], consistent with other d_max fields
        }

        // Terminal >> interior: deamination cross-contamination, not true G→T oxidation.
        if (profile.ox_uniformity_ratio_computed
                && profile.log_ox_uniformity_ratio > std::log(1.5f)) {
            profile.ox_is_artifact = true;
            profile.ox_damage_detected_codon = false;
            profile.ox_damage_detected = false;
        }

        // Interior >> terminal: anomalous suppression at terminus — also artifact.
        if (profile.ox_uniformity_ratio_computed
                && profile.log_ox_uniformity_ratio < std::log(0.7f)) {
            profile.ox_is_artifact = true;
            profile.ox_damage_detected_codon = false;
            profile.ox_damage_detected = false;
        }
    }

    // Channel C3': oxidative G→T stop codon analysis at 3' end.
    {
        double ox_pre_3p = 0, ox_stop_3p = 0;
        double ox_pre_term3 = 0, ox_stop_term3 = 0;
        double ox_pre_mid3  = 0, ox_stop_mid3  = 0;
        for (int p = 0; p < 15; ++p) {
            double pre  = profile.convertible_gag_3prime[p] +
                          profile.convertible_gaa_3prime[p] +
                          profile.convertible_gga_3prime[p];
            double stop = profile.convertible_tag_ox_3prime[p] +
                          profile.convertible_taa_ox_3prime[p] +
                          profile.convertible_tga_ox_3prime[p];
            ox_pre_3p  += pre;  ox_stop_3p  += stop;
            if (p < 5)  { ox_pre_term3 += pre; ox_stop_term3 += stop; }
            else        { ox_pre_mid3  += pre; ox_stop_mid3  += stop; }
        }
        double total_3p = ox_pre_3p + ox_stop_3p;
        if (total_3p >= 200.0) {
            profile.channel_c3_valid = true;
            profile.ox_stop_baseline_3prime = static_cast<float>(ox_stop_3p / total_3p);
        }
        if (ox_pre_term3 + ox_stop_term3 > 50)
            profile.ox_stop_rate_terminal_3prime = static_cast<float>(
                ox_stop_term3 / (ox_pre_term3 + ox_stop_term3));
        if (ox_pre_mid3 + ox_stop_mid3 > 50) {
            profile.ox_stop_rate_interior_3prime = static_cast<float>(
                ox_stop_mid3 / (ox_pre_mid3 + ox_stop_mid3));
            if (profile.ox_stop_rate_interior_3prime > 0.01f && profile.ox_stop_rate_terminal_3prime > 0.0f) {
                profile.ox_uniformity_ratio_3prime =
                    profile.ox_stop_rate_terminal_3prime / profile.ox_stop_rate_interior_3prime;
                profile.log_ox_uniformity_ratio_3prime =
                    std::log((profile.ox_stop_rate_terminal_3prime + 1e-9f) / (profile.ox_stop_rate_interior_3prime + 1e-9f));
                profile.ox_uniformity_ratio_3prime_computed = true;  // C1: real measurement
            }
        }
    }

    // Channel D: Chargaff-balance G→T / C→A oxidation estimate.
    // Reference-free: under no damage, T/(T+G) ≈ A/(A+C) (Chargaff) and terminal ≈ interior.
    // 8-oxoG converts G→T uniformly, raising T/(T+G) above A/(A+C) throughout reads.
    {
        // Interior baseline T/(T+G) and A/(A+C)
        double gt_base_total = profile.baseline_g_total + profile.baseline_g_to_t_count;
        double ca_base_total = profile.baseline_c_ox_total + profile.baseline_c_to_a_count;
        if (gt_base_total >= 500.0) {
            profile.ox_gt_baseline = static_cast<float>(profile.baseline_g_to_t_count / gt_base_total);
        }
        if (ca_base_total >= 500.0) {
            profile.ox_ca_baseline = static_cast<float>(profile.baseline_c_to_a_count / ca_base_total);
        }

        // T/(T+G) at positions 0-4 (terminal) and 5-14 (mid-read)
        double t_term = 0, g_term = 0, t_mid = 0, g_mid = 0;
        for (int p = 0; p < 5; ++p)  { t_term += profile.t_from_g_5prime[p]; g_term += profile.g_count_5prime[p]; }
        for (int p = 5; p < 15; ++p) { t_mid  += profile.t_from_g_5prime[p]; g_mid  += profile.g_count_5prime[p]; }
        if (t_term + g_term >= 200.0) profile.ox_gt_rate_terminal = static_cast<float>(t_term / (t_term + g_term));
        if (t_mid  + g_mid  >= 200.0) {
            profile.ox_gt_rate_interior = static_cast<float>(t_mid / (t_mid + g_mid));
            if (profile.ox_gt_rate_interior > 0.01f && profile.ox_gt_rate_terminal > 0.0f) {
                profile.ox_gt_uniformity = profile.ox_gt_rate_terminal / profile.ox_gt_rate_interior;
                profile.log_ox_gt_uniformity =
                    std::log((profile.ox_gt_rate_terminal + 1e-9f) / (profile.ox_gt_rate_interior + 1e-9f));
                profile.ox_gt_uniformity_computed = true;  // C1: real measurement vs 0.0f default
            }
        }

        // A/(A+C) at positions 0-4 (terminal) and 5-14 (mid-read)
        double a_term = 0, c_term = 0, a_mid = 0, c_mid = 0;
        for (int p = 0; p < 5; ++p)  { a_term += profile.a_from_c_5prime[p]; c_term += profile.c_count_ox_5prime[p]; }
        for (int p = 5; p < 15; ++p) { a_mid  += profile.a_from_c_5prime[p]; c_mid  += profile.c_count_ox_5prime[p]; }
        if (a_term + c_term >= 200.0) profile.ox_ca_rate_terminal = static_cast<float>(a_term / (a_term + c_term));
        if (a_mid  + c_mid  >= 200.0) {
            profile.ox_ca_rate_interior = static_cast<float>(a_mid / (a_mid + c_mid));
            if (profile.ox_ca_rate_interior > 0.01f && profile.ox_ca_rate_terminal > 0.0f) {
                profile.ox_ca_uniformity = profile.ox_ca_rate_terminal / profile.ox_ca_rate_interior;
                profile.log_ox_ca_uniformity =
                    std::log((profile.ox_ca_rate_terminal + 1e-9f) / (profile.ox_ca_rate_interior + 1e-9f));
                profile.ox_ca_uniformity_computed = true;
            }
        }

        // Chargaff asymmetry: excess T/(T+G) over A/(A+C) at interior (= Chargaff deviation)
        // C1: d_computed certifies D / tg_interior / ac_interior / s_gt together;
        // both interior baselines must reach coverage for the contrast to be real.
        if (gt_base_total >= 500.0 && ca_base_total >= 500.0) {
            profile.ox_gt_asymmetry = profile.ox_gt_baseline - profile.ox_ca_baseline;
            profile.d_computed = true;
        }
    }

    // Pre-compute hexamer_excess_tc early so Channel F/G/H adapter skip can use it.
    // Full hexamer analysis runs later; this is just the TC-ratio summary.
    if (profile.hexamer_excess_tc == 0.0f &&
        profile.n_hexamers_5prime >= 1000 && profile.n_hexamers_interior >= 1000) {
        double tt = 0, tc = 0, it = 0, ic = 0;
        for (uint32_t b = 0; b < 1024; ++b) {
            uint32_t c_hex = 0x400 | b, t_hex = 0xC00 | b;
            tc += profile.hexamer_count_5prime[c_hex];
            tt += profile.hexamer_count_5prime[t_hex];
            ic += profile.hexamer_count_interior[c_hex];
            it += profile.hexamer_count_interior[t_hex];
        }
        profile.hexamer_excess_tc = static_cast<float>(
            tt / (tt + tc + 1e-10) - it / (it + ic + 1e-10));
    }

    // Binomial z-score: terminal stop rate vs resolved interior reference.
    // Caller must resolve far-interior vs mid-read fallback before calling.
    // C1: returns NaN (not 0.0f) when the computation cannot run, so 'not computed'
    // is distinguishable from a genuine z==0.
    // C2: the emitted magnitude is clamped to ±kZCap. This z treats correlated
    // reads (PCR dups, admixture, composition) as independent Bernoulli trials, so
    // it scales ~sqrt(N) and is exploratory, NOT a calibrated p-value. Detection
    // gates use z>3 (well inside ±12), so clamping does not change any decision;
    // it only stops absurd hundreds/thousands magnitudes reaching consumers.
    auto binom_z = [](double k_t, double n_t, double k_i, double n_i) -> float {
        return binom_z_clamped(k_t, n_t, k_i, n_i);
    };

    // Shared finalization: fills baseline, terminal rate, interior rate, uniformity.
    // pre_i/stop_i is the already-resolved interior reference (far or mid-read).
    // C3: validity must certify exactly the metrics a consumer reads under it.
    // valid is now atomic: true only when the baseline (total>=200) AND the
    // terminal window (pre_t+stop_t>50) AND the interior reference (pre_i+stop_i>50)
    // are all covered, so term_r / int_r / unif can never stay at their 0.0f
    // default while valid==true.
    auto fin_oxog = [](double pre_tot, double stop_tot,
                       double pre_t,   double stop_t,
                       double pre_i,   double stop_i,
                       float& baseline, float& term_r, float& int_r, float& unif, bool& valid) {
        double total = pre_tot + stop_tot;
        bool have_base = (total >= 200.0);
        bool have_term = (pre_t + stop_t > 50);
        bool have_int  = (pre_i + stop_i > 50);
        if (have_base) baseline = static_cast<float>(stop_tot / total);
        if (have_term) term_r = static_cast<float>(stop_t  / (pre_t  + stop_t));
        if (have_int) {
            int_r = static_cast<float>(stop_i / (pre_i + stop_i));
            if (int_r > 1e-6f && term_r > 0.0f) unif = term_r / int_r;
        }
        valid = have_base && have_term && have_int;
    };

    // Layer-2 producer: one function turns a stop channel's spec + accumulated terminal/interior
    // counts into BOTH the Layer-0 count table AND the emitted z / rates / validity. The emitted z
    // is taken from the count table's own clamped z, so every emitted stop-channel metric traces to
    // exactly one count-table row — there is no longer a separate binom_z call that could silently
    // diverge from the table the gate checks.
    auto compute_stop_channel = [&](const ChannelSpec& spec,
                                    double pre_tot, double stop_tot,
                                    double pre_t, double stop_t, double pre_i, double stop_i,
                                    double shadow_t, double shadow_i,
                                    float& out_z, float& out_or, float& baseline, float& term_r,
                                    float& int_r, float& unif, bool& valid) -> StopChannelCountTable {
        StopChannelCountTable ct = make_stop_count_table(spec, pre_t, stop_t, pre_i, stop_i, shadow_t, shadow_i);
        out_z = static_cast<float>(ct.post_clamp_z);
        // Effect size is selected by the registry: single-stratum channels (G/H) report the 2x2
        // Haldane-Anscombe OR here; F's effect size is the Mantel-Haenszel common OR set separately.
        out_or = (spec.variance == VarianceFamily::PLAIN_2x2_OR)
                 ? static_cast<float>(stop_channel_or_haldane(ct))
                 : std::numeric_limits<float>::quiet_NaN();
        fin_oxog(pre_tot, stop_tot, pre_t, stop_t, pre_i, stop_i, baseline, term_r, int_r, unif, valid);
        return ct;
    };

    // Channels F, G, H: resolve interior reference (far-interior pos 30+ if ≥50 counts,
    // else mid-read pos 5-14), then compute z and fill profile fields.
    // Adapter skip: set to 1 to exclude p=0 from terminal accumulation when a
    // position-0 composition artifact is flagged at either end.
    // TC hexamer bias was previously also gated here for channels F and G, but
    // empirical validation showed TC hexamer drives F/G z-scores *negative* across
    // the full terminal window — excluding p=0 does not prevent this, and the gate
    // removes genuine signal in ancient TC-hexamer libraries. Gate removed.
    const int p0_tc_5 = profile.position_0_artifact_5prime ? 1 : 0;
    const int p0_h5   = profile.position_0_artifact_5prime ? 1 : 0;
    const int p0_3    = profile.position_0_artifact_3prime ? 1 : 0;

    // Channel F: C→A oxidation (bottom-strand 8-oxoG)
    {
        double pre5 = 0, stop5 = 0, pre_t = 0, stop_t = 0, pre_m = 0, stop_m = 0;
        double shadow_t = 0, shadow_m = 0;
        for (int p = 0; p < 15; ++p) {
            double pre    = profile.convertible_tca_5prime[p] + profile.convertible_tcg_5prime[p] +
                            profile.convertible_tac_5prime[p] + profile.convertible_tgc_5prime[p];
            double stop   = profile.convertible_taa_ca_5prime[p] + profile.convertible_tag_ca_5prime[p] +
                            profile.convertible_tga_ca_5prime[p];
            double shadow = profile.ca_deam_shadow_5prime[p];
            pre5 += pre; stop5 += stop;
            if (p >= p0_tc_5 && p < 5) { pre_t += pre; stop_t += stop; shadow_t += shadow; }
            else if (p >= 5)           { pre_m += pre; stop_m += stop; shadow_m += shadow; }
        }
        double pre_far    = (double)profile.ca_pre_interior;
        double stop_far   = (double)profile.ca_stop_interior;
        double shadow_far = (double)profile.ca_deam_shadow_interior;
        bool use_far = (pre_far + stop_far >= 50);
        double pre_i    = use_far ? pre_far    : pre_m;
        double stop_i   = use_far ? stop_far   : stop_m;
        double shadow_i = use_far ? shadow_far : shadow_m;
        float f_or_discard = std::numeric_limits<float>::quiet_NaN();  // F's effect size is the MH common OR, set below
        StopChannelCountTable f_ct = compute_stop_channel(stop_channel_spec('F'),
                 pre5, stop5, pre_t, stop_t, pre_i, stop_i, shadow_t, shadow_i,
                 profile.channel_f_z, f_or_discard, profile.ca_stop_rate_baseline, profile.ca_stop_rate_terminal,
                 profile.ca_stop_rate_interior, profile.ca_uniformity_ratio, profile.channel_f_valid);

        // Mantel-Haenszel stratified test: 3 context strata (TCA+TAC, TCG, TGC)
        // Removes terminal context-composition bias that inflates negative z.
        // Robins-Breslow-Greenland variance of log(MH_OR).
        {
            const double* pre5_arr[3]  = { profile.convertible_tca_5prime.data(), profile.convertible_tcg_5prime.data(), profile.convertible_tgc_5prime.data() };
            const double* pre5b_arr[3] = { profile.convertible_tac_5prime.data(), nullptr, nullptr };
            const double* stop5_arr[3] = { profile.convertible_taa_ca_5prime.data(), profile.convertible_tag_ca_5prime.data(), profile.convertible_tga_ca_5prime.data() };
            const double* sha5_arr[3]  = { profile.ca_shadow_5prime_ctx0.data(), profile.ca_shadow_5prime_ctx1.data(), profile.ca_shadow_5prime_ctx2.data() };
            double mh_R = 0.0, mh_S = 0.0;
            double rbg_PR = 0.0, rbg_QS = 0.0, rbg_PQRS = 0.0;
            for (int k = 0; k < stop_channel_spec('F').n_strata; ++k) {  // registry-driven stratum count (F: 3)
                double a = 0.0, b_val = 0.0;
                for (int p = p0_tc_5; p < 5; ++p) {
                    double pre_k = pre5_arr[k][p] + (pre5b_arr[k] ? pre5b_arr[k][p] : 0.0);
                    b_val += pre_k + sha5_arr[k][p];
                    a     += stop5_arr[k][p];
                }
                double c = profile.ca_stop_interior_by_ctx[k];
                double d = profile.ca_pre_interior_by_ctx[k] + profile.ca_shadow_interior_by_ctx[k];
                f_ct.strata[k] = StopChannelStratum{k, a, b_val, c, d};
                double N = a + b_val + c + d;
                if (N < 8.0) continue;
                double R = a * d / N;
                double S = b_val * c / N;
                mh_R += R;
                mh_S += S;
                double P = (a + d) / N;
                double Q = (b_val + c) / N;
                rbg_PR   += P * R;
                rbg_QS   += Q * S;
                rbg_PQRS += P * S + Q * R;
            }
            if (mh_R > 0.0 && mh_S > 0.0) {
                double common_or  = mh_R / mh_S;
                double log_or     = std::log(common_or);
                double var_log_or = rbg_PR   / (2.0 * mh_R * mh_R)
                                  + rbg_PQRS / (2.0 * mh_R * mh_S)
                                  + rbg_QS   / (2.0 * mh_S * mh_S);
                if (var_log_or > 1e-15) {
                    // C2: RBG variance sums O(N) terms, so this Mantel-Haenszel z is
                    // explosive (|z|≈2000–9000 on the panel) and NOT a calibrated
                    // p-value. Clamp the emitted magnitude to ±kZCap; channel_f_common_or
                    // (the common odds ratio, set below) is the interpretable effect size.
                    double mh_z = log_or / std::sqrt(var_log_or);
                    const double cap = static_cast<double>(SampleDamageProfile::kZCap);
                    profile.channel_f_mh_z = static_cast<float>(std::clamp(mh_z, -cap, cap));
                }
                profile.channel_f_common_or = static_cast<float>(common_or);
            }
        }
        profile.count_tables.push_back(std::move(f_ct));

        double pre3 = 0, stop3 = 0, pre_t3 = 0, stop_t3 = 0, pre_m3 = 0, stop_m3 = 0;
        for (int p = 0; p < 15; ++p) {
            double pre  = profile.convertible_tca_3prime[p] + profile.convertible_tcg_3prime[p] +
                          profile.convertible_tac_3prime[p] + profile.convertible_tgc_3prime[p];
            double stop = profile.convertible_taa_ca_3prime[p] + profile.convertible_tag_ca_3prime[p] +
                          profile.convertible_tga_ca_3prime[p];
            pre3 += pre; stop3 += stop;
            if (p >= p0_3 && p < 5) { pre_t3 += pre; stop_t3 += stop; }
            else if (p >= 5)        { pre_m3 += pre; stop_m3 += stop; }
        }
        fin_oxog(pre3, stop3, pre_t3, stop_t3, pre_m3, stop_m3,
                 profile.ca_stop_rate_baseline_3prime, profile.ca_stop_rate_terminal_3prime,
                 profile.ca_stop_rate_interior_3prime, profile.ca_uniformity_ratio_3prime, profile.channel_f3_valid);
    }

    // Channel G: C->G terminal stop-enrichment (empirical; hydantoins are GUANINE over-oxidation
    // products, not earned by a complement-strand C->G stop observation)
    {
        double pre5 = 0, stop5 = 0, pre_t = 0, stop_t = 0, pre_m = 0, stop_m = 0;
        for (int p = 0; p < 15; ++p) {
            double pre  = profile.convertible_tca_cg_5prime[p] + profile.convertible_tac_cg_5prime[p];
            double stop = profile.convertible_tga_cg_5prime[p] + profile.convertible_tag_cg_5prime[p];
            pre5 += pre; stop5 += stop;
            if (p >= p0_tc_5 && p < 5) { pre_t += pre; stop_t += stop; }
            else if (p >= 5)           { pre_m += pre; stop_m += stop; }
        }
        double pre_far  = (double)profile.cg_pre_interior;
        double stop_far = (double)profile.cg_stop_interior;
        double pre_i  = (pre_far + stop_far >= 50) ? pre_far  : pre_m;
        double stop_i = (pre_far + stop_far >= 50) ? stop_far : stop_m;
        profile.count_tables.push_back(compute_stop_channel(stop_channel_spec('G'),
                 pre5, stop5, pre_t, stop_t, pre_i, stop_i, 0.0, 0.0,
                 profile.channel_g_z, profile.channel_g_or, profile.cg_stop_rate_baseline, profile.cg_stop_rate_terminal,
                 profile.cg_stop_rate_interior, profile.cg_uniformity_ratio, profile.channel_g_valid));

        // Correction 3 (Tier 2): context-stratified Mantel-Haenszel test for G, mirroring F's block.
        // 2 strata = the C→G pre-contexts (ctx0: TCA→TGA, ctx1: TAC→TAG). Removes terminal
        // context-composition bias. No deamination shadow (G is empirical; has_deam_shadow=false).
        {
            const double* pre_arr[2]  = { profile.convertible_tca_cg_5prime.data(), profile.convertible_tac_cg_5prime.data() };
            const double* stop_arr[2] = { profile.convertible_tga_cg_5prime.data(), profile.convertible_tag_cg_5prime.data() };
            double mh_R = 0.0, mh_S = 0.0, rbg_PR = 0.0, rbg_QS = 0.0, rbg_PQRS = 0.0;
            for (int k = 0; k < 2; ++k) {
                double a = 0.0, b_val = 0.0;
                for (int p = p0_tc_5; p < 5; ++p) { a += stop_arr[k][p]; b_val += pre_arr[k][p]; }
                double c = profile.cg_stop_interior_by_ctx[k];
                double d = profile.cg_pre_interior_by_ctx[k];
                double N = a + b_val + c + d;
                if (N < 8.0) continue;
                double R = a * d / N, S = b_val * c / N;
                mh_R += R; mh_S += S;
                double P = (a + d) / N, Q = (b_val + c) / N;
                rbg_PR   += P * R;
                rbg_QS   += Q * S;
                rbg_PQRS += P * S + Q * R;
            }
            if (mh_R > 0.0 && mh_S > 0.0) {
                double common_or = mh_R / mh_S, log_or = std::log(common_or);
                double var_log_or = rbg_PR / (2.0 * mh_R * mh_R) + rbg_PQRS / (2.0 * mh_R * mh_S)
                                  + rbg_QS / (2.0 * mh_S * mh_S);
                if (var_log_or > 1e-15) {
                    double mh_z = log_or / std::sqrt(var_log_or);
                    const double cap = static_cast<double>(SampleDamageProfile::kZCap);
                    profile.channel_g_mh_z = static_cast<float>(std::clamp(mh_z, -cap, cap));
                }
                profile.channel_g_common_or = static_cast<float>(common_or);
            }
        }

        double pre3 = 0, stop3 = 0, pre_t3 = 0, stop_t3 = 0, pre_m3 = 0, stop_m3 = 0;
        for (int p = 0; p < 15; ++p) {
            double pre  = profile.convertible_tca_cg_3prime[p] + profile.convertible_tac_cg_3prime[p];
            double stop = profile.convertible_tga_cg_3prime[p] + profile.convertible_tag_cg_3prime[p];
            pre3 += pre; stop3 += stop;
            if (p >= p0_3 && p < 5) { pre_t3 += pre; stop_t3 += stop; }
            else if (p >= 5)        { pre_m3 += pre; stop_m3 += stop; }
        }
        fin_oxog(pre3, stop3, pre_t3, stop_t3, pre_m3, stop_m3,
                 profile.cg_stop_rate_baseline_3prime, profile.cg_stop_rate_terminal_3prime,
                 profile.cg_stop_rate_interior_3prime, profile.cg_uniformity_ratio_3prime, profile.channel_g3_valid);
    }

    // Channel H: A->T terminal stop-enrichment (empirical; no established direct
    // adenine-oxidation->A->T pathway); z_p2plus excludes p=0 and p=1 as sensitivity check
    {
        double pre5 = 0, stop5 = 0, pre_t = 0, stop_t = 0, pre_m = 0, stop_m = 0;
        double pre_t2 = 0, stop_t2 = 0;
        for (int p = 0; p < 15; ++p) {
            double pre  = profile.convertible_aaa_h_5prime[p] + profile.convertible_aag_h_5prime[p] +
                          profile.convertible_aga_h_5prime[p];
            double stop = profile.convertible_taa_at_5prime[p] + profile.convertible_tag_at_5prime[p] +
                          profile.convertible_tga_at_5prime[p];
            pre5 += pre; stop5 += stop;
            if (p >= p0_h5 && p < 5) { pre_t += pre; stop_t += stop; }
            else if (p >= 5)        { pre_m += pre; stop_m += stop; }
            if (p >= 2 && p < 5)    { pre_t2 += pre; stop_t2 += stop; }
        }
        double pre_far  = (double)profile.at_pre_interior;
        double stop_far = (double)profile.at_stop_interior;
        double pre_i  = (pre_far + stop_far >= 50) ? pre_far  : pre_m;
        double stop_i = (pre_far + stop_far >= 50) ? stop_far : stop_m;
        StopChannelCountTable h_ct = compute_stop_channel(stop_channel_spec('H'),
                 pre5, stop5, pre_t, stop_t, pre_i, stop_i, 0.0, 0.0,
                 profile.channel_h_z, profile.channel_h_or, profile.at_stop_rate_baseline, profile.at_stop_rate_terminal,
                 profile.at_stop_rate_interior, profile.at_uniformity_ratio, profile.channel_h_valid);
        profile.channel_h_z_p2plus = binom_z(stop_t2, pre_t2 + stop_t2, stop_i, pre_i + stop_i);
        // C5: genuine A->T stop-enrichment should produce positive z at BOTH windows. When
        // h_z (incl. p0/p1) and h_z_p2plus (excl. p0/p1) disagree in sign, the signal
        // is artifact-driven at p0/p1. This flag surfaces that contradiction; the
        // detection gate (emitter) ANDs the two windows. NaN (not computed) -> false.
        profile.channel_h_z_consistent =
            std::isfinite(profile.channel_h_z) && std::isfinite(profile.channel_h_z_p2plus) &&
            ((profile.channel_h_z >= 0.0f) == (profile.channel_h_z_p2plus >= 0.0f));
        profile.count_tables.push_back(std::move(h_ct));

        // Correction 3 (Tier 2): context-stratified Mantel-Haenszel test for H, mirroring F's block.
        // 3 strata = the A→T pre-contexts (ctx0: AAA→TAA, ctx1: AAG→TAG, ctx2: AGA→TGA).
        {
            const double* pre_arr[3]  = { profile.convertible_aaa_h_5prime.data(), profile.convertible_aag_h_5prime.data(), profile.convertible_aga_h_5prime.data() };
            const double* stop_arr[3] = { profile.convertible_taa_at_5prime.data(), profile.convertible_tag_at_5prime.data(), profile.convertible_tga_at_5prime.data() };
            double mh_R = 0.0, mh_S = 0.0, rbg_PR = 0.0, rbg_QS = 0.0, rbg_PQRS = 0.0;
            for (int k = 0; k < 3; ++k) {
                double a = 0.0, b_val = 0.0;
                for (int p = p0_h5; p < 5; ++p) { a += stop_arr[k][p]; b_val += pre_arr[k][p]; }
                double c = profile.at_stop_interior_by_ctx[k];
                double d = profile.at_pre_interior_by_ctx[k];
                double N = a + b_val + c + d;
                if (N < 8.0) continue;
                double R = a * d / N, S = b_val * c / N;
                mh_R += R; mh_S += S;
                double P = (a + d) / N, Q = (b_val + c) / N;
                rbg_PR   += P * R;
                rbg_QS   += Q * S;
                rbg_PQRS += P * S + Q * R;
            }
            if (mh_R > 0.0 && mh_S > 0.0) {
                double common_or = mh_R / mh_S, log_or = std::log(common_or);
                double var_log_or = rbg_PR / (2.0 * mh_R * mh_R) + rbg_PQRS / (2.0 * mh_R * mh_S)
                                  + rbg_QS / (2.0 * mh_S * mh_S);
                if (var_log_or > 1e-15) {
                    double mh_z = log_or / std::sqrt(var_log_or);
                    const double cap = static_cast<double>(SampleDamageProfile::kZCap);
                    profile.channel_h_mh_z = static_cast<float>(std::clamp(mh_z, -cap, cap));
                }
                profile.channel_h_common_or = static_cast<float>(common_or);
            }
        }

        double pre3 = 0, stop3 = 0, pre_t3 = 0, stop_t3 = 0, pre_m3 = 0, stop_m3 = 0;
        for (int p = 0; p < 15; ++p) {
            double pre  = profile.convertible_aaa_h_3prime[p] + profile.convertible_aag_h_3prime[p] +
                          profile.convertible_aga_h_3prime[p];
            double stop = profile.convertible_taa_at_3prime[p] + profile.convertible_tag_at_3prime[p] +
                          profile.convertible_tga_at_3prime[p];
            pre3 += pre; stop3 += stop;
            if (p >= p0_3 && p < 5) { pre_t3 += pre; stop_t3 += stop; }
            else if (p >= 5)        { pre_m3 += pre; stop_m3 += stop; }
        }
        fin_oxog(pre3, stop3, pre_t3, stop_t3, pre_m3, stop_m3,
                 profile.at_stop_rate_baseline_3prime, profile.at_stop_rate_terminal_3prime,
                 profile.at_stop_rate_interior_3prime, profile.at_uniformity_ratio_3prime, profile.channel_h3_valid);
    }

    // GT exponential-background fit: GT(p) = A*exp(-mu*p) + B
    // B is the uniform background (8-oxoG rate); A is the terminal artifact.
    // Mirrors the C→T deamination model, extracting the position-independent component.
    // For SS: s_gt = B - ox_ca_baseline is a valid 8-oxoG signal (no Chargaff cancellation).
    // For DS: B and A/(A+C) both increase with 8-oxoG, so s_gt ≈ 0; report B directly.
    {
        float y[15] = {}, w[15] = {};
        float total_w = 0;
        for (int p = 0; p < 15; ++p) {
            double tot = profile.t_from_g_5prime[p] + profile.g_count_5prime[p];
            if (tot > 10.0) {
                w[p] = (float)tot;
                y[p] = (float)(profile.t_from_g_5prime[p] / tot);
            }
            total_w += w[p];
        }

        if (total_w >= 200.0f) {
            static constexpr float mu_grid[] = {0.05f, 0.1f, 0.2f, 0.3f, 0.5f, 0.7f, 1.0f, 1.5f, 2.0f, 3.0f};
            float best_sse = 1e30f, best_A = 0.0f, best_mu = 0.3f, best_B = 0.0f;
            float best_A_raw = 0.0f;
            float best_B_raw = 0.0f;

            for (float mu : mu_grid) {
                float sx2=0, sx=0, sxy=0, sy=0, sw=0;
                for (int p = 0; p < 15; ++p) {
                    if (w[p] == 0) continue;
                    float x = std::exp(-mu * p);
                    sx2 += w[p] * x * x;
                    sx  += w[p] * x;
                    sxy += w[p] * x * y[p];
                    sy  += w[p] * y[p];
                    sw  += w[p];
                }
                float det = sx2*sw - sx*sx;
                if (std::abs(det) < 1e-12f) continue;
                float A_raw = (sxy*sw - sy*sx) / det;
                float B_raw = (sx2*sy - sx*sxy) / det;
                float A = std::max(0.0f, A_raw);
                float B = std::max(0.0f, std::min(0.5f, B_raw));

                float sse = 0;
                for (int p = 0; p < 15; ++p) {
                    if (w[p] == 0) continue;
                    float res = y[p] - A * std::exp(-mu * p) - B;
                    sse += w[p] * res * res;
                }
                if (sse < best_sse) {
                    best_sse = sse;
                    best_A = A; best_mu = mu; best_B = B;
                    best_A_raw = A_raw;
                    best_B_raw = B_raw;
                }
            }

            profile.g_bg_fitted   = best_B;
            profile.g_term_fitted  = best_A;
            profile.g_decay_fitted = best_mu;
            if (profile.ox_ca_baseline > 0.0f)
                profile.s_gt = best_B - profile.ox_ca_baseline;

            // C4: decay rate selected at the grid maximum -> true mu unidentifiable.
            profile.gt_decay_at_upper_boundary = (best_mu == mu_grid[9]);
            // C4: background clamped at the 0.5 ceiling; keep unclamped point estimate.
            profile.gt_bg_at_upper_boundary = (best_B_raw > 0.5f);
            profile.g_bg_fitted_unclamped   = best_B_raw;

            // Interior-mean fallback (positions 5-14)
            {
                float ym = 0.0f, wm = 0.0f;
                for (int p = 5; p < 15; ++p) {
                    if (w[p] > 0) { ym += w[p] * y[p]; wm += w[p]; }
                }
                profile.g_bg_interior_mean = (wm > 0) ? ym / wm : 0.0f;
            }

            profile.g_fit_degenerate = (best_mu <= 0.10f);
            // C4: zero terminal amplitude (inverted/no terminal excess) makes the
            // decay rate meaningless -> null the decay and flag the fit degenerate.
            if (best_A == 0.0f) {
                profile.gt_term_zero_clamped = true;
                profile.g_fit_degenerate     = true;
                profile.g_decay_fitted       = 0.0f;
            }

            // C4: ox_theta (== g_bg_fitted, the 8-oxoG background B) pinned at the
            // 0.5 clamp while the fit is degenerate is a lower bound, not a point
            // estimate. Emitter flags it; g_bg_fitted_unclamped holds the raw magnitude.
            profile.ox_theta_at_clamp = (profile.g_bg_fitted >= 0.5f && profile.g_fit_degenerate);

            // WLS CI95 on B at fixed best_mu.
            // Linear model: y_p = A*x_p + B, x_p = exp(-best_mu*p).
            // Var(B) = sigma2 * sx2 / det  where det = sx2*sw - sx^2.
            {
                float sx2 = 0, sx = 0, sxy = 0, sy = 0, sw2 = 0;
                int   np = 0;
                for (int p = 0; p < 15; ++p) {
                    if (w[p] == 0) continue;
                    float x = std::exp(-best_mu * p);
                    sx2 += w[p] * x * x;
                    sx  += w[p] * x;
                    sxy += w[p] * x * y[p];
                    sy  += w[p] * y[p];
                    sw2 += w[p];
                    ++np;
                }
                float det = sx2 * sw2 - sx * sx;
                if (np > 2 && std::abs(det) > 1e-12f) {
                    float A_hat = (sxy * sw2 - sy * sx) / det;
                    float B_hat = (sx2 * sy - sx * sxy) / det;
                    float sse_ci = 0.0f;
                    for (int p = 0; p < 15; ++p) {
                        if (w[p] == 0) continue;
                        float res = y[p] - A_hat * std::exp(-best_mu * p) - B_hat;
                        sse_ci += w[p] * res * res;
                    }
                    float sigma2 = sse_ci / static_cast<float>(np - 2);
                    float var_B  = sigma2 * sx2 / det;
                    float se_B   = std::sqrt(std::max(0.0f, var_B));
                    // C4: center the CI on the unclamped WLS estimate B_hat (not the
                    // clamped best_B) and drop the 0.5 cap on the upper bound. Centering
                    // on best_B shifted the whole interval down by (B_hat-0.5) and hard-
                    // capped ci_hi at 0.5 whenever B was clamped, making it uninformative.
                    // T/(T+G) is a proportion in [0,1], so >0.5 is valid (GC-rich/oxidised).
                    profile.g_bg_fitted_ci_lo = std::max(0.0f, B_hat - 1.96f * se_B);
                    profile.g_bg_fitted_ci_hi = std::min(1.0f, B_hat + 1.96f * se_B);
                }
            }

            // Update detection: model-based uniform G→T signal replaces codon-based Channel C.
            // For SS: s_gt > threshold (Chargaff contrast valid because no complementary strand).
            // For DS: B elevated above ca_baseline indicates library-level oxidation.
            const bool is_ss_lib = (profile.library_type ==
                                    taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED);
            // s_gt = best_B - ox_ca_baseline (Chargaff contrast) is the valid 8-oxoG signal
            // for both SS and DS: under no OxoG it is ~0; under OxoG it is positive.
            // Using best_B directly for DS compared a raw frequency (~0.43) against a
            // near-zero threshold, firing trivially for any non-degenerate DS library.
            float signal    = profile.s_gt;
            float threshold = is_ss_lib ? 0.004f : 0.006f;

            // Require the terminal component to also be non-trivial (rules out flat noise)
            bool has_data   = total_w >= 500.0f;
            bool elevated   = signal > threshold;
            bool not_inverted = best_A_raw >= 0.0f;

            // C5: the GT model-fit verdict overwrites ox_damage_detected (the legacy
            // field surfaced as channel_c_detected) but does NOT touch ox_d_max, which
            // stays coupled to the codon-block verdict (ox_damage_detected_codon /
            // ox_is_artifact). Record the model verdict explicitly so the two can be
            // reconciled by consumers.
            profile.ox_damage_detected_model = has_data && elevated && not_inverted;
            profile.ox_damage_detected = profile.ox_damage_detected_model;
        }
    }

    // Reference-free oxidation-like signal. Within each length x GC bin, terminal
    // deamination strata are calibrated against the same bin's interior base
    // composition. High-deamination strata are compared to low-deamination strata;
    // a bin contributes only when signal and matched negative-control coverage are
    // both available. The same interior bases never enter both arms.
    {
        struct Arm {
            double sig_t = 0.0, sig_tg = 0.0;
            double sig_a = 0.0, sig_ac = 0.0;
            double ctrl_a = 0.0, ctrl_at = 0.0;
            double ctrl_c = 0.0, ctrl_cg = 0.0;
            double deam_score_reads = 0.0;
            double reads = 0.0;
        };
        struct StratumScore {
            double score = 0.0;
            double z = 0.0;
            bool usable = false;
        };
        struct BinResult {
            double signal = 0.0, signal_var = 0.0;
            double control = 0.0, control_var = 0.0;
            double delta = 0.0, delta_var = 0.0;
        };

        const auto binom_var = [](double p, double n) {
            return (n > 0.0) ? p * (1.0 - p) / n : std::numeric_limits<double>::infinity();
        };
        const auto add_stratum_to_arm = [](Arm& a,
                                           const SampleDamageProfile::OxidationLikeBin::Stratum& s,
                                           double deam_score) {
            a.sig_t += s.sig_t; a.sig_tg += s.sig_tg;
            a.sig_a += s.sig_a; a.sig_ac += s.sig_ac;
            a.ctrl_a += s.ctrl_a; a.ctrl_at += s.ctrl_at;
            a.ctrl_c += s.ctrl_c; a.ctrl_cg += s.ctrl_cg;
            a.deam_score_reads += deam_score * static_cast<double>(s.reads);
            a.reads += static_cast<double>(s.reads);
        };

        std::vector<BinResult> bin_results;
        bin_results.reserve(SampleDamageProfile::N_OX_BINS);
        const bool use_3prime_deam =
            profile.library_type != SampleDamageProfile::LibraryType::SINGLE_STRANDED;

        for (const auto& b : profile.oxidation_like_bins) {
            double bin_int_t = 0.0, bin_int_tc = 0.0;
            double bin_int_a = 0.0, bin_int_ag = 0.0;
            for (const auto& s : b.strata) {
                bin_int_t += s.int_t; bin_int_tc += s.int_tc;
                bin_int_a += s.int_a; bin_int_ag += s.int_ag;
            }
            if (bin_int_tc < 500.0 || bin_int_ag < 500.0) continue;

            std::array<StratumScore, SampleDamageProfile::N_OX_DEAM_STRATA> scores = {};
            double max_score = -std::numeric_limits<double>::infinity();
            for (int i = 0; i < SampleDamageProfile::N_OX_DEAM_STRATA; ++i) {
                const auto& s = b.strata[i];
                double sw = 0.0, sx = 0.0;

                if (s.term_tc5 >= 20.0) {
                    const double p = s.term_t5 / s.term_tc5;
                    const double q = bin_int_t / bin_int_tc;
                    const double var = binom_var(p, s.term_tc5) + binom_var(q, bin_int_tc);
                    if (var > 0.0 && std::isfinite(var)) {
                        const double w = 1.0 / var;
                        sw += w;
                        sx += w * (p - q);
                    }
                }
                if (use_3prime_deam && s.term_ag3 >= 20.0) {
                    const double p = s.term_a3 / s.term_ag3;
                    const double q = bin_int_a / bin_int_ag;
                    const double var = binom_var(p, s.term_ag3) + binom_var(q, bin_int_ag);
                    if (var > 0.0 && std::isfinite(var)) {
                        const double w = 1.0 / var;
                        sw += w;
                        sx += w * (p - q);
                    }
                }

                if (sw > 0.0) {
                    scores[i].score = sx / sw;
                    scores[i].z = scores[i].score / std::sqrt(1.0 / sw);
                    scores[i].usable = true;
                    max_score = std::max(max_score, scores[i].score);
                }
            }
            if (!(max_score > 0.015)) continue;

            Arm anc, bg;
            const double ancient_cut = std::max(0.02, 0.5 * max_score);
            const double bg_cut = std::max(0.005, 0.15 * max_score);
            for (int i = 0; i < SampleDamageProfile::N_OX_DEAM_STRATA; ++i) {
                if (!scores[i].usable) continue;
                const auto& s = b.strata[i];
                if (scores[i].score >= ancient_cut && scores[i].z > 0.5) {
                    add_stratum_to_arm(anc, s, scores[i].score);
                } else if (scores[i].score <= bg_cut || scores[i].z <= 0.0) {
                    add_stratum_to_arm(bg, s, scores[i].score);
                }
            }

            if (anc.reads < 10.0 || bg.reads < 10.0) continue;
            const double anc_deam = anc.deam_score_reads / anc.reads;
            const double bg_deam = bg.deam_score_reads / bg.reads;
            if ((anc_deam - bg_deam) < 0.01) continue;

            if (anc.sig_tg < 100.0 || bg.sig_tg < 100.0 ||
                anc.sig_ac < 100.0 || bg.sig_ac < 100.0 ||
                anc.ctrl_at < 100.0 || bg.ctrl_at < 100.0 ||
                anc.ctrl_cg < 100.0 || bg.ctrl_cg < 100.0) {
                continue;
            }

            const double p_at = anc.sig_t / anc.sig_tg;
            const double p_bt = bg.sig_t / bg.sig_tg;
            const double p_aa = anc.sig_a / anc.sig_ac;
            const double p_ba = bg.sig_a / bg.sig_ac;
            const double gt_delta = p_at - p_bt;
            const double ca_delta = p_aa - p_ba;
            const double signal = 0.5 * (gt_delta + ca_delta);
            const double signal_var = 0.25 * (
                binom_var(p_at, anc.sig_tg) + binom_var(p_bt, bg.sig_tg) +
                binom_var(p_aa, anc.sig_ac) + binom_var(p_ba, bg.sig_ac));

            const double p_ac1 = anc.ctrl_a / anc.ctrl_at;
            const double p_bc1 = bg.ctrl_a / bg.ctrl_at;
            const double p_ac2 = anc.ctrl_c / anc.ctrl_cg;
            const double p_bc2 = bg.ctrl_c / bg.ctrl_cg;
            const double ctrl_1 = p_ac1 - p_bc1;
            const double ctrl_2 = p_ac2 - p_bc2;
            const double control = 0.5 * (std::abs(ctrl_1) + std::abs(ctrl_2));
            const double control_var = 0.25 * (
                binom_var(p_ac1, anc.ctrl_at) + binom_var(p_bc1, bg.ctrl_at) +
                binom_var(p_ac2, anc.ctrl_cg) + binom_var(p_bc2, bg.ctrl_cg));
            const double delta_var = signal_var + control_var;
            if (!(delta_var > 0.0) || !std::isfinite(delta_var)) continue;
            bin_results.push_back({signal, signal_var, control, control_var,
                                   signal - control, delta_var});
        }

        if (!bin_results.empty()) {
            double total_w0 = 0.0;
            for (const auto& r : bin_results) total_w0 += 1.0 / r.delta_var;
            const double cap = total_w0 * 0.25;  // no bin can dominate the meta-estimate

            double sw = 0.0, sw2 = 0.0;
            double signal = 0.0, control = 0.0, delta = 0.0;
            std::vector<double> weights;
            weights.reserve(bin_results.size());
            for (const auto& r : bin_results) {
                const double w = std::min(1.0 / r.delta_var, cap);
                weights.push_back(w);
                sw += w;
                sw2 += w * w;
                signal += w * r.signal;
                control += w * r.control;
                delta += w * r.delta;
            }
            if (sw > 0.0) {
            signal /= sw;
            control /= sw;
            delta /= sw;

            double signal_var = 0.0, control_var = 0.0, delta_var = 0.0;
            double q = 0.0;
            for (size_t i = 0; i < bin_results.size(); ++i) {
                const double a = weights[i] / sw;
                signal_var += a * a * bin_results[i].signal_var;
                control_var += a * a * bin_results[i].control_var;
                delta_var += a * a * bin_results[i].delta_var;
                q += weights[i] * (bin_results[i].delta - delta) *
                     (bin_results[i].delta - delta);
            }
            const double signal_se = std::sqrt(std::max(0.0, signal_var));
            const double control_se = std::sqrt(std::max(0.0, control_var));
            // C2 (behavioral): the fixed-effect binomial delta_var ignores between-bin
            // heterogeneity (I²≈0.99 across the panel), so z scaled ~sqrt(n_reads)
            // (|z|≈300–550). Apply the standard DerSimonian-Laird random-effects floor
            // adjusted_var = max(delta_var, q/sw): when q=0 / single bin it degrades to
            // the binomial variance (no change). This bounds z by the effective number
            // of independent bins, not by read count.
            const double adjusted_var = std::max(delta_var, sw > 0.0 ? q / sw : delta_var);
            const double delta_se = std::sqrt(std::max(0.0, adjusted_var));
            const double excess = std::max(0.0, delta);
            const double z = delta_se > 0.0 ? delta / delta_se : 0.0;
            const double effective_bins = sw2 > 0.0 ? (sw * sw) / sw2 : 0.0;
            const double df = static_cast<double>(bin_results.size() > 0 ? bin_results.size() - 1 : 0);
            const double heterogeneity = (df > 0.0 && q > 0.0)
                ? std::clamp((q - df) / q, 0.0, 1.0)
                : 0.0;

            profile.oxidation_like_signal = static_cast<float>(signal);
            profile.oxidation_like_signal_se = static_cast<float>(signal_se);
            profile.oxidation_like_control = static_cast<float>(control);
            profile.oxidation_like_control_se = static_cast<float>(control_se);
            profile.oxidation_like_adjusted = static_cast<float>(delta);
            profile.oxidation_like_excess = static_cast<float>(excess);
            profile.oxidation_like_se = static_cast<float>(delta_se);
            profile.oxidation_like_z = static_cast<float>(z);
            profile.oxidation_like_bins_used = static_cast<int>(bin_results.size());
            profile.oxidation_like_effective_bins = static_cast<float>(effective_bins);
            profile.oxidation_like_heterogeneity = static_cast<float>(heterogeneity);
            profile.oxidation_like_artifact_suspect =
                control > std::max(0.002, std::abs(signal) * 0.75);

            // C5 (behavioral): reliability_score must measure ESTIMATION QUALITY
            // (data sufficiency), not signal direction. The old formula multiplied by
            // (1-exp(-max(0,z)/5)), which is 0 whenever z<=0 — and z is negative for
            // every non-oxidised library by construction, so reliability_score was
            // structurally always 0. Drop the z_score factor; the signed verdict lives
            // in oxidation_like_z. bin_score / het_penalty / artifact_penalty already
            // capture data quality and are each in [0,1].
            const double bin_score = std::min(1.0, effective_bins / 6.0);
            const double het_penalty = 1.0 - 0.5 * heterogeneity;
            const double artifact_penalty = profile.oxidation_like_artifact_suspect ? 0.5 : 1.0;
            profile.oxidation_like_reliability =
                static_cast<float>(std::clamp(bin_score * het_penalty *
                                              artifact_penalty, 0.0, 1.0));
            }
        }
    }

    // === Correction 4: Chargaff gate for the s_gt strand-asymmetric oxidation contrast ===
    // s_gt = best_B - ox_ca_baseline is a STRAND-ASYMMETRIC oxidation contrast, not total oxidation.
    // It is interpretable only when interior base composition is Chargaff-balanced, so the contrast
    // is not driven by a |G|≠|C| composition skew. Gate (not null): add a validity flag and the
    // measured balance; do NOT alter s_gt itself. baseline_{g,c}_freq are interior (mid-read) fractions
    // by this point; |G-C|/(G+C) is scale-invariant so it reads the same on counts or fractions.
    {
        const double gi = profile.baseline_g_freq;
        const double ci = profile.baseline_c_freq;
        const double gc = gi + ci;
        profile.chargaff_gc_balance = (gc > 0.0) ? static_cast<float>(std::fabs(gi - ci) / gc) : 0.0f;
        profile.s_gt_valid = (profile.chargaff_gc_balance <= 0.02f) && profile.d_computed;
    }
}

} // namespace taph
