#include "damage_estimation_finalize_helpers.hpp"
#include "damage_estimation_finalize_ctx.hpp"
namespace taph {

void finalize_context(SampleDamageProfile& profile, FinalCtx& ctx) {
    // Channel E is computed in finalize_init before raw terminal counts are
    // normalized. Do not recompute it here from mutated rate arrays.

    // Update lambda from fit if reasonable
    // D22: inclusive bounds match the fit's internal clamp range [0.05,0.50]; a steep
    // decay clamped to exactly 0.5 is a real fit result, not the 0.3 init sentinel.
    if (ctx.fit_lambda_5p >= 0.05f && ctx.fit_lambda_5p <= 0.5f) {
        profile.lambda_5prime = ctx.fit_lambda_5p;
        profile.lambda_5prime_fitted = true;
    }
    if (ctx.fit_lambda_3p >= 0.05f && ctx.fit_lambda_3p <= 0.5f) {
        profile.lambda_3prime = ctx.fit_lambda_3p;
        profile.lambda_3prime_fitted = true;
    }

    // Context-split amplitude fits (CpG-like vs non-CpG-like)
    {
        auto fit_cpg = fit_ct5_ctx_amplitude(
            profile.ct_ctx_t_5prime[SampleDamageProfile::CPG_LIKE],
            profile.ct_ctx_total_5prime[SampleDamageProfile::CPG_LIKE],
            profile.ct_ctx_t_interior[SampleDamageProfile::CPG_LIKE],
            profile.ct_ctx_total_interior[SampleDamageProfile::CPG_LIKE],
            profile.lambda_5prime);

        auto fit_ncpg = fit_ct5_ctx_amplitude(
            profile.ct_ctx_t_5prime[SampleDamageProfile::NONCPG_LIKE],
            profile.ct_ctx_total_5prime[SampleDamageProfile::NONCPG_LIKE],
            profile.ct_ctx_t_interior[SampleDamageProfile::NONCPG_LIKE],
            profile.ct_ctx_total_interior[SampleDamageProfile::NONCPG_LIKE],
            profile.lambda_5prime);

        profile.fit_baseline_ct5_cpg_like    = fit_cpg.baseline;
        profile.fit_baseline_ct5_noncpg_like = fit_ncpg.baseline;
        profile.dmax_ct5_cpg_like    = fit_cpg.valid  ? fit_cpg.dmax  : std::numeric_limits<float>::quiet_NaN();
        profile.dmax_ct5_noncpg_like = fit_ncpg.valid ? fit_ncpg.dmax : std::numeric_limits<float>::quiet_NaN();
        profile.cov_ct5_cpg_like_terminal       = fit_cpg.cov_terminal;
        profile.cov_ct5_noncpg_like_terminal    = fit_ncpg.cov_terminal;
        profile.cov_ct5_cpg_like_interior       = fit_cpg.cov_interior;
        profile.cov_ct5_noncpg_like_interior    = fit_ncpg.cov_interior;
        profile.effcov_ct5_cpg_like_terminal    = fit_cpg.effcov_terminal;
        profile.effcov_ct5_noncpg_like_terminal = fit_ncpg.effcov_terminal;
        profile.effcov_ct5_cpg_like_interior    = fit_cpg.effcov_interior;
        profile.effcov_ct5_noncpg_like_interior = fit_ncpg.effcov_interior;
        profile.fit_positions_ct5_cpg_like    = fit_cpg.fit_positions;
        profile.fit_positions_ct5_noncpg_like = fit_ncpg.fit_positions;

        if (fit_cpg.valid && fit_ncpg.valid && fit_ncpg.dmax >= 0.015f) {
            profile.cpg_ratio    = fit_cpg.dmax / fit_ncpg.dmax;
            profile.log2_cpg_ratio = std::log2((fit_cpg.dmax + 1e-6f) / (fit_ncpg.dmax + 1e-6f));
            profile.cpg_ratio_backwards = (profile.log2_cpg_ratio < 0.0f);  // P5 QC: methylated aDNA expects CpG > non-CpG
        }

        // oxoG 16-context finalization
        for (int i = 0; i < SampleDamageProfile::N_OXOG16; ++i) {
            const float t = profile.oxog16_t[i], a = profile.oxog16_a_rc[i];
            const float cov = t + a;
            profile.cov_oxog_16ctx[i] = cov;
            profile.s_oxog_16ctx[i] = (cov >= 500.0f)
                ? (t - a) / cov
                : std::numeric_limits<float>::quiet_NaN();
        }

        // Upstream-context-aware C→T fitting (experimental: AC, CC, GC, TC)
        // Fit each context independently, then compute contrasts
        float dmax_sum = 0.0f;
        int valid_ctx_count = 0;
        for (int uctx = 0; uctx < SampleDamageProfile::N_UPSTREAM_CTX; ++uctx) {
            auto fit = fit_ct5_ctx_amplitude(
                profile.ct5_t_by_upstream[uctx],
                profile.ct5_total_by_upstream[uctx],
                profile.ct5_t_interior_by_upstream[uctx],
                profile.ct5_total_interior_by_upstream[uctx],
                profile.lambda_5prime);

            profile.baseline_ct5_by_upstream[uctx] = fit.baseline;
            profile.dmax_ct5_by_upstream[uctx] = fit.valid ? fit.dmax : std::numeric_limits<float>::quiet_NaN();
            profile.cov_ct5_terminal_by_upstream[uctx] = fit.cov_terminal;
            profile.cov_ct5_interior_by_upstream[uctx] = fit.cov_interior;

            // C4: lower-boundary fit contributes 0 to the contrast (the correct signal-absent value)
            // but still counts toward valid_ctx_count because coverage was sufficient to run the fit.
            const float dmax_for_contrast = fit.valid ? fit.dmax : (fit.lower_boundary ? 0.0f : std::numeric_limits<float>::quiet_NaN());
            if (std::isfinite(dmax_for_contrast)) {
                dmax_sum += dmax_for_contrast;
                ++valid_ctx_count;
            }
        }

        // Compute contrasts if all 4 contexts are valid
        if (valid_ctx_count == 4) {
            const float d_ac = profile.dmax_ct5_by_upstream[SampleDamageProfile::CTX_AC];
            const float d_cc = profile.dmax_ct5_by_upstream[SampleDamageProfile::CTX_CC];
            const float d_gc = profile.dmax_ct5_by_upstream[SampleDamageProfile::CTX_GC];
            const float d_tc = profile.dmax_ct5_by_upstream[SampleDamageProfile::CTX_TC];

            // dipyr_contrast: mean(CC, TC) - mean(AC, GC)
            // Positive = dipyrimidine excess (UV-like)
            profile.dipyr_contrast = 0.5f * (d_cc + d_tc) - 0.5f * (d_ac + d_gc);

            // The upstream-GC arm (d_gc) has no deamination-enhancing chemistry in
            // eukaryotes (upstream-G = GpC, not CpG). The correct CpG/5mC axis is
            // the downstream cpg_like block → log2_cpg_ratio / cpg_z.

            // Chi-squared heterogeneity test: are contexts different from uniform?
            const float mean_d = dmax_sum / 4.0f;
            if (mean_d > 0.001f) {
                float chi2 = 0.0f;
                for (int uctx = 0; uctx < 4; ++uctx) {
                    const float d = profile.dmax_ct5_by_upstream[uctx];
                    const float cov = profile.cov_ct5_terminal_by_upstream[uctx];
                    if (cov > 100.0f) {
                        // Weight by coverage (pseudo chi-squared)
                        chi2 += cov * (d - mean_d) * (d - mean_d) / (mean_d * (1.0f - mean_d) + 1e-6f);
                    }
                }
                // C2: this is a coverage-weighted sum (each term ∝ cov_i), so the
                // statistic grows ~O(total reads) and is NOT a distribution-free
                // chi²(3). The p-value lookup below is therefore only nominal — at
                // high coverage any real between-context difference trivially exceeds
                // 7.81. Treat context_heterogeneity_chi2 as a coverage-scaled index,
                // not a calibrated test statistic (consumers should prefer the
                // effect-size contrast dipyr_contrast). Detection logic unchanged.
                // C2: keep the unclamped coverage-scaled value for audit, but emit the
                // primary statistic clamped to kZCap so consumers do not read an
                // explosive magnitude as a calibrated chi². Detection still gates on the
                // raw value (7.81 < kZCap, so the call is unchanged).
                profile.context_heterogeneity_chi2_raw = chi2;
                profile.context_heterogeneity_chi2 = std::min(chi2, SampleDamageProfile::kZCap);
                // Nominal p-value (chi² with 3 df); not calibrated — see note above.
                profile.context_heterogeneity_detected = (chi2 > 7.81f);
                // Store approximate p-value (very rough)
                if (chi2 < 0.58f) profile.context_heterogeneity_p = 0.9f;
                else if (chi2 < 2.37f) profile.context_heterogeneity_p = 0.5f;
                else if (chi2 < 6.25f) profile.context_heterogeneity_p = 0.1f;
                else if (chi2 < 7.81f) profile.context_heterogeneity_p = 0.05f;
                else if (chi2 < 11.34f) profile.context_heterogeneity_p = 0.01f;
                else profile.context_heterogeneity_p = 0.001f;
                profile.context_heterogeneity_computed = true;
            }
        }
    }

    // Interior clustered C→T finalization
    {
        const auto& acc = profile.interior_ct_cluster;
        auto safe_log2_oe = [](uint64_t obs, double exp) -> float {
            return static_cast<float>(
                std::log((static_cast<double>(obs) + 0.5) / (exp + 0.5)) / std::log(2.0));
        };

        uint64_t obs_ct_s = 0, pairs_ct_s = 0, obs_ag_s = 0;
        double   exp_ct_s = 0.0, var_ct_s = 0.0;
        double   exp_ag_s = 0.0, var_ag_s = 0.0;

        for (int d = 1; d <= 10; ++d) {
            profile.interior_ct_cluster_sep_log2oe[d - 1] =
                safe_log2_oe(acc.obs_ct[d], acc.exp_ct[d]);
            if (d <= 5) {
                obs_ct_s  += acc.obs_ct[d];   pairs_ct_s += acc.pairs_ct[d];
                exp_ct_s  += acc.exp_ct[d];   var_ct_s   += acc.var_ct[d];
                obs_ag_s  += acc.obs_ag[d];
                exp_ag_s  += acc.exp_ag[d];   var_ag_s   += acc.var_ag[d];
            }
        }

        profile.interior_ct_cluster_short_obs          = obs_ct_s;
        profile.interior_ct_cluster_short_pairs        = pairs_ct_s;
        profile.interior_ct_cluster_short_exp          = exp_ct_s;
        profile.interior_ct_cluster_reads_used         = acc.reads_used_ct;
        profile.interior_ct_cluster_reads_used_control = acc.reads_used_ag;
        profile.interior_ct_cluster_reads_skipped      = acc.short_reads_skipped;

        profile.interior_ct_cluster_short_log2oe      = safe_log2_oe(obs_ct_s, exp_ct_s);
        const float ag_log2oe                         = safe_log2_oe(obs_ag_s, exp_ag_s);
        profile.interior_ct_cluster_short_asym_log2oe =
            profile.interior_ct_cluster_short_log2oe - ag_log2oe;

        // C2: short_z is a CT-vs-AG asymmetry contrast whose denominator is a
        // within-read Bernoulli pair variance that underestimates the true variance
        // (T-sites within a read share one deamination rate -> positive correlation),
        // so the magnitude is not a calibrated p-value and is not comparable across
        // ds/ss library types. It measures CT-SPECIFIC asymmetry, not total CT
        // clustering intensity — short_log2oe / short_asym_log2oe are the primary
        // effect-size metrics. Clamp the emitted magnitude to ±kZCap.
        const double denom = std::sqrt(var_ct_s + var_ag_s + 1e-12);
        const double num   = (static_cast<double>(obs_ct_s) - exp_ct_s)
                           - (static_cast<double>(obs_ag_s) - exp_ag_s);
        const double scz_cap = static_cast<double>(SampleDamageProfile::kZCap);
        profile.interior_ct_cluster_short_z =
            static_cast<float>(std::clamp(num / denom, -scz_cap, scz_cap));
    }

    // Step 3: Per-position damage rates using fit baseline
    float fit_baseline_c_frac_5p = 1.0f - profile.fit_baseline_5prime;
    float fit_baseline_g_frac_3p = 1.0f - profile.fit_baseline_3prime;

    for (int i = 0; i < 15; ++i) {
        if (fit_baseline_c_frac_5p > 0.1f) {
            float raw_signal = static_cast<float>(profile.t_freq_5prime[i]) - profile.fit_baseline_5prime;
            profile.damage_rate_5prime[i] = std::max(0.0f, raw_signal / fit_baseline_c_frac_5p);
        } else {
            profile.damage_rate_5prime[i] = 0.0f;
        }

        if (fit_baseline_g_frac_3p > 0.1f) {
            float raw_signal = static_cast<float>(profile.a_freq_3prime[i]) - profile.fit_baseline_3prime;
            profile.damage_rate_3prime[i] = std::max(0.0f, raw_signal / fit_baseline_g_frac_3p);
        } else {
            profile.damage_rate_3prime[i] = 0.0f;
        }
    }

    // Terminal gradients for reporting; inverted_pattern flags set later by hexamer detection
    {
        double terminal_tc_5 = profile.t_freq_5prime[0];
        profile.terminal_gradient_5prime = static_cast<float>(terminal_tc_5 - ctx.baseline_tc);

        double terminal_ag_3 = profile.a_freq_3prime[0];
        profile.terminal_gradient_3prime = static_cast<float>(terminal_ag_3 - ctx.baseline_ag);
    }

    for (int p = 0; p < 3; p++) {
        size_t tc_total = profile.codon_pos_t_count_5prime[p] + profile.codon_pos_c_count_5prime[p];
        if (tc_total > 0) {
            profile.codon_pos_t_fraction_5prime[p] = static_cast<float>(profile.codon_pos_t_count_5prime[p]) / tc_total;
        }

        size_t ag_total = profile.codon_pos_a_count_3prime[p] + profile.codon_pos_g_count_3prime[p];
        if (ag_total > 0) {
            profile.codon_pos_a_fraction_3prime[p] = static_cast<float>(profile.codon_pos_a_count_3prime[p]) / ag_total;
        }
    }

    {
        float pos1_t_rate = profile.codon_pos_t_fraction_5prime[0];
        float pos2_t_rate = profile.codon_pos_t_fraction_5prime[1];
        float pos3_t_rate = profile.codon_pos_t_fraction_5prime[2];  // Wobble position

        float baseline_tc_f = profile.fit_baseline_5prime;
        float baseline_c_f = 1.0f - profile.fit_baseline_5prime;

        if (baseline_c_f > 0.1f) {
            profile.codon_pos1_damage = std::max(0.0f, (pos1_t_rate - baseline_tc_f) / baseline_c_f);
            profile.codon_pos2_damage = std::max(0.0f, (pos2_t_rate - baseline_tc_f) / baseline_c_f);
            profile.codon_pos3_damage = std::max(0.0f, (pos3_t_rate - baseline_tc_f) / baseline_c_f);
        } else {
            profile.codon_pos1_damage = 0.0f;
            profile.codon_pos2_damage = 0.0f;
            profile.codon_pos3_damage = 0.0f;
        }

        // Calculate wobble ratio: pos3 / ((pos1 + pos2) / 2)
        // F4: gate at 0.015 (was 0.005); report on log2 scale; remove the [0.5,3] clamp.
        float avg_pos12_damage = (profile.codon_pos1_damage + profile.codon_pos2_damage) / 2.0f;
        if (avg_pos12_damage > 0.015f) {
            profile.wobble_ratio = profile.codon_pos3_damage / avg_pos12_damage;
            profile.log2_wobble_ratio =
                std::log2((profile.codon_pos3_damage + 1e-6f) / (avg_pos12_damage + 1e-6f));
        } else {
            profile.wobble_ratio = 1.0f;
            profile.log2_wobble_ratio = 0.0f;
        }
    }

    size_t cpg_total = profile.cpg_c_count + profile.cpg_t_count;
    if (cpg_total > 10) {
        profile.cpg_ct_fraction = static_cast<float>(profile.cpg_t_count) / cpg_total;
    }

    size_t non_cpg_total = profile.non_cpg_c_count + profile.non_cpg_t_count;
    if (non_cpg_total > 10) {
        profile.non_cpg_ct_fraction = static_cast<float>(profile.non_cpg_t_count) / non_cpg_total;
    }

    profile.max_damage_5prime = profile.damage_rate_5prime[0];
    profile.max_damage_3prime = profile.damage_rate_3prime[0];

    // Briggs-like damage model parameter estimation (closed-form).
    // start_pos: first position to use (1 when chemistry tag suppresses pos 0).
    // fixed_bg: when >= 0, use this as the background instead of the tail mean
    // computed from rates[10..14]. Tail-anchored bg gives a chemistry-robust
    // baseline that does not trade off with d_max in the fit.
    auto estimate_briggs_params = [](const std::array<float, 15>& rates, float max_rate,
                                     int start_pos, float fixed_bg,
                                     float& delta_s, float& delta_d, float& lambda, float& r_squared,
                                     bool& lambda_fitted) {
        delta_s = 0.0f;
        delta_d = 0.0f;
        lambda = 0.3f;
        r_squared = 0.0f;
        lambda_fitted = false;  // P5: only the converged decay regression below sets this true. This
                                // estimator produces the EMITTED lambda, so the emitter gate
                                // (lambda_*_fitted ? value : null) now nulls the 0.3 bail-out sentinel
                                // instead of leaking it as a real fit.

        if (max_rate < 0.02f) return;
        if (start_pos < 0 || start_pos > 5) start_pos = 0;

        delta_s = rates[start_pos];

        if (fixed_bg >= 0.0f) {
            delta_d = fixed_bg;
        } else {
            float sum_background = 0.0f;
            for (int i = 10; i < 15; ++i) sum_background += rates[i];
            delta_d = sum_background / 5.0f;
        }

        if (delta_s <= delta_d + 0.01f) {
            lambda = 0.3f;
            return;
        }

        // Weighted log-linear regression for lambda estimation.
        // Fit log((rate-bg)/(d_s-bg)) = -slope * (pos - start_pos) over
        // start_pos..14. x is the offset from start_pos so the intercept
        // corresponds to delta_s.
        float sum_xy = 0.0f, sum_x2 = 0.0f, sum_y = 0.0f, sum_x = 0.0f;
        float sum_w  = 0.0f, sum_y2 = 0.0f;

        for (int pos = start_pos; pos < 15; ++pos) {
            float normalized = (rates[pos] - delta_d) / (delta_s - delta_d + 0.001f);
            normalized = std::clamp(normalized, 0.01f, 1.0f);
            float y = std::log(normalized);
            float x = static_cast<float>(pos - start_pos);
            float w = normalized * normalized;
            sum_xy += w * x * y;
            sum_x2 += w * x * x;
            sum_y  += w * y;
            sum_x  += w * x;
            sum_w  += w;
            sum_y2 += w * y * y;
        }

        float denom = sum_w * sum_x2 - sum_x * sum_x;
        if (std::abs(denom) < 0.001f) { lambda = 0.3f; return; }

        float slope = (sum_w * sum_xy - sum_x * sum_y) / denom;
        if (slope >= 0.0f) { lambda = 0.3f; return; }

        lambda = 1.0f - std::exp(slope);
        lambda = std::clamp(lambda, 0.05f, 0.95f);
        lambda_fitted = true;  // P5: real decay fit converged (cleared every 0.3 sentinel guard above)

        float ss_tot = sum_y2 - (sum_y * sum_y) / sum_w;
        float intercept = (sum_y - slope * sum_x) / sum_w;
        float ss_res = 0.0f;
        for (int pos = start_pos; pos < 15; ++pos) {
            float normalized = (rates[pos] - delta_d) / (delta_s - delta_d + 0.001f);
            normalized = std::clamp(normalized, 0.01f, 1.0f);
            float y = std::log(normalized);
            float y_pred = slope * static_cast<float>(pos - start_pos) + intercept;
            float w = normalized * normalized;
            ss_res += w * (y - y_pred) * (y - y_pred);
        }
        r_squared = (ss_tot > 0.001f) ? std::max(0.0f, 1.0f - ss_res / ss_tot) : 0.0f;
    };

    // Tail-anchored background: trimmed mean of per-position C->T (G->A)
    // rates over BG_TAIL_LO..BG_TAIL_HI, requiring denom >= 100. If too few
    // eligible positions, falls through to the legacy joint-fit bg.
    auto compute_anchored_bg = [](const std::array<double, SampleDamageProfile::BG_TAIL_N>& num,
                                  const std::array<double, SampleDamageProfile::BG_TAIL_N>& den,
                                  int min_n_positions, double trim_frac,
                                  float& bg_out, int& n_pos_out, double& denom_out) {
        bg_out = -1.0f;
        n_pos_out = 0;
        denom_out = 0.0;
        std::vector<float> rates;
        rates.reserve(SampleDamageProfile::BG_TAIL_N);
        for (int i = 0; i < SampleDamageProfile::BG_TAIL_N; ++i) {
            if (den[i] >= 100.0) {
                rates.push_back(static_cast<float>(num[i] / den[i]));
                denom_out += den[i];
            }
        }
        n_pos_out = static_cast<int>(rates.size());
        if (n_pos_out < min_n_positions) return;
        std::sort(rates.begin(), rates.end());
        int n_trim = static_cast<int>(std::floor(trim_frac * n_pos_out));
        int lo = n_trim, hi = n_pos_out - n_trim;
        if (hi <= lo) { lo = 0; hi = n_pos_out; }
        double sum = 0.0;
        for (int i = lo; i < hi; ++i) sum += rates[i];
        bg_out = static_cast<float>(sum / static_cast<double>(hi - lo));
    };

    {
        float  bg5 = -1.0f, bg3 = -1.0f;
        int    n5  = 0,     n3  = 0;
        double d5  = 0.0,   d3  = 0.0;
        compute_anchored_bg(profile.tail_t_5prime,  profile.tail_tc_5prime,  10, 0.20, bg5, n5, d5);
        compute_anchored_bg(profile.tail_a_3prime,  profile.tail_ag_3prime,  10, 0.20, bg3, n3, d3);
        profile.bg_5prime_anchored    = (bg5 >= 0.0f) ? bg5 : 0.0f;
        profile.bg_3prime_anchored    = (bg3 >= 0.0f) ? bg3 : 0.0f;
        profile.bg_n_positions_5prime = n5;
        profile.bg_n_positions_3prime = n3;
        profile.bg_denominator_5prime = d5;
        profile.bg_denominator_3prime = d3;

        // Inline chemistry-tag detection: top 5' hexamer matches a curated
        // protocol-tag table at log2fc >= 3.0. Populates profile.protocol_tag_*
        // here (early) so the Briggs fit can mask pos 0; the post-bias-gate
        // rescue further down only mutates library_type if needed.
        if (profile.protocol_tag_5prime.empty()) {
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
                }
            }
        }
        const bool tag5 = !profile.protocol_tag_5prime.empty()
                          && profile.protocol_tag_log2fc >= 3.0f;
        profile.briggs_pos0_masked_5prime = tag5;
        // No 3' chemistry-tag table currently; mirror flag stays false.
        profile.briggs_pos0_masked_3prime = false;

        const int  start5 = profile.briggs_pos0_masked_5prime ? 1 : 0;
        const int  start3 = profile.briggs_pos0_masked_3prime ? 1 : 0;
        const float fb5   = (bg5 >= 0.0f) ? bg5 : -1.0f;
        const float fb3   = (bg3 >= 0.0f) ? bg3 : -1.0f;

        // DELETED-REF-FREE-PATH (kept): estimate_briggs_params fits the per-POSITION λ (wrong axis for the
        // bulk/length law — finalize_tau now owns the correct READ-LENGTH decay τ on the ref-free path).
        // It is NOT deleted because profile.lambda_5prime/3prime it writes here are SHARED with legacy
        // consumers (damage_c_api, library_interpretation_summary, llr_table, finalize_dmax) AND with the
        // ref-free per-read LLR decay (read_ancient_llr.cpp). The ref-free GATE no longer depends on this
        // (finalize_pi gates on tau.state), but the per-read decay constant still does; nulling λ here
        // would silently change every legacy path. Removing it requires migrating read_ancient_llr to a
        // τ-derived per-read decay first — out of scope for the additive τ step.
        estimate_briggs_params(profile.damage_rate_5prime, profile.max_damage_5prime,
                               start5, fb5,
                               profile.delta_s_5prime, profile.delta_d_5prime,
                               profile.lambda_5prime, profile.r_squared_5prime,
                               profile.lambda_5prime_fitted);

        estimate_briggs_params(profile.damage_rate_3prime, profile.max_damage_3prime,
                               start3, fb3,
                               profile.delta_s_3prime, profile.delta_d_3prime,
                               profile.lambda_3prime, profile.r_squared_3prime,
                               profile.lambda_3prime_fitted);

        // Headline area-excess + LR companion vs bg-only null over k..14.
        // Uses RAW per-position T/(T+C) rates (t_freq_5prime / a_freq_3prime,
        // rate-valued post-finalize) and the tail-anchored bg in the same
        // P-space, so area_excess is comparable across libraries without
        // depending on the joint Briggs fit's interior baseline.
        auto compute_headline = [](const std::array<double, 15>& rates,
                                   const std::array<double, 15>& tc_total,
                                   int k, float bg,
                                   float& area_excess, float& lr) {
            area_excess = 0.0f;
            lr = 0.0f;
            if (bg < 0.0f) bg = 0.0f;
            const double bg_safe = std::clamp(static_cast<double>(bg), 1e-6, 1.0 - 1e-6);
            for (int p = k; p < 15; ++p) {
                const double r = rates[p];
                if (r > bg) area_excess += static_cast<float>(r - bg);
                const double n = tc_total[p];
                if (n < 5.0) continue;
                const double k_obs = r * n;
                const double p_alt = std::clamp(r, 1e-6, 1.0 - 1e-6);
                const double ll_alt  = k_obs * std::log(p_alt) + (n - k_obs) * std::log(1.0 - p_alt);
                const double ll_null = k_obs * std::log(bg_safe) + (n - k_obs) * std::log(1.0 - bg_safe);
                if (ll_alt > ll_null) lr += static_cast<float>(ll_alt - ll_null);
            }
        };
        compute_headline(profile.t_freq_5prime, profile.tc_total_5prime,
                         start5, profile.bg_5prime_anchored,
                         profile.damage_5prime_area_excess, profile.damage_5prime_lr);
        compute_headline(profile.a_freq_3prime, profile.ag_total_3prime,
                         start3, profile.bg_3prime_anchored,
                         profile.damage_3prime_area_excess, profile.damage_3prime_lr);
    }

    profile.lambda_5prime = std::clamp(profile.lambda_5prime, 0.1f, 1.0f);
    profile.lambda_3prime = std::clamp(profile.lambda_3prime, 0.1f, 1.0f);

    // Library type: handled below after composition-bias flags are set

    // Flag inversions: terminals depleted relative to interior with statistical support
    // Use both shift threshold AND z-score override for very strong statistical signals
    const bool inversion_5 = (profile.terminal_shift_5prime < -0.01f && profile.terminal_z_5prime < -2.0f)
                           || profile.terminal_z_5prime < -10.0f;  // Very strong depletion signal
    const bool inversion_3 = (profile.terminal_shift_3prime < -0.01f && profile.terminal_z_3prime < -2.0f)
                           || profile.terminal_z_3prime < -10.0f;  // Very strong depletion signal
    profile.terminal_inversion = inversion_5 || inversion_3;

    // Set individual inverted pattern flags (used for asymmetric pattern handling)
    // Note: inverted_pattern_5prime may also be set later by hexamer-based detection
    if (inversion_3) {
        profile.inverted_pattern_3prime = true;
    }
    if (inversion_5) {
        profile.inverted_pattern_5prime = true;
    }

    // Hexamer-based damage: compare C-initial vs T-initial hexamer frequencies at terminal vs interior
    if (profile.n_hexamers_5prime >= 1000 && profile.n_hexamers_interior >= 1000) {
        double llr_sum = 0.0;
        double weight_sum = 0.0;

        for (uint32_t base_code = 0; base_code < 1024; ++base_code) {
            uint32_t c_hex = 0x400 | base_code;
            uint32_t t_hex = 0xC00 | base_code;

            float expected_c = get_hexamer_freq(c_hex);
            float expected_t = get_hexamer_freq(t_hex);

            if (expected_c < 1e-6f || expected_t < 1e-6f) continue;

            double term_c = profile.hexamer_count_5prime[c_hex];
            double term_t = profile.hexamer_count_5prime[t_hex];
            double int_c = profile.hexamer_count_interior[c_hex];
            double int_t = profile.hexamer_count_interior[t_hex];

            double total_term = term_c + term_t;
            double total_int = int_c + int_t;
            if (total_term < 5 || total_int < 5) continue;

            double ratio_term = term_t / total_term;
            double ratio_int = int_t / total_int;
            double excess = ratio_term - ratio_int;

            double weight = (expected_c + expected_t) * std::sqrt(total_term * total_int);
            llr_sum += excess * weight;
            weight_sum += weight;
        }

        if (weight_sum > 0) {
            float raw_llr = static_cast<float>(llr_sum / weight_sum);

            // Store raw hexamer LLR without inversion correction
            // The raw value is the actual terminal-vs-interior hexamer composition difference
            // Positive = T-starting hexamers enriched at terminal (damage pattern)
            // Negative = C-starting hexamers enriched at terminal (AT-rich composition bias)
            profile.hexamer_damage_llr = raw_llr;

            (void)raw_llr;
        }

        // Compute hexamer-based T/(T+C) ratios (more reliable than position 0 or 1 alone)
        double total_term_t = 0, total_term_c = 0;
        double total_int_t = 0, total_int_c = 0;
        for (uint32_t base_code = 0; base_code < 1024; ++base_code) {
            uint32_t c_hex = 0x400 | base_code;
            uint32_t t_hex = 0xC00 | base_code;
            total_term_c += profile.hexamer_count_5prime[c_hex];
            total_term_t += profile.hexamer_count_5prime[t_hex];
            total_int_c += profile.hexamer_count_interior[c_hex];
            total_int_t += profile.hexamer_count_interior[t_hex];
        }
        double term_t_ratio = total_term_t / (total_term_t + total_term_c + 1e-10);
        double int_t_ratio = total_int_t / (total_int_t + total_int_c + 1e-10);

        profile.hexamer_terminal_tc = static_cast<float>(term_t_ratio);
        profile.hexamer_interior_tc = static_cast<float>(int_t_ratio);
        profile.hexamer_excess_tc = static_cast<float>(term_t_ratio - int_t_ratio);

        // Hexamers starting at pos 0 include the pos-0 artifact, so skip inversion
        // detection when a pos-0 artifact is present.
        if (profile.hexamer_damage_llr < -0.02f && !profile.position_0_artifact_5prime) {
            profile.inverted_pattern_5prime = true;
        }
    }

    // Composition bias: flag if |ctrl_shift| >= max(0.005, 0.5 * |damage_shift|)
    float damage_shift_5 = profile.terminal_shift_5prime;
    float ctrl_shift_5 = std::abs(profile.ctrl_shift_5prime);
    float threshold_5 = std::max(0.005f, 0.5f * std::abs(damage_shift_5));
    if (ctrl_shift_5 >= threshold_5 && damage_shift_5 > 0.01f) {
        profile.composition_bias_5prime = true;
    }

    float damage_shift_3 = profile.terminal_shift_3prime;
    float ctrl_shift_3 = std::abs(profile.ctrl_shift_3prime);
    float threshold_3 = std::max(0.005f, 0.5f * std::abs(damage_shift_3));
    if (ctrl_shift_3 >= threshold_3 && damage_shift_3 > 0.01f) {
        profile.composition_bias_3prime = true;
    }


}

} // namespace taph
