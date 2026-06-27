#include "damage_estimation_finalize_helpers.hpp"
#include "damage_estimation_finalize_ctx.hpp"
namespace taph {

void finalize_decay(SampleDamageProfile& profile, FinalCtx& ctx) {
    // Step 2: Fit p(pos) = b + A * exp(-lambda * pos) using middle-of-read baseline
    auto fit_5p = fit_exponential_decay(profile.t_freq_5prime, profile.tc_total_5prime, 0.2f,
                                        static_cast<float>(ctx.baseline_tc));
    profile.fit_baseline_5prime = fit_5p[0];
    profile.fit_amplitude_5prime = fit_5p[1];
    ctx.fit_lambda_5p = fit_5p[2];
    profile.fit_rmse_5prime = fit_5p[3];

    auto fit_3p = fit_exponential_decay(profile.a_freq_3prime, profile.ag_total_3prime, 0.2f,
                                        static_cast<float>(ctx.baseline_ag));
    profile.fit_baseline_3prime = fit_3p[0];
    profile.fit_amplitude_3prime = fit_3p[1];
    ctx.fit_lambda_3p = fit_3p[2];
    profile.fit_rmse_3prime = fit_3p[3];

    // Positive LLR = exponential fits better than constant model (real decay pattern)
    profile.decay_llr_5prime = compute_decay_llr(
        profile.t_freq_5prime, profile.tc_total_5prime,
        profile.fit_baseline_5prime, profile.fit_amplitude_5prime, ctx.fit_lambda_5p);
    profile.decay_llr_3prime = compute_decay_llr(
        profile.a_freq_3prime, profile.ag_total_3prime,
        profile.fit_baseline_3prime, profile.fit_amplitude_3prime, ctx.fit_lambda_3p);

    // Control channel decay LLR: A/(A+G) at 5', T/(T+C) at 3'.
    // If control also shows decay, it's likely composition/trimming artifact, not damage.
    // ctx.ctrl_freq_3p / ctx.ctrl_total_3p are hoisted out of the block so they remain accessible
    // to the library-type BIC classifier further below.
            std::array<double, 15> ctrl_freq_5p  = {};
    std::array<double, 15> ctrl_total_5p = {};
    {

        for (int i = 0; i < 15; ++i) {
            double ag_5p = profile.a_freq_5prime[i] + profile.g_freq_5prime[i];
            ctrl_total_5p[i] = ag_5p;
            ctrl_freq_5p[i] = (ag_5p > 0) ? profile.a_freq_5prime[i] / ag_5p : 0.5;

            double tc_3p = profile.t_freq_3prime[i] + profile.c_freq_3prime[i];
            ctx.ctrl_total_3p[i] = tc_3p;
            ctx.ctrl_freq_3p[i] = (tc_3p > 0) ? profile.t_freq_3prime[i] / tc_3p : 0.5;
            profile.tc_total_3prime[i] = tc_3p;
        }

        double ctrl_baseline_5p = ctx.baseline_ag;
        double ctrl_baseline_3p = ctx.baseline_tc;

        profile.ctrl_decay_llr_5prime = compute_decay_llr(
            ctrl_freq_5p, ctrl_total_5p,
            static_cast<float>(ctrl_baseline_5p), 0.0f, ctx.fit_lambda_5p);
        profile.ctrl_decay_llr_3prime = compute_decay_llr(
            ctx.ctrl_freq_3p, ctx.ctrl_total_3p,
            static_cast<float>(ctrl_baseline_3p), 0.0f, ctx.fit_lambda_3p);

        profile.delta_llr_5prime = profile.decay_llr_5prime - profile.ctrl_decay_llr_5prime;
        profile.delta_llr_3prime = profile.decay_llr_3prime - profile.ctrl_decay_llr_3prime;
    }

    // Channel B: stop codon conversion decay. Real C→T damage must create stops
    // in CAA, CAG, CGA contexts; unlike Channel A, this cannot arise from composition bias.
    {
        double total_pre_interior = profile.convertible_caa_interior +
                                   profile.convertible_cag_interior +
                                   profile.convertible_cga_interior;
        double total_stop_interior = profile.convertible_taa_interior +
                                    profile.convertible_tag_interior +
                                    profile.convertible_tga_interior;
        double total_convertible_interior = total_pre_interior + total_stop_interior;

        if (total_convertible_interior > 100) {
            profile.stop_conversion_rate_baseline = static_cast<float>(
                total_stop_interior / total_convertible_interior);
            profile.channel_b_valid = true;
        }

        std::array<double, 15> stop_rate = {};
        std::array<double, 15> stop_exposure = {};

        for (int p = 0; p < 15; ++p) {
            double pre = profile.convertible_caa_5prime[p] +
                        profile.convertible_cag_5prime[p] +
                        profile.convertible_cga_5prime[p];
            double stop = profile.convertible_taa_5prime[p] +
                         profile.convertible_tag_5prime[p] +
                         profile.convertible_tga_5prime[p];
            stop_exposure[p] = pre + stop;
            if (stop_exposure[p] > 10) {
                stop_rate[p] = stop / stop_exposure[p];
            } else {
                stop_rate[p] = profile.stop_conversion_rate_baseline;  // Use baseline if no data
            }
        }

        // Local baseline from positions 5-14 (same reads, past damage zone)
        double local_pre = 0.0, local_stop = 0.0;
        for (int p = 5; p < 15; ++p) {
            local_pre += profile.convertible_caa_5prime[p] +
                        profile.convertible_cag_5prime[p] +
                        profile.convertible_cga_5prime[p];
            local_stop += profile.convertible_taa_5prime[p] +
                         profile.convertible_tag_5prime[p] +
                         profile.convertible_tga_5prime[p];
        }
        float local_baseline = (local_pre + local_stop > 100)
            ? static_cast<float>(local_stop / (local_pre + local_stop))
            : profile.stop_conversion_rate_baseline;

        if (profile.channel_b_valid) {
            // D25: use the interior 30+ baseline (stop_conversion_rate_baseline)
            // instead of the local positions 5-14 window. Positions 5-14 are still
            // inside the exponential damage decay for ss / broad-lambda ds libraries
            // (up to ~22% of peak amplitude at pos 5), inflating the local baseline
            // above the true interior rate and producing false channel_b_inverted.
            float baseline_b = profile.stop_conversion_rate_baseline;
            float lambda_b = ctx.fit_lambda_5p;

            // D25: when a pos0 adapter-remnant artifact is flagged, start the
            // amplitude/LLR loops at position 1 (mirrors the p0_tc_5 mask used by
            // channels F/G/H) so the depleted-T/enriched-pos1 composition artifact
            // does not corrupt the amplitude estimate or LLR.
            const int p0_b = profile.position_0_artifact_5prime ? 1 : 0;

            // Amplitude from positions 0-4: include pos 0 since steep decay (AT-rich)
            // concentrates signal there and excluding it causes false negatives
            double sum_excess = 0.0, sum_weight = 0.0;
            for (int i = p0_b; i < 5; ++i) {
                if (stop_exposure[i] > 50) {
                    double weight = std::exp(-lambda_b * i);
                    double excess = stop_rate[i] - baseline_b;
                    sum_excess += stop_exposure[i] * excess / weight;
                    sum_weight += stop_exposure[i];
                }
            }
            float amplitude_b = (sum_weight > 0) ? static_cast<float>(sum_excess / sum_weight) : 0.0f;
            profile.stop_amplitude_5prime = std::max(0.0f, amplitude_b);

            // LLR for exponential vs constant (include pos 0 for same reason as above)
            double ll_exp = 0.0, ll_const = 0.0;
            for (int p = p0_b; p < 10; ++p) {  // D25: skip pos0 under adapter artifact
                if (stop_exposure[p] < 50) continue;

                double n = stop_exposure[p];
                double pre = profile.convertible_caa_5prime[p] +
                            profile.convertible_cag_5prime[p] +
                            profile.convertible_cga_5prime[p];
                double stop = profile.convertible_taa_5prime[p] +
                             profile.convertible_tag_5prime[p] +
                             profile.convertible_tga_5prime[p];
                double k = stop;

                double p_exp = baseline_b + amplitude_b * std::exp(-lambda_b * p);
                p_exp = std::clamp(p_exp, 0.001, 0.999);
                ll_exp += binomial_ll(k, n, p_exp);

                double p_const = std::clamp(static_cast<double>(baseline_b), 0.001, 0.999);
                ll_const += binomial_ll(k, n, p_const);
            }

            profile.stop_decay_llr_5prime = static_cast<float>(ll_exp - ll_const);

            // If amplitude is negative (inverted), negate LLR
            if (amplitude_b < 0) {
                profile.stop_decay_llr_5prime = -std::abs(profile.stop_decay_llr_5prime);
            }
            // C2: raw binomial LLR scales with N; clamp the emitted magnitude to
            // ±kLlrCap (exploratory, not a calibrated chi-squared(1) p-value).
            profile.stop_decay_llr_5prime = std::clamp(
                profile.stop_decay_llr_5prime,
                -SampleDamageProfile::kLlrCap, SampleDamageProfile::kLlrCap);

            // Per-stop-type LLRs (same model, but with type-specific counts)
            {
                // Exposure counts for per-type validity gating
                profile.n_convertible_caa = profile.convertible_caa_interior + profile.convertible_taa_interior;
                profile.n_convertible_cag = profile.convertible_cag_interior + profile.convertible_tag_interior;
                profile.n_convertible_cga = profile.convertible_cga_interior + profile.convertible_tga_interior;
                profile.channel_b_valid_tga = (profile.n_convertible_cga >= 50.0);
                // D24: symmetrical per-type validity for TAA/TAG (only valid_tga
                // existed). Lets the emitter null lrt_taa/lrt_tag on insufficient
                // type-specific interior exposure instead of emitting a 0.0 sentinel.
                profile.channel_b_valid_taa = (profile.n_convertible_caa >= 50.0);
                profile.channel_b_valid_tag = (profile.n_convertible_cag >= 50.0);

                struct StopType {
                    const std::array<double,15>* pre;
                    const std::array<double,15>* stop;
                    float* llr_out;
                };
                StopType types[3] = {
                    { &profile.convertible_caa_5prime, &profile.convertible_taa_5prime, &profile.stop_decay_llr_taa_5prime },
                    { &profile.convertible_cag_5prime, &profile.convertible_tag_5prime, &profile.stop_decay_llr_tag_5prime },
                    { &profile.convertible_cga_5prime, &profile.convertible_tga_5prime, &profile.stop_decay_llr_tga_5prime },
                };
                for (auto& t : types) {
                    // Per-type interior baseline
                    double t_pre_int = 0, t_stop_int = 0;
                    for (int p = 5; p < 15; ++p) {
                        t_pre_int  += (*t.pre)[p];
                        t_stop_int += (*t.stop)[p];
                    }
                    double t_baseline = (t_pre_int + t_stop_int > 10)
                        ? t_stop_int / (t_pre_int + t_stop_int)
                        : static_cast<double>(profile.stop_conversion_rate_baseline);  // D25: interior 30+ baseline, not pos5-14 local

                    // Amplitude from positions 0-4
                    double t_sum_excess = 0, t_sum_weight = 0;
                    for (int i = p0_b; i < 5; ++i) {  // D25: skip pos0 under adapter artifact
                        double t_exp = (*t.pre)[i] + (*t.stop)[i];
                        if (t_exp > 10) {
                            double w = std::exp(-lambda_b * i);
                            t_sum_excess += t_exp * ((*t.stop)[i] / t_exp - t_baseline) / w;
                            t_sum_weight += t_exp;
                        }
                    }
                    double t_amp = (t_sum_weight > 0) ? t_sum_excess / t_sum_weight : 0.0;

                    // LLR
                    double t_ll_exp = 0, t_ll_const = 0;
                    for (int p = 0; p < 10; ++p) {
                        double t_n = (*t.pre)[p] + (*t.stop)[p];
                        if (t_n < 20) continue;
                        double t_k = (*t.stop)[p];
                        double t_p_exp = std::clamp(t_baseline + t_amp * std::exp(-lambda_b * p), 0.001, 0.999);
                        double t_p_const = std::clamp(t_baseline, 0.001, 0.999);
                        t_ll_exp   += binomial_ll(t_k, t_n, t_p_exp);
                        t_ll_const += binomial_ll(t_k, t_n, t_p_const);
                    }
                    // C2: clamp per-type LLR magnitude (scales with N; not chi2(1)).
                    *t.llr_out = std::clamp(
                        static_cast<float>(t_amp < 0 ? -(t_ll_exp - t_ll_const) : (t_ll_exp - t_ll_const)),
                        -SampleDamageProfile::kLlrCap, SampleDamageProfile::kLlrCap);
                }
            }

            {
                // WLS fit: r_p = b0 + (1-b0) * d_max * exp(-λp)
                // Solve y_p = a + c*x_p, then b0 = a, d_max = c / (1 - b0)
                const double lambda = std::clamp(static_cast<double>(ctx.fit_lambda_5p), 0.1, 0.5);
                const int N_POSITIONS = 15;

                double S_w = 0, S_x = 0, S_xx = 0, S_y = 0, S_xy = 0, S_yy = 0;
                double total_exposure = 0;
                int    n_fit_pos = 0;

                for (int p = 0; p < N_POSITIONS; ++p) {
                    double x_p = std::exp(-lambda * p);

                    double stops_p = profile.convertible_taa_5prime[p] +
                                    profile.convertible_tag_5prime[p] +
                                    profile.convertible_tga_5prime[p];
                    double pre_p = profile.convertible_caa_5prime[p] +
                                  profile.convertible_cag_5prime[p] +
                                  profile.convertible_cga_5prime[p];
                    double exposure_p = stops_p + pre_p;

                    if (exposure_p < 100) continue;  // Skip low-coverage positions

                    double y_p = stops_p / exposure_p;
                    double w_p = exposure_p;

                    S_w  += w_p;
                    S_x  += w_p * x_p;
                    S_xx += w_p * x_p * x_p;
                    S_y  += w_p * y_p;
                    S_xy += w_p * x_p * y_p;
                    S_yy += w_p * y_p * y_p;
                    total_exposure += exposure_p;
                    ++n_fit_pos;
                }

                double denom = S_w * S_xx - S_x * S_x;

                if (total_exposure > 1000 && std::abs(denom) > 1e-10) {
                    double c = (S_w * S_xy - S_x * S_y) / denom;
                    double a = (S_y - c * S_x) / S_w;

                    profile.channel_b_slope = static_cast<float>(c);

                    if (c > 0) {
                        double b0 = std::clamp(a, 0.01, 0.99);
                        double d_max_b = std::clamp(c / (1.0 - b0), 0.0, 1.0);

                        // B1: WLS slope SE. Weighted SSE = S_yy - a*S_y - c*S_xy (normal-equation
                        // identity); sigma2 = SSE/(n-2); Var(c) = sigma2*S_w/denom. d_max = c/(1-b0)
                        // so SE(d_max) = SE(c)/(1-b0) (delta method, b0 held). Additive only.
                        double sse    = std::max(0.0, S_yy - a * S_y - c * S_xy);
                        double sigma2 = (n_fit_pos > 2) ? sse / (n_fit_pos - 2) : 0.0;
                        double var_c  = sigma2 * S_w / denom;
                        double se_c   = (var_c > 0.0) ? std::sqrt(var_c) : 0.0;
                        profile.d_max_from_channel_b_se =
                            static_cast<float>(se_c / (1.0 - b0));

                        profile.d_max_from_channel_b = static_cast<float>(d_max_b);
                        profile.channel_b_weight = static_cast<float>(S_w);
                        profile.channel_b_quantifiable = true;
                        profile.channel_b_inverted = false;
                    } else {
                        // Inverted pattern: terminal stops LOWER than baseline
                        profile.d_max_from_channel_b = 0.0f;
                        profile.d_max_from_channel_b_se = 0.0f;
                        profile.channel_b_weight = 0.0f;
                        profile.channel_b_quantifiable = false;
                        profile.channel_b_inverted = true;
                    }
                } else {
                    profile.d_max_from_channel_b = 0.0f;
                    profile.d_max_from_channel_b_se = 0.0f;
                    profile.channel_b_weight = 0.0f;
                    profile.channel_b_quantifiable = false;
                    profile.channel_b_slope = 0.0f;
                }
            }

            (void)profile.delta_llr_5prime;
            (void)profile.stop_decay_llr_5prime;
        }
    }


}

namespace {
// Pooled difference-in-differences over 16 flanking contexts for one read end.
// signal_is_pyr: signal = T/(T+C) (C→T), control = A/(A+G) (purine mirror).  Otherwise swapped
// (3' of a ds library, where damage is G→A and the purine ratio is the signal).
// Returns false if no context had >=MIN_CONV convertible coverage on both signal and control.
bool pooled_did(const std::array<uint64_t, 64>& term,
                const std::array<uint64_t, 64>& inte,
                bool signal_is_pyr, float& did_out, float& se_out) {
    constexpr uint64_t MIN_CONV = 50;
    // base codes A=0,C=1,G=2,T=3; trinuc index = prev*16 + mid*4 + next
    const int s_a = signal_is_pyr ? 3 : 0;  // signal "ancient" base: T (pyr) or A (pur)
    const int s_b = signal_is_pyr ? 1 : 2;  // signal "modern" base: C (pyr) or G (pur)
    const int c_a = signal_is_pyr ? 0 : 3;  // control bases (the opposite pair)
    const int c_b = signal_is_pyr ? 2 : 1;
    double num = 0.0, vinv = 0.0;
    auto ratio = [](uint64_t a, uint64_t b, double& p, uint64_t& n) {
        n = a + b; p = n ? static_cast<double>(a) / n : 0.0;
    };
    for (int prev = 0; prev < 4; ++prev) {
        for (int next = 0; next < 4; ++next) {
            const int base = prev * 16 + next;
            auto idx = [&](int mid) { return base + mid * 4; };
            double sp_t, sp_i, cp_t, cp_i; uint64_t nst, nsi, nct, nci;
            ratio(term[idx(s_a)], term[idx(s_b)], sp_t, nst);
            ratio(inte[idx(s_a)], inte[idx(s_b)], sp_i, nsi);
            ratio(term[idx(c_a)], term[idx(c_b)], cp_t, nct);
            ratio(inte[idx(c_a)], inte[idx(c_b)], cp_i, nci);
            if (nst < MIN_CONV || nsi < MIN_CONV || nct < MIN_CONV || nci < MIN_CONV) continue;
            const double did = (sp_t - sp_i) - (cp_t - cp_i);
            const double var = sp_t * (1 - sp_t) / nst + sp_i * (1 - sp_i) / nsi
                             + cp_t * (1 - cp_t) / nct + cp_i * (1 - cp_i) / nci;
            if (var <= 0.0) continue;
            num += did / var; vinv += 1.0 / var;
        }
    }
    if (vinv <= 0.0) return false;
    did_out = static_cast<float>(num / vinv);
    se_out  = static_cast<float>(std::sqrt(1.0 / vinv));
    return true;
}
}  // namespace

void compute_codon_did(SampleDamageProfile& profile) {
    const bool is_ss = (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
    float d5 = 0.0f, se5 = 0.0f, d3 = 0.0f, se3 = 0.0f;
    const bool ok5 = pooled_did(profile.tri_5prime_terminal, profile.tri_5prime_interior,
                                /*signal_is_pyr=*/true, d5, se5);
    // 3' damage is G->A for ds (purine signal) but C->T for ss (pyrimidine signal).
    const bool ok3 = pooled_did(profile.tri_3prime_terminal, profile.tri_3prime_interior,
                                /*signal_is_pyr=*/is_ss, d3, se3);
    profile.did_5prime = d5; profile.did_5prime_se = se5;
    profile.did_3prime = d3; profile.did_3prime_se = se3;
    profile.did_valid = ok5 && ok3;
}

} // namespace taph
