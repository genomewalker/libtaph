#include "damage_estimation_finalize_helpers.hpp"
#include "damage_estimation_finalize_ctx.hpp"
namespace taph {

void finalize_dmax(SampleDamageProfile& profile, const FinalCtx& ctx) {
    // D_max: joint evidence from Channel A (nucleotide frequencies) and Channel B (stop codon conversion).
    // Channel B is the independent validator: T/(T+C) elevation can come from composition OR damage,
    // but stop conversions in CAA/CAG/CGA contexts can only come from real C→T damage.
    {
        // Scan positions 1-5 for peak damage rate when adapter artifact or inversion detected.
        // Using damage_rate[0] fails when pos 0 carries adapter artifact (below baseline);
        // using damage_rate[fit_offset] fails when BIC-best offset points to a below-baseline pos.
        float raw_d_max_5prime, raw_d_max_3prime;
        if (profile.position_0_artifact_5prime || profile.inverted_pattern_5prime
            || profile.briggs_pos0_masked_5prime) {
            float peak = 0.0f;
            for (int p = 1; p <= 5 && p < 15; ++p) {
                if (profile.damage_rate_5prime[p] > peak) peak = profile.damage_rate_5prime[p];
            }
            raw_d_max_5prime = peak;
        } else {
            raw_d_max_5prime = std::clamp(profile.damage_rate_5prime[0], 0.0f, 1.0f);
        }
        // SS damage is C→T at BOTH ends, so the 3' channel is T/(T+C), not A/(A+G).
        // damage_rate_3prime[] is built from the G→A channel (finalize_context) which is
        // the DS 3' signal; for SS recompute the 3' rate from raw t_freq/tc_total_3prime
        // counts (still raw here) against the interior T/(T+C) baseline, mirroring the
        // GA normalization (rate = max(0,(p-b)/(1-b))) and the SS damage_status block below.
        const bool is_ss_3p =
            (profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) ||
            (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
        auto ss_ct3_rate = [&](int p) -> float {
            double n = profile.tc_total_3prime[p];
            if (n <= 0.0) return 0.0f;
            double t_frac = profile.t_freq_3prime[p] / n;
            double denom  = 1.0 - ctx.baseline_tc;
            if (denom <= 0.0) return 0.0f;
            return std::clamp(static_cast<float>((t_frac - ctx.baseline_tc) / denom), 0.0f, 1.0f);
        };
        if (profile.position_0_artifact_3prime || profile.inverted_pattern_3prime
            || profile.briggs_pos0_masked_3prime) {
            float peak = 0.0f;
            for (int p = 1; p <= 5 && p < 15; ++p) {
                float r = is_ss_3p ? ss_ct3_rate(p) : profile.damage_rate_3prime[p];
                if (r > peak) peak = r;
            }
            raw_d_max_3prime = peak;
        } else {
            raw_d_max_3prime = is_ss_3p
                ? ss_ct3_rate(0)
                : std::clamp(profile.damage_rate_3prime[0], 0.0f, 1.0f);
        }


        float d_sum = raw_d_max_5prime + raw_d_max_3prime;
        // Require d_sum>0.04 (4x old 0.01 floor) so the denominator is not near-zero noise.
        // No per-end floor: d5=0.06/d3=0.005 is genuine one-sided asymmetry, not noise.
        if (d_sum > 0.04f) {
            profile.asymmetry = std::abs(raw_d_max_5prime - raw_d_max_3prime) / (d_sum / 2.0f);
            profile.high_asymmetry = (profile.asymmetry > 0.5f);
        } else {
            profile.asymmetry = 0.0f;
            profile.high_asymmetry = false;
        }

        {
            const uint64_t MIN_C_SITES = 10000;  // Minimum C sites for valid per-bin estimate
            float weighted_sum = 0.0f;       // CT-only (Channel A) C-site-weighted accumulator
            float weighted_sum_joint = 0.0f; // max(Channel A, Channel B) C-site-weighted accumulator
            float weight_sum = 0.0f;
            float peak_dmax = 0.0f;
            int peak_bin = -1;

            for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
                auto& b = profile.gc_bins[bin];

                b.c_sites = b.c_interior;
                for (int p = 0; p < 15; ++p) {
                    b.c_sites += b.c_counts[p];
                }

                if (b.n_reads < 1000 || b.c_sites < MIN_C_SITES) {
                    continue;  // Skip bins with insufficient data
                }

                double t_baseline = static_cast<double>(b.t_interior) /
                                   (b.t_interior + b.c_interior + 1);
                double c_baseline = 1.0 - t_baseline;
                double t_terminal = static_cast<double>(b.t_counts[0]) /
                                   (b.t_counts[0] + b.c_counts[0] + 1);

                // For SS libraries the damage signal also appears as C→T at the 3' end;
                // take the max of the two terminal rates to avoid losing signal when the
                // 3' end is dominant (as on typical ssDNA ancient libraries).
                const bool is_ss_bin = (profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) ||
                                       (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
                if (is_ss_bin) {
                    double t_terminal_3 = static_cast<double>(b.t_counts_3prime[0]) /
                                          (b.t_counts_3prime[0] + b.c_counts_3prime[0] + 1);
                    t_terminal = std::max(t_terminal, t_terminal_3);
                }

                if (c_baseline > 0.1) {
                    b.d_max = std::clamp(static_cast<float>((t_terminal - t_baseline) / c_baseline),
                                         0.0f, 1.0f);
                }

                double stop_baseline = static_cast<double>(b.stop_interior) /
                                      (b.stop_interior + b.pre_interior + 1);
                double stop_terminal = static_cast<double>(b.stop_counts[0]) /
                                      (b.stop_counts[0] + b.pre_counts[0] + 1);

                if (stop_baseline < 0.99) {
                    b.d_max_channel_b = std::clamp(
                        static_cast<float>((stop_terminal - stop_baseline) / (1.0 - stop_baseline)),
                        0.0f, 1.0f);
                }

                b.valid = true;

                float bin_dmax = std::max(b.d_max, b.d_max_channel_b);
                float weight = static_cast<float>(b.c_sites);
                // C5: keep gc_stratified_d_max_weighted as the CT-only (Channel A) rate so
                // source="average" cannot exceed the Channel-A d_max_5prime ceiling; the
                // Channel-B-boosted maximum is exposed separately as the joint accumulator.
                weighted_sum += b.d_max * weight;
                weighted_sum_joint += bin_dmax * weight;
                weight_sum += weight;

                if (bin_dmax > peak_dmax) {
                    peak_dmax = bin_dmax;
                    peak_bin = bin;
                }

                int gc_low = bin * 10;
                int gc_high = gc_low + 10;
            }

            if (weight_sum > 0) {
                profile.gc_stratified_d_max_weighted = weighted_sum / weight_sum;
                profile.gc_stratified_d_max_joint    = weighted_sum_joint / weight_sum;
                profile.gc_stratified_d_max_peak = peak_dmax;
                profile.gc_peak_bin = peak_bin;
                profile.gc_stratified_valid = true;


                {
                    constexpr float LLR_THRESHOLD = 10.0f;
                    constexpr float MIN_DMAX_THRESHOLD = 0.01f;

                    uint64_t total_obs = 0;
                    uint64_t damaged_obs = 0;
                    float damaged_weighted_d = 0.0f;
                    uint64_t damaged_weight = 0;
                    int n_damaged = 0;

                    for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
                        auto& b = profile.gc_bins[bin];
                        if (!b.valid) continue;

                        double baseline_tc = static_cast<double>(b.t_interior) /
                                            std::max(1.0, static_cast<double>(b.t_interior + b.c_interior));
                        b.baseline_tc = static_cast<float>(std::clamp(baseline_tc, 0.01, 0.99));

                        uint64_t n_obs = b.n_terminal_obs();
                        total_obs += n_obs;

                        double ll_damaged = 0.0;    // double: at 1e6-1e9 terminal obs the
                        double ll_undamaged = 0.0;  // LLR>10 decision is below float epsilon
                        float lambda = profile.lambda_5prime;
                        const bool is_ss_bin_llr =
                            (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);

                        for (int p = 0; p < 15; ++p) {
                            float decay = std::exp(-lambda * p);
                            float delta_p = b.d_max * decay;
                            float pi_undamaged = b.baseline_tc;
                            float pi_damaged = b.baseline_tc + (1.0f - b.baseline_tc) * delta_p;

                            pi_undamaged = std::clamp(pi_undamaged, 0.001f, 0.999f);
                            pi_damaged = std::clamp(pi_damaged, 0.001f, 0.999f);

                            double k = static_cast<double>(b.t_counts[p]);
                            double n = static_cast<double>(b.t_counts[p] + b.c_counts[p]);
                            if (is_ss_bin_llr && p < static_cast<int>(b.t_counts_3prime.size())) {
                                k += static_cast<double>(b.t_counts_3prime[p]);
                                n += static_cast<double>(b.t_counts_3prime[p] + b.c_counts_3prime[p]);
                            }
                            if (n > 0) {
                                ll_damaged += k * std::log(static_cast<double>(pi_damaged)) + (n - k) * std::log(1.0 - static_cast<double>(pi_damaged));
                                ll_undamaged += k * std::log(static_cast<double>(pi_undamaged)) + (n - k) * std::log(1.0 - static_cast<double>(pi_undamaged));
                            }
                        }

                        b.llr = static_cast<float>(ll_damaged - ll_undamaged);

                        b.classified_damaged = (b.llr > LLR_THRESHOLD) && (b.d_max > MIN_DMAX_THRESHOLD);

                        // Soft probability via logistic on (LLR - threshold)
                        float llr_centered = b.llr - LLR_THRESHOLD;
                        b.p_damaged = 1.0f / (1.0f + std::exp(-0.5f * llr_centered));

                        if (b.d_max < MIN_DMAX_THRESHOLD) {
                            b.p_damaged = 0.0f;
                        }

                        if (b.classified_damaged) {
                            damaged_obs += n_obs;
                            damaged_weighted_d += b.d_max * static_cast<float>(n_obs);
                            damaged_weight += n_obs;
                            ++n_damaged;
                        }
                    }

                    if (total_obs > 0) {
                        profile.pi_damaged = static_cast<float>(damaged_obs) / static_cast<float>(total_obs);
                    }
                    if (damaged_weight > 0) {
                        profile.d_damaged = damaged_weighted_d / static_cast<float>(damaged_weight);
                    }

                    float pop_weighted_d = 0.0f;
                    for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
                        const auto& b = profile.gc_bins[bin];
                        if (!b.valid) continue;
                        pop_weighted_d += b.d_max * static_cast<float>(b.n_terminal_obs());
                    }
                    if (total_obs > 0) {
                        profile.d_population = pop_weighted_d / static_cast<float>(total_obs);
                    }

                    profile.n_damaged_bins = n_damaged;

                }

                // For SS libraries, damage appears as C→T at both ends; feed the
                // combined terminal counts into the GC mixture fit so that bins
                // with 3'-dominant damage (the typical ssDNA case) contribute.
                const bool is_ss_super = (profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) ||
                                         (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
                std::array<SuperRead, N_GC_BINS> super_reads;
                for (int bin = 0; bin < N_GC_BINS; ++bin) {
                    const auto& b = profile.gc_bins[bin];
                    super_reads[bin].gc_bin = bin;
                    super_reads[bin].c_sites = static_cast<double>(b.c_sites);
                    if (is_ss_super) {
                        super_reads[bin].c_sites +=
                            static_cast<double>(b.c_counts_3prime[0] + b.c_counts_3prime[1] + b.c_counts_3prime[2]);
                    }

                    for (int p = 0; p < N_POSITIONS; ++p) {
                        double k = static_cast<double>(b.t_counts[p]);
                        double n = static_cast<double>(b.t_counts[p] + b.c_counts[p]);
                        if (is_ss_super) {
                            k += static_cast<double>(b.t_counts_3prime[p]);
                            n += static_cast<double>(b.t_counts_3prime[p] + b.c_counts_3prime[p]);
                        }
                        super_reads[bin].k_tc[p] = k;
                        super_reads[bin].n_tc[p] = n;
                    }

                    for (int p = 0; p < N_POSITIONS; ++p) {
                        super_reads[bin].k_ag[p] = static_cast<double>(b.a_counts[p]);
                        super_reads[bin].n_ag[p] = static_cast<double>(b.a_counts[p] + b.g_counts[p]);
                    }

                    for (int p = 0; p < N_POSITIONS; ++p) {
                        super_reads[bin].k_stop[p] = static_cast<double>(b.stop_counts[p]);
                        super_reads[bin].n_stop[p] = static_cast<double>(b.stop_counts[p] + b.pre_counts[p]);
                    }

                    super_reads[bin].k_tc_int = static_cast<double>(b.t_interior);
                    super_reads[bin].n_tc_int = static_cast<double>(b.t_interior + b.c_interior);
                    super_reads[bin].k_stop_int = static_cast<double>(b.stop_interior);
                    super_reads[bin].n_stop_int = static_cast<double>(b.stop_interior + b.pre_interior);

                    super_reads[bin].k_ag_int = static_cast<double>(b.a_interior);
                    super_reads[bin].n_ag_int = static_cast<double>(b.a_interior + b.g_interior);
                }

                auto mixture_result = MixtureDamageModel::fit(super_reads);
                profile.mixture_n_components = mixture_result.n_components;
                profile.mixture_d_population = mixture_result.d_population;
                profile.mixture_d_damaged = mixture_result.d_damaged;
                profile.mixture_d_population_highgc = mixture_result.d_population_highgc;
                profile.mixture_pi_damaged = mixture_result.pi_damaged;
                profile.mixture_bic = mixture_result.bic;
                profile.mixture_converged = mixture_result.converged;
                profile.mixture_identifiable = mixture_result.identifiable;

            }
        }

        profile.d_max_5prime = raw_d_max_5prime;
        profile.d_max_3prime = raw_d_max_3prime;

        if (profile.damage_artifact) {
            profile.d_max_5prime = 0.0f;
            profile.d_max_3prime = 0.0f;
            profile.d_max_combined = 0.0f;
            profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
        } else if (profile.damage_validated) {
            // Channel B is NOT verdict-bearing: its terminal-vs-interior stop-codon rate is not
            // composition-immune (modern metagenomes floor at d_max≈0.05, max 0.247 > real ROCS 0.115),
            // so it cannot authenticate. Authenticate on strand-symmetric Channel A only. A position-0
            // artifact no longer diverts to Channel B; samples with genuine symmetric decay (d5≈d3)
            // fall through to joint_bilateral exactly like FLB03mAds3.
            if (profile.inverted_pattern_3prime && !profile.inverted_pattern_5prime) {
                profile.d_max_combined = raw_d_max_5prime;
                profile.d_max_source = SampleDamageProfile::DmaxSource::FIVE_PRIME_ONLY;
            } else if (profile.inverted_pattern_5prime && !profile.inverted_pattern_3prime) {
                profile.d_max_combined = raw_d_max_3prime;
                profile.d_max_source = SampleDamageProfile::DmaxSource::THREE_PRIME_ONLY;
            // d_population_highgc dropped as a d_max source: it is a GC-reweighted mean of
            // the same per-class δ, not an independent estimate (see mixture_damage_model.hpp).
            } else if (profile.gc_stratified_valid) {
                profile.d_max_combined = profile.gc_stratified_d_max_weighted;
                profile.d_max_source = SampleDamageProfile::DmaxSource::JOINT_BILATERAL;
            } else {
                profile.d_max_combined = profile.joint_delta_max;
                profile.d_max_source = SampleDamageProfile::DmaxSource::JOINT_BILATERAL;
            }
        } else if (profile.joint_model_valid && profile.joint_p_damage > 0.5f) {
            if (profile.gc_stratified_valid) {
                profile.d_max_combined = profile.gc_stratified_d_max_weighted;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            } else {
                profile.d_max_combined = profile.joint_delta_max;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            }
        } else if (profile.joint_model_valid) {
            // For single-stranded libraries Channel B is not applicable (ss damage is G→A at 3'),
            // so fall back to Channel A rather than zeroing d_max.
            const bool is_ss = (profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) ||
                               (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
            if (is_ss) {
                profile.d_max_5prime = raw_d_max_5prime;
                profile.d_max_3prime = raw_d_max_3prime;
                profile.d_max_combined = std::max(raw_d_max_5prime, raw_d_max_3prime);
                profile.d_max_source = SampleDamageProfile::DmaxSource::MAX_SS_ASYMMETRY;
            } else {
                // ds libraries: Channel A's terminal C→T carries compositional false
                // positives UNLESS it is length-coupled. Gate on w_length (the validated
                // length-coupling discriminator: genuine deamination falls with read length,
                // a pervasive compositional artifact is flat). When length-coupled, report
                // the validated 5' terminal rate; otherwise leave d_max undetermined (NONE)
                // rather than fabricate one from the GC mixture — which floors at the null
                // (E[d]>0 at H0) and is dominated by the terminal rate (r=0.76 vs 0.90).
                // Threshold 0.6 (not 0.5): null libraries cluster at w_length≈0.5
                // (SOLUTION_pi_delta_dmax.md), and this is the low-confidence regime
                // (not damage_validated, joint p_damage≤0.5) where 'undetermined' is the
                // safe bias over a possibly-compositional terminal rate.
                if (profile.bulk_damage.w_length > 0.6f && raw_d_max_5prime > 0.01f) {
                    profile.d_max_5prime = raw_d_max_5prime;
                    profile.d_max_3prime = raw_d_max_3prime;
                    profile.d_max_combined = raw_d_max_5prime;
                    profile.d_max_source = SampleDamageProfile::DmaxSource::FIVE_PRIME_ONLY;
                } else {
                    profile.d_max_5prime = 0.0f;
                    profile.d_max_3prime = 0.0f;
                    profile.d_max_combined = 0.0f;
                    profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
                }
            }
        } else {
            profile.d_max_5prime = raw_d_max_5prime;
            profile.d_max_3prime = raw_d_max_3prime;

            if (profile.high_asymmetry) {
                const bool is_ss = (profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) ||
                                   (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
                if (is_ss) {
                    // ss libraries have asymmetric damage by design; use max to capture the damaged end
                    profile.d_max_combined = std::max(profile.d_max_5prime, profile.d_max_3prime);
                    profile.d_max_source = SampleDamageProfile::DmaxSource::MAX_SS_ASYMMETRY;
                } else {
                    // ds libraries: high asymmetry suggests artifact - use conservative min
                    profile.d_max_combined = std::min(profile.d_max_5prime, profile.d_max_3prime);
                    profile.d_max_source = SampleDamageProfile::DmaxSource::MIN_ASYMMETRY;
                }
            } else {
                profile.d_max_combined = (profile.d_max_5prime + profile.d_max_3prime) / 2.0f;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            }
        }

        // Skip inversion fallback when Channel B has validated damage (Channel B is ground truth)
        if (!profile.damage_validated) {
            if (profile.inverted_pattern_5prime && !profile.inverted_pattern_3prime) {
                profile.d_max_combined = profile.d_max_3prime;
                profile.d_max_source = SampleDamageProfile::DmaxSource::THREE_PRIME_ONLY;
            } else if (profile.inverted_pattern_3prime && !profile.inverted_pattern_5prime) {
                profile.d_max_combined = profile.d_max_5prime;
                profile.d_max_source = SampleDamageProfile::DmaxSource::FIVE_PRIME_ONLY;
            } else if (profile.inverted_pattern_5prime && profile.inverted_pattern_3prime) {
                // Both ends inverted: normally zero everything (no reliable Channel A signal).
                // Exception: if adapter offset was detected on either end, the scan-for-peak
                // already found the biological damage in pos 1-5. Preserve those values.
                bool has_adapter_5 = profile.position_0_artifact_5prime;
                bool has_adapter_3 = profile.position_0_artifact_3prime;
                if (!has_adapter_5 && !has_adapter_3) {
                    profile.d_max_5prime = 0.0f;
                    profile.d_max_3prime = 0.0f;
                    profile.d_max_combined = 0.0f;
                    profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
                }
            }
        }

    }

    // Per-position flatness diagnostics.
    // A profile is "flat" when max excess rate over positions 0-4 is < 0.005 AND
    // variance over positions 0-9 is < 2.5e-5 (stddev < 0.005 across all positions).
    // This distinguishes genuine-zero damage from attenuated-but-real decay — the
    // latter always leaves a decaying residual above noise at position 0-2.
    {
        auto compute_flatness = [](const std::array<float,15>& rate)
            -> std::tuple<bool,float,float>
        {
            float maxrate = 0.0f;
            for (int p = 0; p <= 4; ++p) maxrate = std::max(maxrate, rate[p]);
            double mean = 0.0;
            for (int p = 0; p <= 9; ++p) mean += rate[p];
            mean /= 10.0;
            double var = 0.0;
            for (int p = 0; p <= 9; ++p) { double d = rate[p] - mean; var += d * d; }
            var /= 10.0;
            bool flat = (maxrate < 0.005f) && (var < 2.5e-5);
            return {flat, maxrate, static_cast<float>(var)};
        };

        auto [flat5, max5, var5] = compute_flatness(profile.damage_rate_5prime);
        auto [flat3, max3, var3] = compute_flatness(profile.damage_rate_3prime);

        profile.d5_profile_flat       = flat5;
        profile.d3_profile_flat       = flat3;
        profile.d5_max_rate_pos0_4    = max5;
        profile.d3_max_rate_pos0_4    = max3;
        profile.d5_profile_var_pos0_9 = var5;
        profile.d3_profile_var_pos0_9 = var3;
        // Blunting suspected: 5' completely flat while 3' shows real decay
        profile.d5_blunting_suspected = flat5 && !flat3 && (profile.d_max_3prime >= 0.05f);
    }

    // Terminal ->G OVERCALL artifact flag (low-abundance cascade degeneracy guard).
    // Raw G-fraction over the FULL base composition, excess over interior, must form a POSITION-FIXED
    // SPIKE at pos0-1 (large peak, NOT continuing a smooth terminal decay into pos2-4). A positive raw-G
    // excess cannot come from G->A deamination (which depletes G), so this axis is artifact-specific and
    // spares the genuine live 3' G->A signal by construction; the non-decay shape test is a second guard.
    {
        const double bg_tot = profile.baseline_a_freq + profile.baseline_g_freq
                            + profile.baseline_t_freq + profile.baseline_c_freq;
        const double g_int  = bg_tot > 0.0 ? profile.baseline_g_freq / bg_tot : 0.0;

        auto g_excess_spike = [&](const std::array<double,15>& a, const std::array<double,15>& g,
                                  const std::array<double,15>& t, const std::array<double,15>& c)
            -> std::pair<bool,float>
        {
            auto gfrac = [&](int p) -> double {
                const double tot = a[p] + g[p] + t[p] + c[p];
                return tot > 0.0 ? g[p] / tot : 0.0;
            };
            // ARTIFACT-SPECIFICITY BY CONSTRUCTION: genuine terminal deamination NEVER raises raw G at
            // either end — 5' C->T leaves G untouched, 3' G->A DEPLETES G — so a POSITIVE raw-G excess
            // cannot be a real decaying channel and never gates out a genuine G-enriched (G->A) live end.
            // Measured over the FULL ACGT composition (raw snapshots), the ->G overcall is a SHARP 2-cycle
            // terminal spike: FLB57md 5' raw-G is 0.48/0.43 at p0/p1 then decays to interior (~0.27) by p2.
            // (The old "sustained plateau" was an artifact of a unit-mixed denominator that excluded T+C.)
            // Require BOTH terminal positions p0 AND p1 elevated >= THR over interior: that is the genuine
            // 2-color overcall shape and rejects a single-position noise blip (e.g. a lone elevated p0).
            // The decay-shaped vs spike discrimination on the OPPOSITE (eligible) end is carried separately
            // by live_end_decay_lrt; this raw-G flag only marks WHICH end is artifact-dead.
            const double spike = std::min(gfrac(0), gfrac(1)) - g_int;  // both terminal positions elevated
            const bool dead = spike >= ARTIFACT_G_SPIKE_THR;
            return {dead, static_cast<float>(std::max(0.0, spike))};
        };

        // Use the raw ACGT snapshots (finalize_init captures them before normalizing
        // 5' t/c and 3' a/g in place); the live arrays mix counts and fractions here.
        auto [sp5, ex5] = g_excess_spike(profile.a_freq_5prime_raw, profile.g_freq_5prime_raw,
                                         profile.t_freq_5prime_raw, profile.c_freq_5prime_raw);
        auto [sp3, ex3] = g_excess_spike(profile.a_freq_3prime_raw, profile.g_freq_3prime_raw,
                                         profile.t_freq_3prime_raw, profile.c_freq_3prime_raw);
        profile.artifact_overcall_5p = sp5;
        profile.artifact_overcall_3p = sp3;
        profile.artifact_g_excess_5p = ex5;
        profile.artifact_g_excess_3p = ex3;
    }

    // d_metamatch: Channel B-anchored estimate.
    // Formula: d_metamatch = d_global + γ × (d_channel_b - d_global)
    // γ from Channel B LLR; blends toward Channel B for validated samples.
    {
        double weighted_sum = 0.0;
        double weight_sum = 0.0;

        for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
            const auto& b = profile.gc_bins[bin];
            if (!b.valid || b.n_reads < 100) continue;

            float gc_center = (bin + 0.5f) / static_cast<float>(SampleDamageProfile::N_GC_BINS);
            float gc_deviation = gc_center - 0.50f;
            float alignability_weight = std::exp(-gc_deviation * gc_deviation / (2.0f * 0.15f * 0.15f));
            double read_weight = std::sqrt(static_cast<double>(b.n_reads));
            double total_weight = alignability_weight * read_weight;

            weighted_sum += b.d_max * total_weight;
            weight_sum += total_weight;
        }

        if (weight_sum > 0) {
            profile.d_alignability_weighted = static_cast<float>(weighted_sum / weight_sum);
        } else {
            profile.d_alignability_weighted = profile.d_max_combined;
        }

        float channel_b_llr = profile.stop_decay_llr_5prime;
        // γ: sigmoid on LLR/50000 (LLR=0 → γ=0.5, LLR=100000 → γ≈0.99)
        float gamma_raw = 1.0f / (1.0f + std::exp(-channel_b_llr / 50000.0f));

        float d_global = profile.d_max_combined;
        float d_channel_b = profile.d_max_from_channel_b;

        if (!profile.damage_validated || profile.damage_artifact) {
            profile.metamatch_gamma = 0.0f;
            profile.d_metamatch = d_global;
        } else if (profile.channel_b_quantifiable && d_channel_b > 0.01f) {
            // Asymmetric: stronger pull toward Channel B when it's higher than d_global,
            // weaker pull when it's lower (avoid under-estimation).
            if (d_channel_b > d_global) {
                profile.metamatch_gamma = gamma_raw;
            } else {
                profile.metamatch_gamma = 0.3f * gamma_raw;
            }
            profile.d_metamatch = d_global + profile.metamatch_gamma * (d_channel_b - d_global);
        } else {
            profile.metamatch_gamma = 0.5f * gamma_raw;
            profile.d_metamatch = d_global + profile.metamatch_gamma * (profile.d_alignability_weighted - d_global);
        }

        profile.d_metamatch = std::clamp(profile.d_metamatch, 0.0f, 1.0f);

        double alignability_total = 0.0;
        double n_total = 0.0;
        for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
            const auto& b = profile.gc_bins[bin];
            if (b.n_reads == 0) continue;

            float gc_center = (bin + 0.5f) / static_cast<float>(SampleDamageProfile::N_GC_BINS);
            float gc_deviation = gc_center - 0.50f;
            float alignability = std::exp(-gc_deviation * gc_deviation / (2.0f * 0.15f * 0.15f));

            alignability_total += alignability * b.n_reads;
            n_total += b.n_reads;
        }
        profile.mean_alignability = (n_total > 0) ? static_cast<float>(alignability_total / n_total) : 0.5f;
    }

    // Damage status: effect-size based, independent of library type.
    // Must run after d_max_5prime/d_max_3prime are finalized (line ~2356+).
    // t_freq_5prime[p] has been normalized to rate at line ~990, so use directly.
    // tc_total_5prime[p] still holds raw T+C counts, used for SE calculation.
    // Uses overdispersed CI (2× SE inflation) so depth alone doesn't drive PRESENT.
    {
        float dmax = std::max(profile.d_max_5prime, profile.d_max_3prime);
        if (dmax >= 0.02f) {
            auto lower95_for_end = [&](const std::array<double,15>& freq,
                                       const std::array<double,15>& total,
                                       double baseline) -> float {
                float exc_max = 0.0f; float n_used = 1.0f;
                for (int p = 1; p <= 4; ++p) {
                    double n = total[p];
                    if (n < 100.0) continue;
                    float exc = static_cast<float>(freq[p] - baseline);
                    if (exc > exc_max) { exc_max = exc; n_used = static_cast<float>(n); }
                }
                if (exc_max <= 0.0f) return -1.0f;
                float p_hat = exc_max + static_cast<float>(baseline);
                float se_od = std::sqrt(p_hat * (1.0f - p_hat) / n_used) * 2.0f;
                return exc_max - 1.96f * se_od;
            };
            float lb5 = lower95_for_end(profile.t_freq_5prime, profile.tc_total_5prime, ctx.baseline_tc);
            // SS 3' damage is C→T (same channel as 5'); DS 3' damage is G→A.
            // t_freq_3prime holds raw counts so compute rate explicitly for SS.
            float lb3;
            if (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) {
                float exc_max = 0.0f; float n_used = 1.0f;
                for (int p = 1; p <= 4; ++p) {
                    double n = profile.tc_total_3prime[p];
                    if (n < 100.0) continue;
                    float exc = static_cast<float>(
                        (n > 0 ? profile.t_freq_3prime[p] / n : 0.0) - ctx.baseline_tc);
                    if (exc > exc_max) { exc_max = exc; n_used = static_cast<float>(n); }
                }
                if (exc_max <= 0.0f) {
                    lb3 = -1.0f;
                } else {
                    float p_hat = exc_max + static_cast<float>(ctx.baseline_tc);
                    float se_od = std::sqrt(p_hat * (1.0f - p_hat) / n_used) * 2.0f;
                    lb3 = exc_max - 1.96f * se_od;
                }
            } else {
                lb3 = lower95_for_end(profile.a_freq_3prime, profile.ag_total_3prime, ctx.baseline_ag);
            }
            float lower95 = std::max(lb5, lb3);
            profile.damage_status = (lower95 >= 0.01f)
                ? SampleDamageProfile::DamageStatus::PRESENT
                : SampleDamageProfile::DamageStatus::WEAK;
            // C5: record WHICH evidence path set the status so consumers can tell an
            // empirical excess-CI confirmation from a fit-amplitude-only WEAK. When
            // both CI lower bounds are -1 (no empirical terminal excess — e.g. an
            // inverted/inward-shifted pattern) but d_max (fit amplitude) cleared 0.02,
            // the verdict rests on the exponential fit alone.
            if (lower95 >= 0.01f) {
                profile.damage_status_basis = (lb5 >= lb3)
                    ? SampleDamageProfile::DamageStatusBasis::CI_5PRIME
                    : SampleDamageProfile::DamageStatusBasis::CI_3PRIME;
            } else if (lb5 < 0.0f && lb3 < 0.0f) {
                profile.damage_status_basis =
                    SampleDamageProfile::DamageStatusBasis::FIT_AMPLITUDE_ONLY;
            } else {
                profile.damage_status_basis = (lb5 >= lb3)
                    ? SampleDamageProfile::DamageStatusBasis::CI_5PRIME
                    : SampleDamageProfile::DamageStatusBasis::CI_3PRIME;
            }
        } else {
            profile.damage_status = SampleDamageProfile::DamageStatus::ABSENT;
            profile.damage_status_basis = SampleDamageProfile::DamageStatusBasis::NONE;
            // Inversion artifacts zero d_max_5/3prime even when BIC fit ct5~=ga3
            // up to ~0.17. Preserve confident BIC verdicts; only fall back to
            // UNKNOWN when the tournament itself was uncertain.
            bool bic_confident =
                profile.library_p_winner >= 0.95f &&
                profile.library_auto_type != SampleDamageProfile::LibraryType::UNKNOWN;
            if (profile.forced_library_type == SampleDamageProfile::LibraryType::UNKNOWN
                && !bic_confident) {
                profile.library_type = SampleDamageProfile::LibraryType::UNKNOWN;
            }
        }
    }

    // === Math-panel relabel: Channel A's d_max estimand (Corrections 1 & 2) ===
    // d_max = A/(1-b) divides out composition, so it consistently estimates the PRODUCT π_dmg·A_b
    // (damaged-molecule fraction × per-ancient terminal C→T amplitude), NOT per-ancient A_b (which is
    // unidentifiable reference-free). terminal_ct_mixture_amp carries d_max_combined under its true
    // estimand label; numerically identical to d_max_combined (byte-for-byte, no recompute).
    profile.terminal_ct_mixture_amp = profile.d_max_combined;
    // Order-statistic selections (max_ss_asymmetry / min_asymmetry pick one end's d_max) are biased as a
    // point estimate of the mixture amplitude; the per-end d_max_5prime/d_max_3prime are the honest objects.
    profile.terminal_ct_mixture_amp_valid_as_point =
        (profile.d_max_source != SampleDamageProfile::DmaxSource::MAX_SS_ASYMMETRY) &&
        (profile.d_max_source != SampleDamageProfile::DmaxSource::MIN_ASYMMETRY);
    // Correction 2 (audit-corrected): the only per-damaged-read amplitude object is a LOWER BOUND.
    // A_b_true = amp / π_dmg, π_dmg ∈ (0,1] ⇒ A_b_true ∈ [amp, ∞); the upper bound is unidentified
    // reference-free (a finite ceiling would assert a hidden π_dmg ≥ threshold prior). Never form amp/w_damaged.
    profile.per_damaged_A_b_lower = profile.terminal_ct_mixture_amp;
    // EM damaged-component weight gate (documents the attenuation; not divided into the amplitude).
    if (profile.mixture_converged && profile.mixture_identifiable)
        profile.w_damaged_gate = SampleDamageProfile::WDamagedGate::IDENTIFIED;
    else if (profile.mixture_converged)
        profile.w_damaged_gate = SampleDamageProfile::WDamagedGate::UNDETERMINED;
    else
        profile.w_damaged_gate = SampleDamageProfile::WDamagedGate::UNAVAILABLE;
}

} // namespace taph
