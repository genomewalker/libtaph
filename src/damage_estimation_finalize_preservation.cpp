#include "damage_estimation_finalize_helpers.hpp"
#include "damage_estimation_finalize_ctx.hpp"
namespace taph {

void finalize_preservation(SampleDamageProfile& profile) {
    // --- Preservation index (evidence × reliability) ---
    {
        static constexpr float EPS = 1e-6f;
        auto sig = [](float x) { return 1.0f / (1.0f + std::exp(-x)); };
        bool is_ss = (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);

        // f5: 5' terminal C→T; half-point at 5% damage (Briggs 2007)
        profile.preservation_f5 = sig((profile.d_max_5prime - 0.05f) / 0.04f);

        // f3: 3' terminal signal; half-point at 5% damage.
        // For ds libraries where d_max_3prime=0 despite real 5' damage, the terminal
        // base is often clipped by the same trimmer artifact that shifts the 5' peak to
        // pos1 — G→A is inward-displaced and absorbed into background. Impute f3 from
        // the 5' signal: ds deamination is symmetric, so f5 is the calibrated estimate.
        bool ds_3prime_censored = !is_ss
            && profile.d_max_5prime > 0.05f
            && profile.d_max_3prime < 0.03f;
        profile.preservation_f3 = ds_3prime_censored
            ? profile.preservation_f5
            : sig((profile.d_max_3prime - 0.05f) / 0.04f);

        // f_coh: mixture coherence — ancient subpop identifiable with real damage
        {
            // Length-coupling coherence from the VALIDATED bulk law (d_max + w_length),
            // not the reference-free-non-identifiable mixture (which floors at H0 and
            // leaked ~0.87-authentic at a true null). A genuine ancient subpopulation
            // carries terminal damage (d_auth) whose mass DECREASES with read length
            // (w_length high); w_length is the load-bearing discriminator — nulls cluster
            // at 0.5, positives push toward 1 (SOLUTION_pi_delta_dmax.md §6). d_max here is
            // already w_length-gated by finalize_dmax, so at a null both factors → ~0.
            float d_auth = sig((profile.d_max_combined - 0.05f) / 0.03f);
            float w_coupling = std::clamp(
                (static_cast<float>(profile.bulk_damage.w_length) - 0.5f) / 0.4f, 0.0f, 1.0f);
            float coh_signal = d_auth * w_coupling;
            // ds: penalise end-asymmetry, but use d_max_combined as floor so that
            // samples where bulk estimator zeroes d3 (noise/rescue limitation) are not
            // wrongly collapsed to near-zero symmetry.
            float d5_sym = std::max(profile.d_max_5prime, profile.d_max_combined);
            // When 3' is censored (trimmer artifact), treat symmetry as perfect —
            // the signal is present but displaced inward, not absent.
            float d3_sym = ds_3prime_censored ? profile.d_max_5prime
                         : std::max(profile.d_max_3prime, profile.d_max_combined);
            float sym = (is_ss || ds_3prime_censored) ? 1.0f
                : std::exp(-std::abs(std::log((d5_sym + EPS)
                                             / (d3_sym + EPS))) / 0.7f);
            profile.preservation_f_coh = std::sqrt(std::max(coh_signal, EPS) * std::max(sym, EPS));
        }

        // f_cpg: CpG age-bias; log2(CpG_dmax/nonCpG_dmax) > 1 = methyl-C deamination enriched.
        // Near-saturation (d_max_combined > 0.20) both contexts saturate and the ratio collapses
        // toward zero — signal is uninformative, fall back to neutral prior.
        if (!std::isnan(profile.log2_cpg_ratio) && profile.dmax_ct5_noncpg_like > 0.01f
                && profile.d_max_combined < 0.20f) {
            // Floor at uninformative prior: absent/inverted CpG signal is not evidence
            // against antiquity (ss libraries, saturated damage, CpG-depleted taxa).
            profile.preservation_f_cpg = std::max(sig((profile.log2_cpg_ratio - 1.0f) / 0.6f), 0.3f);
            profile.cpg_ratio_evaluable = true;  // D21: computed from the CpG ratio
        } else {
            profile.preservation_f_cpg = 0.3f;  // uninformative prior
            profile.cpg_ratio_evaluable = false;  // D21: 0.3 prior, not an evaluated ratio
        }

        // Weighted geometric mean of 4 factors
        float w5, w3, w_coh, w_cpg;
        if (is_ss) { w5 = 0.35f; w3 = 0.20f; w_coh = 0.28f; w_cpg = 0.17f; }
        else        { w5 = 0.27f; w3 = 0.27f; w_coh = 0.28f; w_cpg = 0.18f; }

        profile.preservation_evidence = std::exp(
            w5   * std::log(std::max(profile.preservation_f5,    EPS)) +
            w3   * std::log(std::max(profile.preservation_f3,    EPS)) +
            w_coh* std::log(std::max(profile.preservation_f_coh, EPS)) +
            w_cpg* std::log(std::max(profile.preservation_f_cpg, EPS)));

        // Reliability gates (continuous — no hard cliffs)
        float g_N   = sig((std::log10(static_cast<float>(profile.n_reads) + 1.0f) - 2.7f) / 0.35f);
        // Fit reliability from the VALIDATED bulk GLM (converged + valid), not the
        // reference-free-non-identifiable mixture. The mixture's pure-ancient shortcut
        // (pi>0.90) is dropped — it is exactly the H0 over-confidence we are removing
        // (SOLUTION_pi_delta_dmax.md §6.6). bulk not attempted (<1000 reads) → both false
        // → g_fit=0.15, consistent with g_N already gating low-read libraries.
        float g_fit = (profile.bulk_damage.converged && profile.bulk_damage.valid) ? 1.0f
                    : (profile.bulk_damage.converged || profile.bulk_damage.valid)  ? 0.5f
                    : 0.15f;
        float g_ox  = 1.0f;
        if (profile.ox_is_artifact)
            g_ox = (profile.ox_d_max >= profile.d_max_combined) ? 0.1f : 0.5f;

        profile.preservation_reliability = g_N * g_fit * g_ox;
        profile.preservation_score = std::clamp(
            profile.preservation_evidence * profile.preservation_reliability, 0.0f, 1.0f);

        // Categorical label
        using PL = SampleDamageProfile::PreservationLabel;
        if (g_N < 0.3f)
            profile.preservation_label = PL::INSUFFICIENT;
        else if (profile.ox_is_artifact)
            profile.preservation_label = PL::ARTIFACT_SUSPECTED;
        else if (profile.preservation_score < 0.15f)
            profile.preservation_label = PL::MODERN_LIKE;
        else if (profile.preservation_score < 0.35f)
            profile.preservation_label = PL::WEAK;
        else if (profile.preservation_score < 0.60f)
            profile.preservation_label = PL::MODERATE;
        else if (profile.preservation_score < 0.80f)
            profile.preservation_label = PL::STRONG;
        else
            profile.preservation_label = PL::EXCEPTIONAL;
    }

    // Lifecycle: mark finalized so re-entry is rejected (see SampleDamageProfile::finalized).
    profile.finalized = true;
}

// Bulk damage law (own phase): threshold-free δ(L). Runs BEFORE finalize_dmax so the
// d_max fallback and preservation authenticity can read w_length (the length-coupling
// discriminator). Aggregates the fine fixed length bins into ~equal-read adaptive bins,
// maps terminal counts to damage/control channels (ss/ds-aware), and fits the
// count-level binomial GLM. Reads raw len_bins + final library_type; writes bulk_damage.
void finalize_bulk(SampleDamageProfile& profile, int bulk_fit_threads) {
    {
        const auto& LB = profile.len_bins;
        uint64_t n_total = 0;
        for (const auto& fb : LB) n_total += fb.n_reads;

        constexpr uint64_t MIN_BULK_READS = 1000;   // below this, no meaningful fit
        constexpr uint64_t READS_PER_BIN  = 2000;   // target reads per adaptive length bin
        constexpr int      MAX_LEN_BINS   = 6;

        if (n_total >= MIN_BULK_READS) {
            // greedy quantile grouping of fine bins (ascending length) into ~equal-read bins
            const int K = static_cast<int>(std::clamp<uint64_t>(
                n_total / READS_PER_BIN, uint64_t(1), uint64_t(MAX_LEN_BINS)));
            const uint64_t per = n_total / static_cast<uint64_t>(K);
            std::vector<std::vector<int>> groups;
            std::vector<int> cur;
            uint64_t cum = 0;
            for (int f = 0; f < SampleDamageProfile::N_LEN_FINE; ++f) {
                if (LB[f].n_reads == 0) continue;
                cur.push_back(f);
                cum += LB[f].n_reads;
                if (cum >= per && static_cast<int>(groups.size()) + 1 < K) {
                    groups.push_back(cur);
                    cur.clear();
                    cum = 0;
                }
            }
            if (!cur.empty()) groups.push_back(cur);

            const int L = static_cast<int>(groups.size());
            const bool is_ss =
                profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED ||
                profile.library_type        == SampleDamageProfile::LibraryType::SINGLE_STRANDED;

            BulkDamageSuffStats bs;
            bs.ss = is_ss;
            bs.skip_3p_pos0 = is_ss;          // ss ligation artifact at the 3' terminal base
            // Wave-3: model the 5' terminus as genuine single-strand overhang (r(0)=1) iff this is an
            // ss library AND p0 is not an adapter/composition artifact (identifiable overhang). Otherwise
            // the ss kernel falls back to the ds form and ss_overhang_degenerate is logged below.
            bs.ss_p0_overhang = is_ss
                             && !profile.position_0_artifact_5prime
                             && !profile.briggs_pos0_masked_5prime;
            bs.bin.resize(L);
            bs.k_interior.resize(L);
            bs.n_interior.resize(L);
            bs.median_len.resize(L);
            bs.jstrat.resize(L);
            std::vector<uint64_t> group_reads(L, 0);
            std::vector<int> group_lo(L, 0), group_hi(L, 0);

            // Short-fragment mode: coverage-weighted interior background borrowed from the
            // long (median_len>=50) groups — the least-damaged reads carry the true library
            // interior C→T rate, which the terminal-only short bins have no clean interior for.
            const bool short_mode = profile.short_fragment_floor < 30;
            uint64_t bgTi = 0, bgCi = 0, bgAi = 0, bgGi = 0;

            for (int l = 0; l < L; ++l) {
                std::array<uint64_t, 15> T5{}, C5{}, A5{}, G5{}, T3{}, C3{}, A3{}, G3{};
                uint64_t Ti = 0, Ci = 0, Ai = 0, Gi = 0, nr = 0, lensum = 0;
                JStrat Jds, Jss;
                for (int f : groups[l]) {
                    const auto& fb = LB[f];
                    for (int p = 0; p < 15; ++p) {
                        T5[p] += fb.t_counts[p];        C5[p] += fb.c_counts[p];
                        A5[p] += fb.a_counts[p];        G5[p] += fb.g_counts[p];
                        T3[p] += fb.t_counts_3prime[p]; C3[p] += fb.c_counts_3prime[p];
                        A3[p] += fb.a_counts_3prime[p]; G3[p] += fb.g_counts_3prime[p];
                    }
                    Jds.add(fb.jstrat_ds);
                    Jss.add(fb.jstrat_ss);
                    Ti += fb.t_interior; Ci += fb.c_interior;
                    Ai += fb.a_interior; Gi += fb.g_interior;
                    nr += fb.n_reads;    lensum += fb.len_sum;
                }
                bs.jstrat[l] = is_ss ? Jss : Jds;   // damage 3' channel resolved here
                group_reads[l] = nr;
                group_lo[l] = SampleDamageProfile::LEN_FINE_MIN +
                              groups[l].front() * SampleDamageProfile::LEN_FINE_W;
                group_hi[l] = SampleDamageProfile::LEN_FINE_MIN +
                              (groups[l].back() + 1) * SampleDamageProfile::LEN_FINE_W;
                bs.median_len[l] = nr ? static_cast<double>(lensum) / static_cast<double>(nr) : 0.0;

                auto& dmg5 = bs.bin[l][0][0]; auto& dmg3 = bs.bin[l][0][1];
                auto& ctl5 = bs.bin[l][1][0]; auto& ctl3 = bs.bin[l][1][1];
                for (int p = 0; p < 15; ++p) {
                    // 5' (identical for ss and ds): damage = C→T (T/(T+C)), control = A/(A+G)
                    dmg5.k[p] = T5[p]; dmg5.n[p] = T5[p] + C5[p];
                    ctl5.k[p] = A5[p]; ctl5.n[p] = A5[p] + G5[p];
                    if (is_ss) {
                        // ss: C→T on the 3' end too; control A/(A+G)
                        dmg3.k[p] = T3[p]; dmg3.n[p] = T3[p] + C3[p];
                        ctl3.k[p] = A3[p]; ctl3.n[p] = A3[p] + G3[p];
                    } else {
                        // ds: G→A on the 3' end (A/(A+G)); control T/(T+C)
                        dmg3.k[p] = A3[p]; dmg3.n[p] = A3[p] + G3[p];
                        ctl3.k[p] = T3[p]; ctl3.n[p] = T3[p] + C3[p];
                    }
                }
                // interior baselines (β warm-start only): damage T/(T+C), control A/(A+G)
                bs.k_interior[l][0] = Ti; bs.n_interior[l][0] = Ti + Ci;
                bs.k_interior[l][1] = Ai; bs.n_interior[l][1] = Ai + Gi;

                if (short_mode && bs.median_len[l] >= 50.0) {
                    bgTi += Ti; bgCi += Ci; bgAi += Ai; bgGi += Gi;
                }
            }

            int Ltot = L;
            // Short-fragment injection: prepend one short group ([16,30) reads) whose terminal
            // counts are real but whose interior baseline is BORROWED from the long groups above.
            // Requires (a) short mode, (b) some short reads, (c) a usable borrowed background — if
            // any is missing we simply skip: the estimator abstains rather than fabricate a bin.
            if (short_mode && (bgTi + bgCi) > 0) {
                std::array<uint64_t, 15> T5{}, C5{}, A5{}, G5{}, T3{}, C3{}, A3{}, G3{};
                uint64_t nr = 0, lensum = 0;
                JStrat Jds, Jss;
                for (int b = 0; b < SampleDamageProfile::N_LEN_SHORT; ++b) {
                    const auto& fb = profile.short_len_bins[b];
                    for (int p = 0; p < 15; ++p) {
                        T5[p] += fb.t_counts[p];        C5[p] += fb.c_counts[p];
                        A5[p] += fb.a_counts[p];        G5[p] += fb.g_counts[p];
                        T3[p] += fb.t_counts_3prime[p]; C3[p] += fb.c_counts_3prime[p];
                        A3[p] += fb.a_counts_3prime[p]; G3[p] += fb.g_counts_3prime[p];
                    }
                    Jds.add(fb.jstrat_ds); Jss.add(fb.jstrat_ss);
                    nr += fb.n_reads; lensum += fb.len_sum;
                }
                if (nr > 0) {
                    std::array<std::array<BulkDamageSuffStats::Cell, 2>, 2> gbin{};
                    for (int p = 0; p < 15; ++p) {
                        gbin[0][0].k[p] = T5[p]; gbin[0][0].n[p] = T5[p] + C5[p];  // 5' damage C→T
                        gbin[1][0].k[p] = A5[p]; gbin[1][0].n[p] = A5[p] + G5[p];  // 5' control A/(A+G)
                        if (is_ss) {
                            gbin[0][1].k[p] = T3[p]; gbin[0][1].n[p] = T3[p] + C3[p];
                            gbin[1][1].k[p] = A3[p]; gbin[1][1].n[p] = A3[p] + G3[p];
                        } else {
                            gbin[0][1].k[p] = A3[p]; gbin[0][1].n[p] = A3[p] + G3[p];  // ds 3' G→A
                            gbin[1][1].k[p] = T3[p]; gbin[1][1].n[p] = T3[p] + C3[p];
                        }
                    }
                    bs.bin.insert(bs.bin.begin(), gbin);
                    bs.short_bin_index = 0;   // prepended ⇒ index 0; fit() profiles THIS bin's π CI
                    bs.k_interior.insert(bs.k_interior.begin(), {bgTi, bgAi});
                    bs.n_interior.insert(bs.n_interior.begin(), {bgTi + bgCi, bgAi + bgGi});
                    bs.median_len.insert(bs.median_len.begin(),
                        static_cast<double>(lensum) / static_cast<double>(nr));
                    bs.jstrat.insert(bs.jstrat.begin(), is_ss ? Jss : Jds);
                    group_reads.insert(group_reads.begin(), nr);
                    group_lo.insert(group_lo.begin(), SampleDamageProfile::LEN_SHORT_ANCHOR + 1);
                    group_hi.insert(group_hi.begin(), SampleDamageProfile::LEN_FINE_MIN);
                    Ltot = static_cast<int>(bs.bin.size());
                }
            }

            profile.bulk_attempted = true;
            BulkDamageResult R = BulkDamageModel::fit(bs, bulk_fit_threads);
            // Wave-3: record whether the bulk kernel modeled the 5' ss overhang (r(0)=1), or fell back
            // to the ds exp form. degenerate = ss library whose terminal overhang was not identifiable.
            profile.ss_overhang_modeled    = R.ss_overhang_modeled;
            profile.ss_overhang_degenerate = is_ss && !R.ss_overhang_modeled;

            // fill the per-bin length extents + read counts the solver does not know
            for (int l = 0; l < Ltot && l < static_cast<int>(R.bins.size()); ++l) {
                R.bins[l].length_lo = group_lo[l];
                R.bins[l].length_hi = group_hi[l];
                R.bins[l].n_reads   = static_cast<int64_t>(group_reads[l]);
            }

            // read-weighted headline δ̂ over the valid (live) length bins
            double num = 0.0, den = 0.0;
            for (int l = 0; l < Ltot && l < static_cast<int>(R.bins.size()); ++l) {
                if (!bs.bin_valid(l)) continue;
                num += R.bins[l].delta * static_cast<double>(group_reads[l]);
                den += static_cast<double>(group_reads[l]);
            }
            profile.bulk_headline_delta = den > 0.0 ? num / den : 0.0;
            profile.bulk_damage = R;
        }
    }
}

} // namespace taph
