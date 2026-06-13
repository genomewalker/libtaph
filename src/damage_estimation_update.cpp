// FrameSelector: accumulation, merge, and reset.
#include "taph/frame_selector_decl.hpp"
#include "taph/codon_tables.hpp"
#include "taph/hexamer_tables.hpp"
#include <algorithm>
#include <cmath>
#include <array>
#include <cstring>
#include <stdexcept>
#include <vector>
namespace taph {
void FrameSelector::update_sample_profile(
    SampleDamageProfile& profile,
    std::string_view seq) {

    if (seq.length() < 30) return;  // Too short for reliable statistics

    size_t len = seq.length();

    // Decode entire read once: uppercase and cache so downstream passes
    // pay one & ~0x20u per base instead of one per scan × 8+ scans.
    thread_local std::vector<char> dec_buf;
    if (dec_buf.size() < len) dec_buf.resize(len);
    char* const decoded = dec_buf.data();
    for (size_t i = 0; i < len; ++i)
        decoded[i] = static_cast<char>(static_cast<unsigned char>(seq[i]) & ~0x20u);

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        char base = decoded[i];
        // Damage signal: T/(T+C) for C→T
        if (base == 'T') {
            profile.t_freq_5prime[i]++;
            profile.tc_total_5prime[i]++;
        } else if (base == 'C') {
            profile.c_freq_5prime[i]++;
            profile.tc_total_5prime[i]++;
        }
        // Negative control: A/(A+G) - should NOT be elevated by C→T damage
        if (base == 'A') {
            profile.a_freq_5prime[i]++;
        } else if (base == 'G') {
            profile.g_freq_5prime[i]++;
        }
    }

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        char base = decoded[pos];
        // Damage signal: A/(A+G) for G→A
        if (base == 'A') {
            profile.a_freq_3prime[i]++;
            profile.ag_total_3prime[i]++;
        } else if (base == 'G') {
            profile.g_freq_3prime[i]++;
            profile.ag_total_3prime[i]++;
        }
        // Negative control: T/(T+C) - should NOT be elevated by G→A damage
        if (base == 'T') {
            profile.t_freq_3prime[i]++;
        } else if (base == 'C') {
            profile.c_freq_3prime[i]++;
        }
    }

    // Tail-anchored background sampling: track C->T (G->A) rates at
    // positions BG_TAIL_LO..BG_TAIL_HI from each terminus to provide a
    // chemistry-robust baseline. Only fills positions actually covered by
    // this read.
    {
        const int lo = SampleDamageProfile::BG_TAIL_LO;
        const int hi = SampleDamageProfile::BG_TAIL_HI;
        for (int i = lo; i <= hi && static_cast<size_t>(i) < len; ++i) {
            const int idx = i - lo;
            const char b5 = decoded[i];
            if (b5 == 'T') { profile.tail_t_5prime[idx]++; profile.tail_tc_5prime[idx]++; }
            else if (b5 == 'C') { profile.tail_tc_5prime[idx]++; }

            const size_t pos3 = len - 1 - i;
            const char b3 = decoded[pos3];
            if (b3 == 'A') { profile.tail_a_3prime[idx]++; profile.tail_ag_3prime[idx]++; }
            else if (b3 == 'G') { profile.tail_ag_3prime[idx]++; }
        }
    }

    // Count bases in middle third (undamaged baseline)
    constexpr size_t INTERIOR_TERM_PAD = 15;
    size_t mid_start = len / 3;
    size_t mid_end   = 2 * len / 3;
    if (mid_start < INTERIOR_TERM_PAD)     mid_start = INTERIOR_TERM_PAD;
    if (len > INTERIOR_TERM_PAD && mid_end + INTERIOR_TERM_PAD > len)
        mid_end = len - INTERIOR_TERM_PAD;
    const bool interior_safe = (mid_start < mid_end);
    if (interior_safe) {
        for (size_t i = mid_start; i < mid_end; ++i) {
            char base = decoded[i];
            if (base == 'T') profile.baseline_t_freq++;
            else if (base == 'C') profile.baseline_c_freq++;
            else if (base == 'A') profile.baseline_a_freq++;
            else if (base == 'G') profile.baseline_g_freq++;
        }
    }

    // Reference-free oxidation-like contrast. During the streaming pass we do
    // not assign ancient/background weights directly. Instead, each read is
    // placed into a terminal-deamination-excess stratum using its own interior
    // composition as a null. Finalization calibrates the strata within length
    // x GC bins and only then compares high-deamination to low-deamination
    // strata. This avoids treating raw terminal T-richness as an ancestry score.
    int read_deam_bin = -1;  // per-read deam stratum (0=modern..4=ancient), hoisted for the
                             // stratified trinucleotide spectrum accumulation further below.
    if (interior_safe) {
        double term_t5 = 0.0, term_tc5 = 0.0;
        for (size_t p = 0; p < std::min<size_t>(5, len); ++p) {
            const char b = decoded[p];
            if (b == 'T') { term_t5 += 1.0; term_tc5 += 1.0; }
            else if (b == 'C') { term_tc5 += 1.0; }
        }

        // 3' position 0 is skipped: SS prep can create a ligation artifact at
        // the final base. Positions 1..4 carry the useful G->A deamination axis.
        double term_a3 = 0.0, term_ag3 = 0.0;
        for (size_t off = 1; off < std::min<size_t>(5, len); ++off) {
            const char b = decoded[len - 1 - off];
            if (b == 'A') { term_a3 += 1.0; term_ag3 += 1.0; }
            else if (b == 'G') { term_ag3 += 1.0; }
        }

        double mid_t = 0.0, mid_c = 0.0, mid_a = 0.0, mid_g = 0.0;
        for (size_t i = mid_start; i < mid_end; ++i) {
            const char b = decoded[i];
            if (b == 'T') ++mid_t;
            else if (b == 'C') ++mid_c;
            else if (b == 'A') ++mid_a;
            else if (b == 'G') ++mid_g;
        }

        const double mid_tc = mid_t + mid_c;
        const double mid_ag = mid_a + mid_g;
        const double mid_total = mid_t + mid_c + mid_a + mid_g;
        if (mid_total > 0.0) {
            double score_num = 0.0, score_den = 0.0;
            if (term_tc5 > 0.0 && mid_tc > 0.0) {
                const double term = (term_t5 + 0.5) / (term_tc5 + 1.0);
                const double base = (mid_t + 0.5) / (mid_tc + 1.0);
                score_num += std::max(0.0, term - base) * term_tc5;
                score_den += term_tc5;
            }
            if (term_ag3 > 0.0 && mid_ag > 0.0) {
                const double term = (term_a3 + 0.5) / (term_ag3 + 1.0);
                const double base = (mid_a + 0.5) / (mid_ag + 1.0);
                score_num += std::max(0.0, term - base) * term_ag3;
                score_den += term_ag3;
            }

            const double deam_score = score_den > 0.0 ? score_num / score_den : 0.0;
            int deam_bin = 0;
            if (deam_score > 0.40) deam_bin = 4;
            else if (deam_score > 0.20) deam_bin = 3;
            else if (deam_score > 0.08) deam_bin = 2;
            else if (deam_score > 0.00) deam_bin = 1;
            read_deam_bin = deam_bin;                       // hoist for stratified trinuc spectrum
            ++profile.deam_stratum_reads[deam_bin];

            const double gc_frac = (mid_g + mid_c) / mid_total;
            const int gc_bin = std::clamp(
                static_cast<int>(gc_frac * SampleDamageProfile::N_OX_GC_BINS),
                0, SampleDamageProfile::N_OX_GC_BINS - 1);
            int len_bin = 3;
            if (len <= 50) len_bin = 0;
            else if (len <= 75) len_bin = 1;
            else if (len <= 100) len_bin = 2;

            auto& oxs = profile.oxidation_like_bins[
                len_bin * SampleDamageProfile::N_OX_GC_BINS + gc_bin].strata[deam_bin];
            ++oxs.reads;
            oxs.term_t5 += term_t5; oxs.term_tc5 += term_tc5;
            oxs.term_a3 += term_a3; oxs.term_ag3 += term_ag3;
            oxs.int_t += mid_t; oxs.int_tc += mid_tc;
            oxs.int_a += mid_a; oxs.int_ag += mid_ag;

            for (size_t i = mid_start; i < mid_end; ++i) {
                const char b = decoded[i];
                if (b == 'T' || b == 'G') {
                    oxs.sig_tg += 1.0;
                    if (b == 'T') oxs.sig_t += 1.0;
                }
                if (b == 'A' || b == 'C') {
                    oxs.sig_ac += 1.0;
                    if (b == 'A') oxs.sig_a += 1.0;
                }
                if (b == 'A' || b == 'T') {
                    oxs.ctrl_at += 1.0;
                    if (b == 'A') oxs.ctrl_a += 1.0;
                }
                if (b == 'C' || b == 'G') {
                    oxs.ctrl_cg += 1.0;
                    if (b == 'C') oxs.ctrl_c += 1.0;
                }
            }

            // Two-marker oxidation bins.
            // s1 = C→T 5' proxy: T count at 5' pos 1-3 (reference-free, base observed directly).
            // s2_ga = G→A 3' proxy (DS): A count at 3' pos 1-3 (minus-strand C→T seen as G→A).
            // s2_ct = C→T 3' proxy (SS): T count at 3' pos 1-3 (same-strand deamination).
            // Library type is not yet resolved here; accumulate both panels unconditionally.
            // compute_oxo_two_marker selects the correct panel at regression time via is_ss.
            {
                int s1 = 0;
                for (int p = 1; p <= 3 && static_cast<size_t>(p) < len; ++p) {
                    if (decoded[p] == 'T') ++s1;
                }
                int s2_ga = 0, s2_ct = 0;
                for (int p = 1; p <= 3 && static_cast<size_t>(p) < len; ++p) {
                    if (decoded[len - 1 - p] == 'A') ++s2_ga;
                    if (decoded[len - 1 - p] == 'T') ++s2_ct;
                }
                s1    = std::min(s1,    3);
                s2_ga = std::min(s2_ga, 3);
                s2_ct = std::min(s2_ct, 3);
                int gc_b  = (mid_total > 0)
                    ? std::min(3, static_cast<int>((mid_g + mid_c) * 4.0 / mid_total))
                    : 1;
                int len_b = (len < 45) ? 0 : (len < 70) ? 1 : (len < 110) ? 2 : 3;
                const uint32_t ngt = static_cast<uint32_t>(mid_t + mid_g);
                const uint32_t t   = static_cast<uint32_t>(mid_t);
                const uint32_t nac = static_cast<uint32_t>(mid_a + mid_c);
                const uint32_t a   = static_cast<uint32_t>(mid_a);
                using Bins = SampleDamageProfile::OxoTwoMarkerBins;
                {
                    auto& cell = profile.oxo_two_marker.cells[Bins::idx(s1, s2_ga, gc_b, len_b)];
                    ++cell.n_reads;
                    cell.sum_nGT += ngt; cell.sum_T += t;
                    cell.sum_nAC += nac; cell.sum_A += a;
                }
                {
                    auto& cell = profile.oxo_two_marker_ss.cells[Bins::idx(s1, s2_ct, gc_b, len_b)];
                    ++cell.n_reads;
                    cell.sum_nGT += ngt; cell.sum_T += t;
                    cell.sum_nAC += nac; cell.sum_A += a;
                }
            }
        }
    }

    // Codon-position-aware counting at 5' end (first 15 bases)
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        int codon_pos = i % 3;  // 0, 1, 2
        char base = decoded[i];
        if (base == 'T') profile.codon_pos_t_count_5prime[codon_pos]++;
        else if (base == 'C') profile.codon_pos_c_count_5prime[codon_pos]++;
    }

    // Codon-position-aware counting at 3' end (last 15 bases)
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        int codon_pos = (len - 1 - i) % 3;
        char base = decoded[pos];
        if (base == 'A') profile.codon_pos_a_count_3prime[codon_pos]++;
        else if (base == 'G') profile.codon_pos_g_count_3prime[codon_pos]++;
    }

    // CpG context damage tracking (5' end, first 5 bases)
    for (size_t i = 0; i < std::min(size_t(5), len - 1); ++i) {
        char base = decoded[i];
        char next = decoded[i + 1];

        // Check for CpG context: C followed by G, or T followed by G (damaged CpG)
        if (next == 'G') {
            if (base == 'C') {
                profile.cpg_c_count++;
            } else if (base == 'T') {
                profile.cpg_t_count++;  // Likely C→T in CpG
            }
        } else {
            // Non-CpG context (same position range as CpG)
            if (base == 'C') {
                profile.non_cpg_c_count++;
            } else if (base == 'T') {
                profile.non_cpg_t_count++;
            }
        }
    }

    // Reference-free trinucleotide spectrum (64 contexts).
    // Terminal zone = read positions 1..4 (both flanks available, inside damage zone).
    // Interior zone = read positions 10..14 (null-distribution baseline).
    // Mirror counters at 3' end use positions counted from the 3' terminus.
    {
        auto nuc_idx = [](char c) -> int {
            switch (c) { case 'A': return 0; case 'C': return 1;
                         case 'G': return 2; case 'T': return 3; }
            return -1;
        };
        auto add_ctx = [&](int prev_pos, int mid_pos, int next_pos,
                           std::array<uint64_t, SampleDamageProfile::N_TRINUC>& target) {
            if (prev_pos < 0 || next_pos >= static_cast<int>(len)) return;
            int i0 = nuc_idx(decoded[prev_pos]);
            int i1 = nuc_idx(decoded[mid_pos]);
            int i2 = nuc_idx(decoded[next_pos]);
            if (i0 < 0 || i1 < 0 || i2 < 0) return;
            ++target[i0 * 16 + i1 * 4 + i2];
        };
        // Stratified variant: increments the bulk target AND the per-read deam_bin stratum, so the
        // context-dependent spectrum is split ancient (high-deam) vs modern (low-deam).
        auto add_ctx_s = [&](int prev_pos, int mid_pos, int next_pos,
                             std::array<uint64_t, SampleDamageProfile::N_TRINUC>& bulk,
                             std::array<std::array<uint64_t, SampleDamageProfile::N_TRINUC>,
                                        SampleDamageProfile::N_OX_DEAM_STRATA>& strat) {
            if (prev_pos < 0 || next_pos >= static_cast<int>(len)) return;
            int i0 = nuc_idx(decoded[prev_pos]);
            int i1 = nuc_idx(decoded[mid_pos]);
            int i2 = nuc_idx(decoded[next_pos]);
            if (i0 < 0 || i1 < 0 || i2 < 0) return;
            const int idx = i0 * 16 + i1 * 4 + i2;
            ++bulk[idx];
            if (read_deam_bin >= 0) ++strat[read_deam_bin][idx];
        };
        for (int p = 1; p <= 4 && p + 1 < static_cast<int>(len); ++p)
            add_ctx_s(p - 1, p, p + 1, profile.tri_5prime_terminal, profile.tri_5prime_terminal_by_deam);
        for (int p = 10; p <= 14 && p + 1 < static_cast<int>(len); ++p)
            add_ctx_s(p - 1, p, p + 1, profile.tri_5prime_interior, profile.tri_5prime_interior_by_deam);
        for (int p = 1; p <= 4; ++p) {
            int mid = static_cast<int>(len) - 1 - p;
            add_ctx_s(mid - 1, mid, mid + 1, profile.tri_3prime_terminal, profile.tri_3prime_terminal_by_deam);
        }
        for (int p = 10; p <= 14; ++p) {
            int mid = static_cast<int>(len) - 1 - p;
            add_ctx_s(mid - 1, mid, mid + 1, profile.tri_3prime_interior, profile.tri_3prime_interior_by_deam);
        }
        // Per-position: all positions 1..N_POS_TRI-1 from each end.
        for (int p = 1; p < SampleDamageProfile::N_POS_TRI && p + 1 < static_cast<int>(len); ++p)
            add_ctx(p - 1, p, p + 1, profile.tri_5prime_pos[p]);
        for (int p = 1; p < SampleDamageProfile::N_POS_TRI; ++p) {
            int mid = static_cast<int>(len) - 1 - p;
            if (mid - 1 >= 0 && mid + 1 < static_cast<int>(len))
                add_ctx(mid - 1, mid, mid + 1, profile.tri_3prime_pos[p]);
        }
    }

    // CpG-like context split — 5' terminal positions (all 15)
    // Also accumulate upstream-context-aware bins (AC, CC, GC, TC)
    for (int p = 0; p < SampleDamageProfile::N_POS && (p + 1) < static_cast<int>(len); ++p) {
        const char x = decoded[p];
        const char y = decoded[p + 1];
        if ((x == 'C' || x == 'T') && (y == 'A' || y == 'C' || y == 'G' || y == 'T')) {
            const int ctx = (y == 'G') ? SampleDamageProfile::CPG_LIKE : SampleDamageProfile::NONCPG_LIKE;
            profile.ct_ctx_total_5prime[ctx][p] += 1.0f;
            if (x == 'T') profile.ct_ctx_t_5prime[ctx][p] += 1.0f;
        }
        // Upstream-context-aware: classify by preceding base (for p > 0)
        if (p > 0 && (x == 'C' || x == 'T')) {
            const char u = decoded[p - 1];
            int uctx = -1;
            switch (u) {
                case 'A': uctx = SampleDamageProfile::CTX_AC; break;
                case 'C': uctx = SampleDamageProfile::CTX_CC; break;
                case 'G': uctx = SampleDamageProfile::CTX_GC; break;
                case 'T': uctx = SampleDamageProfile::CTX_TC; break;
            }
            if (uctx >= 0) {
                profile.ct5_total_by_upstream[uctx][p] += 1.0;
                if (x == 'T') profile.ct5_t_by_upstream[uctx][p] += 1.0;
            }
        }
    }

    // Interior baseline + oxoG 16-context (only for reads >= 30 bp, already guarded above)
    {
        constexpr size_t INTERIOR_TERM_PAD_CTX = 15;
        size_t q0s = len / 3;
        size_t q1s = 2 * len / 3;
        if (q0s < INTERIOR_TERM_PAD_CTX) q0s = INTERIOR_TERM_PAD_CTX;
        if (len > INTERIOR_TERM_PAD_CTX && q1s + INTERIOR_TERM_PAD_CTX > len)
            q1s = len - INTERIOR_TERM_PAD_CTX;
        const bool interior_safe_ctx = (q0s < q1s);
        const int q0 = static_cast<int>(q0s), q1 = static_cast<int>(q1s);

        // Context-split interior baseline
        if (interior_safe_ctx)
        for (int q = q0; q < q1 && (q + 1) < static_cast<int>(len); ++q) {
            const char x = decoded[q], y = decoded[q + 1];
            if ((x == 'C' || x == 'T') && (y == 'A' || y == 'C' || y == 'G' || y == 'T')) {
                const int ctx = (y == 'G') ? SampleDamageProfile::CPG_LIKE : SampleDamageProfile::NONCPG_LIKE;
                profile.ct_ctx_total_interior[ctx] += 1.0f;
                if (x == 'T') profile.ct_ctx_t_interior[ctx] += 1.0f;
            }
            // Upstream-context-aware interior baseline
            if (q > 0 && (x == 'C' || x == 'T')) {
                const char u = decoded[q - 1];
                int uctx = -1;
                switch (u) {
                    case 'A': uctx = SampleDamageProfile::CTX_AC; break;
                    case 'C': uctx = SampleDamageProfile::CTX_CC; break;
                    case 'G': uctx = SampleDamageProfile::CTX_GC; break;
                    case 'T': uctx = SampleDamageProfile::CTX_TC; break;
                }
                if (uctx >= 0) {
                    profile.ct5_total_interior_by_upstream[uctx] += 1.0;
                    if (x == 'T') profile.ct5_t_interior_by_upstream[uctx] += 1.0;
                }
            }
        }

        // oxoG 16-context interior panel
        if (interior_safe_ctx)
        for (int q = q0; q < q1; ++q) {
            if (q <= 0 || q >= static_cast<int>(len) - 1) continue;
            const char l = decoded[q-1], b = decoded[q], r = decoded[q+1];
            auto enc = [](char c) -> int {
                switch(c){ case 'A':return 0; case 'C':return 1; case 'G':return 2; case 'T':return 3; default:return -1; }
            };
            auto rc_base = [](char c) -> char {
                switch(c){ case 'A':return 'T'; case 'T':return 'A'; case 'C':return 'G'; case 'G':return 'C'; default:return 'N'; }
            };
            if (b == 'T') {
                int il = enc(l), ir = enc(r);
                if (il >= 0 && ir >= 0) profile.oxog16_t[4*il+ir] += 1.0f;
            } else if (b == 'A') {
                int il = enc(rc_base(r)), ir = enc(rc_base(l));
                if (il >= 0 && ir >= 0) profile.oxog16_a_rc[4*il+ir] += 1.0f;
            }
        }
    }

    if (len >= 18) {
        for (int frame = 0; frame < 3; ++frame) {
            // Scan codons from 5' end up to position 14
            for (size_t k = 0; ; ++k) {
                size_t codon_start = frame + 3 * k;
                if (codon_start + 3 > len || codon_start > 14) break;

                size_t p = codon_start;
                if (p >= 15) break;

                char b0 = decoded[codon_start];
                char b1 = decoded[codon_start + 1];
                char b2 = decoded[codon_start + 2];

                if ((b0 != 'A' && b0 != 'C' && b0 != 'G' && b0 != 'T') ||
                    (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') ||
                    (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T')) {
                    continue;
                }

                profile.total_codons_5prime[p]++;

                if (b1 == 'A' && b2 == 'A') {
                    if (b0 == 'C') profile.convertible_caa_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_taa_5prime[p]++;
                }
                if (b1 == 'A' && b2 == 'G') {
                    if (b0 == 'C') profile.convertible_cag_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tag_5prime[p]++;
                }
                if (b1 == 'G' && b2 == 'A') {
                    if (b0 == 'C') profile.convertible_cga_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tga_5prime[p]++;
                }
            }
        }

        // Track interior convertible codons (positions 30+ from start)
        // This gives us the baseline stop conversion rate
        // Guard: len >= 63 required to have valid interior region [30, len-30)
        // Without this, len - 30 underflows for short reads causing OOB access
        constexpr size_t INTERIOR_MIN_LEN = 63;  // 30 + 3 + 30
        if (len >= INTERIOR_MIN_LEN) {
            for (int frame = 0; frame < 3; ++frame) {
                for (size_t k = 0; ; ++k) {
                    size_t codon_start = frame + 3 * k;
                    // Interior region: 30 to len-30 (away from both ends)
                    if (codon_start < 30 || codon_start + 3 > len - 30) {
                        if (codon_start + 3 > len - 30) break;
                        continue;
                    }

                    char b0 = decoded[codon_start];
                    char b1 = decoded[codon_start + 1];
                    char b2 = decoded[codon_start + 2];

                    if ((b0 != 'A' && b0 != 'C' && b0 != 'G' && b0 != 'T') ||
                        (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') ||
                        (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T')) {
                        continue;
                    }

                    profile.total_codons_interior++;

                    if (b1 == 'A' && b2 == 'A') {
                        if (b0 == 'C') profile.convertible_caa_interior++;
                        else if (b0 == 'T') profile.convertible_taa_interior++;
                    }
                    if (b1 == 'A' && b2 == 'G') {
                        if (b0 == 'C') profile.convertible_cag_interior++;
                        else if (b0 == 'T') profile.convertible_tag_interior++;
                    }
                    if (b1 == 'G' && b2 == 'A') {
                        if (b0 == 'C') profile.convertible_cga_interior++;
                        else if (b0 == 'T') profile.convertible_tga_interior++;
                    }
                }
            }
        }
    }

    // Prefix hexamer codes for adapter-aware F/G/H terminal bucketing.
    // pfx5 = first 6 bases; pfx3 = last 6 bases.  UINT32_MAX when invalid.
    uint32_t pfx5 = UINT32_MAX, pfx3 = UINT32_MAX;
    if (len >= 6) {
        char hbuf[7]; bool ok = true;
        for (int i = 0; i < 6 && ok; ++i) {
            hbuf[i] = static_cast<char>(static_cast<unsigned char>(seq[i]) & ~0x20u);
            if (hbuf[i]!='A'&&hbuf[i]!='C'&&hbuf[i]!='G'&&hbuf[i]!='T') ok = false;
        }
        if (ok) { hbuf[6]='\0'; pfx5 = encode_hexamer(hbuf); }
        ok = true;
        for (int i = 0; i < 6 && ok; ++i) {
            hbuf[i] = static_cast<char>(static_cast<unsigned char>(seq[len-6+i]) & ~0x20u);
            if (hbuf[i]!='A'&&hbuf[i]!='C'&&hbuf[i]!='G'&&hbuf[i]!='T') ok = false;
        }
        if (ok) { hbuf[6]='\0'; pfx3 = encode_hexamer(hbuf); }
    }

    // Fused 5' terminal codon scan: Channels C + F + G + H
    if (len >= 18) {
        for (int frame = 0; frame < 3; ++frame) {
            for (size_t k = 0; ; ++k) {
                size_t codon_start = frame + 3 * k;
                if (codon_start + 3 > len || codon_start > 14) break;
                size_t p = codon_start;

                const char b0 = decoded[codon_start];
                const char b1 = decoded[codon_start + 1];
                const char b2 = decoded[codon_start + 2];

                if ((b0 != 'A' && b0 != 'C' && b0 != 'G' && b0 != 'T') ||
                    (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') ||
                    (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T')) {
                    continue;
                }

                // Channel C: G→T oxidative stop codons
                if (b1 == 'A' && b2 == 'G') {
                    if (b0 == 'G') profile.convertible_gag_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tag_ox_5prime[p]++;
                }
                if (b1 == 'A' && b2 == 'A') {
                    if (b0 == 'G') profile.convertible_gaa_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_taa_ox_5prime[p]++;
                }
                if (b1 == 'G' && b2 == 'A') {
                    if (b0 == 'G') profile.convertible_gga_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tga_ox_5prime[p]++;
                }
                // Channel F: C→A oxidative stop codons
                if (b0 == 'T' && b2 == 'A') {
                    if (b1 == 'C') profile.convertible_tca_5prime[p]++;
                    else if (b1 == 'A') profile.convertible_taa_ca_5prime[p]++;
                }
                if (b0 == 'T' && b2 == 'G') {
                    if (b1 == 'C') profile.convertible_tcg_5prime[p]++;
                    else if (b1 == 'A') profile.convertible_tag_ca_5prime[p]++;
                }
                if (b0 == 'T' && b1 == 'A') {
                    if (b2 == 'C') profile.convertible_tac_5prime[p]++;
                }
                if (b0 == 'T' && b1 == 'G') {
                    if (b2 == 'C') profile.convertible_tgc_5prime[p]++;
                    else if (b2 == 'A') profile.convertible_tga_ca_5prime[p]++;
                }
                // Channel F: deamination shadows — C→T converts TCA→TTA, TCG→TTG, TAC→TAT, TGC→TGT
                if (b0 == 'T') {
                    const bool sha0 = (b1=='T'&&b2=='A') || (b1=='A'&&b2=='T');
                    const bool sha1 = (b1=='T'&&b2=='G');
                    const bool sha2 = (b1=='G'&&b2=='T');
                    if (sha0 || sha1 || sha2) {
                        profile.ca_deam_shadow_5prime[p]++;
                        if      (sha0) profile.ca_shadow_5prime_ctx0[p]++;
                        else if (sha1) profile.ca_shadow_5prime_ctx1[p]++;
                        else           profile.ca_shadow_5prime_ctx2[p]++;
                    }
                }
                // Channel G: C→G oxidative stop codons
                if (b0 == 'T' && b2 == 'A') {
                    if (b1 == 'C') profile.convertible_tca_cg_5prime[p]++;
                    else if (b1 == 'G') profile.convertible_tga_cg_5prime[p]++;
                }
                if (b0 == 'T' && b1 == 'A') {
                    if (b2 == 'C') profile.convertible_tac_cg_5prime[p]++;
                    else if (b2 == 'G') profile.convertible_tag_cg_5prime[p]++;
                }
                // Channel H: A→T oxidative stop codons
                if (b1 == 'A' && b2 == 'A') {
                    if (b0 == 'A') profile.convertible_aaa_h_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_taa_at_5prime[p]++;
                }
                if (b1 == 'A' && b2 == 'G') {
                    if (b0 == 'A') profile.convertible_aag_h_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tag_at_5prime[p]++;
                }
                if (b1 == 'G' && b2 == 'A') {
                    if (b0 == 'A') profile.convertible_aga_h_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tga_at_5prime[p]++;
                }
                // Prefix-conditioned F/G/H terminal sums (p < 5, keyed by first hexamer)
                if (p < 5 && pfx5 < 4096) {
                    const bool t0 = (b0 == 'T'), a0 = (b0 == 'A');
                    // Channel F (C→A): pre = TCA|TCG|TAC|TGC; stop = TAA|TAG|TGA
                    profile.ca_pre_terminal_by_pfx[pfx5] +=
                        t0 && ((b1=='C'&&(b2=='A'||b2=='G')) || (b1=='A'&&b2=='C') || (b1=='G'&&b2=='C'));
                    profile.ca_stop_terminal_by_pfx[pfx5] +=
                        t0 && ((b2=='A'&&b1=='A') || (b2=='G'&&b1=='A') || (b1=='G'&&b2=='A'));
                    profile.ca_deam_shadow_terminal_by_pfx[pfx5] +=
                        t0 && ((b1=='T'&&(b2=='A'||b2=='G')) || (b1=='A'&&b2=='T') || (b1=='G'&&b2=='T'));
                    // Channel G (C→G): pre = TCA|TAC; stop = TGA|TAG
                    profile.cg_pre_terminal_by_pfx[pfx5] +=
                        t0 && ((b1=='C'&&b2=='A') || (b1=='A'&&b2=='C'));
                    profile.cg_stop_terminal_by_pfx[pfx5] +=
                        t0 && ((b1=='G'&&b2=='A') || (b1=='A'&&b2=='G'));
                    // Channel H (A→T): pre = AAA|AAG|AGA; stop = TAA|TAG|TGA
                    profile.at_pre_terminal_by_pfx[pfx5] +=
                        a0 && ((b1=='A'&&b2=='A') || (b1=='A'&&b2=='G') || (b1=='G'&&b2=='A'));
                    profile.at_stop_terminal_by_pfx[pfx5] +=
                        t0 && ((b1=='A'&&b2=='A') || (b1=='A'&&b2=='G') || (b1=='G'&&b2=='A'));
                    if (p >= 2) {
                        profile.at_pre_terminal_p2plus_by_pfx[pfx5] +=
                            a0 && ((b1=='A'&&b2=='A') || (b1=='A'&&b2=='G') || (b1=='G'&&b2=='A'));
                        profile.at_stop_terminal_p2plus_by_pfx[pfx5] +=
                            t0 && ((b1=='A'&&b2=='A') || (b1=='A'&&b2=='G') || (b1=='G'&&b2=='A'));
                    }
                }
            }
        }
    }

    // Fused 3' terminal codon scan: Channels C + F + G + H
    if (len >= 18) {
        for (int frame = 0; frame < 3; ++frame) {
            for (size_t k = 0; ; ++k) {
                size_t p = k * 3 + frame;
                if (p >= 15) break;
                if (len < 3 + k * 3 + frame) break;
                size_t codon_start = len - 3 - k * 3 - frame;
                if (codon_start + 3 > len) break;

                const char b0 = decoded[codon_start];
                const char b1 = decoded[codon_start + 1];
                const char b2 = decoded[codon_start + 2];

                // Channel C
                if (b1 == 'A' && b2 == 'G') {
                    if (b0 == 'G') profile.convertible_gag_3prime[p]++;
                    else if (b0 == 'T') profile.convertible_tag_ox_3prime[p]++;
                }
                if (b1 == 'A' && b2 == 'A') {
                    if (b0 == 'G') profile.convertible_gaa_3prime[p]++;
                    else if (b0 == 'T') profile.convertible_taa_ox_3prime[p]++;
                }
                if (b1 == 'G' && b2 == 'A') {
                    if (b0 == 'G') profile.convertible_gga_3prime[p]++;
                    else if (b0 == 'T') profile.convertible_tga_ox_3prime[p]++;
                }
                // Channel F
                if (b0 == 'T' && b2 == 'A') {
                    if (b1 == 'C') profile.convertible_tca_3prime[p]++;
                    else if (b1 == 'A') profile.convertible_taa_ca_3prime[p]++;
                }
                if (b0 == 'T' && b2 == 'G') {
                    if (b1 == 'C') profile.convertible_tcg_3prime[p]++;
                    else if (b1 == 'A') profile.convertible_tag_ca_3prime[p]++;
                }
                if (b0 == 'T' && b1 == 'A') {
                    if (b2 == 'C') profile.convertible_tac_3prime[p]++;
                }
                if (b0 == 'T' && b1 == 'G') {
                    if (b2 == 'C') profile.convertible_tgc_3prime[p]++;
                    else if (b2 == 'A') profile.convertible_tga_ca_3prime[p]++;
                }
                // Channel F: deamination shadows at 3' end
                if ((b0=='T'&&b2=='A'&&b1=='T') || (b0=='T'&&b2=='G'&&b1=='T') ||
                    (b0=='T'&&b1=='A'&&b2=='T') || (b0=='T'&&b1=='G'&&b2=='T'))
                    profile.ca_deam_shadow_3prime[p]++;
                // Channel G
                if (b0 == 'T' && b2 == 'A') {
                    if (b1 == 'C') profile.convertible_tca_cg_3prime[p]++;
                    else if (b1 == 'G') profile.convertible_tga_cg_3prime[p]++;
                }
                if (b0 == 'T' && b1 == 'A') {
                    if (b2 == 'C') profile.convertible_tac_cg_3prime[p]++;
                    else if (b2 == 'G') profile.convertible_tag_cg_3prime[p]++;
                }
                // Channel H
                if (b1 == 'A' && b2 == 'A') {
                    if (b0 == 'A') profile.convertible_aaa_h_3prime[p]++;
                    else if (b0 == 'T') profile.convertible_taa_at_3prime[p]++;
                }
                if (b1 == 'A' && b2 == 'G') {
                    if (b0 == 'A') profile.convertible_aag_h_3prime[p]++;
                    else if (b0 == 'T') profile.convertible_tag_at_3prime[p]++;
                }
                if (b1 == 'G' && b2 == 'A') {
                    if (b0 == 'A') profile.convertible_aga_h_3prime[p]++;
                    else if (b0 == 'T') profile.convertible_tga_at_3prime[p]++;
                }
                // Prefix-conditioned F/G/H terminal sums (p < 5, keyed by last hexamer)
                if (p < 5 && pfx3 < 4096) {
                    const bool t0 = (b0 == 'T'), a0 = (b0 == 'A');
                    profile.ca_pre_terminal_3p_by_pfx[pfx3] +=
                        t0 && ((b1=='C'&&(b2=='A'||b2=='G')) || (b1=='A'&&b2=='C') || (b1=='G'&&b2=='C'));
                    profile.ca_stop_terminal_3p_by_pfx[pfx3] +=
                        t0 && ((b2=='A'&&b1=='A') || (b2=='G'&&b1=='A') || (b1=='G'&&b2=='A'));
                    profile.cg_pre_terminal_3p_by_pfx[pfx3] +=
                        t0 && ((b1=='C'&&b2=='A') || (b1=='A'&&b2=='C'));
                    profile.cg_stop_terminal_3p_by_pfx[pfx3] +=
                        t0 && ((b1=='G'&&b2=='A') || (b1=='A'&&b2=='G'));
                    profile.at_pre_terminal_3p_by_pfx[pfx3] +=
                        a0 && ((b1=='A'&&b2=='A') || (b1=='A'&&b2=='G') || (b1=='G'&&b2=='A'));
                    profile.at_stop_terminal_3p_by_pfx[pfx3] +=
                        t0 && ((b1=='A'&&b2=='A') || (b1=='A'&&b2=='G') || (b1=='G'&&b2=='A'));
                }
            }
        }
    }

    // Fused interior baseline: Channels C + F + G + H (positions 30 to len-30)
    if (len >= 63) {
        for (int frame = 0; frame < 3; ++frame) {
            for (size_t k = 0; ; ++k) {
                size_t codon_start = frame + 3 * k;
                if (codon_start + 3 > len - 30) break;
                if (codon_start < 30) continue;

                const char b0 = decoded[codon_start];
                const char b1 = decoded[codon_start + 1];
                const char b2 = decoded[codon_start + 2];

                // Channel C interior
                if (b1 == 'A' && b2 == 'G') {
                    if (b0 == 'G') profile.convertible_gag_interior++;
                    else if (b0 == 'T') profile.convertible_tag_ox_interior++;
                }
                if (b1 == 'A' && b2 == 'A') {
                    if (b0 == 'G') profile.convertible_gaa_interior++;
                    else if (b0 == 'T') profile.convertible_taa_ox_interior++;
                }
                if (b1 == 'G' && b2 == 'A') {
                    if (b0 == 'G') profile.convertible_gga_interior++;
                    else if (b0 == 'T') profile.convertible_tga_ox_interior++;
                }
                // Channel F interior
                if (b0 == 'T' && b2 == 'A') {
                    if      (b1 == 'C') profile.convertible_tca_interior++;
                    else if (b1 == 'A') profile.convertible_taa_ca_interior++;
                }
                if (b0 == 'T' && b2 == 'G') {
                    if      (b1 == 'C') profile.convertible_tcg_interior++;
                    else if (b1 == 'A') profile.convertible_tag_ca_interior++;
                }
                if (b0 == 'T' && b1 == 'A' && b2 == 'C') profile.convertible_tac_interior++;
                if (b0 == 'T' && b1 == 'G') {
                    if      (b2 == 'C') profile.convertible_tgc_interior++;
                    else if (b2 == 'A') profile.convertible_tga_ca_interior++;
                }
                // Channel F far-interior (mirrors convertible_tc*/taa_ca/tag_ca/tga_ca)
                if (b0 == 'T') {
                    bool is_pre    = (b1=='C'&&b2=='A') || (b1=='C'&&b2=='G') ||
                                     (b1=='A'&&b2=='C') || (b1=='G'&&b2=='C');
                    bool is_stop   = (b1=='A'&&b2=='A') || (b1=='A'&&b2=='G') || (b1=='G'&&b2=='A');
                    bool is_shadow = (b1=='T'&&(b2=='A'||b2=='G')) ||
                                     (b1=='A'&&b2=='T') || (b1=='G'&&b2=='T');
                    if      (is_pre)    profile.ca_pre_interior++;
                    else if (is_stop)   profile.ca_stop_interior++;
                    else if (is_shadow) profile.ca_deam_shadow_interior++;
                    // Per-context for MH (ctx0=TCA+TAC, ctx1=TCG, ctx2=TGC)
                    if      ((b1=='C'&&b2=='A') || (b1=='A'&&b2=='C')) profile.ca_pre_interior_by_ctx[0]++;
                    else if (b1=='A'&&b2=='A')                          profile.ca_stop_interior_by_ctx[0]++;
                    else if ((b1=='T'&&b2=='A') || (b1=='A'&&b2=='T'))  profile.ca_shadow_interior_by_ctx[0]++;
                    if      (b1=='C'&&b2=='G') profile.ca_pre_interior_by_ctx[1]++;
                    else if (b1=='A'&&b2=='G') profile.ca_stop_interior_by_ctx[1]++;
                    else if (b1=='T'&&b2=='G') profile.ca_shadow_interior_by_ctx[1]++;
                    if      (b1=='G'&&b2=='C') profile.ca_pre_interior_by_ctx[2]++;
                    else if (b1=='G'&&b2=='A') profile.ca_stop_interior_by_ctx[2]++;
                    else if (b1=='G'&&b2=='T') profile.ca_shadow_interior_by_ctx[2]++;
                }
                // Channel G interior (per-context for MH: ctx0=TCA→TGA, ctx1=TAC→TAG)
                if (b0 == 'T' && b2 == 'A') {
                    if      (b1 == 'C') { profile.cg_pre_interior++;  profile.cg_pre_interior_by_ctx[0]++; }
                    else if (b1 == 'G') { profile.cg_stop_interior++; profile.cg_stop_interior_by_ctx[0]++; }
                }
                if (b0 == 'T' && b1 == 'A') {
                    if      (b2 == 'C') { profile.cg_pre_interior++;  profile.cg_pre_interior_by_ctx[1]++; }
                    else if (b2 == 'G') { profile.cg_stop_interior++; profile.cg_stop_interior_by_ctx[1]++; }
                }
                // Channel H interior (per-context: ctx0=AAA→TAA, ctx1=AAG→TAG, ctx2=AGA→TGA)
                if (b1 == 'A' && b2 == 'A') {
                    if      (b0 == 'A') { profile.at_pre_interior++;  profile.at_pre_interior_by_ctx[0]++; }
                    else if (b0 == 'T') { profile.at_stop_interior++; profile.at_stop_interior_by_ctx[0]++; }
                }
                if (b1 == 'A' && b2 == 'G') {
                    if      (b0 == 'A') { profile.at_pre_interior++;  profile.at_pre_interior_by_ctx[1]++; }
                    else if (b0 == 'T') { profile.at_stop_interior++; profile.at_stop_interior_by_ctx[1]++; }
                }
                if (b1 == 'G' && b2 == 'A') {
                    if      (b0 == 'A') { profile.at_pre_interior++;  profile.at_pre_interior_by_ctx[2]++; }
                    else if (b0 == 'T') { profile.at_stop_interior++; profile.at_stop_interior_by_ctx[2]++; }
                }
            }
        }
    }


    // Channel D: G→T and C→A transversion tracking (8-oxoG, Chargaff-balance cross-check).
    // Accumulate raw G, T, C, A counts at 5' terminal positions (0-14).
    // T/(T+G) and A/(A+C) at interior vs terminal positions detect G→T oxidation without alignment.
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        char base = decoded[i];
        if      (base == 'G') profile.g_count_5prime[i]++;
        else if (base == 'T') profile.t_from_g_5prime[i]++;
        if      (base == 'C') profile.c_count_ox_5prime[i]++;
        else if (base == 'A') profile.a_from_c_5prime[i]++;
    }
    // Interior baseline: T/(T+G) and A/(A+C) in middle third (undamaged reference).
    {
        constexpr size_t INTERIOR_TERM_PAD_D = 15;
        size_t mid_s = len / 3, mid_e = 2 * len / 3;
        if (mid_s < INTERIOR_TERM_PAD_D) mid_s = INTERIOR_TERM_PAD_D;
        if (len > INTERIOR_TERM_PAD_D && mid_e + INTERIOR_TERM_PAD_D > len)
            mid_e = len - INTERIOR_TERM_PAD_D;
        if (mid_s < mid_e) {
            for (size_t i = mid_s; i < mid_e; ++i) {
                char base = decoded[i];
                if      (base == 'G') profile.baseline_g_total++;
                else if (base == 'T') profile.baseline_g_to_t_count++;
                if      (base == 'C') profile.baseline_c_ox_total++;
                else if (base == 'A') profile.baseline_c_to_a_count++;
            }
        }
    }

    // Channel E: purine enrichment tracked via a_freq/g_freq; computed in finalize_sample_profile()

    if (len >= 30) {
        // Compute interior GC from positions 5 to end-5 (avoid both terminal regions)
        size_t interior_start = 5;
        size_t interior_end = len > 10 ? len - 5 : len;
        uint64_t gc_count = 0;
        uint64_t at_count = 0;
        for (size_t i = interior_start; i < interior_end; ++i) {
            char base = decoded[i];
            if (base == 'G' || base == 'C') gc_count++;
            else if (base == 'A' || base == 'T') at_count++;
        }

        if (gc_count + at_count > 0) {
            float gc_frac = static_cast<float>(gc_count) / (gc_count + at_count);
            int bin_idx = std::min(static_cast<int>(gc_frac * SampleDamageProfile::N_GC_BINS),
                                   SampleDamageProfile::N_GC_BINS - 1);
            auto& bin = profile.gc_bins[bin_idx];

            // Accumulate Channel A counts (T and C at terminal positions)
            for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
                char base = decoded[i];
                if (base == 'T') bin.t_counts[i]++;
                else if (base == 'C') bin.c_counts[i]++;
                else if (base == 'A') bin.a_counts[i]++;
                else if (base == 'G') bin.g_counts[i]++;
            }

            // Mirror at 3' end: i=0 is last base, i=1 second-to-last, ...
            // Used to recover C→T damage on 3' end for SS libraries and
            // G→A damage on 3' end for DS libraries at joint-mixture time.
            for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
                char base = decoded[len - 1 - i];
                if (base == 'T') bin.t_counts_3prime[i]++;
                else if (base == 'C') bin.c_counts_3prime[i]++;
                else if (base == 'A') bin.a_counts_3prime[i]++;
                else if (base == 'G') bin.g_counts_3prime[i]++;
            }

            // Accumulate interior baselines for both the signal and control channels.
            constexpr size_t INTERIOR_TERM_PAD_GC = 15;
            size_t mid_start = len / 3;
            size_t mid_end = 2 * len / 3;
            if (mid_start < INTERIOR_TERM_PAD_GC) mid_start = INTERIOR_TERM_PAD_GC;
            if (len > INTERIOR_TERM_PAD_GC && mid_end + INTERIOR_TERM_PAD_GC > len)
                mid_end = len - INTERIOR_TERM_PAD_GC;
            if (mid_start < mid_end) {
                for (size_t i = mid_start; i < mid_end; ++i) {
                    char base = decoded[i];
                    if (base == 'T') bin.t_interior++;
                    else if (base == 'C') bin.c_interior++;
                    else if (base == 'A') bin.a_interior++;
                    else if (base == 'G') bin.g_interior++;
                }
            }

            // Accumulate Channel B counts (stop codons at terminal positions)
            for (int frame = 0; frame < 3; ++frame) {
                for (size_t k = 0; ; ++k) {
                    size_t codon_start = frame + 3 * k;
                    if (codon_start + 3 > len || codon_start > 14) break;
                    size_t p = codon_start;

                    char b0 = decoded[codon_start];
                    char b1 = decoded[codon_start + 1];
                    char b2 = decoded[codon_start + 2];

                    // CAA/TAA, CAG/TAG, CGA/TGA
                    if (b1 == 'A' && b2 == 'A') {
                        if (b0 == 'C') bin.pre_counts[p]++;
                        else if (b0 == 'T') bin.stop_counts[p]++;
                    }
                    if (b1 == 'A' && b2 == 'G') {
                        if (b0 == 'C') bin.pre_counts[p]++;
                        else if (b0 == 'T') bin.stop_counts[p]++;
                    }
                    if (b1 == 'G' && b2 == 'A') {
                        if (b0 == 'C') bin.pre_counts[p]++;
                        else if (b0 == 'T') bin.stop_counts[p]++;
                    }
                }
            }

            // Accumulate Channel B interior baseline
            // Guard: len >= 63 required to prevent unsigned underflow in len - 30
            if (len >= 63) {
                for (int frame = 0; frame < 3; ++frame) {
                    for (size_t k = 0; ; ++k) {
                        size_t codon_start = frame + 3 * k;
                        if (codon_start < 30 || codon_start + 3 > len - 30) {
                            if (codon_start + 3 > len - 30) break;
                            continue;
                        }

                        char b0 = decoded[codon_start];
                        char b1 = decoded[codon_start + 1];
                        char b2 = decoded[codon_start + 2];

                        if (b1 == 'A' && b2 == 'A') {
                            if (b0 == 'C') bin.pre_interior++;
                            else if (b0 == 'T') bin.stop_interior++;
                        }
                        if (b1 == 'A' && b2 == 'G') {
                            if (b0 == 'C') bin.pre_interior++;
                            else if (b0 == 'T') bin.stop_interior++;
                        }
                        if (b1 == 'G' && b2 == 'A') {
                            if (b0 == 'C') bin.pre_interior++;
                            else if (b0 == 'T') bin.stop_interior++;
                        }
                    }
                }
            }

            bin.n_reads++;
        }
    }

    // Per-length-bin terminal counts for the bulk damage law (Phase 1). Fine
    // fixed bins here; finalize_sample_profile aggregates them into data-driven
    // adaptive length bins. Every read contributes (no GC guard) so the bulk law
    // sees the full length distribution.
    {
        auto& lbin = profile.len_bins[SampleDamageProfile::len_fine_bin(len)];
        // per-read eligibility-stratified damage co-occurrence (mixture P2): in the first JP terminal
        // positions count damaged sites k and ELIGIBLE sites S per channel (folded into existing scans).
        constexpr size_t JP = static_cast<size_t>(JStrat::JP);
        int t5 = 0, c5 = 0, t3 = 0, c3 = 0, a3 = 0, g3 = 0;
        for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
            char base = decoded[i];
            if (base == 'T') { lbin.t_counts[i]++; if (i < JP) ++t5; }
            else if (base == 'C') { lbin.c_counts[i]++; if (i < JP) ++c5; }
            else if (base == 'A') lbin.a_counts[i]++;
            else if (base == 'G') lbin.g_counts[i]++;
        }
        for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
            char base = decoded[len - 1 - i];
            if (base == 'T') { lbin.t_counts_3prime[i]++; if (i < JP) ++t3; }
            else if (base == 'C') { lbin.c_counts_3prime[i]++; if (i < JP) ++c3; }
            else if (base == 'A') { lbin.a_counts_3prime[i]++; if (i < JP) ++a3; }
            else if (base == 'G') { lbin.g_counts_3prime[i]++; if (i < JP) ++g3; }
        }
        int S5 = t5 + c5, k5 = t5;          // 5' damage C→T: hits=T, eligible sites=T+C (ds & ss)
        int S3d = a3 + g3, k3d = a3;        // ds 3' damage G→A: hits=A, eligible=A+G
        int S3s = t3 + c3, k3s = t3;        // ss 3' damage C→T: hits=T, eligible=T+C
        lbin.jstrat_ds.n[S5][S3d]++;  lbin.jstrat_ds.sk5[S5][S3d] += k5;
        lbin.jstrat_ds.sk3[S5][S3d] += k3d; lbin.jstrat_ds.sk53[S5][S3d] += static_cast<uint64_t>(k5) * k3d;
        lbin.jstrat_ss.n[S5][S3s]++;  lbin.jstrat_ss.sk5[S5][S3s] += k5;
        lbin.jstrat_ss.sk3[S5][S3s] += k3s; lbin.jstrat_ss.sk53[S5][S3s] += static_cast<uint64_t>(k5) * k3s;
        constexpr size_t INTERIOR_TERM_PAD_LEN = 15;
        size_t mid_start = len / 3, mid_end = 2 * len / 3;
        if (mid_start < INTERIOR_TERM_PAD_LEN) mid_start = INTERIOR_TERM_PAD_LEN;
        if (len > INTERIOR_TERM_PAD_LEN && mid_end + INTERIOR_TERM_PAD_LEN > len)
            mid_end = len - INTERIOR_TERM_PAD_LEN;
        if (mid_start < mid_end) {
            for (size_t i = mid_start; i < mid_end; ++i) {
                char base = decoded[i];
                if (base == 'T') lbin.t_interior++;
                else if (base == 'C') lbin.c_interior++;
                else if (base == 'A') lbin.a_interior++;
                else if (base == 'G') lbin.g_interior++;
            }
        }
        lbin.len_sum += len;
        lbin.n_reads++;
    }

    if (len >= 18) {
        // 5' terminal hexamer starting at position 0
        char hex_5prime[7];
        bool valid_5prime = true;
        for (int i = 0; i < 6; ++i) {
            char b = decoded[i];
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') {
                valid_5prime = false;
                break;
            }
            hex_5prime[i] = b;
        }
        hex_5prime[6] = '\0';

        if (valid_5prime) {
            uint32_t code_5prime = encode_hexamer(hex_5prime);
            if (code_5prime < 4096) {
                profile.hexamer_count_5prime[code_5prime] += 1.0;
                profile.n_hexamers_5prime++;
            }
        }

        // Interior hexamer (starting at len/2 - 3, approximately middle)
        size_t interior_start = len / 2 - 3;
        char hex_interior[7];
        bool valid_interior = true;
        for (int i = 0; i < 6; ++i) {
            char b = decoded[interior_start + i];
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') {
                valid_interior = false;
                break;
            }
            hex_interior[i] = b;
        }
        hex_interior[6] = '\0';

        if (valid_interior) {
            uint32_t code_interior = encode_hexamer(hex_interior);
            if (code_interior < 4096) {
                profile.hexamer_count_interior[code_interior] += 1.0;
                profile.n_hexamers_interior++;
            }
        }

        // 3' terminal hexamer (last 6 bases of read)
        char hex_3prime[7];
        bool valid_3prime = true;
        for (int i = 0; i < 6; ++i) {
            char b = decoded[len - 6 + i];
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') { valid_3prime = false; break; }
            hex_3prime[i] = b;
        }
        hex_3prime[6] = '\0';
        if (valid_3prime) {
            uint32_t code_3prime = encode_hexamer(hex_3prime);
            if (code_3prime < 4096) {
                profile.hexamer_count_3prime[code_3prime] += 1.0;
                profile.n_hexamers_3prime++;
            }
        }
    }

    // Interior clustered C→T: excess co-occurrence of T at non-CpG {C,T} sites
    if (len >= 30) {
        const int lo = static_cast<int>(len / 3);
        const int hi = static_cast<int>(len - (len / 3));
        const int wlen = hi - lo;
        if (wlen >= 2) {
            auto& acc = profile.interior_ct_cluster;
            // Build eligible + indicator arrays for CT and AG tracks.
            // Thread-local scratch reused across reads — eliminates 4 heap
            // allocations per read on this hot path. Buffers grow monotonically
            // and are zeroed in the prefix actually used (size_t wlen).
            thread_local std::vector<uint8_t> ct_elig_buf, ct_pos_buf,
                                              ag_elig_buf, ag_pos_buf;
            const size_t W = static_cast<size_t>(wlen);
            if (ct_elig_buf.size() < W) {
                ct_elig_buf.resize(W);
                ct_pos_buf .resize(W);
                ag_elig_buf.resize(W);
                ag_pos_buf .resize(W);
            }
            uint8_t* ct_elig = ct_elig_buf.data();
            uint8_t* ct_pos  = ct_pos_buf.data();
            uint8_t* ag_elig = ag_elig_buf.data();
            uint8_t* ag_pos  = ag_pos_buf.data();
            std::memset(ct_elig, 0, W);
            std::memset(ct_pos,  0, W);
            std::memset(ag_elig, 0, W);
            std::memset(ag_pos,  0, W);
            int n_ct = 0, k_ct = 0, n_ag = 0, k_ag = 0;
            for (int i = lo; i < hi; ++i) {
                const int j = i - lo;
                const char b = decoded[i];
                const char nx = (i + 1 < static_cast<int>(len)) ? decoded[i+1] : 'N';
                const char pv = (i > 0) ? decoded[i-1] : 'N';
                if ((b == 'C' || b == 'T') && nx != 'G') {
                    ct_elig[j] = 1; ++n_ct;
                    if (b == 'T') { ct_pos[j] = 1; ++k_ct; }
                }
                if ((b == 'A' || b == 'G') && pv != 'C') {
                    ag_elig[j] = 1; ++n_ag;
                    if (b == 'A') { ag_pos[j] = 1; ++k_ag; }
                }
            }
            if (wlen <= 64) {
                // Pack eligible/position flags into uint64 bitmaps; replace
                // the O(W)×10 branchy j-loop with 10 shifted-AND-popcount ops.
                uint64_t ct_elig_w = 0, ct_pos_w = 0;
                uint64_t ag_elig_w = 0, ag_pos_w = 0;
                for (int j = 0; j < wlen; ++j) {
                    ct_elig_w |= static_cast<uint64_t>(ct_elig[j]) << j;
                    ct_pos_w  |= static_cast<uint64_t>(ct_pos[j])  << j;
                    ag_elig_w |= static_cast<uint64_t>(ag_elig[j]) << j;
                    ag_pos_w  |= static_cast<uint64_t>(ag_pos[j])  << j;
                }
                if (n_ct >= 2) {
                    ++acc.reads_used_ct;
                    const double q_ct = (k_ct >= 2)
                        ? (static_cast<double>(k_ct) * (k_ct - 1)) /
                          (static_cast<double>(n_ct) * (n_ct - 1))
                        : 0.0;
                    for (int d = 1; d <= 10; ++d) {
                        uint64_t e     = ct_elig_w & (ct_elig_w >> d);
                        uint64_t pairs = static_cast<uint64_t>(__builtin_popcountll(e));
                        uint64_t obs   = static_cast<uint64_t>(__builtin_popcountll(e & ct_pos_w & (ct_pos_w >> d)));
                        acc.pairs_ct[d] += pairs;
                        acc.obs_ct[d]   += obs;
                        acc.exp_ct[d]   += static_cast<double>(pairs) * q_ct;
                        acc.var_ct[d]   += static_cast<double>(pairs) * q_ct * (1.0 - q_ct);
                    }
                }
                if (n_ag >= 2) {
                    ++acc.reads_used_ag;
                    const double q_ag = (k_ag >= 2)
                        ? (static_cast<double>(k_ag) * (k_ag - 1)) /
                          (static_cast<double>(n_ag) * (n_ag - 1))
                        : 0.0;
                    for (int d = 1; d <= 10; ++d) {
                        uint64_t e     = ag_elig_w & (ag_elig_w >> d);
                        uint64_t pairs = static_cast<uint64_t>(__builtin_popcountll(e));
                        uint64_t obs   = static_cast<uint64_t>(__builtin_popcountll(e & ag_pos_w & (ag_pos_w >> d)));
                        acc.pairs_ag[d] += pairs;
                        acc.obs_ag[d]   += obs;
                        acc.exp_ag[d]   += static_cast<double>(pairs) * q_ag;
                        acc.var_ag[d]   += static_cast<double>(pairs) * q_ag * (1.0 - q_ag);
                    }
                }
            } else {
                if (n_ct >= 2) {
                    ++acc.reads_used_ct;
                    const double q_ct = (k_ct >= 2)
                        ? (static_cast<double>(k_ct) * (k_ct - 1)) /
                          (static_cast<double>(n_ct) * (n_ct - 1))
                        : 0.0;
                    for (int d = 1; d <= 10; ++d) {
                        uint64_t pairs = 0, obs = 0;
                        for (int j = 0; j + d < wlen; ++j) {
                            if (!ct_elig[j] || !ct_elig[j + d]) continue;
                            ++pairs;
                            obs += static_cast<uint64_t>(ct_pos[j] & ct_pos[j + d]);
                        }
                        acc.pairs_ct[d] += pairs;
                        acc.obs_ct[d]   += obs;
                        acc.exp_ct[d]   += static_cast<double>(pairs) * q_ct;
                        acc.var_ct[d]   += static_cast<double>(pairs) * q_ct * (1.0 - q_ct);
                    }
                }
                if (n_ag >= 2) {
                    ++acc.reads_used_ag;
                    const double q_ag = (k_ag >= 2)
                        ? (static_cast<double>(k_ag) * (k_ag - 1)) /
                          (static_cast<double>(n_ag) * (n_ag - 1))
                        : 0.0;
                    for (int d = 1; d <= 10; ++d) {
                        uint64_t pairs = 0, obs = 0;
                        for (int j = 0; j + d < wlen; ++j) {
                            if (!ag_elig[j] || !ag_elig[j + d]) continue;
                            ++pairs;
                            obs += static_cast<uint64_t>(ag_pos[j] & ag_pos[j + d]);
                        }
                        acc.pairs_ag[d] += pairs;
                        acc.obs_ag[d]   += obs;
                        acc.exp_ag[d]   += static_cast<double>(pairs) * q_ag;
                        acc.var_ag[d]   += static_cast<double>(pairs) * q_ag * (1.0 - q_ag);
                    }
                }
            }
        } else {
            ++profile.interior_ct_cluster.short_reads_skipped;
        }
    } else {
        ++profile.interior_ct_cluster.short_reads_skipped;
    }

    profile.n_reads++;
}
void FrameSelector::merge_sample_profiles(SampleDamageProfile& dst, const SampleDamageProfile& src) {
    // Lifecycle guard: merge sums raw count arrays. After finalize they hold
    // rates, so summing them produces nonsense. Reject explicitly rather than
    // corrupt downstream silently.
    if (dst.finalized || src.finalized) {
        throw std::logic_error(
            "merge_sample_profiles: profiles must be merged BEFORE "
            "finalize_sample_profile() is called on either side.");
    }
    for (int i = 0; i < 15; ++i) {
        dst.t_freq_5prime[i] += src.t_freq_5prime[i];
        dst.c_freq_5prime[i] += src.c_freq_5prime[i];
        dst.a_freq_3prime[i] += src.a_freq_3prime[i];
        dst.g_freq_3prime[i] += src.g_freq_3prime[i];
        dst.tc_total_5prime[i] += src.tc_total_5prime[i];
        dst.ag_total_3prime[i] += src.ag_total_3prime[i];
        dst.a_freq_5prime[i] += src.a_freq_5prime[i];
        dst.g_freq_5prime[i] += src.g_freq_5prime[i];
        dst.t_freq_3prime[i] += src.t_freq_3prime[i];
        dst.c_freq_3prime[i] += src.c_freq_3prime[i];
        dst.tc_total_3prime[i] += src.tc_total_3prime[i];
    }
    for (int i = 0; i < SampleDamageProfile::BG_TAIL_N; ++i) {
        dst.tail_t_5prime[i]  += src.tail_t_5prime[i];
        dst.tail_tc_5prime[i] += src.tail_tc_5prime[i];
        dst.tail_a_3prime[i]  += src.tail_a_3prime[i];
        dst.tail_ag_3prime[i] += src.tail_ag_3prime[i];
    }
    dst.pe_short_insert_skipped += src.pe_short_insert_skipped;
    if (src.input_mode == SampleDamageProfile::InputMode::PAIRED) {
        dst.input_mode = SampleDamageProfile::InputMode::PAIRED;
    }

    dst.baseline_t_freq += src.baseline_t_freq;
    dst.baseline_c_freq += src.baseline_c_freq;
    dst.baseline_a_freq += src.baseline_a_freq;
    dst.baseline_g_freq += src.baseline_g_freq;

    for (int p = 0; p < 3; ++p) {
        dst.codon_pos_t_count_5prime[p] += src.codon_pos_t_count_5prime[p];
        dst.codon_pos_c_count_5prime[p] += src.codon_pos_c_count_5prime[p];
        dst.codon_pos_a_count_3prime[p] += src.codon_pos_a_count_3prime[p];
        dst.codon_pos_g_count_3prime[p] += src.codon_pos_g_count_3prime[p];
    }

    dst.cpg_c_count += src.cpg_c_count;
    dst.cpg_t_count += src.cpg_t_count;
    dst.non_cpg_c_count += src.non_cpg_c_count;
    dst.non_cpg_t_count += src.non_cpg_t_count;

    for (int ctx = 0; ctx < SampleDamageProfile::N_CT_CTX; ++ctx) {
        for (int p = 0; p < SampleDamageProfile::N_POS; ++p) {
            dst.ct_ctx_t_5prime[ctx][p] += src.ct_ctx_t_5prime[ctx][p];
            dst.ct_ctx_total_5prime[ctx][p] += src.ct_ctx_total_5prime[ctx][p];
        }
        dst.ct_ctx_t_interior[ctx] += src.ct_ctx_t_interior[ctx];
        dst.ct_ctx_total_interior[ctx] += src.ct_ctx_total_interior[ctx];
    }
    for (int i = 0; i < SampleDamageProfile::N_OXOG16; ++i) {
        dst.oxog16_t[i] += src.oxog16_t[i];
        dst.oxog16_a_rc[i] += src.oxog16_a_rc[i];
    }
    for (int i = 0; i < SampleDamageProfile::N_TRINUC; ++i) {
        dst.tri_5prime_terminal[i] += src.tri_5prime_terminal[i];
        dst.tri_5prime_interior[i] += src.tri_5prime_interior[i];
        dst.tri_3prime_terminal[i] += src.tri_3prime_terminal[i];
        dst.tri_3prime_interior[i] += src.tri_3prime_interior[i];
        for (int s = 0; s < SampleDamageProfile::N_OX_DEAM_STRATA; ++s) {
            dst.tri_5prime_terminal_by_deam[s][i] += src.tri_5prime_terminal_by_deam[s][i];
            dst.tri_5prime_interior_by_deam[s][i] += src.tri_5prime_interior_by_deam[s][i];
            dst.tri_3prime_terminal_by_deam[s][i] += src.tri_3prime_terminal_by_deam[s][i];
            dst.tri_3prime_interior_by_deam[s][i] += src.tri_3prime_interior_by_deam[s][i];
        }
    }
    for (int s = 0; s < SampleDamageProfile::N_OX_DEAM_STRATA; ++s)
        dst.deam_stratum_reads[s] += src.deam_stratum_reads[s];
    for (int p = 0; p < SampleDamageProfile::N_POS_TRI; ++p)
        for (int i = 0; i < SampleDamageProfile::N_TRINUC; ++i) {
            dst.tri_5prime_pos[p][i] += src.tri_5prime_pos[p][i];
            dst.tri_3prime_pos[p][i] += src.tri_3prime_pos[p][i];
        }
    // Merge upstream-context-aware accumulators
    for (int uctx = 0; uctx < SampleDamageProfile::N_UPSTREAM_CTX; ++uctx) {
        for (int p = 0; p < SampleDamageProfile::N_POS; ++p) {
            dst.ct5_t_by_upstream[uctx][p] += src.ct5_t_by_upstream[uctx][p];
            dst.ct5_total_by_upstream[uctx][p] += src.ct5_total_by_upstream[uctx][p];
        }
        dst.ct5_t_interior_by_upstream[uctx] += src.ct5_t_interior_by_upstream[uctx];
        dst.ct5_total_interior_by_upstream[uctx] += src.ct5_total_interior_by_upstream[uctx];
    }

    {
        auto& da = dst.interior_ct_cluster;
        const auto& sa = src.interior_ct_cluster;
        da.reads_used_ct        += sa.reads_used_ct;
        da.reads_used_ag        += sa.reads_used_ag;
        da.short_reads_skipped  += sa.short_reads_skipped;
        for (int d = 1; d <= 10; ++d) {
            da.obs_ct[d]   += sa.obs_ct[d];
            da.pairs_ct[d] += sa.pairs_ct[d];
            da.exp_ct[d]   += sa.exp_ct[d];
            da.var_ct[d]   += sa.var_ct[d];
            da.obs_ag[d]   += sa.obs_ag[d];
            da.pairs_ag[d] += sa.pairs_ag[d];
            da.exp_ag[d]   += sa.exp_ag[d];
            da.var_ag[d]   += sa.var_ag[d];
        }
    }

    for (uint32_t i = 0; i < 4096; ++i) {
        dst.hexamer_count_5prime[i] += src.hexamer_count_5prime[i];
        dst.hexamer_count_interior[i] += src.hexamer_count_interior[i];
        dst.hexamer_count_3prime[i] += src.hexamer_count_3prime[i];
        dst.ca_pre_terminal_by_pfx[i]          += src.ca_pre_terminal_by_pfx[i];
        dst.ca_stop_terminal_by_pfx[i]         += src.ca_stop_terminal_by_pfx[i];
        dst.ca_deam_shadow_terminal_by_pfx[i]  += src.ca_deam_shadow_terminal_by_pfx[i];
        dst.cg_pre_terminal_by_pfx[i]    += src.cg_pre_terminal_by_pfx[i];
        dst.cg_stop_terminal_by_pfx[i]   += src.cg_stop_terminal_by_pfx[i];
        dst.at_pre_terminal_by_pfx[i]         += src.at_pre_terminal_by_pfx[i];
        dst.at_stop_terminal_by_pfx[i]        += src.at_stop_terminal_by_pfx[i];
        dst.at_pre_terminal_p2plus_by_pfx[i]  += src.at_pre_terminal_p2plus_by_pfx[i];
        dst.at_stop_terminal_p2plus_by_pfx[i] += src.at_stop_terminal_p2plus_by_pfx[i];
        dst.ca_pre_terminal_3p_by_pfx[i]  += src.ca_pre_terminal_3p_by_pfx[i];
        dst.ca_stop_terminal_3p_by_pfx[i] += src.ca_stop_terminal_3p_by_pfx[i];
        dst.cg_pre_terminal_3p_by_pfx[i]  += src.cg_pre_terminal_3p_by_pfx[i];
        dst.cg_stop_terminal_3p_by_pfx[i] += src.cg_stop_terminal_3p_by_pfx[i];
        dst.at_pre_terminal_3p_by_pfx[i]  += src.at_pre_terminal_3p_by_pfx[i];
        dst.at_stop_terminal_3p_by_pfx[i] += src.at_stop_terminal_3p_by_pfx[i];
    }
    dst.n_hexamers_5prime += src.n_hexamers_5prime;
    dst.n_hexamers_interior += src.n_hexamers_interior;
    dst.n_hexamers_3prime += src.n_hexamers_3prime;

    for (int i = 0; i < 15; ++i) {
        dst.convertible_caa_5prime[i] += src.convertible_caa_5prime[i];
        dst.convertible_taa_5prime[i] += src.convertible_taa_5prime[i];
        dst.convertible_cag_5prime[i] += src.convertible_cag_5prime[i];
        dst.convertible_tag_5prime[i] += src.convertible_tag_5prime[i];
        dst.convertible_cga_5prime[i] += src.convertible_cga_5prime[i];
        dst.convertible_tga_5prime[i] += src.convertible_tga_5prime[i];
        dst.total_codons_5prime[i] += src.total_codons_5prime[i];
    }
    dst.convertible_caa_interior += src.convertible_caa_interior;
    dst.convertible_taa_interior += src.convertible_taa_interior;
    dst.convertible_cag_interior += src.convertible_cag_interior;
    dst.convertible_tag_interior += src.convertible_tag_interior;
    dst.convertible_cga_interior += src.convertible_cga_interior;
    dst.convertible_tga_interior += src.convertible_tga_interior;
    dst.total_codons_interior += src.total_codons_interior;

    for (int i = 0; i < 15; ++i) {
        dst.convertible_gag_5prime[i] += src.convertible_gag_5prime[i];
        dst.convertible_tca_5prime[i]    += src.convertible_tca_5prime[i];
        dst.convertible_tcg_5prime[i]    += src.convertible_tcg_5prime[i];
        dst.convertible_tac_5prime[i]    += src.convertible_tac_5prime[i];
        dst.convertible_tgc_5prime[i]    += src.convertible_tgc_5prime[i];
        dst.convertible_taa_ca_5prime[i] += src.convertible_taa_ca_5prime[i];
        dst.convertible_tag_ca_5prime[i] += src.convertible_tag_ca_5prime[i];
        dst.convertible_tga_ca_5prime[i]  += src.convertible_tga_ca_5prime[i];
        dst.ca_deam_shadow_5prime[i]      += src.ca_deam_shadow_5prime[i];
        dst.convertible_tca_3prime[i]    += src.convertible_tca_3prime[i];
        dst.convertible_tcg_3prime[i]    += src.convertible_tcg_3prime[i];
        dst.convertible_tac_3prime[i]    += src.convertible_tac_3prime[i];
        dst.convertible_tgc_3prime[i]    += src.convertible_tgc_3prime[i];
        dst.convertible_taa_ca_3prime[i] += src.convertible_taa_ca_3prime[i];
        dst.convertible_tag_ca_3prime[i] += src.convertible_tag_ca_3prime[i];
        dst.convertible_tga_ca_3prime[i]  += src.convertible_tga_ca_3prime[i];
        dst.ca_deam_shadow_3prime[i]      += src.ca_deam_shadow_3prime[i];
        dst.convertible_tca_cg_5prime[i]  += src.convertible_tca_cg_5prime[i];
        dst.convertible_tac_cg_5prime[i]  += src.convertible_tac_cg_5prime[i];
        dst.convertible_tga_cg_5prime[i]  += src.convertible_tga_cg_5prime[i];
        dst.convertible_tag_cg_5prime[i]  += src.convertible_tag_cg_5prime[i];
        dst.convertible_tca_cg_3prime[i]  += src.convertible_tca_cg_3prime[i];
        dst.convertible_tac_cg_3prime[i]  += src.convertible_tac_cg_3prime[i];
        dst.convertible_tga_cg_3prime[i]  += src.convertible_tga_cg_3prime[i];
        dst.convertible_tag_cg_3prime[i]  += src.convertible_tag_cg_3prime[i];
        dst.convertible_aaa_h_5prime[i]   += src.convertible_aaa_h_5prime[i];
        dst.convertible_aag_h_5prime[i]   += src.convertible_aag_h_5prime[i];
        dst.convertible_aga_h_5prime[i]   += src.convertible_aga_h_5prime[i];
        dst.convertible_taa_at_5prime[i]  += src.convertible_taa_at_5prime[i];
        dst.convertible_tag_at_5prime[i]  += src.convertible_tag_at_5prime[i];
        dst.convertible_tga_at_5prime[i]  += src.convertible_tga_at_5prime[i];
        dst.convertible_aaa_h_3prime[i]   += src.convertible_aaa_h_3prime[i];
        dst.convertible_aag_h_3prime[i]   += src.convertible_aag_h_3prime[i];
        dst.convertible_aga_h_3prime[i]   += src.convertible_aga_h_3prime[i];
        dst.convertible_taa_at_3prime[i]  += src.convertible_taa_at_3prime[i];
        dst.convertible_tag_at_3prime[i]  += src.convertible_tag_at_3prime[i];
        dst.convertible_tga_at_3prime[i]  += src.convertible_tga_at_3prime[i];
        // Position-INDEPENDENT interior scalars + k-/p-indexed arrays: accumulate
        // ONCE, not once per i. Previously unguarded inside the i<15 loop, which
        // summed each 15× and inflated every interior Channel-F/G/H statistic on the
        // standard multi-shard merge path (thread-count-dependent wrong values).
        if (i == 0) {
            dst.ca_pre_interior          += src.ca_pre_interior;
            dst.ca_stop_interior         += src.ca_stop_interior;
            dst.ca_deam_shadow_interior  += src.ca_deam_shadow_interior;
            for (int k = 0; k < 3; ++k) {
                dst.ca_pre_interior_by_ctx[k]    += src.ca_pre_interior_by_ctx[k];
                dst.ca_stop_interior_by_ctx[k]   += src.ca_stop_interior_by_ctx[k];
                dst.ca_shadow_interior_by_ctx[k] += src.ca_shadow_interior_by_ctx[k];
            }
            for (int p = 0; p < 15; ++p) {
                dst.ca_shadow_5prime_ctx0[p] += src.ca_shadow_5prime_ctx0[p];
                dst.ca_shadow_5prime_ctx1[p] += src.ca_shadow_5prime_ctx1[p];
                dst.ca_shadow_5prime_ctx2[p] += src.ca_shadow_5prime_ctx2[p];
            }
            dst.cg_pre_interior  += src.cg_pre_interior;
            dst.cg_stop_interior += src.cg_stop_interior;
            dst.at_pre_interior  += src.at_pre_interior;
            dst.at_stop_interior += src.at_stop_interior;
            for (int k = 0; k < 2; ++k) {  // Correction 3: G context strata
                dst.cg_pre_interior_by_ctx[k]  += src.cg_pre_interior_by_ctx[k];
                dst.cg_stop_interior_by_ctx[k] += src.cg_stop_interior_by_ctx[k];
            }
            for (int k = 0; k < 3; ++k) {  // Correction 3: H context strata
                dst.at_pre_interior_by_ctx[k]  += src.at_pre_interior_by_ctx[k];
                dst.at_stop_interior_by_ctx[k] += src.at_stop_interior_by_ctx[k];
            }
        }



        dst.convertible_gag_3prime[i]    += src.convertible_gag_3prime[i];
        dst.convertible_tag_ox_3prime[i] += src.convertible_tag_ox_3prime[i];
        dst.convertible_gaa_3prime[i]    += src.convertible_gaa_3prime[i];
        dst.convertible_taa_ox_3prime[i] += src.convertible_taa_ox_3prime[i];
        dst.convertible_gga_3prime[i]    += src.convertible_gga_3prime[i];
        dst.convertible_tga_ox_3prime[i] += src.convertible_tga_ox_3prime[i];

        dst.convertible_tag_ox_5prime[i] += src.convertible_tag_ox_5prime[i];
        dst.convertible_gaa_5prime[i] += src.convertible_gaa_5prime[i];
        dst.convertible_taa_ox_5prime[i] += src.convertible_taa_ox_5prime[i];
        dst.convertible_gga_5prime[i] += src.convertible_gga_5prime[i];
        dst.convertible_tga_ox_5prime[i] += src.convertible_tga_ox_5prime[i];
        dst.g_count_5prime[i]    += src.g_count_5prime[i];
        dst.t_from_g_5prime[i]  += src.t_from_g_5prime[i];
        dst.c_count_ox_5prime[i]+= src.c_count_ox_5prime[i];
        dst.a_from_c_5prime[i]  += src.a_from_c_5prime[i];
    }
    dst.baseline_g_to_t_count  += src.baseline_g_to_t_count;
    dst.baseline_g_total       += src.baseline_g_total;
    dst.baseline_c_to_a_count  += src.baseline_c_to_a_count;
    dst.baseline_c_ox_total    += src.baseline_c_ox_total;

    dst.convertible_gag_interior += src.convertible_gag_interior;
    dst.convertible_tca_interior    += src.convertible_tca_interior;
    dst.convertible_tcg_interior    += src.convertible_tcg_interior;
    dst.convertible_tac_interior    += src.convertible_tac_interior;
    dst.convertible_tgc_interior    += src.convertible_tgc_interior;
    dst.convertible_taa_ca_interior += src.convertible_taa_ca_interior;
    dst.convertible_tag_ca_interior += src.convertible_tag_ca_interior;
    dst.convertible_tga_ca_interior += src.convertible_tga_ca_interior;

    dst.convertible_tag_ox_interior += src.convertible_tag_ox_interior;
    dst.convertible_gaa_interior += src.convertible_gaa_interior;
    dst.convertible_taa_ox_interior += src.convertible_taa_ox_interior;
    dst.convertible_gga_interior += src.convertible_gga_interior;
    dst.convertible_tga_ox_interior += src.convertible_tga_ox_interior;


    for (int i = 0; i < SampleDamageProfile::OxoTwoMarkerBins::TOTAL; ++i) {
        auto merge_cell = [](auto& dc, const auto& sc) {
            dc.n_reads += sc.n_reads;
            dc.sum_nGT += sc.sum_nGT; dc.sum_T += sc.sum_T;
            dc.sum_nAC += sc.sum_nAC; dc.sum_A += sc.sum_A;
        };
        merge_cell(dst.oxo_two_marker.cells[i],    src.oxo_two_marker.cells[i]);
        merge_cell(dst.oxo_two_marker_ss.cells[i], src.oxo_two_marker_ss.cells[i]);
    }

    for (int i = 0; i < SampleDamageProfile::N_OX_BINS; ++i) {
        auto& d = dst.oxidation_like_bins[i];
        const auto& s = src.oxidation_like_bins[i];
        for (int j = 0; j < SampleDamageProfile::N_OX_DEAM_STRATA; ++j) {
            auto& ds = d.strata[j];
            const auto& ss = s.strata[j];
            ds.term_t5 += ss.term_t5; ds.term_tc5 += ss.term_tc5;
            ds.term_a3 += ss.term_a3; ds.term_ag3 += ss.term_ag3;
            ds.int_t += ss.int_t; ds.int_tc += ss.int_tc;
            ds.int_a += ss.int_a; ds.int_ag += ss.int_ag;
            ds.sig_t += ss.sig_t; ds.sig_tg += ss.sig_tg;
            ds.sig_a += ss.sig_a; ds.sig_ac += ss.sig_ac;
            ds.ctrl_a += ss.ctrl_a; ds.ctrl_at += ss.ctrl_at;
            ds.ctrl_c += ss.ctrl_c; ds.ctrl_cg += ss.ctrl_cg;
            ds.reads += ss.reads;
        }
    }

    for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
        auto& db = dst.gc_bins[bin];
        const auto& sb = src.gc_bins[bin];
        for (int p = 0; p < 15; ++p) {
            db.t_counts[p] += sb.t_counts[p];
            db.c_counts[p] += sb.c_counts[p];
            db.a_counts[p] += sb.a_counts[p];
            db.g_counts[p] += sb.g_counts[p];
            db.t_counts_3prime[p] += sb.t_counts_3prime[p];
            db.c_counts_3prime[p] += sb.c_counts_3prime[p];
            db.a_counts_3prime[p] += sb.a_counts_3prime[p];
            db.g_counts_3prime[p] += sb.g_counts_3prime[p];
            db.stop_counts[p] += sb.stop_counts[p];
            db.pre_counts[p] += sb.pre_counts[p];
        }
        db.t_interior += sb.t_interior;
        db.c_interior += sb.c_interior;
        db.a_interior += sb.a_interior;
        db.g_interior += sb.g_interior;
        db.stop_interior += sb.stop_interior;
        db.pre_interior += sb.pre_interior;
        db.n_reads += sb.n_reads;
    }

    for (int b = 0; b < SampleDamageProfile::N_LEN_FINE; ++b) {
        auto& dl = dst.len_bins[b];
        const auto& sl = src.len_bins[b];
        for (int p = 0; p < 15; ++p) {
            dl.t_counts[p] += sl.t_counts[p];
            dl.c_counts[p] += sl.c_counts[p];
            dl.a_counts[p] += sl.a_counts[p];
            dl.g_counts[p] += sl.g_counts[p];
            dl.t_counts_3prime[p] += sl.t_counts_3prime[p];
            dl.c_counts_3prime[p] += sl.c_counts_3prime[p];
            dl.a_counts_3prime[p] += sl.a_counts_3prime[p];
            dl.g_counts_3prime[p] += sl.g_counts_3prime[p];
        }
        dl.t_interior += sl.t_interior;
        dl.c_interior += sl.c_interior;
        dl.a_interior += sl.a_interior;
        dl.g_interior += sl.g_interior;
        dl.n_reads += sl.n_reads;
        dl.len_sum += sl.len_sum;
        dl.jstrat_ds.add(sl.jstrat_ds);
        dl.jstrat_ss.add(sl.jstrat_ss);
    }

    dst.n_reads += src.n_reads;
}

// Paired-end variant. R1 contributes 5'-end counters + interior baseline;
// R2 (complement-mapped) contributes 3'-end counters. Read 3' ends are
// ignored — for inserts shorter than read length, R2's 5' end may read
// through into R1's 5' adapter, contaminating per_pos_3prime. The caller
// (fqdup PE worker) skips short pairs before calling this.
//
// Coverage scope: per-pos C->T and G->A counters at both ends, tail-anchored
// background counters, codon-position counters, hexamer counts at 5', and
// interior baseline. Advanced features filled by single-end update
// (CpG-like ctx, upstream ctx, oxoG 16-ctx, trinuc spectrum, channel D
// transversions, GC bins, channel B stop codons, channel C oxidative codons,
// interior CT cluster) are NOT recomputed here — PE mode is intended for
// raw bilateral 5'/3' damage QA, not full library profiling. Use SE on
// merged reads when those signals are needed.
bool FrameSelector::update_sample_profile_paired(
    SampleDamageProfile& profile,
    std::string_view r1,
    std::string_view r2)
{
    if (r1.length() < 30 || r2.length() < 30) return false;

    // Short-insert detection via R1/R2 overlap. When the molecule (insert)
    // is shorter than the read length, R1 reads through the molecule into
    // adapter A1 and R2 into adapter A2; the per-position damage windows
    // and tail counters then mix molecule and adapter bases, producing an
    // "anti-damage" shape at the 3' end (R2 first 15 bases are largely
    // adapter, complement-mapped into top-strand frame as A-depletion at
    // the 3'-end window).
    //
    // For an insert of length M, R1[M-K..M-1] is the molecule's 3' tail
    // and should reverse-complement to R2[0..K-1] (which reads the
    // molecule's bottom strand from the same end). We scan plausible M
    // values; require K=15 bases of overlap with at most 3 mismatches
    // (allows for sequencing error and aDNA damage at the molecule 3'
    // end). When overlap is detected, the pair is short-insert by
    // definition and belongs to the merged-read SE workflow — skip it.
    // Native PE is intended for true long-insert pairs (insert > read
    // length) where R1 and R2 do not overlap and the per-position
    // windows are clean molecule bases.
    {
        auto rc_base = [](char c) -> char {
            switch (c) { case 'A': return 'T'; case 'T': return 'A';
                         case 'C': return 'G'; case 'G': return 'C'; }
            return 'N';
        };
        constexpr int CHECK_LEN = 15;
        constexpr int MAX_MISMATCH = 3;
        const int max_M = static_cast<int>(std::min(r1.length(), r2.length()));
        bool overlap_found = false;
        for (int M = CHECK_LEN; M <= max_M; ++M) {
            int mismatches = 0;
            for (int i = 0; i < CHECK_LEN; ++i) {
                const char r1b = fast_upper(r1[M - 1 - i]);
                const char r2b = fast_upper(r2[i]);
                if (r1b == 'N' || r2b == 'N' || rc_base(r2b) != r1b) {
                    if (++mismatches > MAX_MISMATCH) break;
                }
            }
            if (mismatches <= MAX_MISMATCH) {
                overlap_found = true;
                break;
            }
        }
        if (overlap_found) {
            profile.pe_short_insert_skipped++;
            return false;
        }
    }

    profile.input_mode = SampleDamageProfile::InputMode::PAIRED;

    const size_t l1 = r1.length();
    const size_t l2 = r2.length();

    // R1 → 5' end counters
    for (size_t i = 0; i < std::min(size_t(15), l1); ++i) {
        char b = fast_upper(r1[i]);
        if (b == 'T')      { profile.t_freq_5prime[i]++; profile.tc_total_5prime[i]++; }
        else if (b == 'C') { profile.c_freq_5prime[i]++; profile.tc_total_5prime[i]++; }
        if (b == 'A')      profile.a_freq_5prime[i]++;
        else if (b == 'G') profile.g_freq_5prime[i]++;
    }

    // R2 → 3' end counters (complement-mapped: R2 reads bottom strand from
    // the molecule 3' inward, so R2[i] = complement(top_strand_at_3prime[i]).
    // The damage signal we want is G->A on top strand 3' end, which appears
    // as C->T on R2. Map R2 base to its complement, then accumulate with the
    // same logic as the SE 3'-end loop.
    for (size_t i = 0; i < std::min(size_t(15), l2); ++i) {
        char b = fast_upper(r2[i]);
        // complement: R2[i] → top_strand_at_3prime[i]
        char top;
        switch (b) {
            case 'A': top = 'T'; break;
            case 'T': top = 'A'; break;
            case 'C': top = 'G'; break;
            case 'G': top = 'C'; break;
            default:  top = 'N'; break;
        }
        if (top == 'A')      { profile.a_freq_3prime[i]++; profile.ag_total_3prime[i]++; }
        else if (top == 'G') { profile.g_freq_3prime[i]++; profile.ag_total_3prime[i]++; }
        if (top == 'T')      profile.t_freq_3prime[i]++;
        else if (top == 'C') profile.c_freq_3prime[i]++;
    }

    // 5' tail from R1, 3' tail from R2 (complement-mapped)
    {
        const int lo = SampleDamageProfile::BG_TAIL_LO;
        const int hi = SampleDamageProfile::BG_TAIL_HI;
        for (int i = lo; i <= hi && static_cast<size_t>(i) < l1; ++i) {
            const int idx = i - lo;
            const char b = fast_upper(r1[i]);
            if (b == 'T')      { profile.tail_t_5prime[idx]++; profile.tail_tc_5prime[idx]++; }
            else if (b == 'C') profile.tail_tc_5prime[idx]++;
        }
        for (int i = lo; i <= hi && static_cast<size_t>(i) < l2; ++i) {
            const int idx = i - lo;
            const char b = fast_upper(r2[i]);
            char top;
            switch (b) {
                case 'A': top = 'T'; break;
                case 'T': top = 'A'; break;
                case 'C': top = 'G'; break;
                case 'G': top = 'C'; break;
                default:  top = 'N'; break;
            }
            if (top == 'A')      { profile.tail_a_3prime[idx]++; profile.tail_ag_3prime[idx]++; }
            else if (top == 'G') profile.tail_ag_3prime[idx]++;
        }
    }

    // Interior baseline from R1's middle third (R2's middle would mostly
    // overlap for short inserts; tracking only R1's middle avoids double-
    // counting and keeps the baseline consistent with SE behavior).
    {
        constexpr size_t INTERIOR_TERM_PAD = 15;
        size_t mid_start = l1 / 3;
        size_t mid_end   = 2 * l1 / 3;
        if (mid_start < INTERIOR_TERM_PAD) mid_start = INTERIOR_TERM_PAD;
        if (l1 > INTERIOR_TERM_PAD && mid_end + INTERIOR_TERM_PAD > l1)
            mid_end = l1 - INTERIOR_TERM_PAD;
        if (mid_start < mid_end) {
            for (size_t i = mid_start; i < mid_end; ++i) {
                char b = fast_upper(r1[i]);
                if (b == 'T') profile.baseline_t_freq++;
                else if (b == 'C') profile.baseline_c_freq++;
                else if (b == 'A') profile.baseline_a_freq++;
                else if (b == 'G') profile.baseline_g_freq++;
            }
        }
    }

    // Codon-position counters: 5' from R1, 3' from R2 (complement-mapped)
    for (size_t i = 0; i < std::min(size_t(15), l1); ++i) {
        char b = fast_upper(r1[i]);
        int cp = i % 3;
        if (b == 'T') profile.codon_pos_t_count_5prime[cp]++;
        else if (b == 'C') profile.codon_pos_c_count_5prime[cp]++;
    }
    for (size_t i = 0; i < std::min(size_t(15), l2); ++i) {
        char b = fast_upper(r2[i]);
        char top;
        switch (b) { case 'A': top='T'; break; case 'T': top='A'; break;
                     case 'C': top='G'; break; case 'G': top='C'; break;
                     default: top='N'; }
        // Codon position in R2: position i in R2 == position (l_mol-1-i) in
        // top strand. Without alignment we can't recover exact codon phase
        // on the top strand, so we use the natural R2 codon phase (which
        // matches the molecule's 3' frame for inserts of length 3k).
        int cp = i % 3;
        if (top == 'A') profile.codon_pos_a_count_3prime[cp]++;
        else if (top == 'G') profile.codon_pos_g_count_3prime[cp]++;
    }

    // 5' hexamer + interior hexamer (from R1)
    if (l1 >= 18) {
        char hex_5prime[7];
        bool valid_5prime = true;
        for (int i = 0; i < 6; ++i) {
            char b = fast_upper(r1[i]);
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') { valid_5prime = false; break; }
            hex_5prime[i] = b;
        }
        hex_5prime[6] = '\0';
        if (valid_5prime) {
            uint32_t code = encode_hexamer(hex_5prime);
            if (code < 4096) {
                profile.hexamer_count_5prime[code] += 1.0;
                profile.n_hexamers_5prime++;
            }
        }

        size_t interior_start = l1 / 2 - 3;
        char hex_interior[7];
        bool valid_interior = true;
        for (int i = 0; i < 6; ++i) {
            char b = fast_upper(r1[interior_start + i]);
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') { valid_interior = false; break; }
            hex_interior[i] = b;
        }
        hex_interior[6] = '\0';
        if (valid_interior) {
            uint32_t code = encode_hexamer(hex_interior);
            if (code < 4096) {
                profile.hexamer_count_interior[code] += 1.0;
                profile.n_hexamers_interior++;
            }
        }

        // 3' terminal hexamer (last 6 bases of read)
        char hex_3prime[7];
        bool valid_3prime = true;
        for (int i = 0; i < 6; ++i) {
            char b = fast_upper(r1[l1 - 6 + i]);
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') { valid_3prime = false; break; }
            hex_3prime[i] = b;
        }
        hex_3prime[6] = '\0';
        if (valid_3prime) {
            uint32_t code = encode_hexamer(hex_3prime);
            if (code < 4096) {
                profile.hexamer_count_3prime[code] += 1.0;
                profile.n_hexamers_3prime++;
            }
        }
    }

    profile.n_reads++;
    return true;
}

void FrameSelector::update_sample_profile_weighted(
    SampleDamageProfile& profile,
    std::string_view seq,
    float weight) {

    if (seq.length() < 30 || weight < 0.001f) return;

    size_t len = seq.length();

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        char base = fast_upper(seq[i]);
        if (base == 'T') {
            profile.t_freq_5prime[i] += weight;
            profile.tc_total_5prime[i] += weight;
        } else if (base == 'C') {
            profile.c_freq_5prime[i] += weight;
            profile.tc_total_5prime[i] += weight;
        }
    }

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        if (base == 'A') {
            profile.a_freq_3prime[i] += weight;
            profile.ag_total_3prime[i] += weight;
        } else if (base == 'G') {
            profile.g_freq_3prime[i] += weight;
            profile.ag_total_3prime[i] += weight;
        }
    }

    // Tail-anchored background sampling (weighted variant)
    {
        const int lo = SampleDamageProfile::BG_TAIL_LO;
        const int hi = SampleDamageProfile::BG_TAIL_HI;
        for (int i = lo; i <= hi && static_cast<size_t>(i) < len; ++i) {
            const int idx = i - lo;
            const char b5 = fast_upper(seq[i]);
            if (b5 == 'T') { profile.tail_t_5prime[idx] += weight; profile.tail_tc_5prime[idx] += weight; }
            else if (b5 == 'C') { profile.tail_tc_5prime[idx] += weight; }

            const size_t pos3 = len - 1 - i;
            const char b3 = fast_upper(seq[pos3]);
            if (b3 == 'A') { profile.tail_a_3prime[idx] += weight; profile.tail_ag_3prime[idx] += weight; }
            else if (b3 == 'G') { profile.tail_ag_3prime[idx] += weight; }
        }
    }

    constexpr size_t INTERIOR_TERM_PAD = 15;
    size_t mid_start = len / 3;
    size_t mid_end   = 2 * len / 3;
    if (mid_start < INTERIOR_TERM_PAD)     mid_start = INTERIOR_TERM_PAD;
    if (len > INTERIOR_TERM_PAD && mid_end + INTERIOR_TERM_PAD > len)
        mid_end = len - INTERIOR_TERM_PAD;
    const bool interior_safe = (mid_start < mid_end);
    if (interior_safe) {
        for (size_t i = mid_start; i < mid_end; ++i) {
            char base = fast_upper(seq[i]);
            if (base == 'T') profile.baseline_t_freq += weight;
            else if (base == 'C') profile.baseline_c_freq += weight;
            else if (base == 'A') profile.baseline_a_freq += weight;
            else if (base == 'G') profile.baseline_g_freq += weight;
        }
    }

    size_t weight_count = std::max(size_t(1), static_cast<size_t>(weight * 10 + 0.5));
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        int codon_pos = i % 3;
        char base = fast_upper(seq[i]);
        if (base == 'T') profile.codon_pos_t_count_5prime[codon_pos] += weight_count;
        else if (base == 'C') profile.codon_pos_c_count_5prime[codon_pos] += weight_count;
    }

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        int codon_pos = (len - 1 - i) % 3;
        char base = fast_upper(seq[pos]);
        if (base == 'A') profile.codon_pos_a_count_3prime[codon_pos] += weight_count;
        else if (base == 'G') profile.codon_pos_g_count_3prime[codon_pos] += weight_count;
    }

    for (size_t i = 0; i < std::min(size_t(5), len - 1); ++i) {
        char base = fast_upper(seq[i]);
        char next = fast_upper(seq[i + 1]);

        if (next == 'G') {
            if (base == 'C') {
                profile.cpg_c_count += weight_count;
            } else if (base == 'T') {
                profile.cpg_t_count += weight_count;
            }
        } else {
            if (base == 'C') {
                profile.non_cpg_c_count += weight_count;
            } else if (base == 'T') {
                profile.non_cpg_t_count += weight_count;
            }
        }
    }

    profile.n_reads++;
}

void FrameSelector::reset_sample_profile(SampleDamageProfile& profile) {
    // Default-construct first: covers every field (including ones the previous
    // hand-rolled reset missed: d_max_*, damage_status, composition_bias_*,
    // inverted_pattern_*, library_bic_*, etc.) and stays in sync as new
    // members are added to SampleDamageProfile. Then restore the few
    // historical non-zero defaults the manual reset relied on (below).
    profile = SampleDamageProfile{};
    for (int i = 0; i < 15; ++i) {
        profile.t_freq_5prime[i] = 0.0;
        profile.c_freq_5prime[i] = 0.0;
        profile.a_freq_5prime[i] = 0.0;
        profile.g_freq_5prime[i] = 0.0;
        profile.a_freq_3prime[i] = 0.0;
        profile.g_freq_3prime[i] = 0.0;
        profile.t_freq_3prime[i] = 0.0;
        profile.c_freq_3prime[i] = 0.0;
        profile.damage_rate_5prime[i] = 0.0f;
        profile.damage_rate_3prime[i] = 0.0f;
        profile.tc_total_5prime[i] = 0.0;
        profile.ag_total_3prime[i] = 0.0;
    }

    profile.baseline_t_freq = 0.0;
    profile.baseline_c_freq = 0.0;
    profile.baseline_a_freq = 0.0;
    profile.baseline_g_freq = 0.0;

    for (int p = 0; p < 3; ++p) {
        profile.codon_pos_t_count_5prime[p] = 0;
        profile.codon_pos_c_count_5prime[p] = 0;
        profile.codon_pos_a_count_3prime[p] = 0;
        profile.codon_pos_g_count_3prime[p] = 0;
        profile.codon_pos_t_fraction_5prime[p] = 0.5f;
        profile.codon_pos_a_fraction_3prime[p] = 0.5f;
    }

    profile.cpg_c_count = 0;
    profile.cpg_t_count = 0;
    profile.non_cpg_c_count = 0;
    profile.non_cpg_t_count = 0;
    profile.cpg_ct_fraction = 0.0f;
    profile.non_cpg_ct_fraction = 0.0f;

    for (int ctx = 0; ctx < SampleDamageProfile::N_CT_CTX; ++ctx) {
        profile.ct_ctx_t_5prime[ctx].fill(0.0f);
        profile.ct_ctx_total_5prime[ctx].fill(0.0f);
        profile.ct_ctx_t_interior[ctx] = 0.0f;
        profile.ct_ctx_total_interior[ctx] = 0.0f;
    }
    profile.fit_baseline_ct5_cpg_like    = std::numeric_limits<float>::quiet_NaN();
    profile.fit_baseline_ct5_noncpg_like = std::numeric_limits<float>::quiet_NaN();
    profile.dmax_ct5_cpg_like    = std::numeric_limits<float>::quiet_NaN();
    profile.dmax_ct5_noncpg_like = std::numeric_limits<float>::quiet_NaN();
    profile.cpg_ratio     = std::numeric_limits<float>::quiet_NaN();
    profile.log2_cpg_ratio = std::numeric_limits<float>::quiet_NaN();

    // Reset upstream-context-aware accumulators
    for (int uctx = 0; uctx < SampleDamageProfile::N_UPSTREAM_CTX; ++uctx) {
        profile.ct5_t_by_upstream[uctx].fill(0.0);
        profile.ct5_total_by_upstream[uctx].fill(0.0);
        profile.ct5_t_interior_by_upstream[uctx] = 0.0;
        profile.ct5_total_interior_by_upstream[uctx] = 0.0;
        profile.dmax_ct5_by_upstream[uctx] = std::numeric_limits<float>::quiet_NaN();
        profile.baseline_ct5_by_upstream[uctx] = std::numeric_limits<float>::quiet_NaN();
        profile.cov_ct5_terminal_by_upstream[uctx] = 0.0f;
        profile.cov_ct5_interior_by_upstream[uctx] = 0.0f;
    }
    profile.dipyr_contrast = std::numeric_limits<float>::quiet_NaN();
    profile.context_heterogeneity_chi2 = 0.0f;
    profile.context_heterogeneity_chi2_raw = 0.0f;
    profile.context_heterogeneity_p = 1.0f;
    profile.context_heterogeneity_detected = false;
    profile.context_heterogeneity_computed = false;
    profile.effcov_ct5_cpg_like_terminal    = 0.0f;
    profile.effcov_ct5_noncpg_like_terminal = 0.0f;
    profile.effcov_ct5_cpg_like_interior    = 0.0f;
    profile.effcov_ct5_noncpg_like_interior = 0.0f;
    profile.cov_ct5_cpg_like_terminal       = 0.0f;
    profile.cov_ct5_noncpg_like_terminal    = 0.0f;
    profile.cov_ct5_cpg_like_interior       = 0.0f;
    profile.cov_ct5_noncpg_like_interior    = 0.0f;
    profile.fit_positions_ct5_cpg_like    = 0;
    profile.fit_positions_ct5_noncpg_like = 0;
    profile.oxog16_t.fill(0.0f);
    profile.oxog16_a_rc.fill(0.0f);
    profile.s_oxog_16ctx.fill(0.0f);
    profile.cov_oxog_16ctx.fill(0.0f);
    profile.tri_5prime_terminal.fill(0.0);
    profile.tri_5prime_interior.fill(0.0);
    profile.tri_3prime_terminal.fill(0.0);
    profile.tri_3prime_interior.fill(0.0);
    for (auto& a : profile.tri_5prime_terminal_by_deam) a.fill(0);
    for (auto& a : profile.tri_5prime_interior_by_deam) a.fill(0);
    for (auto& a : profile.tri_3prime_terminal_by_deam) a.fill(0);
    for (auto& a : profile.tri_3prime_interior_by_deam) a.fill(0);
    profile.deam_stratum_reads.fill(0);
    for (auto& a : profile.tri_5prime_pos) a.fill(0);
    for (auto& a : profile.tri_3prime_pos) a.fill(0);

    profile.max_damage_5prime = 0.0f;
    profile.max_damage_3prime = 0.0f;
    profile.sample_damage_prob = 0.0f;
    profile.lambda_5prime = 0.3f;
    profile.lambda_3prime = 0.3f;
    profile.lambda_5prime_fitted = false;  // D22
    profile.lambda_3prime_fitted = false;  // D22
    profile.terminal_shift_5prime = 0.0f;
    profile.terminal_shift_3prime = 0.0f;
    profile.terminal_z_5prime = 0.0f;
    profile.terminal_z_3prime = 0.0f;
    profile.terminal_inversion = false;
    profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
    profile.library_type_auto_detected = false;
    profile.tc_total_3prime.fill(0.0);

    profile.hexamer_count_5prime.fill(0.0);
    profile.hexamer_count_interior.fill(0.0);
    profile.hexamer_count_3prime.fill(0.0);
    profile.n_hexamers_5prime = 0;
    profile.n_hexamers_interior = 0;
    profile.n_hexamers_3prime = 0;
    profile.hexamer_damage_llr = 0.0f;

    profile.convertible_caa_5prime.fill(0.0);
    profile.convertible_taa_5prime.fill(0.0);
    profile.convertible_cag_5prime.fill(0.0);
    profile.convertible_tag_5prime.fill(0.0);
    profile.convertible_cga_5prime.fill(0.0);
    profile.convertible_tga_5prime.fill(0.0);
    profile.total_codons_5prime.fill(0.0);
    profile.convertible_caa_interior = 0.0;
    profile.convertible_taa_interior = 0.0;
    profile.convertible_cag_interior = 0.0;
    profile.convertible_tag_interior = 0.0;
    profile.convertible_cga_interior = 0.0;
    profile.convertible_tga_interior = 0.0;
    profile.total_codons_interior = 0.0;
    profile.stop_conversion_rate_baseline = 0.0f;
    profile.stop_decay_llr_5prime = 0.0f;
    profile.stop_amplitude_5prime = 0.0f;
    profile.channel_b_valid = false;
    profile.channel_b_valid_tga = false;
    profile.channel_b_valid_taa = false;  // D24
    profile.channel_b_valid_tag = false;  // D24
    profile.damage_validated = false;
    profile.damage_artifact = false;

    profile.joint_delta_max = 0.0f;
    profile.joint_lambda = 0.0f;
    profile.joint_a_max = 0.0f;
    profile.joint_log_lik_m1 = 0.0f;
    profile.joint_log_lik_m0 = 0.0f;
    profile.joint_delta_bic = 0.0f;
    profile.joint_bayes_factor = 0.0f;
    profile.joint_p_damage = 0.0f;
    profile.joint_z_delta = 0.0f;
    profile.joint_z_delta_capped = false;
    profile.joint_n_informative = 0;
    profile.joint_model_valid = false;

    profile.bulk_damage = {};
    profile.bulk_headline_delta = 0.0;
    profile.bulk_attempted = false;

    profile.mixture_K = 0;
    profile.mixture_n_components = 0;
    profile.mixture_d_population = 0.0f;
    profile.mixture_d_ancient = 0.0f;
    profile.mixture_d_population_highgc = 0.0f;
    profile.mixture_pi_ancient = 0.0f;
    profile.mixture_bic = 0.0f;
    profile.mixture_converged = false;
    profile.mixture_identifiable = false;
    profile.gc_histogram.fill(0.0);
    profile.adaptive_gc_threshold = 0.0f;
    profile.gc_threshold_computed = false;
    profile.gc_bins = {};
    profile.len_bins = {};
    profile.gc_stratified_d_max_weighted = 0.0f;
    profile.gc_stratified_d_max_joint = 0.0f;
    profile.gc_stratified_d_max_peak = 0.0f;
    profile.gc_peak_bin = -1;
    profile.gc_stratified_valid = false;
    profile.pi_damaged = 0.0f;
    profile.d_ancient = 0.0f;
    profile.d_population = 0.0f;
    profile.n_damaged_bins = 0;
    profile.n_reads_gc_filtered = 0;
    profile.n_reads_sampled = 0;

    profile.fit_offset_5prime = 1;
    profile.fit_offset_3prime = 1;

    profile.forced_library_type = SampleDamageProfile::LibraryType::UNKNOWN;
    profile.library_type_rescued = false;

    // oxoG / Channel-D accumulators
    profile.g_count_5prime.fill(0.0);
    profile.t_from_g_5prime.fill(0.0);
    profile.baseline_g_to_t_count = 0.0;
    profile.baseline_g_total = 0.0;
    profile.baseline_c_to_a_count = 0.0;
    profile.baseline_c_ox_total = 0.0;
    profile.convertible_gag_5prime.fill(0.0);
    profile.convertible_tca_5prime.fill(0.0);
    profile.convertible_tcg_5prime.fill(0.0);
    profile.convertible_tac_5prime.fill(0.0);
    profile.convertible_tgc_5prime.fill(0.0);
    profile.convertible_taa_ca_5prime.fill(0.0);
    profile.convertible_tag_ca_5prime.fill(0.0);
    profile.convertible_tga_ca_5prime.fill(0.0);
    profile.ca_deam_shadow_5prime.fill(0.0);
    profile.convertible_tca_3prime.fill(0.0);
    profile.convertible_tcg_3prime.fill(0.0);
    profile.convertible_tac_3prime.fill(0.0);
    profile.convertible_tgc_3prime.fill(0.0);
    profile.convertible_taa_ca_3prime.fill(0.0);
    profile.convertible_tag_ca_3prime.fill(0.0);
    profile.convertible_tga_ca_3prime.fill(0.0);
    profile.ca_deam_shadow_3prime.fill(0.0);
    profile.ca_stop_rate_baseline          = 0.0f;
    profile.ca_stop_rate_terminal          = 0.0f;
    profile.ca_stop_rate_interior          = 0.0f;
    profile.channel_f_z                    = std::numeric_limits<float>::quiet_NaN();  // C1: NaN = not computed
    profile.channel_f_mh_z                 = std::numeric_limits<float>::quiet_NaN();  // C1: NaN = not computed
    profile.channel_f_common_or            = 0.0f;
    profile.ca_uniformity_ratio            = 0.0f;
    profile.ca_stop_rate_baseline_3prime        = 0.0f;
    profile.ca_stop_rate_terminal_3prime   = 0.0f;
    profile.ca_stop_rate_interior_3prime   = 0.0f;
    profile.ca_uniformity_ratio_3prime     = 0.0f;
    profile.channel_f_valid                = false;
    profile.channel_f3_valid               = false;
    profile.convertible_tca_cg_5prime.fill(0.0);
    profile.convertible_tac_cg_5prime.fill(0.0);
    profile.convertible_tga_cg_5prime.fill(0.0);
    profile.convertible_tag_cg_5prime.fill(0.0);
    profile.convertible_tca_cg_3prime.fill(0.0);
    profile.convertible_tac_cg_3prime.fill(0.0);
    profile.convertible_tga_cg_3prime.fill(0.0);
    profile.convertible_tag_cg_3prime.fill(0.0);
    profile.cg_stop_rate_terminal          = 0.0f;
    profile.cg_stop_rate_interior          = 0.0f;
    profile.cg_stop_rate_baseline          = 0.0f;
    profile.channel_g_z                    = std::numeric_limits<float>::quiet_NaN();  // C1: NaN = not computed
    profile.channel_g_mh_z                 = std::numeric_limits<float>::quiet_NaN();  // Correction 3
    profile.channel_g_common_or            = 0.0f;
    profile.cg_uniformity_ratio            = 0.0f;
    profile.cg_stop_rate_terminal_3prime   = 0.0f;
    profile.cg_stop_rate_interior_3prime   = 0.0f;
    profile.cg_stop_rate_baseline_3prime   = 0.0f;
    profile.cg_uniformity_ratio_3prime     = 0.0f;
    profile.channel_g_valid                = false;
    profile.channel_g3_valid               = false;
    profile.convertible_aaa_h_5prime.fill(0.0);
    profile.convertible_aag_h_5prime.fill(0.0);
    profile.convertible_aga_h_5prime.fill(0.0);
    profile.convertible_taa_at_5prime.fill(0.0);
    profile.convertible_tag_at_5prime.fill(0.0);
    profile.convertible_tga_at_5prime.fill(0.0);
    profile.convertible_aaa_h_3prime.fill(0.0);
    profile.convertible_aag_h_3prime.fill(0.0);
    profile.convertible_aga_h_3prime.fill(0.0);
    profile.convertible_taa_at_3prime.fill(0.0);
    profile.convertible_tag_at_3prime.fill(0.0);
    profile.convertible_tga_at_3prime.fill(0.0);
    profile.at_stop_rate_terminal          = 0.0f;
    profile.at_stop_rate_interior          = 0.0f;
    profile.at_stop_rate_baseline          = 0.0f;
    profile.channel_h_z                    = std::numeric_limits<float>::quiet_NaN();  // C1: NaN = not computed
    profile.channel_h_z_p2plus             = std::numeric_limits<float>::quiet_NaN();  // C1: NaN = not computed
    profile.channel_h_mh_z                 = std::numeric_limits<float>::quiet_NaN();  // Correction 3
    profile.channel_h_common_or            = 0.0f;
    profile.channel_h_z_consistent         = false;
    profile.at_uniformity_ratio            = 0.0f;
    profile.at_stop_rate_terminal_3prime   = 0.0f;
    profile.at_stop_rate_interior_3prime   = 0.0f;
    profile.at_stop_rate_baseline_3prime   = 0.0f;
    profile.at_uniformity_ratio_3prime     = 0.0f;
    profile.channel_h_valid                = false;
    profile.channel_h3_valid               = false;
    profile.ca_pre_interior  = 0;
    profile.ca_stop_interior        = 0;
    profile.ca_deam_shadow_interior = 0;
    profile.ca_pre_interior_by_ctx.fill(0.0);
    profile.ca_stop_interior_by_ctx.fill(0.0);
    profile.ca_shadow_interior_by_ctx.fill(0.0);
    profile.ca_shadow_5prime_ctx0.fill(0.0);
    profile.ca_shadow_5prime_ctx1.fill(0.0);
    profile.ca_shadow_5prime_ctx2.fill(0.0);
    profile.cg_pre_interior  = 0;
    profile.cg_stop_interior = 0;
    profile.at_pre_interior  = 0;
    profile.at_stop_interior = 0;
    profile.cg_pre_interior_by_ctx.fill(0.0);   // Correction 3
    profile.cg_stop_interior_by_ctx.fill(0.0);
    profile.at_pre_interior_by_ctx.fill(0.0);
    profile.at_stop_interior_by_ctx.fill(0.0);


    profile.convertible_tca_interior = 0;
    profile.convertible_tcg_interior = 0;
    profile.convertible_tac_interior = 0;
    profile.convertible_tgc_interior = 0;
    profile.convertible_taa_ca_interior = 0;
    profile.convertible_tag_ca_interior = 0;
    profile.convertible_tga_ca_interior = 0;

    profile.convertible_gag_3prime.fill(0.0);
    profile.convertible_tag_ox_3prime.fill(0.0);
    profile.convertible_gaa_3prime.fill(0.0);
    profile.convertible_taa_ox_3prime.fill(0.0);
    profile.convertible_gga_3prime.fill(0.0);
    profile.convertible_tga_ox_3prime.fill(0.0);
    profile.ox_stop_rate_terminal_3prime = 0.0f;
    profile.ox_stop_rate_interior_3prime = 0.0f;
    profile.ox_stop_baseline_3prime      = 0.0f;
    profile.ox_uniformity_ratio_3prime   = 0.0f;
    profile.channel_c3_valid             = false;
    profile.ox_uniformity_ratio_3prime_computed = false;

    profile.convertible_gaa_5prime.fill(0.0);
    profile.convertible_gga_5prime.fill(0.0);
    profile.convertible_tag_ox_5prime.fill(0.0);
    profile.convertible_taa_ox_5prime.fill(0.0);
    profile.convertible_tga_ox_5prime.fill(0.0);
    profile.c_count_ox_5prime.fill(0.0);
    profile.a_from_c_5prime.fill(0.0);
    profile.oxidation_like_bins = {};
    profile.oxo_two_marker    = {};
    profile.oxo_two_marker_ss = {};
    profile.oxidation_like_signal = 0.0f;
    profile.oxidation_like_signal_se = 0.0f;
    profile.oxidation_like_control = 0.0f;
    profile.oxidation_like_control_se = 0.0f;
    profile.oxidation_like_adjusted = 0.0f;
    profile.oxidation_like_excess = 0.0f;
    profile.oxidation_like_se = 0.0f;
    profile.oxidation_like_z = 0.0f;
    profile.oxidation_like_reliability = 0.0f;
    profile.oxidation_like_bins_used = 0;
    profile.oxidation_like_effective_bins = 0.0f;
    profile.oxidation_like_heterogeneity = 0.0f;
    profile.oxidation_like_artifact_suspect = false;

    // Interior oxoG codon accumulators (merged in merge_sample_profiles)
    profile.convertible_gag_interior = 0;
    profile.convertible_gaa_interior = 0;
    profile.convertible_gga_interior = 0;
    profile.convertible_tag_ox_interior = 0;
    profile.convertible_taa_ox_interior = 0;
    profile.convertible_tga_ox_interior = 0;

    profile.ox_is_artifact = false;
    profile.ox_d_max = 0.0f;
    profile.ox_damage_detected = false;
    profile.ox_damage_detected_codon = false;
    profile.ox_damage_detected_model = false;

    // Neutral default for oxidation uniformity ratio (real-zero vs not-computed is
    // disambiguated by ox_uniformity_ratio_computed, reset to false here).
    profile.ox_uniformity_ratio = 1.0f;
    profile.ox_uniformity_ratio_computed     = false;
    profile.ox_stop_rate_positional_computed = false;
    profile.ox_gt_uniformity_computed        = false;
    profile.d_computed                       = false;
    profile.channel_e_valid                  = false;

    // GT exponential-background fit boundary/degeneracy flags (C4) + companions.
    profile.gt_decay_at_upper_boundary = false;
    profile.gt_term_zero_clamped       = false;
    profile.gt_bg_at_upper_boundary    = false;
    profile.g_bg_fitted_unclamped      = 0.0f;
    profile.ox_theta_at_clamp          = false;
    // binom_z fields (channel_f_z/f_mh_z/g_z/h_z/h_z_p2plus) are reset to NaN
    // in the per-channel F/G/H reset blocks above.

    profile.n_reads = 0;
}

SampleDamageProfile FrameSelector::compute_sample_profile(
    const std::vector<std::string>& sequences) {

    SampleDamageProfile profile;
    reset_sample_profile(profile);

    for (const auto& seq : sequences) {
        update_sample_profile(profile, seq);
    }

    finalize_sample_profile(profile);
    return profile;
}

} // namespace taph
