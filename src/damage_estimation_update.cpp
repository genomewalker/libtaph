// FrameSelector: accumulation, merge, and reset.
#include "taph/frame_selector_decl.hpp"
#include "taph/codon_tables.hpp"
#include "taph/hexamer_tables.hpp"
#include <algorithm>
#include <cmath>
#include <array>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>
namespace taph {

// Cross-thread-deterministic fixed-point: returns round_half_to_even(num * SCALE / den)
// as __int128, with num,den >= 0 and den > 0. The whole computation is exact integer
// arithmetic, so accumulating these contributions in any order yields a bit-identical
// sum (unlike the non-associative double products they replace). Half-to-even (not
// half-away-from-zero / llround) avoids a one-directional bias that would otherwise
// accumulate over millions of positive addends.
//
// Perf note (2026-07-02): profiled at ~66.5% of update_sample_profile's CPU time on a
// real dataset, with __int128 division/multiplication as the dominant cost within that.
// CLUSTER_FP_SCALE was lowered from 2^40 to 2^20 (see sample_damage_profile.hpp) because
// the result is only ever consumed as a double at finalize (52-bit mantissa) -- the extra
// 20 bits of fixed-point precision were being computed and then discarded. At the smaller
// scale, num*SCALE fits in int64_t for every call site actually exercised (verified: the
// widest is P*A*(D-A) with P<=wlen<=~150, A,D<~22500 for realistic window/read lengths,
// giving num < 2^30, scaled < 2^50 -- comfortably inside int64_t's ~2^63 range). The fast
// path below still verifies this at runtime rather than assuming it, so correctness for
// any future caller with larger inputs is never at risk -- it just falls back to the
// exact __int128 path (identical arithmetic, same result) if the runtime check fails.
// Determinism-verified: a standalone accumulation-order test (matching this exact
// algorithm, both scales) confirmed order-independent bit-identical sums, and the
// double-cast output differs by ~6e-11 relative between scale=2^40 and scale=2^20 --
// far below the double/log2 consumption noise floor this value ultimately feeds.
static inline __int128 fp_round_div(__int128 num, __int128 den) {
    constexpr __int128 kInt64SafeBound = static_cast<__int128>(1) << 62;
    if (num >= 0 && num < kInt64SafeBound && den > 0 && den < kInt64SafeBound) {
        const int64_t n64 = static_cast<int64_t>(num);
        const int64_t d64 = static_cast<int64_t>(den);
        // Re-check the actual scaled product fits int64_t before committing to the
        // fast path -- num alone being < 2^62 does not guarantee num*SCALE does.
        const int64_t scale64 = static_cast<int64_t>(SampleDamageProfile::CLUSTER_FP_SCALE);
        if (n64 <= (std::numeric_limits<int64_t>::max)() / scale64) {
            const int64_t scaled = n64 * scale64;
            int64_t q = scaled / d64;
            int64_t r = scaled - q * d64;
            const int64_t twice = r << 1;
            if (twice > d64) {
                ++q;
            } else if (twice == d64) {
                if (q & 1) ++q;
            }
            return static_cast<__int128>(q);
        }
    }
    // Slow, always-correct path (unchanged from the original implementation).
    const __int128 scaled = num * SampleDamageProfile::CLUSTER_FP_SCALE;
    __int128 q = scaled / den;
    __int128 r = scaled - q * den;            // 0 <= r < den (num,den >= 0)
    const __int128 twice = r << 1;
    if (twice > den) {
        ++q;
    } else if (twice == den) {
        if (q & 1) ++q;                        // round half to even
    }
    return q;
}

// F2: methylation context of a 5' C->T site from its +1 / +2 downstream bases on
// the SAME strand (strand-aware: the read is already oriented to the cytosine).
//   CpG  : +1 == G
//   CHG  : +1 in {A,C,T} and +2 == G
//   CHH  : +1 in {A,C,T} and +2 in {A,C,T}
// Returns -1 when +2 is ambiguous (N / off-end) so the site is dropped, never
// misassigned between CHG and CHH.
static inline int ct_meth_context(char plus1, char plus2) {
    if (plus1 == 'G') return SampleDamageProfile::CpG;
    if (plus1 != 'A' && plus1 != 'C' && plus1 != 'T') return -1;
    if (plus2 == 'G') return SampleDamageProfile::CHG;
    if (plus2 == 'A' || plus2 == 'C' || plus2 == 'T') return SampleDamageProfile::CHH;
    return -1;
}

void FrameSelector::update_sample_profile(
    SampleDamageProfile& profile,
    std::string_view seq) {

    if (seq.length() < 30) return;  // Too short for reliable statistics

    // Resolve the user-forced prep into library_type ONCE, before any accumulation. library_type is
    // the DOUBLE_STRANDED default until finalize_libtype runs; every stream-time consumer (the ss/ds
    // 3' co-occurrence base below, and any future one) must see the real prep, not the default. Without
    // this, a forced-ss library silently streamed as ds and counted the 3' dA-tail (#A) as damage
    // instead of the ss C->T signal (#T). No-op when unforced (UNKNOWN): library_type stays the default
    // and finalize_libtype's BIC auto-detect owns it. Idempotent; finalize_libtype overwrites it anyway.
    if (profile.forced_library_type != SampleDamageProfile::LibraryType::UNKNOWN)
        profile.library_type = profile.forced_library_type;

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

    // Per-position terminal counts for identifiable pi (sample_damage_profile.hpp / finalize_tau).
    // Runs for ALL reads >= 2*P_PI+2 — short, ancient-enriched reads (which interior_safe excludes) included.
    if (len >= 2 * SampleDamageProfile::P_PI + 2) {
        constexpr int P = SampleDamageProfile::P_PI;
        int k5 = 0, cnum = 0;                         // 5' damage-allele count + centroid numerator (for F_class)
        for (int p = 0; p < P; ++p) {
            const char b = decoded[p];
            if (b == 'T') { ++k5; cnum += p; }
        }
        // 3' damage allele is prep-dependent: ss deaminates C->T at BOTH ends (#T), ds shows G->A (#A).
        // UNKNOWN keeps the ds default 'A' (the rest of the pipeline also treats UNKNOWN as ds-ish).
        // library_type is authoritative here: the forced prep was resolved into it at the top of this
        // function (before any accumulation); unforced stays the ds default until finalize resolves it.
        const char dmg3 = (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) ? 'T' : 'A';
        int c3a = 0, c3t = 0;                         // both 3' alleles (QC cross-check, prep-independent)
        for (int p = 0; p < P; ++p) {
            const char b = decoded[len - 1 - p];
            if      (b == 'A') ++c3a;
            else if (b == 'T') ++c3t;
        }
        const int k3a = (dmg3 == 'T') ? c3t : c3a;    // 3' damage-allele count for the FORCED co-occurrence gate
        profile.cooc_qc_term_n += 1; profile.cooc_qc_term_a3 += c3a; profile.cooc_qc_term_t3 += c3t;
        // F_class = 5' C->T decay centroid (overhang proxy): 0=no event, 1=blunt/terminal, 2=broad overhang.
        int C = 0;
        if (k5 > 0)
            C = (static_cast<double>(cnum) / k5 < 0.5 * (P - 1)) ? 1 : 2;
        const int L = SampleDamageProfile::pi_len_bin(static_cast<int>(len));
        auto& a5 = profile.pi_pos_5prime[L][C];
        auto& a3d = profile.pi_pos_3prime_ds[L][C];
        auto& a3s = profile.pi_pos_3prime_ss[L][C];
        for (int p = 0; p < P; ++p) {
            const char b5 = decoded[p];
            if      (b5 == 'T') { ++a5[p].n_elig; ++a5[p].n_deam; }
            else if (b5 == 'C') { ++a5[p].n_elig; }
            const char b3 = decoded[len - 1 - p];
            if      (b3 == 'A') { ++a3d[p].n_elig; ++a3d[p].n_deam; }
            else if (b3 == 'G') { ++a3d[p].n_elig; }
            if      (b3 == 'T') { ++a3s[p].n_elig; ++a3s[p].n_deam; }
            else if (b3 == 'C') { ++a3s[p].n_elig; }
        }
        auto& cc = profile.pi_cooc[L][C];
        cc.n += 1;
        cc.sum_k5   += static_cast<uint64_t>(k5);
        cc.sum_k3   += static_cast<uint64_t>(k3a);
        cc.sum_k5k3 += static_cast<uint64_t>(k5) * k3a;
        cc.hist[k5 * (P + 1) + k3a] += 1;   // joint (k5,k3) dist for Var(Cov)

        // Composition NULL via WITHIN-READ SHUFFLE (route (b) of the fix design). We recount k5(#T in
        // the first P) and k3(#A in the last P) after a DETERMINISTIC permutation of THIS read's bases.
        // The permutation preserves the read's exact base composition (so the between-read compositional
        // covariance that also drives the terminal covT is reproduced identically) while destroying the
        // terminal damage localization, so covI is a clean composition-only pedestal measured in the SAME
        // window geometry as covT. The fixed-budget hypergeometric anti-correlation therefore CANCELS in
        // covT−covI. Interior sub-windows (the previous null) could NOT be used: on a short read there is
        // no room to place two damage-free windows at the terminal end-separation, so their covariance
        // pedestal did not match covT and over-subtracted (covI swung to −0.07..−0.22 at short bins).
        // Seed is derived from read content → fully reproducible regardless of thread scheduling.
        if (len <= 512) {                                          // ceiling: >512 bp read skipped (rare, L=5 long bin)
            uint16_t idx[512];
            const int n = static_cast<int>(len);
            for (int i = 0; i < n; ++i) idx[i] = static_cast<uint16_t>(i);
            uint64_t rng = 0x9E3779B97F4A7C15ULL ^ (static_cast<uint64_t>(len) * 0x100000001B3ULL)
                           ^ (static_cast<uint64_t>(k5) << 32) ^ static_cast<uint64_t>(k3a)
                           ^ (static_cast<uint64_t>(static_cast<uint8_t>(decoded[0])) << 8)
                           ^ static_cast<uint64_t>(static_cast<uint8_t>(decoded[n - 1]));
            auto next = [&rng]() {                                 // splitmix64
                rng += 0x9E3779B97F4A7C15ULL;
                uint64_t z = rng;
                z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
                z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
                return z ^ (z >> 31);
            };
            for (int i = 0; i < 2 * P; ++i) {                      // partial Fisher-Yates: 2P disjoint draws
                const int j = i + static_cast<int>(next() % static_cast<uint64_t>(n - i));
                const uint16_t t = idx[i]; idx[i] = idx[j]; idx[j] = t;
            }
            int k5n = 0, k3nA = 0, k3nT = 0;
            for (int p = 0; p < P; ++p)       if (decoded[idx[p]] == 'T') ++k5n;
            for (int p = P; p < 2 * P; ++p) {
                const char b = decoded[idx[p]];
                if      (b == 'A') ++k3nA;
                else if (b == 'T') ++k3nT;
            }
            const int k3n = (dmg3 == 'T') ? k3nT : k3nA;   // null uses the SAME forced prep base
            profile.cooc_qc_null_n += 1; profile.cooc_qc_null_a3 += k3nA; profile.cooc_qc_null_t3 += k3nT;
            auto& ci = profile.pi_cooc_interior_ds[L];
            ci.n += 1;
            ci.sum_k5   += static_cast<uint64_t>(k5n);
            ci.sum_k3   += static_cast<uint64_t>(k3n);
            ci.sum_k5k3 += static_cast<uint64_t>(k5n) * k3n;
            ci.hist[k5n * (P + 1) + k3n] += 1;
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

    // (Composition NULL for the co-occurrence π now comes from the within-read shuffle above, matched to
    // the terminal read set and window geometry — see the pi_cooc_interior_ds accumulation in the
    // len>=2P+2 block. The old interior sub-window null was removed: it over-subtracted on short reads.)

    // Reference-free oxidation-like contrast. During the streaming pass we do
    // not assign damaged/background weights directly. Instead, each read is
    // placed into a terminal-deamination-excess stratum using its own interior
    // composition as a null. Finalization calibrates the strata within length
    // x GC bins and only then compares high-deamination to low-deamination
    // strata. This avoids treating raw terminal T-richness as an ancestry score.
    int tetra_deam_bin = -1; // cross-end proxy stratum for 5' tetra (ct_5prime_by_deam). See below.
    if (interior_safe) {
        double term_t5 = 0.0, term_tc5 = 0.0, term_a5 = 0.0, term_ag5 = 0.0;
        for (size_t p = 0; p < std::min<size_t>(5, len); ++p) {
            const char b = decoded[p];
            if      (b == 'T') { term_t5 += 1.0; term_tc5 += 1.0; }
            else if (b == 'C') { term_tc5 += 1.0; }
            else if (b == 'A') { term_a5 += 1.0; term_ag5 += 1.0; }
            else if (b == 'G') { term_ag5 += 1.0; }
        }

        // 3' position 0 is skipped: SS prep can create a ligation artifact at
        // the final base. Positions 1..4 carry the useful G->A deamination axis.
        // T/C counts at 3' are also accumulated for the strand-discordant CT3 excess.
        double term_a3 = 0.0, term_ag3 = 0.0, term_t3 = 0.0, term_tc3 = 0.0;
        for (size_t off = 1; off < std::min<size_t>(5, len); ++off) {
            const char b = decoded[len - 1 - off];
            if (b == 'A') { term_a3 += 1.0; term_ag3 += 1.0; }
            else if (b == 'G') { term_ag3 += 1.0; }
            else if (b == 'T') { term_t3 += 1.0; term_tc3 += 1.0; }
            else if (b == 'C') { term_tc3 += 1.0; }
        }

        double mid_t = 0.0, mid_c = 0.0, mid_a = 0.0, mid_g = 0.0;
        for (size_t i = mid_start; i < mid_end; ++i) {
            const char b = decoded[i];
            if (b == 'T') ++mid_t;
            else if (b == 'C') ++mid_c;
            else if (b == 'A') ++mid_a;
            else if (b == 'G') ++mid_g;
        }
        // Interior composition baseline (cross-read accumulator). Folded from a
        // separate interior pass: mid_* are the exact per-read integer counts, so
        // one add == the prior per-base increments, bit-for-bit.
        profile.baseline_t_freq += mid_t;
        profile.baseline_c_freq += mid_c;
        profile.baseline_a_freq += mid_a;
        profile.baseline_g_freq += mid_g;

        const double mid_tc = mid_t + mid_c;
        const double mid_ag = mid_a + mid_g;
        const double mid_total = mid_t + mid_c + mid_a + mid_g;
        if (mid_total > 0.0) {
            double score_num = 0.0, score_den = 0.0;
            double ct5_exc = 0.0, ga3_exc = 0.0;
            double ga5_exc = 0.0;
            // ct3_exc: 3' C→T excess vs interior — strand-discordant control.
            // Computed before deam_score so it is available for cross-end stratification below.
            // Genuine ds-aDNA damage links CT5↔GA3 but not CT5↔CT3; artifacts may link both.
            double ct3_exc = 0.0;
            if (term_tc3 > 0.0 && mid_tc > 0.0)
                ct3_exc = std::max(0.0, (term_t3 + 0.5) / (term_tc3 + 1.0)
                                      - (mid_t  + 0.5) / (mid_tc  + 1.0));
            if (term_tc5 > 0.0 && mid_tc > 0.0) {
                const double term = (term_t5 + 0.5) / (term_tc5 + 1.0);
                const double base = (mid_t + 0.5) / (mid_tc + 1.0);
                ct5_exc    = std::max(0.0, term - base);
                score_num += ct5_exc * term_tc5;
                score_den += term_tc5;
            }
            if (term_ag5 > 0.0 && mid_ag > 0.0) {
                const double term = (term_a5 + 0.5) / (term_ag5 + 1.0);
                const double base = (mid_a + 0.5) / (mid_ag + 1.0);
                ga5_exc = std::max(0.0, term - base);
            }
            if (term_ag3 > 0.0 && mid_ag > 0.0) {
                const double term = (term_a3 + 0.5) / (term_ag3 + 1.0);
                const double base = (mid_a + 0.5) / (mid_ag + 1.0);
                ga3_exc    = std::max(0.0, term - base);
                score_num += ga3_exc * term_ag3;
                score_den += term_ag3;
            }

            // Cross-end stratification for 5' tetra spectrum (ct_5prime_by_deam).
            // Uses max(ga3_exc, ct3_exc): DS reads score via 3' G→A; SS reads via 3' C→T.
            // Axis is orthogonal to the 5' C→T values being reported, eliminating the
            // circular selection bias and fixing SS mis-stratification in one change.
            // Shared terminal-deamination-excess stratum edges (0=none .. 4=heaviest).
            // Used for BOTH the cross-end tetra stratum and the 5' deam_score stratum.
            constexpr double DEAM_BIN_T1 = 0.08, DEAM_BIN_T2 = 0.20, DEAM_BIN_T3 = 0.40;
            auto deam_stratum = [](double excess) -> int {
                if      (excess > DEAM_BIN_T3) return 4;
                else if (excess > DEAM_BIN_T2) return 3;
                else if (excess > DEAM_BIN_T1) return 2;
                else if (excess > 0.00)        return 1;
                return 0;
            };

            tetra_deam_bin = deam_stratum(std::max(ga3_exc, ct3_exc));

            // ---- Cross-fit de-circularized strata ----
            // Split the 3' terminus into two position-folds; key the stratum from one fold,
            // read the misincorporation rate from the held-out fold (both directions pooled).
            // Decoupled interior baselines per fold stop the GC-sort lockstep from re-coupling
            // key and readout. _null reads the same fold under a damage-independent uniform key.
            if (len >= 6) {
                uint64_t h = 1469598103934665603ull;          // FNV-1a over the read
                for (size_t i = 0; i < len; ++i) { h ^= static_cast<uint8_t>(decoded[i]); h *= 1099511628211ull; }
                auto fold = [&](uint64_t x) -> int { return static_cast<int>((h ^ x) & 1ull); };

                // per-fold 3' terminal sums: ga (A success / A+G trial), ct (T / T+C)
                long ga_a[2]={0,0}, ga_n[2]={0,0}, ct_t[2]={0,0}, ct_n[2]={0,0};
                for (size_t off = 1; off < std::min<size_t>(5, len); ++off) {
                    const char b = decoded[len - 1 - off]; const int f = fold(off);
                    if      (b == 'A') { ++ga_a[f]; ++ga_n[f]; }
                    else if (b == 'G') { ++ga_n[f]; }
                    if      (b == 'T') { ++ct_t[f]; ++ct_n[f]; }
                    else if (b == 'C') { ++ct_n[f]; }
                }
                // per-fold interior sums (decoupled baseline) + GC for composition gate
                long iga_a[2]={0,0}, iga_n[2]={0,0}, ict_t[2]={0,0}, ict_n[2]={0,0};
                long igc=0, itot=0;
                for (size_t i = mid_start; i < mid_end; ++i) {
                    const char b = decoded[i]; const int f = fold(i); ++itot;
                    if      (b == 'A') { ++iga_a[f]; ++iga_n[f]; }
                    else if (b == 'G') { ++iga_n[f]; ++igc; }
                    else if (b == 'T') { ++ict_t[f]; ++ict_n[f]; }
                    else if (b == 'C') { ++ict_n[f]; ++igc; }
                }
                auto exc = [](long k, long n, long ik, long in_) -> double {
                    if (n <= 0 || in_ <= 0) return 0.0;
                    return std::max(0.0, (k + 0.5) / (n + 1.0) - (ik + 0.5) / (in_ + 1.0));
                };
                int strat_f[2];
                for (int f = 0; f < 2; ++f)
                    strat_f[f] = deam_stratum(std::max(
                        exc(ga_a[f], ga_n[f], iga_a[f], iga_n[f]),
                        exc(ct_t[f], ct_n[f], ict_t[f], ict_n[f])));
                // direction d keys on fold d, reads out the held-out fold r=1-d
                for (int d = 0; d < 2; ++d) {
                    const int r = 1 - d, S = strat_f[d];
                    profile.cf_ga3.term_k[S] += ga_a[r];  profile.cf_ga3.term_n[S] += ga_n[r];
                    profile.cf_ga3.intr_k[S] += iga_a[r]; profile.cf_ga3.intr_n[S] += iga_n[r];
                    profile.cf_ct3.term_k[S] += ct_t[r];  profile.cf_ct3.term_n[S] += ct_n[r];
                    profile.cf_ct3.intr_k[S] += ict_t[r]; profile.cf_ct3.intr_n[S] += ict_n[r];
                }
                // null: damage-independent uniform key (high bits of h), full readout mass
                const int Sn = static_cast<int>((h >> 17) % SampleDamageProfile::N_OX_DEAM_STRATA);
                profile.cf_ga3.null_term_k[Sn] += ga_a[0]+ga_a[1];  profile.cf_ga3.null_term_n[Sn] += ga_n[0]+ga_n[1];
                profile.cf_ga3.null_intr_k[Sn] += iga_a[0]+iga_a[1]; profile.cf_ga3.null_intr_n[Sn] += iga_n[0]+iga_n[1];
                profile.cf_ct3.null_term_k[Sn] += ct_t[0]+ct_t[1];  profile.cf_ct3.null_term_n[Sn] += ct_n[0]+ct_n[1];
                profile.cf_ct3.null_intr_k[Sn] += ict_t[0]+ict_t[1]; profile.cf_ct3.null_intr_n[Sn] += ict_n[0]+ict_n[1];
                // per-stratum descriptive stats keyed on the full-read stratum
                if (tetra_deam_bin >= 0) {
                    profile.cf_reads[tetra_deam_bin]     += 1;
                    profile.cf_len_sum[tetra_deam_bin]   += len;
                    profile.cf_len_sumsq[tetra_deam_bin] += static_cast<uint64_t>(len) * len;
                    profile.cf_igc_num[tetra_deam_bin]   += igc;
                    profile.cf_igc_den[tetra_deam_bin]   += itot;
                }
            }

            const double deam_score = score_den > 0.0 ? score_num / score_den : 0.0;
            int deam_bin = deam_stratum(deam_score);
            if (deam_score > 0.0) {
                profile.per_read_deam_sum   += deam_score;
                profile.per_read_deam_sumsq += deam_score * deam_score;
                ++profile.per_read_deam_n;
            }
            // Cross-cumulant sufficient stats (all reads, not just score>0).
            // tpg5: strictly TpG at pos0-1 (decoded[0]='T', decoded[1]='G').
            //   Stratified accumulation — κ₂_TpG / n_TpG vs κ₂ / n tests CpG deamination enrichment.
            //   Using only T (not C) avoids mixing original CpG (undamaged) with deaminated CpG.
            const bool tpg5 = (len >= 2) && (decoded[0] == 'T') && (decoded[1] == 'G');
            // Composition-immune channels: g5 = ct5_exc − ga5_exc, g3 = ga3_exc − ct3_exc.
            // Symmetric AT-terminus raises both T/(T+C) and A/(A+G) equally → cancels in diff.
            // Genuine 5' C→T deamination raises only ct5_exc → g5 > 0; composition → g5 ≈ 0.
            const double g5 = ct5_exc - ga5_exc;
            const double g3 = ga3_exc - ct3_exc;
            profile.per_read_ct5_sum    += ct5_exc;
            profile.per_read_ga3_sum    += ga3_exc;
            profile.per_read_ct5ga3     += ct5_exc * ga3_exc;
            profile.per_read_ct5ga3_cpg += tpg5 ? (ct5_exc * ga3_exc) : 0.0;
            profile.per_read_n_tpg      += tpg5 ? 1.0 : 0.0;
            profile.per_read_ct5ct3     += ct5_exc * ct3_exc;
            profile.per_read_g5_sum     += g5;
            profile.per_read_g3_sum     += g3;
            profile.per_read_g5g3       += g5 * g3;
            if (len > 0)
                profile.per_read_score_len += deam_score / static_cast<double>(len);

            const double gc_frac = (mid_g + mid_c) / mid_total;
            const int gc_bin = std::clamp(
                static_cast<int>(gc_frac * SampleDamageProfile::N_OX_GC_BINS),
                0, SampleDamageProfile::N_OX_GC_BINS - 1);
            // Oxidation co-movement: q-weighted per-GC-bin G->T / C->A damaged-vs-non-damaged
            // contrast. q = clamp(deam_score,0,1); middle-third counts are this read's, already
            // strand-oriented. compute_ox_scores derives s_oxog/s_ca/d_oriented at finalize.
            update_ox_bins(profile.ox_comov_bins[gc_bin],
                           std::clamp(deam_score, 0.0, 1.0),
                           static_cast<int>(mid_t), static_cast<int>(mid_g),
                           static_cast<int>(mid_c), static_cast<int>(mid_a));
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

            // sig/ctrl interior sums are pure functions of the mid_* composition
            // counts already scanned above; one add each == the prior per-base
            // increments (exact integer counts in double), bit-for-bit identical.
            oxs.sig_tg  += mid_t + mid_g; oxs.sig_t += mid_t;
            oxs.sig_ac  += mid_a + mid_c; oxs.sig_a += mid_a;
            oxs.ctrl_at += mid_a + mid_t; oxs.ctrl_a += mid_a;
            oxs.ctrl_cg += mid_c + mid_g; oxs.ctrl_c += mid_c;

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
        auto add_ctx_s = [&](int prev_pos, int mid_pos, int next_pos,
                             std::array<uint64_t, SampleDamageProfile::N_TRINUC>& bulk) {
            if (prev_pos < 0 || next_pos >= static_cast<int>(len)) return;
            int i0 = nuc_idx(decoded[prev_pos]);
            int i1 = nuc_idx(decoded[mid_pos]);
            int i2 = nuc_idx(decoded[next_pos]);
            if (i0 < 0 || i1 < 0 || i2 < 0) return;
            ++bulk[i0 * 16 + i1 * 4 + i2];
        };
        // tri gate p+1<len reaches d=1 (the 3′ CpG position tetra's p+2<len cannot).
        for (int p = 1; p <= 4 && p + 1 < static_cast<int>(len); ++p)
            add_ctx_s(p - 1, p, p + 1, profile.tri_5prime_terminal);
        for (int p = 10; p <= 14 && p + 1 < static_cast<int>(len); ++p)
            add_ctx_s(p - 1, p, p + 1, profile.tri_5prime_interior);
        for (int p = 1; p <= 4; ++p) {
            int mid = static_cast<int>(len) - 1 - p;
            add_ctx_s(mid - 1, mid, mid + 1, profile.tri_3prime_terminal);
        }
        for (int p = 10; p <= 14; ++p) {
            int mid = static_cast<int>(len) - 1 - p;
            add_ctx_s(mid - 1, mid, mid + 1, profile.tri_3prime_interior);
        }
        // Per-position: all positions 1..N_POS_TRI-1 from each end.
        for (int p = 1; p < SampleDamageProfile::N_POS_TRI && p + 1 < static_cast<int>(len); ++p)
            add_ctx(p - 1, p, p + 1, profile.tri_5prime_pos[p]);
        for (int p = 1; p < SampleDamageProfile::N_POS_TRI; ++p) {
            int mid = static_cast<int>(len) - 1 - p;
            if (mid - 1 >= 0 && mid + 1 < static_cast<int>(len))
                add_ctx(mid - 1, mid, mid + 1, profile.tri_3prime_pos[p]);
        }
        // Stratified per-position tri: same positions, keyed on tetra_deam_bin (cross-end proxy).
        // Strata 0=non-damaged .. N_OX_DEAM_STRATA-1=damaged. Gate mirrors tetra (bulk==Σstrata).
        if (tetra_deam_bin >= 0) {
            for (int p = 1; p < SampleDamageProfile::N_POS_TRI && p + 1 < static_cast<int>(len); ++p)
                add_ctx(p - 1, p, p + 1, profile.tri_5prime_pos_by_deam[tetra_deam_bin][p]);
            for (int p = 1; p < SampleDamageProfile::N_POS_TRI; ++p) {
                int mid = static_cast<int>(len) - 1 - p;
                if (mid - 1 >= 0 && mid + 1 < static_cast<int>(len))
                    add_ctx(mid - 1, mid, mid + 1, profile.tri_3prime_pos_by_deam[tetra_deam_bin][p]);
            }
        }

        // 4-mer (tetranucleotide) counts for CHG/CHH separation.
        // Requires prev(p-1), mid(p), next1(p+1), next2(p+2) — all in bounds.
        // Encoding: i0*64 + i1*16 + i2*4 + i3  (A=0,C=1,G=2,T=3).
        auto add_ctx4 = [&](int p,
                            std::array<uint64_t, SampleDamageProfile::N_TETRANUC>& bulk,
                            std::array<std::array<uint64_t, SampleDamageProfile::N_TETRANUC>,
                                       SampleDamageProfile::N_OX_DEAM_STRATA>& strat) {
            // Gate on tetra_deam_bin<0: covers both !interior_safe (len==30) and
            // mid_total==0 (all-N interior), ensuring bulk == Σstrata by construction.
            // Uses tetra_deam_bin (cross-end proxy) not read_deam_bin (5' proxy).
            if (tetra_deam_bin < 0 || p < 1 || p + 2 >= static_cast<int>(len)) return;
            int i0 = nuc_idx(decoded[p - 1]);
            int i1 = nuc_idx(decoded[p]);
            int i2 = nuc_idx(decoded[p + 1]);
            int i3 = nuc_idx(decoded[p + 2]);
            if (i0 < 0 || i1 < 0 || i2 < 0 || i3 < 0) return;
            const int idx = i0*64 + i1*16 + i2*4 + i3;
            ++bulk[idx];
            if (tetra_deam_bin >= 0) ++strat[tetra_deam_bin][idx];
        };
        for (int p = 1; p <= 4; ++p)
            add_ctx4(p, profile.tetra_5prime_terminal, profile.tetra_5prime_terminal_by_deam);
        for (int p = 10; p <= 14; ++p)
            add_ctx4(p, profile.tetra_5prime_interior, profile.tetra_5prime_interior_by_deam);
    }

    // CpG-like context split — 5' terminal positions (all 15)
    // Also accumulate upstream-context-aware bins (AC, CC, GC, TC)
    for (int p = 0; p < SampleDamageProfile::N_POS && (p + 1) < static_cast<int>(len); ++p) {
        const char x = decoded[p];
        const char y = decoded[p + 1];
        if ((x == 'C' || x == 'T') && (y == 'A' || y == 'C' || y == 'G' || y == 'T')) {
            // F2: capture +2 to split non-CpG into CHG (+2=G) vs CHH. CpG needs only +1.
            // CHG/CHH require the +2 base in-bounds; if +2 is off the read end the site is
            // ambiguous between CHG and CHH and is dropped (not misassigned).
            const char z = ((p + 2) < static_cast<int>(len)) ? decoded[p + 2] : 'N';
            const int ctx = ct_meth_context(y, z);
            if (ctx >= 0) {
                profile.ct_ctx_total_5prime[ctx][p] += 1.0f;
                if (x == 'T') profile.ct_ctx_t_5prime[ctx][p] += 1.0f;
            }
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
                const char z = ((q + 2) < static_cast<int>(len)) ? decoded[q + 2] : 'N';
                const int ctx = ct_meth_context(y, z);
                if (ctx >= 0) {
                    profile.ct_ctx_total_interior[ctx] += 1.0f;
                    if (x == 'T') profile.ct_ctx_t_interior[ctx] += 1.0f;
                }
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
                if (il >= 0 && ir >= 0) profile.oxog16_t[4*il+ir] += 1;
            } else if (b == 'A') {
                int il = enc(rc_base(r)), ir = enc(rc_base(l));
                if (il >= 0 && ir >= 0) profile.oxog16_a_rc[4*il+ir] += 1;
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
            // F5: gt_elig/gt_pos add the G→T (8-oxoG oxidation) proxy track alongside the
            // existing ct track, so the cross-channel co-occurrence shares this one loop and
            // these buffers — no parallel accumulator and no second pass over the read.
            thread_local std::vector<uint8_t> ct_elig_buf, ct_pos_buf,
                                              ag_elig_buf, ag_pos_buf,
                                              gt_elig_buf, gt_pos_buf;
            const size_t W = static_cast<size_t>(wlen);
            if (ct_elig_buf.size() < W) {
                ct_elig_buf.resize(W);
                ct_pos_buf .resize(W);
                ag_elig_buf.resize(W);
                ag_pos_buf .resize(W);
                gt_elig_buf.resize(W);
                gt_pos_buf .resize(W);
            }
            uint8_t* ct_elig = ct_elig_buf.data();
            uint8_t* ct_pos  = ct_pos_buf.data();
            uint8_t* ag_elig = ag_elig_buf.data();
            uint8_t* ag_pos  = ag_pos_buf.data();
            uint8_t* gt_elig = gt_elig_buf.data();
            uint8_t* gt_pos  = gt_pos_buf.data();
            std::memset(ct_elig, 0, W);
            std::memset(ct_pos,  0, W);
            std::memset(ag_elig, 0, W);
            std::memset(ag_pos,  0, W);
            std::memset(gt_elig, 0, W);
            std::memset(gt_pos,  0, W);
            int n_ct = 0, k_ct = 0, n_ag = 0, k_ag = 0, n_gt = 0, k_gt = 0;
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
                // F5 G→T proxy: eligible at G/T sites, positive when the base reads T (a G→T
                // substitution and a real T are indistinguishable reference-free, so the per-read
                // T-among-{G,T} marginal is the composition-corrected null — same device the ct
                // track uses for C→T). Exclude CpG-context G to mirror the ct non-CpG gate.
                if ((b == 'G' || b == 'T') && pv != 'C') {
                    gt_elig[j] = 1; ++n_gt;
                    if (b == 'T') { gt_pos[j] = 1; ++k_gt; }
                }
            }
            if (wlen <= 64) {
                // Pack eligible/position flags into uint64 bitmaps; replace
                // the O(W)×10 branchy j-loop with 10 shifted-AND-popcount ops.
                uint64_t ct_elig_w = 0, ct_pos_w = 0;
                uint64_t ag_elig_w = 0, ag_pos_w = 0;
                uint64_t gt_elig_w = 0, gt_pos_w = 0;
                for (int j = 0; j < wlen; ++j) {
                    ct_elig_w |= static_cast<uint64_t>(ct_elig[j]) << j;
                    ct_pos_w  |= static_cast<uint64_t>(ct_pos[j])  << j;
                    ag_elig_w |= static_cast<uint64_t>(ag_elig[j]) << j;
                    ag_pos_w  |= static_cast<uint64_t>(ag_pos[j])  << j;
                    gt_elig_w |= static_cast<uint64_t>(gt_elig[j]) << j;
                    gt_pos_w  |= static_cast<uint64_t>(gt_pos[j])  << j;
                }
                if (n_ct >= 2) {
                    ++acc.reads_used_ct;
                    const __int128 A_ct = (k_ct >= 2)
                        ? static_cast<__int128>(k_ct) * (k_ct - 1) : 0;
                    const __int128 D_ct = static_cast<__int128>(n_ct) * (n_ct - 1);
                    for (int d = 1; d <= 10; ++d) {
                        uint64_t e     = ct_elig_w & (ct_elig_w >> d);
                        uint64_t pairs = static_cast<uint64_t>(__builtin_popcountll(e));
                        uint64_t obs   = static_cast<uint64_t>(__builtin_popcountll(e & ct_pos_w & (ct_pos_w >> d)));
                        acc.pairs_ct[d] += pairs;
                        acc.obs_ct[d]   += obs;
                        if (A_ct > 0 && pairs > 0) {
                            const __int128 P = static_cast<__int128>(pairs);
                            acc.exp_ct[d] += fp_round_div(P * A_ct, D_ct);
                            acc.var_ct[d] += fp_round_div(P * A_ct * (D_ct - A_ct), D_ct * D_ct);
                        }
                    }
                }
                if (n_ag >= 2) {
                    ++acc.reads_used_ag;
                    const __int128 A_ag = (k_ag >= 2)
                        ? static_cast<__int128>(k_ag) * (k_ag - 1) : 0;
                    const __int128 D_ag = static_cast<__int128>(n_ag) * (n_ag - 1);
                    for (int d = 1; d <= 10; ++d) {
                        uint64_t e     = ag_elig_w & (ag_elig_w >> d);
                        uint64_t pairs = static_cast<uint64_t>(__builtin_popcountll(e));
                        uint64_t obs   = static_cast<uint64_t>(__builtin_popcountll(e & ag_pos_w & (ag_pos_w >> d)));
                        acc.pairs_ag[d] += pairs;
                        acc.obs_ag[d]   += obs;
                        if (A_ag > 0 && pairs > 0) {
                            const __int128 P = static_cast<__int128>(pairs);
                            acc.exp_ag[d] += fp_round_div(P * A_ag, D_ag);
                            acc.var_ag[d] += fp_round_div(P * A_ag * (D_ag - A_ag), D_ag * D_ag);
                        }
                    }
                }
                // F5 cross-channel (G→T × C→T). Ordered pairs at separation d in BOTH directions
                // (G-site upstream of C-site and vice versa); a positive observation is T at the
                // G-site AND T at the C-site. Poisson-independence expectation uses the two per-read
                // marginals p_g·p_c, so a within-read positive correlation (radiation-track cluster)
                // surfaces as obs − exp > 0. Skip when either marginal is undefined (n<2).
                if (n_gt >= 2 && n_ct >= 2) {
                    ++acc.reads_used_x;
                    // pc = (k_gt/n_gt)*(k_ct/n_ct) = PN/PD as an exact rational.
                    const __int128 PN = static_cast<__int128>(k_gt) * k_ct;
                    const __int128 PD = static_cast<__int128>(n_gt) * n_ct;
                    for (int d = 1; d <= 10; ++d) {
                        // direction 1: G-site at j, C-site at j+d
                        uint64_t e1 = gt_elig_w & (ct_elig_w >> d);
                        uint64_t o1 = e1 & gt_pos_w & (ct_pos_w >> d);
                        // direction 2: C-site at j, G-site at j+d
                        uint64_t e2 = ct_elig_w & (gt_elig_w >> d);
                        uint64_t o2 = e2 & ct_pos_w & (gt_pos_w >> d);
                        uint64_t pairs = static_cast<uint64_t>(__builtin_popcountll(e1))
                                       + static_cast<uint64_t>(__builtin_popcountll(e2));
                        uint64_t obs   = static_cast<uint64_t>(__builtin_popcountll(o1))
                                       + static_cast<uint64_t>(__builtin_popcountll(o2));
                        acc.pairs_x[d] += pairs;
                        acc.obs_x[d]   += obs;
                        if (PN > 0 && pairs > 0) {
                            const __int128 P = static_cast<__int128>(pairs);
                            acc.exp_x[d] += fp_round_div(P * PN, PD);
                            acc.var_x[d] += fp_round_div(P * PN * (PD - PN), PD * PD);
                        }
                    }
                }
            } else {
                if (n_ct >= 2) {
                    ++acc.reads_used_ct;
                    const __int128 A_ct = (k_ct >= 2)
                        ? static_cast<__int128>(k_ct) * (k_ct - 1) : 0;
                    const __int128 D_ct = static_cast<__int128>(n_ct) * (n_ct - 1);
                    for (int d = 1; d <= 10; ++d) {
                        uint64_t pairs = 0, obs = 0;
                        for (int j = 0; j + d < wlen; ++j) {
                            if (!ct_elig[j] || !ct_elig[j + d]) continue;
                            ++pairs;
                            obs += static_cast<uint64_t>(ct_pos[j] & ct_pos[j + d]);
                        }
                        acc.pairs_ct[d] += pairs;
                        acc.obs_ct[d]   += obs;
                        if (A_ct > 0 && pairs > 0) {
                            const __int128 P = static_cast<__int128>(pairs);
                            acc.exp_ct[d] += fp_round_div(P * A_ct, D_ct);
                            acc.var_ct[d] += fp_round_div(P * A_ct * (D_ct - A_ct), D_ct * D_ct);
                        }
                    }
                }
                if (n_ag >= 2) {
                    ++acc.reads_used_ag;
                    const __int128 A_ag = (k_ag >= 2)
                        ? static_cast<__int128>(k_ag) * (k_ag - 1) : 0;
                    const __int128 D_ag = static_cast<__int128>(n_ag) * (n_ag - 1);
                    for (int d = 1; d <= 10; ++d) {
                        uint64_t pairs = 0, obs = 0;
                        for (int j = 0; j + d < wlen; ++j) {
                            if (!ag_elig[j] || !ag_elig[j + d]) continue;
                            ++pairs;
                            obs += static_cast<uint64_t>(ag_pos[j] & ag_pos[j + d]);
                        }
                        acc.pairs_ag[d] += pairs;
                        acc.obs_ag[d]   += obs;
                        if (A_ag > 0 && pairs > 0) {
                            const __int128 P = static_cast<__int128>(pairs);
                            acc.exp_ag[d] += fp_round_div(P * A_ag, D_ag);
                            acc.var_ag[d] += fp_round_div(P * A_ag * (D_ag - A_ag), D_ag * D_ag);
                        }
                    }
                }
                // F5 cross-channel (G→T × C→T) — scalar mirror of the bitmap branch above.
                if (n_gt >= 2 && n_ct >= 2) {
                    ++acc.reads_used_x;
                    const __int128 PN = static_cast<__int128>(k_gt) * k_ct;
                    const __int128 PD = static_cast<__int128>(n_gt) * n_ct;
                    for (int d = 1; d <= 10; ++d) {
                        uint64_t pairs = 0, obs = 0;
                        for (int j = 0; j + d < wlen; ++j) {
                            if (gt_elig[j] && ct_elig[j + d]) {
                                ++pairs;
                                obs += static_cast<uint64_t>(gt_pos[j] & ct_pos[j + d]);
                            }
                            if (ct_elig[j] && gt_elig[j + d]) {
                                ++pairs;
                                obs += static_cast<uint64_t>(ct_pos[j] & gt_pos[j + d]);
                            }
                        }
                        acc.pairs_x[d] += pairs;
                        acc.obs_x[d]   += obs;
                        if (PN > 0 && pairs > 0) {
                            const __int128 P = static_cast<__int128>(pairs);
                            acc.exp_x[d] += fp_round_div(P * PN, PD);
                            acc.var_x[d] += fp_round_div(P * PN * (PD - PN), PD * PD);
                        }
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
    for (int i = 0; i < SampleDamageProfile::N_TETRANUC; ++i) {
        dst.tetra_5prime_terminal[i] += src.tetra_5prime_terminal[i];
        dst.tetra_5prime_interior[i] += src.tetra_5prime_interior[i];
        for (int s = 0; s < SampleDamageProfile::N_OX_DEAM_STRATA; ++s) {
            dst.tetra_5prime_terminal_by_deam[s][i] += src.tetra_5prime_terminal_by_deam[s][i];
            dst.tetra_5prime_interior_by_deam[s][i] += src.tetra_5prime_interior_by_deam[s][i];
        }
    }
    for (int i = 0; i < SampleDamageProfile::N_TRINUC; ++i) {
        dst.tri_5prime_terminal[i] += src.tri_5prime_terminal[i];
        dst.tri_5prime_interior[i] += src.tri_5prime_interior[i];
        dst.tri_3prime_terminal[i] += src.tri_3prime_terminal[i];
        dst.tri_3prime_interior[i] += src.tri_3prime_interior[i];
    }
    for (int p = 0; p < SampleDamageProfile::N_POS_TRI; ++p)
        for (int i = 0; i < SampleDamageProfile::N_TRINUC; ++i) {
            dst.tri_5prime_pos[p][i] += src.tri_5prime_pos[p][i];
            dst.tri_3prime_pos[p][i] += src.tri_3prime_pos[p][i];
            for (int s = 0; s < SampleDamageProfile::N_OX_DEAM_STRATA; ++s) {
                dst.tri_5prime_pos_by_deam[s][p][i] += src.tri_5prime_pos_by_deam[s][p][i];
                dst.tri_3prime_pos_by_deam[s][p][i] += src.tri_3prime_pos_by_deam[s][p][i];
            }
        }
    // Cross-fit de-circularized strata accumulators
    for (int s = 0; s < SampleDamageProfile::N_OX_DEAM_STRATA; ++s) {
        auto merge_cf = [&](SampleDamageProfile::CrossFitStrata& a,
                            const SampleDamageProfile::CrossFitStrata& b) {
            a.term_k[s] += b.term_k[s]; a.term_n[s] += b.term_n[s];
            a.intr_k[s] += b.intr_k[s]; a.intr_n[s] += b.intr_n[s];
            a.null_term_k[s] += b.null_term_k[s]; a.null_term_n[s] += b.null_term_n[s];
            a.null_intr_k[s] += b.null_intr_k[s]; a.null_intr_n[s] += b.null_intr_n[s];
        };
        merge_cf(dst.cf_ga3, src.cf_ga3);
        merge_cf(dst.cf_ct3, src.cf_ct3);
        dst.cf_reads[s]     += src.cf_reads[s];
        dst.cf_len_sum[s]   += src.cf_len_sum[s];
        dst.cf_len_sumsq[s] += src.cf_len_sumsq[s];
        dst.cf_igc_num[s]   += src.cf_igc_num[s];
        dst.cf_igc_den[s]   += src.cf_igc_den[s];
    }
    // Per-position terminal pi counts + co-occurrence scalar
    for (int L = 0; L < SampleDamageProfile::N_PI_LEN; ++L)
      for (int C = 0; C < SampleDamageProfile::N_PI_C; ++C) {
        auto mp = [](SampleDamageProfile::PiPosArr& d, const SampleDamageProfile::PiPosArr& s) {
            for (int p = 0; p < SampleDamageProfile::P_PI; ++p) {
                d[p].n_elig += s[p].n_elig; d[p].n_deam += s[p].n_deam; } };
        mp(dst.pi_pos_5prime[L][C],    src.pi_pos_5prime[L][C]);
        mp(dst.pi_pos_3prime_ds[L][C], src.pi_pos_3prime_ds[L][C]);
        mp(dst.pi_pos_3prime_ss[L][C], src.pi_pos_3prime_ss[L][C]);
        dst.pi_cooc[L][C].n        += src.pi_cooc[L][C].n;
        dst.pi_cooc[L][C].sum_k5    += src.pi_cooc[L][C].sum_k5;
        dst.pi_cooc[L][C].sum_k3    += src.pi_cooc[L][C].sum_k3;
        dst.pi_cooc[L][C].sum_k5k3 += src.pi_cooc[L][C].sum_k5k3;
        for (size_t h = 0; h < dst.pi_cooc[L][C].hist.size(); ++h)
            dst.pi_cooc[L][C].hist[h] += src.pi_cooc[L][C].hist[h];
        if (C == 0) {   // interior null is [L]-only; merge once per L
            dst.pi_cooc_interior_ds[L].n        += src.pi_cooc_interior_ds[L].n;
            dst.pi_cooc_interior_ds[L].sum_k5   += src.pi_cooc_interior_ds[L].sum_k5;
            dst.pi_cooc_interior_ds[L].sum_k3   += src.pi_cooc_interior_ds[L].sum_k3;
            dst.pi_cooc_interior_ds[L].sum_k5k3 += src.pi_cooc_interior_ds[L].sum_k5k3;
            for (size_t h = 0; h < dst.pi_cooc_interior_ds[L].hist.size(); ++h)
                dst.pi_cooc_interior_ds[L].hist[h] += src.pi_cooc_interior_ds[L].hist[h];
        }
      }
      {   // prep QC counters are pooled scalars — merge once, outside the [L][C] loop
        dst.cooc_qc_term_n += src.cooc_qc_term_n; dst.cooc_qc_term_a3 += src.cooc_qc_term_a3; dst.cooc_qc_term_t3 += src.cooc_qc_term_t3;
        dst.cooc_qc_null_n += src.cooc_qc_null_n; dst.cooc_qc_null_a3 += src.cooc_qc_null_a3; dst.cooc_qc_null_t3 += src.cooc_qc_null_t3;
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
        da.reads_used_x         += sa.reads_used_x;
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
            da.obs_x[d]    += sa.obs_x[d];
            da.pairs_x[d]  += sa.pairs_x[d];
            da.exp_x[d]    += sa.exp_x[d];
            da.var_x[d]    += sa.var_x[d];
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
    dst.per_read_deam_sum      += src.per_read_deam_sum;
    dst.per_read_deam_sumsq    += src.per_read_deam_sumsq;
    dst.per_read_deam_n        += src.per_read_deam_n;
    dst.per_read_ct5_sum       += src.per_read_ct5_sum;
    dst.per_read_ga3_sum       += src.per_read_ga3_sum;
    dst.per_read_ct5ga3        += src.per_read_ct5ga3;
    dst.per_read_ct5ga3_cpg    += src.per_read_ct5ga3_cpg;
    dst.per_read_n_tpg         += src.per_read_n_tpg;
    dst.per_read_ct5ct3        += src.per_read_ct5ct3;
    dst.per_read_g5_sum        += src.per_read_g5_sum;
    dst.per_read_g3_sum        += src.per_read_g3_sum;
    dst.per_read_g5g3          += src.per_read_g5g3;
    dst.per_read_score_len     += src.per_read_score_len;

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

    for (int b = 0; b < SampleDamageProfile::N_OX_GC_BINS; ++b)
        merge_ox_bins(dst.ox_comov_bins[b], src.ox_comov_bins[b]);

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
    profile.dmax_ct5_chg         = std::numeric_limits<float>::quiet_NaN();
    profile.dmax_ct5_chh         = std::numeric_limits<float>::quiet_NaN();
    profile.meth_select_excess   = std::numeric_limits<float>::quiet_NaN();
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
    profile.oxog16_t.fill(0);
    profile.oxog16_a_rc.fill(0);
    profile.s_oxog_16ctx.fill(0.0f);
    profile.cov_oxog_16ctx.fill(0);
    profile.tetra_5prime_terminal.fill(0);
    profile.tetra_5prime_interior.fill(0);
    for (auto& a : profile.tetra_5prime_terminal_by_deam) a.fill(0);
    for (auto& a : profile.tetra_5prime_interior_by_deam) a.fill(0);
    profile.tri_5prime_terminal.fill(0.0);
    profile.tri_5prime_interior.fill(0.0);
    profile.tri_3prime_terminal.fill(0.0);
    profile.tri_3prime_interior.fill(0.0);
    for (auto& a : profile.tri_5prime_pos) a.fill(0);
    for (auto& a : profile.tri_3prime_pos) a.fill(0);
    for (auto& strat : profile.tri_5prime_pos_by_deam) for (auto& a : strat) a.fill(0);
    for (auto& strat : profile.tri_3prime_pos_by_deam) for (auto& a : strat) a.fill(0);
    profile.cf_ga3 = {};
    profile.cf_ct3 = {};
    profile.cf_reads.fill(0); profile.cf_len_sum.fill(0); profile.cf_len_sumsq.fill(0);
    profile.cf_igc_num.fill(0); profile.cf_igc_den.fill(0);
    profile.pi_pos_5prime = {};
    profile.pi_pos_3prime_ds = {};
    profile.pi_pos_3prime_ss = {};
    profile.pi_cooc = {};
    profile.pi_cooc_interior_ds = {};
    profile.cooc_qc_term_n = profile.cooc_qc_term_a3 = profile.cooc_qc_term_t3 = 0;
    profile.cooc_qc_null_n = profile.cooc_qc_null_a3 = profile.cooc_qc_null_t3 = 0;
    profile.d_max_cooccurrence = DamageEstimate{};

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

    profile.mixture_n_components = 0;
    profile.mixture_d_population = 0.0f;
    profile.mixture_d_damaged = 0.0f;
    profile.mixture_d_population_highgc = 0.0f;
    profile.mixture_pi_damaged = 0.0f;
    profile.mixture_bic = 0.0f;
    profile.mixture_converged = false;
    profile.mixture_identifiable = false;
    profile.gc_histogram.fill(0.0);
    profile.adaptive_gc_threshold = 0.0f;
    profile.gc_threshold_computed = false;
    profile.gc_bins = {};
    profile.ox_comov_bins = {};
    profile.len_bins = {};
    profile.gc_stratified_d_max_weighted = 0.0f;
    profile.gc_stratified_d_max_joint = 0.0f;
    profile.gc_stratified_d_max_peak = 0.0f;
    profile.gc_peak_bin = -1;
    profile.gc_stratified_valid = false;
    profile.pi_damaged = 0.0f;
    profile.d_damaged = 0.0f;
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
    profile.oxo_two_marker       = {};
    profile.oxo_two_marker_ss    = {};
    profile.per_read_deam_sum    = 0.0;
    profile.per_read_deam_sumsq  = 0.0;
    profile.per_read_deam_n      = 0;
    profile.per_read_ct5_sum     = 0.0;
    profile.per_read_ga3_sum     = 0.0;
    profile.per_read_ct5ga3      = 0.0;
    profile.per_read_ct5ga3_cpg  = 0.0;
    profile.per_read_n_tpg       = 0.0;
    profile.per_read_ct5ct3      = 0.0;
    profile.per_read_g5_sum      = 0.0;
    profile.per_read_g3_sum      = 0.0;
    profile.per_read_g5g3        = 0.0;
    profile.per_read_score_len   = 0.0;
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
