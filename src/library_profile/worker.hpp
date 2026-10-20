#pragma once
// Private worker queue and per-thread state for the library profiler (moved from fqdup).
#include "taph/damage_profile.hpp"
#include "taph/frame_selector_decl.hpp"
#include "taph/sample_damage_profile.hpp"
#include "taph/length_stratified_profile.hpp"
#include <array>
#include <climits>
#include <condition_variable>
#include <mutex>
#include <string>
#include <vector>

namespace taph {

// LSD sig5 encoding:
//   5' position 0: quinary (A=1,G=2,T=3,C=4,other=0) — carries both the T/C
//     deamination readout AND purine/pyrimidine identity, so a joint
//     deamination+depurination per-read LLR classifier can be reconstructed
//     from the histogram without a second FASTQ pass. Positions 1-4 stay
//     ternary (T=1,C=2,other=0) — depurination is a position-0-only signal
//     (matches compute_purine_stats, which only reads position 0).
//     5 (pos0 quinary) * 3^4 (pos1-4 ternary) = 405 values.
//   3' quinary (A=1,G=2,T=3,C=4,other=0) — 5 positions = 5^5 = 3125 values
// ceiling: CpG context not encoded (needs next-base); upgrade via quinary 5' sig
static constexpr int N_LSD_SIG5_5P   = 405;   // 5 * 3^4
static constexpr int N_LSD_SIG5_3R   = 3125;  // 5^5 quinary
static constexpr int N_LSD_SIG5      = N_LSD_SIG5_5P * N_LSD_SIG5_3R;  // 1265625
// ---- bounded work queue ------------------------------------------------
struct WorkQueue {
    std::mutex              mtx;
    std::condition_variable cv_not_empty;
    std::condition_variable cv_not_full;
    std::vector<std::vector<std::string>> batches;
    bool   done      = false;
    int    max_depth;

    explicit WorkQueue(int max_depth) : max_depth(max_depth) {}

    void push(std::vector<std::string>&& batch) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_not_full.wait(lk, [&]{ return (int)batches.size() < max_depth; });
        batches.push_back(std::move(batch));
        cv_not_empty.notify_one();
    }

    bool pop(std::vector<std::string>& batch) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_not_empty.wait(lk, [&]{ return !batches.empty() || done; });
        if (batches.empty()) return false;
        batch = std::move(batches.back());
        batches.pop_back();
        cv_not_full.notify_one();
        return true;
    }

    void set_done() {
        std::unique_lock<std::mutex> lk(mtx);
        done = true;
        cv_not_empty.notify_all();
    }
};

// ---- Single-pass OxoG + LSD accumulation structs -----------------------
// Ternary 10-position signature histogram cell. Second moments allow
// exact SE computation without re-reading the FASTQ.
struct OxSigCell {
    int64_t t=0, g=0, cv=0, av=0, n=0;  // first moments + count
    int64_t t2=0, tm=0, m2=0;            // Σ T^2, Σ T*M, Σ M^2 for SE
};

static constexpr int N_SIG10 = 59049;  // 3^10: ternary 10-position space


// ---- shared signature tables -------------------------------------------
// One dense copy of the OxoG signature histograms and the LSD sig5 counts per worker
// group. Workers buffer their per-read updates and add them in under the lock, so
// memory is one set of tables plus a small buffer per worker, not a set per worker.
// Every cell is an integer sum, so totals do not depend on thread count or order.
struct SigTables {
    std::vector<OxSigCell> ox_fwd, ox_rev;
    std::vector<int32_t>   lsd_cnt;        // n_prov_bins × N_LSD_SIG5 when LSD is on
    std::mutex             mtx;
    SigTables()
        : ox_fwd(taph::SampleDamageProfile::N_GC_BINS * N_SIG10)
        , ox_rev(taph::SampleDamageProfile::N_GC_BINS * N_SIG10) {}
};

// ---- per-thread state --------------------------------------------------
struct WorkerState {
    taph::SampleDamageProfile           profile;
    std::array<uint64_t, LSD_HIST_BINS> lsd_hist{};
    int64_t reads_scanned = 0;
    int64_t reads_skipped = 0;
    int     len_min       = INT_MAX;
    int     len_max       = 0;
    int64_t len_sum       = 0;
    // LSD single-pass: per-(prov_bin, lsd_sig5) read counts go to tables->lsd_cnt.
    taph::LengthBinStats    lbs;
    bool                    lsd_on = false;
    std::vector<int>        lsd_prov_edges; // provisional bin boundaries (copy)
    SigTables*              tables = nullptr;

    // Pending updates, added to `tables` by flush().
    struct OxPending { uint32_t fwd, rev; int32_t t, g, cv, av; };
    std::vector<OxPending> ox_pending;
    std::vector<uint32_t>  lsd_pending;
    static constexpr size_t kFlushAt = 8192;

    void add_ox(uint32_t fwd, uint32_t rev, int t, int g, int cv, int av) {
        ox_pending.push_back({fwd, rev, t, g, cv, av});
        if (ox_pending.size() >= kFlushAt) flush();
    }
    void add_lsd(uint32_t idx) {
        lsd_pending.push_back(idx);
        if (lsd_pending.size() >= kFlushAt) flush();
    }
    void flush() {
        if (ox_pending.empty() && lsd_pending.empty()) return;
        std::lock_guard<std::mutex> lock(tables->mtx);
        for (const OxPending& u : ox_pending) {
            {   // forward orientation
                auto& cell = tables->ox_fwd[u.fwd];
                const int M = u.t + u.g;
                cell.t += u.t; cell.g += u.g; cell.cv += u.cv; cell.av += u.av; ++cell.n;
                cell.t2 += (int64_t)u.t * u.t; cell.tm += (int64_t)u.t * M; cell.m2 += (int64_t)M * M;
            }
            {   // complement-mapped orientation
                const int Tr = u.av, Gr = u.cv, Cvr = u.g, Avr = u.t;
                auto& cell = tables->ox_rev[u.rev];
                const int M = Tr + Gr;
                cell.t += Tr; cell.g += Gr; cell.cv += Cvr; cell.av += Avr; ++cell.n;
                cell.t2 += (int64_t)Tr*Tr; cell.tm += (int64_t)Tr*M; cell.m2 += (int64_t)M*M;
            }
        }
        for (uint32_t idx : lsd_pending) ++tables->lsd_cnt[idx];
        ox_pending.clear();
        lsd_pending.clear();
    }
};

inline static void worker_ox_accumulate(WorkerState& state, const std::string& seq, int L) {
    const int gc  = taph::SampleDamageProfile::get_gc_bin(seq);
    const int beg = L / 3, end_ = L - (L / 3);
    int T = 0, G = 0, Cv = 0, Av = 0;
    for (int j = beg; j < end_; ++j) {
        const char b = seq[j] & ~0x20u;
        if      (b == 'T') ++T;
        else if (b == 'G') ++G;
        else if (b == 'C') ++Cv;
        else if (b == 'A') ++Av;
    }
    // Compute both ternary sigs in one endpoint pass.
    // ALL reads go into BOTH ox_fwd and ox_rev — no proxy orientation decision.
    // is_ss is unknown during pass 1; oxog_from_sig_hist applies only merged_fwd
    // for SS (all reads, correct fwd counts) and both for DS (q-weighting selects
    // the dominant orientation per cell naturally post-fit).
    uint32_t fsig = 0, rsig = 0, fb = 1, rb = 1;
    const int k = (L < 10) ? L : 10;
    for (int i = 0; i < k; ++i) {
        const char fc = seq[i]     & ~0x20u;
        const char rc = seq[L-1-i] & ~0x20u;
        const int  fv = (fc == 'T') ? 1 : (fc == 'C') ? 2 : 0;
        const int  rv = (rc == 'A') ? 1 : (rc == 'G') ? 2 : 0;
        fsig += fv * fb; fb *= 3;
        rsig += rv * rb; rb *= 3;
    }
    // Fwd histogram: original orientation (correct for SS and DS fwd reads);
    // rev histogram: complement-mapped orientation (correct for DS rev reads).
    state.add_ox(static_cast<uint32_t>(gc * N_SIG10 + fsig), static_cast<uint32_t>(gc * N_SIG10 + rsig),
                 T, G, Cv, Av);
    // LSD sig5: joint (5' pos0-quinary+pos1-4-ternary, 3' quinary) for
    // single-pass reconstruction. Computed independently from fsig (which
    // stays pure ternary for OxoG's own 10-position purposes) because
    // position 0 here needs the full A/G/T/C readout for the joint
    // deamination+depurination classifier (see lsd_llr_from_sig).
    if (state.lsd_on) {
        uint32_t sig5f = 0;
        {
            const char c0 = seq[0] & ~0x20u;
            const int v0 = (c0=='A') ? 1 : (c0=='G') ? 2 : (c0=='T') ? 3 : (c0=='C') ? 4 : 0;
            sig5f = static_cast<uint32_t>(v0);
            uint32_t fb2 = 5;
            const int k1 = (L < 5) ? L : 5;
            for (int i = 1; i < k1; ++i) {
                const char ci = seq[i] & ~0x20u;
                const int vi = (ci=='T') ? 1 : (ci=='C') ? 2 : 0;
                sig5f += static_cast<uint32_t>(vi) * fb2;
                fb2 *= 3;
            }
        }
        // sig5r_4: 3' end with 5-way encoding (A=1, G=2, T=3, C=4, other=0)
        uint32_t sig5r_4 = 0, q5 = 1;
        const int k5 = (L < 5) ? L : 5;
        for (int i = 0; i < k5; ++i) {
            const char rc4 = seq[L-1-i] & ~0x20u;
            const int rv4 = (rc4=='A') ? 1 : (rc4=='G') ? 2 : (rc4=='T') ? 3 : (rc4=='C') ? 4 : 0;
            sig5r_4 += rv4 * q5; q5 *= 5;
        }
        const uint32_t lsd_sig = sig5f * N_LSD_SIG5_3R + sig5r_4;
        int cb = 0;
        const int Lcap = (L < LSD_L_MAX) ? L : LSD_L_MAX;
        for (int e : state.lsd_prov_edges) { if (Lcap >= e) ++cb; else break; }
        state.add_lsd(static_cast<uint32_t>(cb * N_LSD_SIG5 + lsd_sig));
    }
}

inline static void worker_fn(WorkQueue& queue, WorkerState& state) {
    std::vector<std::string> batch;
    while (queue.pop(batch)) {
        for (const std::string& seq : batch) {
            int L = static_cast<int>(seq.size());
            if (L < LSD_L_MIN) {
                // Short-fragment mode: route 16-29 bp reads to the terminal-only short path
                // (update_sample_profile handles the floor + abstain internally). They do NOT
                // feed lsd_hist / lbs / oxog / len stats, so those stay >=30-exact.
                if (state.profile.short_fragment_floor < LSD_L_MIN)
                    taph::FrameSelector::update_sample_profile(state.profile, seq);
                ++state.reads_skipped; continue;
            }
            taph::FrameSelector::update_sample_profile(state.profile, seq);
            ++state.lsd_hist[lsd_hist_bin(L)];
            if (L < state.len_min) state.len_min = L;
            if (L > state.len_max) state.len_max = L;
            state.len_sum += L;
            ++state.reads_scanned;
            if (state.lsd_on)
                state.lbs.update(seq, (L < LSD_L_MAX) ? L : LSD_L_MAX);
            worker_ox_accumulate(state, seq, L);
        }
    }
    state.flush();
}

// ---- paired-end work queue + worker -------------------------------------
struct PairedBatch {
    std::vector<std::string> r1;
    std::vector<std::string> r2;
};

struct PairedWorkQueue {
    std::mutex              mtx;
    std::condition_variable cv_not_empty;
    std::condition_variable cv_not_full;
    std::vector<PairedBatch> batches;
    bool   done      = false;
    int    max_depth;

    explicit PairedWorkQueue(int max_depth) : max_depth(max_depth) {}

    void push(PairedBatch&& batch) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_not_full.wait(lk, [&]{ return (int)batches.size() < max_depth; });
        batches.push_back(std::move(batch));
        cv_not_empty.notify_one();
    }
    bool pop(PairedBatch& batch) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_not_empty.wait(lk, [&]{ return !batches.empty() || done; });
        if (batches.empty()) return false;
        batch = std::move(batches.back());
        batches.pop_back();
        cv_not_full.notify_one();
        return true;
    }
    void set_done() {
        std::unique_lock<std::mutex> lk(mtx);
        done = true;
        cv_not_empty.notify_all();
    }
};

inline static void paired_worker_fn(PairedWorkQueue& queue, WorkerState& state) {
    PairedBatch batch;
    while (queue.pop(batch)) {
        const size_t n = std::min(batch.r1.size(), batch.r2.size());
        for (size_t k = 0; k < n; ++k) {
            // Input is trusted clean (fqdup merge strips 3' poly-G run-through upstream).
            std::string_view s1(batch.r1[k]);
            std::string_view s2(batch.r2[k]);
            int L1 = static_cast<int>(s1.size());
            int L2 = static_cast<int>(s2.size());
            if (L1 < LSD_L_MIN || L2 < LSD_L_MIN) { ++state.reads_skipped; continue; }
            const bool ok = taph::FrameSelector::update_sample_profile_paired(
                state.profile, s1, s2);
            if (!ok) { ++state.reads_skipped; continue; }
            // Do NOT feed lsd_hist: R1 read-length != fragment length for the
            // long inserts that fail to merge, and paired reads do not enter LSD
            // until they get their own top stratum. Polluting merged_lsd_hist
            // would shift the merged-channel bin edges.
            if (L1 < state.len_min) state.len_min = L1;
            if (L1 > state.len_max) state.len_max = L1;
            state.len_sum += L1;
            ++state.reads_scanned;
            // Pool into the length-stratified estimator's top stratum so the primary
            // damage estimate reflects ALL reads (matches the all-reads global-metaDMG anchor).
            if (state.lsd_on) {
                state.lbs.update_paired(s1, s2);
                // Feed the unified per-read LSD scorer: 5' terminus from R1, 3' terminus
                // from complement(R2[0..4]) (R2 reads the bottom strand from the molecule
                // 3' inward, so top_strand_3'[i] = complement(R2[i])). Same sig5 encoding
                // as the SE path so paired fragments enter merged_cnt / reconstruct_lsd_llr_accum
                // identically — the primary pi reflects merged + R1 + R2 under one scorer.
                uint32_t sig5f = 0;
                {
                    const char c0 = s1[0] & ~0x20u;
                    const int v0 = (c0=='A')?1:(c0=='G')?2:(c0=='T')?3:(c0=='C')?4:0;
                    sig5f = static_cast<uint32_t>(v0);
                    uint32_t fb2 = 5;
                    const int k1 = (L1 < 5) ? L1 : 5;
                    for (int i = 1; i < k1; ++i) {
                        const char ci = s1[i] & ~0x20u;
                        const int vi = (ci=='T')?1:(ci=='C')?2:0;
                        sig5f += static_cast<uint32_t>(vi) * fb2;
                        fb2 *= 3;
                    }
                }
                uint32_t sig5r_4 = 0, q5 = 1;
                const int k5 = (L2 < 5) ? L2 : 5;
                for (int i = 0; i < k5; ++i) {
                    char b = s2[i] & ~0x20u;                    // complement R2[i] -> top 3'[i]
                    const char top = (b=='A')?'T':(b=='T')?'A':(b=='C')?'G':(b=='G')?'C':'N';
                    const int rv4 = (top=='A')?1:(top=='G')?2:(top=='T')?3:(top=='C')?4:0;
                    sig5r_4 += rv4 * q5; q5 *= 5;
                }
                const uint32_t lsd_sig = sig5f * N_LSD_SIG5_3R + sig5r_4;
                int cb = 0;
                const int Lcap = (L1 < LSD_L_MAX) ? L1 : LSD_L_MAX;
                for (int e : state.lsd_prov_edges) { if (Lcap >= e) ++cb; else break; }
                state.add_lsd(static_cast<uint32_t>(cb * N_LSD_SIG5 + lsd_sig));
            }
        }
    }
    state.flush();
}

// ---- clip pass: re-profile after stripping 5' and/or 3' adapter stubs ------
inline static void clip_worker_fn(WorkQueue& queue, WorkerState& state,
                           const std::vector<std::string>& stubs5,
                           const std::vector<std::string>& stubs3) {
    std::vector<std::string> batch;
    while (queue.pop(batch)) {
        for (std::string& seq : batch) {
            if (!stubs5.empty() && (int)seq.size() >= 6) {
                for (const auto& stub : stubs5) {
                    if (seq.compare(0, 6, stub) == 0) { seq.erase(0, 6); break; }
                }
            }
            if (!stubs3.empty()) {
                bool trimmed;
                do {
                    trimmed = false;
                    int L = static_cast<int>(seq.size());
                    if (L < 12) break;
                    for (const auto& stub : stubs3) {
                        if (seq.compare(L - 6, 6, stub) == 0) {
                            seq.erase(L - 6, 6);
                            trimmed = true;
                            break;
                        }
                    }
                } while (trimmed);
            }
            int L = static_cast<int>(seq.size());
            if (L < LSD_L_MIN) {
                // Short-fragment mode: route 16-29 bp reads to the terminal-only short path
                // (update_sample_profile handles the floor + abstain internally). They do NOT
                // feed lsd_hist / lbs / oxog / len stats, so those stay >=30-exact.
                if (state.profile.short_fragment_floor < LSD_L_MIN)
                    taph::FrameSelector::update_sample_profile(state.profile, seq);
                ++state.reads_skipped; continue;
            }
            taph::FrameSelector::update_sample_profile(state.profile, seq);
            ++state.lsd_hist[lsd_hist_bin(L)];
            if (L < state.len_min) state.len_min = L;
            if (L > state.len_max) state.len_max = L;
            state.len_sum += L;
            ++state.reads_scanned;
            if (state.lsd_on)
                state.lbs.update(seq, (L < LSD_L_MAX) ? L : LSD_L_MAX);
            worker_ox_accumulate(state, seq, L);
        }
    }
    state.flush();
}

} // namespace taph
