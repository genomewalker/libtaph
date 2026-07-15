#include "taph/damage_profiler.hpp"
#include "taph/frame_selector_decl.hpp"
#include <algorithm>
#include <array>
#include <climits>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace taph {

PrescanResult run_prescan(BatchFn          batch_fn,
                          SampleDamageProfile::LibraryType forced_lib,
                          int    min_read_len,
                          size_t prescan_reads) {
    SampleDamageProfile profile;
    profile.forced_library_type = forced_lib;

    std::array<uint32_t, 4096> hex3{};
    uint64_t n_hex3 = 0;
    size_t   done   = 0;

    // prescan_reads == 0 means "scan the entire source" (documented contract on
    // ProfileConfig::prescan_reads and dart --profile-sample). Map it to the max so
    // the bound checks below are unbounded rather than degenerate (0 < 0 == false).
    const size_t cap = prescan_reads ? prescan_reads : ~size_t(0);

    std::vector<std::string> batch;
    bool more = true;
    while (more && done < cap) {
        batch.clear();
        size_t want = std::min(size_t(50000), cap - done);
        more = batch_fn(batch, want);
        for (const auto& seq : batch) {
            if ((int)seq.size() < min_read_len) continue;
            FrameSelector::update_sample_profile(profile, seq);
            int L = (int)seq.size();
            if (L >= 6) {
                int code = encode_hex_at(seq, L - 6);
                if (code >= 0) { ++hex3[code]; ++n_hex3; }
            }
            if (++done >= cap) break;
        }
    }

    FrameSelector::finalize_sample_profile(profile);
    AdapterStubs stubs = detect_adapter_stubs(profile, hex3.data(), n_hex3);
    return {std::move(profile), std::move(stubs)};
}

std::vector<uint32_t> encode_stub_codes(const std::vector<std::string>& stubs) {
    std::vector<uint32_t> codes;
    codes.reserve(stubs.size());
    for (const auto& s : stubs) {
        int code = encode_hex_at(s, 0);
        if (code >= 0) codes.push_back(static_cast<uint32_t>(code));
    }
    return codes;
}

std::string clip_adapters(const std::string&           seq,
                          const std::vector<uint32_t>& stub5_codes,
                          const std::vector<uint32_t>& stub3_codes,
                          int                          min_len) {
    std::string s = seq;

    if (!stub5_codes.empty() && (int)s.size() >= 6) {
        int code = encode_hex_at(s, 0);
        for (uint32_t sc : stub5_codes) {
            if (code == (int)sc) { s = s.substr(6); break; }
        }
    }
    // 3': ITERATIVE terminal hex-code stripping — promoted from fqdup clip_worker_fn so
    // every caller (fqdup full pass, dart Pass-2 clip, profile_reads) strips stacked 3'
    // stubs identically. The L < 12 guard stops a 6 bp remnant being chased below a full
    // hexamer (mirrors fqdup's `L < 12` break).
    if (!stub3_codes.empty()) {
        bool trimmed;
        do {
            trimmed = false;
            int L = (int)s.size();
            if (L < 12) break;
            int code = encode_hex_at(s, L - 6);
            for (uint32_t sc : stub3_codes) {
                if (code == (int)sc) { s = s.substr(0, L - 6); trimmed = true; break; }
            }
        } while (trimmed);
    }

    // min_len <= 0 -> "pure clip": return the clipped bytes even below a floor
    // (profile_reads / the short-fragment path need them). min_len > 0 keeps the legacy
    // drop-to-empty contract for dart's Pass-2 length gate.
    if (min_len > 0 && (int)s.size() < min_len) return {};
    return s;
}

namespace {

// Bounded batch queue — 1:1 with fqdup WorkQueue (damage_worker.hpp).
struct SeqQueue {
    std::mutex              mtx;
    std::condition_variable cv_ne, cv_nf;
    std::vector<std::vector<std::string>> batches;
    bool done = false;
    int  max_depth;
    explicit SeqQueue(int d) : max_depth(d) {}
    void push(std::vector<std::string>&& b) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_nf.wait(lk, [&]{ return (int)batches.size() < max_depth; });
        batches.push_back(std::move(b));
        cv_ne.notify_one();
    }
    bool pop(std::vector<std::string>& b) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_ne.wait(lk, [&]{ return !batches.empty() || done; });
        if (batches.empty()) return false;
        b = std::move(batches.back());
        batches.pop_back();
        cv_nf.notify_one();
        return true;
    }
    void set_done() {
        std::unique_lock<std::mutex> lk(mtx);
        done = true;
        cv_ne.notify_all();
    }
};

struct ThreadAccum {
    int64_t reads_scanned = 0, reads_skipped = 0;
    int     len_min = INT_MAX, len_max = 0;
    int64_t len_sum = 0;
    int64_t h5 = 0, h3 = 0, hc = 0;   // unclipped stub-hit tallies
};

constexpr int BATCH_SZ = 8192;   // matches fqdup profile.cpp BATCH_SZ

} // namespace

ProfileResult profile_reads(const ReaderFactory& make_reader,
                            const ProfileConfig& cfg,
                            const PerReadHook&   hook) {
    const int nt = std::max(1, cfg.threads);
    ProfileResult result;

    // ---- pre-scan (or accept caller-supplied stubs) -----------------------
    if (cfg.precomputed_stubs) {
        result.stubs = *cfg.precomputed_stubs;
    } else {
        BatchFn scan_fn = make_reader();                       // fresh reader @ head
        PrescanResult pr = run_prescan(scan_fn, cfg.forced_lib,
                                       cfg.min_read_len, cfg.prescan_reads);
        result.stubs = std::move(pr.stubs);
    }
    const bool do_clip = result.stubs.adapter_clipped || result.stubs.adapter3_clipped;
    const bool have_stubs = !result.stubs.stubs5.empty() || !result.stubs.stubs3.empty();
    const std::vector<uint32_t> stub5 = encode_stub_codes(result.stubs.stubs5);
    const std::vector<uint32_t> stub3 = encode_stub_codes(result.stubs.stubs3);

    // ---- per-thread profiles ---------------------------------------------
    // LOAD-BEARING: set short_fragment_floor (and forced_library_type) on EACH
    // per-thread profile BEFORE its first update_sample_profile. update reads the
    // floor to decide, on the very first read, whether to run the terminal-only
    // short path or abstain. Setting it after the first update diverges silently.
    std::vector<std::unique_ptr<SampleDamageProfile>> tp(nt);
    std::vector<ThreadAccum> acc(nt);
    for (int t = 0; t < nt; ++t) {
        tp[t] = std::make_unique<SampleDamageProfile>();
        tp[t]->forced_library_type  = cfg.forced_lib;
        tp[t]->short_fragment_floor = cfg.short_fragment_floor;
    }

    // ---- full pass --------------------------------------------------------
    SeqQueue queue(2 * nt);
    std::vector<std::thread> workers;
    workers.reserve(nt);
    for (int t = 0; t < nt; ++t) {
        workers.emplace_back([&, t]() {
            SampleDamageProfile& prof = *tp[t];
            ThreadAccum&         a    = acc[t];
            std::vector<std::string> batch;
            while (queue.pop(batch)) {
                for (std::string& raw : batch) {
                    const int L0 = (int)raw.size();
                    if (cfg.max_read_len > 0 && L0 > cfg.max_read_len) continue;

                    // Unclipped stub-hit tally (drives JSON read fractions). Runs
                    // whenever stubs exist, independent of do_clip. Code-compare is
                    // equivalent to fqdup's seq.compare(0,6,stub) — an N in the 6-mer
                    // gives encode_hex_at == -1 (no match), same as a string mismatch.
                    if (have_stubs && L0 >= 6) {
                        const int c5 = encode_hex_at(raw, 0);
                        const int c3 = encode_hex_at(raw, L0 - 6);
                        bool m5 = false, m3 = false;
                        if (c5 >= 0) for (uint32_t c : stub5) if ((uint32_t)c5 == c) { m5 = true; break; }
                        if (c3 >= 0) for (uint32_t c : stub3) if ((uint32_t)c3 == c) { m3 = true; break; }
                        a.h5 += m5; a.h3 += m3; ++a.hc;
                    }

                    // Clip via the unified clip_adapters. min_len = 0 → "pure clip"
                    // mode: return the clipped bytes even below the floor, so the
                    // short-fragment terminal path still sees them (see Section 3).
                    std::string seq = do_clip ? clip_adapters(raw, stub5, stub3, 0)
                                              : std::move(raw);
                    const int L = (int)seq.size();

                    if (L < cfg.min_read_len) {
                        // Short/empty after clip. Route 16..(min-1) bp reads to the
                        // terminal-only short path ONLY when the floor is lowered;
                        // update_sample_profile enforces the floor + abstain itself.
                        if (prof.short_fragment_floor < cfg.min_read_len)
                            FrameSelector::update_sample_profile(prof, seq);
                        ++a.reads_skipped;
                        continue;
                    }

                    FrameSelector::update_sample_profile(prof, seq);
                    if (L < a.len_min) a.len_min = L;
                    if (L > a.len_max) a.len_max = L;
                    a.len_sum += L;
                    ++a.reads_scanned;
                    if (hook) hook(t, seq, L);     // side channels: OxoG / LSD / etc.
                }
            }
        });
    }

    // ---- producer: one forward pass, batched -----------------------------
    try {
        BatchFn feed = make_reader();              // fresh reader @ head
        std::vector<std::string> src, batch;
        batch.reserve(BATCH_SZ);
        while (feed(src, BATCH_SZ)) {
            for (auto& s : src) {
                batch.push_back(std::move(s));
                if ((int)batch.size() == BATCH_SZ) {
                    queue.push(std::move(batch));
                    batch.clear();
                    batch.reserve(BATCH_SZ);
                }
            }
            src.clear();
        }
        if (!batch.empty()) queue.push(std::move(batch));
    } catch (...) {
        queue.set_done();
        for (auto& w : workers) w.join();
        throw;
    }
    queue.set_done();
    for (auto& w : workers) w.join();

    // ---- merge (thread-index order) + finalize ---------------------------
    result.profile.forced_library_type  = cfg.forced_lib;
    result.profile.short_fragment_floor = cfg.short_fragment_floor;
    for (int t = 0; t < nt; ++t) {
        FrameSelector::merge_sample_profiles(result.profile, *tp[t]);
        const ThreadAccum& a = acc[t];
        result.reads_scanned        += a.reads_scanned;
        result.reads_skipped        += a.reads_skipped;
        result.len_sum              += a.len_sum;
        result.n_stub5_hits         += a.h5;
        result.n_stub3_hits         += a.h3;
        result.n_stub_reads_checked += a.hc;
        if (a.len_min < result.len_min) result.len_min = a.len_min;
        if (a.len_max > result.len_max) result.len_max = a.len_max;
    }
    result.profile.forced_library_type = cfg.forced_lib;   // survive the merge

    FrameSelector::finalize_sample_profile(result.profile, nt);
    return result;
}

} // namespace taph
