#pragma once
#include "taph/sample_damage_profile.hpp"
#include "taph/library_interpretation.hpp"
#include <climits>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace taph {

// ── Pre-scan ──────────────────────────────────────────────────────────────────
//
// Scan up to `prescan_reads` sequences to detect adapter stubs.
// `batch_fn` is called repeatedly; it should append up to `max` sequences to
// `out` and return true, or return false (leaving `out` empty) when exhausted.
// The function reads from the beginning of the source — callers must ensure
// batch_fn restarts at the file head on each invocation.
//
// Returns the finalized pre-scan profile + detected AdapterStubs.

struct PrescanResult {
    SampleDamageProfile profile;
    AdapterStubs        stubs;
};

using BatchFn = std::function<bool(std::vector<std::string>& out, size_t max)>;

PrescanResult run_prescan(BatchFn          batch_fn,
                          SampleDamageProfile::LibraryType forced_lib
                              = SampleDamageProfile::LibraryType::UNKNOWN,
                          int    min_read_len   = 30,
                          size_t prescan_reads  = 100000);

// ── Per-read adapter clipping ─────────────────────────────────────────────────
//
// Encode stubs5/stubs3 vectors into 12-bit hex codes for O(k) per-read lookup.
std::vector<uint32_t> encode_stub_codes(const std::vector<std::string>& stubs);

// Strip matching 5' and/or 3' hex stubs from seq.
// Returns the clipped string, or an empty string if the result is shorter than
// min_len (so callers can skip length check separately).
std::string clip_adapters(const std::string&           seq,
                          const std::vector<uint32_t>& stub5_codes,
                          const std::vector<uint32_t>& stub3_codes,
                          int                          min_len = 30);

// ── Whole-source profiling ─────────────────────────────────────────────────────
//
// profile_reads() runs the complete single-end / merged damage-profiling loop:
//   1. adapter pre-scan (delegates to run_prescan) -> AdapterStubs
//   2. multi-threaded full pass: per-read adapter clip -> update_sample_profile,
//      with an optional caller hook per full-length read for side channels
//   3. merge per-thread profiles -> finalize_sample_profile
// The returned profile is FINALIZED and ready to read.
//
// One engine, one result: fqdup `profile` and dart `predict` self-profiling both
// call this, so the finalized SampleDamageProfile is byte-identical for identical
// reads + identical ProfileConfig (integer-count accumulators + associative merge).

// Factory returning a BatchFn that reads sequences from the START of the source.
// profile_reads() invokes it up to twice: once for the pre-scan (unless stubs are
// supplied in cfg), once for the full pass. Each returned BatchFn MUST start at the
// head of the source and append up to `max` sequences to `out`, returning true, or
// return false (leaving `out` empty) when exhausted.
using ReaderFactory = std::function<BatchFn()>;

// Called once per FULL-LENGTH read (clipped length L >= cfg.min_read_len) AFTER
// adapter clipping, from worker thread `tid` (0 <= tid < cfg.threads). Use it to
// accumulate side channels (OxoG histograms, length-stratified counts, ...) into
// caller-owned, per-thread state indexed by `tid`. Do NOT touch a SampleDamageProfile
// here -- profile_reads owns that. Touch only the slot for `tid` (no locks needed).
using PerReadHook = std::function<void(int tid, const std::string& seq, int L)>;

struct ProfileConfig {
    int    threads              = 1;
    SampleDamageProfile::LibraryType forced_lib
        = SampleDamageProfile::LibraryType::UNKNOWN;
    int    min_read_len         = 30;         // full-length threshold (fqdup LSD_L_MIN)
    int    short_fragment_floor = 30;         // set on EVERY per-thread profile before
                                              // the first update (load-bearing)
    int    max_read_len         = 0;          // 0 = off; drop reads longer than this
    size_t prescan_reads        = 1000000;    // 0 = scan all; qualifying-read cap for
                                              // the internal run_prescan
    // Optional: skip the internal pre-scan and use these stubs verbatim (dart supplies
    // the stubs it already computed for its Pass-2 per-read clip). nullptr = pre-scan.
    const AdapterStubs* precomputed_stubs = nullptr;
};

struct ProfileResult {
    SampleDamageProfile profile;      // merged + FINALIZED
    AdapterStubs        stubs;        // detected (or the supplied precomputed_stubs)
    int64_t reads_scanned = 0;        // clipped L >= min_read_len  (full-length)
    int64_t reads_skipped = 0;        // clipped L <  min_read_len  (short/empty)
    int     len_min = INT_MAX;        // over full-length reads only
    int     len_max = 0;
    int64_t len_sum = 0;
    int64_t n_stub5_hits         = 0; // UNCLIPPED 5' hexamer == a stub5 code
    int64_t n_stub3_hits         = 0; // UNCLIPPED 3' hexamer == a stub3 code
    int64_t n_stub_reads_checked = 0; // reads where the stub tally ran
};

ProfileResult profile_reads(const ReaderFactory& make_reader,
                            const ProfileConfig& cfg,
                            const PerReadHook&   hook = {});

} // namespace taph
