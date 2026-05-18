#pragma once
#include "taph/sample_damage_profile.hpp"
#include "taph/library_interpretation.hpp"
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

} // namespace taph
