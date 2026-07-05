#pragma once

#include "frame_selector_decl.hpp"

#include <array>
#include <cstddef>
#include <string_view>
#include <vector>

namespace taph {

struct LengthBinStats {
    static constexpr std::size_t MAX_BINS = 4;

    std::array<SampleDamageProfile, MAX_BINS> profiles = {};
    std::vector<int> edges;
    std::size_t n_bins = 1;
    SampleDamageProfile::LibraryType forced_library_type =
        SampleDamageProfile::LibraryType::UNKNOWN;

    void configure(const std::vector<int>& new_edges);
    void update(std::string_view seq, int length);
    // Unmerged R1/R2: fragment length unknown (didn't merge => longer than 2x readlen),
    // so route to the top (longest) stratum. R1 gives the 5' terminus, RC(R2) the 3'.
    void update_paired(std::string_view r1, std::string_view r2);
    void merge(const LengthBinStats& other);
    // FIX A: the per-bin BulkDamageModel fits inside finalize_sample_profile are INDEPENDENT (each
    // reads/writes only its own profiles[i]), so run them concurrently. n_threads ≤ 0 ⇒ serial.
    void finalize_all(int n_threads = 1);

    std::size_t bin_index(int length) const;
};

}  // namespace taph
