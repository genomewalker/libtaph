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
    void merge(const LengthBinStats& other);
    // FIX A: the per-bin BulkDamageModel fits inside finalize_sample_profile are INDEPENDENT (each
    // reads/writes only its own profiles[i]), so run them concurrently. n_threads ≤ 0 ⇒ serial.
    void finalize_all(int n_threads = 1);

    std::size_t bin_index(int length) const;
};

}  // namespace taph
