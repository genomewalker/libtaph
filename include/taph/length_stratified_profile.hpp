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
    void finalize_all();

    std::size_t bin_index(int length) const;
};

}  // namespace taph
