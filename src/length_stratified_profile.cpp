#include "taph/length_stratified_profile.hpp"

#include <algorithm>
#include <stdexcept>

namespace taph {

void LengthBinStats::configure(const std::vector<int>& new_edges) {
    if (new_edges.size() >= MAX_BINS) {
        throw std::invalid_argument("LengthBinStats supports at most 4 bins");
    }
    if (!std::is_sorted(new_edges.begin(), new_edges.end())) {
        throw std::invalid_argument("LengthBinStats edges must be sorted");
    }
    for (std::size_t i = 1; i < new_edges.size(); ++i) {
        if (new_edges[i - 1] >= new_edges[i]) {
            throw std::invalid_argument("LengthBinStats edges must be strictly increasing");
        }
    }

    edges = new_edges;
    n_bins = edges.size() + 1;
    profiles = {};
    for (std::size_t i = 0; i < n_bins; ++i) {
        profiles[i].forced_library_type = forced_library_type;
    }
}

std::size_t LengthBinStats::bin_index(int length) const {
    const auto it = std::upper_bound(edges.begin(), edges.end(), length);
    return static_cast<std::size_t>(std::distance(edges.begin(), it));
}

void LengthBinStats::update(std::string_view seq, int length) {
    if (length < 30 || n_bins == 0) {
        return;
    }
    const std::size_t idx = bin_index(length);
    profiles[idx].forced_library_type = forced_library_type;
    FrameSelector::update_sample_profile(profiles[idx], seq);
}

void LengthBinStats::merge(const LengthBinStats& other) {
    if (n_bins != other.n_bins || edges != other.edges) {
        throw std::invalid_argument("LengthBinStats merge requires identical bin edges");
    }
    for (std::size_t i = 0; i < n_bins; ++i) {
        FrameSelector::merge_sample_profiles(profiles[i], other.profiles[i]);
        profiles[i].forced_library_type = forced_library_type;
    }
}

void LengthBinStats::finalize_all() {
    for (std::size_t i = 0; i < n_bins; ++i) {
        profiles[i].forced_library_type = forced_library_type;
        if (profiles[i].n_reads > 0) {
            FrameSelector::finalize_sample_profile(profiles[i]);
        }
    }
}

}  // namespace taph
