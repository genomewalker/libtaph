#include "taph/length_stratified_profile.hpp"

#include <algorithm>
#include <atomic>
#include <new>
#include <stdexcept>
#include <thread>
#include <vector>

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
    // SampleDamageProfile is ~650 KB; any assignment creates a stack temporary
    // of that size (or 4× for the full array). Use placement new to reinitialise
    // in-place — no stack allocation, safe in SSH sessions with reduced ulimit -s.
    for (std::size_t i = 0; i < MAX_BINS; ++i) {
        profiles[i].~SampleDamageProfile();
        ::new (static_cast<void*>(&profiles[i])) SampleDamageProfile{};
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

void LengthBinStats::finalize_all(int n_threads) {
    // forced_library_type must be set before any worker may touch profiles[i] (no data race on it).
    for (std::size_t i = 0; i < n_bins; ++i) {
        profiles[i].forced_library_type = forced_library_type;
    }

    // FIX A: each bin's finalize_sample_profile (incl. its BulkDamageModel MLE fit) reads/writes only
    // profiles[i] — independent, no shared mutable state — so distribute bins across a fixed worker
    // pool. Results land in their own indexed slot, so output order is identical to the serial loop
    // regardless of completion order ⇒ byte-identical. n_threads ≤ 1 (or ≤1 bin) ⇒ serial.
    const int workers = std::min<int>(n_threads < 1 ? 1 : n_threads,
                                      static_cast<int>(n_bins));
    if (workers <= 1) {
        for (std::size_t i = 0; i < n_bins; ++i) {
            if (profiles[i].n_reads > 0) {
                FrameSelector::finalize_sample_profile(profiles[i]);
            }
        }
        return;
    }

    std::atomic<std::size_t> next{0};
    auto run = [&] {
        for (;;) {
            std::size_t i = next.fetch_add(1, std::memory_order_relaxed);
            if (i >= n_bins) break;
            if (profiles[i].n_reads > 0) {
                FrameSelector::finalize_sample_profile(profiles[i]);
            }
        }
    };
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(workers) - 1);
    for (int t = 1; t < workers; ++t) pool.emplace_back(run);
    run();
    for (auto& th : pool) th.join();
}

}  // namespace taph
