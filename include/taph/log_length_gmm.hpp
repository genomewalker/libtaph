#pragma once

// log_length_gmm — 1D Gaussian-mixture natural-bin detection on a
// log-length histogram. Pairs naturally with taph::LengthBinStats to
// auto-detect length bin edges for stratified damage estimation.

#include <cstdint>
#include <vector>

namespace taph {

struct LogLengthGmmResult {
    std::vector<int> edges;
    int n_components = 1;
    double bic = 0.0;
    bool converged = false;
};

LogLengthGmmResult detect_log_length_gmm_edges(
    const std::vector<uint64_t>& histogram,
    double log_min,
    double log_max,
    int min_length,
    int max_length,
    int max_components = 4);

std::vector<int> detect_quantile_length_edges(
    const std::vector<uint64_t>& histogram,
    double log_min,
    double log_max,
    int min_length,
    int max_length,
    int n_bins);

}  // namespace taph
