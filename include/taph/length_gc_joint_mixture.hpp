#pragma once

// length_gc_joint_mixture — shared-component 2-Gaussian mixture fit over
// length × GC cells of a finalized taph::LengthBinStats. Produces one
// d_ancient (shared across all cells) plus per-cell posterior ancient-weight.

#include "taph/joint_damage_model.hpp"
#include "taph/length_stratified_profile.hpp"
#include "taph/sample_damage_profile.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace taph {

struct LengthGcJointMixtureResult {
    static constexpr int N_GC_BINS = SampleDamageProfile::N_GC_BINS;

    double d_ancient    = 0.0;
    double pi_ancient   = 0.0;
    double d_population = 0.0;  // c_sites-weighted mean over all cells
    bool   converged    = false;
    bool   separated    = false;

    // Per-cell posterior ancient-weight w_anc(b,g) = P(damaged | cell).
    // Empty if the fit did not produce a separated damaged component.
    // rows = length bins (size == stats.n_bins), cols = N_GC_BINS.
    // -1.0 marks cells that were invalid / had insufficient C-sites.
    std::vector<std::array<double, N_GC_BINS>> cell_w_ancient;
};

// Fit a 2-component mixture (undamaged μ=0 vs damaged μ=d_ancient) over the
// flat length × GC cell grid of `stats`. Each cell contributes its d_max and
// c_sites. Returns a single shared d_ancient / pi_ancient plus per-cell
// posterior weights. Call after stats.finalize_all().
inline LengthGcJointMixtureResult fit_length_gc_joint_mixture(const LengthBinStats& stats) {
    LengthGcJointMixtureResult out;
    constexpr int  G         = LengthGcJointMixtureResult::N_GC_BINS;
    constexpr std::size_t MAX_CELLS =
        LengthBinStats::MAX_BINS * static_cast<std::size_t>(G);

    // Trust nothing about caller-supplied `n_bins`: clamp to the static
    // MAX_BINS so we never index past stats.profiles[]. Walking past the
    // fixed-size storage on a bad input was the original failure mode.
    const int n_bins = std::clamp(static_cast<int>(stats.n_bins), 0,
                                  static_cast<int>(LengthBinStats::MAX_BINS));
    if (n_bins <= 0) return out;

    std::array<GCBinInput, MAX_CELLS> cells{};
    for (std::size_t i = 0; i < MAX_CELLS; ++i) {
        cells[i].d_max = 0.0f;
        cells[i].c_sites = 0.0f;
        cells[i].valid = false;
    }
    for (int b = 0; b < n_bins; ++b) {
        const auto& pr = stats.profiles[b];
        for (int g = 0; g < G; ++g) {
            std::size_t idx = static_cast<std::size_t>(b) * G + g;
            const auto& gb = pr.gc_bins[g];
            cells[idx].d_max   = gb.valid ? gb.d_max : 0.0f;
            cells[idx].c_sites = static_cast<float>(gb.c_sites);
            cells[idx].valid   = gb.valid && gb.c_sites > 0;
        }
    }

    auto jr = DamageMixtureModel::fit(cells);
    out.d_ancient    = jr.d_ancient;
    out.pi_ancient   = jr.pi_ancient;
    out.d_population = jr.d_mean;
    out.converged    = jr.converged;
    out.separated    = jr.separated;

    if (!jr.separated) return out;

    // Per-cell posterior ancient-weight from the fitted 2-component Gaussian.
    // Shared fixed parameters track DamageMixtureModel's constants.
    constexpr double MU_0    = 0.0;
    constexpr double TAU_0   = 0.01;
    constexpr double TAU_1   = 0.10;
    constexpr double S_FLOOR = 0.02;

    const double mu_1 = jr.d_ancient;
    const double pi_1 = jr.pi_ancient;
    out.cell_w_ancient.assign(n_bins, {});
    for (int b = 0; b < n_bins; ++b) {
        auto& row = out.cell_w_ancient[b];
        for (int g = 0; g < G; ++g) {
            std::size_t idx = static_cast<std::size_t>(b) * G + g;
            const auto& c = cells[idx];
            if (!c.valid || c.c_sites < 1.0f) { row[g] = -1.0; continue; }
            double d = c.d_max;
            double dc = std::clamp(d, 1e-3, 1.0 - 1e-3);
            double sigma_binom = std::sqrt(dc * (1.0 - dc) / c.c_sites);
            double sigma = std::max(S_FLOOR, sigma_binom);
            double var_0 = TAU_0 * TAU_0 + sigma * sigma;
            double var_1 = TAU_1 * TAU_1 + sigma * sigma;
            double ll_0  = -0.5 * std::log(var_0) - 0.5 * (d - MU_0) * (d - MU_0) / var_0;
            double ll_1  = -0.5 * std::log(var_1) - 0.5 * (d - mu_1) * (d - mu_1) / var_1;
            double lp0   = std::log(std::max(1.0 - pi_1, 1e-10)) + ll_0;
            double lp1   = std::log(std::max(pi_1,        1e-10)) + ll_1;
            double lm    = std::max(lp0, lp1);
            double w     = std::exp(lp1 - (lm + std::log(std::exp(lp0 - lm) + std::exp(lp1 - lm))));
            row[g] = w;
        }
    }
    return out;
}

}  // namespace taph
