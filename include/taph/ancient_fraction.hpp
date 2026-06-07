#pragma once
#include "taph/sample_damage_profile.hpp"
#include <cstdint>
#include <type_traits>

namespace taph {

// Dimension constants shared between libtaph and callers.
inline constexpr int TAPH_FRAC_N_POS      = 15;
inline constexpr int TAPH_FRAC_N_SOFT_POS =  4;

// Minimal transfer struct — only the fields compute_ancient_fraction reads.
// Caller populates from their per-bin accumulator (e.g. fqdup LsdLlrBinAccum).
// Sort bins by a stable key (e.g. n_damaged) before passing for deterministic
// float reduction across thread counts.
struct AncientFractionBins {
    int64_t n_damaged, n_undamaged;
    double  sw_sum;
    double  sw_t5_anc[TAPH_FRAC_N_SOFT_POS], sw_tc5_anc[TAPH_FRAC_N_SOFT_POS];
    double  sw_h3_anc[TAPH_FRAC_N_SOFT_POS], sw_n3_anc[TAPH_FRAC_N_SOFT_POS];
    int64_t t_5_anc[TAPH_FRAC_N_POS],  tc_5_anc[TAPH_FRAC_N_POS];
    int64_t h_3_anc[TAPH_FRAC_N_POS],  n_3_anc[TAPH_FRAC_N_POS];
    int64_t t_5_mod[TAPH_FRAC_N_POS],  tc_5_mod[TAPH_FRAC_N_POS];
    int64_t h_3_mod[TAPH_FRAC_N_POS],  n_3_mod[TAPH_FRAC_N_POS];
};
static_assert(std::is_trivially_copyable_v<AncientFractionBins>);

// Returned to caller — not serialized to JSON.
struct AncientFractionResult {
    bool  valid;
    bool  leakage_5prime, leakage_3prime;
    float pi;
    float d_anc5, d_anc3, d_mod5, d_mod3;
};

// Compute the ancient/modern fraction decomposition from pre-merged bin data.
//
// bins[0..n_bins-1]: sorted by a stable key for reproducible float sums.
// bg5/bg3:           interior background rates (from lsd_cls_params or similar).
// position_0_artifact_*: bulk artifact flags from dp (passed explicitly to
//                    avoid needing dp access before this function fills it).
//
// Fills dp.damaged_fraction_* and dp.modern_fraction_* fields for JSON.
// Returns a result struct with validity/leakage flags for the caller.
AncientFractionResult compute_ancient_fraction(
    const AncientFractionBins* bins, int n_bins,
    double bg5, double bg3,
    bool position_0_artifact_5prime,
    bool position_0_artifact_3prime,
    SampleDamageProfile& dp);

} // namespace taph
