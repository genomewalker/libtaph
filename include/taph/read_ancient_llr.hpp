#pragma once

#include <cstdint>
#include <optional>

#include "taph/damage_estimate.hpp"

// Per-read ancient/modern scoring primitive (SOLUTION_pi_delta_dmax.md §6.2). The library owns the
// likelihood; the consumer owns the prior (π) and the split threshold. Additive shadow-mode step 1.

namespace taph {

struct SampleDamageProfile;  // full definition in sample_damage_profile.hpp

// One informative terminal site on a read. `deaminated` = the damage allele was observed:
//   5' end → C→T (T observed);  3' end ds → G→A (A observed);  3' end ss → C→T (T observed).
// The channel is end-determined, so the producer chooses informative sites per end and a single
// `deaminated` flag suffices (no per-site channel tag). `pos` is 0-indexed from the read terminus.
struct TerminalMismatch {
    int  pos = 0;
    bool deaminated = false;
};

// Non-owning view of a read's terminal observations — sufficient statistics only, NEVER a FASTQ/read
// type (that coupling would break the reuse test the whole boundary is built on, §6.2).
struct ReadDamageObs {
    int                     length  = 0;
    const TerminalMismatch* five    = nullptr;
    std::uint32_t           n_five  = 0;
    const TerminalMismatch* three   = nullptr;
    std::uint32_t           n_three = 0;
};

// Prior-free per-read log-likelihood-ratio: log P(obs | ancient) − log P(obs | modern), evaluated under
// the VALIDATED δ/c model — amplitude A = D_MAX_CONSERVED imported (not estimated), λ/baseline from the
// fitted decay; NOT the non-identifiable mixture (§6.5(1)). The consumer forms
//   posterior = sigmoid(*llr + logit(profile.pi.point))
// and thresholds. Returns std::nullopt when profile.pi.state != DETECTED — it refuses rather than emit a
// number that lies at H0 (§6.4).
std::optional<double> read_ancient_llr(const ReadDamageObs& obs, const SampleDamageProfile& profile);

// Reference-free length-decay constant τ of the per-bin amplitude δ(L)≈A·exp(−L/τ). Profiles χ²(τ) over a
// 1-D grid via closed-form WLS on the live length bins (genuine-zero bins censored at a floor) and returns
// τ̂ + its 95% χ²-profile interval + the 3-state verdict. Replaces the wrong-axis per-position Briggs λ on
// the ref-free path; finalize_tau runs before finalize_pi and gates it via state==DETECTED.
DamageEstimate finalize_tau(const SampleDamageProfile& profile, const TauConfig& cfg = TauConfig{});

} // namespace taph
