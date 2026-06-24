#pragma once

#include <cstdint>
#include <optional>

#include "taph/damage_estimate.hpp"

// Per-read damaged/non-damaged scoring primitive (SOLUTION_pi_delta_dmax.md §6.2). The library owns the
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

// PRIOR-FREE per-read both-end, position-resolved log-likelihood-ratio: log P(obs | damaged) − log P(obs |
// non-damaged) under the VALIDATED δ/c model — amplitude A = D_MAX_CONSERVED imported (not estimated),
// λ/baseline from the fitted per-position decay; NOT the non-identifiable mixture (§6.5(1)). This is the
// RANKING primitive: it returns a finite score for EVERY read (no refuse, no prior) so a consumer can sort
// reads by ancient evidence even before the sample-level π gate fires. For single-stranded libraries the two
// ends share one lesion chemistry (5′ C→T and 3′ C→T are the same deamination event read off opposite
// termini of one molecule), so the ends are NOT independent; a dependence-deflation factor SS_DEP_RHO down-
// weights the second end's contribution to avoid double-counting one molecule's overhang (documented
// approximation — a full joint two-end model is the moment-cube path, not this ranking kernel). DS ends ARE
// independent (separate strands), so no deflation. Returns 0.0 when there are no informative sites.
double damage_evidence_llr(const ReadDamageObs& obs, const SampleDamageProfile& profile);

// PRIOR-AWARE per-read posterior: sigmoid(damage_evidence_llr(obs,profile) + logit(profile.pi.point)).
// Returns std::nullopt UNLESS profile.pi.state == DETECTED — it refuses rather than fold in a π that lies
// at H0 (§6.4). The consumer thresholds the returned posterior in (0,1).
std::optional<double> read_ancient_posterior(const ReadDamageObs& obs, const SampleDamageProfile& profile);

// Backward-compat alias (legacy callers): identical contract to the original prior-free LLR — returns the
// raw log-LR when profile.pi.state == DETECTED and the per-position λ was fitted, else std::nullopt.
// New code should call damage_evidence_llr (ranking) or read_ancient_posterior (gated decision) directly.
std::optional<double> read_ancient_llr(const ReadDamageObs& obs, const SampleDamageProfile& profile);

// Reference-free length-decay constant τ of the per-bin amplitude δ(L)≈A·exp(−L/τ). Profiles χ²(τ) over a
// 1-D grid via closed-form WLS on the live length bins (genuine-zero bins censored at a floor) and returns
// τ̂ + its 95% χ²-profile interval + the 3-state verdict. Replaces the wrong-axis per-position Briggs λ on
// the ref-free path; finalize_tau runs before finalize_pi and gates it via state==DETECTED.
DamageEstimate finalize_tau(const SampleDamageProfile& profile, const TauConfig& cfg = TauConfig{});

// Per-position terminal-decay Briggs fit of the 5′ C→T channel from pi_pos_5prime: closed-form
// (Zhao-2025-style) WLS linearisation of r(p)=baseline+(1−baseline)·A·exp(−λ·p) with binomial
// weights, plus a shape LRT (decay vs flat-null, df=2). Runs alongside finalize_tau (separate
// because it needs no length-bin profiling). finalize_pi gates the pi range on PiShapeFit::detected.
PiShapeFit fit_pi_shape(const SampleDamageProfile& profile, const TauConfig& cfg = TauConfig{});

// Reference-free scission rate γ (bp⁻¹): left-truncated exponential MLE over the right tail of the
// fine fragment-length histogram (peak bin onwards). γ tracks molecular degradation; larger γ → shorter
// molecules. CI from Fisher information (Wald 95%). Returns ScissionEstimate::fitted==false when the
// tail is too sparse (< 3 non-empty tail bins or < 50 tail reads) or has zero excess (S ≤ 0).
ScissionEstimate finalize_scission(const SampleDamageProfile& profile);

// Asymmetric leakage control ε from interior Chargaff contrast T/(T+G)−A/(A+C).
// ε ≈ 0 for pooled-DS symmetric damage (oxidation); large |ε_term| identifies C→T contamination.
EpsilonEstimate finalize_epsilon(const SampleDamageProfile& profile);

// Reference-free GC→AT depletion channel σ₀ = T/(T+G)+A/(A+C)−1 from pooled interior counts.
// NOT an oxidation estimator — σ₀ = composition + deamination + oxidation, algebraically fused.
// Emit alongside gc_interior for cross-sample GC-normalisation; length_slope is a confound QC.
GcDepletionEstimate finalize_gc_depletion(const SampleDamageProfile& profile);

// Primary reference-free oxidation readout: deamination-coupled G→T regression over the
// 256-cell (s1×s2×GC×length) panel. Wraps compute_oxo_two_marker; is_ss inferred from
// profile.library_type. Requires oxo_two_marker.cells to be populated by the fqdup damage pass.
OxoTwoMarkerResult finalize_oxidation_comovement(const SampleDamageProfile& profile);

} // namespace taph
