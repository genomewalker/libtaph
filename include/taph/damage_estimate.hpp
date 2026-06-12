#pragma once

// Validated reference-free ancient-fraction contract types (SOLUTION_pi_delta_dmax.md §6).
// Additive shadow-mode step 1: these coexist with the legacy mixture path; no consumer is wired yet.

namespace taph {

// Cohort-conserved per-ancient terminal C→T amplitude (mapping A_b cohort mean; SOLUTION §2).
// Denominator of pi_l = clip(delta_l / D_MAX_CONSERVED, 0, 1) AND the imported amplitude for
// read_ancient_llr (fixed, NOT estimated — keeps the β=π·A confound out of the per-read score, §6.5(1)).
// TODO(later step): expose via CLI/config; documented range 0.34–0.48.
constexpr double D_MAX_CONSERVED = 0.39;

// Length-coupling gate: nulls cluster at w_length≈0.5, positives push toward 1 (SOLUTION §6.5(3)).
// Point-gate at 0.6 (matches finalize_dmax C2); the consumer applies the hysteresis band on w's CI.
constexpr double W_LENGTH_GATE = 0.6;

// 3-state confidence (SOLUTION §6.3/§6.4). NOT_DETECTED is RESERVED: the reference-free bulk estimator
// cannot assert "modern" — mapping recovers dilute ancient that delta=0 misses (§5) — so it emits only
// DETECTED or UNDETERMINED. NOT_DETECTED is for a future reference-anchored path.
enum class DamageConfidence { DETECTED, NOT_DETECTED, UNDETERMINED };

// Gated ancient-fraction estimate: point + 95% interval + state. point/lo/hi == -1 ⇒ not set.
// For profile.tau, point/lo/hi carry τ̂ and its 95% χ²-profile CI; the floor-model decomposition
// {f0, amplitude, overhang_fraction} is populated by finalize_tau alongside τ (default −1 if not fitted).
struct DamageEstimate {
    double point = -1.0;
    double lo    = -1.0;
    double hi    = -1.0;
    DamageConfidence state = DamageConfidence::UNDETERMINED;

    // Floor-model extension: δ(L) = f0 + amplitude·exp(−L/τ). Both non-negative. Populated by
    // finalize_tau; −1 when only 1 live bin available or no data. overhang_fraction = A/(A+f0) is
    // the boundary-robust gate statistic (does not blow up when long-read bins touch δ=0).
    double f0                = -1.0;  // pervasive length-independent C→T floor
    double amplitude         = -1.0;  // overhang component amplitude A
    double overhang_fraction = -1.0;  // A/(A+f0); 1.0 = pure overhang, 0.0 = pure pervasive
    double overhang_lo       = -1.0;  // 95% CI lower (delta method; −1 when projected to boundary)
    double overhang_hi       = -1.0;  // 95% CI upper
};

// Reference-free length-decay constant τ (bp): the e-folding length of the overhang component of
// δ(L) = f0 + A·exp(−L/τ). finalize_tau profiles χ²(τ) over a 1-D grid via closed-form 2-param WLS
// on live bins (f0 and A jointly optimised at each τ, constrained ≥ 0). The hi edge of the 95% CI
// {τ : χ²(τ)−χ²min ≤ 3.84} drives the 3-state gate; overhang_fraction is the boundary-robust supplement.
struct TauConfig {
    double tau_max_detected     = 35.0;   // tau_hi < this → DETECTED
    double tau_max_undetermined = 80.0;
    double a_min                = 0.04;   // amplitude Σδ floor (checked against A+f0)
    int    min_live_bins        = 3;      // bins with CI excluding 0
    double delta_floor          = 0.015;  // censor floor for genuine-zero bins
    double overhang_fraction_min = 0.10;  // minimum overhang_fraction for finalize_pi gate
};

// Reference-free scission rate γ (bp⁻¹): exp(−γ·(L−L_mode)) fit to the right tail of the fragment
// length distribution. L_mode is the mean of the modal fine bin; the right tail is modelled as a
// left-truncated exponential. γ is the molecular nick rate: larger γ → shorter, more degraded DNA.
// MLE: γ = n_tail / Σ n_i·(x̄_i − L_mode).  SE (Fisher): γ / √n_tail.  CI: Wald 95%.
struct ScissionEstimate {
    double gamma        = -1.0;  // nick rate (bp⁻¹); −1 ⇒ not fitted
    double lo           = -1.0;  // 95% CI lower (Wald)
    double hi           = -1.0;  // 95% CI upper (Wald)
    double mean_length  = -1.0;  // read-weighted mean fragment length (bp), all len_bins
    double modal_length = -1.0;  // mean length of the peak fine bin (bp)
    int64_t n_tail      = 0;     // reads in right tail used for γ (peak bin + bins after)
    int64_t n_total     = 0;     // total reads across all fine len_bins
    bool fitted         = false; // false if tail too sparse or S ≤ 0
};

} // namespace taph
