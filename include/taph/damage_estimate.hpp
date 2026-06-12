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
struct DamageEstimate {
    double point = -1.0;
    double lo    = -1.0;
    double hi    = -1.0;
    DamageConfidence state = DamageConfidence::UNDETERMINED;
};

// Reference-free length-decay constant τ (bp): the e-folding length of the per-bin terminal-deamination
// amplitude δ(L) ≈ A·exp(−L/τ). Genuine terminal damage decays fast (small τ); a pervasive per-base
// artifact gives δ flat-or-rising in L (τ→∞). finalize_tau profiles χ²(τ) over a 1-D grid via closed-form
// WLS on the live (identified, positive, non-boundary) length bins, with genuine-zero bins censored at a
// floor. The hi edge of the 95% CI {τ : χ²(τ)−χ²min ≤ 3.84} drives the 3-state gate.
struct TauConfig {
    double tau_max_detected = 35.0;   // hi < this → DETECTED
    double tau_max_undetermined = 80.0;
    double a_min = 0.04;              // amplitude Σδ floor
    int min_live_bins = 3;            // bins with CI excluding 0
    double delta_floor = 0.015;       // censor floor for genuine-zero bins
};

} // namespace taph
