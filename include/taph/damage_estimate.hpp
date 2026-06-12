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

// Asymmetric leakage control: ε = T/(T+G) − A/(A+C) (the DIFFERENCE, which cancels in a pooled
// DS library for symmetric oxidation). ε ≈ 0 is the expected value for 8-oxoG; a large |ε| at
// terminal positions flags C→T contamination or strand-asymmetric processes. Used as a
// null/leakage QC alongside OxidationEstimate, NOT as an oxidation estimator.
struct EpsilonEstimate {
    double epsilon_floor  = -1.0;  // pooled interior ε = T/(T+G)−A/(A+C); −1 ⇒ n/a
    double epsilon_lo     = -1.0;  // 95% CI lower
    double epsilon_hi     = -1.0;  // 95% CI upper
    double epsilon_term   = -1.0;  // ε at pos 0–2 pooled; large negative = C→T contamination
    double phi_share      = -1.0;  // legacy; superseded by OxidationEstimate + BurialFingerprint
    int    n_bins         = 0;
    bool   fitted         = false;
};

// Reference-free GC→AT depletion pressure from interior counts. σ = T/(T+G)+A/(A+C)−1 (the SUM;
// opposed to the cancelled difference ε). σ > 0 for both oxidation AND high AT content; it is
// therefore composition-confounded in absolute terms.
// REQUIRED companion: gc_interior. Consumers MUST residualize σ₀ against gc_interior across
// samples to obtain the oxidation component. σ₀ alone is only interpretable as a relative,
// cross-sample measure (higher burial / more damage → higher σ₀ holding community GC constant).
// length_slope is a confound diagnostic (within-sample slope ≈ 0 expected; large slope = length-
// stratified composition bias or deamination-length coupling — flag, do not gate on it).
struct OxidationEstimate {
    double sigma0          = -1.0;  // pooled interior σ₀ = T/(T+G)+A/(A+C)−1; composition-confounded
    double sigma0_se       = -1.0;  // SE from binomial delta method on pooled counts
    double gc_interior     = -1.0;  // (G+C)/(A+T+G+C) pooled interior; MUST accompany sigma0
    double sigma_term      = -1.0;  // σ at pos 0–2 pooled; σ_term ≫ σ₀ = terminal C→T contamination
    double sigma_long      = -1.0;  // σ of longest contributing bin; approximates intrinsic composition
    double delta_sigma     = -1.0;  // sigma0 − sigma_long: partial composition correction; >0 = damage excess
    double length_slope    = -1.0;  // WLS dσ_l/d(bin_mid) [per bp]; confound QC only
    double length_slope_se = -1.0;
    int    n_bins          = 0;     // bins with min(T+G,A+C)_interior ≥ 200
    uint64_t n_counts      = 0;     // total interior bases across contributing bins
    bool   fitted          = false;
};

// Layer-1 burial fingerprint: dimensionless ratios over existing observables that cancel burial age.
// Θ = ln(γ/f0): backbone-scission vs hydrolytic-deamination rate ratio — mostly a pH proxy.
// overhang_fraction = A/(A+f0): copied from tau.
// phi_share = σ₀/(σ₀+f0): upper bound on oxidation fraction (includes composition baseline; use
// only as a relative measure across samples with similar gc_interior).
// All −1 when constituent estimates are not fitted. Layer-2 (T,pH decode) deferred.
struct BurialFingerprint {
    double theta             = -1.0;  // ln(γ/f0); −1 when γ≤0 or f0≤0
    double theta_lo          = -1.0;  // delta-method CI lower
    double theta_hi          = -1.0;
    double overhang_fraction = -1.0;  // A/(A+f0) from tau
    double tau_hat           = -1.0;  // τ̂ (bp) from tau
    double phi_share         = -1.0;  // σ₀/(σ₀+f0); upper bound; −1 when σ₀≤0 or f0<0
    bool   fitted            = false;
};

} // namespace taph
