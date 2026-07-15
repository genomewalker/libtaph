#pragma once

// Validated reference-free ancient-fraction contract types (SOLUTION_pi_delta_dmax.md §6).
// Additive shadow-mode step 1: these coexist with the legacy mixture path; no consumer is wired yet.

#include <cmath>

namespace taph {

// Length-coupling gate: nulls cluster at w_length≈0.5, positives push toward 1 (SOLUTION §6.5(3)).
// Point-gate at 0.6 (matches finalize_dmax C2); the consumer applies the hysteresis band on w's CI.
constexpr double W_LENGTH_GATE = 0.6;

// 3-state confidence (SOLUTION §6.3/§6.4). NOT_DETECTED is RESERVED: the reference-free bulk estimator
// cannot assert "modern" — mapping recovers dilute ancient that delta=0 misses (§5) — so it emits only
// DETECTED or UNDETERMINED. NOT_DETECTED is for a future reference-anchored path.
// LOW_ABUNDANCE: end-asymmetric recovery. Fired ONLY when channel-5 strand symmetry FAILED because
// one terminus is ->G-overcall artifact-dead (artifact_overcall_Xp) AND the OPPOSITE (live) end passes a
// two-model decay LRT (Briggs exp decay vs position-fixed spike/flat null) on its merged per-position
// histogram. pi is then formed from the live end alone (A = live_end_d_max, denom = the library's own
// co-occurrence-fished damaged-class amplitude, NOT a hardcoded cohort constant).
// Distinct from DETECTED (both ends clean, strand-symmetric) and ABSTAIN (no recoverable signal).
// BELOW_FLOOR: detection floor. Entered on the SAME end-asymmetric path as LOW_ABUNDANCE (exactly one
// terminus artifact-dead) but the surviving end FAILS to authenticate as a live decay channel — either
// its two-model decay LRT does not clear significance, or its terminal-window eligible-site count is
// below PI_FLOOR_MIN_ELIG. This is the honest output for a truly-undetectable / artifact-only library:
// no point estimate is identifiable, so only an UPPER BOUND on pi is emitted (pi.hi set; point/lo null).
// UNIDENTIFIED_INSUFFICIENT_LENGTH_STRUCTURE: a NO-CALL (never zero). The per-library terminal amplitude
// A_b = π_dmg·d_dmg is only identifiable when the fragment-length axis has ≥2 populated strata (a single
// stratum cannot separate a short-fragment damage enrichment from a length-flat composition artifact —
// e.g. the low-input blank with all reads in one length bin). We DECLINE to call, we do not assert
// no-damage.
// EXCESS_PRESENT_DECAY_UNAUTHENTICATED: a NO-CALL (never zero). ≥2 length strata exist and a MATERIAL raw
// terminal excess is present on the biological end, but it does NOT decay short→long (flat or rising with
// fragment length). Reference-free, a non-decaying terminal excess cannot be authenticated as deamination
// vs a length-correlated composition artifact — but neither can it be asserted absent (an alignment-based
// caller may still confirm real damage). We DECLINE to call. Distinct from NOT_DETECTED, which is reserved
// for the case where the raw terminal excess is itself ~0 within noise (a genuine absence → zero).
enum class DamageConfidence { DETECTED, TRACE, NOT_DETECTED, UNDETERMINED, ANCIENT_CPG, LOW_ABUNDANCE,
                              BELOW_FLOOR, UNIDENTIFIED_INSUFFICIENT_LENGTH_STRUCTURE,
                              EXCESS_PRESENT_DECAY_UNAUTHENTICATED };

// Bilateral CpG Δ threshold: min(cpgΔ_5′, cpgΔ_3′GA) ≥ this → ANCIENT_CPG when δ-decay model
// returns zero (sharp-spike or dilute-signal profiles that evade the exponential fit).
constexpr double CPG_BILATERAL_ANCIENT_THR = 0.03;

// Shape-LRT threshold for the per-position terminal-decay fit (PiShapeFit): 2·Δloglik ≥ this ⇒ the
// exponential decay shape is distinguishable from a flat constant-rate null. χ²(df=2) 0.99 quantile
// = 9.21 (conservative: pi is only emitted when the terminal channel unambiguously decays).
constexpr double PI_SHAPE_LRT_THR = 9.21;

// Channel-5 strand-symmetry gate (deam_context_spectrum: ct_5prime_excess vs ga_3prime_rc_excess).
// AUTHENTIC ⇔ amplitude (mean signed excess) > this AND mirror residual < 0.5·amplitude. Genuine
// ancient deamination is strand-symmetric (5′C→T mirrors 3′G→A under RC); pervasive artifacts are not.
// ceiling: 0.015 is a placeholder calibrated on n=2 real libraries (FLB03mAds3 authentic amp~+0.027 vs
// FLB57md artifact amp~-0.017 — the sign cleanly separates; the channel-5 per-context-averaged amp is
// ~5x smaller-scale than a summed-position excess, so the original 0.03 clipped genuine damage). upgrade:
// recalibrate via a gargammel damage titration (plan), and fix the amp↔damaged-class-amplitude scale for pi.
constexpr double PI_CH5_AMP_THR = 0.015;

// SS channel-5 strand-asymmetry ceiling. DS deamination is strand-symmetric by mechanism (5' C->T mirrors
// 3' G->A), so ds uses a PROPORTIONAL residual gate (resid < 0.5*amp). SS deamination is intrinsically
// 5'-dominant (5' >> 3'), so the proportional gate is unmeetable for genuine ss; an ABSOLUTE residual
// ceiling is the ss-appropriate form. EMPIRICAL: real ss FLB resid 0.04-0.05 vs ss-blank composition
// artifact resid 0.27-0.34; 0.15 sits cleanly between (3x real, ~half blank-min). Pairs with amp>THR,
// which already rejects the blanks (amp <= 0.011 < 0.015 < real ss amp 0.020-0.028).
constexpr double PI_SS_RESID_MAX = 0.15;

// Terminal ->G overcall artifact threshold (low-abundance cascade degeneracy guard). The raw G-fraction
// (full ACGT composition) at BOTH terminal positions p0 AND p1 must exceed the interior baseline by at
// least this much to register a dead-end ->G overcall (the dead test is min(g_p0,g_p1)-g_int >= THR).
// Measured on correct units (raw ACGT snapshots, commit 91de2f7), the empirical overcall (FLB57md 5',
// VALIDATED) is a SHARP 2-cycle spike: raw G ~0.48/0.43 at p0/p1 then decays to interior (~0.27) by p2 --
// NOT the sustained plateau the old unit-mixed denominator (G/(A+G), excluding T+C) implied. Requiring
// BOTH terminal positions elevated rejects a single-position noise blip. The flag is artifact-specific BY
// CONSTRUCTION: a positive raw-G excess cannot arise from terminal deamination (5' C->T leaves G untouched;
// 3' G->A DEPLETES G), so it never gates out a genuine G-enriched live (G->A) end. ceiling: 0.05
// placeholder calibrated on FLB57md's dead 5' spike; recalibrate via gargammel titration.
constexpr double ARTIFACT_G_SPIKE_THR = 0.05;

// Detection-floor minimum: the live-end terminal window must carry at least this many eligible (C-bearing,
// event-stratum C>0) sites summed over the fitted positions p=0..P_PI-1 before a LOW_ABUNDANCE point/CI may
// be emitted. Below it the decay LRT is statistically powerless (cov~1 terminal windows go noise-dominated),
// so the honest output is BELOW_FLOOR with an upper bound only. live_end_decay_lrt already requires ≥3
// positions at ≥50 elig each to even fit (≈150 floor on the fit itself); this raises the bar to a window
// total that makes the χ²(df=1) ≥6.63 verdict trustworthy. ceiling: 200 placeholder calibrated against the
// per-position ≥50 fit gate; upgrade: recalibrate the (n_elig, power) curve via a gargammel damage titration.
constexpr double PI_FLOOR_MIN_ELIG = 200.0;

// (Removed PI_LOW_ABUNDANCE_MIN_DMAX: the LOW_ABUNDANCE gate is self-calibrating — it requires the
// d_max-SE-propagated lower CI edge to clear PI_THR, i.e. the recovered interval significantly excludes
// zero, rather than a hardcoded amplitude floor.)

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

// ── Reference-free read-split decision policy (shared: fqdup split, dart, fastst) ─────────────
// Maps the validated pi estimate to a per-read LLR threshold shift. The per-read damage LLR is
// Bayes-thresholded at posterior 0.5 ⇔ llr ≥ −log(pi/(1−pi)); shifting the base cut by the prior
// log-odds is the statistically correct operating point (a low-pi library is not flooded with reads
// that barely cross llr=0). NOTE: the split is ENRICHMENT, not a clean partition — the damaged
// set is pi-calibrated IN EXPECTATION, never molecularly pure.
struct SplitPolicy {
    bool   splittable     = false;  // false ⇒ no damaged stratum (blank/ABSTAIN/BELOW_FLOOR): route all undamaged
    double pi             = -1.0;   // operating-point prior (point estimate)
    double pi_lo          = -1.0;
    double pi_hi          = -1.0;
    double log_prior_odds = 0.0;    // effective cut = base_llr_cut − log_prior_odds
    DamageConfidence state = DamageConfidence::UNDETERMINED;
};

// Split-eligible only when pi is a positive, identifiable point: DETECTED / LOW_ABUNDANCE / TRACE /
// ANCIENT_CPG. ABSTAIN-equivalents (NOT_DETECTED, UNDETERMINED) and the upper-bound-only BELOW_FLOOR
// carry no point ⇒ no damaged stratum, so all reads route undamaged.
inline SplitPolicy split_policy(const DamageEstimate& est) {
    SplitPolicy p;
    p.state = est.state;
    p.pi_lo = est.lo;
    p.pi_hi = est.hi;
    const bool eligible_state =
        est.state == DamageConfidence::DETECTED      ||
        est.state == DamageConfidence::LOW_ABUNDANCE ||
        est.state == DamageConfidence::TRACE         ||
        est.state == DamageConfidence::ANCIENT_CPG;
    if (eligible_state && est.point > 0.0) {
        double pi = est.point;
        if (pi < 1e-6) pi = 1e-6; else if (pi > 1.0 - 1e-6) pi = 1.0 - 1e-6;
        p.splittable     = true;
        p.pi             = pi;
        p.log_prior_odds = std::log(pi / (1.0 - pi));
    }
    return p;
}

// Per-position terminal-decay Briggs fit of the 5′ C→T channel from the pi_pos_5prime counts:
// r(p) = baseline + (1−baseline)·A·exp(−λ·p), p = 0..P_PI−1. Fit by Briggs/Zhao-2025 closed-form
// linearisation (WLS on logit-detrended per-position rates with binomial weights); NOT an EM. The
// shape LRT (binomial deviance of the decay model vs a flat constant-rate null, df = 2) gives the
// DETECTED/ABSTAIN verdict that finalize_pi requires before it may emit a pi range. A is the BULK
// terminal amplitude (averaged over modern+ancient reads); pi = A / (the library's own co-occurrence
// damaged-class amplitude) is formed in finalize_pi. All −1 / false when too few eligible sites or the
// fit did not converge.
struct PiShapeFit {
    double A         = -1.0;  // bulk terminal C→T amplitude at p=0 above baseline (∈[0,1])
    double A_se      = -1.0;  // SE of A (binomial delta-method through the WLS)
    double lambda    = -1.0;  // per-position decay constant (bp⁻¹ along terminal window)
    double baseline  = -1.0;  // background mismatch rate b (interior/asymptotic)
    double lrt       = -1.0;  // 2·(loglik_decay − loglik_flat); χ²(df=2)
    bool   detected  = false; // lrt ≥ shape LRT threshold (decay shape distinguishable from flat)
    bool   fitted    = false; // closed-form fit produced a finite {A,λ,baseline}
    int64_t n_elig   = 0;     // total eligible sites over fitted positions — what the LRT was bought with
                              // (lrt scales ~linearly in n_elig, so it is duplicate-purchasable; report it)
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
// null/leakage QC alongside GcDepletionEstimate, NOT an oxidation estimator.
struct EpsilonEstimate {
    double epsilon_floor  = -1.0;  // pooled interior ε = T/(T+G)−A/(A+C); −1 ⇒ n/a
    double epsilon_lo     = -1.0;  // 95% CI lower
    double epsilon_hi     = -1.0;  // 95% CI upper
    double epsilon_term   = -1.0;  // ε at pos 0–2 pooled; large negative = C→T contamination
    double phi_share      = -1.0;  // legacy; superseded by OxoTwoMarkerResult + BurialFingerprint
    int    n_bins         = 0;
    bool   fitted         = false;
};

// Reference-free GC→AT depletion channel from interior counts: σ = T/(T+G)+A/(A+C)−1.
// NOT an oxidation estimator: σ > 0 for oxidation AND AT-rich composition AND pervasive
// deamination (all three push the same direction). Algebraically inseparable from composition
// at the level of marginal single-base counts. Use as a GC-depletion QC channel only;
// the actual oxidation estimator is OxoTwoMarkerResult (the deamination-coupled G→T regression).
// REQUIRED companion: gc_interior for cross-sample GC-normalisation of σ₀.
struct GcDepletionEstimate {
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

// Deamination-coupled G→T regression over the 256-cell (s1×s2×GC×length) interior panel.
// beta1: s1-coupled D slope (C→T 5' marker); beta2: s2-coupled end marker.
//   DS: s2 = G→A 3' marker (minus-strand deamination seen as G→A on the read).
//   SS: s2 = C→T 3' marker (same-strand deamination; no minus strand present).
// markers_consistent: end-symmetry of deamination→D coupling (library-appropriate test).
// consistency_basis records which marker pairing was used for interpretation.
// This is the PRIMARY reference-free per-sample oxidation readout. Computed by
// finalize_oxidation_comovement (wraps compute_oxo_two_marker from library_interpretation.hpp).
enum class OxoConsistencyBasis {
    DS_STRAND_SYMMETRY,  // beta1≈beta2 via 5'C→T vs 3'G→A; valid for DS libraries
    SS_END_SYMMETRY,     // beta1≈beta2 via 5'C→T vs 3'C→T; valid for SS libraries
};
struct OxoTwoMarkerResult {
    double beta1             = 0.0;  // s1-coupled D slope (C→T 5' marker)
    double beta2             = 0.0;  // s2-coupled D slope (library-type-appropriate 3' marker)
    double beta1_se          = 0.0;
    double beta2_se          = 0.0;
    double beta1_z           = 0.0;  // beta1 / beta1_se
    double beta2_z           = 0.0;  // beta2 / beta2_se
    double alpha             = 0.0;  // intercept (global D floor = prep artifact)
    double sigma2            = 0.0;  // residual variance
    bool   markers_consistent = false;
    OxoConsistencyBasis consistency_basis = OxoConsistencyBasis::DS_STRAND_SYMMETRY;
    double delta_beta        = 0.0;  // beta1 - beta2 (≈0 when consistent)
    // B1: propagated SE of delta_beta = sqrt(beta1_se^2 + beta2_se^2) (beta1,beta2 from a shared
    // 5x5 WLS fit; treated as independent — the cross-covariance is a small same-fit term and the
    // existing DS markers_consistent test already uses this same se_delta). Lets a downstream gate
    // ask "does the oxidation contrast significantly exclude null" via |delta_beta|/delta_beta_se.
    double delta_beta_se     = 0.0;
    int    n_cells_used      = 0;
    bool   valid             = false;
};

// Layer-1 burial fingerprint: dimensionless ratios over existing observables that cancel burial age.
// Θ = ln(γ/f0): backbone-scission vs hydrolytic-deamination rate ratio (hydrolytic fragmentation pressure index).
// Confounds pH, temperature, water activity, and deamination saturation — NOT a specific pH meter.
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
