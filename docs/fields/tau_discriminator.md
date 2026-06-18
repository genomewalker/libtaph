---
type: libtaph-json-field
title: tau_discriminator
tier: standard
estimand: discriminant score separating hydrolytic from oxidative damage regimes
stability: stable
emitted_by: profile_json.cpp
---

### `tau_discriminator`

This block reports the reference-free length-decay constant τ (in base pairs) derived from a floor-decomposition model of length-stratified terminal C→T damage. It appears in every profile.json produced by libtaph; all numeric fields are `null` when the fitting conditions are not met (too few live length bins, or insufficient total damage amplitude). The block is currently marked `"provisional": true`, indicating the τ-based gating logic coexists with the legacy mixture path in shadow mode and does not yet drive any downstream consumer.

#### Background

Terminal C→T excess δ(L) measured across read-length bins carries two separable components. The first is an overhang-deamination component that decays exponentially with read length: δ_overhang(L) = A · exp(−L/τ). Because single-stranded overhangs are short and the longest reads have lost them, this component falls off rapidly; a small τ (≲35 bp) is the hallmark of genuine terminus-concentrated deamination. The second is a length-independent pervasive floor f₀ ≥ 0, arising from in-duplex hydrolytic deamination accumulated during burial — bases modified throughout the fragment, not only at the exposed overhang, so δ(L) remains elevated even in long-read bins.

The full model is δ(L) = f₀ + A · exp(−L/τ). At each candidate τ on a 1-bp grid from 10 to 400, f₀ and A are solved jointly via closed-form, non-negativity-constrained two-parameter WLS. The 95% profile confidence interval on τ is the set {τ : χ²(τ) − χ²_min ≤ 3.84} (1 degree of freedom, f₀ and A profiled out). The upper CI edge drives the three-state gate: if τ_hi < 35 bp the signal is classified DETECTED; if τ_hi < 80 bp it is UNDETERMINED; otherwise NOT_DETECTED. An amplitude guard additionally requires A + f₀ > 0.04 across at least three live bins. A secondary upgrade promotes UNDETERMINED → DETECTED when the overhang fraction A/(A + f₀) ≥ 0.10, because a high-f₀ bias can inflate τ while the decomposition still resolves a dominant overhang component.

The boundary-robust gate statistic `overhang_fraction` = A/(A + f₀) is the key supplement to τ. When f₀ = 0 the fraction saturates at 1 (pure overhang deamination, no pervasive floor); when A = 0 it falls to 0 (pure in-duplex pervasive damage with no detectable terminal excess). This ratio avoids the numerical blow-up that occurs when long-read bins push A toward zero.

| Field | Type | Description |
|-------|------|-------------|
| `provisional` | boolean | Always `true` — block is in shadow mode; gating output does not yet drive downstream consumers |
| `state` | string | Three-state classification: `"DETECTED"`, `"NOT_DETECTED"`, or `"UNDETERMINED"` |
| `point` | number \| null | τ̂ (bp): best-fit e-folding length of the overhang component; null when no live bins are available |
| `ci_lo` | number \| null | Lower edge of the 95% χ²-profile CI on τ̂ (bp) |
| `ci_hi` | number \| null | Upper edge of the 95% χ²-profile CI on τ̂ (bp); state gate is applied to this value |
| `f0` | number \| null | Pervasive length-independent C→T floor at the best-τ solution; null when only one live bin is available |
| `amplitude` | number \| null | Overhang-component amplitude A at the best-τ solution; null when unfitted |
| `overhang_fraction` | number \| null | A/(A + f₀): boundary-robust gate statistic; 1.0 = pure overhang, 0.0 = pure pervasive floor; null when total amplitude is zero |
| `overhang_ci_lo` | number \| null | Delta-method 95% CI lower bound on `overhang_fraction`; null when either A or f₀ is on the non-negativity boundary |
| `overhang_ci_hi` | number \| null | Delta-method 95% CI upper bound on `overhang_fraction`; null when either A or f₀ is on the non-negativity boundary |
| `note` | string | Embedded field-key glossary: "tau = overhang e-fold (bp); f0 = pervasive floor; overhang_fraction = A/(A+f0) boundary-robust gate" |

#### Interpretation

A `state` of `DETECTED` with a small τ̂ (single digits to low tens of bp) and `overhang_fraction` near 1.0 indicates terminal deamination concentrated at short, terminus-proximal positions consistent with single-stranded overhang chemistry — a strong authenticity signal for ancient reads. When `f0` is appreciable relative to `amplitude`, the pervasive floor records accumulated in-duplex hydrolysis rather than purely terminal modification; very high `f0` with low `amplitude` produces large τ̂ and typically yields NOT_DETECTED or UNDETERMINED despite elevated total C→T. Wide or rail-pinned CI edges (`ci_lo` = `ci_hi`) indicate that the likelihood surface is flat — usually caused by too few length bins with resolved damage signal, and the `overhang_fraction` should be checked as the more stable summary. Because the block is provisional, all classification decisions should be cross-referenced against the `bulk_damage` channel and the per-length-bin deamination profiles rather than used as a standalone verdict.

---