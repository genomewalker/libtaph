---
type: libtaph-json-field
title: preservation
tier: standard
estimand: composite taphonomic DNA preservation metrics
stability: stable
emitted_by: profile_json.cpp
---

### preservation

The `preservation` block aggregates all damage-channel evidence into a single composite index of how strongly a library exhibits ancient-DNA damage signatures. It is always present in profile JSON output; there is no null/absent case, but `score` and `evidence` will be near zero and `label` will be `"insufficient"` when the read count is too low for a reliable estimate.

#### Composite scores

| Field | Type | Description |
|-------|------|-------------|
| `score` | float [0, 1] | Final preservation index: `evidence × reliability`, clamped to [0, 1]. The primary single-number summary for ranking or filtering. |
| `evidence` | float [0, 1] | Weighted geometric mean of the four damage factors (`f5`, `f3`, `f_coh`, `f_cpg`). Reflects damage magnitude independent of data-quality penalties. |
| `reliability` | float [0, 1] | Product of three continuous gate factors: `g_N` (read-count sufficiency; half-point at ~500 reads), `g_fit` (bulk-model convergence status), and `g_ox` (oxidation-artifact penalty). Values near 1.0 indicate no quality penalty was applied. |
| `label` | string | Categorical summary derived from `authenticity_eff`, not `score`. Values: `"modern-like"` (< 0.10), `"weak"` (0.10–0.30), `"moderate"` (0.30–0.55), `"ancient"` (≥ 0.55). |

#### Damage factor sub-components

Each factor is a logistic or geometric transform of the underlying damage estimate, designed so 0.5 maps to a half-saturating signal and values accumulate multiplicatively in `evidence`.

| Field | Type | Description |
|-------|------|-------------|
| `f5` | float [0, 1] | 5′ terminal C→T damage factor. Logistic sigmoid of `d_max_5prime` with half-point at 5 % damage and scale 0.04 (i.e., `σ((d5 − 0.05)/0.04)`). |
| `f3` | float [0, 1] | 3′ terminal damage factor. Computed analogously to `f5` from `d_max_3prime`. For DS libraries where 3′ G→A is displaced inward by trimming artifact and `d_max_3prime < 0.03` despite `d_max_5prime > 0.05`, the value is imputed from `f5` (DS deamination is symmetric). |
| `f_coh` | float [0, 1] | Mixture coherence factor: geometric mean of a damage-times-length-coupling signal and a 5′/3′ symmetry term. Penalises libraries where terminal damage is present but the length-coupling discriminator (`bulk_damage.w_length`) is near the null cluster (0.5). For SS libraries or DS libraries with censored 3′ signal, symmetry is treated as perfect. |
| `f_cpg` | float [0, 1] | CpG methylation-context factor. Logistic sigmoid of `log₂(CpG_dmax / nonCpG_dmax)` with half-point at log₂ ratio = 1. Set to the uninformative prior (0.3) when the ratio is not evaluable: near-saturation damage (`d_max_combined > 0.20`), insufficient CpG terminal coverage, or NaN ratio. A value of 0.3 therefore means "uninformative" rather than "evidence against antiquity". |

#### Authenticity summary

| Field | Type | Description |
|-------|------|-------------|
| `authenticity_eff` | float [0, 1] | Weighted average of saturating transforms of `d_max_5prime` and `d_max_3prime` (half-saturation at 5 %). DS libraries weight 5′ at 2× and 3′ at 0.5×; SS libraries weight equally. Replaces the earlier mixture-model `pi` term, which was not reliably identifiable from reference-free data and over-reported authenticity at the null. Determines `label`. |
| `authenticity_evidence` | float [0, 1] | Multi-source evidence aggregate combining z-score evidence from both termini, bulk-model damage magnitude weighted by `w_length`, and CpG terminal-coverage-weighted CpG z-score evidence. Represents how many independent channels converge on the authenticity call. |
| `d5_hexamer_corrected` | float | Present only when adapter artifact suppressed the position-0 C→T rate and a hexamer-based correction was applied. Reports the back-extrapolated 5′ damage estimate (`d5 × exp(λ × fit_offset)`); the uncorrected value is visible as `d_max_5prime` in the deamination block. Absent from output when no correction was applied. |

#### Oxidation-like signal sub-block

The nested `oxidation` object reports a reference-free, length- and GC-stratified composition contrast between deamination-weighted ancient-like reads and background-like reads. It is an empirical proxy for oxidative base modification; a positive excess does not prove a specific lesion without orthogonal controls.

| Field | Type | Description |
|-------|------|-------------|
| `oxidation_eff` | float | Top-level shorthand equal to `oxidation.excess_rate.estimate`. Zero when the contrast is anti-oxidative (clamped). |
| `oxidation_evidence` | float [0, 1] | Equal to `oxidation.reliability_score`, the continuous reliability measure for the oxidation channel. |
| `oxidation.raw_rate` | object | `{estimate, ci95_low, ci95_high}` — the G→T / C→A compositional signal rate in the damage-enriched read stratum (unclamped, can be negative). |
| `oxidation.control_rate` | object | `{estimate, ci95_low, ci95_high}` — matching rate in the background-like stratum used as a composition reference. |
| `oxidation.excess_rate` | object | `{estimate, ci95_low, ci95_high}` — `raw_rate − control_rate`, floor-clamped at zero. Anti-oxidative libraries report `{0, 0, 0}`. |
| `oxidation.z_score` | float | Signed test statistic for the raw-vs-control contrast (unclamped). Negative values indicate anti-oxidative direction and are distinguishable from exact-zero. |
| `oxidation.bins_used` | float | Number of length/deamination-stratum bins that contributed to the estimate. |
| `oxidation.effective_bins` | float | Variance-adjusted effective bin count; lower than `bins_used` when bins are highly heterogeneous. |
| `oxidation.heterogeneity` | float | Measure of between-bin variance relative to within-bin variance. Values near 1.0 indicate high inter-bin variation. |
| `oxidation.reliability_score` | float [0, 1] | Continuous reliability of the oxidation estimate; penalised by low bin count, artifact suspicion, and heterogeneity. |
| `oxidation.reliability` | string | Categorical reliability verdict. `"pass"`: channel operational, positive contrast. `"negative"`: channel operational, but contrast is anti-oxidative (z ≤ 0). `"warning"`: marginal bin coverage (2–5 bins, ≥ 1.5 effective). `"fail"`: insufficient bins for any estimate. |

#### Quality-risk summary

| Field | Type | Description |
|-------|------|-------------|
| `qc_risk_eff` | float [0, 1] | Composite QC risk score: weighted mean of hexamer-artifact risk (penalised when the 5′ hexamer T/C excess is not explained by damage), general damage-artifact flag, and adapter-remnant flags for both termini. Values near zero are desirable. |
| `qc_evidence` | float [0, 1] | Evidence weight for the QC risk estimate; reflects how much data supported the risk call. |

#### Interpretation guidance

`score` is the recommended field for ranking or threshold-based filtering; it is calibrated so that a library with strong terminal damage, good length-coupling, and no quality penalties scores above 0.60. `authenticity_eff` and `label` are more interpretable for reporting: they summarise terminal damage magnitude on a scale that is not deflated by reliability penalties, so a low-coverage but genuinely damaged library can carry a non-zero `label` while `score` appropriately reflects the uncertainty via the `reliability` gate. When `f_cpg` is 0.3 and `authenticity_evidence` is low, interpret the block cautiously — the composite may be driven by a single channel. A `qc_risk_eff` above 0.4 warrants inspection of the `library_qc` block and the deamination offset flags before concluding the preservation signal is biological.

---