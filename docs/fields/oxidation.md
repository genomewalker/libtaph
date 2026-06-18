---
type: libtaph-json-field
title: oxidation
tier: standard
estimand: combined oxidation object (GC-depletion, epsilon, scission, strand asymmetry)
stability: stable
emitted_by: profile_json.cpp
---

### `oxidation` block — GC-depletion channel (composition + deamination + oxidation composite)

The statistic σ = T/(T+G) + A/(A+C) − 1 measures the combined deviation from Chargaff
strand-balance symmetry. Under Chargaff's first rule, complementary base pairing in
double-stranded DNA forces T ≈ A and G ≈ C at the whole-genome level, so T/(T+G) ≈
A/(A+C) and σ ≈ 0. Any process that preferentially converts G or C to another base breaks
this symmetry: cytosine deamination elevates T/(T+G) at read termini; guanine oxidation
(G→T via 8-oxoguanine) elevates T/(T+G) throughout the read; base-compositional bias sets
a non-zero background. The `oxidation` block reports the interior baseline value of σ and
its length and terminal stratifications so that downstream analysis can decompose these
contributions. The block is present for all library types whenever the per-base coverage is
sufficient to fit the length-stratified model (`fitted: true`); it is absent (or
`fitted: false`) when coverage is too sparse.

**Important:** σ is a composite signal. The `note` field embedded in the block documents
this explicitly — it is not a standalone oxidation estimator. Use `oxo_two_marker` (two-marker
regression against ancient-fraction strata) and `oxog_estimate` (uniformity-based 8-oxoG
rate) for quantitative oxidation inference. The `oxidation` block provides the
GC-normalisation baseline: residualizing σ against `gc_interior` across a sample set
removes systematic composition offsets before oxidation signal is interpreted.

| Field | Type | Description |
|-------|------|-------------|
| `fitted` | boolean | `true` when the length-stratified fit converged; `false` when coverage was insufficient |
| `sigma0` | float | Interior σ = T/(T+G) + A/(A+C) − 1; composition + deamination + oxidation baseline estimated from the interior (non-terminal) fraction of reads |
| `sigma0_se` | float | Standard error on `sigma0` |
| `gc_interior` | float | Interior GC fraction (G+C)/(A+C+G+T); used to residualize σ across samples with differing base composition |
| `sigma_term` | float | σ estimated at read termini; elevated above `sigma0` when terminal deamination is present |
| `sigma_long` | float | σ estimated for the longest read-length stratum; used to assess whether σ varies with length, as expected when deamination or oxidation co-varies with fragment size |
| `delta_sigma` | float | `sigma_term` − `sigma0`; positive values indicate terminal enrichment attributable primarily to deamination; negative values may indicate SS library orientation effects |
| `length_slope` | float | Change in σ per unit increase in read length (slope of the length-stratified fit); near-zero for uniform oxidation, non-zero when damage co-varies with fragmentation |
| `length_slope_se` | float | Standard error on `length_slope` |
| `n_bins` | integer | Number of length × GC bins used in the fit |
| `n_counts` | integer | Total base observations contributing to the fit |
| `note` | string | Embedded plain-language description of what σ measures and what it does not; preserved verbatim in the JSON for downstream users |

`sigma0` values near zero indicate a base composition close to Chargaff balance; strongly
positive values in the interior (not attributable to terminal deamination via `delta_sigma`)
may indicate oxidative G→T accumulation, but must be cross-checked against `oxog_estimate`
and `oxo_two_marker` before concluding oxidation is present. For GC-rich or GC-poor
samples, `gc_interior` is the appropriate covariate for normalisation: regress `sigma0` on
`gc_interior` across samples before interpreting cross-sample differences. A large
`delta_sigma` with near-zero `length_slope` is the typical pattern for ancient deamination
without oxidation; large `length_slope` with modest `delta_sigma` warrants further
investigation of fragmentation-correlated composition shifts.


---