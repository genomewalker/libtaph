---
type: libtaph-json-field
title: oxo_two_marker
tier: standard
estimand: dual-marker 5-prime G→T and 3-prime C→A OxoG strand-balance verification
stability: stable
emitted_by: profile_json.cpp
---

### `oxo_two_marker`

The `oxo_two_marker` block tests whether an elevated interior G→T signal is caused by genuine oxidative damage — specifically 8-oxoguanine (8-oxoG) accumulated during burial — or whether it is a structural artefact of the deamination strand architecture. It is present in all profile outputs (both double-stranded and single-stranded library types); it is absent only when the library could not be typed at all.

#### Background

8-oxoG arises from reactive oxygen species attacking guanine at the C8 position and is the dominant oxidative base lesion in aged organic material. The modified base mispairs with adenine during replication, producing G→T transversions. In a reference-free context this signal must be extracted from the interior of reads (away from terminal positions where deamination and strand-end artefacts dominate) and corrected for GC composition, which modulates the baseline G/T ratio independently of damage.

The correction is implemented as a two-predictor ordinary least-squares regression of the interior Chargaff D statistic — the signed asymmetry G→T relative to T→G, composition-normalised — on two deamination-derived markers drawn from the ssDNA-overhang (ancient) stratum: the 5′ C→T rate (`beta1`) and the 3′ G→A rate (`beta2`). The intercept (`alpha`) estimates the GC-composition baseline. `delta_beta` is `beta1 − beta2`, the signed difference that isolates G→T excess after removing the symmetric deamination contribution. Genuine burial oxidation predicts `beta1 ≈ beta2 > 0` and `markers_consistent = true`. When `beta1 > 0` and `beta2 < 0` (opposite signs), the interior G→T signal is tracking the deamination strand polarity, not an independent oxidative process.

| Field | Type | Description |
|-------|------|-------------|
| `valid` | boolean | `true` when sufficient interior counts and both overhang-marker denominators exist for the regression to run |
| `n_cells_used` | integer | Number of (position, context) cells included in the regression; minimum threshold enforced before `valid` is set |
| `alpha` | float | Regression intercept; captures the GC-composition baseline contribution to interior Chargaff D |
| `beta1` | float | Regression coefficient for the 5′ C→T ancient marker; positive under genuine oxidative damage |
| `beta1_se` | float | Standard error of `beta1` |
| `beta1_z` | float | Z-score (`beta1 / beta1_se`); values ≥ 3 are conventionally treated as statistically significant |
| `beta2` | float | Regression coefficient for the 3′ G→A ancient marker; expected to be positive and close to `beta1` under genuine 8-oxoG |
| `beta2_se` | float | Standard error of `beta2` |
| `beta2_z` | float | Z-score (`beta2 / beta2_se`) |
| `delta_beta` | float | `beta1 − beta2`; positive values indicate G→T interior excess above the composition-corrected deamination baseline; the primary quantitative oxidation readout |
| `markers_consistent` | boolean | `true` when `beta1` and `beta2` have the same sign, satisfying the necessary condition for genuine burial oxidation rather than a deamination strand-structure artefact |

#### Interpretation

`delta_beta > 0` with `beta1_z ≥ 3` is the first condition for claiming an interior G→T signal; without `markers_consistent = true` this is insufficient to attribute the signal to oxidative damage. A pattern of `beta1 > 0, beta2 < 0` (opposite signs, `markers_consistent = false`) is strong evidence that the interior asymmetry is driven by the deamination strand polarity and not by an independent oxidative process. Double-stranded libraries in genuinely damaged sediment characteristically show `delta_beta` in the range +0.04–0.06 with `beta1_z` near or above the counting ceiling; single-stranded libraries near zero. A positive correlation between `delta_beta` and a reducing geochemical index, or a negative correlation with an oxidising index, argues against a sample-preparation OxoG artefact (which would run in the opposite direction). When `valid = false` the entire block should be treated as missing.


---