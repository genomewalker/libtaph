---
type: libtaph-json-field
title: oxidation_epsilon
tier: standard
estimand: epsilon parameter: excess G→T magnitude from oxidation model
stability: stable
emitted_by: profile_json.cpp
---

### `oxidation_epsilon`

Oxidative damage to guanine produces 8-oxo-7,8-dihydroguanine (8-oxoG), which
mispairs with adenine during replication, generating G→T transversions in sequenced
reads. Unlike hydrolytic deamination, oxidative lesions accumulate from reactive
oxygen species and are not position-restricted to molecule termini; however, the
8-oxoG signal in ancient molecules is confounded by post-extraction oxidation
introduced during library preparation and sequencing. The `oxidation_epsilon`
block captures the aggregate oxidative damage parameter: an estimate of the excess
G→T substitution rate attributable to genuine 8-oxoG lesions after accounting for
the background transversion rate expected in undamaged sequence. This block appears
only when the library type and base-composition data are sufficient to fit the
oxidation model; it is absent (`null` or not present in the JSON) when coverage or
G-content is inadequate for a reliable estimate.

Object structure depends on library type — see JSON output.

| Field | Type | Description |
|-------|------|-------------|
| *(structure not present in sample JSON)* | — | Object structure depends on library type — see JSON output. |

A non-zero `oxidation_epsilon` indicates measurable G→T excess beyond the background
transversion rate, consistent with 8-oxoG accumulation. Values should be interpreted
alongside the deamination signal: elevated oxidation with low deamination may point
to modern contamination or post-extraction oxidation rather than genuine ancient
damage. When `oxidation_epsilon` is near zero or absent, either the molecule pool
lacks appreciable 8-oxoG or coverage over G-containing positions is too low to
resolve the signal from noise. Cross-check this parameter against the `depurination`
block and the `ancient_fraction` estimate before drawing taphonomic conclusions.


---