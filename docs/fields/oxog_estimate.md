---
type: libtaph-json-field
title: oxog_estimate
tier: standard
estimand: 8-oxoguanine per-context oxidative G→T terminal excess rate
stability: stable
emitted_by: profile_json.cpp
---

### `oxog_estimate`

This block reports a reference-free estimate of terminal oxidative G→T artifact from
8-oxoguanine (8-oxoG). It is present when the library accumulates enough G-bearing reads
for the estimator to converge; it is absent (`null` or missing) for very low-coverage
inputs or when the library type classifier cannot assign a strand model. An absent block
should be treated as uninformative rather than as evidence of absence of oxidation.

#### Background

Guanine is the most electron-rich DNA base and the preferred target of reactive oxygen
species (ROS). The primary oxidation product, 8-oxo-2′-deoxyguanosine (8-oxoG), is
a highly mutagenic lesion: it pairs with adenine rather than cytosine during replication,
producing G→T transversions. In a sequencing library prepared from degraded material,
8-oxoG accumulated post-mortem is read as G→T (or, on the complementary strand, C→A).
Unlike deamination-driven C→T transitions, which are enriched at molecule termini because
hydrolysis attacks the overhanging single-stranded flap produced by exonucleolytic
nick extension, genuine 8-oxoG oxidation arising from bulk chemical damage distributes
relatively uniformly across read positions. A terminal spike in G→T would instead indicate
a library-preparation artifact introduced during end-repair or fill-in chemistry, where
the 3′ overhang provides a transient single-stranded substrate for ROS.

The estimator separates these two scenarios. It quantifies the raw G→T strand asymmetry
parameter θ (the fraction of G/T observations that are G→T rather than T→G) at interior
read positions, then tests whether terminal positions show an excess beyond what interior
θ predicts. A Chargaff-parity or trinucleotide-composition control adjusts for the
local sequence context in which G residues occur.

#### Fields

| Field | Type | Description |
|-------|------|-------------|
| `oxo_schema` | `integer` | Version tag for the estimator output layout; used for forward-compatibility parsing |
| `ox_theta` | `float` | G→T strand-asymmetry parameter at terminal positions: proportion of G/T observations that are G→T; 0.5 = symmetric (no net oxidation); < 0.5 = T→G dominant |
| `ox_theta_ci_lo` | `float` | Lower bound of the 95 % confidence interval on `ox_theta` |
| `ox_theta_ci_hi` | `float` | Upper bound of the 95 % confidence interval on `ox_theta` |
| `ox_theta_interior` | `float` | Same strand-asymmetry parameter estimated from interior (non-terminal) positions; serves as the composition-corrected baseline |
| `ox_like_excess` | `float` | Estimated fractional excess of terminal G→T above the interior baseline; 0.0 = no detectable terminal artifact |
| `ox_like_z` | `float` | Z-score for the terminal G→T excess relative to its standard error; strongly negative values indicate no terminal artifact (the terminal rate is below, not above, the interior) |
| `ox_like_ci_lo` | `float` | Lower bound of the 95 % confidence interval on the excess fraction |
| `ox_like_ci_hi` | `float` | Upper bound of the 95 % confidence interval on the excess fraction |
| `ox_uniformity_ratio` | `float` | Ratio of interior θ to terminal θ; values near 1.0 indicate uniform G→T across the read, consistent with bulk chemical damage; values substantially above 1.0 indicate interior enrichment; values substantially below 1.0 flag terminal artifact |
| `control_mode` | `string` | Method used to compute the composition baseline: `"chargaff"` (Chargaff second-parity correction on the complementary strand) or `"trinuc"` (trinucleotide-context matched control) |
| `fit_degenerate` | `boolean` | `true` when the likelihood surface is flat or the estimator could not separate terminal from interior signal; θ estimates should not be interpreted when this is `true` |
| `gc_skew_warning` | `boolean` | `true` when the library shows pronounced GC skew that may confound the composition control; results remain reported but should be interpreted with caution |

#### Interpretation

`ox_theta_interior` near 0.5 with a strongly negative `ox_like_z` (e.g., below −3) is the
expected pattern for genuine bulk-phase oxidative damage: G→T is distributed uniformly, so
there is no terminal excess to detect. An `ox_like_z` near zero or positive, combined with
`ox_theta` materially above `ox_theta_interior`, would indicate a terminal G→T spike
consistent with library-preparation artifact from 8-oxoG formed during end-repair. The
`ox_uniformity_ratio` provides a complementary view: ratios near 1.0 support uniform
damage, while values well above 1.0 indicate that interior positions are more oxidized
than terminals (which is the profile expected from preserved, pre-extraction oxidation
rather than a preparation artifact). Discard or flag libraries where `fit_degenerate` is
`true`, and cross-check `gc_skew_warning` outputs against the `control_mode` field to
assess whether a trinucleotide control would be more appropriate than a Chargaff one.


---