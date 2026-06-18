---
type: libtaph-json-field
title: fragmentation
tier: standard
estimand: read-length distribution and damage-length coupling (fragmentation proxy)
stability: stable
emitted_by: profile_json.cpp
---

### `fragmentation` block — read-length and damage-length coupling proxy

Reference-free FASTQ can observe fragment-size structure, but it cannot by itself assign a
unique strand-break mechanism. The `fragmentation` block therefore reports read-length
distribution summaries and the already fitted damage-versus-length coupling as an
empirical proxy. It is explicitly separate from `bulk_damage.lambda`, which is a positional
terminal-damage decay parameter.

| Field | Description |
|-------|-------------|
| `valid` | `true` when at least 100 reads contributed to the read-length histogram |
| `observable` | `read_length_distribution_and_damage_length_coupling` |
| `mechanism_status` | `empirical_proxy` |
| `reference_free_identifiability` | Notes that this is a fragmentation/selection proxy, not a causal strand-break rate |
| `not_equivalent_to` | Always `bulk_damage.lambda` |
| `mean_length`, `median_length`, `n50_length` | Reference-free length-distribution summaries |
| `short_fraction_lt_50bp`, `short_fraction_lt_70bp` | Fraction of reads in short-fragment bins |
| `topbin_fraction_ge_225bp` | Fraction in the overflow long-fragment bin |
| `damage_length_coupling_slope`, `damage_length_coupling_weight` | Coupling between terminal-damage excess and read length from the bulk model |
| `length_histogram` | Fixed-width read-length histogram used for the summaries |

Interpret this block alongside extraction method, size selection, trimming, library type,
and deamination evidence. It is the appropriate place for fragmentation-like analyses;
do not use `bulk_damage.lambda` as a fragment-length or strand-break rate.

---