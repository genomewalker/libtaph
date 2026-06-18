---
type: libtaph-json-field
title: oxidative_scission
tier: standard
estimand: oxidative strand scission from C→G and A→T terminal channel enrichment
stability: stable
emitted_by: profile_json.cpp
---

### `oxidative_scission` block — oxidative strand breakage via 8-oxoguanine (Channel G)

Reactive oxygen species attack deoxyguanosine at the C8 position, generating
8-oxo-7,8-dihydroguanine (8-oxoG), the most abundant oxidative base lesion in DNA. In
ancient or degraded material, 8-oxoG accumulates in proportion to cumulative oxidative
stress during deposition, burial, or storage. Unlike deamination, which is driven
primarily by hydrolysis at single-stranded termini, oxidative damage proceeds through
radical chemistry and is not strongly strand-position-biased; however, strand breakage
preferentially occurs at or adjacent to oxidised sites, creating a detectable enrichment
of guanine (and its oxidative mismatch partner thymine) at newly exposed read termini.

The `oxidative_scission` block quantifies this terminal G→T imbalance as evidence of
oxidative fragmentation. The signal is evaluated symmetrically on both 5′ and 3′ ends for
double-stranded libraries and on the appropriate strand for single-stranded library
chemistries. The block is absent (or `null`) when the sequencing data are insufficient to
estimate oxidative rates reliably, or when the library configuration does not expose an
interpretable oxidative channel.

| Field | Type | Description |
|-------|------|-------------|
| Object structure depends on library type — see JSON output. | | |

Interpretation: a `detected: true` result indicates statistically significant terminal
G→T enrichment consistent with oxidative fragmentation; the magnitude of terminal rates
relative to interior background quantifies the excess signal. Low or zero enrichment does
not rule out 8-oxoG accumulation — oxidative lesions that do not cause strand breakage
are invisible to this channel. Co-occurrence of `detected: true` in both the oxidative
and deamination channels supports a genuinely ancient signal, while an isolated oxidative
signal with absent deamination warrants caution and may reflect chemical contamination or
oxidative preservation artefacts rather than age.

---