---
type: libtaph-json-field
title: oxidation_like
tier: standard
estimand: boolean flag: damage pattern matches oxidative profile
stability: stable
emitted_by: profile_json.cpp
---

### `oxidation_like`

This block quantifies evidence of oxidative damage patterns in the library, primarily
G→T transversions arising from 8-oxo-7,8-dihydroguanine (8-oxoG) lesions. 8-oxoG is
formed when reactive oxygen species attack the C8 position of guanine; the modified base
adopts a *syn* conformation that pairs preferentially with adenine rather than cytosine
during replication or end-repair, causing G→T (and complementarily C→A) substitution
signatures. In ancient or degraded samples, oxidative damage accumulates post-mortem
alongside hydrolytic deamination, but the two processes leave chemically and
positionally distinct marks: deamination is enriched at single-stranded overhangs and
produces C→T at read termini, whereas 8-oxoG-driven G→T substitutions are distributed
throughout read bodies and are exacerbated by oxidising preservation environments. The
block is absent (`null` or not present) when the library type or available read depth
does not support a reliable oxidation-like estimate.

| Field | Type | Description |
|-------|------|-------------|
| — | — | Object structure depends on library type — see JSON output. |

The `oxidation_like` block should be interpreted alongside the `deamination` block:
elevated oxidation-like signal with low deamination enrichment suggests chemical
degradation that does not originate from the biological damage expected in genuinely
ancient sequences, and may reflect sample storage, extraction chemistry, or
environmental redox conditions rather than age. Conversely, co-elevation of both signals
is consistent with authentic ancient degradation under oxidising taphonomic conditions.
A `null` or absent block does not indicate the absence of oxidative damage — it indicates
insufficient data or an incompatible library configuration to resolve the signal
reliably; re-evaluation with deeper sequencing or a compatible library type is advisable
before drawing negative conclusions.


---