---
type: libtaph-json-field
title: interior_ct_cluster
tier: standard
estimand: interior C-to-T substitution cluster analysis
stability: stable
emitted_by: profile_json.cpp
---

### `interior_ct_cluster` block — interior co-deamination

Under prolonged exposure conditions, cytosines within single-stranded micro-domains
— short stretches that transiently unpair during storage — can be co-deaminated in the
same hydrolytic event. This produces co-occurring T residues at nearby positions more
often than expected under site independence. The signal is in the read interior (middle
third), away from the terminal overhangs that drive Channel A.

**Chargaff's second parity rule and the AG control.** Chargaff's second (intra-strand)
parity rule states that within a single strand, [A] ≈ [T] and [C] ≈ [G]
[[Rudner et al. 1968](#references)]. This symmetry arises from the evolutionary history
of symmetric mutation processes acting on the genome. It underpins the AG co-occurrence
control used here: if elevated T co-occurrence at adjacent {C,T} sites arose purely from
strand composition rather than damage clustering, the G co-occurrence at adjacent {A,G}
sites should show a matched elevation. The contrast between the CT and AG tracks therefore
isolates genuine deamination clustering from composition effects.

For each pair of non-CpG `{C,T}` sites at separation d = 1–10 bp, the observed co-T
fraction is compared to the expectation under site independence.

| Field | Description |
|-------|-------------|
| `short_asym_log2oe` | log₂(CT obs/exp) − log₂(AG obs/exp); composition-corrected effect size |
| `short_log2oe` | log₂(CT obs/exp) without composition correction |
| `short_z` | Normalised CT-vs-AG contrast statistic, summed over d = 1–5 |
| `short_obs` | Observed co-T pairs at d = 1–5 |
| `short_exp` | Expected co-T pairs at d = 1–5 under independence |
| `reads_used` | Reads contributing ≥ 2 eligible non-CpG {C,T} sites |
| `sep_log2oe[10]` | log₂(obs/exp) at each separation d = 1–10 |

`short_z > 3` with positive `short_asym_log2oe` indicates genuine interior co-deamination
beyond what strand composition predicts. In most aDNA samples this signal is weak; it
becomes informative for very old material or samples with unusual interior damage such as
permafrost specimens.

---