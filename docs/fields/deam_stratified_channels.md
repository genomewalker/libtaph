---
type: libtaph-json-field
title: deam_stratified_channels
tier: standard
estimand: per-deamination-stratum methylation and depurination channel summaries
stability: stable
emitted_by: profile_json.cpp
---

### `deam_stratified_channels` block — damage-stratum-specific methylation and depurination

Reads are binned by the number of terminal C→T events they carry at the 5′ and 3′
positions (the **deamination stratum**). Stratum 0 contains reads with no terminal
C→T at either end (enriched for undamaged or modern molecules); stratum 4 (or the
highest index) contains reads with the maximum counted terminal C→T events at both
ends (enriched for the most-ancient fraction). Three aggregate strata are reported:

| Stratum | Selection | Typical source |
|---------|-----------|----------------|
| `modern` | Reads in stratum 0 (no terminal C→T) | Undamaged / modern contamination |
| `bulk` | All reads | Full library mixture |
| `ancient` | Reads in the highest deamination stratum | Ancient-enriched fraction |

For each stratum, the block reports:

| Field | Description |
|-------|-------------|
| `n_reads` | Number of reads in this stratum |
| `methylation_next_cond_logodds` | Log-odds of the next base being G given a terminal C→T event, relative to a non-CpG baseline; positive values indicate CpG-context enrichment (methylation signal). Computed only over reads with ≥ 1 terminal C→T |
| `depurination_index` | Terminal A+G enrichment (per-read AP-site proxy, Channel E) computed restricted to this stratum |

Comparing `modern` and `ancient` strata isolates deposition-time chemistry from
modern contamination or taphonomic remodelling: if `methylation_next_cond_logodds` is
substantially higher in the `ancient` stratum than in `modern`, the ancient fraction
retains methylation-enhanced CpG deamination consistent with the burial environment.

---