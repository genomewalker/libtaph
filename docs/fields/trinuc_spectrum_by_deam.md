---
type: libtaph-json-field
title: trinuc_spectrum_by_deam
tier: diagnostics
estimand: trinucleotide spectra stratified by deamination (large nested arrays)
stability: experimental
emitted_by: profile_json.cpp
---

### `trinuc_spectrum_by_deam` block — deamination-stratified trinucleotide context counts

The raw trinucleotide context counts underlying `deam_context_spectrum`, stratified by
deamination stratum. This block provides the data necessary to compute per-stratum
C→T excess rates for any context and is the primary input for
stratum-aware ancient-fraction decompositions.

**Encoding.** Contexts are indexed as `prev × 16 + mid × 4 + next` with A=0, C=1, G=2,
T=3 (64 total contexts). Terminal positions are 5′ read positions 1–4; interior
positions are read positions 10–14. The same convention is used at the 3′ end.

| Field | Type | Description |
|-------|------|-------------|
| `n_strata` | int | Number of deamination strata (typically 5: strata 0–4) |
| `deam_stratum_reads` | int[n_strata] | Read count per stratum |
| `tri_5prime_terminal` | int[n_strata][64] | Per-stratum counts of each trinucleotide context at 5′ terminal positions (1–4) |
| `tri_5prime_interior` | int[n_strata][64] | Per-stratum counts at 5′ interior positions (10–14) |
| `tri_3prime_terminal` | int[n_strata][64] | Per-stratum counts at 3′ terminal positions |
| `tri_3prime_interior` | int[n_strata][64] | Per-stratum counts at 3′ interior positions |

To compute a per-stratum C→T excess for context N[C]N (e.g. ACT, index = 0×16+1×4+3 = 7):
`excess_s = terminal[s][idx_C→T] / (terminal[s][idx_C→T] + terminal[s][idx_C→C]) − interior equivalent`
where `idx_C→T` encodes the T-containing trinucleotide and `idx_C→C` the reference
C-containing trinucleotide for the same flanking bases.

Note that this block contains raw counts, not rates. Divide by the sum of the reference
and damaged context counts to obtain rates.

---