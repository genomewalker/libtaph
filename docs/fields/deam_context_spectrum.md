---
type: libtaph-json-field
title: deam_context_spectrum
tier: standard
estimand: 16-channel trinucleotide N[C→T]N excess at 5-prime and 3-prime termini
stability: stable
emitted_by: profile_json.cpp
---

### `deam_context_spectrum` block — trinucleotide C→T context spectrum

Deamination rate depends on the local sequence context because stacking interactions
and base-pairing geometry modulate cytosine accessibility and hydration. The
`deam_context_spectrum` block reports the **16-channel trinucleotide C→T terminal
excess spectrum** — the 5′ terminal C→T rate relative to the interior rate for each
of the 16 N[C→T]N contexts (N ∈ {A,C,G,T}).

Channel index ordering: `ACA ACCi ACG ACT  CCA CCC CCG CCT  GCA GCC GCG GCT  TCA TCC TCG TCT`
(i.e. `prev ∈ ACGT` outer, `next ∈ ACGT` inner). The four CpG channels (next = G,
indices 2, 6, 10, 14: ACG, CCG, GCG, TCG) are driven primarily by 5-methylcytosine
deamination and are mechanistically distinct from the twelve non-CpG channels.

| Field | Type | Description |
|-------|------|-------------|
| `channels` | string[16] | Channel labels in the order used by all 16-element arrays |
| `ct_5prime_excess` | float[16] | 5′ terminal C→T rate minus interior C→T rate, per context |
| `ct_3prime_excess` | float[16] | 3′ terminal C→T rate minus interior (DS complement arm) |
| `ga_3prime_rc_excess` | float[16] | 3′ G→A reverse-complement excess; independent check on the same damage events read from the complementary strand |
| `comparison_vector` | float[16] | Unsigned element-wise maximum of 5′ and 3′ excess — summarises the dominant arm |
| `signed_comparison_vector` | float[16] | Signed version; sign indicates which arm dominates |
| `primary_arm` | string | `"5prime"` or `"3prime"` depending on which terminus carries stronger excess |
| `sum_positive_excess` | float | Sum of positive entries in `ct_5prime_excess` — a scalar summary of total 5′ context-dependent deamination |

CpG channels (ACG, CCG, GCG, TCG) are informative for methylation state: elevated
excess in CpG contexts relative to the matched non-CpG contexts (ACN, CCN, GCN, TCN
for N ≠ G) indicates a methylated source genome. Non-CpG channels reflect unmethylated
hydrolytic deamination and are sensitive to the chemical environment at burial (pH,
organic-nitrogen load, water activity).

---