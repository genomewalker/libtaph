---
type: libtaph-json-field
title: depurination
tier: standard
estimand: depurination (AP-site) damage metrics
stability: stable
emitted_by: profile_json.cpp
---

### `depurination` block — terminal purine enrichment / AP-site proxy (Channel E)

The glycosidic bond linking a purine base (adenine or guanine) to the deoxyribose sugar
is susceptible to acid-catalyzed hydrolysis. The rate under physiological conditions is approximately 3 × 10⁻¹¹ per purine per
second (~2,000–10,000 events per cell per day), substantially faster than cytosine
deamination [[Lindahl & Nyberg 1972](#references); [Lindahl 1993](#references)]. The
resulting apurinic (AP) site is a metastable abasic residue: β-elimination converts the
AP deoxyribose to a strand break, fragmenting the molecule at the site of base loss.
Because purines (A and G) are lost preferentially over pyrimidines (C and T), AP-site
breakage can enrich A+G at newly exposed read starts. Channel E measures this directly as
terminal `(A+G)/(A+C+G+T)` minus the middle-of-read purine fraction. The signal is an
empirical AP-site/depurination proxy independent of deamination; it supports but does not
prove a specific lesion without orthogonal controls [[Dabney et al. 2013](#references)].

| Field | Description |
|-------|-------------|
| `valid` | `true` when the 5′ terminal and interior all-base denominators are adequate |
| `detected` | `true` when valid 5′ purine enrichment is positive and statistically significant; `null` when invalid |
| `rate_terminal_5prime` | A+G fraction at the 5′ read start |
| `rate_terminal_3prime` | A+G fraction at the 3′ read end when evaluable |
| `enrichment_5prime` | Terminal 5′ A+G fraction minus interior A+G fraction |
| `enrichment_3prime` | Terminal 3′ A+G fraction minus interior A+G fraction |
| `rate_interior` | Interior A+G fraction (composition baseline) |
| `purine_z_5prime` | Exploratory two-proportion z-score for the 5′ enrichment |
| `ag_skew_ctrl_shift_5prime` | Deprecated A/(A+G) control shift; not the depurination statistic |

Channel E is reported for sample characterisation and is not used by fqdup for position
masking.

---