---
type: libtaph-json-field
title: bic
tier: diagnostics
estimand: nested model comparison BIC metrics (75 sub-fields)
stability: experimental
emitted_by: profile_json.cpp
---

### `bic` block — library-type classifier

Library preparation determines which damage channels are active, so the BIC model
selection that drives classification also reveals which sub-channels contributed signal.
This block exposes the classifier's evidence for inspecting borderline calls.

Double-stranded (DS) library protocols (e.g., NEBNext, TruSeq [[Meyer & Kircher
2010](#references)]) ligate adapters to both strands. Reads from the original strand show
C→T at the 5′ end; reads from the complementary strand map the same damaged cytosines as
G→A at the 3′ end. Single-stranded (SS) protocols (e.g., Gansauge & Meyer [[2013](#references)];
Santa Cruz [[Kapp et al. 2021](#references)]; SRSLY [[Gansauge et al. 2017](#references)]) circularize and sequence only one strand per
molecule, producing library-type-specific patterns depending on which strand is captured
(see [Library-type classification](#library-type-classification)).

The classifier fits four sub-channel amplitudes simultaneously:

- **ct5** — 5′ C→T exponential decay (deamination on the sequenced strand)
- **ga3** — 3′ G→A exponential decay (deamination on the complementary strand, appearing
  as G→A when read 5′→3′)
- **ga0** — 3′ position-0 G→A spike (ligation-site pattern in SS complement-orientation
  libraries; also produced by certain DS end-repair protocols)
- **ct3** — 3′ C→T decay (SS original-orientation only, where the same strand carries
  damage at both ends)

These are evaluated under seven biological models: a null (all channels flat), two DS
models (symmetric and with end-repair spike), a DS end-repair-only model, and three SS
models (complement, original, and mixed orientations). The model with the lowest BIC wins
[[Schwarz 1978](#references)].

| Field | Description |
|-------|-------------|
| `bias` | BIC of the null model (all channels flat; no damage) |
| `ds` | BIC of the best-fitting DS model |
| `ss` | BIC of the best-fitting SS model |
| `ct5_amp` | Fitted ct5 amplitude |
| `ga3_amp` | Fitted ga3 amplitude |
| `ga0_amp` | Fitted ga0 amplitude (3′ position-0 spike) |
| `ct3_amp` | Fitted ct3 amplitude |

Lower BIC = better fit. `library_type` is `"double-stranded"` when `ds < ss`,
`"single-stranded"` when `ss < ds`, and `"unknown"` when neither beats `bias`.

---