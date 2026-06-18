---
type: libtaph-json-field
title: complement_asymmetry
tier: standard
estimand: G-to-T vs C-to-A strand-balance asymmetry damage proxies
stability: stable
emitted_by: profile_json.cpp
---

### `complement_asymmetry` block — oxidative damage (Channels C, D, F, G, and H)

Guanine is the most easily oxidised DNA base owing to its low one-electron reduction
potential [[Steenken & Jovanovic 1997](#references)]. Reactive oxygen species — primarily hydroxyl radical (•OH) and singlet oxygen (¹O₂) —
attack the C8 position of guanine to form 7,8-dihydro-8-oxoguanine (8-oxoG) [[Cadet et al.
2010](#references)]. The replicative polymerase misreads 8-oxoG as adenine via a
*syn* conformation, inserting dAMP opposite the lesion and producing G→T transversions
in the forward strand and C→A in the reverse strand [[Shibutani et al. 1991](#references)].
Unlike cytosine deamination, which is confined to single-stranded overhangs, oxidative
damage accumulates throughout the molecule with no strong positional preference.

**Chargaff's first rule and the strand asymmetry principle.** Chargaff's first rule
states that in a double-stranded molecule [G] = [C] and [A] = [T] across complementary
strands [[Chargaff 1950](#references)]. As a consequence, if oxidation affects both
strands symmetrically, the G→T rate in forward-strand reads exactly equals the C→A rate
in reverse-strand reads, and they cancel in a pooled DS library; `s_gt` and `D` are
near zero even when oxidation is substantial. These statistics become informative when
oxidation is strand-asymmetric (one strand preferentially oxidised during burial) or when
only one strand is present (SS libraries, where no complementary strand exists to provide
the cancelling C→A signal) [[Mitchell & Bridge 2006](#references)].

#### Rates and strand asymmetry

| Field | Description |
|-------|-------------|
| `tg_interior` | G→T rate in the read interior (background) |
| `ac_interior` | C→A rate in the read interior |
| `tg_terminal` | G→T rate at 5′ terminal positions |
| `ac_terminal` | C→A rate at 5′ terminal positions |
| `gt_bg_fitted` | Fitted uniform background from the G→T model |
| `gt_term_fitted` | Fitted terminal amplitude |
| `gt_decay_fitted` | Fitted decay constant |
| `s_gt` | G→T vs C→A **strand-asymmetric oxidation contrast** (`s_gt_estimand=strand_asymmetric_oxidation_contrast_not_total_oxidation`), not total oxidation; near zero for balanced DS oxidation, informative for SS |
| `s_gt_valid` | `true` only when interior composition is Chargaff-balanced (`chargaff_gc_balance ≤ 0.02`) **and** the C/D coverage gate (`d_computed`) passes. A gate flag — `s_gt`/`D` values are never nulled or altered |
| `chargaff_gc_balance` | `|G_interior − C_interior| / (G_interior + C_interior)` over interior base composition. When large, an `s_gt` contrast may be composition-driven rather than oxidative |
| `D` | Overall 8-oxoG asymmetry index |
| `per_pos_5prime_gt[15]` | Raw G→T fraction at 5′ positions 0–14 |

`ox_gt_uniformity` (terminal / interior ratio) near 1.0 indicates G→T stop-codon
conversions distributed uniformly across the read, consistent with 8-oxoG oxidation
rather than terminal deamination. `ox_gt_asymmetry` encodes the strand contrast
(`ox_gt_baseline − ox_ca_baseline`). The 16-context panel is in `s_oxog_16ctx`.

#### 8-oxoG 16-context panel

G→T asymmetry split by flanking trinucleotide context (N**G**N, all 16 combinations).
The context index is `4 × enc(left) + enc(right)` where `enc(A)=0, C=1, G=2, T=3`.

| Field | Description |
|-------|-------------|
| `s_oxog` | Overall G→T strand asymmetry at interior positions |
| `se_s_oxog` | Standard error of `s_oxog` |
| `s_oxog_16ctx[16]` | Per-context G→T asymmetry; `null` when coverage < 500 |
| `cov_oxog_16ctx[16]` | Coverage per trinucleotide context bin |

In modern oxidative chemistry, guanine oxidation in GG-containing contexts is
preferentially favoured [[Cadet et al. 2010](#references)]. A single-context spike
confined to one or two bins is more consistent with a preparation artefact than with
the broad multi-context enrichment expected from genuine oxidative damage. A spike confined to one or two contexts is more consistent with an
oxidation artefact introduced during library preparation.

---