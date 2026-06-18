---
type: libtaph-json-field
title: tetranuc_damage_rates
tier: diagnostics
estimand: per-tetranucleotide context damage rates (900+ sub-fields)
stability: experimental
emitted_by: profile_json.cpp
---

### `tetranuc_damage_rates` block — 4-mer context C→T and G→T rates

The tetranucleotide (4-mer) context block extends the trinucleotide spectrum to include
one flanking base on each side of the central C or G, providing finer resolution of
sequence-context effects. Four keys are emitted:

| Key | Central base | Damage type | Use case |
|-----|-------------|-------------|----------|
| `ct_5prime` | C | C→T terminal excess | CHG vs CHH vs CpG methylation context separation |
| `gt_5prime` | G | G→T terminal excess | 8-oxoG context specificity; GG-adjacent enrichment |
| `ct_5prime_by_deam` | C | C→T per deamination stratum | Ancient-component C→T rates without EM |
| `gt_5prime_by_deam` | G | G→T per deamination stratum | Ancient-component G→T; stratum-monotonicity test |

**Context encoding.** Keys are 4-character strings `PGNN` where P = prev base, G/C =
central reference base, N1 = next1, N2 = next2 (A=0, C=1, G=2, T=3 in index space).
Example: `"AGCA"` encodes 5′-A–**G**–C–A-3′ with G as the reference base for G→T.

**Per-context fields** (each context key → object):

| Field | Description |
|-------|-------------|
| `terminal_rate` | Fraction of central-base reads showing the substitution at 5′ terminal positions (1–4) |
| `interior_rate` | Same fraction at interior positions (10–14) |
| `excess` | `terminal_rate − interior_rate`; positive = terminal-enriched damage |
| `terminal_n` | Total reads (reference + damaged) at terminal positions for this context |
| `interior_n` | Total reads at interior positions |

**`gt_5prime_by_deam` and stratum monotonicity.** The per-stratum variant is an array of
5 objects (strata 0–4), each with the same per-context structure. Stratum 0 = reads with
no terminal C→T (modern-enriched); stratum 4 = reads with terminal C→T at both ends
(ancient-enriched). For a genuine oxidative signal, the pool-weighted mean G→T excess
should rise monotonically from stratum 1 to stratum 4 because oxidised molecules are
preferentially retained in the ancient fraction. For an OxoG library-preparation artefact
(shearing-induced 8-oxoG), the G→T excess is deposited during library prep and is
uncorrelated with which reads carry terminal C→T — the stratum profile should be flat.

**GG-adjacent enrichment.** Genuine 8-oxoG formation by one-electron guanine oxidation
is strongly favoured at 5′-GG-3′ dinucleotides and GG-containing tracts because
oxidative hole-transport along the DNA helix delocalises to the lowest-potential
guanine (Sugiyama & Saito 1996). The `gt_5prime` context map should therefore show
elevated excess for 4-mer keys where P = G (5′-GG context) or N1 = G (GG 3′-adjacent),
relative to non-GG contexts. A uniform elevation across all 64 G-centred contexts would
be more consistent with a diffuse oxidant (e.g. free hydroxyl radical) or an
artefactual source.

---