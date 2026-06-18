---
type: libtaph-json-field
title: end_motif_enrichment
tier: standard
estimand: terminal nucleotide over-representation from fragmentation sequence preference
stability: stable
emitted_by: profile_json.cpp
---

### `end_motif_enrichment`

This block characterises the per-base composition at the outermost four positions of the 5′ and 3′ read termini, expressed as excess frequency relative to interior read positions. It is emitted for all library types whenever sufficient read coverage is available; a missing or `null` block indicates that coverage was inadequate to compute the profile. The block provides a sequence-context complement to scalar damage metrics, resolving which bases are selectively enriched or depleted at read ends independently of any damage model.

Fragmentation of ancient nucleic acids is non-random. Enzymatic cleavage during library preparation, AP-site β-elimination, and hydrolytic strand scission each leave characteristic sequence preferences at newly exposed termini. Because these preferences are imprinted at the instant of breakage, the terminal base composition encodes information about the biochemical environment at the time of fragmentation. The per-position enrichment vectors capture this signal as a zero-centred excess: a value near +0.3 for cytosine at the 5′ terminal position means cytosine occurs ~30 percentage points more often at that position than it does in interior reads, while negative values indicate depletion. For double-stranded libraries, the 5′ and 3′ motifs are constrained by reverse-complement (RC) symmetry: a base enriched at 5′ position 1 must be matched by its Watson–Crick complement enriched at 3′ position 1. The `rc_symmetry_discordance_pos1` field quantifies the deviation from this constraint at position 1. For double-stranded libraries the discordance is typically near zero (observed median < 0.02); single-stranded libraries systematically violate RC symmetry at position 1 because the two ends derive from different portions of the original molecule rather than complementary strands, yielding discordance values commonly in the range 0.15–0.45.

| Field | Type | Description |
|-------|------|-------------|
| `5prime` | object | Per-base enrichment at positions 1–4 from the 5′ read terminus |
| `5prime.pos1` | object | Base excess frequencies at the first (outermost) 5′ position |
| `5prime.pos1.A` | float | Adenine excess at 5′ position 1 (terminal fraction minus interior fraction) |
| `5prime.pos1.C` | float | Cytosine excess at 5′ position 1 |
| `5prime.pos1.G` | float | Guanine excess at 5′ position 1 |
| `5prime.pos1.T` | float | Thymine excess at 5′ position 1 |
| `5prime.pos2` – `5prime.pos4` | object | Same structure for positions 2–4 from the 5′ end; values decay toward zero as positions move inward |
| `3prime` | object | Per-base enrichment at positions 1–4 from the 3′ read terminus; same structure as `5prime` |
| `3prime.pos1` – `3prime.pos4` | object | Base excess frequencies at positions 1–4 from the 3′ end |
| `rc_symmetry_discordance_pos1` | float | Deviation of the 3′ pos1 vector from the reverse-complement of the 5′ pos1 vector; near 0 for double-stranded libraries, substantially elevated for single-stranded libraries |

Strong cytosine enrichment at 5′ position 1 (pos1.C ≫ 0) accompanied by matching guanine enrichment at 3′ position 1 (pos1.G ≫ 0) is consistent with preferential cleavage 3′ of cytosine-containing di- or trinucleotide contexts, a pattern reported for AP-site and enzymatic fragmentation. Conversely, guanine enrichment at 5′ pos1 with cytosine enrichment at 3′ pos1 suggests cleavage 5′ of guanine. Values at positions 2–4 typically decay in magnitude and are most informative as a gradient: a sharp drop from pos1 to pos2 indicates a motif confined to a single terminal position, whereas extended enrichment across pos1–pos3 suggests a longer recognition context. The `rc_symmetry_discordance_pos1` field serves as an internal consistency check: a high value in a library nominally called double-stranded warrants re-examination of the library type assignment, while a low value in a single-stranded library may indicate hybrid or unexpected library chemistry.

---