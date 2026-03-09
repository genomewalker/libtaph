# Damage types and channels

Ancient DNA carries a record of post-mortem chemistry in its substitution patterns.
libdart-damage tracks five independent biochemical damage channels, each measuring a
distinct degradation process. The channels cross-validate each other, distinguish genuine
ancient damage from modern library-prep artifacts, and together feed the library-type
classifier used by fqdup for damage-aware deduplication.

---

## The five biochemical channels

### Channel A — Cytosine deamination (C→T / G→A)

**Chemistry:** Cytosine spontaneously deaminates to uracil, which the polymerase reads as
thymine. The rate is fastest at single-stranded regions near the ends of fragmented
molecules and decays exponentially toward the interior.

**What is measured:**

| Sub-channel | Signal | Where |
|---|---|---|
| ct5 | T/(T+C) excess above baseline | 5' end, positions 0–14 |
| ga3 | A/(A+G) excess above baseline, smooth decay | 3' end, positions 1–14 |
| ga0 | A/(A+G) excess at position 0 only, spike | 3' position 0 |
| ct3 | T/(T+C) excess above baseline | 3' end, positions 0–14 |

The decay model at each end is:

$$\delta(p) = b + A \cdot e^{-\lambda p}$$

where $b$ is the interior background rate, $A$ is the damage amplitude, $\lambda$ is the
decay constant, and $p$ is the distance from the terminus (0-indexed).

The calibrated damage rate at position 0 is:

$$D_{\max} = \frac{A}{1 - b}$$

**Diagnostic role:** Primary damage quantification channel. Used by fqdup for position
masking. The ct3 sub-channel is a negative control for DS libraries (should be flat) and a
damage signal for SS-original libraries. The ga0 spike at position 0 is a distinct signal
from the smooth ga3 decay (see library classification below).

---

### Channel B — Stop codon conversion (composition-independent C→T validation)

**Chemistry:** The same C→T deamination that creates the ct5 signal also converts specific
sense codons to stop codons. CAA (Gln), CAG (Gln), and CGA (Arg) each become a stop codon
(TAA, TAG, TGA) when the first C is deaminated to T. Because this test uses triplet context
rather than simple base frequency, it is insensitive to GC-composition bias that can
confound Channel A.

**What is measured:** At each terminal position $p$, the stop codon conversion rate:

$$r_B(p) = \frac{\text{TAA}+\text{TAG}+\text{TGA}}{\text{CAA}+\text{CAG}+\text{CGA}+\text{TAA}+\text{TAG}+\text{TGA}}\bigg|_p$$

The interior baseline provides a reference, and WLS fits the same exponential decay model
as Channel A to give a structural $D_{\max}$ estimate from codon transitions alone.

A 3′ variant (**Channel B₃′**) tracks G→A stop codon conversions (TGG→TAG, TGG→TGA) to
validate the 3′ signal in SS libraries.

**Diagnostic role:** Cross-validation of Channel A.

- If Channel A and Channel B both detect signal → `damage_validated = true`.
- If Channel A detects signal but Channel B contradicts it →
  `damage_artifact = true` (likely composition bias or modern contamination passing
  the terminal frequency test).
- If Channel B data are insufficient (low coverage or codon-poor reads) →
  `channel_b_valid = false`, soft suppression only.

---

### Channel C — Oxidative stop codon tracking (G→T uniformity test)

**Chemistry:** 8-Oxoguanine (8-oxoG) forms when guanine is oxidized. The polymerase
misreads 8-oxoG as adenine, producing G→T transversions. Unlike cytosine deamination,
oxidation occurs throughout the molecule, not preferentially at the ends — damage
accumulates wherever guanine is exposed during burial.

**What is measured:** Stop codon conversion via G→T at glutamate and glycine codons
(GAG/GAA/GGA → TAG/TAA/TGA), measured at terminal vs. interior positions.

The key diagnostic statistic is the uniformity ratio:

$$U_C = \frac{\text{stop rate (terminal)}}{\text{stop rate (interior)}}$$

$U_C \approx 1$ means the G→T conversion is evenly distributed along the read — the
signature of genuine oxidative damage. $U_C \gg 1$ means the signal is terminal-enriched,
suggesting co-occurrence with deamination rather than independent oxidation.

**Diagnostic role:** Detects 8-oxoG oxidative damage independently of deamination. Flags
`channel_c_valid` when coverage is sufficient. High terminal enrichment is treated as an
artifact flag, not as additional deamination evidence.

---

### Channel D — Direct G→T / C→A transversion tracking (8-oxoG rate)

**Chemistry:** Same 8-oxoG chemistry as Channel C, but measured directly as a raw
transversion frequency rather than through stop codon context. G→T and its complement C→A
are tracked at terminal positions and compared to the interior baseline.

**What is measured:**

- `ox_gt_rate_terminal` / `ox_gt_rate_interior` — G→T rate at terminals vs. interior.
- `ox_ca_rate_terminal` / `ox_ca_rate_interior` — C→A rate at terminals vs. interior.
- `ox_gt_uniformity` — terminal/interior ratio; ≈ 1 for genuine oxidation.
- `ox_gt_asymmetry` — correlation of G→T with C→A; high asymmetry confirms 8-oxoG
  because C→A is the complementary transversion expected on the opposite strand.

**Diagnostic role:** Direct quantification of oxidative damage rate. `ox_damage_detected`
when both G→T and C→A channels agree. `ox_is_artifact` when the pattern is
terminal-enriched in a way inconsistent with random oxidation (suggests oxidation occurred
during sample preparation rather than in situ).

---

### Channel E — Depurination (AP-site enrichment at strand breaks)

**Chemistry:** Purines (adenine and guanine) are more susceptible than pyrimidines to
hydrolytic base loss (depurination) under acidic and warm conditions. The resulting
apurinic (AP) site is a strand-break precursor; fragmentation preferentially occurs at AP
sites, leaving purines enriched at the newly formed 5′ ends. This produces a purine
enrichment at the 5′ terminus that is distinct from both deamination and oxidation.

**What is measured:**

- `purine_rate_terminal_5prime` — A+G fraction at 5′ terminal positions.
- `purine_rate_interior` — A+G fraction at interior positions (baseline).
- `purine_enrichment_5prime` / `purine_enrichment_3prime` — excess above interior baseline.
- `depurination_detected` — flag when 5′ enrichment is statistically significant.

**Diagnostic role:** Independent evidence that fragmentation occurred at AP sites (genuine
ancient damage), rather than random mechanical shearing. High depurination enrichment
confirms the ancient origin of fragmentation even when deamination is low. Not used
directly by fqdup for masking, but reported for sample characterisation.

---

## Channel summary

| Channel | Measures | Chemistry | Read position | Used for |
|---------|----------|-----------|---------------|----------|
| A (ct5/ga3/ga0/ct3) | C→T and G→A rate | Cytosine deamination | Terminal, exponential decay | Primary damage quantification; position masking in fqdup |
| B / B₃′ | Stop codon C→T / G→A | Deamination, triplet-context | Terminal, exponential decay | Cross-validation of Channel A; artifact detection |
| C | G→T stop codon rate | 8-oxoG oxidation | Uniform across read | Uniformity test; oxidative damage detection |
| D | G→T and C→A rate | 8-oxoG oxidation | Uniform across read | Oxidative damage quantification; artifact flag |
| E | Purine enrichment at 5′ | Depurination / AP-site fragmentation | 5′ terminal | Ancient origin confirmation; sample characterisation |

---

## Library-type classification from Channel A sub-channels

The four Channel A sub-channels (ct5, ga3, ga0, ct3) together determine which of six
biological damage patterns best describes the library. Which sub-channels are active
depends on how the library was prepared.

### DS — double-stranded (symmetric deamination)

**Library prep:** Blunt-end ligation (NEBNext, TruSeq). Both strands of each fragment are
adapter-ligated.

**Pattern:** ct5 ≈ ga3, both smooth exponential decay. The complementary strand read 5'→3'
maps the same C→T deamination onto the 3' end as G→A. Symmetric and equal amplitude.

**Active sub-channels:** ct5 (smooth) + ga3 (smooth).

---

### DS + end-repair artifact

**Library prep:** Same as DS. A prep-associated G→A spike appears at 3' position 0 from
the ligation junction, superimposed on the genuine DS decay.

**Pattern:** ct5 + ga3 (smooth) + ga0 spike at position 0 only.

**Active sub-channels:** ct5 (smooth) + ga3 (smooth) + ga0 (spike).

---

### SS complement-orientation

**Library prep:** Single-stranded protocols (Gansauge & Meyer 2013). In complement
orientation, the 3' ligation junction carries a position-0 spike without a corresponding
smooth decay. The DS mirror logic does not apply because only one strand is processed.

**Pattern:** ga0 spike only; ct5 and ga3 flat.

**Active sub-channels:** ga0 (spike only).

---

### SS original-orientation

**Pattern:** Both the 5' end (ct5) and 3' end (ct3) show C→T excess, because the
same strand is damaged at both termini. No ga3 smooth decay; no ga0 spike.

**Active sub-channels:** ct5 (smooth) + ct3 (elevated).

---

### SS mixed orientations

**Pattern:** Original-orientation reads contribute ct5; complement-orientation reads
contribute ga0. Together: ct5 + ga0 without smooth ga3.

**Active sub-channels:** ct5 (smooth) + ga0 (spike).

---

### UNKNOWN

No Channel A sub-channel produces a signal above the null model. The library has no
detectable deamination, either because it is modern, because damage has not accumulated,
or because coverage is insufficient for any channel to reach significance.

---

## How fqdup uses damage channels

fqdup deduplicates reads by hashing sequences. PCR copies of the same original molecule
that happen to carry deamination at different positions will produce different sequences
and different hashes, inflating the unique read count. Damage-aware hashing corrects this
by masking the affected terminal positions before hashing.

### Which channels drive masking

Only **Channel A** drives position masking. Channels B–E are diagnostic: they validate
that Channel A is genuine ancient damage, flag artifacts, and characterise the sample, but
they do not directly change which positions are masked.

### Revcomp-equivariant masking

fqdup uses a canonical hash — `min(hash(seq), hash(rc(seq)))` — so that a molecule and
its reverse complement always land in the same cluster. The mask must therefore satisfy:

$$\text{mask}(\text{rc}(\text{seq})) = \text{rc}(\text{mask}(\text{seq}))$$

fqdup enforces this by applying a symmetric rule: if position $p$ from the 5′ end is
masked, position $p$ from the 3′ end is also masked. This guarantees the invariant
regardless of which strand is sequenced.

### Masking threshold

A position is masked when the empirically observed background-corrected excess exceeds
a threshold $\tau$ (default 5%):

$$\text{mask position } p \iff \delta(p) - b > \tau$$

The number of masked positions is an emergent property of the data. A slowly decaying
library (low $\lambda$) will mask more positions than a rapidly decaying one.

### Library type and which sub-channels are used

| Library type | Positions masked | Notes |
|---|---|---|
| DS | 5′ pos $p$ and 3′ pos $p$ symmetrically from ct5/ga3 | Genuine bilateral deamination |
| DS + artifact | Same as DS; ga0 spike falls within the masked terminal region | Artifact already covered |
| SS complement | 3′ pos 0 (ga0) only | Single terminal position; paired under revcomp |
| SS original | 5′ and 3′ symmetrically from ct5/ct3 | Both ends carry deamination |
| SS mixed | 5′ from ct5 + 3′ pos 0 from ga0; paired under revcomp | Asymmetric channels, symmetric mask |
| UNKNOWN | None | No trusted channel; standard hashing used |

### Artifact detection

If `damage_artifact = true` (Channel A fires but Channel B contradicts it), fqdup treats
the damage as unvalidated and optionally suppresses masking, depending on user settings.
This prevents modern contamination with GC-composition bias from being incorrectly
masked and over-collapsed.
