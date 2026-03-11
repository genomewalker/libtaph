# Damage types and channels

Ancient DNA carries a record of post-mortem chemistry in its substitution patterns.
libdart-damage tracks five core biochemical damage channels (A–E) plus eight supplementary
analyses: CpG-context split, interior CT clustering, oxoG 16-context panel, hexamer-based
detection, Briggs biophysical model, joint probabilistic model (BIC + Bayes factor),
mixture model (K-component EM), metaDMG-comparable estimate, adapter offset detection,
and two compositional controls (codon-position-aware damage, GC-conditional damage bins). The
channels cross-validate each other, distinguish genuine ancient damage from modern
library-prep artifacts, and together feed the library-type classifier used by fqdup for
damage-aware deduplication.

---

## The five biochemical channels

### Channel A: Cytosine deamination (C→T / G→A)

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

### Channel B: Stop codon conversion (composition-independent C→T validation)

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

### Channel C: Oxidative stop codon tracking (G→T uniformity test)

**Chemistry:** 8-Oxoguanine (8-oxoG) forms when guanine is oxidized. The polymerase
misreads 8-oxoG as adenine, producing G→T transversions. Unlike cytosine deamination,
oxidation occurs throughout the molecule, not preferentially at the ends, damage
accumulates wherever guanine is exposed during burial.

**What is measured:** Stop codon conversion via G→T at glutamate and glycine codons
(GAG/GAA/GGA → TAG/TAA/TGA), measured at terminal vs. interior positions.

The key diagnostic statistic is the uniformity ratio:

$$U_C = \frac{\text{stop rate (terminal)}}{\text{stop rate (interior)}}$$

$U_C \approx 1$ means the G→T conversion is evenly distributed along the read, the
signature of genuine oxidative damage. $U_C \gg 1$ means the signal is terminal-enriched,
suggesting co-occurrence with deamination rather than independent oxidation.

**Diagnostic role:** Detects 8-oxoG oxidative damage independently of deamination. Flags
`channel_c_valid` when coverage is sufficient. High terminal enrichment is treated as an
artifact flag, not as additional deamination evidence.

---

### Channel D: Direct G→T / C→A transversion tracking (8-oxoG rate)

**Chemistry:** Same 8-oxoG chemistry as Channel C, but measured directly as a raw
transversion frequency rather than through stop codon context. G→T and its complement C→A
are tracked at terminal positions and compared to the interior baseline.

**What is measured:**

- `ox_gt_rate_terminal` / `ox_gt_rate_interior`, G→T rate at terminals vs. interior.
- `ox_ca_rate_terminal` / `ox_ca_rate_interior`, C→A rate at terminals vs. interior.
- `ox_gt_uniformity`, terminal/interior ratio; ≈ 1 for genuine oxidation.
- `ox_gt_asymmetry`, correlation of G→T with C→A; high asymmetry confirms 8-oxoG
  because C→A is the complementary transversion expected on the opposite strand.

**Diagnostic role:** Direct quantification of oxidative damage rate. `ox_damage_detected`
when both G→T and C→A channels agree. `ox_is_artifact` when the pattern is
terminal-enriched in a way inconsistent with random oxidation (suggests oxidation occurred
during sample preparation rather than in situ).

---

### Channel E: Depurination (AP-site enrichment at strand breaks)

**Chemistry:** Purines (adenine and guanine) are more susceptible than pyrimidines to
hydrolytic base loss (depurination) under acidic and warm conditions. The resulting
apurinic (AP) site is a strand-break precursor; fragmentation preferentially occurs at AP
sites, leaving purines enriched at the newly formed 5′ ends. This produces a purine
enrichment at the 5′ terminus that is distinct from both deamination and oxidation.

**What is measured:**

- `purine_rate_terminal_5prime`, A+G fraction at 5′ terminal positions.
- `purine_rate_interior`, A+G fraction at interior positions (baseline).
- `purine_enrichment_5prime` / `purine_enrichment_3prime`, excess above interior baseline.
- `depurination_detected`, flag when 5′ enrichment is statistically significant.

**Diagnostic role:** Independent evidence that fragmentation occurred at AP sites (genuine
ancient damage), rather than random mechanical shearing. High depurination enrichment
confirms the ancient origin of fragmentation even when deamination is low. Not used
directly by fqdup for masking, but reported for sample characterisation.

---

## Hexamer-based damage detection

**Background:** Single-position statistics (e.g. T/(T+C) at position 0) can be noisy
for short reads or low coverage. Averaging over a hexamer window (positions 1–6) reduces
single-position artifacts and is less sensitive to the position-0 adapter bias.

**What is measured:** The first and interior hexamers of each read are tallied separately.
Terminal T/(T+C) and interior T/(T+C) are computed from the hexamer counts.

**Key outputs:**

| Field | Description |
|-------|-------------|
| `hexamer_terminal_tc` | T/(T+C) at terminal from hexamer analysis (pos 1–6) |
| `hexamer_interior_tc` | T/(T+C) at interior from hexamer analysis |
| `hexamer_excess_tc` | Terminal − interior excess (negative = inverted) |
| `hexamer_damage_llr` | Hexamer-based damage log-likelihood ratio |

**Interpretation:** `hexamer_excess_tc > 0` with a significant `hexamer_damage_llr`
corroborates Channel A independently of single-position artifacts.

---

## Briggs biophysical model

**Background:** Briggs et al. (2007) parameterised aDNA deamination as a two-state model:
each base is either in a single-stranded overhang (deamination rate δ_s, high) or in a
double-stranded interior (deamination rate δ_d, low). The observed position-dependent
rate is a mixture weighted by the probability of being in the overhang at that distance
from the terminus.

**What is measured:** The same exponential decay profile as Channel A is re-fitted under
the Briggs parameterisation to yield δ_s and δ_d at each end, plus R² goodness of fit.

**Key outputs:**

| Field | Description |
|-------|-------------|
| `delta_s_5prime` / `delta_s_3prime` | Single-stranded deamination rate (overhang) |
| `delta_d_5prime` / `delta_d_3prime` | Double-stranded background deamination rate |
| `r_squared_5prime` / `r_squared_3prime` | Goodness of fit for the Briggs model at each end |

**Interpretation:** High δ_s and low δ_d confirm classic aDNA terminal damage. R² < 0.5
indicates the exponential decay model fits poorly (possibly SS library, composition bias,
or very old heavily degraded material).

---

## Joint probabilistic model (BIC + Bayes factor)

**Background:** The per-channel statistics (Channel A LLR, Channel B LLR, etc.) are
combined into a single Bayesian model comparison: M₁ (damage present) vs M₀ (no damage).
The BIC difference ΔBIC = BIC₀ − BIC₁ quantifies evidence for damage relative to a null
model that expects flat terminal frequencies.

**What is measured:**

| Field | Description |
|-------|-------------|
| `joint_delta_max` | MLE damage rate under M₁ |
| `joint_lambda` | MLE decay constant under M₁ |
| `joint_delta_bic` | ΔBIC = BIC_M0 − BIC_M1 (positive = evidence for damage) |
| `joint_bayes_factor` | BF₁₀ ≈ exp(ΔBIC/2) |
| `joint_p_damage` | Posterior P(damage \| data) |
| `joint_model_valid` | True if sufficient read coverage for the joint model |

**Interpretation:** `joint_delta_bic > 10` is strong evidence for damage (BF > 150).
`joint_p_damage` can be used as a per-sample weight in downstream analyses.

---

## Mixture model (K-component EM over GC bins)

**Background:** A sequencing library may contain a mixture of ancient and modern DNA
molecules. Rather than a single d_max for the whole library, a K-component mixture
fitted over the GC-stratified bins separates high-damage components (ancient) from
low-damage components (modern or contamination).

**What is measured:** EM is run with BIC-guided component selection. Each component
has its own damage rate; reads are assigned soft membership probabilities.

**Key outputs:**

| Field | Description |
|-------|-------------|
| `mixture_K` | Number of components selected by BIC |
| `mixture_d_population` | E[δ] over all C-sites (whole-library average) |
| `mixture_d_ancient` | E[δ \| δ > 5%] — damage rate of the ancient tail |
| `mixture_pi_ancient` | Fraction of C-sites in high-damage components |
| `mixture_d_reference` | E[δ \| GC > 50%] — metaDMG proxy |
| `mixture_bic` | BIC for the selected K |
| `mixture_converged` | Whether EM converged |

**Interpretation:** `mixture_pi_ancient` close to 1 = pure ancient library.
`mixture_pi_ancient` near 0 with elevated `mixture_d_ancient` = small ancient fraction
contaminating a mostly modern library.

---

## metaDMG-comparable estimate (alignability-weighted d_max)

**Background:** Reference-based tools like metaDMG weight damage by read alignability:
highly unique reads (high alignability) contribute more to the d_max estimate than
repetitive reads that map to many locations. libdart-damage approximates this without
a reference by using per-read GC content and sequence complexity as an alignability
proxy, and blending the global and weighted estimates.

**Key outputs:**

| Field | Description |
|-------|-------------|
| `d_metamatch` | Calibrated metaDMG-comparable d_max estimate |
| `d_alignability_weighted` | Raw alignability-proxy-weighted d_max |
| `metamatch_gamma` | Blending coefficient (0 = global, 1 = weighted) |
| `mean_alignability` | Mean alignability proxy score across reads |
| `alignability_damage_corr` | Correlation between alignability and per-read damage |

**Interpretation:** `d_metamatch` is designed to be numerically comparable to metaDMG's
`D_max` output on the same library. It enables cross-tool comparisons without requiring
BAM alignment.

---

## Adapter offset detection

**Background:** Short adapter remnants (1–2 bp) at fragment termini displace the
biological damage signal by one or two positions. If position 0 appears depleted while
position 1 is enriched, the fit is restarted with the signal window shifted accordingly.

**Key outputs:**

| Field | Description |
|-------|-------------|
| `fit_offset_5prime` | 1 = no offset, 2 = 1-bp remnant, 3 = 2-bp remnant |
| `fit_offset_3prime` | Same for 3' end |
| `position_0_artifact_5prime` | True if pos 0 depleted but pos 1 enriched |
| `position_0_artifact_3prime` | Same for 3' end |

**Interpretation:** `fit_offset > 1` does not indicate damage failure — it means the
damage is real but the reported `d_max` is taken from the corrected start position.
Callers should use `d_max_combined` (which already incorporates the offset) rather than
`damage_rate_5prime[0]` directly.

---

## Codon-position-aware damage

**Background:** Not all positions within a codon are equally susceptible to creating a
damaging amino acid change. Position 3 (wobble) changes are often synonymous; positions 1
and 2 (non-wobble) changes are non-synonymous or stop-codon-generating. Tracking C→T
damage by codon position reveals whether terminal deamination disproportionately affects
functionally constrained positions.

**What is measured:** Within the first 15 bases at the 5′ end (and the last 15 at the 3′
end), each cytosine and observed thymine is classified by its position within a codon
(offset modulo 3: position 0, 1, or 2). The T/(T+C) rate is computed per codon position.

**Key outputs:**

| Field | Description |
|-------|-------------|
| `codon_pos_t_rate_5prime[3]` | C→T rate at codon positions 0, 1, 2 at 5′ end |
| `codon_pos_a_rate_3prime[3]` | G→A rate at codon positions 0, 1, 2 at 3′ end |

**Interpretation:** In a strongly damaged ancient library, wobble positions typically show
higher apparent T/(T+C) than non-wobble positions because stop-codon-creating deamination
events at non-wobble positions are under negative selection in the source organism's
evolutionary history. A flat rate across all three positions indicates either very recent
damage or a library with low coding content.

---

## GC-conditional damage bins

**Background:** Read GC content covaries with sequencing depth, mapping efficiency, and
certain library-preparation artefacts. A high-GC composition can create a spurious
T/(T+C) excess even in undamaged libraries, because cytosines are more abundant and
any systematic C-calling error shows up more prominently. Binning reads by GC content
and fitting damage rates separately per bin controls for this composition confound.

**What is measured:** Each read is classified into one of 10 GC-content bins (0–9, where
bin $b$ covers reads with GC fraction in $[b/10, (b+1)/10)$). Within each bin, the
standard exponential decay model is fitted independently, yielding a bin-specific d_max
and baseline.

**Key outputs** (per bin in `gc_bins[10]`):

| Field | Description |
|-------|-------------|
| `d_max` | Fitted D_max for this GC bin |
| `baseline_tc` | Interior T/(T+C) baseline for this GC bin |
| `p_damaged` | Posterior P(damaged) for reads in this bin |
| `valid` | Whether the bin has sufficient coverage for a reliable estimate |

The helper methods `get_gc_bin(seq)`, `get_gc_params(seq)`, and `get_effective_damage(seq,
prob)` provide per-read access to the bin assignments and GC-corrected damage estimates.

**Interpretation:** Consistent d_max across all populated GC bins is strong evidence that
the damage signal is genuine and not a GC-composition artefact. A pattern where only
high-GC or low-GC bins show elevated d_max suggests a composition confound.

---

## Channel summary

| Channel | Measures | Chemistry | Read position | Used for |
|---------|----------|-----------|---------------|----------|
| A (ct5/ga3/ga0/ct3) | C→T and G→A rate | Cytosine deamination | Terminal, exponential decay | Primary damage quantification; position masking in fqdup |
| B / B₃′ | Stop codon C→T / G→A | Deamination, triplet-context | Terminal, exponential decay | Cross-validation of Channel A; artifact detection |
| C | G→T stop codon rate | 8-oxoG oxidation | Uniform across read | Uniformity test; oxidative damage detection |
| D | G→T and C→A rate | 8-oxoG oxidation | Uniform across read | Oxidative damage quantification; artifact flag |
| E | Purine enrichment at 5′ | Depurination / AP-site fragmentation | 5′ terminal | Ancient origin confirmation; sample characterisation |
| CpG split | dmax per CpG/non-CpG context | Methylation-enhanced deamination | 5′ terminal + interior | Methylation signal; cpg_ratio diagnostic |
| Interior clustering | CT co-occurrence at d=1–10 | Clustered deamination | Read interior | Interior damage detection; short_z statistic |
| oxoG 16-ctx | G→T asymmetry per trinucleotide | 8-oxoG specificity | Interior | Context specificity of oxidation |
| Codon-position | C→T rate at codon positions 0/1/2 | Deamination at coding positions | 5′ and 3′ terminal | Wobble vs non-wobble damage; selection bias |
| GC bins | Per-GC-bin d_max and baseline | Composition control | All positions | GC-composition artefact rejection |
| Hexamer | T/(T+C) over hexamer window | Deamination | 5′ terminal (pos 1–6) | Robust terminal enrichment; adapter-bias resistant |
| Briggs model | δ_s, δ_d, R² | Biophysical deamination model | Both termini | Comparable to Briggs 2007; goodness of fit |
| Joint BIC | ΔBIC, Bayes factor, P(damage) | Multi-channel evidence integration | Whole read | Single-number damage evidence score |
| Mixture EM | K, π_ancient, d_ancient | Ancient/modern mixture | Whole library | Ancient fraction estimation |
| metaDMG proxy | d_metamatch, alignability | Alignability-weighted estimate | Whole read | Reference-free metaDMG comparison |
| Adapter offset | fit_offset, pos-0 artifact flag | Adapter remnant detection | Position 0 | Corrects d_max for 1–2 bp adapter remnants |

---

## CpG-like context split

**Background:** Methylated cytosines (5mC) deaminate to thymine several-fold faster than
unmethylated cytosines under the same temperature and pH conditions. In ancient organisms
where CpG dinucleotides are methylated, this produces a distinct C→T excess at positions
followed by G (CpG context) compared to positions where the next base is not G.

**What is measured:** For each of the 15 5′ terminal positions, the T/(C+T) fraction is
tallied separately for positions in a CpG context (next base = G) and non-CpG context
(next base ≠ G). A reference-free interior baseline (middle third of each read) provides
the background T fraction in each context. The 1D MLE then fits the amplitude `d` of an
exponential decay model with the library's existing lambda:

$$\mu(p) = b_{\rm ctx} + (1 - b_{\rm ctx}) \cdot d_{\rm ctx} \cdot e^{-\lambda p}$$

where $b_{\rm ctx}$ is the interior T fraction in that context.

**Key outputs** (`deamination.cpg_like` in JSON):

| Field | Description |
|-------|-------------|
| `dmax_ct5_cpg` | 5′ C→T amplitude in CpG context (1D MLE) |
| `dmax_ct5_noncpg` | 5′ C→T amplitude in non-CpG context (1D MLE) |
| `cpg_ratio` | `dmax_cpg / dmax_noncpg`; >1 indicates methylation-enhanced deamination |
| `log2_cpg_ratio` | log₂ of cpg_ratio; 0 = no context dependence |
| `baseline_cpg` | Interior T/(C+T) at CpG positions (background) |
| `baseline_noncpg` | Interior T/(C+T) at non-CpG positions |

`dmax_ct5_cpg` and `dmax_ct5_noncpg` are reported as `null` when the sample has
insufficient signal (d_max < threshold), or when the MLE converges to the boundary
(flat likelihood), indicating the estimate is uninformative.

**Biological interpretation:** DS aDNA from organisms with methylated CpGs (most
eukaryotes) typically shows `cpg_ratio > 1`. Low values near 1 can indicate unmethylated
CpGs (plants, some invertebrates), modern contamination, or very old samples where
differential context signal is eroded.

---

## Interior C→T clustering

**Background:** Beyond terminal single-stranded overhangs, deamination can occur in
short interior micro-domains where local melting or secondary structure exposes
cytosines. When multiple cytosines in a stretch are co-deaminated in the same event,
adjacent T's in the read are more common than expected from independent site-by-site
deamination. This is distinct from the exponential-decay terminal signal and requires a
pairwise co-occurrence test to detect.

**What is measured:** Within the interior of each read (middle third, read length ≥ 30),
non-CpG `{C,T}` sites are identified. For each pair of eligible sites at separation
d = 1–10 bp, the observed fraction of pairs where both are T is compared to the expected
fraction under site independence (binomial product using the within-read T fraction).
An AG control track (A/(A+G) co-occurrence with the preceding base ≠ C) corrects for
strand-composition and mapping biases.

**Key outputs** (`interior_ct_cluster` in JSON):

| Field | Description |
|-------|-------------|
| `short_z` | Normalised CT-vs-AG contrast statistic (d=1–5 summed) |
| `short_asym_log2oe` | log₂(CT obs/exp) − log₂(AG obs/exp); AG-corrected effect size |
| `short_log2oe` | log₂(CT obs/exp) without AG correction |
| `sep_log2oe[10]` | log₂(obs/exp) for each separation d=1–10 |
| `reads_used` | Reads that contributed ≥ 2 eligible non-CpG {C,T} sites |

**Interpretation:** `short_z > 3` with positive `short_asym_log2oe` indicates excess
CT co-occurrence not explained by strand composition. In practice this signal is weak for
most aDNA samples; its primary value is as a supplementary diagnostic for samples with
unusual interior damage patterns (e.g. very old permafrost material).

---

## 8-oxoG 16-context panel

**Background:** The overall `s_oxog` statistic summarises the G→T strand asymmetry at
interior positions. However, 8-oxoG formation has known sequence-context preferences
(GG-rich contexts are more susceptible). Splitting by flanking dinucleotide reveals
whether oxidation follows the context pattern expected for genuine ancient damage or
is context-uniform (suggesting modern contamination from H₂O₂ exposure, for instance).

**What is measured:** For each interior read position where the base is G (or A on the
complementary strand), the flanking trinucleotide context N**G**N is recorded and the
G→T (or A→C on RC) count is accumulated in one of 16 bins (4 × 4 combinations of
left and right flanking base).

**Key output** (`complement_asymmetry.s_oxog_16ctx[16]` in JSON):

A 16-element float array. Index `4*enc(left) + enc(right)` where `enc(A)=0, enc(C)=1,
enc(G)=2, enc(T)=3`. Elements are `null` when coverage < 500 observations in that bin.

**Interpretation:** Ancient 8-oxoG typically shows broad enrichment across contexts with
slight bias towards GG contexts. A single-context spike suggests a context-specific
artifact. Comparing the 16-context pattern across samples from the same site can
distinguish genuine oxidative damage from prep-introduced oxidation.

---

## Library-type classification from Channel A sub-channels

The four Channel A sub-channels (ct5, ga3, ga0, ct3) together determine which of six
biological damage patterns best describes the library. Which sub-channels are active
depends on how the library was prepared.

### DS: double-stranded (symmetric deamination)

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

fqdup uses a canonical hash, `min(hash(seq), hash(rc(seq)))`, so that a molecule and
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
