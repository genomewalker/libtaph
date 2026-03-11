# Damage types and channels

Ancient DNA carries a record of post-mortem chemistry in its substitution patterns.
libdart-damage tracks five core biochemical damage channels (A–E) plus supplementary
analyses: CpG-context split, interior CT clustering, 8-oxoG 16-context panel, and a
library-type BIC classifier. The channels cross-validate each other, distinguish genuine
ancient damage from modern library-prep artefacts, and together feed the library-type
classifier used by fqdup for damage-aware deduplication.

This document describes the biology of each channel and the JSON output fields they
produce. Fields that are computed internally but not written to JSON (hexamer, Briggs
model, joint BIC, GC bins, mixture model, codon-position tracking) are described briefly
in the [Internal analyses](#internal-analyses) section.

---

## The five biochemical channels

### Channel A: Cytosine deamination (C→T / G→A)

**Chemistry:** Cytosine spontaneously deaminates to uracil, which the polymerase reads as
thymine. The rate is fastest at single-stranded regions near the ends of fragmented
molecules and decays exponentially toward the interior.

**What is measured:**

| Sub-channel | Signal | Where |
|---|---|---|
| ct5 | T/(T+C) excess above baseline | 5′ end, positions 0–14 |
| ga3 | A/(A+G) excess above baseline, smooth decay | 3′ end, positions 1–14 |
| ga0 | A/(A+G) excess at position 0 only, spike | 3′ position 0 |
| ct3 | T/(T+C) excess above baseline | 3′ end, positions 0–14 |

The decay model at each end is:

$$\delta(p) = b + A \cdot e^{-\lambda p}$$

where $b$ is the interior background rate, $A$ is the damage amplitude, $\lambda$ is the
decay constant, and $p$ is the distance from the terminus (0-indexed).

The calibrated damage rate at position 0 is:

$$D_{\max} = \frac{A}{1 - b}$$

**JSON output** (`deamination` block):

| Field | Description |
|-------|-------------|
| `d_max_5prime` | Fitted $D_{\max}$ at 5′ end (background-corrected amplitude) |
| `d_max_3prime` | Fitted $D_{\max}$ at 3′ end |
| `d_max_combined` | Best overall damage estimate — see [source](#d_max-source) |
| `d_metamatch` | Alignability-weighted $D_{\max}$; comparable to metaDMG output |
| `source` | How `d_max_combined` was derived (see below) |
| `lambda_5prime` | Fitted decay constant at 5′ end |
| `lambda_3prime` | Fitted decay constant at 3′ end |
| `bg_5prime` | Interior T/(T+C) baseline at 5′ end |
| `bg_3prime` | Interior A/(A+G) baseline at 3′ end |
| `validated` | `true` when Channel B corroborates Channel A |
| `artifact` | `true` when Channel A fires but Channel B contradicts it |
| `per_pos_5prime_ct[15]` | Raw T/(T+C) at 5′ positions 0–14 |
| `per_pos_3prime[15]` | Raw A/(A+G) at 3′ positions 0–14 (DS) or T/(T+C) (SS-original) |

**`d_max_combined` source values:**

| `source` value | Meaning |
|---|---|
| `"5prime_only"` | 3′ inverted — 5′ value used |
| `"3prime_only"` | 5′ inverted — 3′ value used |
| `"average"` | Both ends valid: average (or mixture/joint model estimate) |
| `"max_ss_asymmetry"` | SS library: `max(d_max_5, d_max_3)` |
| `"min_asymmetry"` | DS with high asymmetry: conservative `min(d_max_5, d_max_3)` |
| `"channel_b_structural"` | Channel B codon estimate used (most reliable) |
| `"none"` | Both ends inverted with no recoverable signal; `d_max_combined = 0` |

**Adapter artifact:** When a 1–2 bp adapter remnant depresses position 0, the fit window
is shifted to start at position 1 or 2. `d_max_combined` already incorporates this
correction; do not use `per_pos_5prime_ct[0]` directly for artifact-affected libraries.

**Diagnostic role:** Primary damage quantification channel. Used by fqdup for position
masking. The ct3 sub-channel is a negative control for DS libraries (should be flat) and a
damage signal for SS-original libraries. The ga0 spike at position 0 is a distinct signal
from the smooth ga3 decay (see [library classification](#library-type-classification)).

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

**JSON output:** Channel B result is reflected in `deamination.validated` and
`deamination.artifact`:

- Both channels detect signal → `validated: true`
- Channel A fires, Channel B contradicts → `artifact: true`
- Insufficient coverage → no flag set (soft suppression only)

---

### Channel C: Oxidative stop codon tracking (G→T uniformity test)

**Chemistry:** 8-Oxoguanine (8-oxoG) forms when guanine is oxidised. The polymerase
misreads 8-oxoG as adenine, producing G→T transversions. Unlike cytosine deamination,
oxidation occurs throughout the molecule, not preferentially at the ends.

**What is measured:** Stop codon conversion via G→T at glutamate and glycine codons
(GAG/GAA/GGA → TAG/TAA/TGA), measured at terminal vs. interior positions.

The key diagnostic statistic is the uniformity ratio:

$$U_C = \frac{\text{stop rate (terminal)}}{\text{stop rate (interior)}}$$

$U_C \approx 1$ → genuine oxidative damage (uniform distribution).
$U_C \gg 1$ → terminal enrichment, suggesting co-occurrence with deamination.

**JSON output** (`complement_asymmetry.channel_c_detected`):
`true` when Channel C detects statistically significant G→T stop codon conversion.

---

### Channel D: Direct G→T / C→A transversion tracking (8-oxoG rate)

**Chemistry:** Same 8-oxoG chemistry as Channel C, measured directly as raw transversion
frequency rather than through stop codon context. G→T and its complement C→A are tracked
at terminal positions and compared to the interior baseline.

**JSON output** (`complement_asymmetry` block):

| JSON field | Description |
|------------|-------------|
| `tg_interior` | G→T background rate in read interior |
| `ac_interior` | C→A background rate in read interior |
| `tg_terminal` | G→T rate at 5′ terminal positions |
| `ac_terminal` | C→A rate at 5′ terminal positions |
| `gt_bg_fitted` | Fitted background for the exponential G→T model |
| `gt_term_fitted` | Fitted terminal amplitude |
| `gt_decay_fitted` | Fitted decay constant |
| `s_gt` | Strand asymmetry score (G→T vs C→A contrast) |
| `D` | Overall 8-oxoG asymmetry index |
| `per_pos_5prime_gt[15]` | Raw G→T fraction at 5′ positions 0–14 |
| `s_oxog` | G→T strand asymmetry at interior positions |
| `se_s_oxog` | Standard error of `s_oxog` |

---

### Channel E: Depurination (AP-site enrichment at strand breaks)

**Chemistry:** Purines (adenine and guanine) are more susceptible than pyrimidines to
hydrolytic base loss (depurination) under acidic and warm conditions. The resulting
apurinic (AP) site is a strand-break precursor; fragmentation preferentially occurs at AP
sites, leaving purines enriched at the newly formed 5′ ends.

**JSON output** (`depurination` block):

| JSON field | Description |
|------------|-------------|
| `detected` | `true` when 5′ purine enrichment is statistically significant |
| `enrichment_5prime` | Purine excess at 5′ terminal vs. interior |
| `enrichment_3prime` | Purine excess at 3′ terminal vs. interior |
| `rate_interior` | Interior A+G fraction (baseline) |

**Diagnostic role:** Independent evidence of genuine ancient fragmentation. Not used by
fqdup for masking; reported for sample characterisation.

---

## Supplementary analyses

### CpG-like context split

**Background:** Methylated cytosines (5mC) deaminate several-fold faster than unmethylated
cytosines. In ancient organisms where CpG dinucleotides are methylated, this produces a
distinct C→T excess at positions followed by G (CpG context) vs. positions where the next
base is not G.

**What is measured:** For each of the 15 5′ terminal positions, T/(C+T) is tallied
separately for CpG (next base = G) and non-CpG contexts. A reference-free interior
baseline (middle third of each read) provides the background T fraction. 1D MLE fits:

$$\mu(p) = b_{\rm ctx} + (1 - b_{\rm ctx}) \cdot d_{\rm ctx} \cdot e^{-\lambda p}$$

**JSON output** (`deamination.cpg_like` block):

| Field | Description |
|-------|-------------|
| `dmax_ct5_cpg` | 5′ C→T amplitude in CpG context |
| `dmax_ct5_noncpg` | 5′ C→T amplitude in non-CpG context |
| `cpg_ratio` | `dmax_cpg / dmax_noncpg`; >1 = methylation-enhanced deamination |
| `log2_cpg_ratio` | log₂(cpg_ratio); 0 = no context dependence |
| `baseline_cpg` | Interior T/(C+T) at CpG positions |
| `baseline_noncpg` | Interior T/(C+T) at non-CpG positions |
| `cov_terminal_cpg` | Total T+C observations at CpG terminal positions |
| `cov_terminal_noncpg` | Total T+C observations at non-CpG terminal positions |
| `effcov_terminal_cpg` | Effective coverage: `cov * (1 − baseline)` |
| `effcov_terminal_noncpg` | Effective coverage: `cov * (1 − baseline)` |

All fields are `null` when the MLE has insufficient signal or converges to its boundary.

**Biological interpretation:** Most eukaryotes with methylated CpGs show `cpg_ratio > 1`.
Values near 1 indicate unmethylated CpGs, modern contamination, or very old samples where
the differential context signal has eroded.

---

### Interior C→T clustering

**Background:** Beyond terminal single-stranded overhangs, deamination can occur in short
interior micro-domains. Co-deamination of adjacent cytosines in the same event produces
more adjacent T's than expected from independent site-by-site deamination.

**What is measured:** Within the interior of each read (middle third, read length ≥ 30),
non-CpG `{C,T}` sites are identified. For each pair at separation d = 1–10 bp, observed
co-T fraction is compared to expected under site independence. An AG control track corrects
for strand-composition bias.

**JSON output** (`interior_ct_cluster` block):

| Field | Description |
|-------|-------------|
| `short_z` | Normalised CT-vs-AG contrast (d=1–5 summed) |
| `short_asym_log2oe` | log₂(CT obs/exp) − log₂(AG obs/exp); AG-corrected effect size |
| `short_log2oe` | log₂(CT obs/exp) without AG correction |
| `short_obs` | Total observed co-T pairs (d=1–5) |
| `short_exp` | Total expected co-T pairs under independence |
| `reads_used` | Reads contributing ≥ 2 eligible non-CpG {C,T} sites |
| `sep_log2oe[10]` | log₂(obs/exp) for each separation d=1–10 |

`short_z > 3` with positive `short_asym_log2oe` indicates excess interior CT
co-occurrence not explained by strand composition.

---

### 8-oxoG 16-context panel

**Background:** 8-oxoG formation has known sequence-context preferences. Splitting the
G→T asymmetry by flanking dinucleotide reveals whether oxidation follows the context
pattern of genuine ancient damage or is context-specific (suggesting a preparation
artefact).

**What is measured:** For each read position where the base is G, the flanking trinucleotide
context N**G**N is recorded and the G→T count is accumulated in one of 16 bins
(4×4 combinations of left and right flanking base). Index = `4*enc(left) + enc(right)`
where `enc(A)=0, enc(C)=1, enc(G)=2, enc(T)=3`.

**JSON output** (`complement_asymmetry` block):

| Field | Description |
|-------|-------------|
| `s_oxog_16ctx[16]` | Per-context G→T strand asymmetry; `null` when coverage < 500 |
| `cov_oxog_16ctx[16]` | Read coverage per context bin |

Ancient 8-oxoG shows broad enrichment across contexts with slight GG-context bias.
A single-context spike suggests a context-specific artefact.

---

### Library-type BIC classifier

**Background:** The four Channel A sub-channels (ct5, ga3, ga0, ct3) are jointly fitted
under seven biological models (M_bias, M_DS_symm, M_DS_spike, M_DS_symm_art, M_SS_comp,
M_SS_orig, M_SS_asym). The model with the lowest joint BIC is selected.

**JSON output** (`bic` block):

| Field | Description |
|-------|-------------|
| `bias` | BIC of the null model (M_bias: all channels flat) |
| `ds` | BIC of the best-fitting DS model |
| `ss` | BIC of the best-fitting SS model |
| `ct5_amp` | Fitted ct5 amplitude used in classification |
| `ga3_amp` | Fitted ga3 amplitude |
| `ga0_amp` | Fitted ga0 amplitude |
| `ct3_amp` | Fitted ct3 amplitude |

Lower BIC = better fit. `library_type` is `"double-stranded"` when `ds < ss`,
`"single-stranded"` when `ss < ds`, and `"unknown"` when neither beats `bias`.

The classification result is at the top level of the JSON: `library_type`,
`library_type_auto`, `library_type_rescued`.

---

## Channel summary

| Channel | Measures | Chemistry | JSON block |
|---------|----------|-----------|-----------|
| A (ct5/ga3/ga0/ct3) | C→T and G→A rate | Cytosine deamination | `deamination` |
| B / B₃′ | Stop codon C→T / G→A | Deamination, triplet-context | `deamination.validated/artifact` |
| C | G→T stop codon uniformity | 8-oxoG oxidation | `complement_asymmetry.channel_c_detected` |
| D | G→T and C→A transversion rate | 8-oxoG oxidation | `complement_asymmetry` |
| E | Purine enrichment at 5′ | Depurination / AP-site fragmentation | `depurination` |
| CpG split | dmax per CpG/non-CpG context | Methylation-enhanced deamination | `deamination.cpg_like` |
| Interior clustering | CT co-occurrence at d=1–10 | Clustered interior deamination | `interior_ct_cluster` |
| 8-oxoG 16-ctx | G→T asymmetry per trinucleotide | 8-oxoG context specificity | `complement_asymmetry.s_oxog_16ctx` |
| Library BIC | BIC per model class | Four-channel joint fit | `bic` |

---

## Library-type classification from Channel A sub-channels

The four Channel A sub-channels (ct5, ga3, ga0, ct3) together determine which of six
biological damage patterns best describes the library.

### DS: double-stranded (symmetric deamination)

**Pattern:** ct5 ≈ ga3, both smooth exponential decay. **Active:** ct5 + ga3.

### DS + end-repair artifact

**Pattern:** ct5 + ga3 (smooth) + ga0 spike at position 0. **Active:** ct5 + ga3 + ga0.

### SS complement-orientation

**Pattern:** ga0 spike only; ct5 and ga3 flat. **Active:** ga0.

### SS original-orientation

**Pattern:** Both 5′ (ct5) and 3′ (ct3) show C→T excess; no ga3/ga0.
**Active:** ct5 + ct3.

### SS mixed orientations

**Pattern:** ct5 from original-orientation reads + ga0 from complement-orientation reads.
**Active:** ct5 + ga0.

### UNKNOWN

No Channel A sub-channel beats the null model. The library has no detectable deamination,
or coverage is insufficient for any channel to reach significance. UNKNOWN is the correct
call for zero-damage libraries — it is not an error.

---

## How fqdup uses damage channels

fqdup deduplicates reads by hashing sequences. Damage-aware hashing corrects for PCR
copies of the same molecule that carry deamination at different positions, by masking
affected terminal positions before hashing.

### Which channels drive masking

Only **Channel A** drives position masking. Channels B–E are diagnostic: they validate
that Channel A is genuine ancient damage, flag artefacts, and characterise the sample, but
do not directly change which positions are masked.

### Revcomp-equivariant masking

fqdup uses a canonical hash `min(hash(seq), hash(rc(seq)))`. The mask must satisfy:

$$\text{mask}(\text{rc}(\text{seq})) = \text{rc}(\text{mask}(\text{seq}))$$

Enforced by symmetric rule: masking position $p$ from 5′ also masks position $p$ from 3′.

### Masking threshold

$$\text{mask position } p \iff \delta(p) - b > \tau \quad (\text{default } \tau = 5\%)$$

### Library type and masking

| Library type | Positions masked |
|---|---|
| DS | 5′ pos $p$ and 3′ pos $p$ symmetrically (ct5/ga3) |
| DS + artifact | Same as DS; ga0 spike within masked region |
| SS complement | 3′ pos 0 (ga0) only |
| SS original | 5′ and 3′ symmetrically (ct5/ct3) |
| SS mixed | 5′ from ct5 + 3′ pos 0 from ga0 |
| UNKNOWN | None — standard hashing |

### Artifact detection

`artifact: true` in the JSON means Channel A fired but Channel B contradicted it
(likely GC-composition bias). fqdup treats damage as unvalidated and optionally suppresses
masking, preventing modern contamination from being over-collapsed.

---

## Internal analyses

The following analyses are computed during damage estimation and used internally (e.g. to
determine `d_max_combined` or detect adapter artefacts) but are **not written to the JSON
output**.

- **Hexamer-based detection:** T/(T+C) averaged over positions 1–6; provides a
  position-0-bias-resistant corroborating statistic for Channel A.
- **Briggs biophysical model:** Re-fits Channel A data as δ_s (single-stranded overhang
  rate) and δ_d (double-stranded background rate) with R² goodness of fit.
- **Joint probabilistic model (BIC + Bayes factor):** Bayesian comparison of M₁ (damage
  present) vs M₀ (no damage); produces `joint_delta_max`, `joint_p_damage`, and Bayes
  factor. Used internally to set `damage_status` and as the source for `d_max_combined`
  when other pathways are unavailable.
- **GC-conditional damage bins:** Independent exponential fits across 10 GC-content bins;
  used to detect composition artefacts and as a GC-adjusted `d_max` for `d_metamatch`.
- **Mixture model (K-component EM over GC bins):** Separates ancient and modern components
  by fitting K damage rates over GC-stratified reads. `mixture_d_reference` (E[δ | GC >
  50%]) feeds into `d_max_combined` when the model converges.
- **Codon-position-aware tracking:** C→T rate split by codon position (0/1/2 mod 3);
  used as supplementary wobble-position diagnostic.
- **Adapter offset detection:** Per-end fit window offset (1/2/3); used internally to
  shift the exponential fit window and set `source` when artefacts are detected.
  `d_max_combined` already incorporates any offset correction.
