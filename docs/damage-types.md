# Damage types and channels

libdart-damage scans raw FASTQ reads for six biochemical damage channels and writes a
structured JSON report. Damage values are fractions in [0, 1]; multiply by 100 for
percentages.

Ancient DNA degradation falls into three broad categories: **hydrolytic damage**
(cytosine deamination; depurination and AP-site fragmentation), **oxidative damage**
(8-oxoguanine formation), and **library-preparation artefacts** that mimic or confound
these signals. Each channel targets one of these processes. Together they cross-validate
that the detected signal is genuine ancient damage rather than a composition or
preparation artefact.

---

## JSON output reference

### Top-level fields

| Field | Type | Description |
|-------|------|-------------|
| `input` | string | Input FASTQ path |
| `n_reads` | int | Total reads scanned |
| `library_type` | string | `"double-stranded"`, `"single-stranded"`, or `"unknown"` |
| `library_type_auto` | bool | `true` if type was auto-detected (not forced by `--library-type`) |
| `library_type_rescued` | bool | `true` if a rescue rule overrode the primary BIC winner |
| `damage_status` | string | `"present"`, `"weak"`, or `"absent"` |

---

### `deamination` block

Cytosine spontaneously deaminates to uracil by hydrolysis of the C4 amino group.
In an intact duplex, uracil-DNA glycosylase removes the lesion efficiently. In the
single-stranded overhangs that form at the ends of nicked or fragmented ancient molecules,
repair is unavailable: the uracil persists and templates thymine during PCR amplification,
producing the characteristic C→T excess in ancient DNA reads [[Briggs et al.
2007](#references); [Lindahl 1993](#references)]. The rate of spontaneous deamination in double-stranded DNA is approximately 10⁻¹³ s⁻¹
per cytosine under physiological conditions (~100–500 events per cell per day), and is
accelerated by low pH and elevated temperature [[Shen et al. 1994](#references); [Lindahl
1993](#references)]; cold, alkaline burial
slows it by orders of magnitude, explaining DNA survival over millennia.

Because only the single-stranded overhang is susceptible, the excess decays exponentially
inward from the fragment end [[Jónsson et al. 2013](#references)]:

$$\delta(p) = b + A \cdot e^{-\lambda p}$$

where $p$ is distance from the terminus (0-indexed), $b$ is the interior background rate,
$A$ is the terminal excess above background, and $\lambda$ is the decay constant. `d_max`
is defined as $A/(1-b)$ — the fraction of terminal cytosine sites that were converted to
thymine, after accounting for the background T proportion.

#### Damage estimates

| Field | Description |
|-------|-------------|
| `d_max_5prime` | C→T excess at the 5′ terminus above the interior background |
| `d_max_3prime` | G→A excess at the 3′ terminus (DS) or C→T excess (SS original-strand) above background |
| `d_max_combined` | Best single damage estimate for this library; `source` records how it was derived |
| `d_metamatch` | GC-composition-weighted estimate comparable to metaDMG's D_max [[Michelsen et al. 2022](#references)] |
| `source` | How `d_max_combined` was derived — see table below |
| `lambda_5prime` | Exponential decay rate at the 5′ end; larger values mean shorter overhang length |
| `lambda_3prime` | Exponential decay rate at the 3′ end |
| `bg_5prime` | Interior T/(T+C) baseline at the 5′ end; reflects sequence composition, not damage |
| `bg_3prime` | Interior A/(A+G) baseline at the 3′ end |

**`source` values:**

| Value | Meaning |
|-------|---------|
| `"5prime_only"` | 3′ signal inverted by adapter artefact; `d_max_combined = d_max_5prime` |
| `"3prime_only"` | 5′ signal inverted; `d_max_combined = d_max_3prime` |
| `"average"` | Both ends valid; weighted average or mixture model |
| `"max_ss_asymmetry"` | SS library with unequal end amplitudes; `max(d_max_5, d_max_3)` |
| `"min_asymmetry"` | DS with high end asymmetry; conservative `min(d_max_5, d_max_3)` |
| `"channel_b_structural"` | Channel B stop-codon estimate; used when more reliable than Channel A |
| `"none"` | Both ends inverted; no recoverable signal; `d_max_combined = 0` |

Adapter remnants of 1–2 bp can depress T/(T+C) at position 0, hiding the biological
signal. When this is detected, the fit window is shifted automatically and `d_max_combined`
already reflects the correction. Do not read `per_pos_5prime_ct[0]` directly for such
libraries.

#### Validation flags

Channel B tests the same C→T deamination through stop-codon context
(CAA/CAG/CGA → TAA/TAG/TGA), which is insensitive to GC-composition bias. It provides a
composition-independent cross-check on Channel A.

| Field | Description |
|-------|-------------|
| `validated` | Channel B corroborates Channel A deamination |
| `artifact` | Channel A signal present but Channel B contradicts it; likely GC-composition artefact rather than genuine damage |

#### Per-position rates

| Field | Description |
|-------|-------------|
| `per_pos_5prime_ct[15]` | Raw T/(T+C) at 5′ positions 0–14 |
| `per_pos_3prime[15]` | A/(A+G) at 3′ positions 0–14 for DS; T/(T+C) for SS original-strand |

Values of `-1` indicate insufficient coverage at that position.

#### `cpg_like` sub-block — methylation signal

When cytosine is methylated at the C5 position (5-methylcytosine, 5mC), its deamination
product is thymine rather than uracil [[Shen et al. 1994](#references)]. Deamination of 5mC in double-stranded DNA occurs at approximately twice the rate of
unmodified cytosine [[Shen et al. 1994](#references)]. An additional factor is that the
thymine product creates a T:G mismatch, which is repaired by thymine-DNA glycosylase (TDG)
less efficiently than the U:G mismatch from cytosine deamination is repaired by
uracil-DNA glycosylase [[Wiebauer & Jiricny 1990](#references); [Neddermann & Jiricny
1994](#references)]. In ancient DNA, this produces an
enhanced C→T rate specifically at positions where the downstream base is G. A `cpg_ratio`
above 1 indicates methylation-dependent deamination, the expected pattern for vertebrates
and most eukaryotes with CpG methylation. Values near 1
indicate either unmethylated CpGs (e.g., plants with low global methylation, many fungi),
modern contamination, or signal too eroded to detect context dependence.

| Field | Description |
|-------|-------------|
| `dmax_ct5_cpg` | 5′ C→T amplitude in CpG context |
| `dmax_ct5_noncpg` | 5′ C→T amplitude in non-CpG context |
| `cpg_ratio` | `dmax_cpg / dmax_noncpg`; values >1 indicate methylation-enhanced deamination |
| `log2_cpg_ratio` | log₂(cpg_ratio); 0 means no context dependence |
| `baseline_cpg` | Interior T/(C+T) at CpG positions |
| `baseline_noncpg` | Interior T/(C+T) at non-CpG positions |
| `cov_terminal_cpg` | Total T+C observations at CpG terminal positions |
| `cov_terminal_noncpg` | Total T+C at non-CpG terminal positions |
| `effcov_terminal_cpg` | Effective coverage: `cov × (1 − baseline)` |
| `effcov_terminal_noncpg` | Effective coverage: `cov × (1 − baseline)` |

All fields are `null` when signal is insufficient or the amplitude estimate hits its
boundary.

---

### `bic` block — library-type classifier

Library preparation determines which damage channels are active, so the BIC model
selection that drives classification also reveals which sub-channels contributed signal.
This block exposes the classifier's evidence for inspecting borderline calls.

Double-stranded (DS) library protocols (e.g., NEBNext, TruSeq [[Meyer & Kircher
2010](#references)]) ligate adapters to both strands. Reads from the original strand show
C→T at the 5′ end; reads from the complementary strand map the same damaged cytosines as
G→A at the 3′ end. Single-stranded (SS) protocols (e.g., Gansauge & Meyer [[2013](#references)];
SRSLY [[Gansauge et al. 2017](#references)]) circularize and sequence only one strand per
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

### `complement_asymmetry` block — oxidative damage (Channels C and D)

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
| `s_gt` | G→T vs C→A strand asymmetry contrast; near zero for balanced DS oxidation, informative for SS |
| `D` | Overall 8-oxoG asymmetry index |
| `per_pos_5prime_gt[15]` | Raw G→T fraction at 5′ positions 0–14 |

`channel_c_detected: true` when G→T stop-codon conversions (GAG/GAA/GGA → stop) are
detected at levels consistent with uniform oxidation throughout the read rather than
terminal enrichment from co-occurring deamination.

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

### `interior_ct_cluster` block — interior co-deamination

Under prolonged exposure conditions, cytosines within single-stranded micro-domains
— short stretches that transiently unpair during storage — can be co-deaminated in the
same hydrolytic event. This produces co-occurring T residues at nearby positions more
often than expected under site independence. The signal is in the read interior (middle
third), away from the terminal overhangs that drive Channel A.

**Chargaff's second parity rule and the AG control.** Chargaff's second (intra-strand)
parity rule states that within a single strand, [A] ≈ [T] and [C] ≈ [G]
[[Rudner et al. 1968](#references)]. This symmetry arises from the evolutionary history
of symmetric mutation processes acting on the genome. It underpins the AG co-occurrence
control used here: if elevated T co-occurrence at adjacent {C,T} sites arose purely from
strand composition rather than damage clustering, the G co-occurrence at adjacent {A,G}
sites should show a matched elevation. The contrast between the CT and AG tracks therefore
isolates genuine deamination clustering from composition effects.

For each pair of non-CpG `{C,T}` sites at separation d = 1–10 bp, the observed co-T
fraction is compared to the expectation under site independence.

| Field | Description |
|-------|-------------|
| `short_asym_log2oe` | log₂(CT obs/exp) − log₂(AG obs/exp); composition-corrected effect size |
| `short_log2oe` | log₂(CT obs/exp) without composition correction |
| `short_z` | Normalised CT-vs-AG contrast statistic, summed over d = 1–5 |
| `short_obs` | Observed co-T pairs at d = 1–5 |
| `short_exp` | Expected co-T pairs at d = 1–5 under independence |
| `reads_used` | Reads contributing ≥ 2 eligible non-CpG {C,T} sites |
| `sep_log2oe[10]` | log₂(obs/exp) at each separation d = 1–10 |

`short_z > 3` with positive `short_asym_log2oe` indicates genuine interior co-deamination
beyond what strand composition predicts. In most aDNA samples this signal is weak; it
becomes informative for very old material or samples with unusual interior damage such as
permafrost specimens.

---

### `depurination` block — AP-site fragmentation (Channel E)

The glycosidic bond linking a purine base (adenine or guanine) to the deoxyribose sugar
is susceptible to acid-catalyzed hydrolysis. The rate under physiological conditions is approximately 3 × 10⁻¹¹ per purine per
second (~2,000–10,000 events per cell per day), substantially faster than cytosine
deamination [[Lindahl & Nyberg 1972](#references); [Lindahl 1993](#references)]. The
resulting apurinic (AP) site is a metastable abasic residue: β-elimination converts the
AP deoxyribose to a strand break, fragmenting the molecule at the site of base loss.
Because purines (A and G) are lost preferentially over pyrimidines (C and T), and because
fragmentation at an AP site leaves the removed purine's position at the 5′ end of the
newly formed fragment, the sequenced 5′ termini are enriched for purines — the signature
detected by Channel E. This signal provides evidence of ancient fragmentation independent
of deamination, and is particularly informative for samples where deamination is low but
preservation conditions are known to be acidic or warm [[Dabney et al. 2013](#references)].

| Field | Description |
|-------|-------------|
| `detected` | `true` when 5′ purine enrichment is statistically significant |
| `enrichment_5prime` | Purine excess at 5′ terminal positions relative to interior |
| `enrichment_3prime` | Purine excess at 3′ terminal positions relative to interior |
| `rate_interior` | Interior A+G fraction (composition baseline) |

Channel E is reported for sample characterisation and is not used by fqdup for position
masking.

---

## Library-type classification

Which damage channels are active depends on library preparation. Double-stranded (DS)
protocols ligate adapters to both strands, so reads from both the original and
complementary strand are sequenced. Single-stranded (SS) protocols circularize and
sequence only one strand per molecule; which strand is captured depends on the protocol
and ligation orientation.

| Library type | Active sub-channels | Pattern |
|---|---|---|
| DS | ct5 + ga3 | Symmetric exponential C→T at 5′ and G→A at 3′ |
| DS + end-repair artifact | ct5 + ga3 + ga0 | DS pattern plus isolated 3′ pos-0 spike |
| SS complement-orientation | ga0 | Strong 3′ pos-0 G→A spike; 5′ flat |
| SS original-orientation | ct5 + ct3 | C→T at both 5′ and 3′ ends; no G→A |
| SS mixed orientations | ct5 + ga0 | 5′ C→T decay plus 3′ pos-0 spike |
| UNKNOWN | — | No sub-channel beats the null model |

In DS libraries, the symmetric ct5 ≈ ga3 pattern is a consequence of Chargaff's first
rule: the original strand contributes C→T at its 5′ end; reads from the complementary
strand map the same damaged positions as G→A at the 3′ end. The shared-amplitude DS model
(M_DS_symm) enforces this symmetry as a model constraint, which also means it penalises
libraries with strongly asymmetric ends.

UNKNOWN is the correct classification for zero-damage or near-zero-damage libraries where
the damage pattern contains insufficient information to distinguish library type. fqdup
falls back to standard hashing for UNKNOWN libraries, which has no effect when d_max is
near zero.

---

## Channel biology reference

| Channel | Measures | Chemistry | JSON output |
|---------|----------|-----------|-------------|
| A (ct5/ga3/ga0/ct3) | C→T and G→A terminal rates | Cytosine deamination | `deamination`, `bic` |
| B / B₃′ | Stop codon C→T / G→A frequency | Deamination in triplet context | `validated`, `artifact` |
| C | G→T stop codon uniformity | 8-oxoG oxidation | `channel_c_detected` |
| D | G→T and C→A transversion rates | 8-oxoG oxidation | `complement_asymmetry` |
| E | Purine 5′ enrichment | AP-site fragmentation | `depurination` |
| CpG split | C→T amplitude by CpG / non-CpG context | Methylation-enhanced deamination | `cpg_like` |
| Interior clustering | Adjacent CT co-occurrence in read interior | Clustered interior deamination | `interior_ct_cluster` |
| 8-oxoG 16-ctx | G→T asymmetry by trinucleotide context | 8-oxoG context specificity | `s_oxog_16ctx` |

Channels B–E cross-validate Channel A without contributing to position masking. If
Channel A detects deamination but Channel B contradicts it, the signal is flagged as a
composition artefact rather than genuine ancient damage.

---

## How fqdup uses damage channels

### Position masking

Only Channel A drives position masking. Position p is masked when the damage rate at
that position, after subtracting the interior baseline b, exceeds the threshold τ (default
5%). The mask is applied symmetrically at both ends so that canonical hashing remains
correct under reverse-complement: `min(hash(seq), hash(rc(seq)))` requires
`mask(rc(seq)) = rc(mask(seq))`.

### Masking by library type

| Library type | Positions masked |
|---|---|
| DS | 5′ pos p and 3′ pos p symmetrically (ct5/ga3) |
| DS + artifact | Same as DS; ga0 spike falls within the masked region |
| SS complement | 3′ pos 0 only (ga0) |
| SS original | 5′ and 3′ symmetrically (ct5/ct3) |
| SS mixed | 5′ from ct5 + 3′ pos 0 from ga0 |
| UNKNOWN | None — standard hashing |

### Artifact suppression

When `artifact: true`, Channel A fired but Channel B contradicted it, indicating
GC-composition bias rather than genuine damage. fqdup can optionally suppress masking in
this case to avoid over-collapsing reads from modern contamination.

---

## Internal analyses (not in JSON)

The following analyses are computed during damage estimation. Their results feed into the
reported fields above but are not written to JSON directly.

- **Hexamer detection:** T/(T+C) averaged over positions 1–6; provides position-0-bias-resistant
  corroboration of Channel A that is unaffected by adapter remnants at position 0.
- **Briggs biophysical model:** Re-fits Channel A as δ_s (single-stranded overhang rate) and
  δ_d (double-stranded background rate), following [[Briggs et al. 2007](#references)]. R²
  assesses fit quality.
- **Joint probabilistic model:** Bayesian comparison of damage-present vs damage-absent using
  both BIC and Bayes factor. Drives `damage_status` and provides a fallback estimate for
  `d_max_combined` when the exponential fit is unreliable.
- **GC-conditional damage bins:** Independent exponential fits across 10 GC-content bins.
  Separates genuine deamination from GC-composition artefacts and produces the adjusted
  estimate that feeds into `d_metamatch`.
- **Mixture model (K-component EM):** EM over GC-stratified reads to separate ancient from
  modern components. When converged, `mixture_d_reference` feeds into `d_max_combined`.
- **Codon-position tracking:** C→T rate at codon positions 0, 1, 2; supplementary wobble
  diagnostic.
- **Adapter offset detection:** Per-end fit window shift of 1–3 bp when adapter remnants
  are detected. Applied before computing `d_max_combined`, which already reflects any
  correction.

---

## References

Neddermann P, Jiricny J (1994) Efficient removal of uracil from G·U mispairs by the
mismatch-specific thymine DNA glycosylase from HeLa cells. *Proc Natl Acad Sci USA*
91:1642–1646.
[DOI: 10.1073/pnas.91.5.1642](https://doi.org/10.1073/pnas.91.5.1642)

Briggs AW, Stenzel U, Johnson PLF, et al. (2007) Patterns of damage in genomic DNA
sequences from a Neandertal. *Proc Natl Acad Sci USA* 104:14616–14621.
[DOI: 10.1073/pnas.0704665104](https://doi.org/10.1073/pnas.0704665104)

Cadet J, Douki T, Ravanat J-L (2010) Oxidatively generated base damage to cellular DNA.
*Free Radic Biol Med* 49:9–21.
[DOI: 10.1016/j.freeradbiomed.2010.03.025](https://doi.org/10.1016/j.freeradbiomed.2010.03.025)

Chargaff E (1950) Chemical specificity of nucleic acids and mechanism of their enzymatic
degradation. *Experientia* 6:201–209.
[DOI: 10.1007/BF02173653](https://doi.org/10.1007/BF02173653)

Dabney J, Knapp M, Glocke I, et al. (2013) Complete mitochondrial genome sequence of a
Middle Pleistocene cave bear reconstructed from ultrashort DNA fragments. *Proc Natl Acad
Sci USA* 110:15758–15763.
[DOI: 10.1073/pnas.1314445110](https://doi.org/10.1073/pnas.1314445110)

Gansauge M-T, Meyer M (2013) Single-stranded DNA library preparation for the sequencing
of ancient or damaged DNA. *Nat Protoc* 8:737–748.
[DOI: 10.1038/nprot.2013.038](https://doi.org/10.1038/nprot.2013.038)

Gansauge M-T, Gerber T, Glocke I, et al. (2017) Single-stranded DNA library preparation
from highly degraded DNA using T4 DNA ligase. *Nucleic Acids Res* 45:e79.
[DOI: 10.1093/nar/gkx033](https://doi.org/10.1093/nar/gkx033)

Jónsson H, Ginolhac A, Schubert M, Johnson PLF, Orlando L (2013) mapDamage2.0: fast
approximate Bayesian estimates of ancient DNA damage parameters. *Bioinformatics*
29:1682–1684.
[DOI: 10.1093/bioinformatics/btt193](https://doi.org/10.1093/bioinformatics/btt193)

Michelsen C, Pedersen MW, Fernandez-Guerra A, et al. (2022) metaDMG: A fast and accurate
ancient DNA damage toolkit for metagenomic data. *bioRxiv*.
[DOI: 10.1101/2022.12.06.519264](https://doi.org/10.1101/2022.12.06.519264)

Lindahl T (1993) Instability and decay of the primary structure of DNA. *Nature*
362:709–715.
[DOI: 10.1038/362709a0](https://doi.org/10.1038/362709a0)

Lindahl T, Nyberg B (1972) Rate of depurination of native deoxyribonucleic acid.
*Biochemistry* 11:3610–3618.
[DOI: 10.1021/bi00769a018](https://doi.org/10.1021/bi00769a018)

Meyer M, Kircher M (2010) Illumina sequencing library preparation for highly multiplexed
target capture and sequencing. *Cold Spring Harb Protoc* 2010:pdb.prot5448.
[DOI: 10.1101/pdb.prot5448](https://doi.org/10.1101/pdb.prot5448)

Mitchell D, Bridge R (2006) A test of Chargaff's second rule. *Biochem Biophys Res
Commun* 340:90–94.
[DOI: 10.1016/j.bbrc.2005.11.160](https://doi.org/10.1016/j.bbrc.2005.11.160)

Rudner R, Karkas JD, Chargaff E (1968) Separation of B. subtilis DNA into complementary
strands. III. Direct analysis. *Proc Natl Acad Sci USA* 60:921–922.
[DOI: 10.1073/pnas.60.3.921](https://doi.org/10.1073/pnas.60.3.921)

Schwarz G (1978) Estimating the dimension of a model. *Ann Stat* 6:461–464.
[DOI: 10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136)

Shen JC, Rideout WM, Jones PA (1994) The rate of hydrolytic deamination of
5-methylcytosine in double-stranded DNA. *Nucleic Acids Res* 22:972–976.
[DOI: 10.1093/nar/22.6.972](https://doi.org/10.1093/nar/22.6.972)

Shibutani S, Takeshita M, Grollman AP (1991) Insertion of specific bases during DNA
synthesis past the oxidation-damaged base 8-oxodG. *Nature* 349:431–434.
[DOI: 10.1038/349431a0](https://doi.org/10.1038/349431a0)

Steenken S, Jovanovic SV (1997) How easily oxidizable is DNA? One-electron reduction
potentials of adenosine and guanosine radicals in aqueous solution. *J Am Chem Soc*
119:617–618.
[DOI: 10.1021/ja962255b](https://doi.org/10.1021/ja962255b)

Wiebauer K, Jiricny J (1990) Mismatch-specific thymine DNA glycosylase and DNA polymerase
β mediate the correction of G·T mispairs in nuclear extracts from human cells.
*Proc Natl Acad Sci USA* 87:5842–5845.
[DOI: 10.1073/pnas.87.15.5842](https://doi.org/10.1073/pnas.87.15.5842)

Willerslev E, Cooper A (2005) Ancient DNA. *Proc Biol Sci* 272:3–16.
[DOI: 10.1098/rspb.2004.2813](https://doi.org/10.1098/rspb.2004.2813)

Orlando L, Gilbert MTP, Willerslev E (2015) Reconstructing ancient genomes and epigenomes.
*Nat Rev Genet* 16:395–408.
[DOI: 10.1038/nrg3935](https://doi.org/10.1038/nrg3935)

Pääbo S, Poinar H, Serre D, et al. (2004) Genetic analyses from ancient DNA. *Annu Rev
Genet* 38:645–679.
[DOI: 10.1146/annurev.genet.37.110801.143214](https://doi.org/10.1146/annurev.genet.37.110801.143214)
