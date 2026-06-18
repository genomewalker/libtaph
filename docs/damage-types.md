# Damage types and channels

libtaph scans raw FASTQ reads for reference-free damage observables and empirical
endpoint/context proxies, then writes a structured JSON report. Damage values are
fractions in [0, 1]; multiply by 100 for percentages.

Ancient DNA degradation falls into three broad categories: **hydrolytic damage**
(cytosine deamination; depurination/AP-site chemistry), **oxidative damage**
(8-oxoguanine formation), and **library-preparation artefacts** that mimic or confound
these signals. Some channels target established lesions; others are empirical
reference-free observables that require cautious mechanism interpretation.

---

## JSON output reference

Full per-field documentation lives in [`docs/fields/`](./fields/index.md).
Each file has YAML frontmatter (`tier`, `estimand`, `stability`) and the full prose description.

<!-- BEGIN FIELD TABLE -->
| Block | File | Tier | Description |
|-------|------|------|-------------|
| `top_level_fields` | [top_level_fields.md](./fields/top_level_fields.md) | summary | scalar top-level fields: d_max, library_type, n_reads, schema_version |
| `bulk_damage` | [bulk_damage.md](./fields/bulk_damage.md) | standard | reference-free bulk deamination rate (Design-B estimator) |
| `burial_fingerprint` | [burial_fingerprint.md](./fields/burial_fingerprint.md) | standard | horizon-specific terminal C→T context fingerprint for stratigraphic identification |
| `channels_fgh` | [channels_fgh.md](./fields/channels_fgh.md) | standard | oxidative terminal enrichment channels F (C→A), G (C→G), H (A→T) |
| `complement_asymmetry` | [complement_asymmetry.md](./fields/complement_asymmetry.md) | standard | G-to-T vs C-to-A strand-balance asymmetry damage proxies |
| `damage_channel_stats` | [damage_channel_stats.md](./fields/damage_channel_stats.md) | standard | 12-channel substitution terminal excess, interior fraction, decay lambda |
| `deam_context_spectrum` | [deam_context_spectrum.md](./fields/deam_context_spectrum.md) | standard | 16-channel trinucleotide N[C→T]N excess at 5-prime and 3-prime termini |
| `deam_stratified_channels` | [deam_stratified_channels.md](./fields/deam_stratified_channels.md) | standard | per-deamination-stratum methylation and depurination channel summaries |
| `deamination` | [deamination.md](./fields/deamination.md) | standard | deamination damage profile across positions and channels |
| `depurination` | [depurination.md](./fields/depurination.md) | standard | depurination (AP-site) damage metrics |
| `depurination_deconvolution` | [depurination_deconvolution.md](./fields/depurination_deconvolution.md) | standard | AP-site C→T pathway signal deconvolved from bulk deamination |
| `end_motif_enrichment` | [end_motif_enrichment.md](./fields/end_motif_enrichment.md) | standard | terminal nucleotide over-representation from fragmentation sequence preference |
| `fragmentation` | [fragmentation.md](./fields/fragmentation.md) | standard | read-length distribution and damage-length coupling (fragmentation proxy) |
| `interior_ct_cluster` | [interior_ct_cluster.md](./fields/interior_ct_cluster.md) | standard | interior C-to-T substitution cluster analysis |
| `library_qc` | [library_qc.md](./fields/library_qc.md) | standard | library quality control metrics and artifact flags |
| `oxidation` | [oxidation.md](./fields/oxidation.md) | standard | combined oxidation object (GC-depletion, epsilon, scission, strand asymmetry) |
| `oxidation_epsilon` | [oxidation_epsilon.md](./fields/oxidation_epsilon.md) | standard | epsilon parameter: excess G→T magnitude from oxidation model |
| `oxidation_like` | [oxidation_like.md](./fields/oxidation_like.md) | standard | boolean flag: damage pattern matches oxidative profile |
| `oxidative_scission` | [oxidative_scission.md](./fields/oxidative_scission.md) | standard | oxidative strand scission from C→G and A→T terminal channel enrichment |
| `oxo_two_marker` | [oxo_two_marker.md](./fields/oxo_two_marker.md) | standard | dual-marker 5-prime G→T and 3-prime C→A OxoG strand-balance verification |
| `oxog_estimate` | [oxog_estimate.md](./fields/oxog_estimate.md) | standard | 8-oxoguanine per-context oxidative G→T terminal excess rate |
| `per_read_overdispersion` | [per_read_overdispersion.md](./fields/per_read_overdispersion.md) | standard | overdispersion of damage across reads; tests for heterogeneous ancient mixture |
| `preservation` | [preservation.md](./fields/preservation.md) | standard | composite taphonomic DNA preservation metrics |
| `scission` | [scission.md](./fields/scission.md) | standard | strand scission rate from exponential fit to fragment-length right tail |
| `tau_discriminator` | [tau_discriminator.md](./fields/tau_discriminator.md) | standard | discriminant score separating hydrolytic from oxidative damage regimes |
| `bic` | [bic.md](./fields/bic.md) | diagnostics | nested model comparison BIC metrics (75 sub-fields) |
| `tetranuc_damage_rates` | [tetranuc_damage_rates.md](./fields/tetranuc_damage_rates.md) | diagnostics | per-tetranucleotide context damage rates (900+ sub-fields) |
| `trinuc_spectrum_by_deam` | [trinuc_spectrum_by_deam.md](./fields/trinuc_spectrum_by_deam.md) | diagnostics | trinucleotide spectra stratified by deamination (large nested arrays) |
<!-- END FIELD TABLE -->


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
| C | G→T stop codon uniformity | 8-oxoG oxidation | `ox_gt_uniformity`, `ox_gt_asymmetry`, `s_oxog_16ctx` |
| D | G→T and C→A transversion rates | 8-oxoG oxidation | `complement_asymmetry` |
| E | Terminal A+G enrichment over all bases | AP-site/depurination proxy | `depurination` |
| Length | Read-length distribution and damage-length coupling | Fragmentation/selection proxy | `fragmentation` |
| F | C→A terminal enrichment (bottom-strand 8-oxoG) | 8-oxoG complement | `complement_asymmetry` |
| G | C→G terminal enrichment | Empirical oxidative-context proxy | `complement_asymmetry` |
| H | A→T terminal enrichment | Empirical; mechanism uncertain | `complement_asymmetry` |
| CpG split | C→T amplitude by CpG / non-CpG context | Methylation-enhanced deamination | `cpg_like` |
| Interior clustering | Adjacent CT co-occurrence in read interior | Clustered interior deamination | `interior_ct_cluster` |
| 8-oxoG 16-ctx | G→T asymmetry by trinucleotide context | 8-oxoG context specificity | `s_oxog_16ctx` |
| Context spectrum | 16-channel trinucleotide C→T excess | Context-dependent deamination / CpG methylation | `deam_context_spectrum` |
| Deam strata | Methylation and depurination stratified by damage level | Ancient-fraction isolation | `deam_stratified_channels` |
| Trinuc × strata | Raw trinucleotide context counts per deamination stratum | Stratum-aware context decomposition | `trinuc_spectrum_by_deam` |
| Channel stats | Per-channel terminal excess, interior fraction, decay λ | Unified substitution-channel summary | `damage_channel_stats` |
| 4-mer CT | 4-mer N[C→T]NN terminal excess + deamination strata | CHG/CHH/CpG separation; ancient C→T rates | `tetranuc_damage_rates.ct_5prime[_by_deam]` |
| 4-mer GT | 4-mer N[G→T]NN terminal excess + deamination strata | 8-oxoG context specificity; stratum monotonicity test | `tetranuc_damage_rates.gt_5prime[_by_deam]` |

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
  modern components. Identifiability requires δ-separated, weight-bearing classes; fits that
  fail this gate report `status: undetermined` and no longer feed `d_max_combined`.
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

Henderson PT, Delaney JC, Gu F, Tannenbaum SR, Essigmann JM (2002) Oxidation of
7,8-dihydro-8-oxoguanine affords lesions that are potent sources of replication errors
in vivo. *Biochemistry* 41:914–921.
[DOI: 10.1021/bi0156355](https://doi.org/10.1021/bi0156355)

Gansauge M-T, Meyer M (2013) Single-stranded DNA library preparation for the sequencing
of ancient or damaged DNA. *Nat Protoc* 8:737–748.
[DOI: 10.1038/nprot.2013.038](https://doi.org/10.1038/nprot.2013.038)

Gansauge M-T, Gerber T, Glocke I, et al. (2017) Single-stranded DNA library preparation
from highly degraded DNA using T4 DNA ligase. *Nucleic Acids Res* 45:e79.
[DOI: 10.1093/nar/gkx033](https://doi.org/10.1093/nar/gkx033)

Kamiya H, Miura H, Murata-Kamiya N, et al. (1995) 8-Hydroxyadenine (7,8-dihydro-8-oxoadenine)
induces misincorporation in in vitro DNA synthesis and mutations in NIH 3T3 cells.
*Nucleic Acids Res* 23:2893–2899.
[DOI: 10.1093/nar/23.15.2893](https://doi.org/10.1093/nar/23.15.2893)

Jónsson H, Ginolhac A, Schubert M, Johnson PLF, Orlando L (2013) mapDamage2.0: fast
approximate Bayesian estimates of ancient DNA damage parameters. *Bioinformatics*
29:1682–1684.
[DOI: 10.1093/bioinformatics/btt193](https://doi.org/10.1093/bioinformatics/btt193)

Kapp JD, Green RE, Shapiro B (2021) A fast and efficient single-stranded genomic library
preparation method optimized for ancient DNA. *J Hered* 112:241–249.
[DOI: 10.1093/jhered/esab012](https://doi.org/10.1093/jhered/esab012)

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

Neeley WL, Essigmann JM (2006) Mechanisms of formation, genotoxicity, and mutation of
guanine oxidation products. *Chem Res Toxicol* 19:491–505.
[DOI: 10.1021/tx0600043](https://doi.org/10.1021/tx0600043)

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
