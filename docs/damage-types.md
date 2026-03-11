# Damage types and channels

libdart-damage scans raw FASTQ reads for six biochemical damage channels and writes a
structured JSON report. Damage values are fractions in [0, 1]; multiply by 100 for
percentages. The report covers cytosine deamination at both termini, library-type
classification, 8-oxoG oxidative damage, interior co-deamination clustering, and AP-site
fragmentation.

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

Cytosine in single-stranded overhangs spontaneously deaminates to uracil, which the
polymerase reads as thymine. The rate is highest at the fragment terminus and decays
exponentially inward: P(p) = d_max × exp(−λ × p) + bg. This block reports the fitted
parameters and the best single-number damage estimate.

#### Damage estimates

| Field | Description |
|-------|-------------|
| `d_max_5prime` | C→T excess at the 5′ terminus above the interior background |
| `d_max_3prime` | G→A excess at the 3′ terminus (DS) or C→T excess (SS original-strand) above background |
| `d_max_combined` | Best single damage estimate for this library; `source` records how it was derived |
| `d_metamatch` | GC-composition-weighted estimate comparable to metaDMG's D_max |
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
| `validated` | `true` — Channel B corroborates Channel A deamination |
| `artifact` | `true` — Channel A fires but Channel B contradicts it; likely a GC-composition artefact rather than genuine damage |

#### Per-position rates

| Field | Description |
|-------|-------------|
| `per_pos_5prime_ct[15]` | Raw T/(T+C) at 5′ positions 0–14 |
| `per_pos_3prime[15]` | A/(A+G) at 3′ positions 0–14 for DS; T/(T+C) for SS original-strand |

Values of `-1` indicate insufficient coverage at that position.

#### `cpg_like` sub-block — methylation signal

Methylated cytosines (5mC) deaminate faster than unmethylated ones at CpG dinucleotides.
This sub-block fits the 5′ C→T amplitude separately for CpG-context positions (next base
= G) and non-CpG positions, using the read interior as a composition baseline. A
`cpg_ratio` above 1 indicates methylation-enhanced deamination — the expected pattern for
organisms with CpG methylation. Values near 1 indicate unmethylated CpGs, modern
contamination, or a signal too eroded to detect context dependence.

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

Library preparation determines which damage channels are active, so the same BIC model
selection that drives classification also reveals which sub-channels contributed signal.
This block exposes the classifier's evidence for inspecting borderline calls.

The classifier fits four sub-channel amplitudes simultaneously:
- **ct5** — 5′ C→T exponential decay
- **ga3** — 3′ G→A exponential decay
- **ga0** — 3′ position-0 G→A spike (ligation-site pattern in SS complement-orientation libraries)
- **ct3** — 3′ C→T decay (SS original-orientation only)

These are evaluated under seven biological models — a null (all channels flat), two DS
models (symmetric and with end-repair spike), a DS end-repair-only model, and three SS
models (complement, original, and mixed orientations). The model with the lowest BIC wins.

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

8-Oxoguanine (8-oxoG) forms when guanine is oxidised. The polymerase misreads 8-oxoG as
adenine, producing G→T transversions. Unlike cytosine deamination, oxidation accumulates
throughout the molecule rather than at the termini, so the diagnostic signal is strand
asymmetry — elevated G→T without a corresponding C→A increase on the complement strand.

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
| `s_gt` | G→T vs C→A strand asymmetry contrast |
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

Genuine ancient 8-oxoG shows broad enrichment across contexts with a slight GG-context
bias. A spike confined to one or two contexts is more consistent with a preparation
artefact.

---

### `interior_ct_cluster` block — interior co-deamination

When cytosines in single-stranded micro-domains are co-deaminated in the same hydrolytic
event, T residues at nearby positions co-occur more often than expected by chance. This
block tests for that excess in the read interior (middle third, reads ≥ 30 bp), using AG
co-occurrence as a composition control.

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
beyond what strand composition alone predicts. In most aDNA samples this signal is weak;
it becomes informative for very old material or samples with unusual interior damage such
as permafrost specimens.

---

### `depurination` block — AP-site fragmentation (Channel E)

Purines (A, G) are lost by hydrolytic base loss more readily than pyrimidines, especially
under acidic or warm burial conditions. The resulting apurinic (AP) sites are strand-break
precursors; when a strand breaks at an AP site, the newly exposed 5′ end is enriched for
purines. This signature provides evidence of ancient fragmentation independent of
deamination.

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
libraries ligate adapters to both strands and sequence both; single-stranded (SS)
libraries circularize and sequence only one strand, and which strand determines the
damage pattern.

| Library type | Active sub-channels | Pattern |
|---|---|---|
| DS | ct5 + ga3 | Symmetric exponential decay at both ends |
| DS + end-repair artifact | ct5 + ga3 + ga0 | DS pattern plus isolated 3′ pos-0 spike |
| SS complement-orientation | ga0 | Strong 3′ pos-0 spike; 5′ flat |
| SS original-orientation | ct5 + ct3 | C→T at both 5′ and 3′ ends; no G→A |
| SS mixed orientations | ct5 + ga0 | 5′ C→T decay plus 3′ pos-0 spike |
| UNKNOWN | — | No sub-channel beats the null model |

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
  δ_d (double-stranded background rate), following Briggs et al. 2007. R² assesses fit quality.
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
