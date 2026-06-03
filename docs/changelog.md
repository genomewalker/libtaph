# Changelog

## main (unreleased)

### Breaking: `cpg_contrast` removed from C++ API and JSON output

The upstream-context contrast (keyed on `decoded[p-1]`, base **5' of** the
cytosine) measured the GpC dinucleotide, not CpG. GpC methylation is absent
in eukaryotes; the field implied chemistry that is not there. The correct
CpG/5mC axis is `log2_cpg_ratio` / `cpg_z` (downstream `cpg_like` block,
keyed on `decoded[p+1]`).

Removed from `SampleDamageProfile`, `LengthBinDamageProfile`, and JSON output.
Update consumers to use `log2_cpg_ratio` / `cpg_z`.

`dipyr_contrast` (CC/TC vs AC/GC, correct CPD dipyrimidine geometry) is
unchanged and its direction has been independently audited.

### Ancient-fraction deamination — soft-EM + independent fitted profiles

#### Added
- `damaged_fraction_valid`, `damaged_fraction_pi`, `damaged_fraction_n`,
  `damaged_fraction_d5`, `damaged_fraction_d3`: soft-EM ancient-fraction
  estimate populated by fqdup's per-read LLR scorer (fused into the oxoG pass).
  `pi` is the posterior-weighted fraction of reads classified as ancient
  (soft-EM E-step), avoiding the selection-on-outcome bias of hard LLR > 0
  thresholding. `d5`/`d3` are derived via mixture identity (`bulk_d / pi`).
- `damaged_fraction_d5_fit` / `damaged_fraction_lambda5` /
  `damaged_fraction_d3_fit` / `damaged_fraction_lambda3`: independently fitted
  exponential decay (log-linear regression over all 15 positions) for the
  ancient-classified read pool. Gives the true ancient d_max + lambda without
  mixture-identity dilution.
- `modern_fraction_d5` / `modern_fraction_d3`: pos-0 peak deamination estimate
  for reads classified as modern.
- `modern_fraction_d5_fit` / `modern_fraction_lambda5` /
  `modern_fraction_d3_fit` / `modern_fraction_lambda3`: fitted exponential for
  the modern-classified pool (typically ≈ 0 for clean libraries).
- All fields serialised in the `ancient_fraction` JSON block under `ancient`
  and `modern` sub-objects.

#### Fixed
- `skip_pos0_5prime` on `LsdClassifyParams`: when `position_0_artifact_5prime`
  is set the per-read 5' LLR scorer starts from position 1, avoiding the
  adapter-blunting depletion at pos 0 from artificially lowering LLR scores
  on genuinely ancient reads.
- **3′ pos-0 artifact detection — DS blunting gap**: `position_0_artifact_3prime`
  previously only fired on the SS ligation spike (pos-0 G→A elevated above
  pos-1–4 mean by > 5 pp). DS libraries with adapter blunting show the
  opposite: pos-0 depleted, pos-1 elevated. Added complementary gap detection
  (`stats_3_pos0 < −0.005 && stats_3_pos1 > 0.01`, or a large pos-0→pos-1
  jump > 3 pp when pos-1 excess > 2%). Using pos-1 alone rather than the
  pos-1–4 average is critical: the four-position average dilutes the pos-1
  signal (shift ~0.008 vs ~0.060 at pos-1 alone), placing it below the
  detection threshold. Fixes `d_max_3prime = 0` in DS libraries affected by
  adapter blunting.

### Channels F/G/H — oxidative terminal enrichment

#### Added
- Channels F (C→A), G (C→G), and H (A→T): binomial z-score comparing terminal
  trinucleotide-context rate (positions 0–4, adapter-adjusted) to a far-interior
  reference (pos ≥ 30, ≥ 50 counts; mid-read fallback pos 5–14).
- `channel_h_z_p2plus`: secondary H score restricted to positions 2–4; robust
  when TC hexamer bias is present without a formal position-0 artifact.
- `ca_pre_interior` / `ca_stop_interior`: pre-aggregated uint64_t counters
  accumulated at read-scan time, giving Channel F the same symmetric interior
  reference as Channels G and H.

#### Fixed
- Channel H adapter skip now uses `p0_h5 = position_0_artifact_5prime` only.
  The shared `p0_tc_5` gate (which also fires on `hexamer_excess_tc < −0.02`)
  was incorrect for Channel H: TC hexamer bias is C→T-specific and does not
  affect A/T-context trinucleotides. Under strong TC bias the old gate silently
  excluded p=0 from Channel H, inflating z by ~3×.
- Renamed `p0_5` → `p0_tc_5`; comment documents scope to Channels F and G only.
- Renamed `ca_stop_baseline_3prime` → `ca_stop_rate_baseline_3prime`.
- Fixed 17 `uint64_t` profile fields initialised with float literal `0.0` → `0`.

### Added
- Native paired-end mode: `update_sample_profile_pe` consumes R1/R2 directly,
  mapping `R2[i]` complement to top-strand 3'-end position `i` (no SE-on-each-mate
  bleed-through). Insert length `M < L` detected by R1/R2 overlap (15 bp window,
  ≤3 mismatches); short-insert pairs are skipped via `pe_short_insert_skipped`
  counter, preventing adapter bases from polluting the C→T / G→A signal.
- Chemistry-aware Briggs fit: `d(i) = d_max·exp(-i/λ) + bg` with tail-anchored
  `bg` from positions 20..49 (interior background, immune to deamination).
  Replaces global-mean `bg` which over-estimated bg in damaged libraries and
  collapsed `d_max`.
- Area-excess statistic: per-channel sum of `(rate - bg)` over the first 10
  positions, plus a likelihood-ratio score against the bg-only null.
  Decouples classification from a single `d_max` point estimate.
- Protocol-tag aware library interpretation; chemistry priors flow through
  `library_interpretation.cpp`.
- `UNKNOWN` library-type category when no model beats M_bias (zero-damage libraries)
- GA0 bilateral rescue: `spike_is_ss && d5 ≤ 0.005 → SS`, validated on 24 DS controls (all d5 ≥ 0.11)
- Channel B₃′: G→A stop codon conversion at 3' end for SS library damage validation

### Fixed
- LV7001884491 PE DS-library misclassified as `SS_orig` in v8: caused by the
  6-base palindrome guard, which only triggered for insert `M=6` and let
  short-insert pairs (≈83% of aDNA pairs at L=101) bleed adapter bases into
  R2[0..14]. After complement-mapping, the adapter's specific A:G ratio
  produced `a_freq_3prime[0]=0.396` (depleted vs interior 0.466) which the
  classifier read as `M_SS_orig`. Replaced with proper R1/R2 overlap detection.
- Bug: missing `best = bic_M_DS_spike` update in spike_is_ss branch of waterfall
- Bug: `best_ss` included M_SS_orig without applying the `ct3.ΔBIC > 0` gate
- Bug: M_DS_spike rescue could fire when M_DS_symm_art won via joint fit with marginal ct5/ga3 ΔBIC ≤ 0, now restricted to `ds_spike_won`

### Validation
- ~80-sample SE regression across heavily to lightly degraded permafrost and sediment libraries:
  77/78 correct in v9 vs 63/78 in v8; all 14 flips in the correct direction (no regressions).
- Dataset 1 (~90 samples, DS+SS mix): 88/91 correct, 3 UNKNOWN (100% on determined calls)
- Dataset 2 (~225 samples, DS+SS mix): 193/196 correct on determined calls, 28 UNKNOWN (98.5%)

---

## Previous milestones

### BIC library-type classifier (initial)
- 4-channel joint BIC model: ct5, ga3, ga0, ct3
- 7 composite BIC models: M_bias, M_DS_symm, M_DS_spike, M_DS_symm_art, M_SS_comp, M_SS_orig, M_SS_asym
- `spike_is_ss` gate: ga0.amplitude ≥ 0.10 routes M_DS_spike to SS model set
- Post-hoc symmetry check: ct5.ΔBIC / ga3.ΔBIC < 0.50 with ga3.ΔBIC > 30,000 → SS
- M_DS_spike rescue using `max_damage_5prime` ≤ 0.005 (d5 bilateral gate)

### Initial release
- Briggs-model exponential decay fitting (WLS, fixed λ)
- GC-stratified damage estimation (10 bins, EM mixture)
- Channel A/B cross-validation
- Channel C (8-oxoG uniformity) and Channel D (G→T / C→A)
- Hexamer-based damage LLR
- Alignability-weighted D_max estimation
