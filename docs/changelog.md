# Changelog

## main (unreleased)

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
- 78-sample SE regression (33 clay91 + 45 ellesmere; 50 ds + 28 ss truth):
  77/78 correct in v9 vs 63/78 in v8; all 14 flips in the correct direction
  (no regressions). Single remaining miss is a real ambiguity
  (`LV7009022725-20S-IS160`: ga3=0.030, ga0=0.047, margin 1.1).
- Dataset 1 (46 DS + 45 SS): 88/91 correct, 3 UNKNOWN (100% on determined calls)
- Dataset 2 (78 DS + 146 SS): 193/196 correct on determined calls, 28 UNKNOWN (98.5%)

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
