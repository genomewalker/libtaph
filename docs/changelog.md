# Changelog

## main (unreleased)

### Added
- `UNKNOWN` library-type category when no model beats M_bias (zero-damage libraries)
- GA0 bilateral rescue: `spike_is_ss && d5 ≤ 0.005 → SS` — validated on 24 DS controls (all d5 ≥ 0.11)
- Channel B₃′: G→A stop codon conversion at 3' end for SS library damage validation

### Fixed
- Bug: missing `best = bic_M_DS_spike` update in spike_is_ss branch of waterfall
- Bug: `best_ss` included M_SS_orig without applying the `ct3.ΔBIC > 0` gate
- Bug: M_DS_spike rescue could fire when M_DS_symm_art won via joint fit with marginal ct5/ga3 ΔBIC ≤ 0 — now restricted to `ds_spike_won`

### Validation
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
