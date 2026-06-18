---
type: libtaph-field-index
title: libtaph JSON field reference
---

# libtaph JSON field reference

Each file documents one JSON block emitted by `fqdup profile --json`.
YAML frontmatter encodes `tier`, `estimand`, `stability` — consumed by `fqdup profile --json-schema`.

| File | Title | Tier | Estimand |
|------|-------|------|----------|
| [top_level_fields.md](./top_level_fields.md) | `top_level_fields` | `summary` | scalar top-level fields: d_max, library_type, n_reads, schema_version |
| [bulk_damage.md](./bulk_damage.md) | `bulk_damage` | `standard` | reference-free bulk deamination rate (Design-B estimator) |
| [burial_fingerprint.md](./burial_fingerprint.md) | `burial_fingerprint` | `standard` | horizon-specific terminal C→T context fingerprint for stratigraphic identification |
| [channels_fgh.md](./channels_fgh.md) | `channels_fgh` | `standard` | oxidative terminal enrichment channels F (C→A), G (C→G), H (A→T) |
| [complement_asymmetry.md](./complement_asymmetry.md) | `complement_asymmetry` | `standard` | G-to-T vs C-to-A strand-balance asymmetry damage proxies |
| [damage_channel_stats.md](./damage_channel_stats.md) | `damage_channel_stats` | `standard` | 12-channel substitution terminal excess, interior fraction, decay lambda |
| [deam_context_spectrum.md](./deam_context_spectrum.md) | `deam_context_spectrum` | `standard` | 16-channel trinucleotide N[C→T]N excess at 5-prime and 3-prime termini |
| [deam_stratified_channels.md](./deam_stratified_channels.md) | `deam_stratified_channels` | `standard` | per-deamination-stratum methylation and depurination channel summaries |
| [deamination.md](./deamination.md) | `deamination` | `standard` | deamination damage profile across positions and channels |
| [depurination.md](./depurination.md) | `depurination` | `standard` | depurination (AP-site) damage metrics |
| [depurination_deconvolution.md](./depurination_deconvolution.md) | `depurination_deconvolution` | `standard` | AP-site C→T pathway signal deconvolved from bulk deamination |
| [end_motif_enrichment.md](./end_motif_enrichment.md) | `end_motif_enrichment` | `standard` | terminal nucleotide over-representation from fragmentation sequence preference |
| [fragmentation.md](./fragmentation.md) | `fragmentation` | `standard` | read-length distribution and damage-length coupling (fragmentation proxy) |
| [interior_ct_cluster.md](./interior_ct_cluster.md) | `interior_ct_cluster` | `standard` | interior C-to-T substitution cluster analysis |
| [library_qc.md](./library_qc.md) | `library_qc` | `standard` | library quality control metrics and artifact flags |
| [oxidation.md](./oxidation.md) | `oxidation` | `standard` | combined oxidation object (GC-depletion, epsilon, scission, strand asymmetry) |
| [oxidation_epsilon.md](./oxidation_epsilon.md) | `oxidation_epsilon` | `standard` | epsilon parameter: excess G→T magnitude from oxidation model |
| [oxidation_like.md](./oxidation_like.md) | `oxidation_like` | `standard` | boolean flag: damage pattern matches oxidative profile |
| [oxidative_scission.md](./oxidative_scission.md) | `oxidative_scission` | `standard` | oxidative strand scission from C→G and A→T terminal channel enrichment |
| [oxo_two_marker.md](./oxo_two_marker.md) | `oxo_two_marker` | `standard` | dual-marker 5-prime G→T and 3-prime C→A OxoG strand-balance verification |
| [oxog_estimate.md](./oxog_estimate.md) | `oxog_estimate` | `standard` | 8-oxoguanine per-context oxidative G→T terminal excess rate |
| [per_read_overdispersion.md](./per_read_overdispersion.md) | `per_read_overdispersion` | `standard` | overdispersion of damage across reads; tests for heterogeneous ancient mixture |
| [preservation.md](./preservation.md) | `preservation` | `standard` | composite taphonomic DNA preservation metrics |
| [scission.md](./scission.md) | `scission` | `standard` | strand scission rate from exponential fit to fragment-length right tail |
| [tau_discriminator.md](./tau_discriminator.md) | `tau_discriminator` | `standard` | discriminant score separating hydrolytic from oxidative damage regimes |
| [bic.md](./bic.md) | `bic` | `diagnostics` | nested model comparison BIC metrics (75 sub-fields) |
| [tetranuc_damage_rates.md](./tetranuc_damage_rates.md) | `tetranuc_damage_rates` | `diagnostics` | per-tetranucleotide context damage rates (900+ sub-fields) |
| [trinuc_spectrum_by_deam.md](./trinuc_spectrum_by_deam.md) | `trinuc_spectrum_by_deam` | `diagnostics` | trinucleotide spectra stratified by deamination (large nested arrays) |

## Tier meanings

| Tier | Emitted by default | Use |
|------|--------------------|-----|
| `summary` | yes | automated pipelines, quick status |
| `standard` | yes | full analysis |
| `diagnostics` | `--json-level full` only | debugging, method development |
| `deprecated` | yes, with warning | legacy; removed in v3 |
