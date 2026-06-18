---
type: libtaph-json-field
title: library_qc
tier: standard
estimand: library quality control metrics and artifact flags
stability: stable
emitted_by: profile_json.cpp
---

### `library_qc`

The `library_qc` block consolidates per-library quality diagnostics that are orthogonal to damage estimation: adapter remnant detection, terminal hexamer composition, position-0 artifact assessment, and downstream output-effect accounting. It is present in every profile regardless of library type or damage status; a `null` or absent block indicates that the read count fell below the minimum threshold required to compute terminal hexamer statistics. The block does not itself validate or invalidate damage estimates — instead it flags specific confounders and records which downstream outputs were corrected or remain suspect.

Nucleotide composition at read termini carries two signals simultaneously: genuine damage (C→T at 5′ for double-stranded libraries, G→A at 3′) and technical artifacts introduced by adapter ligation chemistry, random hexamer priming, and incomplete adapter trimming. Incomplete 5′ adapter removal leaves a residual prefix whose fixed sequence distorts the apparent terminal base composition and inflates apparent C→T rates. Random hexamer priming introduces a well-characterized context bias in the first 6 bases of sequenced reads: certain hexamer sequences are captured at much higher efficiency than their frequency in the template would predict, shifting both nucleotide frequencies and trinucleotide context distributions at the 5′ read start. The `library_qc` block quantifies both of these sources. The hexamer-shift test uses a G-statistic comparing the observed 5′ hexamer frequency distribution against the interior-read baseline; a significant shift (high `hex_shift_z`, low `hex_shift_p`) together with divergent entropy values (`hexamer_entropy_5prime` vs `hexamer_entropy_interior`) indicates hexamer priming bias rather than damage-driven composition change.

#### Top-level fields

| Field | Type | Description |
|-------|------|-------------|
| `adapter_offset_5prime` | integer | Number of bases identified as adapter remnant at the 5′ read start; 0 means no remnant detected, ≥1 means the first N bases are adapter-derived and position-0 base composition is unreliable |
| `adapter_offset_3prime` | integer | Analogous remnant offset at the 3′ end |
| `adapter_prefix_identified` | object\|null | `{"seq": "XXXXXX", "name": "..."}` — the 6-mer prefix detected and its library kit name if known, `null` when no recognisable prefix found |
| `adapter_stubs_clipped` | array of string | 6-mer stubs seen at the 5′ end above threshold |
| `adapter_stubs_clipped_3prime` | array of string | 6-mer stubs seen at the 3′ end above threshold |
| `adapter_stub_reads_checked` | integer\|null | Number of reads inspected for stub detection; `null` when the check was not run |
| `adapter_stub5_read_fraction` | float\|null | Fraction of reads carrying a detectable 5′ adapter stub |
| `adapter_stub3_read_fraction` | float\|null | Fraction of reads carrying a detectable 3′ adapter stub |
| `hexamer_entropy_5prime` | float | Shannon entropy (bits) of the 5′ hexamer distribution; lower than interior entropy signals non-random priming |
| `hexamer_entropy_interior` | float | Shannon entropy of the interior-read hexamer distribution; used as the unbiased baseline |
| `hexamer_terminal_interior_jsd` | float | Jensen–Shannon divergence between 5′ and interior hexamer distributions; values above ~0.05 indicate meaningful terminal composition shift |
| `hex_shift_g` | float | G-statistic (likelihood-ratio test) for the 5′ vs interior hexamer comparison |
| `hex_shift_z` | float | Z-score derived from `hex_shift_g`; large values indicate significant hexamer bias |
| `hex_shift_p` | float | P-value for the hexamer-shift G-test; typically 0.0 at sequencing depth because the test is highly powered |
| `hexamer_excess_tc` | float | Excess T+C fraction at terminal hexamers relative to interior; positive values are consistent with C→T deamination signal but also with primer bias |
| `position0_artifact_5prime` | boolean | `true` when the base at the first read position is dominated by adapter sequence rather than template; downstream damage estimates at position 0 are unreliable |
| `position0_artifact_3prime` | boolean | Analogous flag for the last read position |
| `short_read_frac` | float | Fraction of reads shorter than the minimum length used for interior hexamer estimation; high values reduce the reliability of interior-read baselines |
| `depurination_detected` | boolean | Mirrors the Channel E detection result in this QC context; `false` does not exclude depurination, only indicates the threshold for a positive call was not reached |
| `top_hexamers_5prime` | array of object | Top-5 over-represented hexamers at the 5′ terminus, each with `seq` (string), `log2fc` (float, fold-change over interior baseline), and `damage_consistent` (boolean, whether the hexamer context is compatible with C→T damage) |
| `top_hexamers_3prime` | array of object\|absent | Same structure as `top_hexamers_5prime` but for 3′ terminus; absent in single-stranded original-strand libraries where 3′ composition is structurally asymmetric |
| `flags` | array of string | Machine-readable QC flag list; empty array means no issues detected. Known values: `adapter_remnant_5prime` — unclipped adapter bases confirmed at 5′; `hexamer_terminal_shift` — significant hexamer distribution shift; `hexamer_artifact_bias` — terminal hexamer composition dominated by priming artifact; `ds_3prime_signal_absent` — expected 3′ deamination signal not found in a double-stranded library context |
| `output_effects` | object | Maps downstream output names (e.g., `d_max_3`) to `{"evidence": [...], "status": "corrected"\|"confounded"}` records; `corrected` means the adapter offset was applied and the estimate is usable; `confounded` means the output cannot be trusted |
| `diagnostic_groups` | object | Structured sub-blocks grouping related diagnostics; see below |

#### `diagnostic_groups` sub-blocks

| Sub-block | Fields | Description |
|-----------|--------|-------------|
| `adapter_position_effects` | `adapter_offset_5prime`, `adapter_offset_3prime`, `position0_artifact_5prime`, `position0_artifact_3prime`, `prefix_read_fraction_5prime`, `prefix_read_fraction_3prime`, `corrected_outputs`, `confounded_outputs`, `suspect_outputs`, `residual_outputs`, `selection_biased_outputs` | Detailed breakdown of which outputs were corrected vs remain suspect due to adapter remnants; `residual_outputs` lists signals that survive the correction but retain minor artifact influence |
| `end_damage_profile` | `d5_max_rate_pos0_4`, `d5_amp_over_bg`, `d5_var_pos0_9`, `d5_profile_flat`, `d5_blunting_suspected`, `d3_max_rate_pos0_4`, `d3_amp_over_bg`, `d3_var_pos0_9`, `d3_profile_flat` | Summary statistics of the C→T damage profile shape at each end; `amp_over_bg` is the damage amplitude above the interior background rate; `profile_flat` flags libraries where no position-dependent decay is visible; `d5_blunting_suspected` flags libraries where damage at position 0 is lower than at positions 1–2, consistent with enzymatic end-repair |
| `end_hexamer_asymmetry` | `hexamer_end_fwd_excess_jsd`, `hexamer_end_rc_excess_jsd`, `hexamer_end_rc_excess_jsd_status`, `rc_overlap_topk`, `terminal_excess_mass_5prime`, `terminal_excess_mass_3prime`, `terminal_hexamer_n_5prime`, `terminal_hexamer_n_3prime`, `topk` | Measures asymmetry between forward and reverse-complement hexamer distributions at each end; `hexamer_end_rc_excess_jsd_status` of `ok` means the reverse-strand hexamer excess is within expected range; `terminal_excess_mass` is the fraction of the terminal distribution that exceeds the interior baseline |
| `ss_end_asymmetry` | `d5_suppressed`, `d3_selection_biased`, `no_reliable_estimate`, `recommended_estimate`, `note` | Single-stranded library-specific: flags when 5′ damage is suppressed relative to expectation, when 3′ rates are biased by strand-capture selection, and provides a `recommended_estimate` directive (`none` when both ends are reliable, otherwise names the preferred end) |
| `terminal_hexamer_bias` | `hexamer_excess_tc`, `hexamer_terminal_interior_jsd`, `top_damage_consistent_fraction_5prime`, `top_damage_consistent_fraction_3prime`, `residual_outputs` | Quantifies the fraction of the top-5 over-represented hexamers that are damage-consistent; a low `top_damage_consistent_fraction` alongside a high `hexamer_terminal_interior_jsd` indicates the terminal shift is primarily priming artifact rather than damage |

#### Interpretation

The `flags` array is the fastest triage signal: an empty array means the QC passed without triggering any threshold. When `adapter_remnant_5prime` or `hexamer_artifact_bias` is present, inspect `output_effects` to determine whether the affected downstream estimates (typically `d_max_5`) were corrected or remain confounded before using them. A high `hexamer_terminal_interior_jsd` (above ~0.05–0.10) paired with a low `top_damage_consistent_fraction_5prime` strongly suggests that apparent 5′ enrichment is priming-driven rather than deamination-driven. Conversely, a moderate JSD with a high `top_damage_consistent_fraction` supports genuine damage. The `end_damage_profile.d5_blunting_suspected` field is particularly informative for assessing whether end-repair has attenuated position-0 damage and caused `d_max` to be underestimated.