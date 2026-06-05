# API Reference

All public symbols live in the `taph` namespace. Include `<taph/frame_selector_decl.hpp>`.

---

## Stability tiers

`SampleDamageProfile` fields are grouped into three tiers. The tier of a field is
annotated in `sample_damage_profile.hpp` with a comment block
(`// --- Stable outputs ---`, `// --- Diagnostics ---`,
`// --- Experimental ---`).

### Stable

These fields are committed across minor versions. Their names, types, and
semantics will not change without a major version bump and a changelog entry.

| Field group | Fields |
|---|---|
| Channel A d_max | `d_max_5prime`, `d_max_3prime`, `d_max_combined`, `d_max_source_str()` |
| Channel A decay | `lambda_5prime`, `lambda_3prime` |
| Channel A area | `damage_5prime_area_excess`, `damage_3prime_area_excess`, `damage_5prime_lr`, `damage_3prime_lr` |
| Channel A asymmetry | `asymmetry` |
| Channel A background | `bg_5prime_anchored`, `bg_3prime_anchored` |
| Library type | `library_type`, `library_type_auto_detected`, `library_type_str()` |
| Library BIC | `library_bic_bias`, `library_bic_ds`, `library_bic_ss`, `library_bic_mix` |
| CpG axis | `log2_cpg_ratio`, `cpg_ratio` |
| Context contrasts | `dipyr_contrast`, `gpc_contrast` |
| Channel B headline | `channel_b_valid`, `stop_amplitude_5prime`, `d_max_from_channel_b` |
| Oxidative channels | `channel_f_z`, `channel_g_z`, `channel_h_z`, `channel_h_z_p2plus` |
| Ancient fraction | `damaged_fraction_valid`, `damaged_fraction_pi`, `damaged_fraction_d5_fit`, `damaged_fraction_d3_fit`, `damaged_fraction_lambda5`, `damaged_fraction_lambda3` |
| Validation flags | `damage_validated`, `damage_artifact`, `terminal_inversion`, `is_valid()`, `is_detection_unreliable()` |

### Diagnostics

Useful for debugging and model development but may gain or lose fields in minor
releases. Do not write production code that hard-depends on these names.

- Library BIC internals: `libtype_amp_ct5/ga3/ga0/ct3`, `libtype_dbic_*`, `library_p_ds/ss/bias/winner`, `library_bic_winner_model`, `library_bic_margin`
- Background coverage: `bg_n_positions_*`, `bg_denominator_*`, `briggs_pos0_masked_*`
- Channel B internals: `channel_b_slope`, `channel_b_weight`, `channel_b_quantifiable`, `channel_b_inverted`, `stop_decay_llr_*`, `stop_conversion_rate_baseline`, per-stop-type LLRs (`*_taa_*`, `*_tag_*`, `*_tga_*`)
- Oxidative channel details: `ca_stop_rate_*`, `cg_stop_rate_*`, `at_stop_rate_*`, `*_uniformity_ratio*`, `channel_f/g/h_valid`, `channel_f/g/h3_valid`
- Hexamer chemistry: `hexamer_excess_tc`, `hexamer_damage_llr`
- Fraction d_max (mixture-identity, biased at low π): `damaged_fraction_d5`, `damaged_fraction_d3`, `modern_fraction_*`
- Pattern flags: `composition_bias_*`, `inverted_pattern_*`, `position_0_artifact_*`
- Preservation: `preservation_score`, `preservation_label`, `preservation_evidence`, `preservation_reliability`

### Experimental

No stability guarantee. Schema, semantics, and field names may change between any
two releases. Do not depend on these in downstream tools without pinning the
libtaph version.

- **Raw upstream-context arrays** — `ct5_t_by_upstream[ctx][pos]`, `ct5_total_by_upstream[ctx][pos]`, `ct5_t_interior_by_upstream[ctx]`, `ct5_total_interior_by_upstream[ctx]`, `dmax_ct5_by_upstream[ctx]`, `baseline_ct5_by_upstream[ctx]`, `cov_ct5_*_by_upstream[ctx]`
- **Per-position trinucleotide counts** — `tri_5prime_pos[pos][ctx]`, `tri_3prime_pos[pos][ctx]`, `tri_5prime_terminal`, `tri_5prime_interior`, `tri_3prime_terminal`, `tri_3prime_interior`
- **16-context OxoG panel** — `oxog16_t`, `oxog16_a_rc`, `s_oxog_16ctx`, `cov_oxog_16ctx`
- **Per-hexamer-prefix F/G/H sums** — `ca_pre_terminal_by_pfx`, `ca_stop_terminal_by_pfx`, and all `*_by_pfx` variants
- **OxoTwoMarkerBins struct** and all `ox_two_marker_bins` fields
- **Codon-position arrays** — `codon_pos_t_rate_5prime`, `codon_pos_a_rate_3prime`, `codon_pos_*_count_*`
- **Length×GC joint mixture** (`LengthGcJointMixtureResult`) and `MultiDomainResult` domain ensemble fields

---

## FrameSelector

The main entry point. All methods are `static`.

### `compute_sample_profile`

```cpp
static SampleDamageProfile compute_sample_profile(
    const std::vector<std::string>& sequences);
```

One-shot profile computation from a vector of sequences. Internally calls `reset_sample_profile`, iterates `update_sample_profile`, then `finalize_sample_profile`.

**Parameters**

| Name | Type | Description |
|------|------|-------------|
| `sequences` | `const std::vector<std::string>&` | DNA sequences (ACGT). Reliable estimation requires enough read length and coverage to define both terminal (positions 0–14) and interior (positions 30 to L-30) regions. |

**Returns** A fully populated `SampleDamageProfile`. Call `profile.is_valid()` to check whether enough reads were processed (≥ 1000).

---

### `update_sample_profile`

```cpp
static void update_sample_profile(
    SampleDamageProfile& profile,
    const std::string& seq);
```

Accumulate one read into `profile`.

Safe parallel pattern: update separate `SampleDamageProfile` objects in parallel threads, then merge them on the main thread with `merge_sample_profiles`. Do not call this on the same profile object from multiple threads without external synchronization.

---

### `update_sample_profile_pe`

```cpp
static void update_sample_profile_pe(
    SampleDamageProfile& profile,
    const std::string& r1,
    const std::string& r2);
```

Accumulate one **read pair** into `profile`. Maps `R1[i]` to top-strand 5'-end
position `i` and the complement of `R2[i]` to top-strand 3'-end position `i`,
so a single fragment contributes once to ct5 and once to ga3 (no double-count).

Pairs whose insert length is shorter than the read length read into the
sequencing adapter and would imprint adapter composition onto the terminal
channels. Such pairs are detected by R1/R2 overlap (15 bp window, ≤3
mismatches) and skipped; the count is reported in `pe_short_insert_skipped`.

Sets `profile.input_mode = InputMode::PAIRED`. Same threading contract as
`update_sample_profile` — accumulate into per-thread profiles, then merge.

---

### `update_sample_profile_weighted`

```cpp
static void update_sample_profile_weighted(
    SampleDamageProfile& profile,
    const std::string& seq,
    float weight);
```

Weighted accumulation for alignability-weighted damage estimation.

> **Note:** Only accumulates Channel A (T/(T+C), A/(A+G)), interior baseline, codon-position, and CpG counts. Channel B/C/D/E codon counts, hexamers, and GC-bin stratified data are **not** accumulated. Profiles built with this method will have `channel_b_valid = false` and GC-stratified fields unpopulated after finalization.

---

### `finalize_sample_profile`

```cpp
static void finalize_sample_profile(SampleDamageProfile& profile);
```

Compute all derived statistics from accumulated counts: fits the exponential decay, estimates D_max, runs the BIC library-type classifier. Must be called once after all `update_*` calls.

---

### `merge_sample_profiles`

```cpp
static void merge_sample_profiles(
    SampleDamageProfile& dst,
    const SampleDamageProfile& src);
```

Merge raw counts from `src` into `dst`. Useful for parallel accumulation. Call `finalize_sample_profile(dst)` after merging.

---

### `reset_sample_profile`

```cpp
static void reset_sample_profile(SampleDamageProfile& profile);
```

Zero the main accumulators used by standard damage estimation and library classification.

> **Important:** `reset_sample_profile` does not currently clear every auxiliary diagnostic accumulator — including Channel B₃′/C/D/E, GC-stratified, and alignability-weighted fields. For guaranteed clean reuse, prefer constructing a fresh `SampleDamageProfile{}` (zero-initialized) rather than calling `reset_sample_profile` on a previously populated object.

---

## SampleDamageProfile

All fields are public. Key results after `finalize_sample_profile`:

### Primary damage estimates

| Field | Type | Description |
|-------|------|-------------|
| `d_max_5prime` | `float` | Calibrated D_max at 5' end: $A/(1-b)$ |
| `d_max_3prime` | `float` | Calibrated D_max at 3' end |
| `d_max_combined` | `float` | Final asymmetry-aware D_max estimate |
| `lambda_5prime` | `float` | Fitted decay constant $\lambda$ at 5' |
| `lambda_3prime` | `float` | Fitted decay constant $\lambda$ at 3' |
| `asymmetry` | `float` | `\|d5 - d3\| / mean(d5, d3)`, >0.5 flagged as suspicious |
| `d_max_source_str()` | `const char*` | Source used for `d_max_combined`: `"average"`, `"5prime_only"`, `"3prime_only"`, `"channel_b_structural"`, `"channel_b3_structural"`, `"min_asymmetry"`, `"max_ss_asymmetry"`, `"none"` |

### Chemistry-aware background and area-excess

| Field | Type | Description |
|-------|------|-------------|
| `bg_5prime_anchored` / `bg_3prime_anchored` | `float` | Tail-anchored background: trimmed mean of C→T (G→A) rate over positions 20..49, restricted to positions with denom ≥ 100. Used as `b` in the Briggs fit. |
| `bg_n_positions_5prime` / `bg_n_positions_3prime` | `int` | Number of tail positions that contributed to the anchored bg. |
| `bg_denominator_5prime` / `bg_denominator_3prime` | `double` | Total T+C (A+G) coverage in the bg window. |
| `briggs_pos0_masked_5prime` / `briggs_pos0_masked_3prime` | `bool` | True when position 0 was excluded from the Briggs fit because the protocol-tag fingerprint dominated it (ligation footprint). |
| `damage_5prime_area_excess` / `damage_3prime_area_excess` | `float` | $\sum_{i=k}^{14} \max(0, r_i - b)$ over the first 15 positions (k=1 if pos0 masked, else k=0). Robust headline statistic for cross-library comparison; preferred over `d_max` when chemistry-tag is set. |
| `damage_5prime_lr` / `damage_3prime_lr` | `float` | Log-likelihood ratio of the per-position binomial(rate, bg) model vs the bg-only null over positions k..14. Coverage-scaling companion to area-excess. |

### Input mode

| Field | Type | Description |
|-------|------|-------------|
| `input_mode` | `InputMode` | `SINGLE` (default) or `PAIRED`. Set to `PAIRED` by `update_sample_profile_pe`. |
| `pe_short_insert_skipped` | `uint64_t` | Number of pairs skipped because R1/R2 overlap revealed insert < read length (adapter read-through). |

### Library-type classification

| Field | Type | Description |
|-------|------|-------------|
| `library_type` | `LibraryType` | `DOUBLE_STRANDED`, `SINGLE_STRANDED`, or `UNKNOWN` |
| `library_type_auto_detected` | `bool` | `true` if set by classifier, `false` if user-forced |
| `library_type_str()` | `const char*` | `"double-stranded"`, `"single-stranded"`, `"unknown"` |

`LibraryType::UNKNOWN` means no model beat the null (M_bias): insufficient damage signal for a confident call. Treat as DS for deduplication unless metadata is available.

### BIC model scores (library-type classifier)

| Field | Type | Description |
|-------|------|-------------|
| `library_bic_bias` | `double` | BIC of M_bias (null model) |
| `library_bic_ds` | `double` | Best DS model BIC |
| `library_bic_ss` | `double` | Best SS model BIC |
| `library_bic_mix` | `double` | BIC of M_SS_full (4-channel unconstrained) |

`library_bic_ds - library_bic_ss` > 0 means SS is favoured. Values reach ~10⁹ at high coverage; stored as `double` to preserve precision.

### Per-channel classifier amplitudes and ΔBIC

| Field | Description |
|-------|-------------|
| `libtype_amp_ct5` / `libtype_dbic_ct5` | 5' C→T fitted amplitude and ΔBIC |
| `libtype_amp_ga3` / `libtype_dbic_ga3` | 3' G→A smooth decay (pos 1-10) |
| `libtype_amp_ga0` / `libtype_dbic_ga0` | 3' G→A pos-0 spike |
| `libtype_amp_ct3` / `libtype_dbic_ct3` | 3' C→T (SS original-orientation signal) |

ΔBIC > 0 means the alt (decay) model is preferred over the null for that channel.

### Per-position arrays

All arrays are 15 elements, indexed 0–14 from the read terminus.

| Field | Description |
|-------|-------------|
| `damage_rate_5prime[p]` | C→T excess rate at 5' position `p` |
| `damage_rate_3prime[p]` | G→A excess rate at 3' position `p` |
| `t_freq_5prime[p]` | **Before** `finalize_sample_profile`: raw T count. **After**: T/(T+C) ratio (normalized in-place). |
| `tc_total_5prime[p]` | T+C coverage at 5' position `p` (raw count; not normalized by finalization) |
| `a_freq_3prime[p]` | **Before** `finalize_sample_profile`: raw A count. **After**: A/(A+G) ratio (normalized in-place). |
| `ag_total_3prime[p]` | A+G coverage at 3' position `p` (raw count; not normalized by finalization) |

### GC-stratified mixture model

Available after `finalize_sample_profile`. Requires at least one GC bin with sufficient reads.

| Field | Type | Description |
|-------|------|-------------|
| `mixture_pi_ancient` | `float` | Fraction of C-sites in high-damage components |
| `mixture_d_ancient` | `float` | Expected damage rate among ancient reads (δ > 5%) |
| `mixture_d_population` | `float` | Population-average damage rate across all C-sites |
| `mixture_d_population_highgc` | `float` | High-GC-weighted mean of per-class δ (GC ≥ 50%); diagnostic only, not a reference comparator |
| `mixture_K` | `int` | Number of mixture components selected by BIC |
| `mixture_converged` | `bool` | Whether the EM algorithm converged |
| `gc_stratified_valid` | `bool` | At least one GC bin has a valid estimate |

### Ancient-fraction deamination

Populated by fqdup's per-read LLR scorer when the ancient-fraction pass is
fused into the oxoG worker (`damaged_fraction_valid == true`). Not set by
the standalone libtaph `FrameSelector` API.

| Field | Type | Description |
|-------|------|-------------|
| `damaged_fraction_valid` | `bool` | `true` when the fqdup LLR pass ran and produced a valid estimate |
| `damaged_fraction_pi` | `float` | Soft-EM posterior fraction of reads classified as ancient |
| `damaged_fraction_n` | `int64_t` | Hard-call count of reads with LLR > 0 (informational) |
| `damaged_fraction_d5` | `float` | 5' d_max for the ancient component via mixture identity: `bulk_d5 / pi` |
| `damaged_fraction_d3` | `float` | 3' d_max via mixture identity: `bulk_d3 / pi` |
| `damaged_fraction_d5_fit` | `float` | 5' d_max from independent log-linear regression over all 15 positions of ancient-classified reads |
| `damaged_fraction_lambda5` | `float` | 5' decay constant λ from the independent fit |
| `damaged_fraction_d3_fit` | `float` | 3' d_max from independent fit of ancient-classified reads |
| `damaged_fraction_lambda3` | `float` | 3' decay constant λ from the independent fit |
| `modern_fraction_d5` | `float` | 5' d_max at pos 0 for reads classified as modern |
| `modern_fraction_d3` | `float` | 3' d_max at pos 0 for reads classified as modern |
| `modern_fraction_d5_fit` | `float` | 5' d_max from independent fit of modern-classified reads |
| `modern_fraction_lambda5` | `float` | 5' λ from modern fit |
| `modern_fraction_d3_fit` | `float` | 3' d_max from independent fit of modern-classified reads |
| `modern_fraction_lambda3` | `float` | 3' λ from modern fit |

The mixture-identity fields (`damaged_fraction_d5/d3`) equal `bulk_d / pi`,
valid only under `d_mod ≈ 0` (modern reads carry no deamination) and equal
C-site density between pools. When `d_mod > 0`, the estimate is upward-biased —
the bias is directional, not just noisy — and grows as `pi` decreases. The
fitted fields (`*_fit`) measure the exponential shape directly in the classified
read pools and avoid both assumptions, at the cost of requiring enough ancient
reads for a reliable fit. For high-endogenous samples both approaches agree; for
low-endogenous samples prefer the fitted values. A `damaged_fraction_valid ==
false` result means the LLR pass did not run (paired-mode input, bulk d_max
below the 0.01 gate, or insufficient reads).

Note on `d_max_5prime` (`A/(1−b)`): here `b` is the interior background
deamination **rate** (not a read fraction), estimated from positions 20–49.
The formula gives the model amplitude parameter — the exponential component's
value extrapolated to position 0 — and requires `b < 1`.

---

### Context-aware 5' C→T (upstream base)

Per-context accumulators and shrinkage-fitted amplitudes, indexed by `CTX_AC = 0`, `CTX_CC = 1`, `CTX_GC = 2`, `CTX_TC = 3` (constants under `SampleDamageProfile::UpstreamContext`, with `N_UPSTREAM_CTX = 4`).

| Field | Type | Description |
|-------|------|-------------|
| `ct5_t_by_upstream[ctx][pos]` / `ct5_total_by_upstream[ctx][pos]` | `double[4][N_POS]` | T and T+C counts at 5' position `pos` split by upstream base. Raw before `finalize_sample_profile`. |
| `ct5_t_interior_by_upstream[ctx]` / `ct5_total_interior_by_upstream[ctx]` | `double[4]` | Interior counts per upstream context, used as per-context baseline. |
| `dmax_ct5_by_upstream[ctx]` | `float[4]` | Fitted exponential amplitude per context. `NaN` if insufficient coverage. |
| `baseline_ct5_by_upstream[ctx]` | `float[4]` | Interior baseline rate per context. |
| `cov_ct5_terminal_by_upstream[ctx]` / `cov_ct5_interior_by_upstream[ctx]` | `float[4]` | Coverage used in the per-context fit. |
| `dipyr_contrast` | `float` | `½(d_CC + d_TC) − ½(d_AC + d_GC)`. Positive for dipyrimidine-biased (UV-like) damage. |
| `gpc_contrast` | `float` | `d_GC − ⅓(d_AC + d_CC + d_TC)`. Keyed on the UPSTREAM (5') base (GpC dinucleotide). NOT the CpG/5mC methylation axis — use `log2_cpg_ratio` for that. |
| `context_heterogeneity_chi2` | `float` | Chi-squared statistic (df = 3) for uniformity across the four contexts. |
| `context_heterogeneity_p` | `float` | Ladder-quantized p-value (`0.9, 0.5, 0.1, 0.05, 0.01, 0.001`). |
| `context_heterogeneity_detected` | `bool` | `chi2 > 7.81` (p < 0.05). |

### Validation and reliability flags

| Field | Type | Description |
|-------|------|-------------|
| `damage_validated` | `bool` | Joint model evidence supports genuine terminal deamination |
| `damage_artifact` | `bool` | Channel A-like enrichment is present, but joint evidence indicates composition/artifact rather than real damage |
| `terminal_inversion` | `bool` | Summary flag: terminal damage rate < interior at either end |
| `inverted_pattern_5prime` / `inverted_pattern_3prime` | `bool` | End-specific terminal depletion flags (full window below interior) |
| `position_0_artifact_5prime` / `position_0_artifact_3prime` | `bool` | Sharp position-0 artifact: C→T at p=0 < p=1 (5′) or isolated G→A spike at p=0 (3′). Indicates adapter/ligation artifact at the terminal nucleotide. Independent of `hexamer_excess_tc`, which tracks full-window chemistry bias. |
| `composition_bias_5prime` / `composition_bias_3prime` | `bool` | Control channel rises with the damage channel, suggesting compositional rather than damage-driven enrichment |
| `is_valid()` | `bool` | `n_reads >= 1000` |
| `is_detection_unreliable()` | `bool` | `true` when inversion or composition-bias flags indicate unreliable reference-free detection |

### Utility functions

```cpp
// Validation state: VALIDATED, CONTRADICTED, or UNVALIDATED
taph::DamageValidationState state = taph::get_damage_validation_state(profile);

// Suppression factor for downstream damage-aware deduplication
float factor = taph::get_damage_suppression_factor(profile);
// 1.0 = validated, 0.5 = unvalidated, 0.0 = artifact
```

### Complement-asymmetry oxidative channels (F, G, H)

These fields populate the `complement_asymmetry` JSON block. Each channel
computes a binomial z-score comparing the trinucleotide-context stop rate in
a terminal window (positions p0–4, adapter-adjusted) against a far-interior
reference (pos ≥ 30, ≥ 50 counts; mid-read fallback pos 5–14). libtaph
stores raw z-scores; threshold application is left to the consumer.

See [Damage types — Channels F, G, H](damage-types.md) for the biochemical
background, trinucleotide context lists, and the adapter-skip gate logic.

#### Channel F — C→A (bottom-strand 8-oxoG complement)

Pre-contexts: TCA/TCG/TAC/TGC. Stop contexts: TAA/TAG/TGA.

| Field | Type | Description |
|-------|------|-------------|
| `channel_f_valid` | `bool` | ≥ 200 total trinucleotide counts across pre + stop |
| `channel_f_z` | `float` | Binomial z-score, 5′ terminal vs interior; negative = depletion |
| `ca_stop_rate_terminal` | `float` | Stop rate in terminal window |
| `ca_stop_rate_interior` | `float` | Stop rate in mid-read interior (fallback baseline) |
| `ca_stop_rate_baseline` | `float` | Stop rate in far-interior reference (pos ≥ 30) |
| `ca_uniformity_ratio` | `float` | Terminal / baseline ratio; > 1 indicates enrichment |
| `channel_f3_valid` | `bool` | ≥ 200 counts for 3′ end |
| `ca_stop_rate_terminal_3prime` | `float` | 3′ terminal stop rate |
| `ca_stop_rate_interior_3prime` | `float` | 3′ interior stop rate |
| `ca_stop_rate_baseline_3prime` | `float` | 3′ far-interior baseline |
| `ca_uniformity_ratio_3prime` | `float` | 3′ terminal / baseline ratio |
| `ca_pre_interior` | `uint64_t` | Raw pre-context count, far-interior (pos ≥ 30) |
| `ca_stop_interior` | `uint64_t` | Raw stop-context count, far-interior |

#### Channel G — C→G (further-oxidized guanine: Gh, Sp)

Pre-contexts: TCA/TAC. Stop contexts: TGA (from TCA), TAG (from TAC).

| Field | Type | Description |
|-------|------|-------------|
| `channel_g_valid` | `bool` | ≥ 200 total counts |
| `channel_g_z` | `float` | Binomial z-score, 5′ terminal vs interior |
| `cg_stop_rate_terminal` | `float` | Terminal stop rate |
| `cg_stop_rate_interior` | `float` | Interior stop rate |
| `cg_stop_rate_baseline` | `float` | Far-interior baseline |
| `cg_uniformity_ratio` | `float` | Terminal / baseline ratio |
| `channel_g3_valid` | `bool` | ≥ 200 counts for 3′ end |
| `cg_stop_rate_terminal_3prime` | `float` | 3′ terminal stop rate |
| `cg_stop_rate_interior_3prime` | `float` | 3′ interior stop rate |
| `cg_stop_rate_baseline_3prime` | `float` | 3′ far-interior baseline |
| `cg_uniformity_ratio_3prime` | `float` | 3′ terminal / baseline ratio |
| `cg_pre_interior` | `uint64_t` | Raw pre-context count, far-interior |
| `cg_stop_interior` | `uint64_t` | Raw stop-context count, far-interior |

#### Channel H — A→T (empirical; mechanism uncertain)

Pre-contexts: AAA/AAG/AGA. Stop contexts: TAA/TAG/TGA.

| Field | Type | Description |
|-------|------|-------------|
| `channel_h_valid` | `bool` | ≥ 200 total counts |
| `channel_h_z` | `float` | Binomial z-score, 5′ terminal (positions p0_h5–4) vs interior |
| `channel_h_z_p2plus` | `float` | Same but terminal window restricted to positions 2–4; more robust because position 0 carries a lower background A→T rate |
| `fgh_adapter_prefixes_excluded` | `uint32_t` | Number of 5′ hexamer prefix bins excluded from F/G/H terminal counts; 0 = no adapter stubs detected |
| `at_stop_rate_terminal` | `float` | Terminal stop rate |
| `at_stop_rate_interior` | `float` | Interior stop rate |
| `at_stop_rate_baseline` | `float` | Far-interior baseline |
| `at_uniformity_ratio` | `float` | Terminal / baseline ratio |
| `channel_h3_valid` | `bool` | ≥ 200 counts for 3′ end |
| `at_stop_rate_terminal_3prime` | `float` | 3′ terminal stop rate |
| `at_stop_rate_interior_3prime` | `float` | 3′ interior stop rate |
| `at_stop_rate_baseline_3prime` | `float` | 3′ far-interior baseline |
| `at_uniformity_ratio_3prime` | `float` | 3′ terminal / baseline ratio |
| `at_pre_interior` | `uint64_t` | Raw pre-context count, far-interior |
| `at_stop_interior` | `uint64_t` | Raw stop-context count, far-interior |

---

## Typical usage patterns

### Streaming with thread-parallel accumulation

```cpp
#include <taph/frame_selector_decl.hpp>
#include <thread>

// Per-thread profiles
std::vector<taph::SampleDamageProfile> partial(n_threads);
for (auto& p : partial) taph::FrameSelector::reset_sample_profile(p);

// Fill partial[i] in parallel...

// Merge on main thread
taph::SampleDamageProfile final_profile;
taph::FrameSelector::reset_sample_profile(final_profile);
for (const auto& p : partial)
    taph::FrameSelector::merge_sample_profiles(final_profile, p);

taph::FrameSelector::finalize_sample_profile(final_profile);
```

### Forcing library type

```cpp
taph::SampleDamageProfile profile = /* ... */;
profile.forced_library_type = taph::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
```

When `forced_library_type != UNKNOWN`, downstream code should use the forced type; `library_type_auto_detected` will be `false`.

---

## Damage Interpretation (`taph/library_interpretation.hpp`)

Higher-level functions that operate on a finalized `SampleDamageProfile` to produce scores, flags, and summaries suitable for reporting.

### Hexamer utilities

```cpp
std::array<char,7> taph::decode_hex(int code);
int taph::encode_hex_at(const std::string& s, int pos);
```

Encode/decode 12-bit 6-mer codes over the ACGT alphabet.

### Adapter stub detection

```cpp
taph::AdapterStubs taph::detect_adapter_stubs(
    const SampleDamageProfile& dp,
    const uint32_t hex3_terminal[4096],
    uint64_t n_hex3);
```

Detects 5' and 3' adapter remnants from hexamer enrichment.  Pass the 4096-element terminal hexamer histogram (`hex3_terminal`) from a pre-scan pass alongside `n_hex3` (total counts).  Returns `AdapterStubs` with `stubs5`/`stubs3` (up to 5 each), `adapter_clipped`, `adapter3_clipped`, and `flag_hex_artifact`.

3' detection is gated on `adapter_clipped` (5' stubs must be present first).

```cpp
void taph::FrameSelector::recompute_fgh_excluding_adapter_prefixes(
    SampleDamageProfile& dp,
    const std::vector<uint32_t>& excl5,
    const std::vector<uint32_t>& excl3);
```

Re-sums all F/G/H terminal counts excluding reads whose first hexamer code is in `excl5`
(5′ end) or last hexamer code is in `excl3` (3′ end), then recomputes `channel_f_z`,
`channel_g_z`, `channel_h_z`, and `channel_h_z_p2plus` from the filtered counts. Must
be called after `finalize_sample_profile`. Sets `fgh_adapter_prefixes_excluded` to
`excl5.size()`. Hexamer codes come from `taph::encode_hex_at` or `taph::encode_hexamer`.
Pass empty vectors to recompute with no exclusion (useful for testing the unfiltered
baseline).

### Hexamer statistics

```cpp
taph::HexStats taph::compute_hex_stats(const SampleDamageProfile& dp);
```

Returns Shannon entropy of terminal/interior hexamer distributions, Jensen-Shannon divergence, and a multinomial G-test (statistic, z-score, p-value) comparing terminal vs interior composition.

### Score functions

| Function | Output | Meaning |
|---|---|---|
| `compute_cpg_score(dp)` | `CpgScore{z, p}` | log₂(CpG / non-CpG d_max) / SE |
| `compute_oxog_interior_score(dp)` | `OxogInteriorScore{z, p}` | Chargaff-symmetric G→T excess over 16 trinucleotide contexts |
| `compute_oxog_trinuc(dp)` | `OxogTrinucResult{cosine, n_ctx}` | Cosine similarity of per-context G→T residuals to empirical 8-oxoG reference |
| `compute_depur_score(dp, is_ss)` | `DepurScore{z, p, z5, z3, shift5, shift3}` | Conjunction test on terminal A/(A+G) [5'] and T/(T+C) [3'] enrichment |

### Damage mask

```cpp
taph::DamageMask taph::compute_damage_mask(
    const SampleDamageProfile& dp, bool is_ss,
    double threshold = 0.05, int min_cov = 100);
```

Marks the first `INTERP_N_POS` (15) positions where terminal damage exceeds `threshold` above background.  Returns `DamageMask` with `pos[15]`, `n_masked`, and `masked_str` (comma-separated indices).

### Library QC flags

```cpp
taph::LibraryQcFlags taph::compute_library_qc_flags(
    const SampleDamageProfile& dp, bool is_ss,
    bool flag_hex_artifact, double jsd,
    double h_term, double short_read_frac);
```

Nine boolean flags covering adapter remnants, hexamer biases, short-read spikes, depurination, absent 3' signal (DS), and inward-displaced G→A.

### Preservation summary

```cpp
taph::PreservationSummary taph::compute_preservation_summary(
    const SampleDamageProfile& dp, bool is_ss,
    bool adapter_clipped, bool flag_hex_artifact,
    double cpg_score_z, double oxog_score_z,
    double oxog_context_cosine, double hex_shift_p);
```

Combines evidence from all channels into a single preservation assessment.  Fields: `authenticity_eff`, `authenticity_evidence`, `d5_raw`, `d5_hexamer_corrected`, `d5_was_corrected`, `oxidation_eff`, `oxidation_evidence`, `qc_risk_eff`, `qc_evidence`, `label` (e.g. `"ancient"`, `"weak"`, `"modern-like"`).

### DamageContextProfile

```cpp
struct taph::DamageContextProfile {
    enum class DominantProcess {
        None, LowDamage, CytosineDeamination, CpgEnrichedDeamination,
        DipyrimidineBiased, OxidativeLike,
        FragmentationBias, LibraryArtifactLikely
    };
    static constexpr const char* method = "training_free";
    bool reference_required = false;
    bool alignment_required = false;

    float terminal_deamination_score, cpg_context_score;
    float dipyrimidine_context_score, oxidative_context_score;
    float fragmentation_context_score, library_artifact_score;

    DominantProcess dominant_process;
    std::string     dominant_process_str;
    std::string     interpretation;

    struct Evidence { /* d_max_{5,3}, lambda_{5,3}, log2_cpg_ratio, cpg_z,
                         dipyr_contrast, ox_gt_asymmetry, s_oxog_{mean,max},
                         purine_enrichment_5prime, hex_shift_z,
                         adapter_clipped, adapter3_clipped, flag_hex_artifact,
                         position_0_artifact_{5,3}prime, n_reads */ } evidence;
};

const char* taph::to_string(DamageContextProfile::DominantProcess p);

taph::DamageContextProfile taph::compute_damage_context_profile(
    const SampleDamageProfile& dp,
    double cpg_z, double hex_shift_z,
    bool adapter_clipped, bool adapter3_clipped, bool flag_hex_artifact);
```

Training-free, reference-free summary aggregator. Six scores are in `[0, 1]` (or `NaN` when not evaluable); `dominant_process` is a deterministic rule over the scores. Underlying raw numbers are mirrored into `evidence` for auditing and re-normalisation by downstream tools. See [methods.md](methods.md#damage-context-profile) for score formulas and the rule.

---

## Length-stratified profile (`taph/length_stratified_profile.hpp`)

Up to four per-length `SampleDamageProfile` instances with shared accumulate / merge / finalize semantics.

```cpp
struct taph::LengthBinStats {
    static constexpr std::size_t MAX_BINS = 4;
    std::array<SampleDamageProfile, MAX_BINS> profiles;
    std::vector<int> edges;
    std::size_t n_bins = 1;
    SampleDamageProfile::LibraryType forced_library_type;

    void configure(const std::vector<int>& new_edges);
    void update(std::string_view seq, int length);
    void merge(const LengthBinStats& other);
    void finalize_all();
    std::size_t bin_index(int length) const;
};
```

`configure(edges)` validates that edges are strictly increasing and at most three splits (four bins) and sets up the per-bin profiles. `update` routes a read by `bin_index(length)` and accumulates into that bin's profile. `merge` requires matching edges. `finalize_all` runs `FrameSelector::finalize_sample_profile` on every populated bin; `forced_library_type` propagates to each bin if set.

## Automatic length-bin edges (`taph/log_length_gmm.hpp`)

```cpp
taph::LogLengthGmmResult taph::detect_log_length_gmm_edges(
    const std::vector<uint64_t>& histogram,
    double log_min, double log_max,
    int min_length, int max_length,
    int max_components = 4);

std::vector<int> taph::detect_quantile_length_edges(
    const std::vector<uint64_t>& histogram,
    double log_min, double log_max,
    int min_length, int max_length,
    int n_bins);
```

`detect_log_length_gmm_edges` fits a 1D Gaussian mixture on a log-length histogram, selects `n_components ∈ [1, max_components]` by BIC, and returns split points between adjacent component means. `LogLengthGmmResult` holds `edges`, `n_components`, `bic`, and `converged`. `detect_quantile_length_edges` is the deterministic fallback when the GMM fails to separate modes; it splits the histogram into equal-count quantile bins.

Typical wiring:

```cpp
auto hist = /* log-length histogram from a pre-scan */;
auto gmm = taph::detect_log_length_gmm_edges(hist, log_min, log_max, lmin, lmax);
std::vector<int> edges = gmm.converged && gmm.n_components > 1
    ? gmm.edges
    : taph::detect_quantile_length_edges(hist, log_min, log_max, lmin, lmax, 3);

taph::LengthBinStats stats;
stats.configure(edges);
for (auto& r : reads) stats.update(r.seq, r.length);
stats.finalize_all();
```

## Length × GC joint mixture (`taph/length_gc_joint_mixture.hpp`)

```cpp
taph::LengthGcJointMixtureResult
taph::fit_length_gc_joint_mixture(const LengthBinStats& stats);
```

Shared-component 2-Gaussian mixture over the flat length × GC cell grid of a finalized `LengthBinStats`. Each cell contributes its `d_max` and `c_sites`. Returns:

| Field | Type | Description |
|-------|------|-------------|
| `d_ancient` | `double` | Damaged-component mean, shared across all cells. |
| `pi_ancient` | `double` | Damaged-component mixing weight. |
| `d_population` | `double` | `c_sites`-weighted mean `d_max` across all cells. |
| `converged` | `bool` | EM convergence state. |
| `separated` | `bool` | `true` if the fit produced a distinct damaged component. |
| `cell_w_ancient[b][g]` | `vector<array<double, N_GC_BINS>>` | Posterior `P(damaged | cell)` per (length bin, GC bin). `-1.0` for cells with insufficient C-sites. Empty when `!separated`. |

Fixed priors match `DamageMixtureModel`: `μ₀ = 0`, `τ₀ = 0.01`, `τ₁ = 0.10`, per-cell sigma floor `0.02`. Use `d_ancient` as a single shared amplitude and `cell_w_ancient[b][g]` to recover per-length-bin ancient fractions inside each GC column.

---

## Hexamer tables and domain classifier (`taph/hexamer_tables.hpp`)

Eight pre-computed 4096-bin hexamer frequency tables plus a softmax domain classifier that scores a read against all of them.

```cpp
enum class taph::Domain {
    META, GTDB, FUNGI, PROTOZOA, INVERTEBRATE,
    PLANT, VERTEBRATE_MAMMALIAN, VERTEBRATE_OTHER, VIRAL
};
```

`Domain::GTDB` covers bacteria + archaea. `Domain::META` is the default for ancient metagenomes and selects the ensemble-of-all-domains path. `parse_domain(name)` accepts `"meta" | "metagenome" | "all"` (all three map to `Domain::META`), `"gtdb" | "bacteria" | "archaea"`, `"fungi"`, `"protozoa"`, `"invertebrate"`, `"plant"`, `"mammal" | "vertebrate_mammalian"`, `"vertebrate" | "vertebrate_other"`, `"viral"`, and falls back to `Domain::META` on anything unrecognised. `domain_name(d)` gives the canonical string.

### Encoding and table lookup

```cpp
uint32_t taph::encode_hexamer(const char* seq);                         // 12-bit code, UINT32_MAX on N
float    taph::get_hexamer_freq(uint32_t code);                         // active domain
float    taph::get_hexamer_freq(uint32_t code, Domain d);
float    taph::get_hexamer_freq(const char* seq);                       // active domain
float    taph::get_hexamer_freq(const char* seq, Domain d);
float    taph::get_hexamer_score(const char* seq);                      // log2(freq / uniform)
float    taph::get_hexamer_score(const char* seq, Domain d);
```

`get_hexamer_score` takes a sequence pointer (there is no `uint32_t` overload), returns `log2(freq / (1/4096))`, and floors at `-10.0` for unseen hexamers.

### Active domain and ensemble mode

```cpp
void          taph::set_active_domain(Domain d);
taph::Domain& taph::active_domain();             // reference to the global
taph::Domain  taph::get_active_domain();         // value accessor
void          taph::set_ensemble_mode(bool enabled);
bool          taph::get_ensemble_mode();
void          taph::set_domain_probs(const MultiDomainResult& r);
taph::MultiDomainResult& taph::get_domain_probs();
float         taph::get_ensemble_hexamer_freq(uint32_t code);
```

`set_active_domain` writes a global `Domain` (not thread-local) that gates the default-table overloads; read it via `get_active_domain()`. The ensemble flag (`ensemble_mode_enabled()`) and the posterior buffer (`current_domain_probs()`) are thread-local, so per-thread pipelines can carry their own blend. With ensemble mode off, `get_ensemble_hexamer_freq` falls through to the active-domain table. With ensemble mode on, it returns the posterior-weighted average of the eight tables using the probabilities written by `set_domain_probs`.

### Domain classifier

```cpp
struct taph::MultiDomainResult {
    float gtdb_prob, fungi_prob, protozoa_prob, invertebrate_prob;
    float plant_prob, vertebrate_mammalian_prob, vertebrate_other_prob, viral_prob;
    Domain best_domain;
    float  best_score;
};

taph::MultiDomainResult taph::score_all_domains(const std::string& seq, int frame);
```

Sums `log2(freq[h] * 4096 + 1e-10)` across every in-frame hexamer for each of the eight domains, then normalises with temperature-scaled softmax (`exp(0.1 * (score - max))`). `best_domain` is the argmax and `best_score` is the pre-softmax log-sum for that domain.

Typical ensemble flow for mixed-community samples:

```cpp
auto probs = taph::score_all_domains(read, /*frame=*/0);
taph::set_domain_probs(probs);
taph::set_ensemble_mode(true);
// Subsequent get_ensemble_hexamer_freq(code) returns the posterior-weighted
// background rather than committing to one domain table.
```

### Hexamer fields on `SampleDamageProfile`

Raw 4096-bin terminal and interior histograms populated during accumulation and consumed by the interpretation functions in `taph/library_interpretation.hpp`.

| Field | Type | Description |
|-------|------|-------------|
| `hexamer_count_5prime[4096]` | `double[4096]` | Hexamer counts at 5' positions 0-5. |
| `hexamer_count_interior[4096]` | `double[4096]` | Hexamer counts at the interior (middle third). |
| `n_hexamers_5prime` | `size_t` | Total counted terminal hexamers. |
| `n_hexamers_interior` | `size_t` | Total counted interior hexamers. |
| `hexamer_terminal_tc` | `float` | T/(T+C) ratio at 5′ terminal positions across all C/T hexamer pairs. |
| `hexamer_interior_tc` | `float` | T/(T+C) ratio at interior positions (same hexamer pairs). |
| `hexamer_excess_tc` | `float` | `hexamer_terminal_tc − hexamer_interior_tc`. **Library chemistry / priming bias metric — not an adapter remnant signal.** A negative value means the library prep preferentially produces C-starting read starts (e.g. SS ligation junction composition), which enriches C-context trinucleotides in the channels F and G pre-count denominator and suppresses their z-scores toward negative values. Adapter remnants are separately tracked by `adapter_clipped`, `flag_hex_artifact`, and `position_0_artifact_*`. |
| `hexamer_damage_llr` | `float` | Log-likelihood ratio of T-enrichment over C-enrichment at 5′ terminal hexamers vs interior. Positive = T-starting hexamers enriched (consistent with deamination). Used to detect `inverted_pattern_5prime` when < −0.02 and no p0 artifact. |

Paired with the 3' terminal hexamer histogram produced by the pre-scan (passed to `detect_adapter_stubs` as `hex3_terminal[4096]` with `n_hex3`). The 3' histogram is kept out of `SampleDamageProfile` so it does not bloat per-thread state in accumulators that do not need it.
