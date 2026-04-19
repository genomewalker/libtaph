# API Reference

All public symbols live in the `taph` namespace. Include `<dart/frame_selector_decl.hpp>`.

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
| `mixture_d_reference` | `float` | Damage rate in GC > 50% bins (metaDMG proxy) |
| `mixture_K` | `int` | Number of mixture components selected by BIC |
| `mixture_converged` | `bool` | Whether the EM algorithm converged |
| `gc_stratified_valid` | `bool` | At least one GC bin has a valid estimate |

### Validation and reliability flags

| Field | Type | Description |
|-------|------|-------------|
| `damage_validated` | `bool` | Joint model evidence supports genuine terminal deamination |
| `damage_artifact` | `bool` | Channel A-like enrichment is present, but joint evidence indicates composition/artifact rather than real damage |
| `terminal_inversion` | `bool` | Summary flag: terminal damage rate < interior at either end |
| `inverted_pattern_5prime` / `inverted_pattern_3prime` | `bool` | End-specific terminal depletion flags |
| `position_0_artifact_5prime` / `position_0_artifact_3prime` | `bool` | Position-0 depletion with downstream enrichment; likely adapter/ligation artifact |
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
    double oxog_trinuc_cosine, double hex_shift_p);
```

Combines evidence from all channels into a single preservation assessment.  Fields: `authenticity_eff`, `authenticity_evidence`, `d5_raw`, `d5_hexamer_corrected`, `d5_was_corrected`, `oxidation_eff`, `oxidation_evidence`, `qc_risk_eff`, `qc_evidence`, `label` (e.g. `"ancient"`, `"weak"`, `"modern-like"`).
