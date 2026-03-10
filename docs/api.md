# API Reference

All public symbols live in the `dart` namespace. Include `<dart/frame_selector_decl.hpp>`.

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
| `sequences` | `const std::vector<std::string>&` | DNA sequences (ACGT, any length ≥ 20) |

**Returns** A fully populated `SampleDamageProfile`. Call `profile.is_valid()` to check whether enough reads were processed (≥ 1000).

---

### `update_sample_profile`

```cpp
static void update_sample_profile(
    SampleDamageProfile& profile,
    const std::string& seq);
```

Accumulate one read into `profile`. Thread-safe per profile object (not across profiles).

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

Zero all counters and reset Channel A, Channel B (5'), joint model, and mixture model fields. Use before reusing a profile object for a new sample.

> **Note:** Channel B₃′, Channel C, Channel D, GC histogram, and alignability accumulator arrays are **not** explicitly reset. Reusing a profile via `reset_sample_profile` followed by `update_sample_profile` calls will leave those fields with stale values. For a guaranteed clean state, construct a new `SampleDamageProfile()` (zero-initialized) instead.

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
| `damage_validated` | `bool` | Channel A and B agree |
| `damage_artifact` | `bool` | Channel A fires, Channel B contradicts |
| `terminal_inversion` | `bool` | Terminal damage rate < interior (detection unreliable) |
| `composition_bias_5prime/3prime` | `bool` | Control channel comparable to damage channel, likely bias |
| `is_valid()` | `bool` | `n_reads >= 1000` |
| `is_detection_unreliable()` | `bool` | Any inversion or composition bias flag is set |

### Utility functions

```cpp
// Validation state: VALIDATED, CONTRADICTED, or UNVALIDATED
dart::DamageValidationState state = dart::get_damage_validation_state(profile);

// Suppression factor for downstream damage-aware deduplication
float factor = dart::get_damage_suppression_factor(profile);
// 1.0 = validated, 0.5 = unvalidated, 0.0 = artifact
```

---

## Typical usage patterns

### Streaming with thread-parallel accumulation

```cpp
#include <dart/frame_selector_decl.hpp>
#include <thread>

// Per-thread profiles
std::vector<dart::SampleDamageProfile> partial(n_threads);
for (auto& p : partial) dart::FrameSelector::reset_sample_profile(p);

// Fill partial[i] in parallel...

// Merge on main thread
dart::SampleDamageProfile final_profile;
dart::FrameSelector::reset_sample_profile(final_profile);
for (const auto& p : partial)
    dart::FrameSelector::merge_sample_profiles(final_profile, p);

dart::FrameSelector::finalize_sample_profile(final_profile);
```

### Forcing library type

```cpp
dart::SampleDamageProfile profile = /* ... */;
profile.forced_library_type = dart::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
```

When `forced_library_type != UNKNOWN`, downstream code should use the forced type; `library_type_auto_detected` will be `false`.
