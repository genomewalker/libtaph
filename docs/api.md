# API Reference

All public symbols live in the `taph` namespace. Include `<taph/frame_selector_decl.hpp>`.

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
| `mixture_d_reference` | `float` | Damage rate in GC > 50% bins (metaDMG proxy) |
| `mixture_K` | `int` | Number of mixture components selected by BIC |
| `mixture_converged` | `bool` | Whether the EM algorithm converged |
| `gc_stratified_valid` | `bool` | At least one GC bin has a valid estimate |

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
| `cpg_contrast` | `float` | `d_GC − ⅓(d_AC + d_CC + d_TC)`. Positive for CpG-dominant (methylation-driven) damage. |
| `context_heterogeneity_chi2` | `float` | Chi-squared statistic (df = 3) for uniformity across the four contexts. |
| `context_heterogeneity_p` | `float` | Ladder-quantized p-value (`0.9, 0.5, 0.1, 0.05, 0.01, 0.001`). |
| `context_heterogeneity_detected` | `bool` | `chi2 > 7.81` (p < 0.05). |

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

Paired with the 3' terminal hexamer histogram produced by the pre-scan (passed to `detect_adapter_stubs` as `hex3_terminal[4096]` with `n_hex3`). The 3' histogram is kept out of `SampleDamageProfile` so it does not bloat per-thread state in accumulators that do not need it.
