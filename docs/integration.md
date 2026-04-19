# Classifier Integration Guide

libtaph ships a C-compatible API (`dart/damage_c_api.h`) that lets any
classifier link the library regardless of its own C++ standard.  The
implementation compiles at C++17; the header is plain C and can be included
from C, C++14, or any later standard.

---

## CMake setup

```cmake
# In your top-level or src/ CMakeLists.txt
add_subdirectory(path/to/libtaph)
target_link_libraries(your_classifier PRIVATE taph)
```

`taph` sets its own `CXX_STANDARD 17` — it does not affect your
target's standard.

---

## Two-pass workflow

### Pass 1 — estimate sample damage

Feed all reads through the accumulator before classification begins.

```c
#include <taph/damage_c_api.h"

taph_profile_t *profile = taph_profile_create();
if (!profile) { /* OOM */ }

/* For each read in every input file: */
taph_profile_add_read(profile, seq, len);

taph_profile_finalize(profile);    /* fit the model — call exactly once */
```

Inspect the result:

```c
float dmax      = taph_profile_dmax(profile);      /* 0–1 */
int   lib_type  = taph_profile_library_type(profile); /* 0=UNKNOWN 1=DS 2=SS */
int   validated = taph_profile_damage_validated(profile);
int   artifact  = taph_profile_damage_artifact(profile);
int   reliable  = taph_profile_is_reliable(profile);
```

Only activate correction when the signal is genuine:

```c
if (dmax >= 0.05f && !artifact && reliable) {
    /* enable Pass 2 */
}
```

### Pass 2 — correct or mask each read

Two functions are available. Both use the same per-position damage probabilities
(empirical `damage_rate_5prime/3prime` arrays for the first 15 positions,
exponential extrapolation `d_max × e^{-λ × pos}` beyond that).

**`taph_correct_read`** — back-convert damaged bases to their inferred originals:

```c
char corrected[MAX_READ_LEN + 1];
size_t n_fixes = taph_correct_read(profile,
                                   seq, len,
                                   corrected,
                                   0.30f);   /* confidence threshold */
/* use corrected instead of seq */
```

Reverts **T → C** at 5′ positions and **A → G** at 3′ positions where the
damage probability ≥ threshold.  Appropriate when you want to recover the
likely original sequence before alignment or classification.

**`taph_mask_read`** — replace damaged positions with a mask character:

```c
char masked[MAX_READ_LEN + 1];
size_t n_masked = taph_mask_read(profile,
                                 seq, len,
                                 masked,
                                 0.30f,   /* confidence threshold */
                                 'N');    /* mask character */
/* use masked instead of seq */
```

Writes `mask_char` (typically `'N'`) at any position where the C→T or G→A
damage probability ≥ threshold, leaving all other bases unchanged.  Useful
for k-mer classifiers (e.g. Metabuli) that skip k-mers containing `'N'`:
damaged positions are excluded from k-mer extraction without assuming the
original base.  Prefer masking over correction when the downstream tool
handles ambiguous bases natively.

### Cleanup

```c
taph_profile_destroy(profile);
```

---

## Threading

`taph_profile_add_read` is **not thread-safe**.  For multi-threaded Pass 1,
create one profile per thread and merge before finalizing:

```cpp
#include <taph/frame_selector_decl.hpp"

std::vector<taph_profile_t *> per_thread(n_threads);
for (auto &p : per_thread) p = taph_profile_create();

/* ... fill per_thread[tid] from each thread's reads ... */

/* Merge into thread 0's profile (single-threaded): */
for (int i = 1; i < n_threads; ++i)
    taph::FrameSelector::merge_sample_profiles(
        per_thread[0]->profile, per_thread[i]->profile);

taph_profile_finalize(per_thread[0]);
/* use per_thread[0] for Pass 2 */
```

`taph_correct_read` is **thread-safe** — it only reads the finalized profile.

---

## API reference

### `taph_profile_create`
```c
taph_profile_t *taph_profile_create(void);
```
Allocate a new profile accumulator.  Returns `NULL` on OOM.

---

### `taph_profile_destroy`
```c
void taph_profile_destroy(taph_profile_t *p);
```
Free all resources.  Safe to call with `NULL`.

---

### `taph_profile_add_read`
```c
void taph_profile_add_read(taph_profile_t *p, const char *seq, size_t len);
```
Accumulate one DNA read into the profile.  Silently ignored after
`taph_profile_finalize`.  Not thread-safe (see Threading above).

---

### `taph_profile_finalize`
```c
void taph_profile_finalize(taph_profile_t *p);
```
Fit the exponential decay model.  Must be called exactly once before any
getter or `taph_correct_read` / `taph_mask_read`.  Subsequent calls are no-ops.

---

### Getters (require finalized profile)

| Function | Returns |
|---|---|
| `taph_profile_dmax(p)` | Combined D_max (0–1) |
| `taph_profile_lambda5(p)` | 5′ decay constant λ |
| `taph_profile_lambda3(p)` | 3′ decay constant λ |
| `taph_profile_library_type(p)` | 0=UNKNOWN, 1=DS, 2=SS |
| `taph_profile_damage_validated(p)` | 1 if both channels agree |
| `taph_profile_damage_artifact(p)` | 1 if signal looks like adapter bias |
| `taph_profile_is_reliable(p)` | 1 if estimate is statistically reliable |

All return safe defaults (0 / 0.0) before finalization.

---

### `taph_correct_read`
```c
size_t taph_correct_read(const taph_profile_t *p,
                          const char *seq, size_t len,
                          char *out_buf,
                          float confidence_threshold);
```
Back-convert damaged bases (T→C at 5′, A→G at 3′) where the position-dependent
damage probability ≥ `confidence_threshold`.  `out_buf` must be at least
`len + 1` bytes.  Returns the number of bases corrected.  Returns 0 if the
profile is not finalized.

---

### `taph_mask_read`
```c
size_t taph_mask_read(const taph_profile_t *p,
                      const char *seq, size_t len,
                      char *out_buf,
                      float confidence_threshold,
                      char mask_char);
```
Replace damaged positions with `mask_char` (e.g. `'N'`) where the
position-dependent damage probability ≥ `confidence_threshold`.  `out_buf`
must be at least `len + 1` bytes.  Returns the number of bases masked.
Returns 0 if the profile is not finalized.
