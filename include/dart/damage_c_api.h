/**
 * damage_c_api.h — C-compatible API for libdart-damage
 *
 * Exposes two-pass aDNA damage estimation and read correction through
 * opaque handles with extern "C" linkage.  Any C or C++ project can
 * include this header regardless of the C++ standard it uses.
 *
 * Usage:
 *   // Pass 1 – accumulate reads, estimate sample damage
 *   dart_profile_t *p = dart_profile_create();
 *   while (read_fastq(&seq, &len))
 *       dart_profile_add_read(p, seq, len);
 *   dart_profile_finalize(p);
 *
 *   // Inspect results
 *   if (dart_profile_dmax(p) >= threshold && dart_profile_is_reliable(p)) {
 *       // Pass 2 – correct each read before k-mer extraction
 *       dart_correct_read(p, seq, len, buf, confidence);
 *   }
 *   dart_profile_destroy(p);
 */

#ifndef DART_DAMAGE_C_API_H
#define DART_DAMAGE_C_API_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle to a SampleDamageProfile + finalized model state. */
typedef struct dart_profile_t dart_profile_t;

/* ── Pass 1: estimation ──────────────────────────────────────────────────── */

/** Allocate a new profile accumulator.  Returns NULL on OOM. */
dart_profile_t *dart_profile_create(void);

/** Release all memory owned by the handle. */
void dart_profile_destroy(dart_profile_t *p);

/**
 * Feed one raw DNA read (no quality scores needed).
 * Not thread-safe: use one handle per thread and merge with
 * dart::FrameSelector::merge_sample_profiles() before finalizing.
 * Calls after dart_profile_finalize() are silently ignored.
 */
void dart_profile_add_read(dart_profile_t *p, const char *seq, size_t len);

/**
 * Finalize the profile after all reads have been added.
 * Must be called exactly once before any getter or dart_correct_read.
 */
void dart_profile_finalize(dart_profile_t *p);

/* ── Result accessors (call after finalize) ─────────────────────────────── */

/** Combined D_max (0-1): maximum terminal damage rate across both ends. */
float dart_profile_dmax(const dart_profile_t *p);

/** Exponential decay constant for C→T damage from the 5' end. */
float dart_profile_lambda5(const dart_profile_t *p);

/** Exponential decay constant for G→A damage from the 3' end. */
float dart_profile_lambda3(const dart_profile_t *p);

/**
 * Inferred library type:
 *   0 = UNKNOWN
 *   1 = DOUBLE_STRANDED
 *   2 = SINGLE_STRANDED
 */
int dart_profile_library_type(const dart_profile_t *p);

/** 1 if both C→T and G→A channels agree on damage, 0 otherwise. */
int dart_profile_damage_validated(const dart_profile_t *p);

/** 1 if terminal signal looks like adapter/composition artifact, not real damage. */
int dart_profile_damage_artifact(const dart_profile_t *p);

/**
 * 1 if the estimate is statistically reliable (enough reads, consistent signal).
 * When 0, treat d_max as uninformative.
 */
int dart_profile_is_reliable(const dart_profile_t *p);

/* ── Pass 2: per-read correction ────────────────────────────────────────── */

/**
 * Correct damage in one DNA read using the finalized profile.
 *
 * Reverts T→C at 5' positions and A→G at 3' positions where the
 * position-dependent damage probability exceeds confidence_threshold.
 *
 * @param p                   Finalized profile handle.
 * @param seq                 Input sequence (ASCII, may be upper or lower case).
 * @param len                 Sequence length in bytes.
 * @param out_buf             Caller-allocated buffer of at least len+1 bytes.
 *                            Receives the corrected sequence (NUL-terminated).
 * @param confidence_threshold  Minimum damage probability to trigger correction
 *                              (0–1; typical values 0.2–0.5).
 * @return Number of bases corrected.
 */
size_t dart_correct_read(const dart_profile_t *p,
                         const char           *seq,
                         size_t                len,
                         char                 *out_buf,
                         float                 confidence_threshold);

#ifdef __cplusplus
}
#endif

#endif /* DART_DAMAGE_C_API_H */
