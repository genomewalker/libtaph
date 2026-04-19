/**
 * damage_c_api.h — C-compatible API for libtaph
 *
 * Exposes two-pass aDNA damage estimation, read correction, and read
 * masking through opaque handles with extern "C" linkage.  Any C or
 * C++ project can include this header regardless of the C++ standard
 * it uses.
 *
 * Usage (correction):
 *   taph_profile_t *p = taph_profile_create();
 *   while (read_fastq(&seq, &len))
 *       taph_profile_add_read(p, seq, len);
 *   taph_profile_finalize(p);
 *
 *   if (taph_profile_dmax(p) >= threshold && taph_profile_is_reliable(p)) {
 *       // Back-convert damaged T→C and A→G before k-mer extraction:
 *       taph_correct_read(p, seq, len, buf, confidence);
 *       // Or replace damaged positions with 'N' instead:
 *       taph_mask_read(p, seq, len, buf, confidence, 'N');
 *   }
 *   taph_profile_destroy(p);
 */

#ifndef TAPH_DAMAGE_C_API_H
#define TAPH_DAMAGE_C_API_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle to a SampleDamageProfile + finalized model state. */
typedef struct taph_profile_t taph_profile_t;

/* ── Pass 1: estimation ──────────────────────────────────────────────────── */

/** Allocate a new profile accumulator.  Returns NULL on OOM. */
taph_profile_t *taph_profile_create(void);

/** Release all memory owned by the handle. */
void taph_profile_destroy(taph_profile_t *p);

/**
 * Feed one raw DNA read (no quality scores needed).
 * Not thread-safe: use one handle per thread and merge with
 * taph::FrameSelector::merge_sample_profiles() before finalizing.
 * Calls after taph_profile_finalize() are silently ignored.
 */
void taph_profile_add_read(taph_profile_t *p, const char *seq, size_t len);

/**
 * Finalize the profile after all reads have been added.
 * Must be called exactly once before any getter or taph_correct_read.
 */
void taph_profile_finalize(taph_profile_t *p);

/* ── Result accessors (call after finalize) ─────────────────────────────── */

/** Combined D_max (0-1): maximum terminal damage rate across both ends. */
float taph_profile_dmax(const taph_profile_t *p);

/** Exponential decay constant for C→T damage from the 5' end. */
float taph_profile_lambda5(const taph_profile_t *p);

/** Exponential decay constant for G→A damage from the 3' end. */
float taph_profile_lambda3(const taph_profile_t *p);

/**
 * Inferred library type:
 *   0 = UNKNOWN
 *   1 = DOUBLE_STRANDED
 *   2 = SINGLE_STRANDED
 */
int taph_profile_library_type(const taph_profile_t *p);

/** 1 if both C→T and G→A channels agree on damage, 0 otherwise. */
int taph_profile_damage_validated(const taph_profile_t *p);

/** 1 if terminal signal looks like adapter/composition artifact, not real damage. */
int taph_profile_damage_artifact(const taph_profile_t *p);

/**
 * 1 if the estimate is statistically reliable (enough reads, consistent signal).
 * When 0, treat d_max as uninformative.
 */
int taph_profile_is_reliable(const taph_profile_t *p);

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
size_t taph_correct_read(const taph_profile_t *p,
                         const char           *seq,
                         size_t                len,
                         char                 *out_buf,
                         float                 confidence_threshold);

/**
 * Mask damage in one DNA read using the finalized profile.
 *
 * Replaces T at 5' positions and A at 3' positions where the
 * position-dependent damage probability exceeds confidence_threshold
 * with mask_char (typically 'N').  k-mer extractors that skip k-mers
 * containing the mask character will naturally avoid damage-affected
 * positions without assuming the original base.
 *
 * @param p                   Finalized profile handle.
 * @param seq                 Input sequence (ASCII, upper or lower case).
 * @param len                 Sequence length in bytes.
 * @param out_buf             Caller-allocated buffer of at least len+1 bytes.
 *                            Receives the masked sequence (NUL-terminated).
 * @param confidence_threshold  Minimum damage probability to trigger masking
 *                              (0–1; typical values 0.2–0.5).
 * @param mask_char           Character to write at masked positions (e.g. 'N').
 * @return Number of bases masked.
 */
size_t taph_mask_read(const taph_profile_t *p,
                      const char           *seq,
                      size_t                len,
                      char                 *out_buf,
                      float                 confidence_threshold,
                      char                  mask_char);

#ifdef __cplusplus
}
#endif

#endif /* TAPH_DAMAGE_C_API_H */
