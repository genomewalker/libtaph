# Deamination d_max — the bulk↔enriched RANGE contract (reference-free)

## Problem
A single d_max scalar cannot be both an all-reads population estimate and a
damaged-subpopulation (mapping/metaDMG-A_b-comparable) estimate. The pooled
terminal fit COLLAPSES under length-mixing (ancient-short + modern-long → one
diluted decay → ~0), which is the ss near-zero-d5 bug (Task #27). The fix is not
a better scalar; it is to emit an INTERVAL whose two endpoints are two different,
each-honest estimands, plus a confidence that is bound to authentication.

## Contract (uniform across channels)
`{bulk, enriched, range:{lo,hi}, confidence}` per channel. One estimand per
endpoint (bulk always bulk, enriched always enriched — never sorted). Only
deamination and oxidation populate `enriched`; all other channels emit
`enriched:null, mechanism:none`. `confidence` widens/flags when the enriched
endpoint is depth-biased or authentication fails.

### bulk endpoint — length-recombined pooled (collapse cure)
Read-weighted mean of the PER-BIN terminal d5 (each length bin fit separately, so
no ancient-short/modern-long dilution), recombined by n_reads:
`bulk = Σ n_i·d5_i / Σ n_i`. Cures the spurious 0 honestly (no constant):
3S-205A-ss .035→.077, 10S-11-ss .166→.232, 10S-11-ds .210→.399; stays ~0 for
blanks (11S-NC-ds .002) and low-damage (21S-IS100-ss .001). Never overshoots.

### enriched endpoint — T2, the prior-free LSD classifier (0.39-free)
`enriched_upper = d_max_5_damaged` of the SHORTEST/purest length bin (exploits
d_anc length-invariance → least attenuated, closest to truth). Companion
`enriched_wt` = n_damaged-weighted mean across bins (conservative/robust). CI on
the top-bin point from that bin's damaged read count (binomial); when the CI is
wide (low effective N) the interval widens and confidence drops to LOW — the
noise shows up as WIDTH, not a false-precise truth-reaching point. No N-threshold
constant.

Provenance (call-graph of the exact shipped field, verified 2026-07):
`profile_json.cpp:937 d_max_5_damaged` ← fqdup
`damage_profile_lengthstrat.cpp:569 pick_dmax(per_pos_5prime_ct_damaged, llr_bg_5)`
(empirical terminal C→T rate of classified-ancient reads, counts − bg) ← counts
from libtaph `lsd_accumulate` over reads where `lsd_classify_read = lsd_llr_score>0`
(`lsd_accumulator.cpp:79 → :36`). The classifier uses `p.d_anc` (data-driven =
mixture_d_damaged if converged else max(d5,d3), `lsd_accumulator.cpp:14-16`),
NEVER `p.d_anc_contract`. The chain does not enter `read_ancient_llr.cpp`.

### Why T2 over the T1 length×GC joint mixture
`by_length_joint.d_damaged` (T1, `fit_length_gc_joint_mixture`) is also 0.39-free
and monotone across the truth gradient, but it is MORE attenuated than T2 (it is a
population mixture mean under the ss identifiability ceiling below). T2 reaches
82–103% of truth where T1 reaches ~55–90%. T1 is retained as the fallback/diagnostic.

## REJECTED path — the per-read POSTERIOR (do not re-litigate)
`read_ancient_posterior` / `damage_evidence_llr` (`read_ancient_llr.cpp`) were
rejected as the enriched estimator for TWO dispositive reasons:
1. It REFUSES reference-free: returns `nullopt` unless `pi.state==DETECTED`
   (Briggs λ unfitted in ref-free mode, `profile_json.cpp:3007`).
2. It is 0.39-ANCHORED: amplitude `A = D_MAX_CONSERVED` (`read_ancient_llr.cpp:480`,
   denom at :212/:323) — a fixed imported amplitude, so its "reach" would be partly
   the constant, re-opening circularity.
This is a genuinely DIFFERENT function from the LSD classifier that produces the
shipped T2. Conflating the two ("per-read posterior" vs "per-read LSD classifier")
caused an internal contradiction during design; the classifier is 0.39-free, the
posterior is not, and only the classifier ships.

## The 0.39 constant (D_MAX_CONSERVED, damage_estimate.hpp:14) — status
NOT removed. Two residues, neither on a deam-range endpoint:
- DEAD write: `lsd_accumulator.cpp:31` assigns `d_anc_contract = D_MAX_CONSERVED`,
  never read (exactly 2 occurrences in libtaph: this assign + the struct decl).
  Deletable in the constant-purge (Task #26).
- LIVE reads: `read_ancient_llr.cpp:212/323/480` inside the rejected per-read
  posterior, `nullopt` in ref-free → feeds no deam endpoint.
All THREE deam-range endpoints (bulk wmean, T1 by_length_joint, T2 d_max_5_damaged)
are 0.39-free.

## The reference-free ss identifiability CEILING (scientific result)
ss `pi_damaged` saturates 0.90–0.99 with `joint_separated=false`: reference-free,
there is NO clean non-damaged reference class within an ss library, so the mixture
`d_damaged` degrades toward the attenuated POPULATION mean rather than the damaged
amplitude. Hence the ss enriched endpoint sits BELOW truth by construction:

| lib (baseline, no-fix) | truth | bulk_raw | bulk_wmean | T2_top | T2_top/truth |
|---|---|---|---|---|---|
| 3S-205A-ss | .414 | .035 | .077 | .340 | 82% |
| 10S-11-ss  | .455 | .166 | .232 | .354 | 78% |
| 1S-198A-ss | .326 | .013 | .011 | .303 | 93% |
| 3S-205A-ds | .418 | .239 | .270 | .402 | 96% |
| 10S-11-ds  | .490 | .210 | .399 | .399 | 81% |
| 1S-198A-ds | .371 | .130 | .123 | .359 | 97% |

DS ESCAPES the ceiling via the independent 3' strand (an internal non-damaged
handle), so ds enriched reaches ~81–97% while ss reaches ~78–93% and never brackets
truth for the highest-damage ss libs. This is the honest answer to "why is ss
enriched below truth" — an identifiability limit of reference-free ss, not a bug.

## The ss BLANK — reference-free CANNOT separate it (honest limit)
π-saturation does NOT separate blank from real: 11S-NC-ss (blank) π≈1.0 AND
10S-11-ss (real) π .987 both saturate. On the controlled no-fix baseline:
- ds blank (11S-NC-ds): confidence NONE — `source→NONE` via the w_length gate
  (Status as of 2026-10-20: that gate never passed in any library, since w_length reads
  ~0.5 everywhere; finalize_dmax now sets NONE in that regime unconditionally.)
  (commit 6afc8e7, `finalize_dmax.cpp:372`); correctly suppressed.
- ss blank (11S-NC-ss): confidence HIGH, T2_top .350 — auth_eff .898,
  dominant_process cytosine_deamination, decay_llr_5prime +75149>0, source≠NONE.
  EVERY authentication gate reads "real damage."
Conclusion: reference-free CANNOT distinguish the ss blank from real ss damage on
amplitude, π, decay-sign, OR the auth ensemble. `range_flag`/`confidence=LOW` is a
DISCLOSURE of this limit, not a fix; we do NOT hard-suppress with a threshold
(would break calibration-free). The auth ensemble (auth_eff→0 AND
dominant_process∈{library_artifact_likely,low_damage}) catches ds blanks and
low-damage libs but NOT the ss blank, which presents as fully authentic.

## confidence wiring (constant-free)
- NONE: `source==NONE` OR `decay_llr_5prime ≤ 0` — not evidence of ancient damage;
  the enriched field must be explicitly flagged, never left as a bare number a
  consumer could read as damage.
- LOW: blank ensemble (auth_eff<0.10 AND dom∈{library_artifact_likely,low_damage}),
  OR top-bin CI wide relative to the point (diluted-lib guard).
- HIGH: otherwise.
Drivers auth_eff and dominant_process are themselves emitted without any 0.39/fixed
damage cutoff (producers `library_interpretation_summary.cpp`, `profile_json.cpp`).

## Known overshoot (track)
18S-C-ds T2_top .404 vs truth .270 (also flagged in the original REVIEW_BRIEF,
.42 vs .27). Mid-gradient ds; CI tight (n_damaged 78M) so the width guard does not
catch it — an amplitude overshoot, not a precision issue. Whether the metaDMG A_b
under-reads 18S-C or the terminal rate genuinely exceeds it is unresolved.
