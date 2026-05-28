# Interpreting libtaph output

This page is a practical guide for researchers reading a `SampleDamageProfile`
or its JSON serialisation. It does not document every field — the
[API reference](api.md) does that — but explains what the headline numbers mean,
what counts as a strong or weak signal, and how to combine fields into a
defensible call about a sample. For the underlying biology and statistics see
[Damage types](damage-types.md) and [Methods](methods.md).

A useful default rule of thumb: **never look at one number in isolation.**
libtaph emits a dozen partly redundant statistics so that compositional
artefacts and genuine damage can be distinguished. A high `d_max_5prime` on its
own is suggestive; the same number with a positive Channel B validation, an
elevated F/G/H trio, and a clean QC flag panel is a confident call.

## 1. Core deamination metrics

### `d_max_5prime` and `d_max_3prime`

These are the calibrated terminal C→T excess at the 5′ end and the calibrated
terminal G→A excess at the 3′ end (or terminal C→T at the 3′ end for
single-stranded original-orientation libraries). They are computed as
`A / (1 − b)` where `A` is the fitted exponential amplitude above the interior
background and `b` is the anchored background rate. Physically they answer:
"by how much does the cytosine-deamination–driven substitution rate at the very
end of the read exceed what the interior of the same library shows?"

A working interpretation scale, for chemistry-untreated libraries:

| `d_max` value | Interpretation |
|---|---|
| ≥ 0.20 | Strong signal — well-preserved aDNA, often Pleistocene |
| 0.08 – 0.20 | Moderate — typical Holocene aDNA |
| 0.02 – 0.08 | Weak — borderline; needs cross-channel support |
| < 0.02 | Effectively absent — modern, USER-treated, or contaminated |

These thresholds are heuristics, not hard cutoffs. When the library has been
chemically treated (e.g. USER, single-stranded), a moderate `d_max` may still
be meaningful — see the source string returned by `d_max_source_str()`.

### `lambda_5prime` and `lambda_3prime`

`lambda` is the decay constant of the fitted exponential. Larger λ means the
elevated rate falls back to background more steeply, i.e. damage is concentrated
in the first one or two terminal positions. Lower λ means a slower decay with
damage extending several bases inward, reflecting longer single-stranded overhangs
(Briggs 2007; Jónsson 2013), which tend to accompany older material. As a
guideline, λ ≈ 0.2–0.4 is typical of well-preserved aDNA; λ < 0.1 indicates long
damage tails from extended single-stranded overhangs; λ > 0.5 indicates damage
concentrated in the first 1–2 positions, consistent with short overhangs, chemical
treatment, or near-modern material (but not exclusively). These bands are heuristic
— they overlap and do not replace cross-channel validation.

### `d_max_combined` versus the per-end values

`d_max_combined` is the asymmetry-aware single-number summary; the string
returned by `d_max_source_str()` records how it was derived. Use the
combined value when reporting a single damage number per library. Use the
individual `d_max_5prime` / `d_max_3prime` when comparing the two ends, when
diagnosing library type, or when one end is known to be compromised by
adapter remnants or end-repair artefacts.

If `d_max_source_str()` returns `"5prime_only"`, `"3prime_only"`, or
`"none"`, libtaph has decided that one or both ends are inverted (terminal
rate below interior) and should not be averaged. Treat the combined value
with care in those cases.

### `bg_5prime` and `bg_3prime` (interior baseline)

These are the interior C→T (5′) and G→A (3′) rates used as the baseline `b`
in the exponential fit (the field referred to as "fit baseline" in some
downstream tools). They are *not* damage measurements — they reflect the
ambient sequencing-error and sequence-composition rate. A baseline of
0.002–0.005 is typical for clean Illumina data; baselines above ~0.01 suggest
elevated sequencing error, mapping artefacts, or — in the context of an
aDNA workflow — a substantial fraction of modern contamination diluting the
ancient signal. Use `bg_5prime` and `bg_3prime` as a modern-contamination
proxy alongside the d_max values: high `d_max` with low baseline is a clean
ancient call, while moderate `d_max` with elevated baseline is suspicious.

## 2. Library type classification (DS vs SS)

The `library_type` field takes values `DOUBLE_STRANDED`, `SINGLE_STRANDED`, or
`UNKNOWN`. The classifier examines which terminal sub-channels (ct5, ga3,
ct3, ga0) are active under the seven-model BIC selection described in
[Damage types](damage-types.md#libraryresolution-block).

In a **double-stranded** library (NEBNext, TruSeq), strand-overhang
deamination produces C→T at the 5′ end and G→A at the 3′ end. The 3′ end
shows G→A rather than C→T because the complementary strand of a damaged
C-overhang carries a G, and that strand is what is sequenced from the 3′
direction.

In a **single-stranded** library (Gansauge & Meyer, Santa Cruz), both ends
of the original strand survive intact through circularisation, so the
underlying C→T deamination is visible at *both* the 5′ and 3′ ends.
Consequently `d_max_5prime` and `d_max_3prime` should be of similar magnitude
in an SS library, whereas a DS library only shows the 5′ C→T directly
(the 3′ value reports G→A on the opposite strand).

`UNKNOWN` means no sub-channel beat the null model. This can happen
legitimately in modern libraries with no damage, but in a sample that should
be ancient it is a warning. Check `is_valid()` first: if `n_reads < 1000`
the classifier has not had enough evidence to commit. Below that threshold
the call is unreliable and downstream tools should treat the type as missing
rather than as `UNKNOWN`.

The flag `library_type_auto_detected` distinguishes a classifier call from a
user-forced override. If you have set `forced_library_type` on the profile,
honour it; otherwise inspect the auto-detected value.

## 3. Damage validation through cross-channel evidence

libtaph cross-validates Channel A (terminal C→T excess) against Channel B
(stop-codon-context C→T) to separate genuine deamination from
composition-driven enrichment.

- `damage_validated == true` — Channel B independently confirms the
  Channel A signal. The terminal enrichment is consistent with cytosine
  deamination acting in stop-codon-bearing contexts.
- `damage_artifact == true` — Channel A enrichment is present, but Channel B
  contradicts it. The terminal excess is driven by sequence composition
  (often GC-rich genomes or mapping artefacts), not deamination. **A sample
  with `damage_artifact == true` should not be reported as ancient on the
  strength of `d_max` alone.** Re-examine read lengths, mapping quality,
  and the F/G/H channels.
- Neither flag set — there is too little Channel B coverage to corroborate
  or refute, or the Channel A signal is too weak to test. The result is
  inconclusive and downstream interpretation rests on Channel A alone.

`composition_bias_5prime` / `composition_bias_3prime` fire when the control
channel rises in lockstep with the damage channel as the read terminus is
approached. That co-elevation is the signature of a compositional artefact
(GC-rich material driving both). When set, the corresponding terminal
damage call is unreliable even if `d_max` looks reasonable.

`position_0_artifact_5prime` / `position_0_artifact_3prime` indicate a
depleted position-0 with downstream enrichment, the typical fingerprint of
adapter or ligation residue. Genuine deamination in a DS library produces a
peak *at* position 0, not a depletion there; a depleted p0 with rising p1
through p4 is almost always a clipping or ligation artefact. When this flag
is set, libtaph automatically shifts the fit window inward and (in DS
libraries with `flag_hex_artifact`) reports a hexamer-corrected
`d5_hexamer_corrected` extrapolated back to position 0.

## 4. Oxidative damage channels (F, G, H)

Channels F, G, and H detect oxidative damage through terminal trinucleotide
context z-scores against the interior baseline. All three are reported as
binomial z-scores (`channel_f_z`, `channel_g_z`, `channel_h_z`).

A practical scale for the 5′ z-scores:

| `channel_*_z` | Interpretation |
|---|---|
| ~0 | No detectable enrichment over interior |
| 1–2 | Borderline; coverage-limited noise |
| 3–5 | Detectable signal — the conventional `z > 3` threshold |
| 5–10 | Strong signal — typical of oxidatively damaged ancient material |
| > 10 | Very strong, but check for compositional drivers |

A z-score of 5 for Channel F combined with Channel H near zero is consistent
with GC-rich bacterial DNA driving the C-context channels through
composition rather than oxidation. A z-score of 5 for F *and* H indicates
the enrichment extends to AT-context channels too, which composition alone
cannot explain.

**The co-elevation pattern F+ G+ H+ is the diagnostic signature of authentic
oxidative aDNA.** All three channels rise together because 8-oxoG-mediated
miscoding affects many trinucleotide contexts. Channel H (A→T) is
particularly informative: its pre-contexts are AT-rich (AAA, AAG, AGA), so
unlike F and G it cannot be inflated by GC-rich composition. A pattern of
**F+ G+ with H ≈ 0** is therefore diagnostic of compositional bias, most
commonly GC-rich bacterial contamination (e.g. *Pseudomonas*) rather than
genuine oxidative damage on the host material.

For Channel H specifically, prefer `channel_h_z_p2plus` over `channel_h_z`
when the position-0 background A→T rate is depressed. The p2plus variant
restricts the terminal window to positions 2–4 and avoids dilution by the
low p0 rate; downstream consumers commonly use
`max(channel_h_z, channel_h_z_p2plus)`.

The per-channel `*_uniformity_ratio` (`ca_uniformity_ratio`,
`cg_uniformity_ratio`, `at_uniformity_ratio`) is the terminal-to-baseline
ratio. A ratio above 1 indicates terminal enrichment over interior; below 1
indicates depletion. Values near 1 mean no terminal effect and should be
read together with the corresponding z-score before drawing a conclusion.

The 3′ counterparts (`channel_f3_z`, `channel_g3_z`, `channel_h3_z`,
together with their `*_3prime` uniformity ratios) are typically weaker than
their 5′ counterparts. This is partly because 3′ end-repair processes mask
some terminal positions, and partly because adapter trimming on the 3′ side
removes some of the terminal context. Treat 3′ z-scores below the 5′ values
by a factor of 2–3 as normal; only worry when the 3′ signal is *higher*
than the 5′ signal, which often signals an inversion or trimming artefact.

## 5. Channel C/D — 8-oxoG context diagnostics

Channels C and D summarise the strand-asymmetry of G→T versus C→A at
interior positions, the canonical 8-oxoG signature.

`ox_gt_uniformity` is the terminal-to-interior ratio of the G→T rate.
Values **near 1.0** indicate that the G→T enrichment is interior-dominated:
the signal is interpretable as 8-oxoG rather than as a terminal deamination
artefact. Values **substantially above 1.0** mean the G→T rate spikes at
the terminus, which is more consistent with a terminal compositional
effect than with bulk oxidative damage.

`ox_gt_asymmetry` is the strand contrast `ox_gt_baseline − ox_ca_baseline`.
**Positive values** (G→T > C→A) are the expected sign for genuine 8-oxoG,
because oxidation of guanine on one strand is read as G→T and has no
corresponding C→A counterpart on the other strand of a pooled DS library.
**Negative values** point the other way and usually indicate a sequencing
or composition artefact rather than oxidation.

`s_oxog_16ctx` reports the per-context G→T asymmetry across the 16
trinucleotide contexts NGN. It is informative when coverage is high
(reported as non-`null` when the per-context coverage exceeds 500). A
broad enrichment across most contexts is consistent with genuine oxidative
damage. A spike confined to one or two contexts (typically GG-flanked) is
characteristic of modern oxidative chemistry rather than diffuse aDNA
oxidation; it can also be noise when coverage is borderline. When
`s_oxog_16ctx` is `null` or only a couple of contexts are populated, treat
the panel as uninformative and rely on the F/G/H z-scores instead.

## 6. Channel E — depurination and AP-site fragmentation

Channel E captures the AP-site-mediated fragmentation pattern: ancient DNA
preferentially breaks 3′ of purine residues (A and G), so reads from
fragmented aDNA have purines enriched at the position immediately 5′ of the
fragmentation site.

`purine_enrichment_5prime` is the headline statistic — the excess of A+G in
the −1 / 5′-flanking position over what a uniform composition would predict.

- **Positive enrichment** (e.g. > 0.05) is the expected aDNA signature;
  values above 0.15 are strong evidence of AP-site-mediated fragmentation.
- **Near-zero enrichment** is consistent with random shearing, USER
  treatment, or modern DNA where AP-site chemistry has not had time to
  operate.
- **Negative enrichment** is unusual and typically indicates a
  sonication-driven or otherwise non-biological fragmentation profile.

Channel E is independent of the deamination channels. A sample can show
strong fragmentation enrichment with weak `d_max` (USER-treated aDNA, or
material where deamination has saturated and lost its terminal asymmetry).
In the dominant-process rule (see below), a strong Channel E with weak
deamination is what produces the `fragmentation_bias` label.

## 7. Damage-context profile and `dominant_process`

The dominant-process classifier emits a single label summarising what
process best describes the library, alongside six normalised scores in
`[0, 1]` (or `NaN` when unevaluable). The label is assigned by a
deterministic rule top-to-bottom; the first match wins. See
[Methods](methods.md#damage-context-profile) for the exact thresholds.

| Label | Implication |
|---|---|
| `cytosine_deamination` | The default ancient-DNA call: terminal C→T deamination dominates with no contradictory evidence. |
| `cpg_enriched_deamination` | Deamination is amplified at CpG sites; consistent with vertebrate methylation. Suggests host-mammalian aDNA. |
| `dipyrimidine_biased` | Damage is enriched at CC/TC dipyrimidine contexts; a UV-photoproduct signature seen in surface-exposed material (skin, sun-exposed bone fragments). |
| `oxidative_like` | Oxidative-channel score is high relative to terminal deamination. May reflect long-term storage oxidation, post-mortem oxygen exposure, or — if F+G+ but H≈0 — bacterial-composition contamination masquerading as oxidation. |
| `fragmentation_bias` | Strong purine enrichment with weak terminal deamination. Often seen in chemically-treated libraries or very old material where deamination has saturated. |
| `library_artifact_likely` | Adapter / hexamer / position-0 evidence is present *and* terminal deamination is below the threshold for a damage call. The library is unreliable; investigate trimming and re-prep. |
| `low_damage` | Terminal deamination is essentially zero. Either modern DNA, a USER-treated library, or contamination so heavy that ancient signal is diluted. |
| `none` | Insufficient reads (< 1000) or no finite `d_max` at either end. The classifier abstains. |

The six scores let downstream tools re-normalise without rescanning. A
library with `terminal_deamination_score = 0.85` and
`library_artifact_score = 0.05` is a confident damage call; the reverse
ratio is a confident artefact call.

## 8. GC-stratified mixture model

For samples that are mixtures of ancient and modern DNA, the per-read
exponential summary in Channel A is a population average and can mislead.
The two-component mixture over the GC × length cell grid disentangles them.

- `mixture_pi_ancient` — the fraction of C-sites in the high-damage
  component. A value near 1.0 means almost all reads carry damage; near
  0.0 means almost none do; intermediate values indicate a true mixture.
- `mixture_d_ancient` — the expected damage rate among the reads assigned
  to the ancient component (i.e. the damage rate conditional on being
  ancient). This is the number to compare against `d_max` from clean
  libraries.

When the mixture model and `d_max` agree, the library is uniform (ancient
or modern; no mixing). When `mixture_pi_ancient` is intermediate (say
0.2 – 0.6) and `mixture_d_ancient` is high while `d_max_combined` is
moderate, the library is a mixture: there is a real ancient component, but
it is diluted by modern reads. In that case `mixture_d_ancient` is the
better number to report for the ancient fraction, and `mixture_pi_ancient`
is the proxy for the contamination level.

## 9. Ancient-fraction deamination (fqdup integration)

When fqdup runs the per-read LLR scorer fused into the oxoG pass, it
populates a set of fields that unmix the ancient and modern read pools and
estimate deamination for each fraction independently. These are only available
when `damaged_fraction_valid == true`.

### `damaged_fraction_pi` — ancient fraction

The posterior-weighted fraction of reads classified as ancient. This is the
fraction to compare against metaDMG "damaged fraction" values. A value near 1
means almost all reads are ancient; near 0 means the ancient signal is diluted
by a large modern background.

The estimate uses a single posterior-weighted E-step: π is initialised from the
mixture identity (bulk_d5 / d_anc), per-read posteriors are computed as
`sigmoid(LLR + log(π/(1−π)))` under a two-class model with reads conditionally
independent, and π is then re-estimated as the mean posterior. This is **not**
iterated to convergence (full EM); the one-step estimate is generally not a fixed
point of the EM operator. In practice the initialisation is close enough that a
single step suffices, but the result can differ from a converged estimate when π
is far from 0.5. Compared to hard LLR > 0 thresholding, the posterior weighting
avoids the selection-on-outcome floor that pins hard-call pi at ~0.25–0.35 in
low-endogenous libraries.

### `damaged_fraction_d5` / `d3` — mixture-identity estimate

These divide the bulk `d_max` by `pi`:

```
d_anc = bulk_d_max / pi
```

derived from `bulk_d = π × d_anc + (1−π) × d_mod`, which reduces to the above
only if `d_mod ≈ 0` (modern reads carry no deamination). **This assumption is
not always met.** If the modern component carries background deamination
(d_mod > 0), then bulk_d/π **overestimates** d_anc — the bias is directional
(upward), not just noisy, and it grows as π decreases. At low endogenous
fractions (<5%), even a small d_mod (e.g. 0.01) inflates the estimate
substantially.

Additionally, the identity weights components by C-site count, not read count.
`pi` is a read fraction; these coincide only if ancient and modern reads have
similar GC content and length. If ancient reads are GC-poorer or shorter, the
mixture identity mis-scales d_anc. The independent fit fields (`*_fit`) avoid
both assumptions.

### `damaged_fraction_d5_fit` / `lambda5` — independent fitted estimate

These fit an exponential decay directly to the 15-position per-read profile
of reads classified as ancient, using log-linear regression (OLS on
log-transformed excess):

```
log(T/TC(p) − bg) = log(d_max) − λ × p
```

fitted over positions p = 0..14 where TC(p) ≥ 50 and excess > 0. The fit stops
at the first position whose coverage falls below 50 (regardless of later
positions); individual positions with non-positive excess before that point are
skipped. OLS in log space assumes homoscedastic log-residuals, which binomial
sampling does not guarantee — positions near the background threshold carry the
most noise and disproportionate leverage on λ, so the fitted λ is approximate.

This does not depend on the d_mod ≈ 0 or C-site/read-fraction assumptions of
the mixture-identity estimate and gives both `d_max` and `λ` directly from the
ancient-classified read pool. It is the preferred estimate when `pi` is uncertain
or when you need the shape of the decay.

Under d_mod ≈ 0 and clean classification, d5_fit typically meets or exceeds d5
(mixture identity). The inequality can reverse if the modern pool carries
background deamination (inflating bulk_d/π) or if classification leakage at very
low π dilutes the ancient pool and pulls d5_fit down.

### `modern_fraction_*` — modern component

The same fields for reads classified as modern. For a clean ancient library
these will be near zero. A non-zero `modern_fraction_d5_fit` can indicate
modern contamination carrying genuine (but low-level) background deamination,
or — at very low ancient fractions — false-positive classification leaking
ancient reads into the modern pool.

### Reading the trio together

| Signal | Interpretation |
|---|---|
| `pi ≈ 1`, `d5_fit ≈ d5` | High endogenous; both estimates converge toward bulk d_max — use either, provided both use the same background/normalization convention |
| `pi` moderate (0.2–0.7), `d5_fit > d5` | Mixed library under d_mod ≈ 0; `d5_fit` better describes the ancient component |
| `d5_fit` and `lambda5` both finite, `modern_d5_fit ≈ 0` | Likely clean separation — but at very low π, classification leakage can make `modern_d5_fit` unreliable as a cleanliness indicator |
| `valid == false` | Paired-mode input, bulk d_max < 0.01, or insufficient reads — field unavailable |

---

## 10. QC flags and when to distrust results

Several flags should override an otherwise positive damage call.

- `terminal_inversion` — terminal damage rate is *below* interior at one
  or both ends. The fit has been forced into a fallback configuration.
  Likely causes: read-orientation bug upstream, mis-trimmed adapters, or a
  reverse-complemented BAM. Investigate before reporting any `d_max` from
  this profile.
- `hexamer_excess_tc` — a **library chemistry / priming bias** metric, not an
  adapter remnant signal (see `flag_hex_artifact` and `adapter_clipped` for
  those). It measures T/(T+C) at read starts minus T/(T+C) in the read interior
  across all matched C/T hexamer pairs. A negative value means the library prep
  preferentially produces C-starting read starts (e.g. ligation junction
  composition in SS protocols), enriching C-context trinucleotides in the
  channels F and G pre-count denominator across the full 5′ terminal window.
  The effect is false **negatives** in F/G (suppression toward large negative
  z-scores), not false positives in Channel A. Interpret F/G z-scores with
  caution when `hexamer_excess_tc < −0.05`; the raw value is always reported.
- `short_read_frac` — the fraction of reads below the short-read threshold.
  High values indicate library-prep size selection failure or extreme
  fragmentation; both reduce the reliability of terminal statistics
  because end-effects dominate the per-read information budget.
- `flag_hex_artifact` / `adapter_clipped` / `adapter3_clipped` — direct
  evidence of adapter remnants. The hexamer-corrected `d5` is computed
  for DS libraries when these conditions co-occur with a position-0
  artefact, but the safer call is to re-trim and re-run.
- `is_detection_unreliable()` — true when terminal inversion or
  composition-bias flags fire. **Treat every damage call from such a
  profile with caution**; report values only with their flag context.

## 11. Common patterns in practice

| Pattern | Likely interpretation |
|---|---|
| High `d_max`, `damage_validated`, F+G+H+ all elevated, DS library | Authentic ancient DNA with concurrent oxidative damage. The standard "good ancient sample" call. |
| High `d_max`, `damage_validated`, F+G+ with H ≈ 0, DS library | Authentic ancient DNA with minimal oxidative damage; or, in unsupervised data, GC-rich-composition contamination — disambiguate with `damage_validated` and host coverage. |
| High `d_max` in a negative control with F+G+ H ≈ 0 | GC-rich bacterial contamination driving Channel A through composition. Not ancient damage. |
| Moderate `d_max` with `damage_artifact == true` | Composition artefact — Channel B contradicts Channel A. Report as unreliable. |
| `d_max ≈ 0` and `dominant_process == low_damage` | Modern DNA, USER-treated library, or contamination so heavy that ancient signal is diluted. |
| SS library with `d_max_5prime ≈ d_max_3prime` | Single-stranded library showing the expected symmetric C→T at both ends. |
| `dominant_process == dipyrimidine_biased` | UV-driven photoproduct signature. Common in skin, hair, or surface-exposed bone fragments. |
| `dominant_process == fragmentation_bias` | AP-site-mediated fragmentation dominates over terminal deamination. Typical of USER-treated or heavily chemically processed libraries. |
| `terminal_inversion == true` with otherwise high `d_max` | Likely orientation or trimming bug upstream — investigate before trusting the call. |
| Mixture model: `mixture_pi_ancient = 0.3`, `mixture_d_ancient = 0.18`, `d_max_combined = 0.06` | A genuine ancient signal diluted by ~70% modern DNA (0.3 × 0.18 ≈ 0.054 ≈ 0.06, assuming d_mod ≈ 0). Note: `mixture_pi_ancient` is the fraction of **C-sites** in the high-damage component, not necessarily the read fraction — they differ when ancient and modern reads have unequal GC or length. Report `mixture_d_ancient` as the ancient damage rate. |

When in doubt, look at all three layers — the d_max value, the
cross-channel validation flags, and the dominant-process label — and
require them to agree. libtaph is built so that no single number can be
the whole story.
