---
type: libtaph-json-field
title: channels_fgh
tier: standard
estimand: oxidative terminal enrichment channels F (C→A), G (C→G), H (A→T)
stability: stable
emitted_by: profile_json.cpp
---

### Channels F, G, and H — complement-asymmetry oxidative channels

Channels F, G, and H detect oxidative damage through **terminal trinucleotide-context
enrichment**: the rate of a substitution type at the 5′ or 3′ terminal window (positions
0–4, adapter-adjusted) is compared to an interior reference using a binomial z-score.

**Counting model.** For each channel, reads at each position are classified into
two mutually exclusive categories: `stop` (the trinucleotide represents the substituted
outcome, e.g. TAA/TAG/TGA for C→A in Channel F) and `pre` (convertible source contexts
that would produce `stop` under the target substitution, e.g. TCA/TCG/TAC/TGC for
Channel F). The event rate in any window is `stop / (pre + stop)`.

**Detection statistic.** For a terminal window accumulating `k_t` stop events out of
`n_t = pre_t + stop_t` total, and an interior reference with `k_i` out of `n_i`:

```
p_pool = (k_t + k_i) / (n_t + n_i)
z = (k_t/n_t − k_i/n_i) / sqrt(p_pool·(1−p_pool)·(1/n_t + 1/n_i))
```

This is a descriptive two-proportion z-statistic, not a fully calibrated hypothesis
test. Read counts within a sample are correlated (duplicate molecules, taxon admixture,
trinucleotide composition biases), so z-scores of 100+ are possible in heavily damaged or
contaminated samples and should not be interpreted as independent-site p-values.

**Interior reference resolution.** The far-interior region (position ≥ 30, minimum 50
counts) is preferred; when insufficient, the mid-read region (positions 5–14) is used
as fallback. This avoids leakage of terminal damage signal into the baseline.

**Adapter skip (p0 gate).** Position 0 is excluded from the terminal window when a
position-0 composition artifact is detected (`position_0_artifact_5prime` at the 5′
end, `position_0_artifact_3prime` at the 3′ end). The 5′ flag fires when the C→T rate
at p=0 falls below p=1 — the signature of an adapter ligation junction suppressing
the true damage signal at the terminal nucleotide. The 3′ flag fires when the G→A
spike at p=0 is isolated without smooth decay through p=1–4. All three channels (F,
G, H) apply the same gate on both ends; no channel applies a TC-hexamer gate at p=0.

| Variable | Trigger condition | Channels affected |
|----------|-------------------|-------------------|
| `p0_tc_5` | `position_0_artifact_5prime` | F, G (C-context) |
| `p0_h5`   | `position_0_artifact_5prime` | H (A/T-context) |
| `p0_3`    | `position_0_artifact_3prime` | F, G, H (3′ end) |

**Adapter remnants and TC hexamer bias are tracked by two separate systems.**

*Adapter remnant detection* uses `position_0_artifact_5prime/3prime` (sharp p=0 spike
without downstream decay) and hexamer-based stub enrichment (`adapter_clipped`,
`adapter3_clipped`, `flag_hex_artifact` — a non-T-starting hexamer enriched > 3× over
interior at the terminal). `inverted_pattern_5prime/3prime` fires when the full terminal
window is depressed below interior. When adapter flags fire, Channel A shifts its
exponential fit start point (`fit_offset_5prime`) to skip the suppressed positions and
extrapolates back to p=0 for reporting.

*TC hexamer bias* (`hexamer_excess_tc`) is a **library chemistry signal, not an adapter
remnant**. It computes T/(T+C) at read starts minus T/(T+C) in the interior, across all
matched C/T hexamer pairs. A negative value means the library prep preferentially produces
C-starting read starts (e.g. ligation junction composition in SS protocols), not that
adapter sequence is present in the reads. The two systems fire on different evidence and
may fire independently or together.

**TC hexamer bias and F/G z-scores.** TC hexamer priming (`hexamer_excess_tc < 0`)
enriches C-context trinucleotides (TCA/TCG/TAC/TGC) across the full terminal window
(positions 0–4), increasing the pre-context denominator without a proportional increase
in stop contexts. This suppresses the terminal stop rate relative to interior, driving
F and G z-scores **negative**. TC hexamer bias therefore cannot produce false-positive
F/G signals; it only suppresses genuine signal. The `hexamer_excess_tc` value is
reported so consumers can flag or down-weight F/G z-scores when strongly negative
(empirically, values below −0.05 are sufficient to drive F_z below −100 at n > 500M reads).

**Channel F — C→A terminal enrichment (8-oxoG bottom-strand complement).** Pre-contexts
are TCA/TCG/TAC/TGC; stop contexts are TAA/TAG/TGA. Elevated terminal C→A rate indicates
8-oxoG on the opposite strand read as C→A in the sequenced strand.

**Channel G — C→G terminal enrichment (empirical oxidative-context proxy).**
Pre-contexts are TCA/TAC; stop contexts are TGA (from TCA) and TAG (from TAC).
Further-oxidized guanine products such as guanidinohydantoin (Gh) and
spiroiminodihydantoin (Sp) provide one possible biochemical route to C→G-like
readouts [[Henderson et al. 2002](#references); [Neeley & Essigmann 2006](#references)],
but the reference-free statistic itself is the observed terminal C→G enrichment,
not a lesion-specific call.

**Channel H — A→T terminal enrichment (empirical; mechanism uncertain).** Pre-contexts
are AAA/AAG/AGA; stop contexts are TAA/TAG/TGA. The biological basis has not been
conclusively established: 8-oxoadenine (8-oxoA), the most studied adenine oxidation
product, primarily causes A→C rather than A→T transversions [[Kamiya et al. 1995](#references)],
so the lesion driving the observed terminal A→T enrichment in ancient DNA remains to be
identified. Position 0 carries a lower background A→T rate than positions 2–4
independently of any artifact flag; accordingly a secondary score `channel_h_z_p2plus`
is computed from positions 2–4 only. Detection fires when either score exceeds the
threshold.

**Distinguishing genuine aDNA oxidation from composition bias.** Genuine ancient
oxidative damage produces co-elevation of F, G, and H. A sample showing F+/G+ with H
near zero is more consistent with GC-rich bacterial contamination (e.g. *Burkholderia*,
*Pseudomonas*) driving composition asymmetry than with authentic oxidative aDNA damage.
The mechanism is sequence-compositional: reads from high-GC organisms carry a higher
density of C-context trinucleotides (TCA/TCG/TAC/TGC), which are the pre- and
stop-contexts for channels F and G. When reads are fragmented preferentially at certain
positions, the terminal trinucleotide composition differs from the interior in a way that
cannot be normalised away by the binomial test, inflating z-scores for F and G. A/T
context trinucleotides (AAA/AAG/AGA) are unaffected by GC-content asymmetry, so Channel
H remains near zero.

| Field | Description |
|-------|-------------|
| `channel_f_valid` | `true` when ≥ 200 total trinucleotide counts for baseline |
| `channel_f_z` | Binomial z-score, 5′ C→A terminal vs interior; negative = depletion |
| `channel_f3_valid` | `true` when 3′ end has sufficient counts |
| `channel_g_valid` | `true` when ≥ 200 total counts for G baseline |
| `channel_g_z` | Pooled binomial z-score, 5′ C→G terminal vs interior. **Legacy/descriptive** (`channel_g_z_inference=descriptive_not_calibrated_p_value`); the MH z below is the corrected statistic |
| `channel_g_mh_z` | **Corrected** context-stratified Mantel-Haenszel z for G (2 strata: TCA→TGA, TAC→TAG), mirroring F's stratification. Removes terminal context-composition bias |
| `channel_g_common_or` | MH common odds ratio for G |
| `channel_g3_valid` | `true` when 3′ end has sufficient counts |
| `channel_h_valid` | `true` when ≥ 200 total counts for H baseline |
| `channel_h_z` | Pooled binomial z-score, 5′ A→T terminal (positions p0_h5–4) vs interior. **Legacy/descriptive**; the MH z below is the corrected statistic |
| `channel_h_mh_z` | **Corrected** context-stratified MH z for H (3 strata: AAA→TAA, AAG→TAG, AGA→TGA) |
| `channel_h_common_or` | MH common odds ratio for H |
| `channel_h_z_p2plus` | Same as `channel_h_z` but terminal window restricted to positions 2–4; more robust because position 0 carries a lower background A→T rate that dilutes the signal when included |
| `channel_h3_valid` | `true` when 3′ end has sufficient counts |
| `fgh_adapter_prefixes_excluded` | Number of 5′ hexamer prefixes excluded from the F/G/H terminal counts (set by `recompute_fgh_excluding_adapter_prefixes`; 0 if no adapter stubs were detected) |

libtaph stores the raw z-scores; it does not apply a detection threshold internally.
Downstream consumers commonly use z > 3.0 as a heuristic. For Channel H, take
`max(channel_h_z, channel_h_z_p2plus)` to avoid the position-0 dilution effect.

**Pooled z is descriptive, not a calibrated p-value.** The pooled `channel_{f,g,h}_z`
scale ~√N on correlated reads and carry `channel_{f,g,h}_z_inference =
descriptive_not_calibrated_p_value`. The context-stratified Mantel-Haenszel z
(`channel_{f,g,h}_mh_z`) is the corrected inferential statistic; F already used it, and
G and H now extend the same per-context stratification (estimands and assumptions are
sourced from the channel registry and emitted in the `damage_types` legend as `estimand`
/ `assumptions`). The emitter's F/G/H **detection gate** requires the pooled z to agree in
sign with the corrected MH z before firing — a library whose pooled and MH z disagree in
sign (Simpson's paradox from terminal context composition) is no longer called detected.

**Adapter prefix exclusion.** During read accumulation each terminal F/G/H codon count
is bucketed by the read's first hexamer (4096 bins). After `detect_adapter_stubs`
identifies enriched adapter hexamers, `recompute_fgh_excluding_adapter_prefixes` re-sums
all terminal counts omitting those bins and recomputes all three z-scores (including
`channel_h_z_p2plus`) from the clean subset. This corrects for adapter remnants whose
sequence happens to contain stop codons (e.g. TGA/TAA/TAG) within the first five
nucleotides of the read, which would otherwise inflate F/G/H z-scores. T-starting
hexamers are never flagged as adapters — genuine aDNA C→T damage naturally produces
T-starting read starts. `fgh_adapter_prefixes_excluded` records how many prefixes were
excluded; 0 means the pre-scan found no non-T-starting enriched stubs.

---