# Methods

libdart-damage builds a sample profile in two phases. During `update_sample_profile`, it accumulates terminal and interior counts from raw reads. During `finalize_sample_profile`, it converts those counts into derived estimates in this order: joint damage validation, library-type classification, GC-stratified mixture modelling, and final `d_max` selection. The sections below follow that execution order while keeping the underlying mathematical model explicit.

## Per-position accumulation and baseline estimation

libdart-damage scans each read and accumulates base counts at the first and last 15 positions. For each position $p$:

- **5' damage channel (ct5)**: $r_p = T_p / (T_p + C_p)$, excess T where C expected
- **3' damage channel (ga3)**: $r_p = A_p / (A_p + G_p)$, excess A where G expected
- **5' control channel**: $A_p / (A_p + G_p)$ at 5' end, should be flat if ct5 is real
- **3' control channel (ct3)**: $T_p / (T_p + C_p)$ at 3' end, carries SS original-orientation signal

The middle third of each read (positions 30 to $L-30$) provides the background baseline $b$.

WLS amplitude estimation uses the fixed decay rate $\hat\lambda$ estimated from these accumulated terminal profiles:

$$\hat{A} = \max\!\left(0,\; \frac{\sum_p n_p \cdot e^{-\hat\lambda p} \cdot (r_p - b)}{\sum_p n_p \cdot e^{-2\hat\lambda p}}\right)$$

where $n_p$ is coverage (T+C count) at position $p$.

---

## Damage model

With terminal and interior counts in hand, libdart-damage models ancient-DNA deamination as an exponential decay from the read terminus. In double-stranded libraries this produces C→T substitutions at the 5' end and G→A substitutions at the 3' end (the complement of C→T on the opposite strand):

$$\delta(p) = b + A \cdot e^{-\lambda p}$$

where $p$ is the distance from the read terminus (0-indexed), $b$ is the background (interior) substitution rate, $A$ is the damage amplitude, and $\lambda$ is the decay constant.

D_max is the calibrated damage estimate at position 0:

$$D_{\max} = \frac{A}{1 - b}$$

This matches the definition used by metaDMG and mapDamage.

---

## Multi-channel validation

The first major inference step in `finalize_sample_profile` asks whether apparent terminal C→T enrichment is genuine damage or a compositional artefact. Conceptually, libdart-damage validates terminal C→T damage by asking whether Channel A (terminal nucleotide excess) is supported by Channel B (composition-robust stop-codon conversion). Operationally, the implementation fits a joint model over three observed 5' signals — Channel A (T/(T+C)), the 5' control channel (A/(A+G)), and Channel B (stop/(pre+stop)) — that simultaneously estimates genuine deamination (`delta_max`) and composition artefact (`a_max`).

Six biochemical channels are tracked in total:

| Channel | Signal | Notes |
|---------|--------|-------|
| A | C→T rate per position | Primary deamination channel |
| B | Stop codon conversion (CAA/CAG/CGA → TAA/TAG/TGA) at 5' | Sequence-composition-independent |
| B₃′ | Stop codon conversion via G→A at 3' (TGG → TAG/TGA) | Validates SS 3' damage |
| C | G→T transversions (8-oxoG) via stop codon uniformity | Uniform across read; non-exponential distinguishes from deamination |
| D | G→T / C→A direct transversion rates | Terminal/interior ratio confirms genuine 8-oxoG vs. artefact |
| E | Purine enrichment at 5' termini (depurination) | AP-site fragmentation; evidence of ancient origin independent of deamination |

At a high level, `damage_validated` means the primary nucleotide signal is supported by composition-robust evidence, whereas `damage_artifact` means terminal enrichment is present in Channel A but contradicted by the stop-codon signal.

### Channel B: frame-agnostic codon scanning

Channel B requires no reference alignment, no ORF prediction, and no frame or strand selection. It works by scanning all three forward reading frames simultaneously at every position in every read:

- For each frame offset $f \in \{0, 1, 2\}$ and each codon starting at position $p$, classify the codon as:
  - **pre-image** (CAA, CAG, CGA) — a codon that becomes a stop codon if C→T damage occurs
  - **stop** (TAA, TAG, TGA) — could be genuine stop or damage-converted

- Accumulate counts at each distance from the 5' terminus (positions 0–14) and in the interior (positions 30 to $L-30$) independently.

The key insight is that frame composition bias cancels in the ratio. Because the same mixture of three reading frames is scanned at both terminal and interior positions, any frame-specific enrichment of stop codons is identical at both locations. The only thing that can cause the **ratio** $r_p = \text{stop} / (\text{pre} + \text{stop})$ to be higher at terminal positions than in the interior is real C→T damage converting pre-image codons into stops.

### Joint damage model (JointDamageModel)

Rather than evaluating Channel A and Channel B independently, libdart-damage fits them jointly via `JointDamageModel::fit()`. This separates three signals that are otherwise conflated:

| Signal | Formula | What it captures |
|--------|---------|-----------------|
| Channel A (ct5) | $\pi_{TC}(p) = (b_{tc} + a(p)) + (1 - b_{tc} - a(p)) \cdot \delta_{\max} \cdot e^{-\lambda p}$ | Deamination + compositional artefact |
| Control (5' AG) | $\pi_{AG}(p) = b_{ag} + a(p)$ | Compositional artefact only |
| Channel B (stop) | $\pi_{stop}(p) = b_{stop} + (1 - b_{stop}) \cdot \delta_{\max} \cdot e^{-\lambda p}$ | Deamination only |

The artefact term $a(p) = a_{\max} \cdot e^{-\lambda p}$ is shared by Channel A and the control channel, but is **absent from Channel B**. This means:

- If T/(T+C) rises terminally because of GC composition bias, both Channel A and the control rise together and $a_{\max} > 0$ absorbs it; Channel B stays flat.
- If T/(T+C) rises because of genuine deamination, Channel A rises but the control stays flat; simultaneously Channel B rises. Only real damage drives both $\delta_{\max} > 0$ and a Channel B excess.

The fitted model has three free parameters: $\delta_{\max}$, $\lambda$, and $a_{\max}$. The background rates $b_{tc}$, $b_{ag}$, and $b_{stop}$ are estimated directly from the interior counts and then treated as fixed. The implementation searches a grid over $(\lambda, \delta_{\max})$ and uses a nested golden-section optimization for $a_{\max}$ at each grid point.

`damage_validated` is set when the joint damage model ($M_1$: $\delta_{\max} > 0$) is favored over the no-damage model ($M_0$: $\delta_{\max} = 0$) by posterior probability ($p_{\text{damage}} > 0.95$), by BIC evidence ($\Delta\text{BIC} > 10$), or by a degenerate Bayes factor.

`damage_artifact` is set when terminal enrichment is present in Channel A but the stop-codon channel (Channel B) remains flat or inverted, indicating the excess is better explained by sequence composition than by genuine deamination. When `damage_artifact = true`, `d_max_combined` is set to zero.

---

## Library-type classifier

Once damage validation has been assessed, `finalize_sample_profile` fits a separate four-channel BIC classifier to distinguish double-stranded and single-stranded library signatures. See [Damage types](damage-types.md) for a full description of the biological origin of each class, which substitution channels distinguish them, and how fqdup uses the classification for damage-aware hashing.

### Biological basis

| Library prep | Typical observed pattern |
|-------------|--------------------------|
| Double-stranded (DS) | Symmetric ct5 and ga3 exponential decay |
| DS + end-repair artifact | DS pattern plus a ga0 spike |
| Single-stranded (SS), complement orientation | ga0-dominated 3' G→A signal; often without smooth ga3 |
| SS, original orientation | ct5 + ct3; no ga3 |
| SS, mixed orientations | ct5 + ga0; weak residual ga3 may be present depending on protocol |

### Four channels and the GA0 spike

In addition to the smooth exponential channels, a single-position model fits position 0 of the 3' end (**ga0**). This captures:

- DS end-repair artefacts: bilateral, both 5' pos-0 CT and 3' pos-0 GA elevated
- SS complement-orientation reads: unilateral, 3' GA0 spike only, no 5' counterpart

### BIC models

Seven composite BIC models partition the four channels into active (alt) and inactive (null):

| Model | ct5 | ga3 | ga0 | ct3 | Interpretation |
|-------|-----|-----|-----|-----|----------------|
| M_bias | null | null | null | null | No signal |
| M_DS_symm | joint↑ | joint↑ | null | null | DS symmetric deamination |
| M_DS_spike | null | null | alt | null | DS end-repair artifact only |
| M_DS_symm_art | joint↑ | joint↑ | alt | null | DS deamination + end-repair |
| M_SS_comp | null | alt | alt | null | SS complement-orientation only |
| M_SS_orig | alt | null | null | alt | SS original-orientation only |
| M_SS_asym | alt | null | alt | null | SS both orientations, no smooth ga3 |

↑ `joint` = ct5 and ga3 share one amplitude parameter (M_DS_symm constraint). This is 1 free parameter for 2 channels; if ct5 ≠ ga3, the joint fit is penalised and a pure SS model wins.

Each composite BIC is the sum of component BICs:

$$\text{BIC}(M) = \sum_{\text{channels}} \text{BIC}_{\text{component}}(\text{channel}, \text{active/null})$$

Lower BIC wins. The null model has 0 free parameters; each active channel contributes 1 free parameter (amplitude) with penalty $\ln(n_{\text{trials}})$.

In addition to the candidate models used in the waterfall, the implementation also evaluates an unconstrained four-channel SS model (`M_SS_full`: ct5=alt, ga3=alt, ga0=alt, ct3=alt), reported as `library_bic_mix` in `SampleDamageProfile`. This score is available for diagnostic purposes and for constructing likelihood ratios but is not itself a winner in the cascade.

### Cascade

The classifier runs a waterfall over models in this order:

1. Start: `best = BIC(M_bias)`, type = DS
2. M_DS_symm, DS
3. M_DS_spike (only when `spike_is_ss = false`), DS
4. M_DS_symm_art, DS
5. M_SS_comp, SS
6. M_SS_orig (only when ct3 ΔBIC > 0), SS
7. M_DS_spike as SS (only when `spike_is_ss = true`), SS
8. M_SS_asym (only when `spike_is_ss = true`), SS

`spike_is_ss = (ga0.amplitude ≥ 0.10)`: a GA0 spike above 10% amplitude is too large to be a DS end-repair artifact and enters the SS model set.

### Post-hoc symmetry check

If DS wins but the ct5/ga3 amplitude ratio is highly asymmetric, the DS symmetry assumption is violated:

```
if DS wins AND ga3.ΔBIC > 30,000 AND ct5.ΔBIC / ga3.ΔBIC < 0.50 → SS
```

Applied when SS libraries have original-orientation reads dominating ct5 and complement reads dominating ga3, creating apparent asymmetry that breaks M_DS_symm.

### Rescue rules

**M_DS_spike rescue** (for `spike_is_ss = false`): when M_DS_spike won but d5 ≈ 0, the pos-0 spike is unilateral (complement-orientation SS), not bilateral DS end-repair:

```
if DS wins AND ds_spike_won AND ct5.ΔBIC ≤ 0 AND ga3.ΔBIC ≤ 0
   AND ga0.ΔBIC > 0 AND ga0.amplitude > 0.02 AND d5 ≤ 0.005 → SS
```

**GA0 bilateral rescue** (for `spike_is_ss = true`): when a large GA0 spike (≥ 0.10) drives M_DS_symm_art to win, d5 discriminates DS bilateral from SS complement-only:

```
if DS wins AND spike_is_ss AND ga0.ΔBIC > 0 AND d5 ≤ 0.005 → SS
```

Validated on 24 DS controls with ga0_amp ≥ 0.10: all have d5 ≥ 0.11. All SS mixed failures with ga0_amp ≥ 0.10 have d5 = 0.00. The gap is >20× the threshold.

### UNKNOWN category

If no model beats M_bias (`best == BIC(M_bias)` exactly), the library has no detectable damage in any channel. The result is `UNKNOWN` rather than a default DS call. The caller should use `--library-type ds` or provide metadata.

---

## GC-stratified estimation

After the global damage and library-type calls are established, reads are binned by their interior GC content (10 bins, 0–100%). Within each bin, damage is estimated independently, enabling separation of high-damage ancient DNA from low-damage modern contamination in mixed samples.

Two mixture models operate over GC bins:

**`DamageMixtureModel`** (2-component Gaussian, used internally by `JointDamageModel`): fits a simple undamaged (μ=0, σ=0.01) and damaged (μ estimated, σ=0.10) component over per-bin d_max values. Runs standard EM up to 50 iterations.

**`MixtureDamageModel`** (K=2–4 categorical, full GC-aware model): tries K=2, 3, 4 components with 5 random restarts each, selects by BIC. Per-component δ_max is updated by grid search (61 points, 0–0.60); shared a_max by golden-section search. Reports:

- `mixture_pi_ancient`: fraction of C-sites in high-damage components (δ > 5%)
- `mixture_d_ancient`: expected damage rate among ancient reads
- `mixture_d_population`: population-average damage rate
- `mixture_d_reference`: damage rate in GC ≥ 50% bins (metaDMG proxy)
- `mixture_K`: number of components selected by BIC

---

## D_max combination

The final step combines the end-specific estimates into `d_max_combined` using an asymmetry-aware decision rule:

| Condition | Strategy |
|-----------|----------|
| Both ends valid, low asymmetry | Average of d5 and d3 |
| High asymmetry (SS mixed) | Use whichever end is more reliable |
| Only one end valid | Use that end alone |
| Channel B₃′ available | Structural estimate from stop codon decay |
