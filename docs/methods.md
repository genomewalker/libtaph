# Methods

libtaph estimates ancient-DNA damage directly from raw reads, without alignment or a reference genome. Working directly on raw reads makes the method applicable early in a workflow, before any alignment step. The cost is that terminal base composition, library preparation, and genuine post-mortem damage can all produce superficially similar sequence patterns, making it difficult to attribute a terminal signal to any one cause. The method therefore uses multiple independent signals to reduce that ambiguity — each targeting a different biochemical process — and validates agreement between them before reporting a damage estimate.

The processing pipeline has four stages: per-position accumulation, joint damage validation, library-type classification, GC-stratified mixture modelling, and final `d_max` selection. The sections below follow that order. The emphasis is not only on what is computed, but on why each stage is needed and what ambiguity it resolves.

## Per-position accumulation and baseline estimation

The raw material for all downstream inference is a set of terminal and interior counts accumulated directly from reads. For each read, libtaph inspects the first and last 15 positions and records the base combinations relevant to deamination and its controls. The goal at this stage is deliberately simple: count what is observed at the termini, count what is observed away from the termini, and postpone interpretation until enough evidence has accumulated.

For each terminal position $p$:

- **5' damage channel (ct5)**: $r_p = T_p / (T_p + C_p)$, excess T where C is expected
- **3' damage channel (ga3)**: $r_p = A_p / (A_p + G_p)$, excess A where G is expected
- **5' control channel**: $A_p / (A_p + G_p)$ at the 5' end, which should remain flat if the 5' C→T signal is genuine deamination rather than composition bias
- **3' control channel (ct3)**: $T_p / (T_p + C_p)$ at the 3' end, which acts as a negative control in DS libraries but becomes a biologically meaningful SS signal when the original strand is sequenced

The interior of the read provides the undamaged reference state. libtaph uses the middle third of each read (positions 30 to $L-30$) to estimate the background rate $b$. This region is far enough from the termini to be largely insensitive to terminal overhang damage, but still reflects the sample's intrinsic sequence composition. In a reference-free setting this interior baseline is essential: the method does not ask whether a specific read differs from a reference genome, but whether terminal positions deviate from that sample's own interior composition.

These per-position ratios are later summarized with a weighted least-squares (WLS) amplitude estimate using a fixed decay rate $\hat\lambda$:

$$\hat{A} = \max\!\left(0,\; \frac{\sum_p n_p \cdot e^{-\hat\lambda p} \cdot (r_p - b)}{\sum_p n_p \cdot e^{-2\hat\lambda p}}\right)$$

where $n_p$ is the coverage at terminal position $p$. This expression asks whether terminal positions that should carry the strongest damage signal are systematically elevated above the interior baseline, weighted by how many reads cover each position. The WLS estimate is not the whole method; it is one summary of a richer terminal pattern that is cross-checked later against composition controls and codon-based evidence.

---

## Paired-end mode

`update_sample_profile_pe(R1, R2)` consumes a read pair directly rather than
running SE accumulation on each mate separately. The mapping is biological:
`R1[i]` measures the top strand 5'-end at position `i`, while `R2[i]`'s
complement measures the **same fragment's** top strand 3'-end at position `i`.
Treating R2 as an independent SE read would attribute its 5'-end signal to a
second 5' damage process that does not exist, and would let G→A events on the
bottom strand contribute to ct5 rather than ga3.

When the insert length `M` is shorter than the read length `L`, both mates
read into the sequencing adapter. Adapter bases have a fixed composition that
differs from the genomic background, and after complement-mapping they imprint
a deterministic pattern on the terminal channels — in the worst case enough to
reroute a DS library through the SS branch of the classifier. libtaph detects
short-insert pairs by R1/R2 overlap (a 15 bp window with at most 3 mismatches)
and skips them, recording the count under `pe_short_insert_skipped`. Pairs
with `M ≥ L` enter the accumulator unchanged.

---

## Damage model

Once terminal and interior counts have been accumulated, libtaph models ancient-DNA deamination as a terminal process that decays inward from each end of the fragment. In double-stranded libraries this appears as C→T excess at the 5' end and G→A excess at the 3' end, because the complementary strand carries the same damaged cytosines as G→A when read in the opposite orientation. This characteristic pattern was first systematically described by [Briggs et al. (2007)](https://www.pnas.org/doi/10.1073/pnas.0704665104) and subsequently formalized as an exponential decay model by [Jónsson et al. (2013)](https://doi.org/10.1093/bioinformatics/btt193).

The terminal signal is modelled as

$$\delta(p) = b + A \cdot e^{-\lambda p}$$

where:

- $p$ is distance from the terminus (0-indexed)
- $b$ is the interior background rate
- $A$ is the terminal excess above that background
- $\lambda$ is the decay constant that controls how rapidly the excess disappears with distance from the end

The role of this equation is not just descriptive. It encodes the central biological assumption of the method: genuine terminal damage should be strongest at the edge of the fragment and should decay smoothly inward, whereas many confounders either remain flat or fail to reproduce the same decay in independent channels.

To report a quantity comparable to [mapDamage2.0 (Jónsson et al. 2013)](https://doi.org/10.1093/bioinformatics/btt193) and [metaDMG (Michelsen et al. 2022)](https://doi.org/10.1101/2022.12.06.519264), libtaph converts the fitted amplitude into a calibrated terminal damage rate,

$$D_{\max} = \frac{A}{1 - b}$$

which estimates the fraction of terminal C sites that were converted to T, after accounting for the background T proportion. $A$ is the absolute excess above that background; $D_{\max}$ converts it to the fraction of damage-susceptible terminal cytosines that were actually damaged.

The background $b$ used in the fit is **tail-anchored**: it is estimated from
positions 20..49 of the terminal accumulator rather than from the global mean
of the per-position rates. This matters because the per-position rates near
the terminus carry the deamination signal we are trying to fit, so a global
mean leaks signal into the baseline and biases $A$ downward. Tail positions
20..49 sit far enough inward that the exponential damage component has decayed
by 1–2 orders of magnitude (for typical $\lambda$), so they are effectively a
chemistry-aware baseline immune to the deamination signal.

Alongside $D_{\max}$ libtaph reports an **area-excess** statistic per channel:
the sum of $(r_p - b)$ over the first 10 positions, together with a
likelihood-ratio score against the bg-only null. Area-excess is a non-parametric
companion to $D_{\max}$: it does not depend on the exponential functional form
and is robust when the decay is slower or noisier than the model assumes. The
classifier consumes both, so a library that fits the model poorly but has
clear terminal excess can still be classified correctly.

---

## Multi-channel validation

The first major inference step in `finalize_sample_profile` asks a question that is specific to reference-free damage estimation: is terminal T enrichment actually ancient-DNA damage, or is it simply a compositional feature of the reads? A single nucleotide-frequency channel cannot answer that reliably, because elevated terminal T/(T+C) can arise from real C→T deamination or from a terminal sequence bias unrelated to damage.

libtaph addresses this by evaluating several independent biochemical channels. Some channels are primary evidence for deamination, while others act as controls or orthogonal diagnostics.

| Channel | Signal | Notes |
|---------|--------|-------|
| A | C→T rate per position | Primary deamination channel |
| B | Stop codon conversion (CAA/CAG/CGA → TAA/TAG/TGA) at 5' | Sequence-composition-independent validation |
| B₃′ | Stop codon conversion via G→A at 3' (TGG → TAG/TGA) | Validates SS 3' damage |
| C | G→T transversions (8-oxoG) via stop codon uniformity | Uniform across read; distinguishes oxidation from terminal deamination |
| D | G→T / C→A direct transversion rates | Terminal/interior contrast for oxidative damage |
| E | Purine enrichment at 5' termini (depurination) | Independent evidence of fragmentation at AP sites |

In practice, the decision about whether 5' deamination is real is driven by a joint model built from Channel A, its 5' control channel, and Channel B. The other channels remain informative diagnostics, especially in asymmetric or unusual libraries, but they do not replace the core logic: a valid deamination call should be supported by an independent signal that composition alone cannot mimic. Channels D and E are described in detail in [Damage types](damage-types.md); the biochemical basis for Channel D (8-oxoG G→T misincorporation) is reviewed in [Shibutani et al. (1991)](#ref-shibutani1991) and [Neeley & Essigmann (2006)](#ref-neeley2006), and for Channel E (depurination fragmentation and AP-site strand cleavage) in [Lindahl (1993)](#ref-lindahl1993).

At a high level, `damage_validated` means the primary nucleotide signal is supported by composition-robust evidence. Conversely, `damage_artifact` means terminal enrichment is present in Channel A but contradicted by the stop-codon signal, so the observed excess is more plausibly explained by composition than by genuine damage.

### Channel B: frame-agnostic codon scanning

Channel B is designed to answer exactly the question that Channel A cannot: if terminal C→T damage is real, do we also see the specific codon changes that such damage should create? The key convertible codons are CAA, CAG, and CGA, which become the stop codons TAA, TAG, and TGA after a single C→T event at the first position.

Crucially, this test does **not** require reference alignment, ORF prediction, strand selection, or choosing a single reading frame. Instead, libtaph scans all three forward reading frames simultaneously at every position in every read:

- For each frame offset $f \in \{0, 1, 2\}$ and each codon starting at position $p$, classify the codon as:
  - **pre-image** (CAA, CAG, CGA): a codon that would become a stop after C→T damage
  - **stop** (TAA, TAG, TGA): either a genuine stop or a damage-converted pre-image
- Accumulate those counts separately at each terminal distance (positions 0-14) and in the interior (positions 30 to $L-30$)

Why does this work without frame selection? Because the method does not need to know which individual codons are biologically translated. It only needs the **relative excess** of stop codons near the terminus compared with the interior. Any bias caused by the mixture of frames is present at both locations, because the same three-frame scan is applied everywhere. As a result, the ratio

$$r_p = \frac{\text{stop}}{\text{pre} + \text{stop}}$$

is largely insensitive to frame composition itself. What can change this ratio specifically at the termini is real C→T damage converting pre-image codons into stops. This is why Channel B is useful as a composition-robust validator rather than merely another way of counting terminal bases.

The trade-off is statistical rather than conceptual. Channel B is more specific than Channel A, but it uses a smaller subset of informative sites, so it needs more data before it becomes decisive. In low-complexity or low-coverage samples, Channel A may show a suggestive terminal pattern while Channel B remains underpowered.

### Joint damage model (`JointDamageModel`)

Rather than interpreting Channel A, the control channel, and Channel B separately, libtaph fits them jointly with `JointDamageModel::fit()`. The observed 5' terminal profile can be a mixture of at least two effects.

1. **Real deamination**, which should elevate Channel A and Channel B but not the control channel.
2. **Compositional artifact**, which can elevate Channel A and the control channel together, but has no reason to elevate Channel B.

The joint model makes that distinction explicit:

| Signal | Formula | What it captures |
|--------|---------|-----------------|
| Channel A (ct5) | $\pi_{TC}(p) = (b_{tc} + a(p)) + (1 - b_{tc} - a(p)) \cdot \delta_{\max} \cdot e^{-\lambda p}$ | Deamination + compositional artifact |
| Control (5' AG) | $\pi_{AG}(p) = b_{ag} + a(p)$ | Compositional artifact only |
| Channel B (stop) | $\pi_{stop}(p) = b_{stop} + (1 - b_{stop}) \cdot \delta_{\max} \cdot e^{-\lambda p}$ | Deamination only |

Here $b_{tc}$, $b_{ag}$, and $b_{stop}$ are interior baselines estimated directly from the accumulated counts. The term

$$a(p) = a_{\max} \cdot e^{-\lambda p}$$

represents a terminal compositional shift that affects Channel A and the control in parallel. The damage term, parameterized by $\delta_{\max}$, affects Channel A and Channel B. This structure is what makes the model identifiable in practice:

- If the apparent terminal excess is due to composition, Channel A and the control rise together and $a_{\max}$ absorbs that pattern, while Channel B remains flat.
- If the excess is due to genuine deamination, Channel A rises, the control stays comparatively flat, and Channel B rises in parallel with the damage term.

This is the core reason the joint model is more robust than thresholding Channel A alone. It does not merely ask whether the terminal frequency is high; it asks which of two mechanistic explanations better accounts for the joint behavior of all three observed channels.

The fitted model has three free parameters: $\delta_{\max}$, $\lambda$, and $a_{\max}$. The baselines are fixed from the interior counts. The implementation searches a grid over $(\lambda, \delta_{\max})$ and, at each grid point, optimizes $a_{\max}$ by a nested golden-section search. It then compares the best damage model ($M_1$: $\delta_{\max} > 0$) to a no-damage model ($M_0$: $\delta_{\max} = 0$) using both BIC and an approximate Bayes factor.

`damage_validated` is set when the joint damage model is favored strongly enough by the implemented criteria: posterior probability $p_{\text{damage}} > 0.95$, BIC evidence $\Delta\text{BIC} > 10$ ([Kass & Raftery 1995](#ref-kass1995)), or a degenerate Bayes factor in extreme cases. `damage_artifact` is set when Channel A shows terminal enrichment but the stop-codon channel remains flat or inverted, implying that the enrichment is not supported by composition-independent evidence. When `damage_artifact = true`, `d_max_combined` is forced to zero so that downstream tools do not treat a compositional bias as authentic damage.

---

## Library-type classifier

After establishing whether terminal damage is genuine, libtaph asks a different question: **which library architecture produced the observed terminal asymmetry?** This is a separate inference problem. Damage validation is about distinguishing authentic deamination from artifact; library classification is about determining which ends and strands are expected to carry that damage signal.

This separation matters. A sample can have real damage but still show very different terminal patterns depending on whether the library is double-stranded, single-stranded in the original orientation, single-stranded in the complement orientation, or a mixture. The classifier therefore works with four observable channels rather than only the damage-validation trio.

See [Damage types](damage-types.md) for the biological background of these categories and how fqdup uses the resulting call for damage-aware hashing.

### Biological basis

| Library prep | Typical observed pattern |
|-------------|--------------------------|
| Double-stranded (DS) | Symmetric ct5 and ga3 exponential decay |
| DS + end-repair artifact | DS pattern plus a ga0 spike |
| Single-stranded (SS), complement orientation | ga0-dominated 3' G→A signal; often without smooth ga3 |
| SS, original orientation | ct5 + ct3; no ga3 |
| SS, mixed orientations | ct5 + ga0; weak residual ga3 may be present depending on protocol |

Single-stranded library protocols include [Gansauge & Meyer (2013)](https://doi.org/10.1038/nprot.2013.038), the Santa Cruz method ([Kapp et al. 2021](#ref-kapp2021)), and SRSLY (Claret Bioscience; [Gansauge et al. 2017](#ref-gansauge2017)).

The main difficulty is that these patterns are not perfectly separable by any single statistic. For example, a strong `ga0` spike can reflect either a DS end-repair artifact or a genuinely single-stranded complement-orientation library. Likewise, asymmetric mixtures of SS orientations can partially mimic a DS signal if one end dominates the evidence. The classifier therefore uses a model-comparison approach rather than a fixed decision tree.

### Four channels and the GA0 spike

In addition to the smooth exponential channels, the classifier treats position 0 of the 3' end as a separate single-position feature (**ga0**). This is necessary because some library-preparation artifacts and some SS signatures are localized almost entirely at the ligation junction rather than forming a smooth inward decay.

The `ga0` term is particularly informative because it captures two biologically distinct cases:

- **DS end-repair artifacts**: bilateral behavior, with elevated terminal substitutions caused by enzymatic processing rather than unilateral SS damage
- **SS complement-orientation reads**: unilateral 3' G→A spikes without the 5' counterpart expected under DS symmetry

Treating `ga0` as its own channel prevents the classifier from forcing a spike into an exponential model that does not actually describe the data.

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

In `M_DS_symm`, ct5 and ga3 share one amplitude parameter (`joint↑`). This reflects a direct consequence of Chargaff's first rule: in double-stranded DNA, [C] = [G] between strands by complementarity ([Chargaff et al. 1950](#ref-chargaff1950)). Deaminated cytosines on the original strand appear as C→T when that strand is sequenced; the same damage events on the complementary strand appear as G→A when it is sequenced. With unbiased adapter ligation capturing both strands in equal proportion, as expected in standard DS library preparation, the two terminal signals should have equal amplitude. The shared parameter improves efficiency when that assumption holds, but it also deliberately penalizes a DS interpretation when the ends are strongly asymmetric — in effect, the model is using symmetry itself as evidence. The post-hoc symmetry check below makes that penalty explicit when asymmetry is extreme.

Each composite BIC is the sum of channel-level BIC terms:

$$\text{BIC}(M) = \sum_{\text{channels}} \text{BIC}_{\text{component}}(\text{channel}, \text{active/null})$$

Lower BIC wins ([Schwarz 1978](#ref-schwarz1978)). The null model contributes no amplitude parameter; each active channel contributes one amplitude parameter with a complexity penalty of $\ln(n_{\text{trials}})$. This formulation makes the trade-off explicit: a more complex library explanation is only preferred if it improves fit enough to justify the additional flexibility.

In addition to the waterfall candidates, the implementation also evaluates an unconstrained four-channel SS model (`M_SS_full`: ct5=alt, ga3=alt, ga0=alt, ct3=alt), reported as `library_bic_mix` in `SampleDamageProfile`. This model is useful diagnostically and for likelihood-based summaries, even though it is not itself the winner in the cascade.

### Cascade

The classifier evaluates candidate models in the following order:

1. Start with `best = BIC(M_bias)`, type = DS
2. M_DS_symm, DS
3. M_DS_spike (only when `spike_is_ss = false`), DS
4. M_DS_symm_art, DS
5. M_SS_comp, SS
6. M_SS_orig (only when ct3 ΔBIC > 0), SS
7. M_DS_spike as SS (only when `spike_is_ss = true`), SS
8. M_SS_asym (only when `spike_is_ss = true`), SS

Here `spike_is_ss = (ga0.amplitude >= 0.10)`. A very large `ga0` spike is difficult to reconcile with a simple DS end-repair artifact, so once it crosses that threshold the classifier explicitly allows single-stranded interpretations that would otherwise be disfavored.

The cascade structure is a practical compromise between exhaustive model comparison and biological prior knowledge. It keeps the classifier interpretable while still allowing asymmetric or spike-dominated libraries to escape an overly rigid DS default.

### Post-hoc symmetry check

If a DS model wins but the evidence is strongly asymmetric between ct5 and ga3, the symmetry assumption underlying the DS call is revisited:

```
if DS wins AND ga3.ΔBIC > 30,000 AND ct5.ΔBIC / ga3.ΔBIC < 0.50 -> SS
```

This rule is motivated by a known failure mode: an SS library can contain a mixture of original-orientation and complement-orientation reads such that one end dominates ct5 and the other dominates ga3, creating a superficially DS-like pattern without true bilateral symmetry. The post-hoc check prevents the shared-amplitude DS model from winning purely because both ends have signal, even when that signal is clearly imbalanced.

### Rescue rules

The rescue rules exist because real libraries do not always occupy cleanly separated regions of model space. In particular, a localized `ga0` spike can dominate the likelihood and pull the classifier toward a DS artifact explanation even when the underlying biology is unilateral SS damage.

**M_DS_spike rescue** (for `spike_is_ss = false`): when `M_DS_spike` wins but `d5` is effectively absent, the 3' spike is more consistent with complement-orientation SS than with bilateral DS end repair:

```
if DS wins AND ds_spike_won AND ct5.ΔBIC <= 0 AND ga3.ΔBIC <= 0
   AND ga0.ΔBIC > 0 AND ga0.amplitude > 0.02 AND d5 <= 0.005 -> SS
```

**GA0 bilateral rescue** (for `spike_is_ss = true`): when a very large `ga0` spike drives `M_DS_symm_art` to win, `d5` is used to distinguish true bilateral DS behavior from complement-only SS behavior:

```
if DS wins AND spike_is_ss AND ga0.ΔBIC > 0 AND d5 <= 0.005 -> SS
```

These rescue rules are not ad hoc patches in the pejorative sense; they encode biological asymmetries that the base model family only approximates. Without them, the classifier would be too willing to interpret unilateral junction spikes as DS artifacts simply because those spikes are statistically strong. The trade-off is that the rules introduce a small amount of thresholded expert knowledge, but they do so precisely to avoid a larger systematic bias.

Validated on 24 DS controls with `ga0_amp >= 0.10`, all had `d5 >= 0.11`. In contrast, the observed SS mixed failures with `ga0_amp >= 0.10` had `d5 = 0.00`. That separation is large enough to support a deterministic rescue threshold.

### UNKNOWN category

If no model beats `M_bias` exactly, the library has no detectable damage in any channel and the result is `UNKNOWN` rather than a default DS call. This is an important design choice: in a genuinely low-damage or modern library, the correct answer is often that library type cannot be inferred from sequence alone. Reporting `UNKNOWN` makes that uncertainty explicit instead of converting absence of evidence into a potentially misleading structural call.

---

## GC-stratified estimation

After the global damage call and library-type classification are established, libtaph asks whether the sample is compositionally heterogeneous. This matters because environmental and metagenomic ancient-DNA datasets often contain mixtures of genuinely ancient molecules and lower-damage background DNA, and those components can differ systematically in GC content. A single global `d_max` is still useful, but it can blur together distinct populations.

To expose that heterogeneity, reads are binned by **interior** GC content into 10 bins spanning 0-100%. Interior GC is used rather than whole-read GC so that the binning itself is less sensitive to terminal damage. Within each bin, the method re-estimates damage using terminal and interior counts aggregated across reads in that bin.

This introduces a deliberate trade-off. Stratifying by GC can reveal mixtures that a global estimate would average away, but each bin has fewer observations than the full sample. The implementation therefore regularizes low-information bins and treats bins with insufficient support as invalid rather than forcing unstable estimates.

Two mixture models operate over the GC bins:

**`DamageMixtureModel`** is a simple 2-component Gaussian model over per-bin `d_max` values. It treats one component as effectively undamaged ($\mu = 0$, $\sigma = 0.01$) and one as damaged ($\mu$ estimated, $\sigma = 0.10$), and uses the EM algorithm ([Dempster, Laird & Rubin 1977](#ref-dempster1977)) to estimate the fraction and mean damage of the damaged component. This model provides a compact summary of whether the GC bins separate into low-damage and high-damage groups.

**`MixtureDamageModel`** is the full GC-aware model. It fits $K = 2, 3, 4$ latent classes with multiple random restarts and selects the best solution by BIC. Each class has its own `delta_max`, while `lambda` and `a_max` are shared. The model also learns a GC distribution for each class, so it can represent situations in which high-damage and low-damage molecules occupy different GC ranges rather than simply different damage levels.

In implementation terms, the full model operates on per-bin "super-reads": for each GC bin it aggregates Channel A, the control channel, and Channel B counts, then evaluates how well each latent class explains those summaries. Per-class `delta_max` values are updated by grid search, shared `a_max` is updated by golden-section search, and per-bin baselines are shrunk toward global baselines to stabilize sparse bins. This regularization is important: without it, a few low-count bins could dominate the mixture fit through sampling noise alone.

The reported outputs summarize different biological questions:

- `mixture_pi_ancient`: fraction of C-sites assigned to high-damage classes (`delta_max > 5%`)
- `mixture_d_ancient`: expected damage rate among the high-damage component only
- `mixture_d_population`: population-average damage across all C-sites
- `mixture_d_reference`: expected damage restricted to GC >= 50% classes, intended as a rough proxy for the GC-rich fraction that reference-based tools often emphasize
- `mixture_K`: number of classes selected by BIC

The full GC-stratified model is therefore not just a refinement of the global call. It is an attempt to separate "how damaged are the ancient molecules?" from "what fraction of the sample is ancient at all?" Those are distinct questions in mixed environmental samples, and collapsing them into one global `d_max` can be misleading.

---

## D_max combination

The final step converts several end-specific and channel-specific estimates into a single reported `d_max_combined`. A naive strategy would simply average the 5' and 3' estimates, but that is not robust across library types. In DS libraries, averaging is sensible because both ends should reflect the same terminal damage process. In SS, mixed-orientation, or artifact-affected libraries, however, one end can be informative while the other is structurally absent, biased by adapter effects, or driven by a different biochemical process.

libtaph therefore uses an asymmetry-aware combination rule rather than unconditional averaging. The selected strategy is recorded in `d_max_source`.

| Condition | Strategy |
|-----------|----------|
| Both ends valid, low asymmetry | Average of d5 and d3 |
| High asymmetry (SS mixed) | Use whichever end is more reliable |
| Only one end valid | Use that end alone |
| Channel B₃′ available | Structural estimate from stop codon decay |

The rationale is as follows:

- **When both ends agree**, averaging reduces variance and yields the most stable estimate.
- **When one end is clearly less reliable**, averaging would dilute the biologically meaningful signal with structural noise or protocol-specific asymmetry.
- **When only one end is supported**, forcing a two-end summary would create a false sense of symmetry.
- **When Channel B₃′ is informative**, its stop-codon conversion signal can provide a structural 3' estimate even when a simple terminal nucleotide fit is weak or confounded.

This final combination stage is where several earlier diagnostics matter operationally. Position-0 artifacts, inverted terminal patterns, strong library asymmetry, and Channel B/B₃′ structural estimates all influence which end should be trusted. The output `d_max_combined` is therefore not merely the result of one exponential fit; it is the end product of the full evidence cascade described above.

---

## Hexamer signals and domain classifier

Ancient-DNA libraries rarely fail in a uniform way. Adapter stubs, ligation artifacts, and residual protocol sequence all manifest as over-represented 6-mers at the read terminus, often with upstream bases that damage alone would not explain. libtaph accumulates a full 4096-bin hexamer histogram at both the 5' terminus (`hexamer_count_5prime[4096]`, hexamers at read positions 0-5) and the interior (`hexamer_count_interior[4096]`, middle third of the read), plus running totals `n_hexamers_5prime` and `n_hexamers_interior`. These counts drive three independent analyses.

**Terminal-vs-interior enrichment.** `compute_hex_enriched_5prime` reports every hexamer whose terminal frequency exceeds its interior frequency by `log2fc > 1.5` (3x enrichment) with a minimum of 20 observations on each side. The function also flags whether the hexamer starts with T, which is the expected leading base for a C-to-T deamination product. An enriched terminal hexamer whose first base is not T is suspicious because genuine deamination cannot produce arbitrary 5' base changes.

**Adapter stub detection.** `detect_adapter_stubs` keeps up to five enriched 5' stubs (`log2fc > 3.0`, first base != T) and up to five enriched 3' stubs (`log2fc > 3.0`, last base != A, computed from the pre-scan 3' terminal hexamer histogram passed in as `hex3_terminal[4096]` with total `n_hex3`). 3' detection is gated on 5' stubs being present first, reflecting the biology: if the 5' end is clean, a smaller 3' enrichment is unlikely to be a real adapter. The boolean `flag_hex_artifact` fires when the single top-enriched 5' hexamer has `log2fc > 1.5` and a first base other than T. Genuine deamination cannot lead with any base except T, so a non-T top hexamer above that threshold is evidence of a composition or protocol artifact rather than damage.

**Hexamer distribution statistics.** `compute_hex_stats` summarises the terminal vs interior 6-mer distributions with Shannon entropy of each, Jensen-Shannon divergence (`jsd`), and a multinomial G-test on the two empirical distributions. The statistic is reported alongside a z-score against a chi-squared null and a p-value (`shift_p`). Large `jsd` values mean the terminal composition is not merely noisy: it is drawn from a different distribution than the interior, which is the expected fingerprint of a protocol-driven terminal bias.

**Hexamer-corrected 5' d_max.** In DS libraries where an adapter artifact suppresses the position-0 C-to-T rate, the WLS exponential fit is shifted inward: `fit_offset_5prime >= 1` indicates the fit started at position 1 or later. When `position_0_artifact_5prime && flag_hex_artifact && !adapter_clipped && lambda_5prime > 0`, the preservation summary extrapolates the decay back to position 0 via `d5_corr = d5_raw * exp(lambda_5prime * fit_offset_5prime)`, capped at 0.50 so a fit on a small offset cannot blow up. The correction is recorded in `d5_hexamer_corrected` and `d5_was_corrected` for auditability; the raw value is never overwritten.

**Domain-aware hexamer scoring.** libtaph ships eight pre-computed 4096-bin hexamer frequency tables covering GTDB (bacteria + archaea), fungi, protozoa, invertebrate, plant, vertebrate mammalian, vertebrate other, and viral sequences. `score_all_domains(seq, frame)` sums `log2(freq[h] * 4096 + epsilon)` across every in-frame hexamer in a read for each domain, then normalises to a posterior with temperature-scaled softmax (`exp(0.1 * (score - max))`). The result is an eight-element `MultiDomainResult` with per-domain probabilities, a `best_domain`, and a `best_score`. Two mechanisms consume that posterior:

- `set_active_domain(d)` sets a global selected domain (read via `get_active_domain()`). It controls which table `get_hexamer_freq(code)` and `get_hexamer_score(seq)` fall through to when the caller does not specify a domain. The active domain is not thread-local: set it once for the pass.
- `set_ensemble_mode(true)` and `set_domain_probs(probs)` switch hexamer lookups through `get_ensemble_hexamer_freq(code)`, which returns the posterior-weighted average of the eight tables rather than committing to one. Ensemble state and the posterior (`current_domain_probs`) are thread-local, so per-thread pipelines can carry their own blend.

The domain classifier is not a taxonomic identifier in the metagenomic sense. It provides a calibrated hexamer background for Channel B codon scanning and library-composition diagnostics that is better than a single universal model when samples are known to be non-bacterial.

---

## Handling sparse and noisy data

Every stage of the pipeline makes explicit decisions about what to do when evidence is thin, because collapsing sparse counts into a point estimate silently fabricates structure.

**Coverage thresholds.** Per-position baseline and shift estimators (`src/damage_estimation.cpp`) require `total[pos] >= MIN_COVERAGE = 100` before using a position. Positions below the threshold are skipped, not imputed. `compute_damage_mask` uses the same threshold (`min_cov = 100`) by default: positions that do not clear `min_cov` are never flagged as masked regardless of apparent rate. At the sample level, `is_valid()` returns `true` only when `n_reads >= 1000`; consumers that want a conservative gate check `is_detection_unreliable()` which also folds in the inversion and composition-bias flags.

**Shrinkage baselines.** In the GC-stratified fits, per-bin baselines are shrunk toward the global baseline before any likelihood is evaluated. Low-count bins therefore pull toward the overall sample rather than producing a noise-driven `d_max`, while high-count bins dominate their own posterior. The GC categorical component of `MixtureDamageModel` applies Dirichlet smoothing with `alpha = 1` (`p_gc[k][b] = w[b][k] * c_sites + 1`), which prevents a single unrepresented GC bin from zeroing out a class's likelihood.

**Per-cell invalidity markers.** `fit_length_gc_joint_mixture` marks cells with insufficient C-sites (`c_sites < 1`) or that were never valid after `finalize_all` with `cell_w_ancient[b][g] = -1.0`. Downstream code distinguishes "low ancient fraction" (value near 0) from "no evidence" (value -1), so consumers do not silently treat missing cells as uncontaminated.

**Uninformative priors as floors.** Several late-stage signals fall back to explicit priors rather than aggressive extrapolation. Preservation-layer `preservation_f_cpg` defaults to `0.3` (uninformative) when the CpG signal is absent or inverted, which prevents a dropped signal from being interpreted as evidence of no CpG methylation. The sample-level `preservation_f_cpg` fallback fires for `total < 1000`, returning `0.40`. Joint damage model priors follow the same philosophy: when the posterior collapses to the no-damage model, `p_damage = BF / (1 + BF)` under a `0.5` prior, so the classifier reports "not detected" rather than "no damage" when the data are uninformative.

**Validation flags and rescue rules.** The cascade in `finalize_sample_profile` is explicitly layered to avoid false positives from composition or artifacts while still rescuing partially valid signals:

- `composition_bias_5prime` / `composition_bias_3prime` fire when the control channel rises in lockstep with the damage channel, meaning the apparent excess is compositional. These flags feed `is_detection_unreliable()`.
- `position_0_artifact_5prime` / `position_0_artifact_3prime` flag a depleted position 0 with downstream enrichment, which is the signature of an adapter / ligation artifact rather than a terminal decay. When this fires together with `flag_hex_artifact` in a DS library, the hexamer-corrected `d5` is computed as described above.
- `inverted_pattern_5prime` / `inverted_pattern_3prime` flag terminal rates below the interior; the inversion fallback in the library classifier is suppressed when Channel B has validated damage, because the stop-codon signal is treated as ground truth when available.
- `damage_artifact = true` forces `d_max_combined = 0.0` so downstream consumers never apply an artifact-driven mask.

**Numerical guards.** Inputs near the probability boundary are clamped before use (`std::clamp(d, 1e-3, 1.0 - 1e-3)` in the per-cell posterior for the length x GC fit; `std::max(S_FLOOR, sigma_binom)` with `S_FLOOR = 0.02` for per-cell variances; `std::max(p, 1e-300)` before `-log10(p)`). Log-sum-exp is computed with the standard max-subtract trick so posterior weights remain finite when either component has very low likelihood.

The combined effect is that under-powered samples do not produce confident-looking point estimates: they produce flagged estimates, quantized p-values, or explicit fallbacks. Every path that ends with a number records why that number is believed, and the preservation summary distinguishes evidence strength (`_evidence`) from the estimate itself (`_eff`).

---

## Length-stratified estimation

Aggregating every read into one profile treats short, damage-rich molecules and longer, less-damaged molecules as a single population. Ancient samples are usually mixtures along the length axis, so libtaph can stratify the profile by read length.

`taph::LengthBinStats` holds up to four independent `SampleDamageProfile` instances plus an `edges` vector that defines the length bins. A read of length `L` is routed to the bin whose half-open interval contains `L`. The same accumulate / merge / finalize contract as the unstratified profile applies per bin:

```cpp
taph::LengthBinStats stats;
stats.configure({35, 60, 90});          // bins: [min,35), [35,60), [60,90), [90,max]
for (auto& read : reads)
    stats.update(read.seq, read.length);
stats.finalize_all();                    // finalizes every bin's SampleDamageProfile
```

Bin edges can be supplied explicitly or derived from the empirical length distribution. `taph::detect_log_length_gmm_edges` fits a 1D Gaussian mixture on a log-length histogram, selects the number of components by BIC up to `max_components`, and returns the density valleys between adjacent component means as split points. Separation filters are applied before the valleys are computed: the histogram must contain at least 256 counts, components with weight below `3%` are dropped, and adjacent means within `0.25` in log space (roughly 30% length change) are merged. If fewer than two components survive, `n_components` is reported as 1 and no edges are returned. `taph::detect_quantile_length_edges` is a separate helper that splits the histogram into equal-count quantile bins; libtaph does not fall through to it automatically, the caller chooses when to use it.

On top of a finalized `LengthBinStats`, `fit_length_gc_joint_mixture` runs a shared-component 2-Gaussian mixture over the flat length × GC cell grid. The damaged component has a single `d_ancient` shared across all cells; the undamaged component is fixed at `μ=0`. Each cell contributes its `d_max` and `c_sites`. The fit returns:

- `d_ancient`, `pi_ancient`: shared mean and mixing weight of the ancient component.
- `d_population`: `c_sites`-weighted mean `d_max` across all cells.
- `cell_w_ancient[b][g]`: posterior `P(damaged | cell)` for each (length bin, GC bin). Cells with insufficient C-sites are marked `-1.0`.
- `converged`, `separated`: numerical state of the EM fit.

The design reuses `DamageMixtureModel`'s fixed priors (undamaged width `τ₀ = 0.01`, damaged width `τ₁ = 0.10`, minimum per-cell sigma `0.02`) so length × GC cells live in the same calibration as the GC-only fit. Per-cell posteriors give downstream code a per-length-bin ancient fraction rather than a single global number, which matters when the short end of the distribution is almost entirely ancient and the long end is dominated by a modern contaminant.

---

## Context-aware terminal damage

C→T rates at the 5' terminus average over four upstream bases (A, C, G, T). Different biological processes imprint different context distributions, and a single `d_max` erases that structure. libtaph accumulates the 5' C→T channel separately for each upstream base and exposes fitted amplitudes per context.

The counts live on `SampleDamageProfile`:

- `ct5_t_by_upstream[ctx][pos]`, `ct5_total_by_upstream[ctx][pos]` for `ctx ∈ {CTX_AC, CTX_CC, CTX_GC, CTX_TC}` and `pos ∈ [0, N_POS)`.
- `ct5_t_interior_by_upstream[ctx]`, `ct5_total_interior_by_upstream[ctx]` as per-context interior baselines.

Finalization fits a fixed-lambda binomial likelihood per context via golden-section search on `d` (`fit_ct5_ctx_amplitude`), using the interior T/(T+C) as the baseline and clamping `mu` to `[1e-6, 1-1e-6]`. Coverage gates are explicit: the interior must have `total >= 500` with effective coverage `>= 100`, each fit position must have `>= 50` observations, at least three positions must clear that gate, and the cumulative terminal effective coverage must be `>= 50`. A fit that lands at `d >= 0.98` is reported as invalid (the optimizer hit the upper wall and the signal is indistinguishable from saturation). Results are stored as `dmax_ct5_by_upstream[ctx]`, `baseline_ct5_by_upstream[ctx]`, and the terminal / interior coverage arrays. Two derived contrasts summarize the pattern:

- `dipyr_contrast = ½(d_CC + d_TC) − ½(d_AC + d_GC)`: positive when dipyrimidine contexts (CC, TC) carry more damage than non-dipyrimidine, as expected for UV-driven 6-4 photoproduct / CPD-linked deamination.
- `cpg_contrast = d_GC − ⅓(d_AC + d_CC + d_TC)`: positive when 5-methylcytosine at CpG sites dominates, as expected for methylation-driven hydrolytic deamination.

A coverage-weighted pseudo-chi-squared statistic on the four fitted `dmax` values tests the null that the context amplitudes are equal. With three degrees of freedom, `context_heterogeneity_chi2 > 7.81` sets `context_heterogeneity_detected = true`. The p-value is read off a coarse critical-value ladder (`0.001, 0.01, 0.05, 0.1, 0.5, 0.9`) rather than a full chi-squared CDF, since the statistic is a pseudo-chi-squared and a precise p-value would imply more calibration than the construction warrants.

In practice the contrast pair separates three regimes we observe in real samples: ancient DNA with a flat spectrum across contexts (cytosine hydrolysis everywhere), methylation-driven deamination with a CpG-dominant signal (`cpg_contrast > 0`, other contexts near baseline), and UV-damaged material where dipyrimidine contexts lead (`dipyr_contrast > 0`).

---

## Damage-context profile

The damage-context profile is a training-free, reference-free, alignment-free summary of the per-process signals already present in `SampleDamageProfile`. It is a pure function of the finalized profile plus the already-computed `cpg_z` and `hex_shift_z`; it introduces no new scan pass, no fitted model, and no external reference panel.

Six scores are emitted in `[0, 1]`, with `NaN` when the underlying signal is not evaluable. The normalization constants are tunable and live as named `constexpr` values at the top of `library_interpretation.cpp`.

| Score | Definition | Underlying fields |
|---|---|---|
| `terminal_deamination_score` | `1 − exp(−max(d_max_5, d_max_3) / 0.10)` | `d_max_5prime`, `d_max_3prime` |
| `cpg_context_score` | `sigmoid(cpg_z)` from `compute_cpg_score` | `log2_cpg_ratio`, `effcov_ct5_cpg_like_*` |
| `dipyrimidine_context_score` | `clamp(dp.dipyr_contrast / 0.05, 0, 1)` where `dipyr_contrast = 0.5·(d_CC + d_TC) − 0.5·(d_AC + d_GC)` | `dipyr_contrast`, `dmax_ct5_by_upstream[AC,CC,GC,TC]` |
| `oxidative_context_score` | `clamp(max(|ox_gt_asymmetry|, mean(s_oxog_16ctx)) / 0.05, 0, 1)` | `ox_gt_asymmetry`, `s_oxog_16ctx` |
| `fragmentation_context_score` | `clamp(purine_enrichment_5prime / 0.15, 0, 1)` | `purine_enrichment_5prime` |
| `library_artifact_score` | `max(indicator(flag_hex_artifact ∨ adapter_clipped ∨ adapter3_clipped ∨ pos0 artifact ∨ fit_offset_{5,3}prime > 1), sigmoid(hex_shift_z − 4))` | `flag_hex_artifact`, `adapter_clipped`, `hex_shift_z`, `position_0_artifact_*`, `fit_offset_{5,3}prime` |

A single `dominant_process` label is assigned by a deterministic rule over the six scores. The rule is evaluated top-to-bottom and stops at the first match:

1. `n_reads < 1000` → `none` (insufficient coverage). The six scores are populated where the underlying signals are evaluable; fields whose source signal is `NaN` remain `NaN` in the output.
2. Terminal deamination score `NaN` (neither end has a finite `d_max`) → `none`.
3. A boolean artifact flag is set (`flag_hex_artifact`, `adapter_clipped`, `adapter3_clipped`, a position-0 artifact on either end, or `fit_offset_{5,3}prime > 1`) **and** `terminal_deamination_score < 0.5` → `library_artifact_likely`. The categorical label fires only when adapter/hexamer evidence is present *and* genuine deamination is not dominating; strong terminal damage alongside adapter contamination keeps a damage label (the artifact booleans remain visible in `evidence`). A high `hex_shift_z` alone is reported in the score but cannot trigger the label on its own, because clean libraries can reach z ≈ 10–20 purely from compositional variance.
4. `fragmentation_context_score > 0.5` and `terminal_deamination_score < 0.3` → `fragmentation_bias`.
5. `terminal_deamination_score < 0.10` → `low_damage`.
6. `cpg_context_score > 0.7`, `terminal_deamination_score > 0.3`, and `log2_cpg_ratio > 0.15` → `cpg_enriched_deamination`. The effect-size floor on `log2_cpg_ratio` prevents large-sample z-scores from assigning the label on negligible CpG vs non-CpG differences.
7. `oxidative_context_score > 0.5` and `terminal_deamination_score < 0.5` → `oxidative_like`.
8. `dipyrimidine_context_score > 0.4` → `dipyrimidine_biased`.
9. Otherwise → `cytosine_deamination`.

Threshold comparisons require a finite score; `NaN` scores never satisfy a `>` / `<` branch. Only a `NaN` `terminal_deamination_score` forces `none`; if terminal deamination is finite and all other scores are `NaN`, the rule falls through to `cytosine_deamination` as the catch-all.

The `evidence` block in the JSON output mirrors the raw underlying numbers (d_max, λ, log2 CpG ratio, dipyr contrast, `ox_gt_asymmetry`, `s_oxog_{mean,max}`, purine enrichment, `hex_shift_z`, adapter and position-0 flags, `fit_offset_{5,3}prime`, `n_reads`). Downstream tools can therefore re-normalize scores or replace the rule without rescanning.

**Design notes.** The six scores only summarise mechanisms that are independently measurable in the existing `SampleDamageProfile`. The profile deliberately avoids latent-process decomposition (for example `NMF` over a trinucleotide matrix against a reference panel), since ancient-DNA damage mechanisms are small in number, well characterised, and already directly observable here. The rule is intentionally legible so a human can check which signals drove each label.

---

## References

<a id="ref-briggs2007"></a>**Briggs AW, Stenzel U, Johnson PLF, Green RE, Kelso J, Prüfer K, Meyer M, Krause J, Ronan MT, Lachmann M, Pääbo S** (2007) Patterns of damage in genomic DNA sequences from a Neandertal. *Proc Natl Acad Sci USA* **104**:14616–14621. [doi:10.1073/pnas.0704665104](https://doi.org/10.1073/pnas.0704665104)

<a id="ref-chargaff1950"></a>**Chargaff E, Zamenhof S, Green C** (1950) Composition of human desoxypentose nucleic acid. *Nature* **165**:756–757. [doi:10.1038/165756b0](https://doi.org/10.1038/165756b0)

<a id="ref-dempster1977"></a>**Dempster AP, Laird NM, Rubin DB** (1977) Maximum likelihood from incomplete data via the EM algorithm. *J R Stat Soc Ser B* **39**:1–22. [doi:10.1111/j.2517-6161.1977.tb01600.x](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x)

<a id="ref-gansauge2013"></a>**Gansauge M-T, Meyer M** (2013) Single-stranded DNA library preparation for the sequencing of ancient or damaged DNA. *Nat Protoc* **8**:737–748. [doi:10.1038/nprot.2013.038](https://doi.org/10.1038/nprot.2013.038)

<a id="ref-gansauge2017"></a>**Gansauge M-T, Gerber T, Glocke I, Korlević P, Lippik L, Nagel S, Riehl LM, Schmidt A, Meyer M** (2017) Single-stranded DNA library preparation from highly degraded DNA using T4 DNA ligase. *Nucleic Acids Res* **45**:e79. [doi:10.1093/nar/gkx033](https://doi.org/10.1093/nar/gkx033)

<a id="ref-jonsson2013"></a>**Jónsson H, Ginolhac A, Schubert M, Johnson PLF, Orlando L** (2013) mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters. *Bioinformatics* **29**:1682–1684. [doi:10.1093/bioinformatics/btt193](https://doi.org/10.1093/bioinformatics/btt193)

<a id="ref-kapp2021"></a>**Kapp JD, Green RE, Shapiro B** (2021) A fast and efficient single-stranded genomic library preparation method optimized for ancient DNA. *J Hered* **112**:241–249. [doi:10.1093/jhered/esab012](https://doi.org/10.1093/jhered/esab012)

<a id="ref-kass1995"></a>**Kass RE, Raftery AE** (1995) Bayes factors. *J Am Stat Assoc* **90**:773–795. [doi:10.1080/01621459.1995.10476572](https://doi.org/10.1080/01621459.1995.10476572)

<a id="ref-lindahl1993"></a>**Lindahl T** (1993) Instability and decay of the primary structure of DNA. *Nature* **362**:709–715. [doi:10.1038/362709a0](https://doi.org/10.1038/362709a0)

<a id="ref-michelsen2022"></a>**Michelsen C, Fortunato G, Warinner C, Jonsson H, Schroeder H, Orlando L, Renaud G, Librado P** (2022) metaDMG: a fast and accurate ancient DNA damage toolkit for metagenomic sequencing data. *bioRxiv*. [doi:10.1101/2022.12.06.519264](https://doi.org/10.1101/2022.12.06.519264)

<a id="ref-neeley2006"></a>**Neeley WL, Essigmann JM** (2006) Mechanisms of formation, genotoxicity, and mutation of guanine oxidation products. *Chem Res Toxicol* **19**:491–505. [doi:10.1021/tx0600043](https://doi.org/10.1021/tx0600043)

<a id="ref-schwarz1978"></a>**Schwarz G** (1978) Estimating the dimension of a model. *Ann Stat* **6**:461–464. [doi:10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136)

<a id="ref-shibutani1991"></a>**Shibutani S, Takeshita M, Grollman AP** (1991) Insertion of specific bases during DNA synthesis past the oxidation-damaged base 8-oxodG. *Nature* **349**:431–434. [doi:10.1038/349431a0](https://doi.org/10.1038/349431a0)
