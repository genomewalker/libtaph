---
type: libtaph-json-field
title: bulk_damage
tier: standard
estimand: reference-free bulk deamination rate (Design-B estimator)
stability: stable
emitted_by: profile_json.cpp
---

### `bulk_damage`

This block reports the reference-free, length-stratified terminal C→T damage estimator (design B). It is present whenever `attempted` is `true`; if the solver did not run (e.g. insufficient reads to fill at least two length bins), the entire block is absent from the JSON. A `null` block or `"valid": false` means the fitting attempt produced no usable output.

The estimator partitions reads into quantile-spaced length bins and fits a joint binomial model across damage and control channels simultaneously. For each length bin *l* and terminal position *p*, the expected damage-channel rate is

$$\mu_{l,p} = \eta_{l,p} + (1 - \eta_{l,p}) \cdot \delta_l \cdot e^{-\lambda p}$$

where $\eta_{l,p} = \sigma(\beta_l + \kappa_e \cdot a_{e,p})$ is the composition baseline plus a per-end artifact curve $a_{e,p}$ (with $a_{e,14} = 0$ as the interior anchor), and the control channel provides the same $\eta$ with no damage component. The per-bin damage amplitudes $\delta_l$ are constrained to be isotonically non-increasing in read length (shorter reads carry more terminal damage signal than longer reads in genuinely fragmented ancient material), with the longest-bin amplitude anchored at zero.

### Short-fragment fold-in (opt-in; `short_fragment_floor < 30`)

By default the estimator fits reads ≥ 30 bp only. When short-fragment mode is enabled (the driver lowers `short_fragment_floor` below 30, e.g. to 16), reads in the `[floor, 30)` band are folded into the **same** joint length-isotonic fit rather than estimated separately — so `headline_delta`, `tau_discriminator`, `burial_fingerprint` and `preservation` then reflect **all** lengths ≥ floor, not purely the ≥ 30 band. The fold-in is **not additive**: you cannot difference the short-mode and default headlines. These fields are emitted **only** in short mode, so the default (≥ 30) JSON stays byte-identical for downstream consumers.

- **`short_fragment_joint_refit`** (`true` when present) — flags that the short band was injected into the joint fit and the headline quantities above are the refit values.
- **`short_fragment_delta_metric`** (`"cooccurrence_composition_controlled"`) — names the estimand of the short-band headline delta: the profile's native per-C-site terminal C→T excess, de-confounded by the `(S5, S3)` eligibility conditioning. It is **not** the raw shuffle-null co-occurrence lift reported on the reference-free thread (same order of magnitude, not directly comparable — do not difference the two).
- **`short_fragment_bins`** — one entry per short length bin (anchor 15, width 5, 3 bins → `[15,20)`, `[20,25)`, `[25,30)`; the active subset depends on `floor`). Each carries `length_lo`, `length_hi`, `n_reads`, `mean_len`, and `per_pos_5prime_ct`: the 15-position 5′ terminal $T/(T+C)$ profile for that bin, so the recovered short-band deamination is inspectable per bin.

Two internal consequences of the fold-in (not separate JSON fields): the injected short bin's π interval comes from its **own** per-bin profile-likelihood δ curve, $\{\delta : \ell(\delta) - \hat\ell \ge -\tfrac12\chi^2_{1,0.95}\}$ scaled by the point $d_{\max}$ (replacing the shared-$d_{\max}$ Wald band, which under-covers on thin short-band data); and the GC-mixture `pi_damaged` / damaged-fraction population absorbs the short band as one un-GC-stratified group (short reads have no clean interior to stratify GC over). At `floor = 30` (no injection) both paths are byte-identical to the default.


---