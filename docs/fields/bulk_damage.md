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

where $\eta_{l,p} = \sigma(\beta_l + \kappa_e \cdot a_{e,p})$ is the composition baseline plus a per-end artifact curve $a_{e,p}$ (with $a_{e,14} = 0$ as the interior anchor), and the control channel provides the same $\eta$ with no damage component. The per-bin damage amplitudes $\delta_l$ are constrained to be isotonically non-increasing in read length (shorter reads carry more terminal damage signal than longer reads in genuinely fragmented ancient material), with the longest-bin amplitude anchored at zer


---