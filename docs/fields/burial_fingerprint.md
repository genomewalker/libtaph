---
type: libtaph-json-field
title: burial_fingerprint
tier: standard
estimand: horizon-specific terminal C→T context fingerprint for stratigraphic identification
stability: stable
emitted_by: profile_json.cpp
---

### `burial_fingerprint` block — hydrolytic fragmentation pressure and overhang architecture

This block encodes a joint summary of the two dominant hydrolytic processes shaping ancient DNA fragment architecture: AP-site–driven strand scission (depurination followed by β-elimination) and cytosine deamination. Together these processes leave a coupled imprint on read-terminal damage rates and fragment-length distributions that can be decomposed into a scission-to-deamination pressure ratio and an overhang geometry estimate without reference to any genome. The block is present in every profile but fields are set to sentinel value `-1.0` when the joint model did not converge (`fitted: false`); this occurs when read count, damage amplitude, or fragment-length signal is insufficient to separate the two processes.

Strand scission at apurinic sites and cytosine deamination share a common environmental driver — hydrolysis — but with different rate sensitivities: scission is accelerated by acidity and responds rapidly to pH, while deamination is more sensitive to temperature and sustained water activity. The ratio of the two rates therefore encodes information about the integrated chemical environment experienced by the DNA during burial, independent of any external calibration or reference sequence. The overhang geometry term captures how much of the terminal damage signal is concentrated in single-stranded overhangs versus a pervasive, position-independent floor, which reflects the balance between initial nick density and subsequent exonucleolytic processing or molecular degradation during burial.

| Field | Type | Description |
|-------|------|-------------|
| `fitted` | bool | `true` when the joint scission/deamination model converged; `false` when read count or signal amplitude was insufficient |
| `theta` | float | Hydrolytic fragmentation pressure index: `ln(γ/f₀)`, the log-ratio of AP-site scission rate `γ` to pervasive deamination floor `f₀`. Positive values indicate scission-dominant chemistry; negative values indicate deamination-dominant chemistry. Set to `-1.0` when not fitted. |
| `theta_ci_lo` | float | Lower bound of the 95 % confidence interval on `theta`. Set to `-1.0` when not fitted. |
| `theta_ci_hi` | float | Upper bound of the 95 % confidence interval on `theta`. Set to `-1.0` when not fitted. |
| `tau_hat` | float | Estimated overhang e-fold length in bp: the exponential decay constant of the terminal damage signal. Larger values correspond to longer single-stranded overhangs and less exonucleolytic trimming. Reported even when `fitted` is `false` as a preliminary length estimate. |
| `overhang_fraction` | float | Boundary-robust overhang gate: `A/(A + f₀)`, the fraction of the terminal signal attributable to the single-stranded overhang amplitude `A` relative to the pervasive floor `f₀`. Values near 1 indicate an overhang-dominated signal; values near 0 indicate that terminal damage is spread uniformly across all positions. |
| `phi_share` | float | Upper bound on the oxidative contribution: `σ₀/(σ₀ + f₀)`, where `σ₀` is the oxidation-like component extracted during decomposition. This is a ceiling, not a direct oxidation rate — it constrains how much of the pervasive floor could plausibly be explained by oxidation rather than deamination. Set to `-1.0` when not separately estimable. |
| `note` | string | Embedded interpretation reminder documenting the semantics of `theta` and `phi_share`; present in every output for inline reference. |

`theta` is a chemical-environment proxy, not a direct pH readout: it integrates the effects of pH, temperature, water activity, and deamination saturation over the entire burial history. A high `theta` is consistent with acidic, wet, or warm conditions that accelerate AP-site formation relative to deamination, but the same value can arise from different combinations of these factors. `tau_hat` should be interpreted alongside `deamination.lambda_5prime`: when both are large the overhang is long and well-preserved; when `tau_hat` is large but `lambda_5prime` is small, the fragment length distribution and the damage decay rate are decoupled, which may indicate post-depositional fragmentation without overhang loss. `phi_share` values well below 1 do not exclude oxidative damage — they only constrain the pervasive floor component; the dedicated `oxidation` and `oxog_estimate` blocks should be consulted for oxidation assessment. When `fitted` is `false`, only `tau_hat` and `overhang_fraction` carry partial information; `theta`, `theta_ci_lo`, `theta_ci_hi`, and `phi_share` are set to `-1.0` and must not be interpreted.


---