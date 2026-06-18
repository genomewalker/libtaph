---
type: libtaph-json-field
title: deamination
tier: standard
estimand: deamination damage profile across positions and channels
stability: stable
emitted_by: profile_json.cpp
---

### `deamination` block

Cytosine spontaneously deaminates to uracil by hydrolysis of the C4 amino group.
In an intact duplex, uracil-DNA glycosylase removes the lesion efficiently. In the
single-stranded overhangs that form at the ends of nicked or fragmented ancient molecules,
repair is unavailable: the uracil persists and templates thymine during PCR amplification,
producing the characteristic C→T excess in ancient DNA reads [[Briggs et al.
2007](#references); [Lindahl 1993](#references)]. The rate of spontaneous deamination in double-stranded DNA is approximately 7×10⁻¹³ s⁻¹
per cytosine under physiological conditions (~100–500 events per cell per day), and is
accelerated by low pH and elevated temperature [[Lindahl
1993](#references)]; cold, alkaline burial
slows it by orders of magnitude, explaining DNA survival over millennia.

Because only the single-stranded overhang is susceptible, the excess decays exponentially
inward from the fragment end [[Jónsson et al. 2013](#references)]:

$$\delta(p) = b + A \cdot e^{-\lambda p}$$

where $p$ is distance from the terminus (0-indexed), $b$ is the interior background rate,
$A$ is the terminal excess above background, and $\lambda$ is the decay constant. `d_max`
is defined as $A/(1-b)$ — the fraction of terminal cytosine sites that were converted to
thymine, after accounting for the background T proportion.

#### Damage estimates

| Field | Description |
|-------|-------------|
| `d_max_5prime` | C→T excess at the 5′ terminus above the interior background |
| `d_max_3prime` | G→A excess at the 3′ terminus (DS) or C→T excess (SS original-strand) above background |
| `d_max_combined` | Best single damage estimate for this library; `source` records how it was derived |
| `terminal_ct_mixture_amp` | **Relabel of `d_max_combined` to its true estimand.** Channel A's `d_max = A/(1−b)` divides out composition, so it consistently estimates the *product* π_dmg·A_b (damaged-molecule fraction × per-ancient terminal C→T amplitude), **not** per-ancient A_b — which is unidentifiable reference-free. Numerically byte-identical to `d_max_combined`; the `terminal_ct_estimand` object states what it measures. It is an attenuated **lower bound** on the per-ancient amplitude, **not comparable** to mapped metaDMG D_max. |
| `terminal_ct_mixture_amp_valid_as_point` | `false` when `source ∈ {max_ss_asymmetry, min_asymmetry}` (the scalar is order-statistic biased); use `d_max_5prime`/`d_max_3prime` as the honest per-end objects in that case |
| `terminal_ct_estimand` | Static metadata object: `estimand=pi_dmg*A_b_true`, `is_lower_bound_on=A_b_true`, `not_comparable_to=mapped_metaDMG_Dmax`, plus the assumptions under which it holds |
| `w_ancient` | EM ancient-component weight (= `mixture_pi_ancient`) or `null` if EM did not run. Documents the attenuation only — **never** divided into the amplitude to form a per-ancient point estimate |
| `w_ancient_gate_status` | `identified` / `undetermined` / `unavailable` from the mixture-EM status |
| `per_ancient_A_b_lower` | **Lower bound only** on the per-ancient terminal C→T amplitude (= `terminal_ct_mixture_amp`). Since `A_b_true = amp / π_dmg` with `π_dmg ∈ (0,1]`, the true amplitude lies in `[amp, ∞)`: `amp` is a rigorous lower bound. |
| `per_ancient_A_b_upper` | Always `null` — the upper bound is **unidentified reference-free**. A finite ceiling (e.g. `2×amp`) would silently assert a `π_dmg ≥ 0.5` prior, breaking study-independence (at `π_dmg ≈ 0.2`, true `A_b ≈ 5×amp`). No π_dmg assumption is made. |
| `per_ancient_A_b_note` | States the rigorous fact: `upper bound unidentified reference-free (= amp/pi_dmg, pi_dmg unknown)` |
| `d_metamatch` | GC-composition-weighted estimate comparable to metaDMG's D_max [[Michelsen et al. 2022](#references)] |
| `source` | How `d_max_combined` was derived — see table below |
| `lambda_5prime` | Exponential decay rate at the 5′ end; larger values mean shorter overhang length |
| `lambda_3prime` | Exponential decay rate at the 3′ end |
| `bg_5prime` | Interior T/(T+C) baseline at the 5′ end; reflects sequence composition, not damage |
| `bg_3prime` | Interior A/(A+G) baseline at the 3′ end |

**`source` values:**

| Value | Meaning |
|-------|---------|
| `"5prime_only"` | 3′ signal inverted by adapter artefact; `d_max_combined = d_max_5prime` |
| `"3prime_only"` | 5′ signal inverted; `d_max_combined = d_max_3prime` |
| `"average"` | Both ends valid; weighted average or mixture model |
| `"max_ss_asymmetry"` | SS library with unequal end amplitudes; `max(d_max_5, d_max_3)` |
| `"min_asymmetry"` | DS with high end asymmetry; conservative `min(d_max_5, d_max_3)` |
| `"channel_b_structural"` | Channel B stop-codon estimate; used when more reliable than Channel A |
| `"none"` | Both ends inverted; no recoverable signal; `d_max_combined = 0` |

Adapter remnants of 1–2 bp can depress T/(T+C) at position 0, hiding the biological
signal. When this is detected, the fit window is shifted automatically and `d_max_combined`
already reflects the correction. Do not read `per_pos_5prime_ct[0]` directly for such
libraries.

#### Validation flags

Channel B tests the same C→T deamination through stop-codon context
(CAA/CAG/CGA → TAA/TAG/TGA), which is insensitive to GC-composition bias. It provides a
composition-independent cross-check on Channel A.

| Field | Description |
|-------|-------------|
| `validated` | Channel B corroborates Channel A deamination |
| `artifact` | Channel A signal present but Channel B contradicts it; likely GC-composition artefact rather than genuine damage |

#### Per-position rates

| Field | Description |
|-------|-------------|
| `per_pos_5prime_ct[15]` | Raw T/(T+C) at 5′ positions 0–14 |
| `per_pos_3prime[15]` | A/(A+G) at 3′ positions 0–14 for DS; T/(T+C) for SS original-strand |

Values of `-1` indicate insufficient coverage at that position.

#### `cpg_like` sub-block — methylation signal

When cytosine is methylated at the C5 position (5-methylcytosine, 5mC), its deamination
product is thymine rather than uracil [[Shen et al. 1994](#references)]. Deamination of 5mC in double-stranded DNA occurs at approximately twice the rate of
unmodified cytosine [[Shen et al. 1994](#references)]. An additional factor is that the
thymine product creates a T:G mismatch, which is repaired by thymine-DNA glycosylase (TDG)
less efficiently than the U:G mismatch from cytosine deamination is repaired by
uracil-DNA glycosylase [[Wiebauer & Jiricny 1990](#references); [Neddermann & Jiricny
1994](#references)]. In ancient DNA, this produces an
enhanced C→T rate specifically at positions where the downstream base is G. A `cpg_ratio`
above 1 indicates methylation-dependent deamination, the expected pattern for vertebrates
and most eukaryotes with CpG methylation. Values near 1
indicate either unmethylated CpGs (e.g., plants with low global methylation, many fungi),
modern contamination, or signal too eroded to detect context dependence.

| Field | Description |
|-------|-------------|
| `dmax_ct5_cpg` | 5′ C→T amplitude in CpG context |
| `dmax_ct5_noncpg` | 5′ C→T amplitude in non-CpG context |
| `cpg_ratio` | `dmax_cpg / dmax_noncpg`; values >1 indicate methylation-enhanced deamination |
| `log2_cpg_ratio` | log₂(cpg_ratio); 0 means no context dependence |
| `baseline_cpg` | Interior T/(C+T) at CpG positions |
| `baseline_noncpg` | Interior T/(C+T) at non-CpG positions |
| `cov_terminal_cpg` | Total T+C observations at CpG terminal positions |
| `cov_terminal_noncpg` | Total T+C at non-CpG terminal positions |
| `effcov_terminal_cpg` | Effective coverage: `cov × (1 − baseline)` |
| `effcov_terminal_noncpg` | Effective coverage: `cov × (1 − baseline)` |

All fields are `null` when signal is insufficient or the amplitude estimate hits its
boundary.

---