---
type: libtaph-json-field
title: scission
tier: standard
estimand: strand scission rate from exponential fit to fragment-length right tail
stability: stable
emitted_by: profile_json.cpp
---

### `scission`

The `scission` block quantifies the rate at which DNA fragments terminate per unit length, estimated as a single exponential decay parameter (γ, bp⁻¹) fit to the right tail of the fragment-length distribution. It is always present in the profile JSON; when the quality gates are not met, `fitted` is `false` and `gamma`, `ci_lo`, and `ci_hi` are `null`.

#### Method

Fragment lengths are accumulated into a fine-resolution histogram across all reads in the library. The modal length L_mode is taken as the mean length of the most-populated bin. The right tail — the modal bin plus all longer bins — is modelled as a left-truncated exponential with rate γ:

P(L | L ≥ L_mode) = γ · exp(−γ · (L − L_mode))

The sufficient statistic for the MLE is S = Σ nᵢ · (x̄ᵢ − L_mode) over tail bins, giving γ̂ = n_tail / S with standard error γ̂ / √n_tail (Fisher information of the left-truncated exponential). A 95% Wald confidence interval is computed as (max(0, γ̂ − 1.96 SE), γ̂ + 1.96 SE). The fit is attempted only when there are at least 3 non-empty bins in the tail, at least 50 tail reads, and S > 0.

| Field | Type | Description |
|-------|------|-------------|
| `fitted` | `bool` | `true` when quality gates passed and γ was estimated; `false` otherwise |
| `gamma` | `number \| null` | MLE scission rate in bp⁻¹; `null` when not fitted |
| `ci_lo` | `number \| null` | Lower bound of 95% Wald CI for γ; floored at 0; `null` when not fitted |
| `ci_hi` | `number \| null` | Upper bound of 95% Wald CI for γ; `null` when not fitted |
| `mean_length` | `number \| null` | Mean read length across all reads (bp); `null` when no reads |
| `modal_length` | `number \| null` | Mean length within the modal (most-frequent) length bin, used as L_mode (bp); `null` when not fitted |
| `n_tail` | `integer` | Number of reads in the right tail (modal bin and longer) used in the MLE |
| `n_total` | `integer` | Total number of reads in the fine length histogram |

#### Interpretation

Larger γ indicates steeper exponential fall-off to the right of the mode, corresponding to shorter, more fragmented DNA. In practice, γ ≈ 0.02–0.10 bp⁻¹ is typical for moderately degraded ancient material; values above ~0.15 bp⁻¹ reflect severe size truncation, often from aggressive size selection or extreme degradation. When `fitted` is `false`, the library either contains too few reads, the length distribution is too narrow to populate the tail, or is dominated by very short reads that collapse the modal bin. The `mean_length` and `modal_length` fields are always informative even when the fit fails, and should be inspected alongside γ when comparing libraries. Because this estimate is derived entirely from read lengths in the FASTQ, it reflects both biological fragmentation and any upstream size selection or trimming; it cannot be interpreted as an absolute strand-break rate without those processing details.

---