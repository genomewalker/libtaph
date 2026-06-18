---
type: libtaph-json-field
title: per_read_overdispersion
tier: standard
estimand: overdispersion of damage across reads; tests for heterogeneous ancient mixture
stability: stable
emitted_by: profile_json.cpp
---

### `per_read_overdispersion`

This block characterises the between-read spread of deamination damage beyond what a
binomial model predicts. In a purely modern library every read draws its terminal C→T
substitution independently with the same probability; the resulting read-level damage
score follows a near-binomial distribution and its variance equals its mean (coefficient
of variation squared ≈ 0). When the library is a mixture of genuinely ancient and modern
molecules the damaged sub-population introduces an additional between-class variance
component: ancient reads carry systematically more terminal damage than modern reads, so
the marginal variance exceeds the binomial expectation by a term proportional to
π(1 − π)d²\_max, where π is the ancient fraction and d\_max is the per-read damage
amplitude. The cross-cumulant fields — the empirical second joint cumulants κ₂ of the
5′ C→T rate paired with the 3′ G→A rate, or with the 3′ C→T rate, or with 5′ and 3′
oxidative G rates — provide complementary handles on the same between-read heterogeneity
that are channel-specific and partially independent of read-length variation. The block
is emitted for every run regardless of library type; there is no null/absent state.

| Field | Type | Description |
|-------|------|-------------|
| `n_damaged_reads` | integer | Number of reads that contributed to the overdispersion statistics (reads passing the minimum-quality filter used to compute the deamination score) |
| `mean_deam_score` | float | Mean per-read deamination score across all qualifying reads; the score combines 5′ C→T and 3′ G→A terminal counts normalised by read length |
| `variance` | float | Observed variance of the per-read deamination score; includes both within-read binomial noise and between-read heterogeneity |
| `cv2` | float | Squared coefficient of variation of the per-read deamination score (`variance / mean_deam_score²`); values substantially above 0 indicate between-read heterogeneity consistent with a mixture |
| `ct5_mean` | float | Mean 5′ terminal C→T fraction across reads (first-position C→T rate averaged over qualifying reads) |
| `ga3_mean` | float | Mean 3′ terminal G→A fraction across reads (last-position G→A rate, the complement-strand deamination signal) |
| `k2_ct5_ga3` | float | Empirical second joint cumulant κ₂(CT5, GA3): sample covariance of the 5′ C→T rate and the 3′ G→A rate across reads; the primary cross-cumulant estimator for mixture overdispersion — expected near zero for homogeneous modern libraries and positive when ancient reads co-express both terminal signals |
| `k2_ct5_ct3` | float | Empirical second joint cumulant κ₂(CT5, CT3): covariance of 5′ and 3′ C→T rates; non-zero only in single-stranded or biased libraries where 3′ C→T is present; provides a complementary same-strand cross-cumulant |
| `k2_ct5_ga3_corr` | float | Residual cross-cumulant after correcting for length-composition coupling; removes the fraction of κ₂(CT5, GA3) attributable to read-length variation shared between the two terminal channels |
| `k2_tpg` | float | Second joint cumulant involving the T-at-5′ and G-at-3′ terminal fractions; a control cross-cumulant using non-deamination terminal bases to probe composition covariance not driven by cytosine deamination |
| `n_tpg_reads` | float | Number of reads entering the T-at-5′/G-at-3′ cumulant calculation |
| `g5_mean` | float | Mean 5′ terminal G→T fraction across reads; the 5′ oxidative signal (8-oxoguanine → G→T at the 5′ read terminus) |
| `g3_mean` | float | Mean 3′ terminal G→T fraction across reads; the 3′ oxidative signal |
| `k2_g5g3` | float | Empirical second joint cumulant κ₂(G5, G3): covariance of 5′ and 3′ G→T rates across reads; expected negative under 8-oxoguanine damage because oxidative damage is strand-specific and the two terminals report opposite strands |
| `score_len_cov` | float | Covariance between the per-read deamination score and read length; positive values indicate that longer reads carry more damage on average, consistent with terminal-decay models where additional positions accumulate signal; used to partition length-driven variance from true between-molecule heterogeneity |

**Interpretation.** A `cv2` substantially above 0 (typically > 0.2 in libraries with
high ancient fractions) combined with a positive `k2_ct5_ga3` is the primary signature
of a heterogeneous mixture of ancient and modern molecules. The cross-cumulant
`k2_ct5_ga3` is more specific than `cv2` because it requires co-expression of both
terminal deamination signals within the same read, which a per-base technical artifact
or simple composition noise does not produce. Values of `k2_ct5_ga3_corr` near zero
after correction confirm that the raw `k2_ct5_ga3` signal is genuine rather than a
read-length confound. The oxidative channel `k2_g5g3` typically takes a negative value
because 8-oxoguanine on the sequenced strand appears preferentially at 3′ termini of
short reads when the damage is strand-asymmetric; a near-zero or positive `k2_g5g3`
alongside elevated `g5_mean` or `g3_mean` warrants comparison with the `oxog_estimate`
and `oxidation` blocks. When `n_damaged_reads` is small (< ~50 000), cross-cumulant
estimates are noisy and should be interpreted with caution; the `score_len_cov` field
is useful for diagnosing libraries where apparent overdispersion is length-driven rather
than mixture-driven.


---