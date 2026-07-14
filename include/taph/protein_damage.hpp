#pragma once

// Protein-level ancient-damage model: Beta-Binomial LRT with self-calibrated
// null (H0: theta = pi0, the library's own technical floor), Laplace Bayes
// factor, and profile-likelihood CI. libtaph is the single damage engine;
// DART/fqdup are thin consumers that build per-read observations and call this.
// Migrated verbatim from dart/damage_stats.hpp (fit_protein_damage cluster).

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace taph {

inline double log_gamma(double x) {
    // Lanczos coefficients for g=7
    static constexpr double LANCZOS_G = 7.0;
    static constexpr double LANCZOS_COEFF[9] = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };

    if (x < 0.5) {
        // Reflection formula
        return std::log(M_PI / std::sin(M_PI * x)) - log_gamma(1.0 - x);
    }

    x -= 1.0;
    double a = LANCZOS_COEFF[0];
    for (int i = 1; i < 9; ++i) {
        a += LANCZOS_COEFF[i] / (x + i);
    }

    double t = x + LANCZOS_G + 0.5;
    return 0.5 * std::log(2 * M_PI) + (x + 0.5) * std::log(t) - t + std::log(a);
}

inline double log_beta(double a, double b) {
    return log_gamma(a) + log_gamma(b) - log_gamma(a + b);
}

inline double digamma(double x) {
    double result = 0.0;

    // Use recurrence for small x
    while (x < 6.0) {
        result -= 1.0 / x;
        x += 1.0;
    }

    // Asymptotic series for large x
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;
    result += std::log(x) - 0.5 * inv_x
              - inv_x2 * (1.0/12.0 - inv_x2 * (1.0/120.0 - inv_x2 / 252.0));

    return result;
}

inline double trigamma(double x) {
    double result = 0.0;

    // Use recurrence for small x
    while (x < 6.0) {
        result += 1.0 / (x * x);
        x += 1.0;
    }

    // Asymptotic series for large x
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;
    result += inv_x + 0.5 * inv_x2
              + inv_x2 * inv_x * (1.0/6.0 - inv_x2 * (1.0/30.0 - inv_x2 / 42.0));

    return result;
}

struct ProteinDamageResult {
    // MLE estimates
    double delta_max = 0.0;         // Damage rate at terminal (MLE)
    double lambda = 0.3;            // Decay constant (fixed or estimated)
    double phi = 50.0;              // Concentration parameter (controls overdispersion)

    // Log-likelihoods
    double log_lik_m1 = 0.0;        // Log-likelihood under damage model (H1: δ > 0)
    double log_lik_m0 = 0.0;        // Log-likelihood under null (H0: δ = 0)

    // Test statistics
    double likelihood_ratio = 0.0;  // 2 * (log_lik_m1 - log_lik_m0)
    double p_value = 1.0;           // From χ²(1) under null

    // Bayesian posterior
    double log_bayes_factor = 0.0;  // log BF₁₀ (positive favors damage)
    double p_damaged = 0.0;         // P(damaged | data)

    // Confidence intervals (from profile likelihood)
    double ci_lower = 0.0;          // 95% CI lower bound for delta_max
    double ci_upper = 0.0;          // 95% CI upper bound for delta_max
    double se_delta = 0.0;          // Standard error of delta_max

    // Summary statistics from data
    size_t n_reads = 0;             // Number of reads
    size_t n_damaged = 0;           // Reads with detected damage
    size_t total_sites = 0;         // Total damage sites
    double mean_p_damaged = 0.0;    // Mean per-read damage posterior

    // Clustering diagnostics (Binomial fit). A pooled (K,N) tail is anti-conservative
    // when N is concentrated in a few correlated molecules; these bound that failure
    // so a downstream sensitivity filter can flag/collapse without a per-protein phi.
    double max_ni_frac = 0.0;       // max single-read n_i / n_eff (1.0 = all N from one read)
    size_t n_informative_reads = 0; // reads contributing n_eligible > 0

    bool is_significant(double alpha = 0.05) const {
        return p_value < alpha && delta_max > 0.01;
    }
};

// Observation from a single read (for protein aggregation)
struct ProteinReadObs {
    float p_damaged;                // Per-read posterior P(damaged|read)
    float log_lr;                   // Log-likelihood ratio for this read
    float info;                     // Informativeness (sum of pC*D at terminals)
    int n_sites;                    // Number of damage sites detected
    bool is_damaged;                // Hard classification threshold
    int n_eligible;                 // #2: C→T/G→A eligible opportunities on this read (per-opportunity denominator)
    int k_damaged;                  // #2: damage-consistent mismatches on this read (per-opportunity numerator)

    ProteinReadObs() : p_damaged(0), log_lr(0), info(0), n_sites(0), is_damaged(false),
                       n_eligible(0), k_damaged(0) {}
    ProteinReadObs(float p, float lr, float i, int s, bool d)
        : p_damaged(p), log_lr(lr), info(i), n_sites(s), is_damaged(d),
          n_eligible(0), k_damaged(0) {}
};

inline double log_beta_binomial(double k, double n, double theta, double phi) {
    if (n <= 0.0) return 0.0;
    k = std::clamp(k, 0.0, n);

    // Reparameterize: α = θφ, β = (1-θ)φ
    theta = std::clamp(theta, 1e-10, 1.0 - 1e-10);
    phi = std::max(phi, 1.0);

    double alpha = theta * phi;
    double beta = (1.0 - theta) * phi;

    // Beta-Binomial PMF:
    // P(k|n,α,β) = C(n,k) * B(k+α, n-k+β) / B(α, β)
    // log P = log C(n,k) + log B(k+α, n-k+β) - log B(α, β)

    double log_binom = log_gamma(n + 1.0) - log_gamma(k + 1.0) - log_gamma(n - k + 1.0);
    double log_beta_num = log_beta(k + alpha, n - k + beta);
    double log_beta_denom = log_beta(alpha, beta);

    return log_binom + log_beta_num - log_beta_denom;
}

inline ProteinDamageResult fit_protein_damage(
    const std::vector<ProteinReadObs>& obs,
    double prior_theta = 0.5,
    double prior_strength = 2.0,
    double null_theta = 0.05)
{
    ProteinDamageResult result;
    result.n_reads = obs.size();

    if (obs.empty()) {
        // No reads carry damage evidence -> not a damage call (no data). Returning a
        // 0.5 "prior" here produced a spurious spike at p_protein_damaged==0.5.
        result.p_damaged = 0.0;
        return result;
    }

    // Aggregate statistics
    double sum_p = 0.0, sum_p2 = 0.0;
    size_t n_damaged = 0, total_sites = 0;

    for (const auto& o : obs) {
        sum_p += o.p_damaged;
        sum_p2 += o.p_damaged * o.p_damaged;
        if (o.is_damaged) ++n_damaged;
        total_sites += o.n_sites;
    }

    result.n_damaged = n_damaged;
    result.total_sites = total_sites;
    result.mean_p_damaged = sum_p / obs.size();

    const size_t n = obs.size();
    const double n_eff = static_cast<double>(n);
    const double k_eff = sum_p;  // Soft count from per-read posteriors
    const double prior_theta_clamped = std::clamp(prior_theta, 1e-6, 1.0 - 1e-6);
    const double prior_strength_pos = std::max(prior_strength, 1e-6);
    const double prior_alpha = prior_theta_clamped * prior_strength_pos;
    const double prior_beta = (1.0 - prior_theta_clamped) * prior_strength_pos;

    // =========================================================================
    // Method of moments for initial estimates
    // =========================================================================

    // Sample mean and variance of per-read posteriors
    double mean_p = sum_p / n;
    double var_p = (sum_p2 / n) - (mean_p * mean_p);
    var_p = std::max(var_p, 1e-6);  // Avoid division by zero

    // Method of moments for φ:
    // Var(X) = nθ(1-θ)(n+φ)/(1+φ) for Beta-Binomial
    // Solving for φ: φ = (nθ(1-θ) - Var(X)) / (Var(X) - θ(1-θ))
    // But we're using continuous posteriors, so approximate with:
    // φ ≈ θ(1-θ)/var_p - 1

    double theta_init = (k_eff + prior_alpha) / (n_eff + prior_alpha + prior_beta);
    theta_init = std::clamp(theta_init, 0.01, 0.99);
    double phi_init = theta_init * (1.0 - theta_init) / var_p - 1.0;
    phi_init = std::clamp(phi_init, 2.0, 1000.0);  // Reasonable range

    // =========================================================================
    // Null model: θ = π₀, the per-sample TECHNICAL FLOOR (library's own non-damage
    // substitution / p_damaged baseline), NOT θ = 0. Testing against θ≈0 made ANY
    // nonzero soft count k_eff astronomically unlikely under H0 -> log BF saturated
    // -> p_damaged pinned at 1.0 (the bimodal {0.5,1.0} spike). Against θ=π₀ the LRT
    // asks "is this protein's damage rate ABOVE the sample floor" -> graded posterior.
    // The test is one-sided: proteins whose MLE rate sits at/below the floor are
    // forced to p_damaged=0 by the directional guard after the BF (see below).
    // =========================================================================
    const double null_theta_eff = std::clamp(null_theta, 1e-4, 0.5);
    result.log_lik_m0 = log_beta_binomial(k_eff, n_eff, null_theta_eff, phi_init);

    // =========================================================================
    // Damage model: θ > 0 with Newton-Raphson MLE
    // =========================================================================

    double theta = theta_init;
    double phi = phi_init;

    // Newton iterations for theta while holding phi fixed
    for (int iter = 0; iter < 20; ++iter) {
        double alpha = theta * phi;
        double beta = (1.0 - theta) * phi;

        // Score: d/dθ log L
        // = φ * (digamma(k + α) - digamma(n - k + β) - digamma(α) + digamma(β))
        double score = phi * (digamma(k_eff + alpha) - digamma(n_eff - k_eff + beta)
                              - digamma(alpha) + digamma(beta));

        // Fisher information (expected)
        double info = phi * phi * (trigamma(k_eff + alpha) + trigamma(n_eff - k_eff + beta)
                                   - trigamma(alpha) - trigamma(beta));
        info = std::max(std::abs(info), 1e-6);

        double step = score / info;
        theta = std::clamp(theta + step, 1e-4, 1.0 - 1e-4);

        if (std::abs(step) < 1e-6) break;
    }

    result.delta_max = theta;
    result.phi = phi;
    result.log_lik_m1 = log_beta_binomial(k_eff, n_eff, theta, phi);

    // =========================================================================
    // Likelihood Ratio Test
    // =========================================================================
    result.likelihood_ratio = 2.0 * (result.log_lik_m1 - result.log_lik_m0);
    result.likelihood_ratio = std::max(result.likelihood_ratio, 0.0);

    // P-value from χ²(1) distribution using survival function approximation
    // For LRT with one constrained parameter
    if (result.likelihood_ratio > 0.0) {
        // Quick approximation: P(χ² > x) ≈ erfc(sqrt(x/2))
        double x = result.likelihood_ratio;
        result.p_value = std::erfc(std::sqrt(x / 2.0));
    } else {
        result.p_value = 1.0;
    }

    // =========================================================================
    // Bayes Factor via Laplace approximation
    // =========================================================================
    // BF₁₀ ≈ exp(log_lik_m1 - log_lik_m0) * sqrt(2π/I(θ̂)) * p(θ̂) / ∫p(θ)dθ
    // With uniform prior on [0,1], this simplifies
    double alpha_hat = theta * phi;
    double beta_hat = (1.0 - theta) * phi;
    double fisher_info = phi * phi * (trigamma(k_eff + alpha_hat) + trigamma(n_eff - k_eff + beta_hat)
                                      - trigamma(alpha_hat) - trigamma(beta_hat));
    fisher_info = std::max(std::abs(fisher_info), 1e-6);

    result.se_delta = 1.0 / std::sqrt(fisher_info);

    // Laplace approximation to Bayes factor
    double log_bf = result.log_lik_m1 - result.log_lik_m0;
    log_bf += 0.5 * std::log(2 * M_PI / fisher_info);  // Laplace correction
    result.log_bayes_factor = log_bf;

    // Posterior probability with uniform prior
    // P(H1|data) = BF10 / (1 + BF10)
    double bf = std::exp(std::clamp(log_bf, -100.0, 100.0));
    result.p_damaged = bf / (1.0 + bf);

    // Directional: the damage test is one-sided (H1: θ > π₀). If the MLE damage rate
    // sits at or below the sample floor there is no evidence for damage ABOVE the
    // floor, whatever the Laplace/BF arithmetic says -> not a damage call.
    if (result.delta_max <= null_theta_eff) {
        result.p_damaged = 0.0;
        result.likelihood_ratio = 0.0;
        result.p_value = 1.0;
    }

    // =========================================================================
    // Profile Likelihood Confidence Interval
    // =========================================================================
    // 95% CI: {θ : 2(log L(θ̂) - log L(θ)) ≤ χ²₀.₉₅(1) = 3.84}

    double chi2_crit = 3.84;  // χ²(1) at 95%
    double ll_max = result.log_lik_m1;
    double ll_threshold = ll_max - chi2_crit / 2.0;

    // Search for lower bound
    double lower = 0.0;
    for (double t = theta; t >= 1e-4; t -= 0.01) {
        double ll = log_beta_binomial(k_eff, n_eff, t, phi);
        if (ll < ll_threshold) {
            lower = t + 0.01;
            break;
        }
        lower = t;
    }

    // Search for upper bound
    double upper = 1.0;
    for (double t = theta; t <= 1.0 - 1e-4; t += 0.01) {
        double ll = log_beta_binomial(k_eff, n_eff, t, phi);
        if (ll < ll_threshold) {
            upper = t - 0.01;
            break;
        }
        upper = t;
    }

    result.ci_lower = std::max(lower, 0.0);
    result.ci_upper = std::min(upper, 1.0);

    return result;
}

// Binomial log-likelihood of k successes in n Bernoulli(theta) trials, up to the
// n-choose-k constant (which is independent of theta and so cancels in every LR /
// Bayes-factor difference below). Kept separate from log_beta_binomial: the
// per-opportunity model here is a pure Binomial, not overdispersed (see note).
inline double log_binomial_no_const(double k, double n, double theta) {
    theta = std::clamp(theta, 1e-12, 1.0 - 1e-12);
    return k * std::log(theta) + (n - k) * std::log(1.0 - theta);
}

// ── #2: per-OPPORTUNITY one-sided Binomial LRT (post-migration deliberate correction) ─
// fit_protein_damage above aggregates SOFT per-read posteriors (k_eff=Σp_damaged,
// n_eff=n_reads): every eligible-clean read (damage opportunities present, none damaged)
// is credited a fabricated prior-derived soft count, diluting/biasing the protein rate.
// This variant uses the calibration-consistent object (the M1 per-opportunity grain):
//   k = damage-consistent mismatches, n = C→T/G→A eligible opportunities.
// Info-free reads (n_eligible==0) contribute nothing; eligible-clean reads (k=0, n>0)
// contribute legitimate 0/n negative evidence. null_theta is the per-OPPORTUNITY
// technical floor (q_modern: baseline substitution rate at eligible sites), NOT the
// per-read p_damaged floor.
//
// Model choice — why Binomial, not Beta-Binomial:
//   The per-read opportunity counts are dominated by single-opportunity reads (n_i=1),
//   so a per-protein overdispersion φ is NOT identifiable: estimating φ from the moments
//   of p_i=k_i/n_i just reads back Bernoulli variance as overdispersion and, held fixed
//   at the moment estimate, lets H0 "explain" a genuine excess (a prior review showed a
//   single hit at rate == null scoring p≈0.66). We therefore model the opportunities as
//   exchangeable Bernoulli(theta) within a protein — for which the pooled (k_eff, n_eff)
//   are the sufficient statistic, so aggregation across the protein's reads is exact —
//   and test one-sided H1: theta > null. Between-read heterogeneity, if it needs to be
//   accounted for, belongs at the DATASET level (a shared/calibrated dispersion or
//   variance-inflation on the LR), not re-estimated per protein; that hook is deliberately
//   left out here pending calibration. fit_protein_damage is left byte-for-byte unchanged
//   so the migration bit-exactness proof still targets it directly.
inline ProteinDamageResult fit_protein_damage_binomial(
    const std::vector<ProteinReadObs>& obs,
    double null_theta,
    double prior_theta = 0.5,
    double prior_strength = 2.0)
{
    // prior_theta/prior_strength are retained for call-site/signature compatibility;
    // a Beta pseudo-count would bias an all-clean protein's rate ABOVE null and trip a
    // false call, so the one-sided test uses the unregularized MLE.
    (void)prior_theta;
    (void)prior_strength;

    ProteinDamageResult result;
    result.n_reads = obs.size();
    result.phi = 0.0;  // no overdispersion parameter in the Binomial model

    double k_eff = 0.0, n_eff = 0.0, max_ni = 0.0;
    size_t n_damaged = 0, total_sites = 0, n_informative = 0;
    for (const auto& o : obs) {
        if (o.is_damaged) ++n_damaged;
        total_sites += static_cast<size_t>(o.n_sites);
        if (o.n_eligible <= 0) continue;   // info-free: contributes n=0 (nothing)
        ++n_informative;
        const double ni = static_cast<double>(o.n_eligible);
        if (ni > max_ni) max_ni = ni;
        k_eff += static_cast<double>(o.k_damaged);
        n_eff += ni;
    }
    result.n_damaged = n_damaged;
    result.total_sites = total_sites;
    result.n_informative_reads = n_informative;
    result.max_ni_frac = (n_eff > 0.0) ? max_ni / n_eff : 0.0;

    if (n_informative == 0 || n_eff <= 0.0) {
        // No eligible opportunities anywhere -> no data, no damage call.
        result.p_damaged = 0.0;
        return result;
    }

    const double null_eff = std::clamp(null_theta, 1e-4, 0.5);
    const double theta_hat = k_eff / n_eff;          // unregularized Binomial MLE
    result.mean_p_damaged = theta_hat;
    result.delta_max = theta_hat;

    // One-sided: an observed rate at or below the technical floor is no evidence for
    // damage. All-clean proteins (k_eff=0 -> theta_hat=0) land here by construction.
    if (theta_hat <= null_eff) {
        result.log_lik_m0 = log_binomial_no_const(k_eff, n_eff, null_eff);
        result.log_lik_m1 = result.log_lik_m0;
        result.se_delta = std::sqrt(std::max(theta_hat * (1.0 - theta_hat), 1e-12) / n_eff);
        result.ci_lower = 0.0;
        result.ci_upper = theta_hat + 1.96 * result.se_delta;
        result.p_damaged = 0.0;
        result.likelihood_ratio = 0.0;
        result.p_value = 1.0;
        return result;
    }

    const double ll1 = log_binomial_no_const(k_eff, n_eff, theta_hat);
    const double ll0 = log_binomial_no_const(k_eff, n_eff, null_eff);
    result.log_lik_m1 = ll1;
    result.log_lik_m0 = ll0;

    result.likelihood_ratio = std::max(2.0 * (ll1 - ll0), 0.0);
    // One-sided LRT p-value. For H0: theta = null vs H1: theta > null (null interior, MLE in
    // H1 here since theta_hat>null is guaranteed above), the LR statistic — 0 when
    // theta_hat<=null — has the Self–Liang/Chernoff null distribution 0.5·χ²₀ + 0.5·χ²₁, so
    // the tail is 0.5·P(χ²₁ > LR) = 0.5·erfc(sqrt(LR/2)). The bare erfc is the FULL χ²₁ tail
    // and double-counts the down side (over-conservative + mislabeled).
    result.p_value = (result.likelihood_ratio > 0.0)
        ? 0.5 * std::erfc(std::sqrt(result.likelihood_ratio / 2.0)) : 1.0;

    // PRIMARY reported quantity: LRT-based detection confidence that theta exceeds the null.
    // Kept as p_damaged (the emitted p_protein_damaged) so downstream "higher = more damage"
    // thresholds hold, and — unlike the Laplace BF below — it is well-behaved at the K=N
    // boundary MLE (where the Gaussian-integral approximation degenerates).
    result.p_damaged = 1.0 - result.p_value;

    // Laplace Bayes factor: retained as a SECONDARY diagnostic only. Fisher info for a
    // Binomial rate at the MLE is I = n/(theta(1-theta)); at a boundary MLE (theta_hat→1,
    // i.e. K=N) the variance → 0, the Laplace approximation is invalid and the volume term
    // blows up, so guard it rather than emitting a numerically arbitrary value.
    const double var_hat = theta_hat * (1.0 - theta_hat);
    if (var_hat > 1e-9) {
        const double fisher_info = n_eff / var_hat;
        result.se_delta = 1.0 / std::sqrt(fisher_info);
        result.log_bayes_factor = (ll1 - ll0) + 0.5 * std::log(2.0 * M_PI / fisher_info);
    } else {
        // Boundary MLE (K=N): Laplace volume undefined; report the half-LR as the BF proxy.
        result.se_delta = 0.0;
        result.log_bayes_factor = 0.5 * result.likelihood_ratio;
    }

    // Profile-likelihood 95% CI on the Binomial rate.
    const double chi2_crit = 3.84;
    const double ll_threshold = ll1 - chi2_crit / 2.0;
    double lower = theta_hat;
    for (double t = theta_hat; t >= 1e-4; t -= 0.01) {
        if (log_binomial_no_const(k_eff, n_eff, t) < ll_threshold) { lower = t + 0.01; break; }
        lower = t;
    }
    double upper = theta_hat;
    for (double t = theta_hat; t <= 1.0 - 1e-4; t += 0.01) {
        if (log_binomial_no_const(k_eff, n_eff, t) < ll_threshold) { upper = t - 0.01; break; }
        upper = t;
    }
    result.ci_lower = std::max(lower, 0.0);
    result.ci_upper = std::min(upper, 1.0);

    return result;
}

// Self-calibrated protein null pi0: mean per-read p_damaged in the sample's
// low tail, clamped to [floor, cap]; floor fallback when the tail is empty.
// DART supplies the accumulated low-tail sum/count (data); the H0 floor/cap
// semantics live with the model (engine).
inline float self_calibrated_pi0(double low_tail_sum, std::size_t low_tail_n,
                                 float floor = 0.02f, float cap = 0.5f) {
    float pi0 = (low_tail_n > 0)
        ? static_cast<float>(low_tail_sum / static_cast<double>(low_tail_n))
        : floor;
    return std::clamp(pi0, floor, cap);
}

} // namespace taph
