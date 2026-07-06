#pragma once

#include <array>
#include <cmath>
#include <algorithm>

namespace taph {

struct JointDamageResult {
    // Fitted parameters
    float delta_max = 0.0f;      // Maximum damage rate at terminal (δ_max)
    float lambda = 0.0f;         // Decay constant
    float a_max = 0.0f;          // Artifact amplitude (can be negative for inverted)

    // Baselines (fixed from interior data)
    float b_tc = 0.0f;           // Interior T/(T+C) baseline
    float b_ag = 0.0f;           // Interior A/(A+G) baseline
    float b_stop = 0.0f;         // Interior stop/(pre+stop) baseline

    // Model comparison
    // double: at 1e8+ binomial trials the LL/BIC magnitudes reach ~1e9 where float
    // epsilon (~64-256) swamps the ~10-200-unit ΔBIC that distinguishes the models.
    double log_lik_m1 = 0.0;     // Log-likelihood for M1 (δ_max > 0)
    double log_lik_m0 = 0.0;     // Log-likelihood for M0 (δ_max = 0)
    double bic_m1 = 0.0;         // BIC for M1
    double bic_m0 = 0.0;         // BIC for M0
    double delta_bic = 0.0;      // BIC_M0 - BIC_M1 (positive favors damage)
    // ΔBIC accumulates over n_trials reads, so it scales ~O(N) for any fixed rate
    // difference: exploratory, not a calibrated Bayes factor on correlated reads.
    float delta_bic_normalized = 0.0f; // ΔBIC / log(n_trials): per-parameter, N-independent
    bool  delta_bic_saturated = false; // |ΔBIC| beyond interpretable evidence scale (|·|>200)
    float bayes_factor = 0.0f;   // BF_10 ≈ exp(ΔBIC/2)
    float p_damage = 0.0f;       // P(damage | data) = BF / (1 + BF)

    // Diagnostics
    float rmse = 0.0f;           // Root mean squared error of fit
    int n_positions = 15;        // Number of positions used
    uint64_t n_trials = 0;       // Total binomial trials (for BIC)
    // Wald z for δ̂ from observed Fisher information. The Fisher info sums over all
    // reads at 15 positions treated as independent Bernoulli trials, so it inflates
    // ~O(N) and z ~O(sqrt(N)·δ) on correlated reads (PCR dups, admixture): NOT a
    // calibrated test statistic. Use delta_max as the N-independent effect size.
    float z_delta = 0.0f;
    bool  z_delta_capped = false; // |z_delta| exceeded Z_CAP (exploratory magnitude clamped)
    float se_delta = 0.0f;       // se(δ̂) from the same curvature
    int n_informative = 0;       // 5' Channel-A positions with cov≥100 & excess > 2·noise
    bool valid = false;          // Sufficient data for estimation
};

struct JointDamageSuffStats {
    // Channel A: T/(T+C) at 5' end
    std::array<uint64_t, 15> k_tc = {};   // T counts (successes)
    std::array<uint64_t, 15> n_tc = {};   // T+C counts (trials)

    // Control: A/(A+G) at 5' end (same artifact, no damage)
    std::array<uint64_t, 15> k_ag = {};   // A counts
    std::array<uint64_t, 15> n_ag = {};   // A+G counts

    // Channel B: stop/(pre+stop) conversions
    std::array<uint64_t, 15> k_stop = {}; // Stop codon counts
    std::array<uint64_t, 15> n_stop = {}; // Pre-stop + stop counts

    // Interior baselines
    uint64_t k_tc_interior = 0, n_tc_interior = 0;
    uint64_t k_ag_interior = 0, n_ag_interior = 0;
    uint64_t k_stop_interior = 0, n_stop_interior = 0;

    // Compute baselines from interior
    float baseline_tc() const {
        return n_tc_interior > 0 ? static_cast<float>(k_tc_interior) / n_tc_interior : 0.5f;
    }
    float baseline_ag() const {
        return n_ag_interior > 0 ? static_cast<float>(k_ag_interior) / n_ag_interior : 0.5f;
    }
    float baseline_stop() const {
        return n_stop_interior > 0 ? static_cast<float>(k_stop_interior) / n_stop_interior : 0.05f;
    }

    uint64_t total_trials() const {
        uint64_t total = 0;
        for (int p = 0; p < 15; ++p) {
            total += n_tc[p] + n_ag[p] + n_stop[p];
        }
        return total;
    }

    bool is_valid() const {
        // Coverage floor over DATA-BEARING positions only. pos-0 artifact masking
        // zeroes n_tc[0] and short-fragment libraries zero distal positions; a min
        // over all 15 would wrongly invalidate exactly the heavily-damaged / short
        // libraries this model targets (n_tc[0]==0 < 100 ⇒ valid=false). Require a
        // minimum count of informative positions instead.
        uint64_t min_tc = ~0ULL, min_ag = ~0ULL;
        int n_pos_tc = 0, n_pos_ag = 0;
        for (int p = 0; p < 15; ++p) {
            if (n_tc[p] > 0) { min_tc = std::min(min_tc, n_tc[p]); ++n_pos_tc; }
            if (n_ag[p] > 0) { min_ag = std::min(min_ag, n_ag[p]); ++n_pos_ag; }
        }
        return n_pos_tc >= 10 && n_pos_ag >= 10 &&
               min_tc >= 100 && min_ag >= 100 && n_tc_interior >= 1000;
    }
};

class JointDamageModel {
public:
    // Grid parameters
    static constexpr int N_LAMBDA = 10;
    static constexpr int N_DELTA = 100;
    static constexpr float LAMBDA_MIN = 0.05f;
    static constexpr float LAMBDA_MAX = 0.50f;
    static constexpr float DELTA_MAX_LIMIT = 0.60f;
    // Reporting clamps for the exploratory, N-inflating diagnostics (z_delta, ΔBIC).
    static constexpr float Z_CAP = 12.0f;          // emitted Wald-z magnitude bound
    static constexpr float DELTA_BIC_CAP = 200.0f; // emitted ΔBIC saturation bound

    static JointDamageResult fit(const JointDamageSuffStats& stats);

    static double log_likelihood(
        const JointDamageSuffStats& stats,
        float delta_max, float lambda, float a_max,
        float b_tc, float b_ag, float b_stop);

    static float optimize_a_max(
        const JointDamageSuffStats& stats,
        float delta_max, float lambda,
        float b_tc, float b_ag);

private:
    // Binomial log-likelihood: k * log(p) + (n-k) * log(1-p)
    static double binom_ll(uint64_t k, uint64_t n, float p) {
        if (n == 0) return 0.0;
        p = std::clamp(p, 1e-10f, 1.0f - 1e-10f);
        const double pd = static_cast<double>(p);
        return static_cast<double>(k) * std::log(pd) +
               static_cast<double>(n - k) * std::log(1.0 - pd);
    }
};

inline double JointDamageModel::log_likelihood(
    const JointDamageSuffStats& stats,
    float delta_max, float lambda, float a_max,
    float b_tc, float b_ag, float b_stop)
{
    double ll = 0.0;

    for (int p = 0; p < 15; ++p) {
        float decay = std::exp(-lambda * p);
        float delta_p = delta_max * decay;
        float a_p = a_max * decay;

        // Channel A: π_TC(p) = b_tc + a(p) + (1 - b_tc - a(p)) · δ(p)
        float base_tc = std::clamp(b_tc + a_p, 0.01f, 0.99f);
        float pi_tc = base_tc + (1.0f - base_tc) * delta_p;
        ll += binom_ll(stats.k_tc[p], stats.n_tc[p], pi_tc);

        // Control: π_AG(p) = b_ag + a(p)
        float pi_ag = std::clamp(b_ag + a_p, 0.01f, 0.99f);
        ll += binom_ll(stats.k_ag[p], stats.n_ag[p], pi_ag);

        // Channel B: π_stop(p) = b_stop + a(p) + (1 − b_stop − a(p)) · δ(p)
        // Same compositional artifact correction as Channel A: AT-rich sequences
        // have elevated T-starting codon rates at termini that are not damage-induced.
        float base_stop = std::clamp(b_stop + a_p, 0.001f, 0.999f);
        float pi_stop   = std::clamp(base_stop + (1.0f - base_stop) * delta_p, 0.001f, 0.999f);
        ll += binom_ll(stats.k_stop[p], stats.n_stop[p], pi_stop);
    }

    return ll;
}

inline float JointDamageModel::optimize_a_max(
    const JointDamageSuffStats& stats,
    float delta_max, float lambda,
    float b_tc, float b_ag)
{
    // Golden section search for a_max in [-0.3, 0.3]
    const float phi = 0.618033988749895f;
    float a = -0.3f, b = 0.3f;
    float c = b - phi * (b - a);
    float d = a + phi * (b - a);

    auto eval = [&](float a_max) {
        double ll = 0.0;
        for (int p = 0; p < 15; ++p) {
            float decay = std::exp(-lambda * p);
            float delta_p = delta_max * decay;
            float a_p = a_max * decay;

            // Channel A
            float base_tc = std::clamp(b_tc + a_p, 0.01f, 0.99f);
            float pi_tc = base_tc + (1.0f - base_tc) * delta_p;
            ll += binom_ll(stats.k_tc[p], stats.n_tc[p], pi_tc);

            // Control
            float pi_ag = std::clamp(b_ag + a_p, 0.01f, 0.99f);
            ll += binom_ll(stats.k_ag[p], stats.n_ag[p], pi_ag);
        }
        return ll;
    };

    for (int iter = 0; iter < 30; ++iter) {
        if (eval(c) > eval(d)) {
            b = d;
            d = c;
            c = b - phi * (b - a);
        } else {
            a = c;
            c = d;
            d = a + phi * (b - a);
        }
    }

    return (a + b) / 2.0f;
}

inline JointDamageResult JointDamageModel::fit(const JointDamageSuffStats& stats) {
    JointDamageResult result;

    // Get baselines from interior
    result.b_tc = stats.baseline_tc();
    result.b_ag = stats.baseline_ag();
    result.b_stop = stats.baseline_stop();
    result.n_trials = stats.total_trials();
    result.valid = stats.is_valid();

    if (!result.valid) {
        return result;
    }

    // Grid search over (λ, δ_max)
    double best_ll = -1e300;
    float best_delta = 0.0f;
    float best_lambda = 0.2f;
    float best_a = 0.0f;

    for (int i_lambda = 0; i_lambda < N_LAMBDA; ++i_lambda) {
        float lambda = LAMBDA_MIN + (LAMBDA_MAX - LAMBDA_MIN) * i_lambda / (N_LAMBDA - 1);

        for (int i_delta = 0; i_delta <= N_DELTA; ++i_delta) {
            float delta_max = DELTA_MAX_LIMIT * i_delta / N_DELTA;

            // Optimize a_max for this (λ, δ_max)
            float a_max = optimize_a_max(stats, delta_max, lambda,
                                         result.b_tc, result.b_ag);

            double ll = log_likelihood(stats, delta_max, lambda, a_max,
                                      result.b_tc, result.b_ag, result.b_stop);

            if (ll > best_ll) {
                best_ll = ll;
                best_delta = delta_max;
                best_lambda = lambda;
                best_a = a_max;
            }
        }
    }

    result.delta_max = best_delta;
    result.lambda = best_lambda;
    result.a_max = best_a;
    result.log_lik_m1 = best_ll;

    // Compute M0 likelihood (δ_max = 0)
    float a_max_m0 = optimize_a_max(stats, 0.0f, best_lambda,
                                    result.b_tc, result.b_ag);
    result.log_lik_m0 = log_likelihood(stats, 0.0f, best_lambda, a_max_m0,
                                       result.b_tc, result.b_ag, result.b_stop);

    // BIC comparison
    // M1 has 3 parameters (δ_max, λ, a_max), M0 has 2 (λ, a_max)
    // BIC = -2 * log_lik + k * log(N)
    double log_n = std::log(static_cast<double>(result.n_trials));
    result.bic_m1 = -2.0 * result.log_lik_m1 + 3.0 * log_n;
    result.bic_m0 = -2.0 * result.log_lik_m0 + 2.0 * log_n;

    // ΔBIC = BIC_M0 - BIC_M1 (positive favors M1 = damage)
    result.delta_bic = result.bic_m0 - result.bic_m1;

    // ΔBIC accumulates over n_trials reads → ~O(N), uninterpretable as a Bayes factor.
    // Report the per-parameter, N-independent contribution and flag saturation; the raw
    // value is kept (exploratory only). p_damage below uses a stable clamped logistic.
    result.delta_bic_normalized = log_n > 0.0f ? result.delta_bic / log_n : 0.0f;
    result.delta_bic_saturated  = std::abs(result.delta_bic) > DELTA_BIC_CAP;

    // Bayes factor approximation: BF_10 ≈ exp(ΔBIC/2). Cap exponent to keep
    // BF / p_damage finite under extreme separations (avoid inf/inf = NaN).
    float half_dbic = static_cast<float>(std::clamp(result.delta_bic / 2.0, -80.0, 80.0));
    result.bayes_factor = std::exp(half_dbic);

    // P(damage) computed via stable logistic on ΔBIC/2:
    //   p = 1 / (1 + exp(-ΔBIC/2))
    // Equivalent to BF/(1+BF) with prior 0.5, but never NaN.
    result.p_damage = 1.0f / (1.0f + std::exp(-half_dbic));

    // Compute RMSE for diagnostics
    float sse = 0.0f;
    int n_obs = 0;
    for (int p = 0; p < 15; ++p) {
        if (stats.n_tc[p] > 0) {
            float decay = std::exp(-result.lambda * p);
            float base_tc = std::clamp(result.b_tc + result.a_max * decay, 0.01f, 0.99f);
            float pi_tc = base_tc + (1.0f - base_tc) * result.delta_max * decay;
            float obs_tc = static_cast<float>(stats.k_tc[p]) / stats.n_tc[p];
            sse += (pi_tc - obs_tc) * (pi_tc - obs_tc);
            ++n_obs;
        }
    }
    result.rmse = n_obs > 0 ? std::sqrt(sse / n_obs) : 0.0f;
    result.n_positions = n_obs;

    // Wald z for δ̂ from the observed Fisher information I(δ̂) = -d²ℓ/dδ², via a
    // boundary-safe central second difference with λ, a_max and baselines held at
    // their fitted optima. δ̂ pinned to the lower boundary (0) or a non-concave
    // neighbourhood (I ≤ 0) yields z = 0 — no Wald signal to report.
    if (result.delta_max > 1e-6f) {
        const float h = DELTA_MAX_LIMIT / N_DELTA;
        float dl = std::clamp(result.delta_max - h, 0.0f, DELTA_MAX_LIMIT);
        float dr = std::clamp(result.delta_max + h, 0.0f, DELTA_MAX_LIMIT);
        float hl = result.delta_max - dl;
        float hr = dr - result.delta_max;
        if (hl > 1e-9f && hr > 1e-9f) {
            double ll_l = log_likelihood(stats, dl, result.lambda, result.a_max,
                                        result.b_tc, result.b_ag, result.b_stop);
            double ll_r = log_likelihood(stats, dr, result.lambda, result.a_max,
                                        result.b_tc, result.b_ag, result.b_stop);
            double fisher = -2.0 * (static_cast<double>(hr) * ll_l
                            - (static_cast<double>(hl) + hr) * result.log_lik_m1
                            + static_cast<double>(hl) * ll_r) /
                           (static_cast<double>(hl) * hr * (hl + hr));
            if (fisher > 0.0) {
                result.se_delta = 1.0f / std::sqrt(fisher);
                // exploratory; clamped, not a calibrated p-value (correlated reads)
                float z_raw = result.delta_max / result.se_delta;
                result.z_delta_capped = std::abs(z_raw) > Z_CAP;
                result.z_delta = std::clamp(z_raw, -Z_CAP, Z_CAP);
            }
        }
    }

    // Informative positions: 5' Channel-A coverage ≥ 100 AND observed terminal
    // excess over the interior baseline exceeds 2× its binomial sampling noise.
    for (int p = 0; p < 15; ++p) {
        if (stats.n_tc[p] < 100) continue;
        float obs = static_cast<float>(stats.k_tc[p]) / stats.n_tc[p];
        float se_pos = std::sqrt(std::max(result.b_tc * (1.0f - result.b_tc), 1e-8f) /
                                 static_cast<float>(stats.n_tc[p]));
        if (obs - result.b_tc > 2.0f * se_pos) ++result.n_informative;
    }

    return result;
}

struct DamageMixtureResult {
    float d_mean = 0.0f;         // GC-weighted average (reference-free population avg)
    float d_ancient = 0.0f;      // Damaged component mean (μ_1)
    float pi_ancient = 0.0f;     // Fraction of C-sites in damaged component
    float tau_ancient = 0.0f;    // Damaged component std dev
    int n_iterations = 0;        // EM iterations to converge
    bool converged = false;      // Did EM converge?
    bool separated = false;      // Are components well-separated? (d_ancient > 0.02)
};

struct GCBinInput {
    float d_max;      // Estimated damage for this bin
    float c_sites;    // Weight (number of C sites)
    bool valid;       // Is this bin valid?
};

class DamageMixtureModel {
public:
    // Fixed parameters
    static constexpr float MU_0 = 0.0f;       // Undamaged component mean
    static constexpr float TAU_0 = 0.01f;     // Undamaged component std dev
    static constexpr float TAU_1 = 0.10f;     // Damaged component std dev (fixed)
    static constexpr float SIGMA_FLOOR = 0.02f; // Minimum observation noise
    static constexpr int MAX_ITER = 50;
    static constexpr float CONVERGENCE_TOL = 1e-6f;

    template<size_t N>
    static DamageMixtureResult fit(const std::array<GCBinInput, N>& bins) {
        DamageMixtureResult result;

        // Compute weighted mean (d_mean)
        float total_weight = 0.0f;
        float weighted_sum = 0.0f;
        int n_valid = 0;
        float max_d = 0.0f;

        for (size_t i = 0; i < N; ++i) {
            if (bins[i].valid && bins[i].c_sites > 0) {
                weighted_sum += bins[i].d_max * bins[i].c_sites;
                total_weight += bins[i].c_sites;
                max_d = std::max(max_d, bins[i].d_max);
                ++n_valid;
            }
        }

        if (total_weight < 1.0f || n_valid < 2) {
            return result;  // Not enough data
        }

        result.d_mean = weighted_sum / total_weight;

        // Initialize EM
        float pi_1 = 0.2f;  // Initial mixing proportion for damaged
        float mu_1 = std::max(0.05f, max_d * 0.8f);  // Initialize near max

        // Precompute observation noise σ_i for each bin
        std::array<float, N> sigma;
        for (size_t i = 0; i < N; ++i) {
            if (bins[i].valid && bins[i].c_sites > 0) {
                float d = std::clamp(bins[i].d_max, 0.001f, 0.999f);
                float sigma_binom = std::sqrt(d * (1.0f - d) / bins[i].c_sites);
                sigma[i] = std::max(SIGMA_FLOOR, sigma_binom);
            } else {
                sigma[i] = SIGMA_FLOOR;
            }
        }

        // EM iterations
        std::array<float, N> r;  // Responsibilities (posterior prob of damaged)
        float prev_ll = -1e30f;
        float last_sum_p = 0.0f;   // B6: final-iteration Sum prec_i, for the mu_1 SE gate

        for (int iter = 0; iter < MAX_ITER; ++iter) {
            // E-step: compute responsibilities using log-space for stability
            float ll = 0.0f;
            for (size_t i = 0; i < N; ++i) {
                if (!bins[i].valid || bins[i].c_sites < 1.0f) {
                    r[i] = 0.0f;
                    continue;
                }

                float d = bins[i].d_max;
                float w = bins[i].c_sites;

                // Variance for each component: τ_k² + σ_i²
                float var_0 = TAU_0 * TAU_0 + sigma[i] * sigma[i];
                float var_1 = TAU_1 * TAU_1 + sigma[i] * sigma[i];

                // Log-likelihood under each component
                float ll_0 = -0.5f * std::log(var_0) - 0.5f * (d - MU_0) * (d - MU_0) / var_0;
                float ll_1 = -0.5f * std::log(var_1) - 0.5f * (d - mu_1) * (d - mu_1) / var_1;

                // Log-posterior (unnormalized)
                float log_p0 = std::log(1.0f - pi_1 + 1e-10f) + ll_0;
                float log_p1 = std::log(pi_1 + 1e-10f) + ll_1;

                // Normalize using log-sum-exp
                float log_sum = log_p0;
                if (log_p1 > log_p0) {
                    log_sum = log_p1 + std::log(1.0f + std::exp(log_p0 - log_p1));
                } else {
                    log_sum = log_p0 + std::log(1.0f + std::exp(log_p1 - log_p0));
                }

                r[i] = std::exp(log_p1 - log_sum);
                ll += w * log_sum;
            }

            // Check convergence
            if (std::abs(ll - prev_ll) < CONVERGENCE_TOL * std::abs(prev_ll)) {
                result.converged = true;
                result.n_iterations = iter + 1;
                break;
            }
            prev_ll = ll;

            // M-step: pi_1 by C-site prevalence; mu_1 by IVW (Gaussian MLE).
            // Weighting mu_1 by 1/(TAU_1²+σ_i²) caps over-influence from deep bins:
            // once c_sites > d(1−d)/SIGMA_FLOOR² ≈ 500, sigma[i]→SIGMA_FLOOR and extra
            // depth no longer buys extra weight.  Consistent with the E-step variances.
            static constexpr float PI_FLOOR   = 1e-3f;
            static constexpr float MIN_MU_EFF = 0.5f;  // ≥0.5 "effective bins" required

            float sum_w = 0.0f, sum_w_r = 0.0f;        // prevalence (c_sites-weighted)
            float sum_p = 0.0f, sum_p_d = 0.0f;        // precision-weighted for mu_1
            float tau1_var = TAU_1 * TAU_1;

            for (size_t i = 0; i < N; ++i) {
                if (bins[i].valid && bins[i].c_sites > 0) {
                    float c    = static_cast<float>(bins[i].c_sites);
                    float var  = tau1_var + sigma[i] * sigma[i];
                    float prec = r[i] / var;
                    sum_w += c;  sum_w_r += c * r[i];
                    sum_p += prec;  sum_p_d += prec * bins[i].d_max;
                }
            }
            last_sum_p = sum_p;   // B6: capture final-iteration precision sum for the mu_1 SE gate
            // Effective-mass check: sum_p * (TAU_1²+SIGMA_FLOOR²) ≥ MIN_MU_EFF
            float eff_mass = sum_p * (tau1_var + SIGMA_FLOOR * SIGMA_FLOOR);

            if (sum_w > 1e-10f)
                pi_1 = std::clamp(sum_w_r / sum_w, PI_FLOOR, 1.0f - PI_FLOOR);
            if (eff_mass >= MIN_MU_EFF)
                mu_1 = std::clamp(sum_p_d / sum_p, 0.0f, 1.0f);

            if (sum_w_r < 1e-10f) {
                // No mass in damaged component - declare no separation
                result.converged = true;
                result.n_iterations = iter + 1;
                result.separated = false;
                result.d_ancient = 0.0f;
                result.pi_ancient = 0.0f;
                return result;
            }

            result.n_iterations = iter + 1;
        }

        // Store results
        result.d_ancient = mu_1;
        result.pi_ancient = pi_1;
        result.tau_ancient = TAU_1;
        // B6: per-library separation gate. mu_1 is the IVW mean (M-step above),
        // SE = sqrt(1/sum_p) (sum_p = Sum r_i/(TAU_1^2+sigma_i^2) already carries
        // precision units). Require mu_1 > 2*SE above zero -- one-sided, so
        // mu_1 > 0.0f guards the sign that mu_1^2*sum_p > 4 would otherwise pass --
        // instead of the fixed 0.02 floor, so a faint lib no longer flips
        // 'separated' on where EM happened to initialise mu_1.
        result.separated = (mu_1 > 0.0f && last_sum_p > 0.0f
                            && mu_1 * mu_1 * last_sum_p > 4.0f && pi_1 > 0.01f);

        return result;
    }
};

} // namespace taph
