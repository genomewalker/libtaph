#pragma once
// Static helpers shared by finalize_sample_profile split TUs.
// FrameSelector::finalize_sample_profile and supporting decay/ctx helpers.
#include "taph/frame_selector_decl.hpp"
#include "taph/codon_tables.hpp"
#include "taph/hexamer_tables.hpp"
#include "taph/library_interpretation.hpp"
#include "taph/channel_count_table.hpp"
#include "taph/channel_registry.hpp"
#include "damage_estimation_detail.hpp"
#include <algorithm>
#include <cmath>
#include <array>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>
#include <array>
namespace taph {
// Negative return value signals inverted pattern (terminal lower than interior).
static float compute_decay_llr(
    const std::array<double, 15>& freq,
    const std::array<double, 15>& total,
    float baseline,
    float /*amplitude*/,
    float lambda) {

    const double MIN_COVERAGE = 100.0;

    // Best-fit amplitude without clamping: allows detecting inverted (negative) patterns
    double num = 0.0, den = 0.0;
    for (int i = 1; i < 10; ++i) {
        if (total[i] < MIN_COVERAGE) continue;
        double x = std::exp(-lambda * i);
        double excess = freq[i] - baseline;
        num += total[i] * x * excess;
        den += total[i] * x * x;
    }
    double raw_amplitude = (den > 0.0) ? (num / den) : 0.0;

    double ll_exp = 0.0;
    double ll_const = 0.0;

    for (int i = 1; i < 10; ++i) {
        if (total[i] < MIN_COVERAGE) continue;

        double n = total[i];
        double k = freq[i] * n;

        double p_exp = baseline + raw_amplitude * std::exp(-lambda * i);
        p_exp = std::clamp(p_exp, 0.001, 0.999);
        ll_exp += binomial_ll(k, n, p_exp);

        double p_const = std::clamp(static_cast<double>(baseline), 0.001, 0.999);
        ll_const += binomial_ll(k, n, p_const);
    }

    float llr = static_cast<float>(ll_exp - ll_const);

    // Negate LLR for inverted patterns so callers see negative values for T-depletion
    if (raw_amplitude < 0) {
        return -std::abs(llr);
    }
    return llr;
}

// Fit p(pos) = b + A * exp(-lambda * pos) via weighted least squares.
// Returns {b, A, lambda, rmse}. If external_baseline >= 0, use it instead of
// estimating from positions 10-14 (middle-of-read baseline is more reliable).
static std::array<float, 4> fit_exponential_decay(
    const std::array<double, 15>& freq,
    const std::array<double, 15>& coverage,
    float lambda_init = 0.2f,
    float external_baseline = -1.0f) {

    const double MIN_COVERAGE = 100.0;

    // Fallback to positions 10-14 when no external baseline is supplied, but
    // those positions still carry end-composition artifacts; b=0 would inflate
    // damage to the raw T/(T+C) ratio causing downstream failure modes.
    float b;
    if (external_baseline >= 0.0f) {
        b = std::clamp(external_baseline, 0.001f, 0.999f);
    } else {
        double baseline_sum = 0.0, baseline_weight = 0.0;
        for (int i = 10; i < 15; ++i) {
            if (coverage[i] >= MIN_COVERAGE) {
                baseline_sum += freq[i] * coverage[i];
                baseline_weight += coverage[i];
            }
        }
        b = (baseline_weight > 0) ?
            static_cast<float>(baseline_sum / baseline_weight) : 0.5f;
        b = std::clamp(b, 0.001f, 0.999f);
    }

    int n_valid = 0;
    for (int i = 1; i < 10; ++i) {
        if (coverage[i] >= MIN_COVERAGE) n_valid++;
    }

    // Amplitude from position 1, not 0: pos 0 is prone to ligation/trimming artifacts
    float A = (coverage[1] >= MIN_COVERAGE) ?
              std::max(0.0f, static_cast<float>(freq[1]) - b) : 0.0f;

    if (n_valid < 5) {
        return {b, A, lambda_init, 1.0f};
    }

    // Refine lambda via log-linear regression on positions 1-9 (exclude pos 0)
    float lambda = lambda_init;
    if (A > 0.01f) {
        double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0, sum_w = 0.0;
        for (int i = 1; i < 10; ++i) {
            if (coverage[i] < MIN_COVERAGE) continue;
            double excess = freq[i] - b;
            if (excess > 0.005) {
                double w = coverage[i];
                double y = std::log(excess);
                sum_x += w * i;
                sum_y += w * y;
                sum_xy += w * i * y;
                sum_xx += w * i * i;
                sum_w += w;
            }
        }
        if (sum_w > 0 && (sum_w * sum_xx - sum_x * sum_x) > 1e-6) {
            double slope = (sum_w * sum_xy - sum_x * sum_y) /
                          (sum_w * sum_xx - sum_x * sum_x);
            lambda = std::max(0.05f, std::min(0.5f, static_cast<float>(-slope)));
            double intercept = (sum_y - slope * sum_x) / sum_w;
            float A_new = static_cast<float>(std::exp(intercept));
            if (A_new > 0.0f && A_new < 1.0f - b) {
                A = A_new;
            }
        }
    }

    A = std::max(0.0f, std::min(A, 1.0f - b - 0.001f));

    double sse = 0.0, weight_sum = 0.0;
    for (int i = 1; i < 15; ++i) {
        if (coverage[i] >= MIN_COVERAGE) {
            double pred = b + A * std::exp(-lambda * i);
            double resid = freq[i] - pred;
            sse += resid * resid * coverage[i];
            weight_sum += coverage[i];
        }
    }
    float rmse = (weight_sum > 0) ?
                 static_cast<float>(std::sqrt(sse / weight_sum)) : 1.0f;

    return {b, A, lambda, rmse};
}

// Fit p(pos) = baseline + amplitude * exp(-lambda * pos) with baseline and lambda fixed.
// Used for 3' library-type classification over positions 1-10 (pos-0 excluded as artifact).
// Amplitude is constrained >= 0 (damage only adds signal, never removes it).
// Returns BIC of the amplitude model vs the flat null (bias-only) model.
struct ChannelDecayFit {
    float    baseline   = 0.0f;
    float    amplitude  = 0.0f;
    float    lambda     = 0.0f;
    // BIC stored in double: at 100M+ reads the values reach ~1e9 and float loses
    // the ~20-unit differences that distinguish models (float epsilon ~256 at 1e9).
    double   log_lik_alt  = 0.0;
    double   log_lik_null = 0.0;
    double   bic_alt   = 0.0;
    double   bic_null  = 0.0;
    double   delta_bic = 0.0;  // bic_null - bic_alt; positive favours decay model
    uint64_t n_trials  = 0;
    bool     valid     = false;
};

static ChannelDecayFit fit_decay_fixed_lambda(
        const std::array<double, 15>& freq,
        const std::array<double, 15>& coverage,
        float baseline,
        float lambda,
        int start_pos = 1,
        int end_pos   = 10,
        int min_valid = 3) {

    ChannelDecayFit fit;
    fit.baseline = std::clamp(baseline, 0.001f, 0.999f);
    // Allow up to 10.0 to capture steep single-position spikes (e.g. adapter-truncated profiles).
    fit.lambda   = std::clamp(lambda,   0.05f,  10.0f);

    // WLS amplitude: A_hat = max(0, Σ n·x·(y - b) / Σ n·x²)
    // Basis is offset-normalised: x = exp(-lambda*(p-start_pos)) so that A is the amplitude
    // at the first fitted position and large lambdas remain numerically well-conditioned.
    double numer = 0.0, denom = 0.0;
    int n_valid = 0;
    for (int p = start_pos; p <= end_pos && p < 15; ++p) {
        double n = coverage[p];
        if (n < 100.0) continue;
        double x = std::exp(-fit.lambda * (p - start_pos));
        numer += n * x * (freq[p] - fit.baseline);
        denom += n * x * x;
        fit.n_trials += static_cast<uint64_t>(n);
        ++n_valid;
    }
    if (n_valid < min_valid || denom <= 0.0 || fit.n_trials == 0) return fit;

    fit.amplitude = std::clamp(static_cast<float>(numer / denom),
                               0.0f, 1.0f - fit.baseline - 0.001f);

    // Binomial log-likelihoods for alt (decay) and null (flat baseline) models
    for (int p = start_pos; p <= end_pos && p < 15; ++p) {
        double n = coverage[p];
        if (n < 100.0) continue;
        double k     = freq[p] * n;
        double x     = std::exp(-fit.lambda * (p - start_pos));
        double p_alt = std::clamp(static_cast<double>(fit.baseline + fit.amplitude * x), 0.001, 0.999);
        double p_null = std::clamp(static_cast<double>(fit.baseline), 0.001, 0.999);
        fit.log_lik_alt  += binomial_ll(k, n, p_alt);
        fit.log_lik_null += binomial_ll(k, n, p_null);
    }

    // BIC: alt has 1 free parameter (amplitude); null has 0.
    // Keep in double — at 100M+ reads BIC values reach ~1e9, losing float precision.
    double log_n = std::log(static_cast<double>(fit.n_trials));
    fit.bic_alt  = -2.0 * fit.log_lik_alt  + log_n;
    fit.bic_null = -2.0 * fit.log_lik_null;
    fit.delta_bic = fit.bic_null - fit.bic_alt;
    fit.valid = true;
    return fit;
}

// Joint two-channel WLS amplitude fit: ct5 + ga3 share ONE amplitude parameter.
// Biologically, DS deamination produces ct5 ≈ ga3 (same exponential decay, symmetric).
// When ct5 << ga3 (SS_comp: complement-orientation reads only, residual original reads),
// A_joint is a poor compromise for both channels → LL(A_joint) < LL_null for ct5 and
// < LL_alt_optimal for ga3 → combined BIC_alt exceeds independent fits → M_DS_symm rejected.
//
// Returns ChannelDecayFit where:
//   amplitude = A_joint (single value fit to both channels simultaneously)
//   bic_alt   = -2*LL_joint(A_joint) + log(n1+n2)   [1 free parameter for both channels]
//   bic_null  = ch1.bic_null + ch2.bic_null          [= -2*LL_joint(0)]
//   delta_bic = bic_null - bic_alt
// valid iff both channels separately had >= min_valid positions with coverage >= 100.
static ChannelDecayFit fit_decay_fixed_lambda_joint(
        const std::array<double, 15>& freq1,
        const std::array<double, 15>& cov1,
        float baseline1,
        const std::array<double, 15>& freq2,
        const std::array<double, 15>& cov2,
        float baseline2,
        float lambda,
        int start_pos = 1,
        int end_pos   = 10,
        int min_valid = 3) {

    ChannelDecayFit fit;
    fit.baseline = std::clamp(baseline1, 0.001f, 0.999f);
    fit.lambda   = std::clamp(lambda, 0.05f, 10.0f);

    // WLS joint amplitude summed over both channels (offset-normalised basis)
    double numer = 0.0, denom = 0.0;
    int n_valid1 = 0, n_valid2 = 0;
    for (int p = start_pos; p <= end_pos && p < 15; ++p) {
        double x = std::exp(-fit.lambda * (p - start_pos));
        double n1 = cov1[p];
        if (n1 >= 100.0) {
            numer += n1 * x * (freq1[p] - baseline1);
            denom += n1 * x * x;
            fit.n_trials += static_cast<uint64_t>(n1);
            ++n_valid1;
        }
        double n2 = cov2[p];
        if (n2 >= 100.0) {
            numer += n2 * x * (freq2[p] - baseline2);
            denom += n2 * x * x;
            fit.n_trials += static_cast<uint64_t>(n2);
            ++n_valid2;
        }
    }
    if (n_valid1 < min_valid || n_valid2 < min_valid || denom <= 0.0 || fit.n_trials == 0) return fit;

    double max_base = std::max(static_cast<double>(baseline1), static_cast<double>(baseline2));
    fit.amplitude = std::clamp(static_cast<float>(numer / denom),
                               0.0f, static_cast<float>(1.0 - max_base - 0.001));
    double A = fit.amplitude;

    // Joint log-likelihoods at A_joint and at null (A=0) for both channels
    for (int p = start_pos; p <= end_pos && p < 15; ++p) {
        double x = std::exp(-fit.lambda * (p - start_pos));
        double n1 = cov1[p];
        if (n1 >= 100.0) {
            double k1    = freq1[p] * n1;
            double p_alt = std::clamp(baseline1 + A * x, 0.001, 0.999);
            double p_null = std::clamp(static_cast<double>(baseline1), 0.001, 0.999);
            fit.log_lik_alt  += binomial_ll(k1, n1, p_alt);
            fit.log_lik_null += binomial_ll(k1, n1, p_null);
        }
        double n2 = cov2[p];
        if (n2 >= 100.0) {
            double k2    = freq2[p] * n2;
            double p_alt = std::clamp(baseline2 + A * x, 0.001, 0.999);
            double p_null = std::clamp(static_cast<double>(baseline2), 0.001, 0.999);
            fit.log_lik_alt  += binomial_ll(k2, n2, p_alt);
            fit.log_lik_null += binomial_ll(k2, n2, p_null);
        }
    }

    // BIC: 1 free parameter (shared amplitude); n = n1 + n2 combined
    double log_n = std::log(static_cast<double>(fit.n_trials));
    fit.bic_alt   = -2.0 * fit.log_lik_alt  + log_n;
    fit.bic_null  = -2.0 * fit.log_lik_null;
    fit.delta_bic = fit.bic_null - fit.bic_alt;
    fit.valid = true;
    return fit;
}

// Try start_pos in {1,2,3} × lambda in {caller_lambda, 2.0, 5.0, 10.0}.
// Returns the fit with the highest delta_bic (best BIC improvement over null).
// High-lambda templates capture single-position spikes caused by adapter remnants.
// GA0 is always at pos 0 and must NOT use this helper.
static std::pair<ChannelDecayFit, int> fit_decay_best_offset(
        const std::array<double, 15>& freq,
        const std::array<double, 15>& coverage,
        float baseline, float lambda,
        int end_pos = 10, int min_valid = 3, int min_sp = 1) {
    static constexpr float lambda_extra[] = {2.0f, 5.0f, 10.0f};
    ChannelDecayFit best;
    int best_offset = min_sp;
    for (int sp = min_sp; sp <= std::max(min_sp, 3); ++sp) {
        // Try caller's lambda
        ChannelDecayFit f = fit_decay_fixed_lambda(freq, coverage, baseline, lambda, sp, end_pos, min_valid);
        if (f.valid && (!best.valid || f.delta_bic > best.delta_bic)) { best = f; best_offset = sp; }
        // Try steep-decay templates
        for (float lam : lambda_extra) {
            f = fit_decay_fixed_lambda(freq, coverage, baseline, lam, sp, end_pos, min_valid);
            if (f.valid && (!best.valid || f.delta_bic > best.delta_bic)) { best = f; best_offset = sp; }
        }
    }
    if (!best.valid) {
        best = fit_decay_fixed_lambda(freq, coverage, baseline, lambda, min_sp, end_pos, min_valid);
        best_offset = min_sp;
    }
    return {best, best_offset};
}

// Joint offset search for the DS symmetric model (CT5 + GA3 share one amplitude).
// Also iterates over high-lambda templates for consistency with the single-channel search.
static std::pair<ChannelDecayFit, int> fit_decay_joint_best_offset(
        const std::array<double, 15>& freq1, const std::array<double, 15>& cov1, float baseline1,
        const std::array<double, 15>& freq2, const std::array<double, 15>& cov2, float baseline2,
        float lambda, int end_pos = 10, int min_valid = 3, bool restrict_high_lambda = false,
        int min_sp = 1) {
    static constexpr float lambda_extra[] = {2.0f, 5.0f, 10.0f};
    ChannelDecayFit best;
    int best_offset = min_sp;
    for (int sp = min_sp; sp <= std::max(min_sp, 3); ++sp) {
        ChannelDecayFit f = fit_decay_fixed_lambda_joint(
            freq1, cov1, baseline1, freq2, cov2, baseline2, lambda, sp, end_pos, min_valid);
        if (f.valid && (!best.valid || f.delta_bic > best.delta_bic)) { best = f; best_offset = sp; }
        if (!restrict_high_lambda) {
            for (float lam : lambda_extra) {
                f = fit_decay_fixed_lambda_joint(
                    freq1, cov1, baseline1, freq2, cov2, baseline2, lam, sp, end_pos, min_valid);
                if (f.valid && (!best.valid || f.delta_bic > best.delta_bic)) { best = f; best_offset = sp; }
            }
        }
    }
    if (!best.valid) {
        best = fit_decay_fixed_lambda_joint(
            freq1, cov1, baseline1, freq2, cov2, baseline2, lambda, min_sp, end_pos, min_valid);
        best_offset = min_sp;
    }
    return {best, best_offset};
}



// Context-split 1D amplitude fit for CpG-like / non-CpG-like C→T damage.
// Fixed lambda (from global fit), golden-section search over d in [0,1].
struct CtCtxFit {
    float baseline = std::numeric_limits<float>::quiet_NaN();
    float dmax     = std::numeric_limits<float>::quiet_NaN();
    // Coverage/effective-coverage are SUMS of per-position double counts. At 100M+
    // libraries these exceed float32's 2^24 exact-integer limit (rounded by ±16 ULP).
    // Keep in double: computed once at finalize from exact double accumulators, so
    // double is both exact and deterministic. (dmax/baseline are bounded rates -> float.)
    double cov_terminal = 0.0, cov_interior = 0.0;
    double effcov_terminal = 0.0, effcov_interior = 0.0;
    int   fit_positions = 0;
    bool  valid = false;
    bool  lower_boundary = false;  // C4: fit ran but dmax clamped to lower wall; valid=false for individual emission, but contributes 0 to contrasts
};

static CtCtxFit fit_ct5_ctx_amplitude(
    const std::array<double, SampleDamageProfile::N_POS>& t_counts,
    const std::array<double, SampleDamageProfile::N_POS>& total_counts,
    double t_interior, double total_interior,
    float lambda)
{
    constexpr float MIN_POS_COV  = 50.0f;
    constexpr float MIN_INT_COV  = 500.0f;
    constexpr float MIN_EFF_INT  = 100.0f;
    constexpr float MIN_EFF_TERM = 50.0f;
    constexpr float EPS = 1e-6f;

    CtCtxFit fit;
    fit.cov_interior = total_interior;
    if (total_interior < MIN_INT_COV) return fit;

    const double b_d = std::clamp(t_interior / total_interior, (double)EPS, 1.0 - (double)EPS);
    const float b = static_cast<float>(b_d);
    fit.baseline = b;
    fit.effcov_interior = total_interior * (1.0 - b_d);
    if (fit.effcov_interior < MIN_EFF_INT) return fit;

    double cov_term = 0.0, effcov_term = 0.0;
    int n_fit_pos = 0;
    for (int p = 0; p < SampleDamageProfile::N_POS; ++p) {
        const double n = total_counts[p];
        cov_term += n;
        if (n >= MIN_POS_COV) { ++n_fit_pos; effcov_term += n * (1.0 - b_d); }
    }
    fit.cov_terminal = cov_term;
    fit.effcov_terminal = effcov_term;
    fit.fit_positions = n_fit_pos;
    if (n_fit_pos < 3 || effcov_term < MIN_EFF_TERM) return fit;

    const float lam = std::clamp(lambda, 0.01f, 2.0f);

    auto neg_ll = [&](double d) -> double {
        d = std::clamp(d, 0.0, 1.0);
        double nll = 0.0;
        for (int p = 0; p < SampleDamageProfile::N_POS; ++p) {
            const double n = static_cast<double>(total_counts[p]);
            if (n < MIN_POS_COV) continue;
            const double k = static_cast<double>(t_counts[p]);
            const double mu = std::clamp(
                static_cast<double>(b) + (1.0 - b) * d * std::exp(-lam * p),
                1e-6, 1.0 - 1e-6);
            nll -= k * std::log(mu) + (n - k) * std::log(1.0 - mu);
        }
        return nll;
    };

    // Golden section search
    double lo = 0.0, hi = 1.0;
    const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
    double c = hi - gr * (hi - lo), d_pt = lo + gr * (hi - lo);
    for (int iter = 0; iter < 60; ++iter) {
        if (neg_ll(c) < neg_ll(d_pt)) { hi = d_pt; d_pt = c; c = hi - gr * (hi - lo); }
        else                           { lo = c;    c = d_pt; d_pt = lo + gr * (hi - lo); }
        if (hi - lo < 1e-7) break;
    }
    fit.dmax = static_cast<float>((lo + hi) / 2.0);
    // Boundary check: d near 1.0 means the optimizer hit the upper wall —
    // the true optimum is outside [0,1] or the signal is indistinguishable
    // from noise. Report as invalid rather than a misleading saturated value.
    // C4: also reject the LOWER wall — with no terminal elevation in this context
    // the search converges to dmax≈0.0 with valid=true, emitting a numeric 0.0 that
    // is indistinguishable from a genuine zero. Treat a lower-boundary fit as
    // undetectable signal so the consumer sees null (the dmax==NaN guard downstream).
    fit.lower_boundary = (fit.dmax <= 1e-4f);
    fit.valid = (!fit.lower_boundary && fit.dmax < 0.98f);
    return fit;
}

// Empirical GG-breakpoint proxy: GG-vs-A endpoint double-difference (reference-free,
// composition-internal). Trinuc idx = prev*16 + mid*4 + next (A=0,C=1,G=2,T=3). 5'-of-GG run =
// mid=G,next=G  => idx%16 == 10. A-context comparator = mid=A => (idx/4)%4 == 0.
static float compute_oxscission_delta(const std::array<uint64_t, 64>& term,
                                      const std::array<uint64_t, 64>& inter) {
    auto frac = [](const std::array<uint64_t, 64>& s, bool gg) -> double {
        uint64_t num = 0, den = 0;
        for (int i = 0; i < 64; ++i) {
            den += s[i];
            const bool hit = gg ? ((i % 16) == 10) : (((i / 4) % 4) == 0);
            if (hit) num += s[i];
        }
        return den ? static_cast<double>(num) / static_cast<double>(den) : 0.0;
    };
    const double gg = frac(term, true)  - frac(inter, true);
    const double a  = frac(term, false) - frac(inter, false);
    return static_cast<float>(gg - a);
}

} // namespace taph
