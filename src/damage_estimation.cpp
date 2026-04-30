// Sample-level damage profile management

#include "taph/frame_selector_decl.hpp"
#include "taph/codon_tables.hpp"
#include "taph/hexamer_tables.hpp"
#include "taph/library_interpretation.hpp"
#include <algorithm>
#include <cmath>
#include <array>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace taph {

// Compute binomial log-likelihood for a single observation
// k = successes (e.g., T count), n = total trials (e.g., T+C count), p = probability
static inline double binomial_ll(double k, double n, double p) {
    if (n < 1 || p <= 0 || p >= 1) return 0.0;
    // Log-likelihood: k*log(p) + (n-k)*log(1-p) + constant (ignored for LLR)
    return k * std::log(p) + (n - k) * std::log(1.0 - p);
}

// LLR of exponential decay vs constant model over positions 1-9 (excludes pos 0 artifacts).
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
        int end_pos = 10, int min_valid = 3) {
    static constexpr float lambda_extra[] = {2.0f, 5.0f, 10.0f};
    ChannelDecayFit best;
    int best_offset = 1;
    for (int sp = 1; sp <= 3; ++sp) {
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
        best = fit_decay_fixed_lambda(freq, coverage, baseline, lambda, 1, end_pos, min_valid);
        best_offset = 1;
    }
    return {best, best_offset};
}

// Joint offset search for the DS symmetric model (CT5 + GA3 share one amplitude).
// Also iterates over high-lambda templates for consistency with the single-channel search.
static std::pair<ChannelDecayFit, int> fit_decay_joint_best_offset(
        const std::array<double, 15>& freq1, const std::array<double, 15>& cov1, float baseline1,
        const std::array<double, 15>& freq2, const std::array<double, 15>& cov2, float baseline2,
        float lambda, int end_pos = 10, int min_valid = 3, bool restrict_high_lambda = false) {
    static constexpr float lambda_extra[] = {2.0f, 5.0f, 10.0f};
    ChannelDecayFit best;
    int best_offset = 1;
    for (int sp = 1; sp <= 3; ++sp) {
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
            freq1, cov1, baseline1, freq2, cov2, baseline2, lambda, 1, end_pos, min_valid);
        best_offset = 1;
    }
    return {best, best_offset};
}

void FrameSelector::update_sample_profile(
    SampleDamageProfile& profile,
    std::string_view seq) {

    if (seq.length() < 30) return;  // Too short for reliable statistics

    size_t len = seq.length();

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        char base = fast_upper(seq[i]);
        // Damage signal: T/(T+C) for C→T
        if (base == 'T') {
            profile.t_freq_5prime[i]++;
            profile.tc_total_5prime[i]++;
        } else if (base == 'C') {
            profile.c_freq_5prime[i]++;
            profile.tc_total_5prime[i]++;
        }
        // Negative control: A/(A+G) - should NOT be elevated by C→T damage
        if (base == 'A') {
            profile.a_freq_5prime[i]++;
        } else if (base == 'G') {
            profile.g_freq_5prime[i]++;
        }
    }

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        // Damage signal: A/(A+G) for G→A
        if (base == 'A') {
            profile.a_freq_3prime[i]++;
            profile.ag_total_3prime[i]++;
        } else if (base == 'G') {
            profile.g_freq_3prime[i]++;
            profile.ag_total_3prime[i]++;
        }
        // Negative control: T/(T+C) - should NOT be elevated by G→A damage
        if (base == 'T') {
            profile.t_freq_3prime[i]++;
        } else if (base == 'C') {
            profile.c_freq_3prime[i]++;
        }
    }

    // Tail-anchored background sampling: track C->T (G->A) rates at
    // positions BG_TAIL_LO..BG_TAIL_HI from each terminus to provide a
    // chemistry-robust baseline. Only fills positions actually covered by
    // this read.
    {
        const int lo = SampleDamageProfile::BG_TAIL_LO;
        const int hi = SampleDamageProfile::BG_TAIL_HI;
        for (int i = lo; i <= hi && static_cast<size_t>(i) < len; ++i) {
            const int idx = i - lo;
            const char b5 = fast_upper(seq[i]);
            if (b5 == 'T') { profile.tail_t_5prime[idx]++; profile.tail_tc_5prime[idx]++; }
            else if (b5 == 'C') { profile.tail_tc_5prime[idx]++; }

            const size_t pos3 = len - 1 - i;
            const char b3 = fast_upper(seq[pos3]);
            if (b3 == 'A') { profile.tail_a_3prime[idx]++; profile.tail_ag_3prime[idx]++; }
            else if (b3 == 'G') { profile.tail_ag_3prime[idx]++; }
        }
    }

    // Count bases in middle third (undamaged baseline)
    constexpr size_t INTERIOR_TERM_PAD = 15;
    size_t mid_start = len / 3;
    size_t mid_end   = 2 * len / 3;
    if (mid_start < INTERIOR_TERM_PAD)     mid_start = INTERIOR_TERM_PAD;
    if (len > INTERIOR_TERM_PAD && mid_end + INTERIOR_TERM_PAD > len)
        mid_end = len - INTERIOR_TERM_PAD;
    const bool interior_safe = (mid_start < mid_end);
    if (interior_safe) {
        for (size_t i = mid_start; i < mid_end; ++i) {
            char base = fast_upper(seq[i]);
            if (base == 'T') profile.baseline_t_freq++;
            else if (base == 'C') profile.baseline_c_freq++;
            else if (base == 'A') profile.baseline_a_freq++;
            else if (base == 'G') profile.baseline_g_freq++;
        }
    }

    // Codon-position-aware counting at 5' end (first 15 bases)
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        int codon_pos = i % 3;  // 0, 1, 2
        char base = fast_upper(seq[i]);
        if (base == 'T') profile.codon_pos_t_count_5prime[codon_pos]++;
        else if (base == 'C') profile.codon_pos_c_count_5prime[codon_pos]++;
    }

    // Codon-position-aware counting at 3' end (last 15 bases)
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        int codon_pos = (len - 1 - i) % 3;
        char base = fast_upper(seq[pos]);
        if (base == 'A') profile.codon_pos_a_count_3prime[codon_pos]++;
        else if (base == 'G') profile.codon_pos_g_count_3prime[codon_pos]++;
    }

    // CpG context damage tracking (5' end, first 5 bases)
    for (size_t i = 0; i < std::min(size_t(5), len - 1); ++i) {
        char base = fast_upper(seq[i]);
        char next = fast_upper(seq[i + 1]);

        // Check for CpG context: C followed by G, or T followed by G (damaged CpG)
        if (next == 'G') {
            if (base == 'C') {
                profile.cpg_c_count++;
            } else if (base == 'T') {
                profile.cpg_t_count++;  // Likely C→T in CpG
            }
        } else {
            // Non-CpG context (same position range as CpG)
            if (base == 'C') {
                profile.non_cpg_c_count++;
            } else if (base == 'T') {
                profile.non_cpg_t_count++;
            }
        }
    }

    // Reference-free trinucleotide spectrum (64 contexts).
    // Terminal zone = read positions 1..4 (both flanks available, inside damage zone).
    // Interior zone = read positions 10..14 (null-distribution baseline).
    // Mirror counters at 3' end use positions counted from the 3' terminus.
    {
        auto nuc_idx = [](char c) -> int {
            switch (c) { case 'A': return 0; case 'C': return 1;
                         case 'G': return 2; case 'T': return 3; }
            return -1;
        };
        auto add_ctx = [&](int prev_pos, int mid_pos, int next_pos,
                           std::array<uint64_t, SampleDamageProfile::N_TRINUC>& target) {
            if (prev_pos < 0 || next_pos >= static_cast<int>(len)) return;
            int i0 = nuc_idx(fast_upper(seq[prev_pos]));
            int i1 = nuc_idx(fast_upper(seq[mid_pos]));
            int i2 = nuc_idx(fast_upper(seq[next_pos]));
            if (i0 < 0 || i1 < 0 || i2 < 0) return;
            ++target[i0 * 16 + i1 * 4 + i2];
        };
        for (int p = 1; p <= 4 && p + 1 < static_cast<int>(len); ++p)
            add_ctx(p - 1, p, p + 1, profile.tri_5prime_terminal);
        for (int p = 10; p <= 14 && p + 1 < static_cast<int>(len); ++p)
            add_ctx(p - 1, p, p + 1, profile.tri_5prime_interior);
        for (int p = 1; p <= 4; ++p) {
            int mid = static_cast<int>(len) - 1 - p;
            add_ctx(mid - 1, mid, mid + 1, profile.tri_3prime_terminal);
        }
        for (int p = 10; p <= 14; ++p) {
            int mid = static_cast<int>(len) - 1 - p;
            add_ctx(mid - 1, mid, mid + 1, profile.tri_3prime_interior);
        }
    }

    // CpG-like context split — 5' terminal positions (all 15)
    // Also accumulate upstream-context-aware bins (AC, CC, GC, TC)
    for (int p = 0; p < SampleDamageProfile::N_POS && (p + 1) < static_cast<int>(len); ++p) {
        const char x = fast_upper(seq[p]);
        const char y = fast_upper(seq[p + 1]);
        if ((x == 'C' || x == 'T') && (y == 'A' || y == 'C' || y == 'G' || y == 'T')) {
            const int ctx = (y == 'G') ? SampleDamageProfile::CPG_LIKE : SampleDamageProfile::NONCPG_LIKE;
            profile.ct_ctx_total_5prime[ctx][p] += 1.0f;
            if (x == 'T') profile.ct_ctx_t_5prime[ctx][p] += 1.0f;
        }
        // Upstream-context-aware: classify by preceding base (for p > 0)
        if (p > 0 && (x == 'C' || x == 'T')) {
            const char u = fast_upper(seq[p - 1]);  // upstream base
            int uctx = -1;
            switch (u) {
                case 'A': uctx = SampleDamageProfile::CTX_AC; break;
                case 'C': uctx = SampleDamageProfile::CTX_CC; break;
                case 'G': uctx = SampleDamageProfile::CTX_GC; break;
                case 'T': uctx = SampleDamageProfile::CTX_TC; break;
            }
            if (uctx >= 0) {
                profile.ct5_total_by_upstream[uctx][p] += 1.0;
                if (x == 'T') profile.ct5_t_by_upstream[uctx][p] += 1.0;
            }
        }
    }

    // Interior baseline + oxoG 16-context (only for reads >= 30 bp, already guarded above)
    {
        constexpr size_t INTERIOR_TERM_PAD_CTX = 15;
        size_t q0s = len / 3;
        size_t q1s = 2 * len / 3;
        if (q0s < INTERIOR_TERM_PAD_CTX) q0s = INTERIOR_TERM_PAD_CTX;
        if (len > INTERIOR_TERM_PAD_CTX && q1s + INTERIOR_TERM_PAD_CTX > len)
            q1s = len - INTERIOR_TERM_PAD_CTX;
        const bool interior_safe_ctx = (q0s < q1s);
        const int q0 = static_cast<int>(q0s), q1 = static_cast<int>(q1s);

        // Context-split interior baseline
        if (interior_safe_ctx)
        for (int q = q0; q < q1 && (q + 1) < static_cast<int>(len); ++q) {
            const char x = fast_upper(seq[q]), y = fast_upper(seq[q + 1]);
            if ((x == 'C' || x == 'T') && (y == 'A' || y == 'C' || y == 'G' || y == 'T')) {
                const int ctx = (y == 'G') ? SampleDamageProfile::CPG_LIKE : SampleDamageProfile::NONCPG_LIKE;
                profile.ct_ctx_total_interior[ctx] += 1.0f;
                if (x == 'T') profile.ct_ctx_t_interior[ctx] += 1.0f;
            }
            // Upstream-context-aware interior baseline
            if (q > 0 && (x == 'C' || x == 'T')) {
                const char u = fast_upper(seq[q - 1]);
                int uctx = -1;
                switch (u) {
                    case 'A': uctx = SampleDamageProfile::CTX_AC; break;
                    case 'C': uctx = SampleDamageProfile::CTX_CC; break;
                    case 'G': uctx = SampleDamageProfile::CTX_GC; break;
                    case 'T': uctx = SampleDamageProfile::CTX_TC; break;
                }
                if (uctx >= 0) {
                    profile.ct5_total_interior_by_upstream[uctx] += 1.0;
                    if (x == 'T') profile.ct5_t_interior_by_upstream[uctx] += 1.0;
                }
            }
        }

        // oxoG 16-context interior panel
        if (interior_safe_ctx)
        for (int q = q0; q < q1; ++q) {
            if (q <= 0 || q >= static_cast<int>(len) - 1) continue;
            const char l = fast_upper(seq[q-1]), b = fast_upper(seq[q]), r = fast_upper(seq[q+1]);
            auto enc = [](char c) -> int {
                switch(c){ case 'A':return 0; case 'C':return 1; case 'G':return 2; case 'T':return 3; default:return -1; }
            };
            auto rc_base = [](char c) -> char {
                switch(c){ case 'A':return 'T'; case 'T':return 'A'; case 'C':return 'G'; case 'G':return 'C'; default:return 'N'; }
            };
            if (b == 'T') {
                int il = enc(l), ir = enc(r);
                if (il >= 0 && ir >= 0) profile.oxog16_t[4*il+ir] += 1.0f;
            } else if (b == 'A') {
                int il = enc(rc_base(r)), ir = enc(rc_base(l));
                if (il >= 0 && ir >= 0) profile.oxog16_a_rc[4*il+ir] += 1.0f;
            }
        }
    }

    if (len >= 18) {
        for (int frame = 0; frame < 3; ++frame) {
            // Scan codons from 5' end up to position 14
            for (size_t k = 0; ; ++k) {
                size_t codon_start = frame + 3 * k;
                if (codon_start + 3 > len || codon_start > 14) break;

                size_t p = codon_start;
                if (p >= 15) break;

                char b0 = fast_upper(seq[codon_start]);
                char b1 = fast_upper(seq[codon_start + 1]);
                char b2 = fast_upper(seq[codon_start + 2]);

                if ((b0 != 'A' && b0 != 'C' && b0 != 'G' && b0 != 'T') ||
                    (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') ||
                    (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T')) {
                    continue;
                }

                profile.total_codons_5prime[p]++;

                if (b1 == 'A' && b2 == 'A') {
                    if (b0 == 'C') profile.convertible_caa_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_taa_5prime[p]++;
                }
                if (b1 == 'A' && b2 == 'G') {
                    if (b0 == 'C') profile.convertible_cag_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tag_5prime[p]++;
                }
                if (b1 == 'G' && b2 == 'A') {
                    if (b0 == 'C') profile.convertible_cga_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tga_5prime[p]++;
                }
            }
        }

        // Track interior convertible codons (positions 30+ from start)
        // This gives us the baseline stop conversion rate
        // Guard: len >= 63 required to have valid interior region [30, len-30)
        // Without this, len - 30 underflows for short reads causing OOB access
        constexpr size_t INTERIOR_MIN_LEN = 63;  // 30 + 3 + 30
        if (len >= INTERIOR_MIN_LEN) {
            for (int frame = 0; frame < 3; ++frame) {
                for (size_t k = 0; ; ++k) {
                    size_t codon_start = frame + 3 * k;
                    // Interior region: 30 to len-30 (away from both ends)
                    if (codon_start < 30 || codon_start + 3 > len - 30) {
                        if (codon_start + 3 > len - 30) break;
                        continue;
                    }

                    char b0 = fast_upper(seq[codon_start]);
                    char b1 = fast_upper(seq[codon_start + 1]);
                    char b2 = fast_upper(seq[codon_start + 2]);

                    if ((b0 != 'A' && b0 != 'C' && b0 != 'G' && b0 != 'T') ||
                        (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') ||
                        (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T')) {
                        continue;
                    }

                    profile.total_codons_interior++;

                    if (b1 == 'A' && b2 == 'A') {
                        if (b0 == 'C') profile.convertible_caa_interior++;
                        else if (b0 == 'T') profile.convertible_taa_interior++;
                    }
                    if (b1 == 'A' && b2 == 'G') {
                        if (b0 == 'C') profile.convertible_cag_interior++;
                        else if (b0 == 'T') profile.convertible_tag_interior++;
                    }
                    if (b1 == 'G' && b2 == 'A') {
                        if (b0 == 'C') profile.convertible_cga_interior++;
                        else if (b0 == 'T') profile.convertible_tga_interior++;
                    }
                }
            }
        }
    }

    // Channel C: oxidative G→T stop codon tracking (GAG→TAG, GAA→TAA, GGA→TGA)
    if (len >= 18) {
        for (int frame = 0; frame < 3; ++frame) {
            for (size_t k = 0; ; ++k) {
                size_t codon_start = frame + 3 * k;
                if (codon_start + 3 > len || codon_start > 14) break;

                size_t p = codon_start;
                if (p >= 15) break;

                char b0 = fast_upper(seq[codon_start]);
                char b1 = fast_upper(seq[codon_start + 1]);
                char b2 = fast_upper(seq[codon_start + 2]);

                if ((b0 != 'A' && b0 != 'C' && b0 != 'G' && b0 != 'T') ||
                    (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') ||
                    (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T')) {
                    continue;
                }

                if (b1 == 'A' && b2 == 'G') {
                    if (b0 == 'G') profile.convertible_gag_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tag_ox_5prime[p]++;
                }
                if (b1 == 'A' && b2 == 'A') {
                    if (b0 == 'G') profile.convertible_gaa_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_taa_ox_5prime[p]++;
                }
                if (b1 == 'G' && b2 == 'A') {
                    if (b0 == 'G') profile.convertible_gga_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tga_ox_5prime[p]++;
                }
            }
        }

        // Track interior oxidative convertible codons (positions 30+)
        constexpr size_t INTERIOR_MIN_LEN_OX = 63;
        if (len >= INTERIOR_MIN_LEN_OX) {
            for (int frame = 0; frame < 3; ++frame) {
                for (size_t k = 0; ; ++k) {
                    size_t codon_start = frame + 3 * k;
                    if (codon_start < 30 || codon_start + 3 > len - 30) {
                        if (codon_start + 3 > len - 30) break;
                        continue;
                    }

                    char b0 = fast_upper(seq[codon_start]);
                    char b1 = fast_upper(seq[codon_start + 1]);
                    char b2 = fast_upper(seq[codon_start + 2]);

                    if ((b0 != 'A' && b0 != 'C' && b0 != 'G' && b0 != 'T') ||
                        (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') ||
                        (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T')) {
                        continue;
                    }

                    if (b1 == 'A' && b2 == 'G') {
                        if (b0 == 'G') profile.convertible_gag_interior++;
                        else if (b0 == 'T') profile.convertible_tag_ox_interior++;
                    }
                    if (b1 == 'A' && b2 == 'A') {
                        if (b0 == 'G') profile.convertible_gaa_interior++;
                        else if (b0 == 'T') profile.convertible_taa_ox_interior++;
                    }
                    if (b1 == 'G' && b2 == 'A') {
                        if (b0 == 'G') profile.convertible_gga_interior++;
                        else if (b0 == 'T') profile.convertible_tga_ox_interior++;
                    }
                }
            }
        }
    }

    // Channel D: G→T and C→A transversion tracking (8-oxoG, Chargaff-balance cross-check).
    // Accumulate raw G, T, C, A counts at 5' terminal positions (0-14).
    // T/(T+G) and A/(A+C) at interior vs terminal positions detect G→T oxidation without alignment.
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        char base = fast_upper(seq[i]);
        if      (base == 'G') profile.g_count_5prime[i]++;
        else if (base == 'T') profile.t_from_g_5prime[i]++;
        if      (base == 'C') profile.c_count_ox_5prime[i]++;
        else if (base == 'A') profile.a_from_c_5prime[i]++;
    }
    // Interior baseline: T/(T+G) and A/(A+C) in middle third (undamaged reference).
    {
        constexpr size_t INTERIOR_TERM_PAD_D = 15;
        size_t mid_s = len / 3, mid_e = 2 * len / 3;
        if (mid_s < INTERIOR_TERM_PAD_D) mid_s = INTERIOR_TERM_PAD_D;
        if (len > INTERIOR_TERM_PAD_D && mid_e + INTERIOR_TERM_PAD_D > len)
            mid_e = len - INTERIOR_TERM_PAD_D;
        if (mid_s < mid_e) {
            for (size_t i = mid_s; i < mid_e; ++i) {
                char base = fast_upper(seq[i]);
                if      (base == 'G') profile.baseline_g_total++;
                else if (base == 'T') profile.baseline_g_to_t_count++;
                if      (base == 'C') profile.baseline_c_ox_total++;
                else if (base == 'A') profile.baseline_c_to_a_count++;
            }
        }
    }

    // Channel E: purine enrichment tracked via a_freq/g_freq; computed in finalize_sample_profile()

    if (len >= 30) {
        // Compute interior GC from positions 5 to end-5 (avoid both terminal regions)
        size_t interior_start = 5;
        size_t interior_end = len > 10 ? len - 5 : len;
        uint64_t gc_count = 0;
        uint64_t at_count = 0;
        for (size_t i = interior_start; i < interior_end; ++i) {
            char base = fast_upper(seq[i]);
            if (base == 'G' || base == 'C') gc_count++;
            else if (base == 'A' || base == 'T') at_count++;
        }

        if (gc_count + at_count > 0) {
            float gc_frac = static_cast<float>(gc_count) / (gc_count + at_count);
            int bin_idx = std::min(static_cast<int>(gc_frac * SampleDamageProfile::N_GC_BINS),
                                   SampleDamageProfile::N_GC_BINS - 1);
            auto& bin = profile.gc_bins[bin_idx];

            // Accumulate Channel A counts (T and C at terminal positions)
            for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
                char base = fast_upper(seq[i]);
                if (base == 'T') bin.t_counts[i]++;
                else if (base == 'C') bin.c_counts[i]++;
                else if (base == 'A') bin.a_counts[i]++;
                else if (base == 'G') bin.g_counts[i]++;
            }

            // Mirror at 3' end: i=0 is last base, i=1 second-to-last, ...
            // Used to recover C→T damage on 3' end for SS libraries and
            // G→A damage on 3' end for DS libraries at joint-mixture time.
            for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
                char base = fast_upper(seq[len - 1 - i]);
                if (base == 'T') bin.t_counts_3prime[i]++;
                else if (base == 'C') bin.c_counts_3prime[i]++;
                else if (base == 'A') bin.a_counts_3prime[i]++;
                else if (base == 'G') bin.g_counts_3prime[i]++;
            }

            // Accumulate interior baselines for both the signal and control channels.
            constexpr size_t INTERIOR_TERM_PAD_GC = 15;
            size_t mid_start = len / 3;
            size_t mid_end = 2 * len / 3;
            if (mid_start < INTERIOR_TERM_PAD_GC) mid_start = INTERIOR_TERM_PAD_GC;
            if (len > INTERIOR_TERM_PAD_GC && mid_end + INTERIOR_TERM_PAD_GC > len)
                mid_end = len - INTERIOR_TERM_PAD_GC;
            if (mid_start < mid_end) {
                for (size_t i = mid_start; i < mid_end; ++i) {
                    char base = fast_upper(seq[i]);
                    if (base == 'T') bin.t_interior++;
                    else if (base == 'C') bin.c_interior++;
                    else if (base == 'A') bin.a_interior++;
                    else if (base == 'G') bin.g_interior++;
                }
            }

            // Accumulate Channel B counts (stop codons at terminal positions)
            for (int frame = 0; frame < 3; ++frame) {
                for (size_t k = 0; ; ++k) {
                    size_t codon_start = frame + 3 * k;
                    if (codon_start + 3 > len || codon_start > 14) break;
                    size_t p = codon_start;

                    char b0 = fast_upper(seq[codon_start]);
                    char b1 = fast_upper(seq[codon_start + 1]);
                    char b2 = fast_upper(seq[codon_start + 2]);

                    // CAA/TAA, CAG/TAG, CGA/TGA
                    if (b1 == 'A' && b2 == 'A') {
                        if (b0 == 'C') bin.pre_counts[p]++;
                        else if (b0 == 'T') bin.stop_counts[p]++;
                    }
                    if (b1 == 'A' && b2 == 'G') {
                        if (b0 == 'C') bin.pre_counts[p]++;
                        else if (b0 == 'T') bin.stop_counts[p]++;
                    }
                    if (b1 == 'G' && b2 == 'A') {
                        if (b0 == 'C') bin.pre_counts[p]++;
                        else if (b0 == 'T') bin.stop_counts[p]++;
                    }
                }
            }

            // Accumulate Channel B interior baseline
            // Guard: len >= 63 required to prevent unsigned underflow in len - 30
            if (len >= 63) {
                for (int frame = 0; frame < 3; ++frame) {
                    for (size_t k = 0; ; ++k) {
                        size_t codon_start = frame + 3 * k;
                        if (codon_start < 30 || codon_start + 3 > len - 30) {
                            if (codon_start + 3 > len - 30) break;
                            continue;
                        }

                        char b0 = fast_upper(seq[codon_start]);
                        char b1 = fast_upper(seq[codon_start + 1]);
                        char b2 = fast_upper(seq[codon_start + 2]);

                        if (b1 == 'A' && b2 == 'A') {
                            if (b0 == 'C') bin.pre_interior++;
                            else if (b0 == 'T') bin.stop_interior++;
                        }
                        if (b1 == 'A' && b2 == 'G') {
                            if (b0 == 'C') bin.pre_interior++;
                            else if (b0 == 'T') bin.stop_interior++;
                        }
                        if (b1 == 'G' && b2 == 'A') {
                            if (b0 == 'C') bin.pre_interior++;
                            else if (b0 == 'T') bin.stop_interior++;
                        }
                    }
                }
            }

            bin.n_reads++;
        }
    }

    if (len >= 18) {
        // 5' terminal hexamer starting at position 0
        char hex_5prime[7];
        bool valid_5prime = true;
        for (int i = 0; i < 6; ++i) {
            char b = fast_upper(seq[i]);
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') {
                valid_5prime = false;
                break;
            }
            hex_5prime[i] = b;
        }
        hex_5prime[6] = '\0';

        if (valid_5prime) {
            uint32_t code_5prime = encode_hexamer(hex_5prime);
            if (code_5prime < 4096) {
                profile.hexamer_count_5prime[code_5prime] += 1.0;
                profile.n_hexamers_5prime++;
            }
        }

        // Interior hexamer (starting at len/2 - 3, approximately middle)
        size_t interior_start = len / 2 - 3;
        char hex_interior[7];
        bool valid_interior = true;
        for (int i = 0; i < 6; ++i) {
            char b = fast_upper(seq[interior_start + i]);
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') {
                valid_interior = false;
                break;
            }
            hex_interior[i] = b;
        }
        hex_interior[6] = '\0';

        if (valid_interior) {
            uint32_t code_interior = encode_hexamer(hex_interior);
            if (code_interior < 4096) {
                profile.hexamer_count_interior[code_interior] += 1.0;
                profile.n_hexamers_interior++;
            }
        }
    }

    // Interior clustered C→T: excess co-occurrence of T at non-CpG {C,T} sites
    if (len >= 30) {
        const int lo = static_cast<int>(len / 3);
        const int hi = static_cast<int>(len - (len / 3));
        const int wlen = hi - lo;
        if (wlen >= 2) {
            auto& acc = profile.interior_ct_cluster;
            // Build eligible + indicator arrays for CT and AG tracks.
            // Thread-local scratch reused across reads — eliminates 4 heap
            // allocations per read on this hot path. Buffers grow monotonically
            // and are zeroed in the prefix actually used (size_t wlen).
            thread_local std::vector<uint8_t> ct_elig_buf, ct_pos_buf,
                                              ag_elig_buf, ag_pos_buf;
            const size_t W = static_cast<size_t>(wlen);
            if (ct_elig_buf.size() < W) {
                ct_elig_buf.resize(W);
                ct_pos_buf .resize(W);
                ag_elig_buf.resize(W);
                ag_pos_buf .resize(W);
            }
            uint8_t* ct_elig = ct_elig_buf.data();
            uint8_t* ct_pos  = ct_pos_buf.data();
            uint8_t* ag_elig = ag_elig_buf.data();
            uint8_t* ag_pos  = ag_pos_buf.data();
            std::memset(ct_elig, 0, W);
            std::memset(ct_pos,  0, W);
            std::memset(ag_elig, 0, W);
            std::memset(ag_pos,  0, W);
            int n_ct = 0, k_ct = 0, n_ag = 0, k_ag = 0;
            for (int i = lo; i < hi; ++i) {
                const int j = i - lo;
                const char b = fast_upper(seq[i]);
                const char nx = (i + 1 < static_cast<int>(len)) ? fast_upper(seq[i+1]) : 'N';
                const char pv = (i > 0) ? fast_upper(seq[i-1]) : 'N';
                if ((b == 'C' || b == 'T') && nx != 'G') {
                    ct_elig[j] = 1; ++n_ct;
                    if (b == 'T') { ct_pos[j] = 1; ++k_ct; }
                }
                if ((b == 'A' || b == 'G') && pv != 'C') {
                    ag_elig[j] = 1; ++n_ag;
                    if (b == 'A') { ag_pos[j] = 1; ++k_ag; }
                }
            }
            if (n_ct >= 2) {
                ++acc.reads_used_ct;
                const double q_ct = (k_ct >= 2)
                    ? (static_cast<double>(k_ct) * (k_ct - 1)) /
                      (static_cast<double>(n_ct) * (n_ct - 1))
                    : 0.0;
                for (int d = 1; d <= 10; ++d) {
                    uint64_t pairs = 0, obs = 0;
                    for (int j = 0; j + d < wlen; ++j) {
                        if (!ct_elig[j] || !ct_elig[j + d]) continue;
                        ++pairs;
                        obs += static_cast<uint64_t>(ct_pos[j] & ct_pos[j + d]);
                    }
                    acc.pairs_ct[d] += pairs;
                    acc.obs_ct[d]   += obs;
                    acc.exp_ct[d]   += static_cast<double>(pairs) * q_ct;
                    acc.var_ct[d]   += static_cast<double>(pairs) * q_ct * (1.0 - q_ct);
                }
            }
            if (n_ag >= 2) {
                ++acc.reads_used_ag;
                const double q_ag = (k_ag >= 2)
                    ? (static_cast<double>(k_ag) * (k_ag - 1)) /
                      (static_cast<double>(n_ag) * (n_ag - 1))
                    : 0.0;
                for (int d = 1; d <= 10; ++d) {
                    uint64_t pairs = 0, obs = 0;
                    for (int j = 0; j + d < wlen; ++j) {
                        if (!ag_elig[j] || !ag_elig[j + d]) continue;
                        ++pairs;
                        obs += static_cast<uint64_t>(ag_pos[j] & ag_pos[j + d]);
                    }
                    acc.pairs_ag[d] += pairs;
                    acc.obs_ag[d]   += obs;
                    acc.exp_ag[d]   += static_cast<double>(pairs) * q_ag;
                    acc.var_ag[d]   += static_cast<double>(pairs) * q_ag * (1.0 - q_ag);
                }
            }
        } else {
            ++profile.interior_ct_cluster.short_reads_skipped;
        }
    } else {
        ++profile.interior_ct_cluster.short_reads_skipped;
    }

    profile.n_reads++;
}

// Context-split 1D amplitude fit for CpG-like / non-CpG-like C→T damage.
// Fixed lambda (from global fit), golden-section search over d in [0,1].
struct CtCtxFit {
    float baseline = std::numeric_limits<float>::quiet_NaN();
    float dmax     = std::numeric_limits<float>::quiet_NaN();
    float cov_terminal = 0.0f, cov_interior = 0.0f;
    float effcov_terminal = 0.0f, effcov_interior = 0.0f;
    int   fit_positions = 0;
    bool  valid = false;
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
    fit.cov_interior = static_cast<float>(total_interior);
    if (total_interior < MIN_INT_COV) return fit;

    const double b_d = std::clamp(t_interior / total_interior, (double)EPS, 1.0 - (double)EPS);
    const float b = static_cast<float>(b_d);
    fit.baseline = b;
    fit.effcov_interior = static_cast<float>(total_interior * (1.0 - b_d));
    if (fit.effcov_interior < MIN_EFF_INT) return fit;

    float cov_term = 0.0f, effcov_term = 0.0f;
    int n_fit_pos = 0;
    for (int p = 0; p < SampleDamageProfile::N_POS; ++p) {
        const double n = total_counts[p];
        cov_term += static_cast<float>(n);
        if (n >= MIN_POS_COV) { ++n_fit_pos; effcov_term += static_cast<float>(n * (1.0 - b_d)); }
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
    fit.valid = (fit.dmax < 0.98f);
    return fit;
}

void FrameSelector::finalize_sample_profile(SampleDamageProfile& profile) {
    if (profile.n_reads == 0) return;
    // Lifecycle guard (see SampleDamageProfile::finalized): finalize mutates
    // raw counts into rates in place — calling it twice would normalize the
    // already-normalized arrays again. Throw rather than silently corrupt.
    if (profile.finalized) {
        throw std::logic_error(
            "finalize_sample_profile: profile already finalized; "
            "reset_sample_profile() must be called before re-running update/finalize.");
    }

    // Capture raw counts before normalization for statistical tests
    double base_tc_total = profile.baseline_t_freq + profile.baseline_c_freq;
    double base_ag_total = profile.baseline_a_freq + profile.baseline_g_freq;

    auto compute_terminal_stats = [](double term_t, double term_c, double base_t, double base_c) {
        double term_total = term_t + term_c;
        double base_total = base_t + base_c;
        if (term_total < 1.0 || base_total < 1.0) {
            return std::pair<float, float>{0.0f, 0.0f};
        }
        double p_term = term_t / term_total;
        double p_base = base_t / base_total;
        double p_pool = (term_t + base_t) / (term_total + base_total);
        double se = std::sqrt(std::max(1e-9, p_pool * (1.0 - p_pool) * (1.0 / term_total + 1.0 / base_total)));
        double z = se > 0.0 ? (p_term - p_base) / se : 0.0;
        return std::pair<float, float>{static_cast<float>(p_term - p_base), static_cast<float>(z)};
    };

    // Terminal enrichment/depletion statistics
    // Check for position-0 artifact: if pos0 shows depletion but pos1 shows enrichment,
    // use pos1 instead (common adapter ligation bias pattern)
    auto stats_5_pos0 = compute_terminal_stats(profile.t_freq_5prime[0], profile.c_freq_5prime[0],
                                               profile.baseline_t_freq, profile.baseline_c_freq);
    auto stats_5_pos1 = compute_terminal_stats(profile.t_freq_5prime[1], profile.c_freq_5prime[1],
                                               profile.baseline_t_freq, profile.baseline_c_freq);

    // Detect position-0 artifact: two criteria
    // 1. Classic: pos0 depleted but pos1 enriched (relative to baseline)
    // 2. Jump pattern: Large pos0-to-pos1 difference (>3%) with pos1 enriched (>2%)
    //    This catches cases where baseline elevation masks the pos0 depletion
    bool pos0_artifact_5 = false;
    if (stats_5_pos0.first < -0.005f && stats_5_pos1.first > 0.01f) {
        // Classic pattern: pos0 below baseline, pos1 above
        pos0_artifact_5 = true;
    } else if (stats_5_pos1.first > 0.02f && (stats_5_pos1.first - stats_5_pos0.first) > 0.03f) {
        // Large pos0-to-pos1 jump pattern: pos1 clearly elevated but pos0 isn't
        // This indicates adapter artifact masking damage at pos0
        pos0_artifact_5 = true;
    }

    if (pos0_artifact_5) {
        profile.terminal_shift_5prime = stats_5_pos1.first;
        profile.terminal_z_5prime = stats_5_pos1.second;
        profile.position_0_artifact_5prime = true;
    } else {
        profile.terminal_shift_5prime = stats_5_pos0.first;
        profile.terminal_z_5prime = stats_5_pos0.second;
        profile.position_0_artifact_5prime = false;
    }

    // Always use pooled positions 1-4 for the 3' terminal G→A signal.
    // Position 0 is excluded: SS library prep introduces an adapter ligation artifact
    // at 3' pos-0 (elevated A/(A+G)) that is unrelated to aDNA damage. At high coverage
    // this artifact drives terminal_z_3prime into the hundreds and falsely blocks SS
    // auto-detection. Report pos-0 separately as position_0_artifact_3prime.
    double sum_a_3_1_4 = 0.0, sum_g_3_1_4 = 0.0;
    for (int p = 1; p <= 4; ++p) {
        sum_a_3_1_4 += profile.a_freq_3prime[p];
        sum_g_3_1_4 += profile.g_freq_3prime[p];
    }
    auto stats_3_pos1_4 = compute_terminal_stats(sum_a_3_1_4, sum_g_3_1_4,
                                                 4.0 * profile.baseline_a_freq,
                                                 4.0 * profile.baseline_g_freq);
    profile.terminal_shift_3prime = stats_3_pos1_4.first;
    profile.terminal_z_3prime     = stats_3_pos1_4.second;

    // Flag pos-0 as an artifact when it is substantially elevated above the pos1-4 mean
    {
        double tot0   = profile.a_freq_3prime[0] + profile.g_freq_3prime[0];
        double tot1_4 = sum_a_3_1_4 + sum_g_3_1_4;
        double ga0    = tot0   > 0.0 ? profile.a_freq_3prime[0] / tot0   : 0.0;
        double ga1_4  = tot1_4 > 0.0 ? sum_a_3_1_4              / tot1_4 : 0.0;
        profile.position_0_artifact_3prime = (ga0 - ga1_4) > 0.05;
    }

    // Negative control statistics (should NOT show enrichment if damage is real)
    // 5' control: A/(A+G) at 5' end vs interior
    auto ctrl_5 = compute_terminal_stats(profile.a_freq_5prime[0], profile.g_freq_5prime[0],
                                         profile.baseline_a_freq, profile.baseline_g_freq);
    profile.ctrl_shift_5prime = ctrl_5.first;
    profile.ctrl_z_5prime = ctrl_5.second;

    // 3' control: T/(T+C) at 3' end vs interior
    auto ctrl_3 = compute_terminal_stats(profile.t_freq_3prime[0], profile.c_freq_3prime[0],
                                         profile.baseline_t_freq, profile.baseline_c_freq);
    profile.ctrl_shift_3prime = ctrl_3.first;
    profile.ctrl_z_3prime = ctrl_3.second;

    // Channel divergence: difference between damage shift and control shift
    // High divergence = real damage (only damage channel elevated)
    // Low divergence = composition bias (both channels elevated together)
    profile.channel_divergence_5prime = std::abs(profile.terminal_shift_5prime - profile.ctrl_shift_5prime);
    profile.channel_divergence_3prime = std::abs(profile.terminal_shift_3prime - profile.ctrl_shift_3prime);

    // Normalize baseline frequencies (using double to preserve precision)
    double mid_total = profile.baseline_t_freq + profile.baseline_c_freq +
                       profile.baseline_a_freq + profile.baseline_g_freq;
    if (mid_total > 0) {
        profile.baseline_t_freq /= mid_total;
        profile.baseline_c_freq /= mid_total;
        profile.baseline_a_freq /= mid_total;
        profile.baseline_g_freq /= mid_total;
    }

    // Compute baseline T/C and A/G ratios from the middle of reads (for inverted pattern detection)
    double baseline_tc = profile.baseline_t_freq /
                        (profile.baseline_t_freq + profile.baseline_c_freq + 0.001);
    double baseline_ag = profile.baseline_a_freq /
                        (profile.baseline_a_freq + profile.baseline_g_freq + 0.001);

    {
        JointDamageSuffStats jstats;

        // Channel A: T/(T+C) at 5' end (raw counts still in t_freq_5prime, c_freq_5prime)
        for (int p = 0; p < 15; ++p) {
            jstats.k_tc[p] = static_cast<uint64_t>(profile.t_freq_5prime[p]);
            jstats.n_tc[p] = static_cast<uint64_t>(profile.t_freq_5prime[p] + profile.c_freq_5prime[p]);
        }

        // Control: A/(A+G) at 5' end
        for (int p = 0; p < 15; ++p) {
            jstats.k_ag[p] = static_cast<uint64_t>(profile.a_freq_5prime[p]);
            jstats.n_ag[p] = static_cast<uint64_t>(profile.a_freq_5prime[p] + profile.g_freq_5prime[p]);
        }

        // Channel B: stop codon conversions (TAA+TAG+TGA) / (CAA+CAG+CGA + TAA+TAG+TGA)
        for (int p = 0; p < 15; ++p) {
            double stops = profile.convertible_taa_5prime[p] +
                          profile.convertible_tag_5prime[p] +
                          profile.convertible_tga_5prime[p];
            double pre = profile.convertible_caa_5prime[p] +
                        profile.convertible_cag_5prime[p] +
                        profile.convertible_cga_5prime[p];
            jstats.k_stop[p] = static_cast<uint64_t>(stops);
            jstats.n_stop[p] = static_cast<uint64_t>(pre + stops);
        }

        // Interior baselines (from middle of reads)
        jstats.k_tc_interior = static_cast<uint64_t>(base_tc_total * baseline_tc);
        jstats.n_tc_interior = static_cast<uint64_t>(base_tc_total);
        jstats.k_ag_interior = static_cast<uint64_t>(base_ag_total * baseline_ag);
        jstats.n_ag_interior = static_cast<uint64_t>(base_ag_total);

        // Channel B interior
        double stop_int = profile.convertible_taa_interior +
                         profile.convertible_tag_interior +
                         profile.convertible_tga_interior;
        double pre_int = profile.convertible_caa_interior +
                        profile.convertible_cag_interior +
                        profile.convertible_cga_interior;
        jstats.k_stop_interior = static_cast<uint64_t>(stop_int);
        jstats.n_stop_interior = static_cast<uint64_t>(pre_int + stop_int);

        // Fit the joint model
        JointDamageResult jresult = JointDamageModel::fit(jstats);

        // Store results in profile
        profile.joint_delta_max = jresult.delta_max;
        profile.joint_lambda = jresult.lambda;
        profile.joint_a_max = jresult.a_max;
        profile.joint_log_lik_m1 = jresult.log_lik_m1;
        profile.joint_log_lik_m0 = jresult.log_lik_m0;
        profile.joint_delta_bic = jresult.delta_bic;
        profile.joint_bayes_factor = jresult.bayes_factor;
        profile.joint_p_damage = jresult.p_damage;
        profile.joint_model_valid = jresult.valid;

        // Use joint model for damage decision
        if (jresult.valid) {
            // Check for strong evidence of damage
            // P(damage) > 0.95 OR ΔBIC > 10 (handles overflow when BF=inf)
            bool strong_damage_evidence = jresult.p_damage > 0.95f ||
                                          jresult.delta_bic > 10.0f ||
                                          std::isinf(jresult.bayes_factor);

            if (strong_damage_evidence) {
                profile.damage_validated = true;
                profile.damage_artifact = false;
            } else if (jresult.p_damage < 0.05f && jresult.a_max > 0.02f) {
                // Strong evidence for artifact (composition bias detected)
                profile.damage_validated = false;
                profile.damage_artifact = true;
            } else if (jresult.delta_bic < -10.0f) {
                // Strong evidence against damage (M0 fits much better)
                profile.damage_validated = false;
                profile.damage_artifact = (jresult.a_max > 0.02f);
            } else {
                // Uncertain - could go either way
                profile.damage_validated = false;
                profile.damage_artifact = false;
            }
        }
    }

    // Step 1: Normalize position-specific frequencies (but don't compute damage rates yet)
    for (int i = 0; i < 15; ++i) {
        double tc_total = profile.t_freq_5prime[i] + profile.c_freq_5prime[i];
        if (tc_total > 0) {
            profile.t_freq_5prime[i] = profile.t_freq_5prime[i] / tc_total;
            profile.c_freq_5prime[i] = 1.0 - profile.t_freq_5prime[i];
        }
        double ag_total = profile.a_freq_3prime[i] + profile.g_freq_3prime[i];
        if (ag_total > 0) {
            profile.a_freq_3prime[i] = profile.a_freq_3prime[i] / ag_total;
            profile.g_freq_3prime[i] = 1.0 - profile.a_freq_3prime[i];
        }
    }

    // Step 2: Fit p(pos) = b + A * exp(-lambda * pos) using middle-of-read baseline
    auto fit_5p = fit_exponential_decay(profile.t_freq_5prime, profile.tc_total_5prime, 0.2f,
                                        static_cast<float>(baseline_tc));
    profile.fit_baseline_5prime = fit_5p[0];
    profile.fit_amplitude_5prime = fit_5p[1];
    float fit_lambda_5p = fit_5p[2];
    profile.fit_rmse_5prime = fit_5p[3];

    auto fit_3p = fit_exponential_decay(profile.a_freq_3prime, profile.ag_total_3prime, 0.2f,
                                        static_cast<float>(baseline_ag));
    profile.fit_baseline_3prime = fit_3p[0];
    profile.fit_amplitude_3prime = fit_3p[1];
    float fit_lambda_3p = fit_3p[2];
    profile.fit_rmse_3prime = fit_3p[3];

    // Positive LLR = exponential fits better than constant model (real decay pattern)
    profile.decay_llr_5prime = compute_decay_llr(
        profile.t_freq_5prime, profile.tc_total_5prime,
        profile.fit_baseline_5prime, profile.fit_amplitude_5prime, fit_lambda_5p);
    profile.decay_llr_3prime = compute_decay_llr(
        profile.a_freq_3prime, profile.ag_total_3prime,
        profile.fit_baseline_3prime, profile.fit_amplitude_3prime, fit_lambda_3p);

    // Control channel decay LLR: A/(A+G) at 5', T/(T+C) at 3'.
    // If control also shows decay, it's likely composition/trimming artifact, not damage.
    // ctrl_freq_3p / ctrl_total_3p are hoisted out of the block so they remain accessible
    // to the library-type BIC classifier further below.
    std::array<double, 15> ctrl_freq_3p  = {};
    std::array<double, 15> ctrl_total_3p = {};
    std::array<double, 15> ctrl_freq_5p  = {};
    std::array<double, 15> ctrl_total_5p = {};
    {

        for (int i = 0; i < 15; ++i) {
            double ag_5p = profile.a_freq_5prime[i] + profile.g_freq_5prime[i];
            ctrl_total_5p[i] = ag_5p;
            ctrl_freq_5p[i] = (ag_5p > 0) ? profile.a_freq_5prime[i] / ag_5p : 0.5;

            double tc_3p = profile.t_freq_3prime[i] + profile.c_freq_3prime[i];
            ctrl_total_3p[i] = tc_3p;
            ctrl_freq_3p[i] = (tc_3p > 0) ? profile.t_freq_3prime[i] / tc_3p : 0.5;
            profile.tc_total_3prime[i] = tc_3p;
        }

        double ctrl_baseline_5p = baseline_ag;
        double ctrl_baseline_3p = baseline_tc;

        profile.ctrl_decay_llr_5prime = compute_decay_llr(
            ctrl_freq_5p, ctrl_total_5p,
            static_cast<float>(ctrl_baseline_5p), 0.0f, fit_lambda_5p);
        profile.ctrl_decay_llr_3prime = compute_decay_llr(
            ctrl_freq_3p, ctrl_total_3p,
            static_cast<float>(ctrl_baseline_3p), 0.0f, fit_lambda_3p);

        profile.delta_llr_5prime = profile.decay_llr_5prime - profile.ctrl_decay_llr_5prime;
        profile.delta_llr_3prime = profile.decay_llr_3prime - profile.ctrl_decay_llr_3prime;
    }

    // Channel B: stop codon conversion decay. Real C→T damage must create stops
    // in CAA, CAG, CGA contexts; unlike Channel A, this cannot arise from composition bias.
    {
        double total_pre_interior = profile.convertible_caa_interior +
                                   profile.convertible_cag_interior +
                                   profile.convertible_cga_interior;
        double total_stop_interior = profile.convertible_taa_interior +
                                    profile.convertible_tag_interior +
                                    profile.convertible_tga_interior;
        double total_convertible_interior = total_pre_interior + total_stop_interior;

        if (total_convertible_interior > 100) {
            profile.stop_conversion_rate_baseline = static_cast<float>(
                total_stop_interior / total_convertible_interior);
            profile.channel_b_valid = true;
        }

        std::array<double, 15> stop_rate = {};
        std::array<double, 15> stop_exposure = {};

        for (int p = 0; p < 15; ++p) {
            double pre = profile.convertible_caa_5prime[p] +
                        profile.convertible_cag_5prime[p] +
                        profile.convertible_cga_5prime[p];
            double stop = profile.convertible_taa_5prime[p] +
                         profile.convertible_tag_5prime[p] +
                         profile.convertible_tga_5prime[p];
            stop_exposure[p] = pre + stop;
            if (stop_exposure[p] > 10) {
                stop_rate[p] = stop / stop_exposure[p];
            } else {
                stop_rate[p] = profile.stop_conversion_rate_baseline;  // Use baseline if no data
            }
        }

        // Local baseline from positions 5-14 (same reads, past damage zone)
        double local_pre = 0.0, local_stop = 0.0;
        for (int p = 5; p < 15; ++p) {
            local_pre += profile.convertible_caa_5prime[p] +
                        profile.convertible_cag_5prime[p] +
                        profile.convertible_cga_5prime[p];
            local_stop += profile.convertible_taa_5prime[p] +
                         profile.convertible_tag_5prime[p] +
                         profile.convertible_tga_5prime[p];
        }
        float local_baseline = (local_pre + local_stop > 100)
            ? static_cast<float>(local_stop / (local_pre + local_stop))
            : profile.stop_conversion_rate_baseline;

        if (profile.channel_b_valid) {
            float baseline_b = local_baseline;
            float lambda_b = fit_lambda_5p;

            // Amplitude from positions 0-4: include pos 0 since steep decay (AT-rich)
            // concentrates signal there and excluding it causes false negatives
            double sum_excess = 0.0, sum_weight = 0.0;
            for (int i = 0; i < 5; ++i) {
                if (stop_exposure[i] > 50) {
                    double weight = std::exp(-lambda_b * i);
                    double excess = stop_rate[i] - baseline_b;
                    sum_excess += stop_exposure[i] * excess / weight;
                    sum_weight += stop_exposure[i];
                }
            }
            float amplitude_b = (sum_weight > 0) ? static_cast<float>(sum_excess / sum_weight) : 0.0f;
            profile.stop_amplitude_5prime = std::max(0.0f, amplitude_b);

            // LLR for exponential vs constant (include pos 0 for same reason as above)
            double ll_exp = 0.0, ll_const = 0.0;
            for (int p = 0; p < 10; ++p) {
                if (stop_exposure[p] < 50) continue;

                double n = stop_exposure[p];
                double pre = profile.convertible_caa_5prime[p] +
                            profile.convertible_cag_5prime[p] +
                            profile.convertible_cga_5prime[p];
                double stop = profile.convertible_taa_5prime[p] +
                             profile.convertible_tag_5prime[p] +
                             profile.convertible_tga_5prime[p];
                double k = stop;

                double p_exp = baseline_b + amplitude_b * std::exp(-lambda_b * p);
                p_exp = std::clamp(p_exp, 0.001, 0.999);
                ll_exp += binomial_ll(k, n, p_exp);

                double p_const = std::clamp(static_cast<double>(baseline_b), 0.001, 0.999);
                ll_const += binomial_ll(k, n, p_const);
            }

            profile.stop_decay_llr_5prime = static_cast<float>(ll_exp - ll_const);

            // If amplitude is negative (inverted), negate LLR
            if (amplitude_b < 0) {
                profile.stop_decay_llr_5prime = -std::abs(profile.stop_decay_llr_5prime);
            }

            {
                // WLS fit: r_p = b0 + (1-b0) * d_max * exp(-λp)
                // Solve y_p = a + c*x_p, then b0 = a, d_max = c / (1 - b0)
                const double lambda = std::clamp(static_cast<double>(fit_lambda_5p), 0.1, 0.5);
                const int N_POSITIONS = 15;

                double S_w = 0, S_x = 0, S_xx = 0, S_y = 0, S_xy = 0;
                double total_exposure = 0;

                for (int p = 0; p < N_POSITIONS; ++p) {
                    double x_p = std::exp(-lambda * p);

                    double stops_p = profile.convertible_taa_5prime[p] +
                                    profile.convertible_tag_5prime[p] +
                                    profile.convertible_tga_5prime[p];
                    double pre_p = profile.convertible_caa_5prime[p] +
                                  profile.convertible_cag_5prime[p] +
                                  profile.convertible_cga_5prime[p];
                    double exposure_p = stops_p + pre_p;

                    if (exposure_p < 100) continue;  // Skip low-coverage positions

                    double y_p = stops_p / exposure_p;
                    double w_p = exposure_p;

                    S_w  += w_p;
                    S_x  += w_p * x_p;
                    S_xx += w_p * x_p * x_p;
                    S_y  += w_p * y_p;
                    S_xy += w_p * x_p * y_p;
                    total_exposure += exposure_p;
                }

                double denom = S_w * S_xx - S_x * S_x;

                if (total_exposure > 1000 && std::abs(denom) > 1e-10) {
                    double c = (S_w * S_xy - S_x * S_y) / denom;
                    double a = (S_y - c * S_x) / S_w;

                    profile.channel_b_slope = static_cast<float>(c);

                    if (c > 0) {
                        double b0 = std::clamp(a, 0.01, 0.99);
                        double d_max_b = std::clamp(c / (1.0 - b0), 0.0, 1.0);

                        profile.d_max_from_channel_b = static_cast<float>(d_max_b);
                        profile.channel_b_weight = static_cast<float>(S_w);
                        profile.channel_b_quantifiable = true;
                        profile.channel_b_inverted = false;
                    } else {
                        // Inverted pattern: terminal stops LOWER than baseline
                        profile.d_max_from_channel_b = 0.0f;
                        profile.channel_b_weight = 0.0f;
                        profile.channel_b_quantifiable = false;
                        profile.channel_b_inverted = true;
                    }
                } else {
                    profile.d_max_from_channel_b = 0.0f;
                    profile.channel_b_weight = 0.0f;
                    profile.channel_b_quantifiable = false;
                    profile.channel_b_slope = 0.0f;
                }
            }

            (void)profile.delta_llr_5prime;
            (void)profile.stop_decay_llr_5prime;
        }
    }

    // Channel C: oxidative stop codon analysis.
    // Unlike deamination, oxidation is uniform across reads (terminal ≈ interior rate).
    {
        double ox_pre_interior = profile.convertible_gag_interior +
                                 profile.convertible_gaa_interior +
                                 profile.convertible_gga_interior;
        double ox_stop_interior = profile.convertible_tag_ox_interior +
                                  profile.convertible_taa_ox_interior +
                                  profile.convertible_tga_ox_interior;
        double ox_total_interior = ox_pre_interior + ox_stop_interior;

        if (ox_total_interior > 100) {
            profile.ox_stop_conversion_rate_baseline = static_cast<float>(
                ox_stop_interior / ox_total_interior);
            profile.channel_c_valid = true;
        }

        double ox_pre_terminal = 0.0, ox_stop_terminal = 0.0;
        double ox_pre_mid = 0.0, ox_stop_mid = 0.0;

        for (int p = 0; p < 5; ++p) {
            ox_pre_terminal += profile.convertible_gag_5prime[p] +
                              profile.convertible_gaa_5prime[p] +
                              profile.convertible_gga_5prime[p];
            ox_stop_terminal += profile.convertible_tag_ox_5prime[p] +
                               profile.convertible_taa_ox_5prime[p] +
                               profile.convertible_tga_ox_5prime[p];
        }
        for (int p = 5; p < 15; ++p) {
            ox_pre_mid += profile.convertible_gag_5prime[p] +
                         profile.convertible_gaa_5prime[p] +
                         profile.convertible_gga_5prime[p];
            ox_stop_mid += profile.convertible_tag_ox_5prime[p] +
                          profile.convertible_taa_ox_5prime[p] +
                          profile.convertible_tga_ox_5prime[p];
        }

        double ox_total_terminal = ox_pre_terminal + ox_stop_terminal;
        double ox_total_mid = ox_pre_mid + ox_stop_mid;

        profile.ox_uniformity_ratio = 1.0f;
        if (ox_total_terminal > 50 && ox_total_mid > 50) {
            profile.ox_stop_rate_terminal = static_cast<float>(ox_stop_terminal / ox_total_terminal);
            profile.ox_stop_rate_interior = static_cast<float>(ox_stop_mid / ox_total_mid);

            // Uniformity ratio: terminal/interior (≈1 for uniform oxidation)
            if (profile.ox_stop_rate_interior > 0.001f) {
                profile.ox_uniformity_ratio = profile.ox_stop_rate_terminal / profile.ox_stop_rate_interior;
            }
        }

        float ox_stop_excess = profile.ox_stop_rate_terminal - profile.ox_stop_conversion_rate_baseline;

        // Use stricter threshold if deamination is present (correlated damage)
        // Deamination enriches ALL terminal damage, including G→T via adjacent effects
        float threshold = 0.02f;  // 2% excess required (was 0.5%)
        if (profile.d_max_combined > 0.05f || profile.damage_validated) {
            threshold = 0.05f;  // 5% excess if deamination present
        }

        bool elevated = ox_stop_excess > threshold;
        bool uniform = profile.ox_uniformity_ratio > 0.85f && profile.ox_uniformity_ratio < 1.15f;

        if (profile.channel_c_valid && elevated && uniform) {
            profile.ox_damage_detected = true;
            profile.ox_d_max = std::max(0.0f, ox_stop_excess);  // fraction [0,1], consistent with other d_max fields
        }

        // Terminal much higher than interior suggests deamination cross-contamination, not true G→T
        if (profile.ox_uniformity_ratio > 1.5f) {
            profile.ox_is_artifact = true;
            profile.ox_damage_detected = false;  // Override - not true oxidation
        }

        if (profile.ox_uniformity_ratio < 0.7f) {
            profile.ox_is_artifact = true;
            profile.ox_damage_detected = false;
        }
    }

    // Channel D: Chargaff-balance G→T / C→A oxidation estimate.
    // Reference-free: under no damage, T/(T+G) ≈ A/(A+C) (Chargaff) and terminal ≈ interior.
    // 8-oxoG converts G→T uniformly, raising T/(T+G) above A/(A+C) throughout reads.
    {
        // Interior baseline T/(T+G) and A/(A+C)
        double gt_base_total = profile.baseline_g_total + profile.baseline_g_to_t_count;
        double ca_base_total = profile.baseline_c_ox_total + profile.baseline_c_to_a_count;
        if (gt_base_total >= 500.0) {
            profile.ox_gt_baseline = static_cast<float>(profile.baseline_g_to_t_count / gt_base_total);
        }
        if (ca_base_total >= 500.0) {
            profile.ox_ca_baseline = static_cast<float>(profile.baseline_c_to_a_count / ca_base_total);
        }

        // T/(T+G) at positions 0-4 (terminal) and 5-14 (mid-read)
        double t_term = 0, g_term = 0, t_mid = 0, g_mid = 0;
        for (int p = 0; p < 5; ++p)  { t_term += profile.t_from_g_5prime[p]; g_term += profile.g_count_5prime[p]; }
        for (int p = 5; p < 15; ++p) { t_mid  += profile.t_from_g_5prime[p]; g_mid  += profile.g_count_5prime[p]; }
        if (t_term + g_term >= 200.0) profile.ox_gt_rate_terminal = static_cast<float>(t_term / (t_term + g_term));
        if (t_mid  + g_mid  >= 200.0) {
            profile.ox_gt_rate_interior = static_cast<float>(t_mid / (t_mid + g_mid));
            if (profile.ox_gt_rate_interior > 0.001f && profile.ox_gt_rate_terminal > 0.0f)
                profile.ox_gt_uniformity = profile.ox_gt_rate_terminal / profile.ox_gt_rate_interior;
        }

        // A/(A+C) at positions 0-4 (terminal) and 5-14 (mid-read)
        double a_term = 0, c_term = 0, a_mid = 0, c_mid = 0;
        for (int p = 0; p < 5; ++p)  { a_term += profile.a_from_c_5prime[p]; c_term += profile.c_count_ox_5prime[p]; }
        for (int p = 5; p < 15; ++p) { a_mid  += profile.a_from_c_5prime[p]; c_mid  += profile.c_count_ox_5prime[p]; }
        if (a_term + c_term >= 200.0) profile.ox_ca_rate_terminal = static_cast<float>(a_term / (a_term + c_term));
        if (a_mid  + c_mid  >= 200.0) {
            profile.ox_ca_rate_interior = static_cast<float>(a_mid / (a_mid + c_mid));
            if (profile.ox_ca_rate_interior > 0.001f && profile.ox_ca_rate_terminal > 0.0f)
                profile.ox_ca_uniformity = profile.ox_ca_rate_terminal / profile.ox_ca_rate_interior;
        }

        // Chargaff asymmetry: excess T/(T+G) over A/(A+C) at interior (= Chargaff deviation)
        if (gt_base_total >= 500.0 && ca_base_total >= 500.0)
            profile.ox_gt_asymmetry = profile.ox_gt_baseline - profile.ox_ca_baseline;
    }

    // GT exponential-background fit: GT(p) = A*exp(-mu*p) + B
    // B is the uniform background (8-oxoG rate); A is the terminal artifact.
    // Mirrors the C→T deamination model, extracting the position-independent component.
    // For SS: s_gt = B - ox_ca_baseline is a valid 8-oxoG signal (no Chargaff cancellation).
    // For DS: B and A/(A+C) both increase with 8-oxoG, so s_gt ≈ 0; report B directly.
    {
        float y[15] = {}, w[15] = {};
        float total_w = 0;
        for (int p = 0; p < 15; ++p) {
            double tot = profile.t_from_g_5prime[p] + profile.g_count_5prime[p];
            if (tot > 10.0) {
                w[p] = (float)tot;
                y[p] = (float)(profile.t_from_g_5prime[p] / tot);
            }
            total_w += w[p];
        }

        if (total_w >= 200.0f) {
            static constexpr float mu_grid[] = {0.05f, 0.1f, 0.2f, 0.3f, 0.5f, 0.7f, 1.0f, 1.5f, 2.0f, 3.0f};
            float best_sse = 1e30f, best_A = 0.0f, best_mu = 0.3f, best_B = 0.0f;
            float best_A_raw = 0.0f;

            for (float mu : mu_grid) {
                float sx2=0, sx=0, sxy=0, sy=0, sw=0;
                for (int p = 0; p < 15; ++p) {
                    if (w[p] == 0) continue;
                    float x = std::exp(-mu * p);
                    sx2 += w[p] * x * x;
                    sx  += w[p] * x;
                    sxy += w[p] * x * y[p];
                    sy  += w[p] * y[p];
                    sw  += w[p];
                }
                float det = sx2*sw - sx*sx;
                if (std::abs(det) < 1e-12f) continue;
                float A_raw = (sxy*sw - sy*sx) / det;
                float B_raw = (sx2*sy - sx*sxy) / det;
                float A = std::max(0.0f, A_raw);
                float B = std::max(0.0f, std::min(0.5f, B_raw));

                float sse = 0;
                for (int p = 0; p < 15; ++p) {
                    if (w[p] == 0) continue;
                    float res = y[p] - A * std::exp(-mu * p) - B;
                    sse += w[p] * res * res;
                }
                if (sse < best_sse) {
                    best_sse = sse;
                    best_A = A; best_mu = mu; best_B = B;
                    best_A_raw = A_raw;
                }
            }

            profile.g_bg_fitted   = best_B;
            profile.g_term_fitted  = best_A;
            profile.g_decay_fitted = best_mu;
            if (profile.ox_ca_baseline > 0.0f)
                profile.s_gt = best_B - profile.ox_ca_baseline;

            // Update detection: model-based uniform G→T signal replaces codon-based Channel C.
            // For SS: s_gt > threshold (Chargaff contrast valid because no complementary strand).
            // For DS: B elevated above ca_baseline indicates library-level oxidation.
            const bool is_ss_lib = (profile.library_type ==
                                    taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED);
            // SS: Chargaff contrast (s_gt = B - ca_baseline) is the valid signal.
            // DS: B and ca_baseline rise together so s_gt cancels; use best_B directly.
            float signal    = is_ss_lib ? profile.s_gt : best_B;
            float threshold = is_ss_lib ? 0.004f : 0.006f;

            // Require the terminal component to also be non-trivial (rules out flat noise)
            bool has_data   = total_w >= 500.0f;
            bool elevated   = signal > threshold;
            bool not_inverted = best_A_raw >= 0.0f;

            profile.ox_damage_detected = has_data && elevated && not_inverted;
        }
    }

    // Channel E: depurination detection (terminal purine enrichment)
    {
        double purine_terminal = 0.0, total_terminal = 0.0;
        for (int p = 0; p < 5; ++p) {
            // a_freq_5prime/g_freq_5prime are normalized; reconstruct purine count via tc_total
            purine_terminal += profile.tc_total_5prime[p] * (1.0 - profile.t_freq_5prime[p] - profile.c_freq_5prime[p]);
            total_terminal += profile.tc_total_5prime[p];
        }

        double purine_baseline = profile.baseline_a_freq + profile.baseline_g_freq;

        if (total_terminal > 100 && purine_baseline > 0.01) {
            profile.purine_rate_interior = static_cast<float>(purine_baseline);

            profile.purine_enrichment_5prime = profile.ctrl_shift_5prime;
            profile.purine_enrichment_3prime = profile.ctrl_shift_3prime;

            if (profile.purine_enrichment_5prime > 0.02f &&
                profile.channel_divergence_5prime > 0.01f) {
                profile.depurination_detected = true;
            }
        }
    }

    // Update lambda from fit if reasonable
    if (fit_lambda_5p > 0.05f && fit_lambda_5p < 0.5f) {
        profile.lambda_5prime = fit_lambda_5p;
    }
    if (fit_lambda_3p > 0.05f && fit_lambda_3p < 0.5f) {
        profile.lambda_3prime = fit_lambda_3p;
    }

    // Context-split amplitude fits (CpG-like vs non-CpG-like)
    {
        auto fit_cpg = fit_ct5_ctx_amplitude(
            profile.ct_ctx_t_5prime[SampleDamageProfile::CPG_LIKE],
            profile.ct_ctx_total_5prime[SampleDamageProfile::CPG_LIKE],
            profile.ct_ctx_t_interior[SampleDamageProfile::CPG_LIKE],
            profile.ct_ctx_total_interior[SampleDamageProfile::CPG_LIKE],
            profile.lambda_5prime);

        auto fit_ncpg = fit_ct5_ctx_amplitude(
            profile.ct_ctx_t_5prime[SampleDamageProfile::NONCPG_LIKE],
            profile.ct_ctx_total_5prime[SampleDamageProfile::NONCPG_LIKE],
            profile.ct_ctx_t_interior[SampleDamageProfile::NONCPG_LIKE],
            profile.ct_ctx_total_interior[SampleDamageProfile::NONCPG_LIKE],
            profile.lambda_5prime);

        profile.fit_baseline_ct5_cpg_like    = fit_cpg.baseline;
        profile.fit_baseline_ct5_noncpg_like = fit_ncpg.baseline;
        profile.dmax_ct5_cpg_like    = fit_cpg.valid  ? fit_cpg.dmax  : std::numeric_limits<float>::quiet_NaN();
        profile.dmax_ct5_noncpg_like = fit_ncpg.valid ? fit_ncpg.dmax : std::numeric_limits<float>::quiet_NaN();
        profile.cov_ct5_cpg_like_terminal       = fit_cpg.cov_terminal;
        profile.cov_ct5_noncpg_like_terminal    = fit_ncpg.cov_terminal;
        profile.cov_ct5_cpg_like_interior       = fit_cpg.cov_interior;
        profile.cov_ct5_noncpg_like_interior    = fit_ncpg.cov_interior;
        profile.effcov_ct5_cpg_like_terminal    = fit_cpg.effcov_terminal;
        profile.effcov_ct5_noncpg_like_terminal = fit_ncpg.effcov_terminal;
        profile.effcov_ct5_cpg_like_interior    = fit_cpg.effcov_interior;
        profile.effcov_ct5_noncpg_like_interior = fit_ncpg.effcov_interior;
        profile.fit_positions_ct5_cpg_like    = fit_cpg.fit_positions;
        profile.fit_positions_ct5_noncpg_like = fit_ncpg.fit_positions;

        if (fit_cpg.valid && fit_ncpg.valid && fit_ncpg.dmax >= 0.005f) {
            profile.cpg_ratio    = fit_cpg.dmax / fit_ncpg.dmax;
            profile.log2_cpg_ratio = std::log2((fit_cpg.dmax + 1e-6f) / (fit_ncpg.dmax + 1e-6f));
        }

        // oxoG 16-context finalization
        for (int i = 0; i < SampleDamageProfile::N_OXOG16; ++i) {
            const float t = profile.oxog16_t[i], a = profile.oxog16_a_rc[i];
            const float cov = t + a;
            profile.cov_oxog_16ctx[i] = cov;
            profile.s_oxog_16ctx[i] = (cov >= 500.0f)
                ? (t - a) / cov
                : std::numeric_limits<float>::quiet_NaN();
        }

        // Upstream-context-aware C→T fitting (experimental: AC, CC, GC, TC)
        // Fit each context independently, then compute contrasts
        float dmax_sum = 0.0f;
        int valid_ctx_count = 0;
        for (int uctx = 0; uctx < SampleDamageProfile::N_UPSTREAM_CTX; ++uctx) {
            auto fit = fit_ct5_ctx_amplitude(
                profile.ct5_t_by_upstream[uctx],
                profile.ct5_total_by_upstream[uctx],
                profile.ct5_t_interior_by_upstream[uctx],
                profile.ct5_total_interior_by_upstream[uctx],
                profile.lambda_5prime);

            profile.baseline_ct5_by_upstream[uctx] = fit.baseline;
            profile.dmax_ct5_by_upstream[uctx] = fit.valid ? fit.dmax : std::numeric_limits<float>::quiet_NaN();
            profile.cov_ct5_terminal_by_upstream[uctx] = fit.cov_terminal;
            profile.cov_ct5_interior_by_upstream[uctx] = fit.cov_interior;

            if (fit.valid && fit.dmax >= 0.0f) {
                dmax_sum += fit.dmax;
                ++valid_ctx_count;
            }
        }

        // Compute contrasts if all 4 contexts are valid
        if (valid_ctx_count == 4) {
            const float d_ac = profile.dmax_ct5_by_upstream[SampleDamageProfile::CTX_AC];
            const float d_cc = profile.dmax_ct5_by_upstream[SampleDamageProfile::CTX_CC];
            const float d_gc = profile.dmax_ct5_by_upstream[SampleDamageProfile::CTX_GC];
            const float d_tc = profile.dmax_ct5_by_upstream[SampleDamageProfile::CTX_TC];

            // dipyr_contrast: mean(CC, TC) - mean(AC, GC)
            // Positive = dipyrimidine excess (UV-like)
            profile.dipyr_contrast = 0.5f * (d_cc + d_tc) - 0.5f * (d_ac + d_gc);

            // cpg_contrast: GC - mean(AC, CC, TC)
            // Positive = CpG excess (methylation-enhanced deamination)
            profile.cpg_contrast = d_gc - (d_ac + d_cc + d_tc) / 3.0f;

            // Chi-squared heterogeneity test: are contexts different from uniform?
            const float mean_d = dmax_sum / 4.0f;
            if (mean_d > 0.001f) {
                float chi2 = 0.0f;
                for (int uctx = 0; uctx < 4; ++uctx) {
                    const float d = profile.dmax_ct5_by_upstream[uctx];
                    const float cov = profile.cov_ct5_terminal_by_upstream[uctx];
                    if (cov > 100.0f) {
                        // Weight by coverage (pseudo chi-squared)
                        chi2 += cov * (d - mean_d) * (d - mean_d) / (mean_d * (1.0f - mean_d) + 1e-6f);
                    }
                }
                profile.context_heterogeneity_chi2 = chi2;
                // Approximate p-value (chi2 with 3 df)
                // Using simple threshold: chi2 > 7.81 → p < 0.05
                profile.context_heterogeneity_detected = (chi2 > 7.81f);
                // Store approximate p-value (very rough)
                if (chi2 < 0.58f) profile.context_heterogeneity_p = 0.9f;
                else if (chi2 < 2.37f) profile.context_heterogeneity_p = 0.5f;
                else if (chi2 < 6.25f) profile.context_heterogeneity_p = 0.1f;
                else if (chi2 < 7.81f) profile.context_heterogeneity_p = 0.05f;
                else if (chi2 < 11.34f) profile.context_heterogeneity_p = 0.01f;
                else profile.context_heterogeneity_p = 0.001f;
            }
        }
    }

    // Interior clustered C→T finalization
    {
        const auto& acc = profile.interior_ct_cluster;
        auto safe_log2_oe = [](uint64_t obs, double exp) -> float {
            return static_cast<float>(
                std::log((static_cast<double>(obs) + 0.5) / (exp + 0.5)) / std::log(2.0));
        };

        uint64_t obs_ct_s = 0, pairs_ct_s = 0, obs_ag_s = 0;
        double   exp_ct_s = 0.0, var_ct_s = 0.0;
        double   exp_ag_s = 0.0, var_ag_s = 0.0;

        for (int d = 1; d <= 10; ++d) {
            profile.interior_ct_cluster_sep_log2oe[d - 1] =
                safe_log2_oe(acc.obs_ct[d], acc.exp_ct[d]);
            if (d <= 5) {
                obs_ct_s  += acc.obs_ct[d];   pairs_ct_s += acc.pairs_ct[d];
                exp_ct_s  += acc.exp_ct[d];   var_ct_s   += acc.var_ct[d];
                obs_ag_s  += acc.obs_ag[d];
                exp_ag_s  += acc.exp_ag[d];   var_ag_s   += acc.var_ag[d];
            }
        }

        profile.interior_ct_cluster_short_obs          = obs_ct_s;
        profile.interior_ct_cluster_short_pairs        = pairs_ct_s;
        profile.interior_ct_cluster_short_exp          = exp_ct_s;
        profile.interior_ct_cluster_reads_used         = acc.reads_used_ct;
        profile.interior_ct_cluster_reads_used_control = acc.reads_used_ag;
        profile.interior_ct_cluster_reads_skipped      = acc.short_reads_skipped;

        profile.interior_ct_cluster_short_log2oe      = safe_log2_oe(obs_ct_s, exp_ct_s);
        const float ag_log2oe                         = safe_log2_oe(obs_ag_s, exp_ag_s);
        profile.interior_ct_cluster_short_asym_log2oe =
            profile.interior_ct_cluster_short_log2oe - ag_log2oe;

        const double denom = std::sqrt(var_ct_s + var_ag_s + 1e-12);
        const double num   = (static_cast<double>(obs_ct_s) - exp_ct_s)
                           - (static_cast<double>(obs_ag_s) - exp_ag_s);
        profile.interior_ct_cluster_short_z = static_cast<float>(num / denom);
    }

    // Step 3: Per-position damage rates using fit baseline
    float fit_baseline_c_frac_5p = 1.0f - profile.fit_baseline_5prime;
    float fit_baseline_g_frac_3p = 1.0f - profile.fit_baseline_3prime;

    for (int i = 0; i < 15; ++i) {
        if (fit_baseline_c_frac_5p > 0.1f) {
            float raw_signal = static_cast<float>(profile.t_freq_5prime[i]) - profile.fit_baseline_5prime;
            profile.damage_rate_5prime[i] = std::max(0.0f, raw_signal / fit_baseline_c_frac_5p);
        } else {
            profile.damage_rate_5prime[i] = 0.0f;
        }

        if (fit_baseline_g_frac_3p > 0.1f) {
            float raw_signal = static_cast<float>(profile.a_freq_3prime[i]) - profile.fit_baseline_3prime;
            profile.damage_rate_3prime[i] = std::max(0.0f, raw_signal / fit_baseline_g_frac_3p);
        } else {
            profile.damage_rate_3prime[i] = 0.0f;
        }
    }

    // Terminal gradients for reporting; inverted_pattern flags set later by hexamer detection
    {
        double terminal_tc_5 = profile.t_freq_5prime[0];
        profile.terminal_gradient_5prime = static_cast<float>(terminal_tc_5 - baseline_tc);

        double terminal_ag_3 = profile.a_freq_3prime[0];
        profile.terminal_gradient_3prime = static_cast<float>(terminal_ag_3 - baseline_ag);
    }

    for (int p = 0; p < 3; p++) {
        size_t tc_total = profile.codon_pos_t_count_5prime[p] + profile.codon_pos_c_count_5prime[p];
        if (tc_total > 0) {
            profile.codon_pos_t_rate_5prime[p] = static_cast<float>(profile.codon_pos_t_count_5prime[p]) / tc_total;
        }

        size_t ag_total = profile.codon_pos_a_count_3prime[p] + profile.codon_pos_g_count_3prime[p];
        if (ag_total > 0) {
            profile.codon_pos_a_rate_3prime[p] = static_cast<float>(profile.codon_pos_a_count_3prime[p]) / ag_total;
        }
    }

    {
        float pos1_t_rate = profile.codon_pos_t_rate_5prime[0];
        float pos2_t_rate = profile.codon_pos_t_rate_5prime[1];
        float pos3_t_rate = profile.codon_pos_t_rate_5prime[2];  // Wobble position

        float baseline_tc_f = profile.fit_baseline_5prime;
        float baseline_c_f = 1.0f - profile.fit_baseline_5prime;

        if (baseline_c_f > 0.1f) {
            profile.codon_pos1_damage = std::max(0.0f, (pos1_t_rate - baseline_tc_f) / baseline_c_f);
            profile.codon_pos2_damage = std::max(0.0f, (pos2_t_rate - baseline_tc_f) / baseline_c_f);
            profile.codon_pos3_damage = std::max(0.0f, (pos3_t_rate - baseline_tc_f) / baseline_c_f);
        } else {
            profile.codon_pos1_damage = 0.0f;
            profile.codon_pos2_damage = 0.0f;
            profile.codon_pos3_damage = 0.0f;
        }

        // Calculate wobble ratio: pos3 / ((pos1 + pos2) / 2)
        float avg_pos12_damage = (profile.codon_pos1_damage + profile.codon_pos2_damage) / 2.0f;
        if (avg_pos12_damage > 0.005f) {
            profile.wobble_ratio = profile.codon_pos3_damage / avg_pos12_damage;
        } else if (profile.codon_pos3_damage > 0.005f) {
            profile.wobble_ratio = 2.0f;  // Cap at 2x
        } else {
            profile.wobble_ratio = 1.0f;
        }

        profile.wobble_ratio = std::clamp(profile.wobble_ratio, 0.5f, 3.0f);
    }

    size_t cpg_total = profile.cpg_c_count + profile.cpg_t_count;
    if (cpg_total > 10) {
        profile.cpg_damage_rate = static_cast<float>(profile.cpg_t_count) / cpg_total;
    }

    size_t non_cpg_total = profile.non_cpg_c_count + profile.non_cpg_t_count;
    if (non_cpg_total > 10) {
        profile.non_cpg_damage_rate = static_cast<float>(profile.non_cpg_t_count) / non_cpg_total;
    }

    profile.max_damage_5prime = profile.damage_rate_5prime[0];
    profile.max_damage_3prime = profile.damage_rate_3prime[0];

    // Briggs-like damage model parameter estimation (closed-form).
    // start_pos: first position to use (1 when chemistry tag suppresses pos 0).
    // fixed_bg: when >= 0, use this as the background instead of the tail mean
    // computed from rates[10..14]. Tail-anchored bg gives a chemistry-robust
    // baseline that does not trade off with d_max in the fit.
    auto estimate_briggs_params = [](const std::array<float, 15>& rates, float max_rate,
                                     int start_pos, float fixed_bg,
                                     float& delta_s, float& delta_d, float& lambda, float& r_squared) {
        delta_s = 0.0f;
        delta_d = 0.0f;
        lambda = 0.3f;
        r_squared = 0.0f;

        if (max_rate < 0.02f) return;
        if (start_pos < 0 || start_pos > 5) start_pos = 0;

        delta_s = rates[start_pos];

        if (fixed_bg >= 0.0f) {
            delta_d = fixed_bg;
        } else {
            float sum_background = 0.0f;
            for (int i = 10; i < 15; ++i) sum_background += rates[i];
            delta_d = sum_background / 5.0f;
        }

        if (delta_s <= delta_d + 0.01f) {
            lambda = 0.3f;
            return;
        }

        // Weighted log-linear regression for lambda estimation.
        // Fit log((rate-bg)/(d_s-bg)) = -slope * (pos - start_pos) over
        // start_pos..14. x is the offset from start_pos so the intercept
        // corresponds to delta_s.
        float sum_xy = 0.0f, sum_x2 = 0.0f, sum_y = 0.0f, sum_x = 0.0f;
        float sum_w  = 0.0f, sum_y2 = 0.0f;

        for (int pos = start_pos; pos < 15; ++pos) {
            float normalized = (rates[pos] - delta_d) / (delta_s - delta_d + 0.001f);
            normalized = std::clamp(normalized, 0.01f, 1.0f);
            float y = std::log(normalized);
            float x = static_cast<float>(pos - start_pos);
            float w = normalized * normalized;
            sum_xy += w * x * y;
            sum_x2 += w * x * x;
            sum_y  += w * y;
            sum_x  += w * x;
            sum_w  += w;
            sum_y2 += w * y * y;
        }

        float denom = sum_w * sum_x2 - sum_x * sum_x;
        if (std::abs(denom) < 0.001f) { lambda = 0.3f; return; }

        float slope = (sum_w * sum_xy - sum_x * sum_y) / denom;
        if (slope >= 0.0f) { lambda = 0.3f; return; }

        lambda = 1.0f - std::exp(slope);
        lambda = std::clamp(lambda, 0.05f, 0.95f);

        float ss_tot = sum_y2 - (sum_y * sum_y) / sum_w;
        float intercept = (sum_y - slope * sum_x) / sum_w;
        float ss_res = 0.0f;
        for (int pos = start_pos; pos < 15; ++pos) {
            float normalized = (rates[pos] - delta_d) / (delta_s - delta_d + 0.001f);
            normalized = std::clamp(normalized, 0.01f, 1.0f);
            float y = std::log(normalized);
            float y_pred = slope * static_cast<float>(pos - start_pos) + intercept;
            float w = normalized * normalized;
            ss_res += w * (y - y_pred) * (y - y_pred);
        }
        r_squared = (ss_tot > 0.001f) ? std::max(0.0f, 1.0f - ss_res / ss_tot) : 0.0f;
    };

    // Tail-anchored background: trimmed mean of per-position C->T (G->A)
    // rates over BG_TAIL_LO..BG_TAIL_HI, requiring denom >= 100. If too few
    // eligible positions, falls through to the legacy joint-fit bg.
    auto compute_anchored_bg = [](const std::array<double, SampleDamageProfile::BG_TAIL_N>& num,
                                  const std::array<double, SampleDamageProfile::BG_TAIL_N>& den,
                                  int min_n_positions, double trim_frac,
                                  float& bg_out, int& n_pos_out, double& denom_out) {
        bg_out = -1.0f;
        n_pos_out = 0;
        denom_out = 0.0;
        std::vector<float> rates;
        rates.reserve(SampleDamageProfile::BG_TAIL_N);
        for (int i = 0; i < SampleDamageProfile::BG_TAIL_N; ++i) {
            if (den[i] >= 100.0) {
                rates.push_back(static_cast<float>(num[i] / den[i]));
                denom_out += den[i];
            }
        }
        n_pos_out = static_cast<int>(rates.size());
        if (n_pos_out < min_n_positions) return;
        std::sort(rates.begin(), rates.end());
        int n_trim = static_cast<int>(std::floor(trim_frac * n_pos_out));
        int lo = n_trim, hi = n_pos_out - n_trim;
        if (hi <= lo) { lo = 0; hi = n_pos_out; }
        double sum = 0.0;
        for (int i = lo; i < hi; ++i) sum += rates[i];
        bg_out = static_cast<float>(sum / static_cast<double>(hi - lo));
    };

    {
        float  bg5 = -1.0f, bg3 = -1.0f;
        int    n5  = 0,     n3  = 0;
        double d5  = 0.0,   d3  = 0.0;
        compute_anchored_bg(profile.tail_t_5prime,  profile.tail_tc_5prime,  10, 0.20, bg5, n5, d5);
        compute_anchored_bg(profile.tail_a_3prime,  profile.tail_ag_3prime,  10, 0.20, bg3, n3, d3);
        profile.bg_5prime_anchored    = (bg5 >= 0.0f) ? bg5 : 0.0f;
        profile.bg_3prime_anchored    = (bg3 >= 0.0f) ? bg3 : 0.0f;
        profile.bg_n_positions_5prime = n5;
        profile.bg_n_positions_3prime = n3;
        profile.bg_denominator_5prime = d5;
        profile.bg_denominator_3prime = d3;

        // Inline chemistry-tag detection: top 5' hexamer matches a curated
        // protocol-tag table at log2fc >= 3.0. Populates profile.protocol_tag_*
        // here (early) so the Briggs fit can mask pos 0; the post-bias-gate
        // rescue further down only mutates library_type if needed.
        if (profile.protocol_tag_5prime.empty()) {
            auto enriched = compute_hex_enriched_5prime(profile, 3.0f);
            if (!enriched.empty()) {
                auto seq = decode_hex(enriched[0].idx);
                const ProtocolTag* tag = lookup_protocol_tag(seq.data());
                if (tag) {
                    profile.protocol_tag_5prime   = std::string(seq.data(), 6);
                    profile.protocol_tag_protocol = tag->protocol;
                    profile.protocol_tag_class    = tag->klass;
                    profile.protocol_tag_log2fc   = static_cast<float>(enriched[0].log2fc);
                    profile.protocol_tag_log_lr   = tag->log_lr;
                }
            }
        }
        const bool tag5 = !profile.protocol_tag_5prime.empty()
                          && profile.protocol_tag_log2fc >= 3.0f;
        profile.briggs_pos0_masked_5prime = tag5;
        // No 3' chemistry-tag table currently; mirror flag stays false.
        profile.briggs_pos0_masked_3prime = false;

        const int  start5 = profile.briggs_pos0_masked_5prime ? 1 : 0;
        const int  start3 = profile.briggs_pos0_masked_3prime ? 1 : 0;
        const float fb5   = (bg5 >= 0.0f) ? bg5 : -1.0f;
        const float fb3   = (bg3 >= 0.0f) ? bg3 : -1.0f;

        estimate_briggs_params(profile.damage_rate_5prime, profile.max_damage_5prime,
                               start5, fb5,
                               profile.delta_s_5prime, profile.delta_d_5prime,
                               profile.lambda_5prime, profile.r_squared_5prime);

        estimate_briggs_params(profile.damage_rate_3prime, profile.max_damage_3prime,
                               start3, fb3,
                               profile.delta_s_3prime, profile.delta_d_3prime,
                               profile.lambda_3prime, profile.r_squared_3prime);

        // Headline area-excess + LR companion vs bg-only null over k..14.
        // Uses RAW per-position T/(T+C) rates (t_freq_5prime / a_freq_3prime,
        // rate-valued post-finalize) and the tail-anchored bg in the same
        // P-space, so area_excess is comparable across libraries without
        // depending on the joint Briggs fit's interior baseline.
        auto compute_headline = [](const std::array<double, 15>& rates,
                                   const std::array<double, 15>& tc_total,
                                   int k, float bg,
                                   float& area_excess, float& lr) {
            area_excess = 0.0f;
            lr = 0.0f;
            if (bg < 0.0f) bg = 0.0f;
            const double bg_safe = std::clamp(static_cast<double>(bg), 1e-6, 1.0 - 1e-6);
            for (int p = k; p < 15; ++p) {
                const double r = rates[p];
                if (r > bg) area_excess += static_cast<float>(r - bg);
                const double n = tc_total[p];
                if (n < 5.0) continue;
                const double k_obs = r * n;
                const double p_alt = std::clamp(r, 1e-6, 1.0 - 1e-6);
                const double ll_alt  = k_obs * std::log(p_alt) + (n - k_obs) * std::log(1.0 - p_alt);
                const double ll_null = k_obs * std::log(bg_safe) + (n - k_obs) * std::log(1.0 - bg_safe);
                if (ll_alt > ll_null) lr += static_cast<float>(ll_alt - ll_null);
            }
        };
        compute_headline(profile.t_freq_5prime, profile.tc_total_5prime,
                         start5, profile.bg_5prime_anchored,
                         profile.damage_5prime_area_excess, profile.damage_5prime_lr);
        compute_headline(profile.a_freq_3prime, profile.ag_total_3prime,
                         start3, profile.bg_3prime_anchored,
                         profile.damage_3prime_area_excess, profile.damage_3prime_lr);
    }

    profile.lambda_5prime = std::clamp(profile.lambda_5prime, 0.1f, 1.0f);
    profile.lambda_3prime = std::clamp(profile.lambda_3prime, 0.1f, 1.0f);

    // Library type: handled below after composition-bias flags are set

    // Flag inversions: terminals depleted relative to interior with statistical support
    // Use both shift threshold AND z-score override for very strong statistical signals
    const bool inversion_5 = (profile.terminal_shift_5prime < -0.01f && profile.terminal_z_5prime < -2.0f)
                           || profile.terminal_z_5prime < -10.0f;  // Very strong depletion signal
    const bool inversion_3 = (profile.terminal_shift_3prime < -0.01f && profile.terminal_z_3prime < -2.0f)
                           || profile.terminal_z_3prime < -10.0f;  // Very strong depletion signal
    profile.terminal_inversion = inversion_5 || inversion_3;

    // Set individual inverted pattern flags (used for asymmetric pattern handling)
    // Note: inverted_pattern_5prime may also be set later by hexamer-based detection
    if (inversion_3) {
        profile.inverted_pattern_3prime = true;
    }
    if (inversion_5) {
        profile.inverted_pattern_5prime = true;
    }

    // Hexamer-based damage: compare C-initial vs T-initial hexamer frequencies at terminal vs interior
    if (profile.n_hexamers_5prime >= 1000 && profile.n_hexamers_interior >= 1000) {
        double llr_sum = 0.0;
        double weight_sum = 0.0;

        for (uint32_t base_code = 0; base_code < 1024; ++base_code) {
            uint32_t c_hex = 0x400 | base_code;
            uint32_t t_hex = 0xC00 | base_code;

            float expected_c = get_hexamer_freq(c_hex);
            float expected_t = get_hexamer_freq(t_hex);

            if (expected_c < 1e-6f || expected_t < 1e-6f) continue;

            double term_c = profile.hexamer_count_5prime[c_hex];
            double term_t = profile.hexamer_count_5prime[t_hex];
            double int_c = profile.hexamer_count_interior[c_hex];
            double int_t = profile.hexamer_count_interior[t_hex];

            double total_term = term_c + term_t;
            double total_int = int_c + int_t;
            if (total_term < 5 || total_int < 5) continue;

            double ratio_term = term_t / total_term;
            double ratio_int = int_t / total_int;
            double excess = ratio_term - ratio_int;

            double weight = (expected_c + expected_t) * std::sqrt(total_term * total_int);
            llr_sum += excess * weight;
            weight_sum += weight;
        }

        if (weight_sum > 0) {
            float raw_llr = static_cast<float>(llr_sum / weight_sum);

            // Store raw hexamer LLR without inversion correction
            // The raw value is the actual terminal-vs-interior hexamer composition difference
            // Positive = T-starting hexamers enriched at terminal (damage pattern)
            // Negative = C-starting hexamers enriched at terminal (AT-rich composition bias)
            profile.hexamer_damage_llr = raw_llr;

            (void)raw_llr;
        }

        // Compute hexamer-based T/(T+C) ratios (more reliable than position 0 or 1 alone)
        double total_term_t = 0, total_term_c = 0;
        double total_int_t = 0, total_int_c = 0;
        for (uint32_t base_code = 0; base_code < 1024; ++base_code) {
            uint32_t c_hex = 0x400 | base_code;
            uint32_t t_hex = 0xC00 | base_code;
            total_term_c += profile.hexamer_count_5prime[c_hex];
            total_term_t += profile.hexamer_count_5prime[t_hex];
            total_int_c += profile.hexamer_count_interior[c_hex];
            total_int_t += profile.hexamer_count_interior[t_hex];
        }
        double term_t_ratio = total_term_t / (total_term_t + total_term_c + 1e-10);
        double int_t_ratio = total_int_t / (total_int_t + total_int_c + 1e-10);

        profile.hexamer_terminal_tc = static_cast<float>(term_t_ratio);
        profile.hexamer_interior_tc = static_cast<float>(int_t_ratio);
        profile.hexamer_excess_tc = static_cast<float>(term_t_ratio - int_t_ratio);

        // Hexamers starting at pos 0 include the pos-0 artifact, so skip inversion
        // detection when a pos-0 artifact is present.
        if (profile.hexamer_damage_llr < -0.02f && !profile.position_0_artifact_5prime) {
            profile.inverted_pattern_5prime = true;
        }
    }

    // Composition bias: flag if |ctrl_shift| >= max(0.005, 0.5 * |damage_shift|)
    float damage_shift_5 = profile.terminal_shift_5prime;
    float ctrl_shift_5 = std::abs(profile.ctrl_shift_5prime);
    float threshold_5 = std::max(0.005f, 0.5f * std::abs(damage_shift_5));
    if (ctrl_shift_5 >= threshold_5 && damage_shift_5 > 0.01f) {
        profile.composition_bias_5prime = true;
    }

    float damage_shift_3 = profile.terminal_shift_3prime;
    float ctrl_shift_3 = std::abs(profile.ctrl_shift_3prime);
    float threshold_3 = std::max(0.005f, 0.5f * std::abs(damage_shift_3));
    if (ctrl_shift_3 >= threshold_3 && damage_shift_3 > 0.01f) {
        profile.composition_bias_3prime = true;
    }

    // Library type detection.
    // DS: C→T at 5' + G→A at 3' (terminal_shift_3prime elevated).
    // SS: C→T at both ends — 3' shows elevated T/(T+C) (ctrl_shift_3prime) with no G→A.
    //
    // Library type detection via 4-model BIC comparison on the 3' end (positions 1-10).
    //
    // Four competing models:
    //   M_bias: no 3' decay in either channel (composition/ligation bias only)
    //   M_DS:   G→A decay at 3' only  (classic double-stranded aDNA)
    //   M_SS:   C→T decay at 3' only  (single-stranded: damage at both ends)
    //   M_mix:  both channels show decay (ambiguous / mixed library)
    //
    // Lambda is fixed to the fitted 5' C→T decay, so coverage inflating z-scores
    // does not affect the classification — only whether the decay shape fits.
    // Position 0 is excluded entirely (known adapter ligation artifact at 3' in SS prep).
    // P0-1: BIC tournament runs UNCONDITIONALLY. The forced-library override is
    // applied AFTER, by overwriting profile.library_type only — diagnostic
    // BIC fields, posteriors, and submodel scores remain populated so KapK-forced
    // runs produce identical numbers to KapK-auto runs.
    {
        float lambda_lib = std::clamp(fit_lambda_5p, 0.05f, 0.50f);

        // Four channels, all with lambda fixed to the fitted 5' C→T decay rate.
        //
        // ct5: T/(T+C) at 5' pos 1-10  — 5' C→T damage (DS reads and SS-original reads).
        // ga3: A/(A+G) at 3' pos 1-10  — smooth 3' G→A decay (DS reads and SS-complement reads).
        // ga0: A/(A+G) at 3' pos 0     — pos-0 G→A spike (SS-complement ligation artifact).
        // ct3: T/(T+C) at 3' pos 1-10  — 3' C→T (SS-original reads only; absent in DS).
        //
        // Seven competing joint models:
        //
        //   M_bias         = ct5_null + ga3_null  + ga0_null + ct3_null  (no damage anywhere)
        //   M_DS_symm      = ds_symm  + ga0_null  + ct3_null  (DS: ct5+ga3 joint fit, symmetric decay)
        //   M_DS_spike     = ct5_null + ga3_null  + ga0_alt  + ct3_null  (DS: only pos-0 end-repair spike)
        //   M_DS_symm_art  = ds_symm  + ga0_alt   + ct3_null  (DS: symmetric decay + pos-0 artifact)
        //   M_SS_comp      = ct5_null + ga3_alt   + ga0_alt  + ct3_null  (SS complement-orientation reads)
        //   M_SS_orig      = ct5_alt  + ga3_null  + ga0_null + ct3_alt   (SS original-orientation reads)
        //   M_SS_full      = ct5_alt  + ga3_alt   + ga0_alt  + ct3_alt   (SS both orientations)
        //
        // M_DS_symm, M_DS_spike, M_DS_symm_art classify as DS; M_SS_* as SS.
        //
        // M_DS_symm: ct5 and ga3 are jointly fitted with a SINGLE shared amplitude. For genuine DS,
        // ct5 ≈ ga3 and A_joint is near-optimal for both → low BIC. For SS_comp (ct5 << ga3), A_joint
        // is too high for ct5 (over-predicts → LL_ct5 < LL_null) and too low for ga3 → combined BIC
        // exceeds M_SS_comp → M_DS_symm rejected. This replaces the earlier threshold-based approach.
        //
        // M_DS_spike covers DS libraries with composition bias at 5' where positions 1-10 show no
        // usable ct5/ga3 signal (amplitudes clamped to 0), yet a genuine chemical pos-0 signal exists
        // from end-repair. Without M_DS_spike, M_SS_comp would win for such samples.
        //
        // Joint BIC = sum of per-channel BICs (bic_alt for alt channels, bic_null for null).
        // M_DS_symm uses 1 free parameter for both ct5+ga3; all other models use per-channel params.
        // Winner = lowest joint BIC. No post-hoc thresholds.
        // Offset search: only applied when an adapter artifact or terminal inversion was detected
        // (position_0_artifact or inverted_pattern). For normal samples always use start_pos=1
        // to avoid perturbing the BIC scores used in library-type classification.
        const bool artifact_5 = profile.position_0_artifact_5prime || profile.inverted_pattern_5prime;
        const bool artifact_3 = profile.position_0_artifact_3prime || profile.inverted_pattern_3prime;
        ChannelDecayFit ct5;  int ct5_offset = 1;
        ChannelDecayFit ga3;  int ga3_offset = 1;
        ChannelDecayFit ct3;  int ct3_offset = 1;
        ChannelDecayFit ds_symm;  int ds_symm_offset = 1;
        if (artifact_5) {
            std::tie(ct5, ct5_offset) = fit_decay_best_offset(
                profile.t_freq_5prime, profile.tc_total_5prime,
                static_cast<float>(baseline_tc), lambda_lib, 10);
        } else {
            ct5 = fit_decay_fixed_lambda(profile.t_freq_5prime, profile.tc_total_5prime,
                static_cast<float>(baseline_tc), lambda_lib, 1, 10);
        }
        if (artifact_3) {
            std::tie(ga3, ga3_offset) = fit_decay_best_offset(
                profile.a_freq_3prime, profile.ag_total_3prime,
                static_cast<float>(baseline_ag), lambda_lib, 10);
            // ct3 does NOT use offset search even when artifact_3 is true: the 3' adapter
            // artifact suppresses GA3 (complement-strand G→A) but not CT3 (original-strand
            // C→T). Offset search for ct3 finds spurious C→T signals in adapter-affected DS
            // libraries and drives false SS-original classifications.
            ct3 = fit_decay_fixed_lambda(ctrl_freq_3p, ctrl_total_3p,
                static_cast<float>(baseline_tc), lambda_lib, 1, 10);
        } else {
            ga3 = fit_decay_fixed_lambda(profile.a_freq_3prime, profile.ag_total_3prime,
                static_cast<float>(baseline_ag), lambda_lib, 1, 10);
            ct3 = fit_decay_fixed_lambda(ctrl_freq_3p, ctrl_total_3p,
                static_cast<float>(baseline_tc), lambda_lib, 1, 10);
        }
        // ga0: single-position fit at 3' pos-0; no offset search (it IS the pos-0 spike signal)
        ChannelDecayFit ga0 = fit_decay_fixed_lambda(
            profile.a_freq_3prime, profile.ag_total_3prime,
            static_cast<float>(baseline_ag), lambda_lib, 0, 0, 1);
        // Symmetric DS model: offset search only when either end has artifact.
        // Restrict high-lambda candidates when ga0 dominates ct5 (SS-complement indicator):
        // the high-lambda CT5 fit should not inflate the DS joint BIC when the sample already
        // shows a large GA0 spike (spike_is_ss-like pattern). ga0 is fit before this point.
        const bool restrict_joint_lambda = (ga0.amplitude >= 0.10f && ga0.delta_bic > ct5.delta_bic);
        // P2: spike-gate diagnostic — record the joint-lambda gating decision
        profile.library_joint_lambda_restricted = restrict_joint_lambda;
        if (artifact_5 || artifact_3) {
            std::tie(ds_symm, ds_symm_offset) = fit_decay_joint_best_offset(
                profile.t_freq_5prime, profile.tc_total_5prime, static_cast<float>(baseline_tc),
                profile.a_freq_3prime, profile.ag_total_3prime, static_cast<float>(baseline_ag),
                lambda_lib, 10, 3, restrict_joint_lambda);
        } else {
            ds_symm = fit_decay_fixed_lambda_joint(
                profile.t_freq_5prime, profile.tc_total_5prime, static_cast<float>(baseline_tc),
                profile.a_freq_3prime, profile.ag_total_3prime, static_cast<float>(baseline_ag),
                lambda_lib, 1, 10);
        }

        // Store detected offsets; use the single-channel offsets as the canonical per-end values
        // (ds_symm_offset is derived from the joint fit and may differ from individual channels)
        profile.fit_offset_5prime = ct5_offset;
        profile.fit_offset_3prime = ga3_offset;

        // When adapter artifact or terminal inversion is detected, scan positions 1-5 for the peak
        // damage rate. Using damage_rate[fit_offset] is unreliable because the BIC-best offset may
        // point to a position still below baseline (e.g. BPN103cm: fit_offset=3 but pos3 < bg,
        // while actual biological damage is at pos1). Scan-for-peak handles all shift lengths.
        if (profile.position_0_artifact_5prime || profile.inverted_pattern_5prime) {
            float peak = 0.0f;
            for (int p = 1; p <= 5 && p < 15; ++p) {
                if (profile.damage_rate_5prime[p] > peak) peak = profile.damage_rate_5prime[p];
            }
            if (peak > 0.0f) profile.max_damage_5prime = peak;
        }
        if (profile.position_0_artifact_3prime || profile.inverted_pattern_3prime) {
            float peak = 0.0f;
            for (int p = 1; p <= 5 && p < 15; ++p) {
                if (profile.damage_rate_3prime[p] > peak) peak = profile.damage_rate_3prime[p];
            }
            if (peak > 0.0f) profile.max_damage_3prime = peak;
        }
        (void)ct3_offset; (void)ds_symm_offset;

        profile.libtype_amp_ct5  = ct5.amplitude;
        profile.libtype_amp_ga3  = ga3.amplitude;
        profile.libtype_amp_ga0  = ga0.amplitude;
        profile.libtype_amp_ct3  = ct3.amplitude;
        profile.libtype_dbic_ct5 = ct5.delta_bic;
        profile.libtype_dbic_ga3 = ga3.delta_bic;
        profile.libtype_dbic_ga0 = ga0.delta_bic;
        profile.libtype_dbic_ct3 = ct3.delta_bic;

        if (ct5.valid && ga3.valid && ga0.valid && ct3.valid && ds_symm.valid) {
            // Hard biological gates (no tunable parameters):
            //   M_DS_symm:     DS damage is symmetric — requires ga3 to have real
            //                  signal. Use the standard Bayes-factor "no evidence"
            //                  threshold: ga3.delta_bic > log(2) ≈ 0.693
            //                  (Kass & Raftery 1995). If ga3 is indistinguishable
            //                  from null, there is nothing for ct5 to be symmetric
            //                  with — the joint ds_symm fit is fitting noise.
            //   M_DS_symm_art: 3' GA spike is the complementary-strand reflection
            //                  of 5' CT damage. A reflection cannot exceed its
            //                  source — require ga0 ≤ ct5.
            // When violated, the model is invalid for this sample and excluded
            // from the tournament (BIC = +inf).
            const bool ds_symm_valid     = (ga3.delta_bic > std::log(2.0));
            const bool ds_symm_art_valid = (ga0.amplitude <= ct5.amplitude);
            constexpr double kInvalidBIC = 1e300;

            const double bic_M_bias        = ct5.bic_null   + ga3.bic_null + ga0.bic_null + ct3.bic_null;
            const double bic_M_DS_symm     = ds_symm_valid     ? (ds_symm.bic_alt + ga0.bic_null + ct3.bic_null) : kInvalidBIC;
            const double bic_M_DS_spike    = ct5.bic_null   + ga3.bic_null + ga0.bic_alt  + ct3.bic_null;
            const double bic_M_DS_symm_art = ds_symm_art_valid ? (ds_symm.bic_alt + ga0.bic_alt  + ct3.bic_null) : kInvalidBIC;
            const double bic_M_SS_comp     = ct5.bic_null   + ga3.bic_alt  + ga0.bic_alt  + ct3.bic_null;
            const double bic_M_SS_orig     = ct5.bic_alt    + ga3.bic_null + ga0.bic_null + ct3.bic_alt;
            const double bic_M_SS_full     = ct5.bic_alt    + ga3.bic_alt  + ga0.bic_alt  + ct3.bic_alt;
            // S1 telemetry only: asymmetric DS counterfactual (independent ct5/ga3 amps,
            // GA0 artifact, no CT3). Never enters cascade/softmax.
            const double bic_M_DS_asym_art = ct5.bic_alt    + ga3.bic_alt  + ga0.bic_alt  + ct3.bic_null;
            // S1 invariant probe: M_DS_symm_art rebuilt from a forced no-offset joint fit
            // (start_pos=1, no joint best-offset search) regardless of artifact_5/artifact_3.
            // Catches future regressions where joint best-offset search inflates ds_symm.
            ChannelDecayFit ds_symm_no_off = fit_decay_fixed_lambda_joint(
                profile.t_freq_5prime, profile.tc_total_5prime, static_cast<float>(baseline_tc),
                profile.a_freq_3prime, profile.ag_total_3prime, static_cast<float>(baseline_ag),
                lambda_lib, 1, 10);
            const double bic_M_DS_symm_art_no_offset = (ds_symm_no_off.valid
                ? ds_symm_no_off.bic_alt + ga0.bic_alt + ct3.bic_null
                : 0.0);
            // SS asymmetric: original-orientation CT5 + complement-orientation GA0 spike; no GA3 smooth decay.
            // Wins over M_DS_symm_art when ga3.delta_bic < log(2), i.e. ga3 has no real signal.
            // Only considered when spike_is_ss=true so it cannot compete with M_DS_spike in the DS-only path.
            const double bic_M_SS_asym     = ct5.bic_alt    + ga3.bic_null + ga0.bic_alt  + ct3.bic_null;

            // ga0 amplitude distinguishes DS end-repair artifact (<0.10) from SS ligation spike (>=0.10).
            // Exception: if the channel-B structural analysis (computed before BIC, line ~1200)
            // confirmed bilateral symmetric damage with d_max_from_channel_b > 0.10, then the
            // 3' pos-0 GA excess is the mirror of the bilateral 5' C→T on the bottom strand —
            // not an SS ligation artifact.  Without this guard, highly damaged DS libraries
            // (d_max > 0.20, steep lambda) produce large GA0_ΔBIC from the 3' bilateral
            // reflection and are wrongly classified as SS.
            const bool structural_bilateral = profile.channel_b_quantifiable
                                           && (profile.d_max_from_channel_b > 0.10f);
            const bool spike_is_ss = (ga0.amplitude >= 0.10f) && !structural_bilateral;
            // P2: spike-gate diagnostics — record the gating inputs and decision
            profile.library_spike_is_ss                    = spike_is_ss;
            profile.library_spike_gate_ga0_amp             = ga0.amplitude;
            profile.library_spike_gate_structural_bilateral = structural_bilateral;
            // P1-B: store all 7 submodel BICs and active-flags so callers can audit the tournament
            profile.library_bic_M_bias        = bic_M_bias;
            profile.library_bic_M_DS_symm     = bic_M_DS_symm;
            profile.library_bic_M_DS_spike    = bic_M_DS_spike;
            profile.library_bic_M_DS_symm_art = bic_M_DS_symm_art;
            profile.library_bic_M_SS_comp     = bic_M_SS_comp;
            profile.library_bic_M_SS_orig     = bic_M_SS_orig;
            profile.library_bic_M_SS_asym     = bic_M_SS_asym;
            profile.library_bic_M_SS_full     = bic_M_SS_full;
            profile.library_M_SS_orig_active  = (ct3.delta_bic > 0.0);
            profile.library_M_SS_asym_active  = spike_is_ss;
            // S1 telemetry: counterfactual + invariant probe BICs (never enter cascade)
            profile.library_bic_M_DS_asym_art          = bic_M_DS_asym_art;
            profile.library_bic_M_DS_symm_art_no_offset = bic_M_DS_symm_art_no_offset;
            // S1 telemetry: gate inputs / per-channel offsets / ds_symm joint diagnostics
            profile.library_gate_artifact_5                   = artifact_5;
            profile.library_gate_artifact_3                   = artifact_3;
            profile.library_gate_position_0_artifact_5prime   = profile.position_0_artifact_5prime;
            profile.library_gate_position_0_artifact_3prime   = profile.position_0_artifact_3prime;
            profile.library_gate_inverted_pattern_5prime      = profile.inverted_pattern_5prime;
            profile.library_gate_inverted_pattern_3prime      = profile.inverted_pattern_3prime;
            profile.library_gate_max_damage_5prime            = profile.max_damage_5prime;
            profile.library_gate_structural_bilateral         = structural_bilateral;
            profile.library_gate_ga0_dominates_ct5            =
                (ga0.amplitude >= 0.10f && ga0.delta_bic > ct5.delta_bic);
            profile.library_ct3_offset      = ct3_offset;
            profile.library_ds_symm_offset  = ds_symm_offset;
            profile.library_ds_symm_lambda_used = ds_symm.lambda;
            profile.library_ds_symm_amp     = ds_symm.amplitude;
            profile.library_ds_symm_ct5_resid = ct5.amplitude - ds_symm.amplitude;
            profile.library_ds_symm_ga3_resid = ga3.amplitude - ds_symm.amplitude;
            const double best_ds = std::min({bic_M_DS_symm,
                                             bic_M_DS_symm_art,
                                             spike_is_ss ? 1e300 : bic_M_DS_spike});
            // M_SS_full excluded: 4-param model unfairly defeats M_DS_symm_art for asymmetric DS.
            // M_SS_asym only enters the SS set when spike_is_ss, preventing competition with M_DS_spike
            // in non-spike samples (which would misclassify one-sided DS libraries).
            const double best_ss = std::min({bic_M_SS_comp,
                                             ct3.delta_bic > 0.0 ? bic_M_SS_orig : 1e300,
                                             spike_is_ss ? bic_M_DS_spike : 1e300,
                                             spike_is_ss ? bic_M_SS_asym  : 1e300});
            profile.library_bic_bias = bic_M_bias;
            profile.library_bic_ds   = best_ds;
            profile.library_bic_ss   = best_ss;
            profile.library_bic_mix  = bic_M_SS_full;

            // Softmax posterior probabilities over all 7 candidate models:
            //   P(M_i | data) ∝ exp(-BIC_i/2)
            // Subtracting min(BIC) before exp() keeps everything finite under
            // the typical BIC range we see (Δ ~ 1e3–1e6 between models). Sum
            // by class — DS includes M_DS_symm/spike/symm_art; SS includes
            // M_SS_comp/orig/full/asym (the latter two only when active).
            // M_DS_spike contributes to SS instead when spike_is_ss=true,
            // matching the cascade's interpretation of bilateral pos-0 GA.
            double best = bic_M_bias;
            profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            bool ds_spike_won = false;  // tracks whether M_DS_spike is the current winning model

            if (bic_M_DS_symm < best) {
                best = bic_M_DS_symm;
                profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            }
            if (!spike_is_ss && bic_M_DS_spike < best) {
                best = bic_M_DS_spike;
                profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
                ds_spike_won = true;
            }
            if (bic_M_DS_symm_art < best) {
                best = bic_M_DS_symm_art;
                profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
                ds_spike_won = false;
            }
            if (bic_M_SS_comp < best) {
                best = bic_M_SS_comp;
                profile.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
                ds_spike_won = false;
            }
            // M_SS_orig requires ct3 signal: SS original-orientation reads produce CT3 whenever
            // they produce CT5. Without CT3, a one-sided DS pattern is more likely.
            if (bic_M_SS_orig < best && ct3.delta_bic > 0.0) {
                best = bic_M_SS_orig;
                profile.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
                ds_spike_won = false;
            }
            if (spike_is_ss && bic_M_DS_spike < best) {
                best = bic_M_DS_spike;
                profile.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
                ds_spike_won = false;
            }
            // M_SS_asym: SS with CT5 from original-orientation + GA0 spike from complement-orientation,
            // but no detectable GA3 smooth decay (ga3.delta_bic ≈ 0). Analytically beats M_DS_symm_art
            // iff ga3.delta_bic < log(2) ≈ 0.693, which is the gap the joint-fit BIC penalty creates.
            if (spike_is_ss && bic_M_SS_asym < best) {
                best = bic_M_SS_asym;
                profile.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            }
            // Post-hoc symmetry check: DS_symm constrains ct5_amp ≈ ga3_amp.
            // If DS wins but CT5 ΔBIC / GA3 ΔBIC < 0.50, the winning model's own symmetry
            // assumption is violated → reclassify as SS. Guard ga3.delta_bic > 3e4 to avoid
            // misfiring on low-damage DS libraries where small asymmetry is noise.
            // Symmetry veto: DS won BIC, but GA3 >> CT5 asymmetry suggests SS.
            // Only apply when there is positive biological evidence for SS orientation:
            //   ga0.delta_bic > 1e4 — meaningful GA0 ligation spike (amplitude-independent,
            //                         avoids false negatives from samples near the 0.10 amplitude
            //                         threshold that spike_is_ss uses)
            //   inverted_5prime     — 5' CT profile depressed below interior (SS adapter suppression)
            // Without either, the asymmetry likely reflects 3'-adapter inflation of GA3_ΔBIC in a DS
            // library (e.g., high-λ fit capturing sharp pos2+ GA3 decay after 3' adapter suppression).
            const bool ss_orientation_evidence = (ga0.delta_bic > 1e4)
                                              || profile.inverted_pattern_5prime;
            if (profile.library_type == SampleDamageProfile::LibraryType::DOUBLE_STRANDED &&
                ss_orientation_evidence &&
                ga3.delta_bic > 3e4 &&
                ct5.delta_bic / ga3.delta_bic < 0.50) {
                profile.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
                ds_spike_won = false;
            }
            // M_DS_spike rescue: a GA0 pos-0 spike with no CT5 and no GA3 smooth decay could be
            // M_DS_spike (DS end-repair, bilateral: both 5' pos0 and 3' pos0 elevated) or
            // complement-orientation SS (3' GA0 spike only, no 5' pos-0 counterpart).
            // Discriminate with max_damage_5prime (background-corrected excess CT at 5' pos-0).
            // Requires ga0.amplitude > 0.02 to exclude near-zero noise spikes (e.g. tiny end-repair
            // artifacts with ga0_amp ≈ 0.001 that also have negligible d_max_5 ≈ 0.005).
            // DS bilateral artifacts have d_max_5 ≈ ga0_amp >> 0.005; SS complement-only has d_max_5 ≈ 0.
            // Restricted to ds_spike_won: M_DS_symm_art can win via joint fit even with marginal
            // ct5/ga3 delta_bic ≤ 0; rescue only fires when M_DS_spike was the actual winner.
            if (profile.library_type == SampleDamageProfile::LibraryType::DOUBLE_STRANDED &&
                ds_spike_won &&
                !spike_is_ss &&
                ct5.delta_bic <= 0.0 &&
                ga3.delta_bic <= 0.0 &&
                ga0.delta_bic > 0.0 &&
                ga0.amplitude > 0.02f &&
                profile.max_damage_5prime <= 0.005f) {
                profile.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            }
            // GA0 bilateral rescue (spike_is_ss path): when ga0.amplitude >= 0.10 and DS wins,
            // discriminate DS end-repair (bilateral: both 5' pos-0 CT and 3' pos-0 GA elevated)
            // from SS complement-orientation reads (3' GA0 spike only, no 5' CT0 counterpart).
            // Validated on 24 DS controls with ga0_amp >= 0.10: all have d5 >= 0.11.
            // All SS mixed failures with ga0_amp >= 0.10 have d5 = 0. Gap is >20x the threshold.
            if (profile.library_type == SampleDamageProfile::LibraryType::DOUBLE_STRANDED &&
                spike_is_ss &&
                ga0.delta_bic > 0.0 &&
                profile.max_damage_5prime <= 0.005f) {
                profile.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            }
            // Channel-B structural DS rescue: if channel-B (stop codon conversion, immune to
            // composition bias) confirmed bilateral symmetric damage > 0.10, AND the 3' ligation
            // spike (GA0) dominates the smooth GA3 decay by ≥5×, override SS→DS.
            // Handles highly damaged DS libraries where steep lambda suppresses CT5 in the BIC
            // fit, causing M_SS_comp to win by capturing GA3+GA0 independently. Channel B is
            // only quantifiable when real C→T bilateral damage exists at 5'; SS-complement
            // libraries have ch_b_quant=False and are not affected.
            // The GA0 dominance guard (ga3 * 5 < ga0) excludes SS-orig libraries: they have
            // real bilateral-looking channels (ch_b_quant=True, large GA3) but their GA3 decay
            // is comparable to GA0 rather than being dwarfed by it.
            const bool ga_spike_dominant = (ga3.delta_bic * 5.0 < ga0.delta_bic);
            if (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED
                && structural_bilateral && ga_spike_dominant
                && ct3.delta_bic <= 0.0) {
                profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
                profile.library_type_rescued = true;
            }
            // GA0-spike DS-symm veto: when the 3' ga0 ligation spike unambiguously
            // dominates (ga0 ≥ 0.10, ga0_dominates_ct5, ga0 > both ct5 and ga3),
            // an apparent CT5/GA3 symmetry absorbed by M_DS_symm is artifact, not
            // real DS damage — true DS has ct5≈ga3 ≫ ga0. Restricted to ds_symm
            // winners; independent of max_damage_5prime (the artifact inflates it).
            if (profile.library_type == SampleDamageProfile::LibraryType::DOUBLE_STRANDED &&
                profile.library_gate_ga0_dominates_ct5 &&
                ga0.amplitude >= 0.10f &&
                ga0.amplitude > std::max(ct5.amplitude, ga3.amplitude) &&
                (best == bic_M_DS_symm || best == bic_M_DS_symm_art)) {
                profile.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            }
            // Channel-B DS rescue from M_SS_comp: structural_bilateral confirms
            // bilateral symmetric damage (channel B is artifact-immune); ct5≈ga3
            // (within 30%, both ≥ 0.03) confirms the symmetric DS pair. M_SS_comp
            // wins by capturing ct5+ga3+ga0 as 3 independent bumps; M_DS_symm has
            // no ga0 term and loses on raw BIC despite being the correct model
            // (the ga0 here is end-repair / ligation residue co-occurring with
            // real DS damage, not the SS complement-orientation signature).
            if (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED &&
                structural_bilateral &&
                best == bic_M_SS_comp &&
                ct3.delta_bic <= 0.0 &&
                ct5.amplitude >= 0.03f &&
                ga3.amplitude >= 0.03f &&
                std::abs(ct5.amplitude - ga3.amplitude) <=
                    0.30f * std::max(ct5.amplitude, ga3.amplitude)) {
                profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
                profile.library_type_rescued = true;
            }
            // Low-amp symmetric DS rescue: zero 3' damage collapses the per-end
            // DS amplitudes so M_SS_full / M_SS_orig wins raw BIC, but the joint
            // M_DS_symm fit is still valid (ct5 ~= ga3, residuals ~= 0). ga0 < 0.005
            // and !spike_is_ss exclude SS-complement orientation.
            if (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED &&
                !spike_is_ss &&
                ds_symm.amplitude >= 0.005f &&
                ct5.amplitude >= 0.005f && ga3.amplitude >= 0.005f &&
                ga0.amplitude < 0.005f &&
                std::abs(ct5.amplitude - ga3.amplitude) <=
                    0.50f * std::max(ct5.amplitude, ga3.amplitude) &&
                profile.max_damage_3prime < 0.005f) {
                profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
                profile.library_type_rescued = true;
            }
            // S1 telemetry: capture cascade-derived gate outputs.
            profile.library_gate_ss_orientation_evidence = ss_orientation_evidence;
            profile.library_gate_ga_spike_dominant       = ga_spike_dominant;
            profile.library_gate_ds_spike_won            = ds_spike_won;
            // Uninformative: if nothing beat M_bias (best unchanged), no damage channel
            // provided evidence for either DS or SS. Use exact equality — best is only
            // updated via assignment from another BIC value, so equality is safe here.
            // Does NOT affect low-damage DS libraries where M_DS_symm still beats M_bias.
            if (best == bic_M_bias) {
                profile.library_type = SampleDamageProfile::LibraryType::UNKNOWN;
            }
            // Protocol-tag rescue: top 5' hexamer is a chemistry fingerprint
            // (deterministic per library prep, independent of damage shape).
            // Runs AFTER the bias gate so chemistry evidence overrides UNKNOWN
            // for low-damage libraries where BIC has nothing to fit. Only fires
            // when log2fc >= 3.0 (8x enrichment) and the hex matches the curated
            // table, so noise hexamers cannot trip it.
            {
                auto enriched = compute_hex_enriched_5prime(profile, 3.0f);
                if (!enriched.empty()) {
                    auto seq = decode_hex(enriched[0].idx);
                    const ProtocolTag* tag = lookup_protocol_tag(seq.data());
                    if (tag) {
                        profile.protocol_tag_5prime   = std::string(seq.data(), 6);
                        profile.protocol_tag_protocol = tag->protocol;
                        profile.protocol_tag_class    = tag->klass;
                        profile.protocol_tag_log2fc   = static_cast<float>(enriched[0].log2fc);
                        profile.protocol_tag_log_lr   = tag->log_lr;
                        if (profile.library_type != tag->klass) {
                            profile.library_type        = tag->klass;
                            profile.library_type_rescued = true;
                            profile.protocol_tag_applied = true;
                            profile.library_artifact_reasons.push_back("protocol_tag_5prime");
                        }
                    }
                }
            }
            // Posterior class probabilities P(class | data) ∝ exp(-BIC/2),
            // computed AFTER the cascade so M_DS_spike (and M_SS_asym) are
            // routed to whichever class the cascade's domain rescues actually
            // assigned. Keeps p_ds/p_ss consistent with library_type.
            {
                // Hypothesis set must match the cascade's competition exactly:
                //   M_bias, M_DS_symm, M_DS_spike, M_DS_symm_art, M_SS_comp,
                //   M_SS_orig (only when ct3.delta_bic > 0),
                //   M_SS_asym (only when spike_is_ss).
                // M_SS_full is intentionally NOT in the cascade (its 4-param fit
                // unfairly defeats M_DS_symm_art on asymmetric DS) so it is
                // excluded from the softmax too — otherwise its mass would
                // appear in p_ss without ever being a valid winner.
                const bool inc_orig = (ct3.delta_bic > 0.0);
                const bool inc_asym = spike_is_ss;
                const double BIG = 1e300;
                const double bics[7] = {
                    bic_M_bias,
                    bic_M_DS_symm,
                    bic_M_DS_spike,
                    bic_M_DS_symm_art,
                    bic_M_SS_comp,
                    inc_orig ? bic_M_SS_orig : BIG,
                    inc_asym ? bic_M_SS_asym : BIG
                };
                double bmin = bics[0];
                for (int i = 1; i < 7; ++i) if (bics[i] < bmin) bmin = bics[i];
                double w[7];
                double Z = 0.0;
                for (int i = 0; i < 7; ++i) {
                    double half = (bmin - bics[i]) * 0.5;
                    if (half < -80.0) half = -80.0;
                    w[i] = std::exp(half);
                    Z += w[i];
                }
                if (Z > 0.0) {
                    for (int i = 0; i < 7; ++i) w[i] /= Z;
                    const double p_bias = w[0];
                    double p_ds = w[1] + w[3];                 // M_DS_symm + M_DS_symm_art
                    double p_ss = w[4] + w[5] + w[6];          // M_SS_comp + M_SS_orig (gated) + M_SS_asym (gated)
                    // Route ambiguous M_DS_spike (w[2]) to the class the cascade
                    // (with all rescues applied) ended up choosing. UNKNOWN/bias
                    // falls back to the spike_is_ss heuristic.
                    const auto LT = profile.library_type;
                    if (LT == SampleDamageProfile::LibraryType::SINGLE_STRANDED) {
                        p_ss += w[2];
                    } else if (LT == SampleDamageProfile::LibraryType::DOUBLE_STRANDED) {
                        p_ds += w[2];
                    } else {
                        if (spike_is_ss) p_ss += w[2]; else p_ds += w[2];
                    }
                    double sum = p_bias + p_ds + p_ss;
                    if (sum > 0.0) {
                        profile.library_p_bias   = static_cast<float>(p_bias / sum);
                        profile.library_p_ds     = static_cast<float>(p_ds   / sum);
                        profile.library_p_ss     = static_cast<float>(p_ss   / sum);
                        profile.library_p_winner = std::max({profile.library_p_bias,
                                                             profile.library_p_ds,
                                                             profile.library_p_ss});
                        profile.library_type_evaluable = true;
                    }
                }
            }
            // P1-A: winner / second-best model + margin and class-min softmax.
            // The candidate set MUST match the cascade exactly (same gating as
            // the existing posterior block above): always-on M_bias/M_DS_symm/
            // M_DS_spike/M_DS_symm_art/M_SS_comp; M_SS_orig only when ct3
            // shows positive evidence; M_SS_asym only when spike_is_ss.
            {
                struct CandM { const char* name; double bic; bool active; int klass; };
                // klass: 0=bias, 1=DS, 2=SS. M_DS_spike routes per cascade outcome.
                int spike_klass = 1;  // default DS
                if (spike_is_ss) {
                    spike_klass = 2;
                } else if (profile.library_type ==
                           SampleDamageProfile::LibraryType::SINGLE_STRANDED) {
                    spike_klass = 2;
                }
                CandM cands[7] = {
                    {"M_bias",        bic_M_bias,        true,                       0},
                    {"M_DS_symm",     bic_M_DS_symm,     true,                       1},
                    {"M_DS_spike",    bic_M_DS_spike,    true,                       spike_klass},
                    {"M_DS_symm_art", bic_M_DS_symm_art, true,                       1},
                    {"M_SS_comp",     bic_M_SS_comp,     true,                       2},
                    {"M_SS_orig",     bic_M_SS_orig,     ct3.delta_bic > 0.0,        2},
                    {"M_SS_asym",     bic_M_SS_asym,     spike_is_ss,                2},
                };
                int win_i = -1, sec_i = -1;
                for (int i = 0; i < 7; ++i) {
                    if (!cands[i].active) continue;
                    if (win_i < 0 || cands[i].bic < cands[win_i].bic) {
                        sec_i = win_i; win_i = i;
                    } else if (sec_i < 0 || cands[i].bic < cands[sec_i].bic) {
                        sec_i = i;
                    }
                }
                if (win_i >= 0) {
                    profile.library_bic_winner_model = cands[win_i].name;
                    if (sec_i >= 0) {
                        profile.library_bic_second_model = cands[sec_i].name;
                        profile.library_bic_margin = cands[sec_i].bic - cands[win_i].bic;
                    }
                }
                // Class-min softmax: best-BIC representative per class only.
                // 3-way over { best DS, best SS, M_bias }, subtract-min trick.
                double best_per_class[3] = {1e300, 1e300, 1e300};
                for (int i = 0; i < 7; ++i) {
                    if (!cands[i].active) continue;
                    int k = cands[i].klass;
                    if (cands[i].bic < best_per_class[k]) best_per_class[k] = cands[i].bic;
                }
                double bmin_c = std::min({best_per_class[0], best_per_class[1], best_per_class[2]});
                double w_c[3] = {0.0, 0.0, 0.0};
                double Z_c = 0.0;
                for (int k = 0; k < 3; ++k) {
                    if (best_per_class[k] >= 1e299) continue;
                    double half = (bmin_c - best_per_class[k]) * 0.5;
                    if (half < -80.0) half = -80.0;
                    w_c[k] = std::exp(half);
                    Z_c += w_c[k];
                }
                if (Z_c > 0.0) {
                    profile.library_p_bias_class_min = static_cast<float>(w_c[0] / Z_c);
                    profile.library_p_ds_class_min   = static_cast<float>(w_c[1] / Z_c);
                    profile.library_p_ss_class_min   = static_cast<float>(w_c[2] / Z_c);
                }
            }
            // S1 telemetry: raw 9-candidate ranking — ignores cascade gating, exposes
            // the absolute-best fit so callers can audit cascade exclusions.
            // 8 cascade contenders + M_DS_asym_art (telemetry-only counterfactual).
            {
                struct Raw { const char* name; double bic; const char* klass; bool in_cascade; };
                const Raw raw[9] = {
                    {"M_bias",        bic_M_bias,        "bias", true},
                    {"M_DS_symm",     bic_M_DS_symm,     "ds",   true},
                    {"M_DS_spike",    bic_M_DS_spike,    "ds",   true},
                    {"M_DS_symm_art", bic_M_DS_symm_art, "ds",   true},
                    {"M_SS_comp",     bic_M_SS_comp,     "ss",   true},
                    {"M_SS_orig",     bic_M_SS_orig,     "ss",   true},
                    {"M_SS_asym",     bic_M_SS_asym,     "ss",   true},
                    {"M_SS_full",     bic_M_SS_full,     "ss",   false},  // structurally excluded from cascade
                    {"M_DS_asym_art", bic_M_DS_asym_art, "ds",   false},  // telemetry-only counterfactual
                };
                int rwin = 0;
                for (int i = 1; i < 9; ++i) if (raw[i].bic < raw[rwin].bic) rwin = i;
                int rsec = -1;
                for (int i = 0; i < 9; ++i) {
                    if (i == rwin) continue;
                    if (rsec < 0 || raw[i].bic < raw[rsec].bic) rsec = i;
                }
                profile.library_bic_raw_winner_model = raw[rwin].name;
                profile.library_bic_raw_winner_class = raw[rwin].klass;
                profile.library_bic_raw_winner_in_cascade = raw[rwin].in_cascade;
                if (rsec >= 0) {
                    profile.library_bic_raw_second_model = raw[rsec].name;
                    profile.library_bic_raw_margin       = raw[rsec].bic - raw[rwin].bic;
                }
                // Cascade-exclusion booleans: which gating reasons could have
                // suppressed the raw winner from cascade competition. Multi-valued.
                profile.library_bic_excl_in_cascade            = raw[rwin].in_cascade;
                profile.library_bic_excl_M_SS_full_hardcoded   = (std::string(raw[rwin].name) == "M_SS_full");
                profile.library_bic_excl_ct3_zero              = (std::string(raw[rwin].name) == "M_SS_orig"
                                                                  && ct3.delta_bic <= 0.0);
                profile.library_bic_excl_spike_is_ss           = ((std::string(raw[rwin].name) == "M_SS_asym" && !spike_is_ss)
                                                                  || (std::string(raw[rwin].name) == "M_DS_spike" && spike_is_ss));
                profile.library_bic_excl_structural_bilateral  = (std::string(raw[rwin].name) == "M_DS_asym_art"
                                                                  && !structural_bilateral);
            }
            // Low-confidence override: if the posterior winner sits below the
            // confidence threshold AND no domain rescue rule fired, fall back
            // to UNKNOWN so downstream consumers (fqdup) use a neutral damage
            // prior instead of committing to a marginal call.
            if (!profile.library_type_rescued &&
                profile.library_p_winner > 0.0f &&
                profile.library_p_winner <
                    SampleDamageProfile::kLibraryTypeConfidenceThreshold) {
                profile.library_type = SampleDamageProfile::LibraryType::UNKNOWN;
            }
        } else {
            profile.library_bic_bias = 0.0;
            profile.library_bic_ds   = 0.0;
            profile.library_bic_ss   = 0.0;
            profile.library_bic_mix  = 0.0;
            profile.library_p_bias   = 0.0f;
            profile.library_p_ds     = 0.0f;
            profile.library_p_ss     = 0.0f;
            profile.library_p_winner = 0.0f;
            profile.library_type_evaluable = false;
            profile.library_type = SampleDamageProfile::LibraryType::UNKNOWN;
        }

        // Rescue rule: when BIC classification returns UNKNOWN due to known model
        // misspecification (adapter artifact / composition bias) but empirical CT5 signal
        // clearly shows ancient damage and SS-specific channels (GA0, CT3) are flat,
        // classify as DS. Uses raw per-position rates (not BIC-fitted amplitudes) to
        // avoid circularity. Effect-size thresholds (not z-scores) for depth robustness.
        if (profile.library_type == SampleDamageProfile::LibraryType::UNKNOWN) {
            const bool misspec = profile.composition_bias_5prime
                              || profile.position_0_artifact_5prime
                              || (profile.fit_offset_5prime > 1 && profile.position_0_artifact_5prime);
            // BUG FIX: rescue ran before d_max_5prime was assigned (that
            // happens further below at the raw_d_max_5prime block). Use
            // max_damage_5prime which is set up-stream at line ~1994.
            if (misspec && profile.max_damage_5prime >= 0.03f) {
                // CT5 empirical excess: max T/(T+C) at pos 1-4 over baseline
                float ct5_exc = 0.0f;
                float n_ct5   = 1.0f;
                for (int p = 1; p <= 4; ++p) {
                    double ntc = profile.tc_total_5prime[p];
                    if (ntc < 100.0) continue;
                    float exc = static_cast<float>(profile.t_freq_5prime[p] / ntc)
                              - static_cast<float>(baseline_tc);
                    if (exc > ct5_exc) { ct5_exc = exc; n_ct5 = static_cast<float>(ntc); }
                }
                // GA0 empirical excess at 3' pos-0 (SS complement-orientation indicator)
                float ga0_exc = 0.0f;
                if (profile.ag_total_3prime[0] >= 100.0) {
                    float ga0_rate = static_cast<float>(
                        profile.a_freq_3prime[0] / profile.ag_total_3prime[0]);
                    ga0_exc = std::max(0.0f, ga0_rate - static_cast<float>(baseline_ag));
                }
                // CT3 empirical excess at 3' pos 1-4 (SS original-orientation indicator)
                float ct3_exc = 0.0f;
                for (int p = 1; p <= 4; ++p) {
                    double ntc3 = profile.tc_total_3prime[p];
                    if (ntc3 < 100.0) continue;
                    float exc = static_cast<float>(ctrl_freq_3p[p])
                              - static_cast<float>(baseline_tc);
                    if (exc > ct3_exc) ct3_exc = exc;
                }
                // Overdispersed lower-95 CI for ct5_exc: SE × 2× inflation factor.
                // Avoids z-score dependence (which is depth-driven at 10M+ reads).
                float p_hat    = ct5_exc + static_cast<float>(baseline_tc);
                float se_od    = std::sqrt(p_hat * (1.0f - p_hat) / n_ct5) * 2.0f;
                float lower95  = ct5_exc - 1.96f * se_od;

                if (ct5_exc  >= 0.025f
                    && lower95  >= 0.01f
                    && ga0_exc  <= 0.01f
                    && ct3_exc  <= 0.01f
                    && std::max(ga0_exc, ct3_exc) <= 0.4f * ct5_exc) {
                    profile.library_type         = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
                    profile.library_type_rescued = true;
                }
            }
        }

        profile.library_type_auto_detected = true;

        // S1 telemetry: final winner after all post-hoc rescues / vetoes / UNKNOWN
        // overrides. final_library_bic_winner_model echoes the cascade-tournament
        // winner; override_reason is non-empty when post-hoc logic moved
        // library_type to a class that disagrees with the winner model's class.
        {
            const std::string& cw = profile.library_bic_winner_model;
            std::string cw_class;
            if (cw == "M_bias")                                    cw_class = "bias";
            else if (cw == "M_DS_symm" || cw == "M_DS_symm_art")   cw_class = "ds";
            else if (cw == "M_DS_spike")                           cw_class = profile.library_spike_is_ss ? "ss" : "ds";
            else                                                   cw_class = "ss";  // M_SS_*
            std::string lt_class;
            switch (profile.library_type) {
                case SampleDamageProfile::LibraryType::DOUBLE_STRANDED: lt_class = "ds"; break;
                case SampleDamageProfile::LibraryType::SINGLE_STRANDED: lt_class = "ss"; break;
                default:                                                lt_class = "unknown"; break;
            }
            profile.final_library_bic_winner_model = cw;
            if (cw_class == lt_class) {
                profile.final_library_bic_override_reason.clear();
            } else if (lt_class == "unknown") {
                profile.final_library_bic_override_reason =
                    profile.library_type_rescued
                        ? "post_hoc_rescue_to_unknown"
                        : "low_confidence_override_to_unknown";
            } else if (profile.library_type_rescued) {
                profile.final_library_bic_override_reason =
                    "post_hoc_rescue_" + cw_class + "_to_" + lt_class;
            } else {
                profile.final_library_bic_override_reason =
                    "post_hoc_veto_" + cw_class + "_to_" + lt_class;
            }
            // F3: post-veto final probabilities. One-hot when override fired
            // (veto/rescue/UNKNOWN), mirror raw probs when class survived.
            if (profile.final_library_bic_override_reason.empty()) {
                profile.library_p_ds_final     = profile.library_p_ds;
                profile.library_p_ss_final     = profile.library_p_ss;
                profile.library_p_bias_final   = profile.library_p_bias;
                profile.library_p_winner_final = profile.library_p_winner;
            } else {
                profile.library_p_ds_final     = (lt_class == "ds")   ? 1.0f : 0.0f;
                profile.library_p_ss_final     = (lt_class == "ss")   ? 1.0f : 0.0f;
                profile.library_p_bias_final   = (lt_class == "bias") ? 1.0f : 0.0f;
                profile.library_p_winner_final = (lt_class == "unknown") ? 0.0f : 1.0f;
            }
        }

        // P0-1: capture the BIC tournament's verdict regardless of any forced override.
        profile.library_auto_type      = profile.library_type;
        profile.library_auto_evaluable = profile.library_type_evaluable;

        // P0-1: forced-library override — applied AFTER the tournament so all
        // diagnostic BICs / posteriors remain populated. Only library_type and
        // the auto-detect flag change; library_auto_type preserves the auto call.
        profile.library_forced_type = profile.forced_library_type;
        if (profile.forced_library_type != SampleDamageProfile::LibraryType::UNKNOWN) {
            profile.library_type               = profile.forced_library_type;
            profile.library_type_auto_detected = false;
        }
    }

    float damage_signal = (profile.max_damage_5prime + profile.max_damage_3prime) / 2.0f;

    float hexamer_boost = 0.0f;
    if (profile.hexamer_damage_llr > 0.02f && !profile.terminal_inversion) {
        // Normal sample with clear hexamer signal
        hexamer_boost = profile.hexamer_damage_llr * 8.0f;
    } else if (profile.terminal_inversion) {
        // Inverted sample: z-score asymmetry distinguishes real damage from composition bias.
        // Real damage: 3' G→A relatively stronger (z_ratio < 1.2).
        // Composition bias: 5' T dominates (z_ratio > 1.5).
        float z5_abs = std::abs(profile.terminal_z_5prime);
        float z3_abs = std::abs(profile.terminal_z_3prime);
        float z_ratio = (z3_abs > 0) ? z5_abs / z3_abs : 10.0f;

        // Use absolute hexamer LLR magnitude for inverted samples
        float abs_llr = std::abs(profile.hexamer_damage_llr);

        if (z_ratio < 1.2f && abs_llr > 0.02f) {
            hexamer_boost = abs_llr * 8.0f;
        } else if (z_ratio > 1.5f) {
            // 5' signal dominates → likely composition bias
        } else {
            hexamer_boost = abs_llr * 4.0f;
        }

        if (hexamer_boost > 0.01f && damage_signal < 0.01f) {
            damage_signal = hexamer_boost;
            profile.max_damage_5prime = hexamer_boost;
            profile.max_damage_3prime = hexamer_boost;
        }
    }

    float cpg_boost = 0.0f;
    float wobble_boost = 0.0f;

    if (damage_signal > 0.01f) {
        if (profile.cpg_damage_rate > profile.non_cpg_damage_rate + 0.05f) {
            cpg_boost = 0.03f;
        }

        float wobble_enrichment = profile.codon_pos_t_rate_5prime[2] -
                                  (profile.codon_pos_t_rate_5prime[0] + profile.codon_pos_t_rate_5prime[1]) / 2.0f;
        if (wobble_enrichment > 0.05f) {
            wobble_boost = 0.02f;
        }
    }

    float total_signal = damage_signal + cpg_boost + wobble_boost + hexamer_boost;

    if (total_signal > 0.12f) {
        profile.sample_damage_prob = 0.95f;
    } else if (total_signal > 0.06f) {
        profile.sample_damage_prob = 0.80f;
    } else if (total_signal > 0.03f) {
        profile.sample_damage_prob = 0.50f;
    } else {
        profile.sample_damage_prob = 0.20f;
    }

    // D_max: joint evidence from Channel A (nucleotide frequencies) and Channel B (stop codon conversion).
    // Channel B is the independent validator: T/(T+C) elevation can come from composition OR damage,
    // but stop conversions in CAA/CAG/CGA contexts can only come from real C→T damage.
    {
        // Scan positions 1-5 for peak damage rate when adapter artifact or inversion detected.
        // Using damage_rate[0] fails when pos 0 carries adapter artifact (below baseline);
        // using damage_rate[fit_offset] fails when BIC-best offset points to a below-baseline pos.
        float raw_d_max_5prime, raw_d_max_3prime;
        if (profile.position_0_artifact_5prime || profile.inverted_pattern_5prime
            || profile.briggs_pos0_masked_5prime) {
            float peak = 0.0f;
            for (int p = 1; p <= 5 && p < 15; ++p) {
                if (profile.damage_rate_5prime[p] > peak) peak = profile.damage_rate_5prime[p];
            }
            raw_d_max_5prime = peak;
        } else {
            raw_d_max_5prime = std::clamp(profile.damage_rate_5prime[0], 0.0f, 1.0f);
        }
        if (profile.position_0_artifact_3prime || profile.inverted_pattern_3prime
            || profile.briggs_pos0_masked_3prime) {
            float peak = 0.0f;
            for (int p = 1; p <= 5 && p < 15; ++p) {
                if (profile.damage_rate_3prime[p] > peak) peak = profile.damage_rate_3prime[p];
            }
            raw_d_max_3prime = peak;
        } else {
            raw_d_max_3prime = std::clamp(profile.damage_rate_3prime[0], 0.0f, 1.0f);
        }


        float d_sum = raw_d_max_5prime + raw_d_max_3prime;
        if (d_sum > 0.01f) {
            profile.asymmetry = std::abs(raw_d_max_5prime - raw_d_max_3prime) / (d_sum / 2.0f);
        } else {
            profile.asymmetry = 0.0f;
        }
        profile.high_asymmetry = (profile.asymmetry > 0.5f);

        {
            const uint64_t MIN_C_SITES = 10000;  // Minimum C sites for valid per-bin estimate
            float weighted_sum = 0.0f;
            float weight_sum = 0.0f;
            float peak_dmax = 0.0f;
            int peak_bin = -1;

            for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
                auto& b = profile.gc_bins[bin];

                b.c_sites = b.c_interior;
                for (int p = 0; p < 15; ++p) {
                    b.c_sites += b.c_counts[p];
                }

                if (b.n_reads < 1000 || b.c_sites < MIN_C_SITES) {
                    continue;  // Skip bins with insufficient data
                }

                double t_baseline = static_cast<double>(b.t_interior) /
                                   (b.t_interior + b.c_interior + 1);
                double c_baseline = 1.0 - t_baseline;
                double t_terminal = static_cast<double>(b.t_counts[0]) /
                                   (b.t_counts[0] + b.c_counts[0] + 1);

                // For SS libraries the damage signal also appears as C→T at the 3' end;
                // take the max of the two terminal rates to avoid losing signal when the
                // 3' end is dominant (as on typical ssDNA ancient libraries).
                const bool is_ss_bin = (profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) ||
                                       (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
                if (is_ss_bin) {
                    double t_terminal_3 = static_cast<double>(b.t_counts_3prime[0]) /
                                          (b.t_counts_3prime[0] + b.c_counts_3prime[0] + 1);
                    t_terminal = std::max(t_terminal, t_terminal_3);
                }

                if (c_baseline > 0.1) {
                    b.d_max = std::clamp(static_cast<float>((t_terminal - t_baseline) / c_baseline),
                                         0.0f, 1.0f);
                }

                double stop_baseline = static_cast<double>(b.stop_interior) /
                                      (b.stop_interior + b.pre_interior + 1);
                double stop_terminal = static_cast<double>(b.stop_counts[0]) /
                                      (b.stop_counts[0] + b.pre_counts[0] + 1);

                if (stop_baseline < 0.99) {
                    b.d_max_channel_b = std::clamp(
                        static_cast<float>((stop_terminal - stop_baseline) / (1.0 - stop_baseline)),
                        0.0f, 1.0f);
                }

                b.valid = true;

                float bin_dmax = std::max(b.d_max, b.d_max_channel_b);
                float weight = static_cast<float>(b.c_sites);
                weighted_sum += bin_dmax * weight;
                weight_sum += weight;

                if (bin_dmax > peak_dmax) {
                    peak_dmax = bin_dmax;
                    peak_bin = bin;
                }

                int gc_low = bin * 10;
                int gc_high = gc_low + 10;
            }

            if (weight_sum > 0) {
                profile.gc_stratified_d_max_weighted = weighted_sum / weight_sum;
                profile.gc_stratified_d_max_peak = peak_dmax;
                profile.gc_peak_bin = peak_bin;
                profile.gc_stratified_valid = true;


                {
                    constexpr float LLR_THRESHOLD = 10.0f;
                    constexpr float MIN_DMAX_THRESHOLD = 0.01f;

                    uint64_t total_obs = 0;
                    uint64_t damaged_obs = 0;
                    float damaged_weighted_d = 0.0f;
                    uint64_t damaged_weight = 0;
                    int n_damaged = 0;

                    for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
                        auto& b = profile.gc_bins[bin];
                        if (!b.valid) continue;

                        double baseline_tc = static_cast<double>(b.t_interior) /
                                            std::max(1.0, static_cast<double>(b.t_interior + b.c_interior));
                        b.baseline_tc = static_cast<float>(std::clamp(baseline_tc, 0.01, 0.99));

                        uint64_t n_obs = b.n_terminal_obs();
                        total_obs += n_obs;

                        float ll_damaged = 0.0f;
                        float ll_undamaged = 0.0f;
                        float lambda = profile.lambda_5prime;
                        const bool is_ss_bin_llr =
                            (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);

                        for (int p = 0; p < 15; ++p) {
                            float decay = std::exp(-lambda * p);
                            float delta_p = b.d_max * decay;
                            float pi_undamaged = b.baseline_tc;
                            float pi_damaged = b.baseline_tc + (1.0f - b.baseline_tc) * delta_p;

                            pi_undamaged = std::clamp(pi_undamaged, 0.001f, 0.999f);
                            pi_damaged = std::clamp(pi_damaged, 0.001f, 0.999f);

                            double k = static_cast<double>(b.t_counts[p]);
                            double n = static_cast<double>(b.t_counts[p] + b.c_counts[p]);
                            if (is_ss_bin_llr && p < static_cast<int>(b.t_counts_3prime.size())) {
                                k += static_cast<double>(b.t_counts_3prime[p]);
                                n += static_cast<double>(b.t_counts_3prime[p] + b.c_counts_3prime[p]);
                            }
                            if (n > 0) {
                                ll_damaged += static_cast<float>(k * std::log(pi_damaged) + (n - k) * std::log(1.0f - pi_damaged));
                                ll_undamaged += static_cast<float>(k * std::log(pi_undamaged) + (n - k) * std::log(1.0f - pi_undamaged));
                            }
                        }

                        b.llr = ll_damaged - ll_undamaged;

                        b.classified_damaged = (b.llr > LLR_THRESHOLD) && (b.d_max > MIN_DMAX_THRESHOLD);

                        // Soft probability via logistic on (LLR - threshold)
                        float llr_centered = b.llr - LLR_THRESHOLD;
                        b.p_damaged = 1.0f / (1.0f + std::exp(-0.5f * llr_centered));

                        if (b.d_max < MIN_DMAX_THRESHOLD) {
                            b.p_damaged = 0.0f;
                        }

                        if (b.classified_damaged) {
                            damaged_obs += n_obs;
                            damaged_weighted_d += b.d_max * static_cast<float>(n_obs);
                            damaged_weight += n_obs;
                            ++n_damaged;
                        }
                    }

                    if (total_obs > 0) {
                        profile.pi_damaged = static_cast<float>(damaged_obs) / static_cast<float>(total_obs);
                    }
                    if (damaged_weight > 0) {
                        profile.d_ancient = damaged_weighted_d / static_cast<float>(damaged_weight);
                    }

                    float pop_weighted_d = 0.0f;
                    for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
                        const auto& b = profile.gc_bins[bin];
                        if (!b.valid) continue;
                        pop_weighted_d += b.d_max * static_cast<float>(b.n_terminal_obs());
                    }
                    if (total_obs > 0) {
                        profile.d_population = pop_weighted_d / static_cast<float>(total_obs);
                    }

                    profile.n_damaged_bins = n_damaged;

                }

                // For SS libraries, damage appears as C→T at both ends; feed the
                // combined terminal counts into the GC mixture fit so that bins
                // with 3'-dominant damage (the typical ssDNA case) contribute.
                const bool is_ss_super = (profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) ||
                                         (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
                std::array<SuperRead, N_GC_BINS> super_reads;
                for (int bin = 0; bin < N_GC_BINS; ++bin) {
                    const auto& b = profile.gc_bins[bin];
                    super_reads[bin].gc_bin = bin;
                    super_reads[bin].c_sites = static_cast<double>(b.c_sites);
                    if (is_ss_super) {
                        super_reads[bin].c_sites +=
                            static_cast<double>(b.c_counts_3prime[0] + b.c_counts_3prime[1] + b.c_counts_3prime[2]);
                    }

                    for (int p = 0; p < N_POSITIONS; ++p) {
                        double k = static_cast<double>(b.t_counts[p]);
                        double n = static_cast<double>(b.t_counts[p] + b.c_counts[p]);
                        if (is_ss_super) {
                            k += static_cast<double>(b.t_counts_3prime[p]);
                            n += static_cast<double>(b.t_counts_3prime[p] + b.c_counts_3prime[p]);
                        }
                        super_reads[bin].k_tc[p] = k;
                        super_reads[bin].n_tc[p] = n;
                    }

                    for (int p = 0; p < N_POSITIONS; ++p) {
                        super_reads[bin].k_ag[p] = static_cast<double>(b.a_counts[p]);
                        super_reads[bin].n_ag[p] = static_cast<double>(b.a_counts[p] + b.g_counts[p]);
                    }

                    for (int p = 0; p < N_POSITIONS; ++p) {
                        super_reads[bin].k_stop[p] = static_cast<double>(b.stop_counts[p]);
                        super_reads[bin].n_stop[p] = static_cast<double>(b.stop_counts[p] + b.pre_counts[p]);
                    }

                    super_reads[bin].k_tc_int = static_cast<double>(b.t_interior);
                    super_reads[bin].n_tc_int = static_cast<double>(b.t_interior + b.c_interior);
                    super_reads[bin].k_stop_int = static_cast<double>(b.stop_interior);
                    super_reads[bin].n_stop_int = static_cast<double>(b.stop_interior + b.pre_interior);

                    super_reads[bin].k_ag_int = static_cast<double>(b.a_interior);
                    super_reads[bin].n_ag_int = static_cast<double>(b.a_interior + b.g_interior);
                }

                auto mixture_result = MixtureDamageModel::fit(super_reads);
                profile.mixture_K = mixture_result.K;
                profile.mixture_n_components = mixture_result.n_components;
                profile.mixture_d_population = mixture_result.d_population;
                profile.mixture_d_ancient = mixture_result.d_ancient;
                profile.mixture_d_reference = mixture_result.d_reference;
                profile.mixture_pi_ancient = mixture_result.pi_ancient;
                profile.mixture_bic = mixture_result.bic;
                profile.mixture_converged = mixture_result.converged;
                profile.mixture_identifiable = mixture_result.identifiable;

            }
        }

        profile.d_max_5prime = raw_d_max_5prime;
        profile.d_max_3prime = raw_d_max_3prime;

        if (profile.damage_artifact) {
            profile.d_max_5prime = 0.0f;
            profile.d_max_3prime = 0.0f;
            profile.d_max_combined = 0.0f;
            profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
        } else if (profile.damage_validated) {
            bool channel_a_unreliable = (profile.inverted_pattern_5prime && profile.inverted_pattern_3prime)
                                      || profile.position_0_artifact_5prime
                                      || profile.position_0_artifact_3prime;

            if (channel_a_unreliable && profile.channel_b_quantifiable && profile.d_max_from_channel_b > 0.01f) {
                profile.d_max_combined = profile.d_max_from_channel_b;
                profile.d_max_source = SampleDamageProfile::DmaxSource::CHANNEL_B_STRUCTURAL;
            } else if (profile.inverted_pattern_3prime && !profile.inverted_pattern_5prime) {
                profile.d_max_combined = raw_d_max_5prime;
                profile.d_max_source = SampleDamageProfile::DmaxSource::FIVE_PRIME_ONLY;
            } else if (profile.inverted_pattern_5prime && !profile.inverted_pattern_3prime) {
                profile.d_max_combined = raw_d_max_3prime;
                profile.d_max_source = SampleDamageProfile::DmaxSource::THREE_PRIME_ONLY;
            } else if (profile.mixture_converged && profile.mixture_d_reference > 0.01f) {
                profile.d_max_combined = profile.mixture_d_reference;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            } else if (profile.gc_stratified_valid) {
                profile.d_max_combined = profile.gc_stratified_d_max_weighted;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            } else {
                profile.d_max_combined = profile.joint_delta_max;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            }
        } else if (profile.joint_model_valid && profile.joint_p_damage > 0.5f) {
            if (profile.mixture_converged && profile.mixture_d_reference > 0.01f) {
                profile.d_max_combined = profile.mixture_d_reference;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            } else if (profile.gc_stratified_valid) {
                profile.d_max_combined = profile.gc_stratified_d_max_weighted;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            } else {
                profile.d_max_combined = profile.joint_delta_max;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            }
        } else if (profile.joint_model_valid) {
            // For single-stranded libraries Channel B is not applicable (ss damage is G→A at 3'),
            // so fall back to Channel A rather than zeroing d_max.
            const bool is_ss = (profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) ||
                               (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
            if (is_ss) {
                profile.d_max_5prime = raw_d_max_5prime;
                profile.d_max_3prime = raw_d_max_3prime;
                profile.d_max_combined = std::max(raw_d_max_5prime, raw_d_max_3prime);
                profile.d_max_source = SampleDamageProfile::DmaxSource::MAX_SS_ASYMMETRY;
            } else {
                // ds libraries: avoid Channel A's compositional false positives.
                // If the GC-stratified mixture is identifiable and finds a real
                // ancient subpopulation, prefer its estimate over zeroing.
                if (profile.mixture_converged && profile.mixture_identifiable &&
                    profile.mixture_d_ancient > 0.02f) {
                    profile.d_max_5prime = raw_d_max_5prime;
                    profile.d_max_3prime = raw_d_max_3prime;
                    profile.d_max_combined = profile.mixture_d_population;
                    profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
                } else {
                    profile.d_max_5prime = 0.0f;
                    profile.d_max_3prime = 0.0f;
                    profile.d_max_combined = 0.0f;
                    profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
                }
            }
        } else if (profile.channel_b_quantifiable) {
            profile.d_max_5prime = raw_d_max_5prime;
            profile.d_max_3prime = raw_d_max_3prime;
            profile.d_max_combined = profile.d_max_from_channel_b;
            profile.d_max_source = SampleDamageProfile::DmaxSource::CHANNEL_B_STRUCTURAL;
        } else {
            profile.d_max_5prime = raw_d_max_5prime;
            profile.d_max_3prime = raw_d_max_3prime;

            if (profile.high_asymmetry) {
                const bool is_ss = (profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) ||
                                   (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
                if (is_ss) {
                    // ss libraries have asymmetric damage by design; use max to capture the damaged end
                    profile.d_max_combined = std::max(profile.d_max_5prime, profile.d_max_3prime);
                    profile.d_max_source = SampleDamageProfile::DmaxSource::MAX_SS_ASYMMETRY;
                } else {
                    // ds libraries: high asymmetry suggests artifact - use conservative min
                    profile.d_max_combined = std::min(profile.d_max_5prime, profile.d_max_3prime);
                    profile.d_max_source = SampleDamageProfile::DmaxSource::MIN_ASYMMETRY;
                }
            } else {
                profile.d_max_combined = (profile.d_max_5prime + profile.d_max_3prime) / 2.0f;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            }
        }

        // Skip inversion fallback when Channel B has validated damage (Channel B is ground truth)
        if (!profile.damage_validated) {
            if (profile.inverted_pattern_5prime && !profile.inverted_pattern_3prime) {
                profile.d_max_combined = profile.d_max_3prime;
                profile.d_max_source = SampleDamageProfile::DmaxSource::THREE_PRIME_ONLY;
            } else if (profile.inverted_pattern_3prime && !profile.inverted_pattern_5prime) {
                profile.d_max_combined = profile.d_max_5prime;
                profile.d_max_source = SampleDamageProfile::DmaxSource::FIVE_PRIME_ONLY;
            } else if (profile.inverted_pattern_5prime && profile.inverted_pattern_3prime) {
                // Both ends inverted: normally zero everything (no reliable Channel A signal).
                // Exception: if adapter offset was detected on either end, the scan-for-peak
                // already found the biological damage in pos 1-5. Preserve those values.
                bool has_adapter_5 = profile.position_0_artifact_5prime;
                bool has_adapter_3 = profile.position_0_artifact_3prime;
                if (!has_adapter_5 && !has_adapter_3) {
                    profile.d_max_5prime = 0.0f;
                    profile.d_max_3prime = 0.0f;
                    profile.d_max_combined = 0.0f;
                    profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
                }
            }
        }

    }

    // d_metamatch: Channel B-anchored estimate.
    // Formula: d_metamatch = d_global + γ × (d_channel_b - d_global)
    // γ from Channel B LLR; blends toward Channel B for validated samples.
    {
        double weighted_sum = 0.0;
        double weight_sum = 0.0;

        for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
            const auto& b = profile.gc_bins[bin];
            if (!b.valid || b.n_reads < 100) continue;

            float gc_center = (bin + 0.5f) / static_cast<float>(SampleDamageProfile::N_GC_BINS);
            float gc_deviation = gc_center - 0.50f;
            float alignability_weight = std::exp(-gc_deviation * gc_deviation / (2.0f * 0.15f * 0.15f));
            double read_weight = std::sqrt(static_cast<double>(b.n_reads));
            double total_weight = alignability_weight * read_weight;

            weighted_sum += b.d_max * total_weight;
            weight_sum += total_weight;
        }

        if (weight_sum > 0) {
            profile.d_alignability_weighted = static_cast<float>(weighted_sum / weight_sum);
        } else {
            profile.d_alignability_weighted = profile.d_max_combined;
        }

        float channel_b_llr = profile.stop_decay_llr_5prime;
        // γ: sigmoid on LLR/50000 (LLR=0 → γ=0.5, LLR=100000 → γ≈0.99)
        float gamma_raw = 1.0f / (1.0f + std::exp(-channel_b_llr / 50000.0f));

        float d_global = profile.d_max_combined;
        float d_channel_b = profile.d_max_from_channel_b;

        if (!profile.damage_validated || profile.damage_artifact) {
            profile.metamatch_gamma = 0.0f;
            profile.d_metamatch = d_global;
        } else if (profile.channel_b_quantifiable && d_channel_b > 0.01f) {
            // Asymmetric: stronger pull toward Channel B when it's higher than d_global,
            // weaker pull when it's lower (avoid under-estimation).
            if (d_channel_b > d_global) {
                profile.metamatch_gamma = gamma_raw;
            } else {
                profile.metamatch_gamma = 0.3f * gamma_raw;
            }
            profile.d_metamatch = d_global + profile.metamatch_gamma * (d_channel_b - d_global);
        } else {
            profile.metamatch_gamma = 0.5f * gamma_raw;
            profile.d_metamatch = d_global + profile.metamatch_gamma * (profile.d_alignability_weighted - d_global);
        }

        profile.d_metamatch = std::clamp(profile.d_metamatch, 0.0f, 1.0f);

        double alignability_total = 0.0;
        double n_total = 0.0;
        for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
            const auto& b = profile.gc_bins[bin];
            if (b.n_reads == 0) continue;

            float gc_center = (bin + 0.5f) / static_cast<float>(SampleDamageProfile::N_GC_BINS);
            float gc_deviation = gc_center - 0.50f;
            float alignability = std::exp(-gc_deviation * gc_deviation / (2.0f * 0.15f * 0.15f));

            alignability_total += alignability * b.n_reads;
            n_total += b.n_reads;
        }
        profile.mean_alignability = (n_total > 0) ? static_cast<float>(alignability_total / n_total) : 0.5f;
    }

    // Damage status: effect-size based, independent of library type.
    // Must run after d_max_5prime/d_max_3prime are finalized (line ~2356+).
    // t_freq_5prime[p] has been normalized to rate at line ~990, so use directly.
    // tc_total_5prime[p] still holds raw T+C counts, used for SE calculation.
    // Uses overdispersed CI (2× SE inflation) so depth alone doesn't drive PRESENT.
    {
        float dmax = std::max(profile.d_max_5prime, profile.d_max_3prime);
        if (dmax >= 0.02f) {
            auto lower95_for_end = [&](const std::array<double,15>& freq,
                                       const std::array<double,15>& total,
                                       double baseline) -> float {
                float exc_max = 0.0f; float n_used = 1.0f;
                for (int p = 1; p <= 4; ++p) {
                    double n = total[p];
                    if (n < 100.0) continue;
                    float exc = static_cast<float>(freq[p] - baseline);
                    if (exc > exc_max) { exc_max = exc; n_used = static_cast<float>(n); }
                }
                if (exc_max <= 0.0f) return -1.0f;
                float p_hat = exc_max + static_cast<float>(baseline);
                float se_od = std::sqrt(p_hat * (1.0f - p_hat) / n_used) * 2.0f;
                return exc_max - 1.96f * se_od;
            };
            float lb5 = lower95_for_end(profile.t_freq_5prime, profile.tc_total_5prime, baseline_tc);
            // SS 3' damage is C→T (same channel as 5'); DS 3' damage is G→A.
            // t_freq_3prime holds raw counts so compute rate explicitly for SS.
            float lb3;
            if (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) {
                float exc_max = 0.0f; float n_used = 1.0f;
                for (int p = 1; p <= 4; ++p) {
                    double n = profile.tc_total_3prime[p];
                    if (n < 100.0) continue;
                    float exc = static_cast<float>(
                        (n > 0 ? profile.t_freq_3prime[p] / n : 0.0) - baseline_tc);
                    if (exc > exc_max) { exc_max = exc; n_used = static_cast<float>(n); }
                }
                if (exc_max <= 0.0f) {
                    lb3 = -1.0f;
                } else {
                    float p_hat = exc_max + static_cast<float>(baseline_tc);
                    float se_od = std::sqrt(p_hat * (1.0f - p_hat) / n_used) * 2.0f;
                    lb3 = exc_max - 1.96f * se_od;
                }
            } else {
                lb3 = lower95_for_end(profile.a_freq_3prime, profile.ag_total_3prime, baseline_ag);
            }
            float lower95 = std::max(lb5, lb3);
            profile.damage_status = (lower95 >= 0.01f)
                ? SampleDamageProfile::DamageStatus::PRESENT
                : SampleDamageProfile::DamageStatus::WEAK;
        } else {
            profile.damage_status = SampleDamageProfile::DamageStatus::ABSENT;
            // Inversion artifacts zero d_max_5/3prime even when BIC fit ct5~=ga3
            // up to ~0.17. Preserve confident BIC verdicts; only fall back to
            // UNKNOWN when the tournament itself was uncertain.
            bool bic_confident =
                profile.library_p_winner >= 0.95f &&
                profile.library_auto_type != SampleDamageProfile::LibraryType::UNKNOWN;
            if (profile.forced_library_type == SampleDamageProfile::LibraryType::UNKNOWN
                && !bic_confident) {
                profile.library_type = SampleDamageProfile::LibraryType::UNKNOWN;
            }
        }
    }

    // --- Preservation index (evidence × reliability) ---
    {
        static constexpr float EPS = 1e-6f;
        auto sig = [](float x) { return 1.0f / (1.0f + std::exp(-x)); };
        bool is_ss = (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);

        // f5: 5' terminal C→T; half-point at 5% damage (Briggs 2007)
        profile.preservation_f5 = sig((profile.d_max_5prime - 0.05f) / 0.04f);

        // f3: 3' terminal signal; half-point at 5% damage.
        // For ds libraries where d_max_3prime=0 despite real 5' damage, the terminal
        // base is often clipped by the same trimmer artifact that shifts the 5' peak to
        // pos1 — G→A is inward-displaced and absorbed into background. Impute f3 from
        // the 5' signal: ds deamination is symmetric, so f5 is the calibrated estimate.
        bool ds_3prime_censored = !is_ss
            && profile.d_max_5prime > 0.05f
            && profile.d_max_3prime < 0.03f;
        profile.preservation_f3 = ds_3prime_censored
            ? profile.preservation_f5
            : sig((profile.d_max_3prime - 0.05f) / 0.04f);

        // f_coh: mixture coherence — ancient subpop identifiable with real damage
        {
            float mix_signal = 0.0f;
            float trust = (profile.mixture_identifiable && profile.mixture_converged) ? 1.0f
                        : (profile.mixture_identifiable || profile.mixture_converged)  ? 0.5f
                        : 0.1f;
            mix_signal = trust * profile.mixture_pi_ancient
                       * sig((profile.mixture_d_ancient - 0.05f) / 0.03f);
            // ds: penalise end-asymmetry, but use d_max_combined as floor so that
            // samples where bulk estimator zeroes d3 (noise/rescue limitation) are not
            // wrongly collapsed to near-zero symmetry.
            float d5_sym = std::max(profile.d_max_5prime, profile.d_max_combined);
            // When 3' is censored (trimmer artifact), treat symmetry as perfect —
            // the signal is present but displaced inward, not absent.
            float d3_sym = ds_3prime_censored ? profile.d_max_5prime
                         : std::max(profile.d_max_3prime, profile.d_max_combined);
            float sym = (is_ss || ds_3prime_censored) ? 1.0f
                : std::exp(-std::abs(std::log((d5_sym + EPS)
                                             / (d3_sym + EPS))) / 0.7f);
            profile.preservation_f_coh = std::sqrt(std::max(mix_signal, EPS) * std::max(sym, EPS));
        }

        // f_cpg: CpG age-bias; log2(CpG_dmax/nonCpG_dmax) > 1 = methyl-C deamination enriched.
        // Near-saturation (d_max_combined > 0.20) both contexts saturate and the ratio collapses
        // toward zero — signal is uninformative, fall back to neutral prior.
        if (!std::isnan(profile.log2_cpg_ratio) && profile.dmax_ct5_noncpg_like > 0.01f
                && profile.d_max_combined < 0.20f) {
            // Floor at uninformative prior: absent/inverted CpG signal is not evidence
            // against antiquity (ss libraries, saturated damage, CpG-depleted taxa).
            profile.preservation_f_cpg = std::max(sig((profile.log2_cpg_ratio - 1.0f) / 0.6f), 0.3f);
        } else {
            profile.preservation_f_cpg = 0.3f;  // uninformative prior
        }

        // Weighted geometric mean of 4 factors
        float w5, w3, w_coh, w_cpg;
        if (is_ss) { w5 = 0.35f; w3 = 0.20f; w_coh = 0.28f; w_cpg = 0.17f; }
        else        { w5 = 0.27f; w3 = 0.27f; w_coh = 0.28f; w_cpg = 0.18f; }

        profile.preservation_evidence = std::exp(
            w5   * std::log(std::max(profile.preservation_f5,    EPS)) +
            w3   * std::log(std::max(profile.preservation_f3,    EPS)) +
            w_coh* std::log(std::max(profile.preservation_f_coh, EPS)) +
            w_cpg* std::log(std::max(profile.preservation_f_cpg, EPS)));

        // Reliability gates (continuous — no hard cliffs)
        float g_N   = sig((std::log10(static_cast<float>(profile.n_reads) + 1.0f) - 2.7f) / 0.35f);
        // When pi_ancient > 0.90 and EM converged, the mixture is essentially pure-ancient:
        // identifiability criterion doesn't apply (no modern class to separate from).
        bool effectively_pure_ancient = profile.mixture_converged && profile.mixture_pi_ancient > 0.90f;
        float g_fit = (effectively_pure_ancient ||
                       (profile.mixture_identifiable && profile.mixture_converged)) ? 1.0f
                    : (profile.mixture_identifiable || profile.mixture_converged)    ? 0.5f
                    : 0.15f;
        float g_ox  = 1.0f;
        if (profile.ox_is_artifact)
            g_ox = (profile.ox_d_max >= profile.d_max_combined) ? 0.1f : 0.5f;

        profile.preservation_reliability = g_N * g_fit * g_ox;
        profile.preservation_score = std::clamp(
            profile.preservation_evidence * profile.preservation_reliability, 0.0f, 1.0f);

        // Categorical label
        using PL = SampleDamageProfile::PreservationLabel;
        if (g_N < 0.3f)
            profile.preservation_label = PL::INSUFFICIENT;
        else if (profile.ox_is_artifact)
            profile.preservation_label = PL::ARTIFACT_SUSPECTED;
        else if (profile.preservation_score < 0.15f)
            profile.preservation_label = PL::MODERN_LIKE;
        else if (profile.preservation_score < 0.35f)
            profile.preservation_label = PL::WEAK;
        else if (profile.preservation_score < 0.60f)
            profile.preservation_label = PL::MODERATE;
        else if (profile.preservation_score < 0.80f)
            profile.preservation_label = PL::STRONG;
        else
            profile.preservation_label = PL::EXCEPTIONAL;
    }
    // Lifecycle: mark finalized so re-entry is rejected (see SampleDamageProfile::finalized).
    profile.finalized = true;
}

void FrameSelector::merge_sample_profiles(SampleDamageProfile& dst, const SampleDamageProfile& src) {
    // Lifecycle guard: merge sums raw count arrays. After finalize they hold
    // rates, so summing them produces nonsense. Reject explicitly rather than
    // corrupt downstream silently.
    if (dst.finalized || src.finalized) {
        throw std::logic_error(
            "merge_sample_profiles: profiles must be merged BEFORE "
            "finalize_sample_profile() is called on either side.");
    }
    for (int i = 0; i < 15; ++i) {
        dst.t_freq_5prime[i] += src.t_freq_5prime[i];
        dst.c_freq_5prime[i] += src.c_freq_5prime[i];
        dst.a_freq_3prime[i] += src.a_freq_3prime[i];
        dst.g_freq_3prime[i] += src.g_freq_3prime[i];
        dst.tc_total_5prime[i] += src.tc_total_5prime[i];
        dst.ag_total_3prime[i] += src.ag_total_3prime[i];
        dst.a_freq_5prime[i] += src.a_freq_5prime[i];
        dst.g_freq_5prime[i] += src.g_freq_5prime[i];
        dst.t_freq_3prime[i] += src.t_freq_3prime[i];
        dst.c_freq_3prime[i] += src.c_freq_3prime[i];
        dst.tc_total_3prime[i] += src.tc_total_3prime[i];
    }
    for (int i = 0; i < SampleDamageProfile::BG_TAIL_N; ++i) {
        dst.tail_t_5prime[i]  += src.tail_t_5prime[i];
        dst.tail_tc_5prime[i] += src.tail_tc_5prime[i];
        dst.tail_a_3prime[i]  += src.tail_a_3prime[i];
        dst.tail_ag_3prime[i] += src.tail_ag_3prime[i];
    }
    dst.pe_short_insert_skipped += src.pe_short_insert_skipped;
    if (src.input_mode == SampleDamageProfile::InputMode::PAIRED) {
        dst.input_mode = SampleDamageProfile::InputMode::PAIRED;
    }

    dst.baseline_t_freq += src.baseline_t_freq;
    dst.baseline_c_freq += src.baseline_c_freq;
    dst.baseline_a_freq += src.baseline_a_freq;
    dst.baseline_g_freq += src.baseline_g_freq;

    for (int p = 0; p < 3; ++p) {
        dst.codon_pos_t_count_5prime[p] += src.codon_pos_t_count_5prime[p];
        dst.codon_pos_c_count_5prime[p] += src.codon_pos_c_count_5prime[p];
        dst.codon_pos_a_count_3prime[p] += src.codon_pos_a_count_3prime[p];
        dst.codon_pos_g_count_3prime[p] += src.codon_pos_g_count_3prime[p];
    }

    dst.cpg_c_count += src.cpg_c_count;
    dst.cpg_t_count += src.cpg_t_count;
    dst.non_cpg_c_count += src.non_cpg_c_count;
    dst.non_cpg_t_count += src.non_cpg_t_count;

    for (int ctx = 0; ctx < SampleDamageProfile::N_CT_CTX; ++ctx) {
        for (int p = 0; p < SampleDamageProfile::N_POS; ++p) {
            dst.ct_ctx_t_5prime[ctx][p] += src.ct_ctx_t_5prime[ctx][p];
            dst.ct_ctx_total_5prime[ctx][p] += src.ct_ctx_total_5prime[ctx][p];
        }
        dst.ct_ctx_t_interior[ctx] += src.ct_ctx_t_interior[ctx];
        dst.ct_ctx_total_interior[ctx] += src.ct_ctx_total_interior[ctx];
    }
    for (int i = 0; i < SampleDamageProfile::N_OXOG16; ++i) {
        dst.oxog16_t[i] += src.oxog16_t[i];
        dst.oxog16_a_rc[i] += src.oxog16_a_rc[i];
    }
    for (int i = 0; i < SampleDamageProfile::N_TRINUC; ++i) {
        dst.tri_5prime_terminal[i] += src.tri_5prime_terminal[i];
        dst.tri_5prime_interior[i] += src.tri_5prime_interior[i];
        dst.tri_3prime_terminal[i] += src.tri_3prime_terminal[i];
        dst.tri_3prime_interior[i] += src.tri_3prime_interior[i];
    }
    // Merge upstream-context-aware accumulators
    for (int uctx = 0; uctx < SampleDamageProfile::N_UPSTREAM_CTX; ++uctx) {
        for (int p = 0; p < SampleDamageProfile::N_POS; ++p) {
            dst.ct5_t_by_upstream[uctx][p] += src.ct5_t_by_upstream[uctx][p];
            dst.ct5_total_by_upstream[uctx][p] += src.ct5_total_by_upstream[uctx][p];
        }
        dst.ct5_t_interior_by_upstream[uctx] += src.ct5_t_interior_by_upstream[uctx];
        dst.ct5_total_interior_by_upstream[uctx] += src.ct5_total_interior_by_upstream[uctx];
    }

    {
        auto& da = dst.interior_ct_cluster;
        const auto& sa = src.interior_ct_cluster;
        da.reads_used_ct        += sa.reads_used_ct;
        da.reads_used_ag        += sa.reads_used_ag;
        da.short_reads_skipped  += sa.short_reads_skipped;
        for (int d = 1; d <= 10; ++d) {
            da.obs_ct[d]   += sa.obs_ct[d];
            da.pairs_ct[d] += sa.pairs_ct[d];
            da.exp_ct[d]   += sa.exp_ct[d];
            da.var_ct[d]   += sa.var_ct[d];
            da.obs_ag[d]   += sa.obs_ag[d];
            da.pairs_ag[d] += sa.pairs_ag[d];
            da.exp_ag[d]   += sa.exp_ag[d];
            da.var_ag[d]   += sa.var_ag[d];
        }
    }

    for (uint32_t i = 0; i < 4096; ++i) {
        dst.hexamer_count_5prime[i] += src.hexamer_count_5prime[i];
        dst.hexamer_count_interior[i] += src.hexamer_count_interior[i];
    }
    dst.n_hexamers_5prime += src.n_hexamers_5prime;
    dst.n_hexamers_interior += src.n_hexamers_interior;

    for (int i = 0; i < 15; ++i) {
        dst.convertible_caa_5prime[i] += src.convertible_caa_5prime[i];
        dst.convertible_taa_5prime[i] += src.convertible_taa_5prime[i];
        dst.convertible_cag_5prime[i] += src.convertible_cag_5prime[i];
        dst.convertible_tag_5prime[i] += src.convertible_tag_5prime[i];
        dst.convertible_cga_5prime[i] += src.convertible_cga_5prime[i];
        dst.convertible_tga_5prime[i] += src.convertible_tga_5prime[i];
        dst.total_codons_5prime[i] += src.total_codons_5prime[i];
    }
    dst.convertible_caa_interior += src.convertible_caa_interior;
    dst.convertible_taa_interior += src.convertible_taa_interior;
    dst.convertible_cag_interior += src.convertible_cag_interior;
    dst.convertible_tag_interior += src.convertible_tag_interior;
    dst.convertible_cga_interior += src.convertible_cga_interior;
    dst.convertible_tga_interior += src.convertible_tga_interior;
    dst.total_codons_interior += src.total_codons_interior;

    for (int i = 0; i < 15; ++i) {
        dst.convertible_gag_5prime[i] += src.convertible_gag_5prime[i];
        dst.convertible_tag_ox_5prime[i] += src.convertible_tag_ox_5prime[i];
        dst.convertible_gaa_5prime[i] += src.convertible_gaa_5prime[i];
        dst.convertible_taa_ox_5prime[i] += src.convertible_taa_ox_5prime[i];
        dst.convertible_gga_5prime[i] += src.convertible_gga_5prime[i];
        dst.convertible_tga_ox_5prime[i] += src.convertible_tga_ox_5prime[i];
        dst.g_count_5prime[i]    += src.g_count_5prime[i];
        dst.t_from_g_5prime[i]  += src.t_from_g_5prime[i];
        dst.c_count_ox_5prime[i]+= src.c_count_ox_5prime[i];
        dst.a_from_c_5prime[i]  += src.a_from_c_5prime[i];
    }
    dst.baseline_g_to_t_count  += src.baseline_g_to_t_count;
    dst.baseline_g_total       += src.baseline_g_total;
    dst.baseline_c_to_a_count  += src.baseline_c_to_a_count;
    dst.baseline_c_ox_total    += src.baseline_c_ox_total;

    dst.convertible_gag_interior += src.convertible_gag_interior;
    dst.convertible_tag_ox_interior += src.convertible_tag_ox_interior;
    dst.convertible_gaa_interior += src.convertible_gaa_interior;
    dst.convertible_taa_ox_interior += src.convertible_taa_ox_interior;
    dst.convertible_gga_interior += src.convertible_gga_interior;
    dst.convertible_tga_ox_interior += src.convertible_tga_ox_interior;

    for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
        auto& db = dst.gc_bins[bin];
        const auto& sb = src.gc_bins[bin];
        for (int p = 0; p < 15; ++p) {
            db.t_counts[p] += sb.t_counts[p];
            db.c_counts[p] += sb.c_counts[p];
            db.a_counts[p] += sb.a_counts[p];
            db.g_counts[p] += sb.g_counts[p];
            db.t_counts_3prime[p] += sb.t_counts_3prime[p];
            db.c_counts_3prime[p] += sb.c_counts_3prime[p];
            db.a_counts_3prime[p] += sb.a_counts_3prime[p];
            db.g_counts_3prime[p] += sb.g_counts_3prime[p];
            db.stop_counts[p] += sb.stop_counts[p];
            db.pre_counts[p] += sb.pre_counts[p];
        }
        db.t_interior += sb.t_interior;
        db.c_interior += sb.c_interior;
        db.a_interior += sb.a_interior;
        db.g_interior += sb.g_interior;
        db.stop_interior += sb.stop_interior;
        db.pre_interior += sb.pre_interior;
        db.n_reads += sb.n_reads;
    }

    dst.n_reads += src.n_reads;
}

// Paired-end variant. R1 contributes 5'-end counters + interior baseline;
// R2 (complement-mapped) contributes 3'-end counters. Read 3' ends are
// ignored — for inserts shorter than read length, R2's 5' end may read
// through into R1's 5' adapter, contaminating per_pos_3prime. The caller
// (fqdup PE worker) skips short pairs before calling this.
//
// Coverage scope: per-pos C->T and G->A counters at both ends, tail-anchored
// background counters, codon-position counters, hexamer counts at 5', and
// interior baseline. Advanced features filled by single-end update
// (CpG-like ctx, upstream ctx, oxoG 16-ctx, trinuc spectrum, channel D
// transversions, GC bins, channel B stop codons, channel C oxidative codons,
// interior CT cluster) are NOT recomputed here — PE mode is intended for
// raw bilateral 5'/3' damage QA, not full library profiling. Use SE on
// merged reads when those signals are needed.
bool FrameSelector::update_sample_profile_paired(
    SampleDamageProfile& profile,
    std::string_view r1,
    std::string_view r2)
{
    if (r1.length() < 30 || r2.length() < 30) return false;

    // Short-insert detection via R1/R2 overlap. When the molecule (insert)
    // is shorter than the read length, R1 reads through the molecule into
    // adapter A1 and R2 into adapter A2; the per-position damage windows
    // and tail counters then mix molecule and adapter bases, producing an
    // "anti-damage" shape at the 3' end (R2 first 15 bases are largely
    // adapter, complement-mapped into top-strand frame as A-depletion at
    // the 3'-end window).
    //
    // For an insert of length M, R1[M-K..M-1] is the molecule's 3' tail
    // and should reverse-complement to R2[0..K-1] (which reads the
    // molecule's bottom strand from the same end). We scan plausible M
    // values; require K=15 bases of overlap with at most 3 mismatches
    // (allows for sequencing error and aDNA damage at the molecule 3'
    // end). When overlap is detected, the pair is short-insert by
    // definition and belongs to the merged-read SE workflow — skip it.
    // Native PE is intended for true long-insert pairs (insert > read
    // length) where R1 and R2 do not overlap and the per-position
    // windows are clean molecule bases.
    {
        auto rc_base = [](char c) -> char {
            switch (c) { case 'A': return 'T'; case 'T': return 'A';
                         case 'C': return 'G'; case 'G': return 'C'; }
            return 'N';
        };
        constexpr int CHECK_LEN = 15;
        constexpr int MAX_MISMATCH = 3;
        const int max_M = static_cast<int>(std::min(r1.length(), r2.length()));
        bool overlap_found = false;
        for (int M = CHECK_LEN; M <= max_M; ++M) {
            int mismatches = 0;
            for (int i = 0; i < CHECK_LEN; ++i) {
                const char r1b = fast_upper(r1[M - 1 - i]);
                const char r2b = fast_upper(r2[i]);
                if (r1b == 'N' || r2b == 'N' || rc_base(r2b) != r1b) {
                    if (++mismatches > MAX_MISMATCH) break;
                }
            }
            if (mismatches <= MAX_MISMATCH) {
                overlap_found = true;
                break;
            }
        }
        if (overlap_found) {
            profile.pe_short_insert_skipped++;
            return false;
        }
    }

    profile.input_mode = SampleDamageProfile::InputMode::PAIRED;

    const size_t l1 = r1.length();
    const size_t l2 = r2.length();

    // R1 → 5' end counters
    for (size_t i = 0; i < std::min(size_t(15), l1); ++i) {
        char b = fast_upper(r1[i]);
        if (b == 'T')      { profile.t_freq_5prime[i]++; profile.tc_total_5prime[i]++; }
        else if (b == 'C') { profile.c_freq_5prime[i]++; profile.tc_total_5prime[i]++; }
        if (b == 'A')      profile.a_freq_5prime[i]++;
        else if (b == 'G') profile.g_freq_5prime[i]++;
    }

    // R2 → 3' end counters (complement-mapped: R2 reads bottom strand from
    // the molecule 3' inward, so R2[i] = complement(top_strand_at_3prime[i]).
    // The damage signal we want is G->A on top strand 3' end, which appears
    // as C->T on R2. Map R2 base to its complement, then accumulate with the
    // same logic as the SE 3'-end loop.
    for (size_t i = 0; i < std::min(size_t(15), l2); ++i) {
        char b = fast_upper(r2[i]);
        // complement: R2[i] → top_strand_at_3prime[i]
        char top;
        switch (b) {
            case 'A': top = 'T'; break;
            case 'T': top = 'A'; break;
            case 'C': top = 'G'; break;
            case 'G': top = 'C'; break;
            default:  top = 'N'; break;
        }
        if (top == 'A')      { profile.a_freq_3prime[i]++; profile.ag_total_3prime[i]++; }
        else if (top == 'G') { profile.g_freq_3prime[i]++; profile.ag_total_3prime[i]++; }
        if (top == 'T')      profile.t_freq_3prime[i]++;
        else if (top == 'C') profile.c_freq_3prime[i]++;
    }

    // 5' tail from R1, 3' tail from R2 (complement-mapped)
    {
        const int lo = SampleDamageProfile::BG_TAIL_LO;
        const int hi = SampleDamageProfile::BG_TAIL_HI;
        for (int i = lo; i <= hi && static_cast<size_t>(i) < l1; ++i) {
            const int idx = i - lo;
            const char b = fast_upper(r1[i]);
            if (b == 'T')      { profile.tail_t_5prime[idx]++; profile.tail_tc_5prime[idx]++; }
            else if (b == 'C') profile.tail_tc_5prime[idx]++;
        }
        for (int i = lo; i <= hi && static_cast<size_t>(i) < l2; ++i) {
            const int idx = i - lo;
            const char b = fast_upper(r2[i]);
            char top;
            switch (b) {
                case 'A': top = 'T'; break;
                case 'T': top = 'A'; break;
                case 'C': top = 'G'; break;
                case 'G': top = 'C'; break;
                default:  top = 'N'; break;
            }
            if (top == 'A')      { profile.tail_a_3prime[idx]++; profile.tail_ag_3prime[idx]++; }
            else if (top == 'G') profile.tail_ag_3prime[idx]++;
        }
    }

    // Interior baseline from R1's middle third (R2's middle would mostly
    // overlap for short inserts; tracking only R1's middle avoids double-
    // counting and keeps the baseline consistent with SE behavior).
    {
        constexpr size_t INTERIOR_TERM_PAD = 15;
        size_t mid_start = l1 / 3;
        size_t mid_end   = 2 * l1 / 3;
        if (mid_start < INTERIOR_TERM_PAD) mid_start = INTERIOR_TERM_PAD;
        if (l1 > INTERIOR_TERM_PAD && mid_end + INTERIOR_TERM_PAD > l1)
            mid_end = l1 - INTERIOR_TERM_PAD;
        if (mid_start < mid_end) {
            for (size_t i = mid_start; i < mid_end; ++i) {
                char b = fast_upper(r1[i]);
                if (b == 'T') profile.baseline_t_freq++;
                else if (b == 'C') profile.baseline_c_freq++;
                else if (b == 'A') profile.baseline_a_freq++;
                else if (b == 'G') profile.baseline_g_freq++;
            }
        }
    }

    // Codon-position counters: 5' from R1, 3' from R2 (complement-mapped)
    for (size_t i = 0; i < std::min(size_t(15), l1); ++i) {
        char b = fast_upper(r1[i]);
        int cp = i % 3;
        if (b == 'T') profile.codon_pos_t_count_5prime[cp]++;
        else if (b == 'C') profile.codon_pos_c_count_5prime[cp]++;
    }
    for (size_t i = 0; i < std::min(size_t(15), l2); ++i) {
        char b = fast_upper(r2[i]);
        char top;
        switch (b) { case 'A': top='T'; break; case 'T': top='A'; break;
                     case 'C': top='G'; break; case 'G': top='C'; break;
                     default: top='N'; }
        // Codon position in R2: position i in R2 == position (l_mol-1-i) in
        // top strand. Without alignment we can't recover exact codon phase
        // on the top strand, so we use the natural R2 codon phase (which
        // matches the molecule's 3' frame for inserts of length 3k).
        int cp = i % 3;
        if (top == 'A') profile.codon_pos_a_count_3prime[cp]++;
        else if (top == 'G') profile.codon_pos_g_count_3prime[cp]++;
    }

    // 5' hexamer + interior hexamer (from R1)
    if (l1 >= 18) {
        char hex_5prime[7];
        bool valid_5prime = true;
        for (int i = 0; i < 6; ++i) {
            char b = fast_upper(r1[i]);
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') { valid_5prime = false; break; }
            hex_5prime[i] = b;
        }
        hex_5prime[6] = '\0';
        if (valid_5prime) {
            uint32_t code = encode_hexamer(hex_5prime);
            if (code < 4096) {
                profile.hexamer_count_5prime[code] += 1.0;
                profile.n_hexamers_5prime++;
            }
        }

        size_t interior_start = l1 / 2 - 3;
        char hex_interior[7];
        bool valid_interior = true;
        for (int i = 0; i < 6; ++i) {
            char b = fast_upper(r1[interior_start + i]);
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') { valid_interior = false; break; }
            hex_interior[i] = b;
        }
        hex_interior[6] = '\0';
        if (valid_interior) {
            uint32_t code = encode_hexamer(hex_interior);
            if (code < 4096) {
                profile.hexamer_count_interior[code] += 1.0;
                profile.n_hexamers_interior++;
            }
        }
    }

    profile.n_reads++;
    return true;
}

void FrameSelector::update_sample_profile_weighted(
    SampleDamageProfile& profile,
    std::string_view seq,
    float weight) {

    if (seq.length() < 30 || weight < 0.001f) return;

    size_t len = seq.length();

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        char base = fast_upper(seq[i]);
        if (base == 'T') {
            profile.t_freq_5prime[i] += weight;
            profile.tc_total_5prime[i] += weight;
        } else if (base == 'C') {
            profile.c_freq_5prime[i] += weight;
            profile.tc_total_5prime[i] += weight;
        }
    }

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        if (base == 'A') {
            profile.a_freq_3prime[i] += weight;
            profile.ag_total_3prime[i] += weight;
        } else if (base == 'G') {
            profile.g_freq_3prime[i] += weight;
            profile.ag_total_3prime[i] += weight;
        }
    }

    // Tail-anchored background sampling (weighted variant)
    {
        const int lo = SampleDamageProfile::BG_TAIL_LO;
        const int hi = SampleDamageProfile::BG_TAIL_HI;
        for (int i = lo; i <= hi && static_cast<size_t>(i) < len; ++i) {
            const int idx = i - lo;
            const char b5 = fast_upper(seq[i]);
            if (b5 == 'T') { profile.tail_t_5prime[idx] += weight; profile.tail_tc_5prime[idx] += weight; }
            else if (b5 == 'C') { profile.tail_tc_5prime[idx] += weight; }

            const size_t pos3 = len - 1 - i;
            const char b3 = fast_upper(seq[pos3]);
            if (b3 == 'A') { profile.tail_a_3prime[idx] += weight; profile.tail_ag_3prime[idx] += weight; }
            else if (b3 == 'G') { profile.tail_ag_3prime[idx] += weight; }
        }
    }

    constexpr size_t INTERIOR_TERM_PAD = 15;
    size_t mid_start = len / 3;
    size_t mid_end   = 2 * len / 3;
    if (mid_start < INTERIOR_TERM_PAD)     mid_start = INTERIOR_TERM_PAD;
    if (len > INTERIOR_TERM_PAD && mid_end + INTERIOR_TERM_PAD > len)
        mid_end = len - INTERIOR_TERM_PAD;
    const bool interior_safe = (mid_start < mid_end);
    if (interior_safe) {
        for (size_t i = mid_start; i < mid_end; ++i) {
            char base = fast_upper(seq[i]);
            if (base == 'T') profile.baseline_t_freq += weight;
            else if (base == 'C') profile.baseline_c_freq += weight;
            else if (base == 'A') profile.baseline_a_freq += weight;
            else if (base == 'G') profile.baseline_g_freq += weight;
        }
    }

    size_t weight_count = std::max(size_t(1), static_cast<size_t>(weight * 10 + 0.5));
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        int codon_pos = i % 3;
        char base = fast_upper(seq[i]);
        if (base == 'T') profile.codon_pos_t_count_5prime[codon_pos] += weight_count;
        else if (base == 'C') profile.codon_pos_c_count_5prime[codon_pos] += weight_count;
    }

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        int codon_pos = (len - 1 - i) % 3;
        char base = fast_upper(seq[pos]);
        if (base == 'A') profile.codon_pos_a_count_3prime[codon_pos] += weight_count;
        else if (base == 'G') profile.codon_pos_g_count_3prime[codon_pos] += weight_count;
    }

    for (size_t i = 0; i < std::min(size_t(5), len - 1); ++i) {
        char base = fast_upper(seq[i]);
        char next = fast_upper(seq[i + 1]);

        if (next == 'G') {
            if (base == 'C') {
                profile.cpg_c_count += weight_count;
            } else if (base == 'T') {
                profile.cpg_t_count += weight_count;
            }
        } else {
            if (base == 'C') {
                profile.non_cpg_c_count += weight_count;
            } else if (base == 'T') {
                profile.non_cpg_t_count += weight_count;
            }
        }
    }

    profile.n_reads++;
}

void FrameSelector::reset_sample_profile(SampleDamageProfile& profile) {
    // Default-construct first: covers every field (including ones the previous
    // hand-rolled reset missed: d_max_*, damage_status, composition_bias_*,
    // inverted_pattern_*, library_bic_*, etc.) and stays in sync as new
    // members are added to SampleDamageProfile. Then restore the few
    // historical non-zero defaults the manual reset relied on (below).
    profile = SampleDamageProfile{};
    for (int i = 0; i < 15; ++i) {
        profile.t_freq_5prime[i] = 0.0;
        profile.c_freq_5prime[i] = 0.0;
        profile.a_freq_5prime[i] = 0.0;
        profile.g_freq_5prime[i] = 0.0;
        profile.a_freq_3prime[i] = 0.0;
        profile.g_freq_3prime[i] = 0.0;
        profile.t_freq_3prime[i] = 0.0;
        profile.c_freq_3prime[i] = 0.0;
        profile.damage_rate_5prime[i] = 0.0f;
        profile.damage_rate_3prime[i] = 0.0f;
        profile.tc_total_5prime[i] = 0.0;
        profile.ag_total_3prime[i] = 0.0;
    }

    profile.baseline_t_freq = 0.0;
    profile.baseline_c_freq = 0.0;
    profile.baseline_a_freq = 0.0;
    profile.baseline_g_freq = 0.0;

    for (int p = 0; p < 3; ++p) {
        profile.codon_pos_t_count_5prime[p] = 0;
        profile.codon_pos_c_count_5prime[p] = 0;
        profile.codon_pos_a_count_3prime[p] = 0;
        profile.codon_pos_g_count_3prime[p] = 0;
        profile.codon_pos_t_rate_5prime[p] = 0.5f;
        profile.codon_pos_a_rate_3prime[p] = 0.5f;
    }

    profile.cpg_c_count = 0;
    profile.cpg_t_count = 0;
    profile.non_cpg_c_count = 0;
    profile.non_cpg_t_count = 0;
    profile.cpg_damage_rate = 0.0f;
    profile.non_cpg_damage_rate = 0.0f;

    for (int ctx = 0; ctx < SampleDamageProfile::N_CT_CTX; ++ctx) {
        profile.ct_ctx_t_5prime[ctx].fill(0.0f);
        profile.ct_ctx_total_5prime[ctx].fill(0.0f);
        profile.ct_ctx_t_interior[ctx] = 0.0f;
        profile.ct_ctx_total_interior[ctx] = 0.0f;
    }
    profile.fit_baseline_ct5_cpg_like    = std::numeric_limits<float>::quiet_NaN();
    profile.fit_baseline_ct5_noncpg_like = std::numeric_limits<float>::quiet_NaN();
    profile.dmax_ct5_cpg_like    = std::numeric_limits<float>::quiet_NaN();
    profile.dmax_ct5_noncpg_like = std::numeric_limits<float>::quiet_NaN();
    profile.cpg_ratio     = std::numeric_limits<float>::quiet_NaN();
    profile.log2_cpg_ratio = std::numeric_limits<float>::quiet_NaN();

    // Reset upstream-context-aware accumulators
    for (int uctx = 0; uctx < SampleDamageProfile::N_UPSTREAM_CTX; ++uctx) {
        profile.ct5_t_by_upstream[uctx].fill(0.0);
        profile.ct5_total_by_upstream[uctx].fill(0.0);
        profile.ct5_t_interior_by_upstream[uctx] = 0.0;
        profile.ct5_total_interior_by_upstream[uctx] = 0.0;
        profile.dmax_ct5_by_upstream[uctx] = std::numeric_limits<float>::quiet_NaN();
        profile.baseline_ct5_by_upstream[uctx] = std::numeric_limits<float>::quiet_NaN();
        profile.cov_ct5_terminal_by_upstream[uctx] = 0.0f;
        profile.cov_ct5_interior_by_upstream[uctx] = 0.0f;
    }
    profile.dipyr_contrast = std::numeric_limits<float>::quiet_NaN();
    profile.cpg_contrast = std::numeric_limits<float>::quiet_NaN();
    profile.context_heterogeneity_chi2 = 0.0f;
    profile.context_heterogeneity_p = 1.0f;
    profile.context_heterogeneity_detected = false;
    profile.effcov_ct5_cpg_like_terminal    = 0.0f;
    profile.effcov_ct5_noncpg_like_terminal = 0.0f;
    profile.effcov_ct5_cpg_like_interior    = 0.0f;
    profile.effcov_ct5_noncpg_like_interior = 0.0f;
    profile.cov_ct5_cpg_like_terminal       = 0.0f;
    profile.cov_ct5_noncpg_like_terminal    = 0.0f;
    profile.cov_ct5_cpg_like_interior       = 0.0f;
    profile.cov_ct5_noncpg_like_interior    = 0.0f;
    profile.fit_positions_ct5_cpg_like    = 0;
    profile.fit_positions_ct5_noncpg_like = 0;
    profile.oxog16_t.fill(0.0f);
    profile.oxog16_a_rc.fill(0.0f);
    profile.s_oxog_16ctx.fill(0.0f);
    profile.cov_oxog_16ctx.fill(0.0f);
    profile.tri_5prime_terminal.fill(0);
    profile.tri_5prime_interior.fill(0);
    profile.tri_3prime_terminal.fill(0);
    profile.tri_3prime_interior.fill(0);

    profile.max_damage_5prime = 0.0f;
    profile.max_damage_3prime = 0.0f;
    profile.sample_damage_prob = 0.0f;
    profile.lambda_5prime = 0.3f;
    profile.lambda_3prime = 0.3f;
    profile.terminal_shift_5prime = 0.0f;
    profile.terminal_shift_3prime = 0.0f;
    profile.terminal_z_5prime = 0.0f;
    profile.terminal_z_3prime = 0.0f;
    profile.terminal_inversion = false;
    profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
    profile.library_type_auto_detected = false;
    profile.tc_total_3prime.fill(0.0);

    profile.hexamer_count_5prime.fill(0.0);
    profile.hexamer_count_interior.fill(0.0);
    profile.n_hexamers_5prime = 0;
    profile.n_hexamers_interior = 0;
    profile.hexamer_damage_llr = 0.0f;

    profile.convertible_caa_5prime.fill(0.0);
    profile.convertible_taa_5prime.fill(0.0);
    profile.convertible_cag_5prime.fill(0.0);
    profile.convertible_tag_5prime.fill(0.0);
    profile.convertible_cga_5prime.fill(0.0);
    profile.convertible_tga_5prime.fill(0.0);
    profile.total_codons_5prime.fill(0.0);
    profile.convertible_caa_interior = 0.0;
    profile.convertible_taa_interior = 0.0;
    profile.convertible_cag_interior = 0.0;
    profile.convertible_tag_interior = 0.0;
    profile.convertible_cga_interior = 0.0;
    profile.convertible_tga_interior = 0.0;
    profile.total_codons_interior = 0.0;
    profile.stop_conversion_rate_baseline = 0.0f;
    profile.stop_decay_llr_5prime = 0.0f;
    profile.stop_amplitude_5prime = 0.0f;
    profile.channel_b_valid = false;
    profile.damage_validated = false;
    profile.damage_artifact = false;

    profile.joint_delta_max = 0.0f;
    profile.joint_lambda = 0.0f;
    profile.joint_a_max = 0.0f;
    profile.joint_log_lik_m1 = 0.0f;
    profile.joint_log_lik_m0 = 0.0f;
    profile.joint_delta_bic = 0.0f;
    profile.joint_bayes_factor = 0.0f;
    profile.joint_p_damage = 0.0f;
    profile.joint_model_valid = false;

    profile.mixture_K = 0;
    profile.mixture_n_components = 0;
    profile.mixture_d_population = 0.0f;
    profile.mixture_d_ancient = 0.0f;
    profile.mixture_d_reference = 0.0f;
    profile.mixture_pi_ancient = 0.0f;
    profile.mixture_bic = 0.0f;
    profile.mixture_converged = false;
    profile.mixture_identifiable = false;
    profile.gc_histogram.fill(0);
    profile.adaptive_gc_threshold = 0.0f;
    profile.gc_threshold_computed = false;
    profile.gc_bins = {};
    profile.gc_stratified_d_max_weighted = 0.0f;
    profile.gc_stratified_d_max_peak = 0.0f;
    profile.gc_peak_bin = -1;
    profile.gc_stratified_valid = false;
    profile.pi_damaged = 0.0f;
    profile.d_ancient = 0.0f;
    profile.d_population = 0.0f;
    profile.n_damaged_bins = 0;
    profile.n_reads_gc_filtered = 0;
    profile.n_reads_sampled = 0;

    profile.fit_offset_5prime = 1;
    profile.fit_offset_3prime = 1;

    profile.forced_library_type = SampleDamageProfile::LibraryType::UNKNOWN;
    profile.library_type_rescued = false;

    // oxoG / Channel-D accumulators
    profile.g_count_5prime.fill(0.0);
    profile.t_from_g_5prime.fill(0.0);
    profile.baseline_g_to_t_count = 0.0;
    profile.baseline_g_total = 0.0;
    profile.baseline_c_to_a_count = 0.0;
    profile.baseline_c_ox_total = 0.0;
    profile.convertible_gag_5prime.fill(0.0);
    profile.convertible_gaa_5prime.fill(0.0);
    profile.convertible_gga_5prime.fill(0.0);
    profile.convertible_tag_ox_5prime.fill(0.0);
    profile.convertible_taa_ox_5prime.fill(0.0);
    profile.convertible_tga_ox_5prime.fill(0.0);
    profile.c_count_ox_5prime.fill(0.0);
    profile.a_from_c_5prime.fill(0.0);

    // Interior oxoG codon accumulators (merged in merge_sample_profiles)
    profile.convertible_gag_interior = 0.0;
    profile.convertible_gaa_interior = 0.0;
    profile.convertible_gga_interior = 0.0;
    profile.convertible_tag_ox_interior = 0.0;
    profile.convertible_taa_ox_interior = 0.0;
    profile.convertible_tga_ox_interior = 0.0;

    profile.ox_is_artifact = false;
    profile.ox_d_max = 0.0f;

    // Neutral default for oxidation uniformity ratio
    profile.ox_uniformity_ratio = 1.0f;

    profile.n_reads = 0;
}

SampleDamageProfile FrameSelector::compute_sample_profile(
    const std::vector<std::string>& sequences) {

    SampleDamageProfile profile;
    reset_sample_profile(profile);

    for (const auto& seq : sequences) {
        update_sample_profile(profile, seq);
    }

    finalize_sample_profile(profile);
    return profile;
}

} // namespace taph
