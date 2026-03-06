// Sample-level damage profile management

#include "dart/frame_selector_decl.hpp"
#include "dart/codon_tables.hpp"
#include "dart/hexamer_tables.hpp"
#include <algorithm>
#include <cmath>
#include <array>

namespace dart {

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
    double sum_signal = 0.0, sum_weight = 0.0;
    for (int i = 1; i < 10; ++i) {
        if (total[i] < MIN_COVERAGE) continue;
        double weight = std::exp(-lambda * i);
        double excess = freq[i] - baseline;
        sum_signal += total[i] * excess / weight;
        sum_weight += total[i];
    }
    double raw_amplitude = (sum_weight > 0) ? sum_signal / sum_weight : 0.0;

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
    for (int i = 0; i < 15; ++i) {
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
    fit.lambda   = std::clamp(lambda,   0.05f,  0.50f);

    // WLS amplitude: A_hat = max(0, Σ n·x·(y - b) / Σ n·x²)
    double numer = 0.0, denom = 0.0;
    int n_valid = 0;
    for (int p = start_pos; p <= end_pos && p < 15; ++p) {
        double n = coverage[p];
        if (n < 100.0) continue;
        double x = std::exp(-fit.lambda * p);
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
        double x     = std::exp(-fit.lambda * p);
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

void FrameSelector::update_sample_profile(
    SampleDamageProfile& profile,
    const std::string& seq) {

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

    // Count bases in middle third (undamaged baseline)
    size_t mid_start = len / 3;
    size_t mid_end = 2 * len / 3;
    for (size_t i = mid_start; i < mid_end; ++i) {
        char base = fast_upper(seq[i]);
        if (base == 'T') profile.baseline_t_freq++;
        else if (base == 'C') profile.baseline_c_freq++;
        else if (base == 'A') profile.baseline_a_freq++;
        else if (base == 'G') profile.baseline_g_freq++;
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

    // Channel D: G count at 5' end for G→T asymmetry tracking
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        char base = fast_upper(seq[i]);
        if (base == 'G') {
            profile.g_count_5prime[i]++;
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
            }

            // Accumulate Channel A interior baseline (middle third)
            size_t mid_start = len / 3;
            size_t mid_end = 2 * len / 3;
            for (size_t i = mid_start; i < mid_end; ++i) {
                char base = fast_upper(seq[i]);
                if (base == 'T') bin.t_interior++;
                else if (base == 'C') bin.c_interior++;
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

    profile.n_reads++;
}

void FrameSelector::finalize_sample_profile(SampleDamageProfile& profile) {
    if (profile.n_reads == 0) return;

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
        if (profile.d_max_combined > 5.0f || profile.damage_validated) {
            threshold = 0.05f;  // 5% excess if deamination present
        }

        bool elevated = ox_stop_excess > threshold;
        bool uniform = profile.ox_uniformity_ratio > 0.85f && profile.ox_uniformity_ratio < 1.15f;

        if (profile.channel_c_valid && elevated && uniform) {
            profile.ox_damage_detected = true;
            profile.ox_d_max = std::max(0.0f, ox_stop_excess * 100.0f);  // Convert to percentage
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

    // Channel D: G→T asymmetry inferred from Channel C stop codon data
    // (reference-free; cannot directly measure G→T vs T→G without alignment)
    {
        if (profile.channel_c_valid) {
            profile.ox_gt_asymmetry = profile.ox_uniformity_ratio;
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

    // Briggs-like damage model parameter estimation (closed-form)
    auto estimate_briggs_params = [](const std::array<float, 15>& rates, float max_rate,
                                     float& delta_s, float& delta_d, float& lambda, float& r_squared) {
        delta_s = 0.0f;
        delta_d = 0.0f;
        lambda = 0.3f;
        r_squared = 0.0f;

        if (max_rate < 0.02f) return;

        delta_s = rates[0];

        float sum_background = 0.0f;
        for (int i = 10; i < 15; ++i) {
            sum_background += rates[i];
        }
        delta_d = sum_background / 5.0f;

        if (delta_s <= delta_d + 0.01f) {
            lambda = 0.3f;
            return;
        }

        // Weighted log-linear regression for lambda estimation
        float sum_xy = 0.0f;
        float sum_x2 = 0.0f;
        float sum_y = 0.0f;
        float sum_x = 0.0f;
        float sum_w = 0.0f;
        float sum_y2 = 0.0f;

        for (int pos = 0; pos < 15; ++pos) {
            float normalized = (rates[pos] - delta_d) / (delta_s - delta_d + 0.001f);
            normalized = std::clamp(normalized, 0.01f, 1.0f);

            float y = std::log(normalized);
            float x = static_cast<float>(pos);
            float w = normalized * normalized;

            sum_xy += w * x * y;
            sum_x2 += w * x * x;
            sum_y += w * y;
            sum_x += w * x;
            sum_w += w;
            sum_y2 += w * y * y;
        }

        float denom = sum_w * sum_x2 - sum_x * sum_x;
        if (std::abs(denom) < 0.001f) {
            lambda = 0.3f;
            return;
        }

        float slope = (sum_w * sum_xy - sum_x * sum_y) / denom;

        if (slope >= 0.0f) {
            lambda = 0.3f;
            return;
        }

        lambda = 1.0f - std::exp(slope);
        lambda = std::clamp(lambda, 0.05f, 0.95f);

        // Calculate R²
        float ss_tot = sum_y2 - (sum_y * sum_y) / sum_w;
        float intercept = (sum_y - slope * sum_x) / sum_w;
        float ss_res = 0.0f;
        for (int pos = 0; pos < 15; ++pos) {
            float normalized = (rates[pos] - delta_d) / (delta_s - delta_d + 0.001f);
            normalized = std::clamp(normalized, 0.01f, 1.0f);
            float y = std::log(normalized);
            float y_pred = slope * static_cast<float>(pos) + intercept;
            float w = normalized * normalized;
            ss_res += w * (y - y_pred) * (y - y_pred);
        }
        r_squared = (ss_tot > 0.001f) ? std::max(0.0f, 1.0f - ss_res / ss_tot) : 0.0f;
    };

    estimate_briggs_params(profile.damage_rate_5prime, profile.max_damage_5prime,
                           profile.delta_s_5prime, profile.delta_d_5prime,
                           profile.lambda_5prime, profile.r_squared_5prime);

    estimate_briggs_params(profile.damage_rate_3prime, profile.max_damage_3prime,
                           profile.delta_s_3prime, profile.delta_d_3prime,
                           profile.lambda_3prime, profile.r_squared_3prime);

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
    if (profile.forced_library_type != SampleDamageProfile::LibraryType::UNKNOWN) {
        profile.library_type = profile.forced_library_type;
        profile.library_type_auto_detected = false;
    } else {
        float lambda_lib = std::clamp(fit_lambda_5p, 0.05f, 0.50f);

        // ga3: A/(A+G) at 3' of read, positions 1-10 — smooth exponential decay (DS signal).
        //   For DS (BEST, top-strand read): bottom strand 3' C→T → G→A at read 3' pos 1-10.
        //   For SS (SCR, complement-strand read): biological 5' C→T → G→A at read 3' pos 1-10.
        //   Both library types can show ga3 signal; it is NOT the discriminating feature alone.
        //
        // ga0: A/(A+G) at 3' of read, position 0 only — single-position spike (SS signal).
        //   For SS (SCR complement strand): biological 5' C→T → complement G→A at read 3' pos-0.
        //     Observed: 3'_GA pos-0 ≈ 0.73 in SS environmental samples (far above baseline 0.46).
        //   For DS (BEST): 3' ligation chemistry DEPRESSES A at pos-0 (3'_GA pos-0 ≈ 0.36 < baseline).
        //     pos-0 spike is absent or negative → ga0_delta_bic ≤ 0.
        //
        // 4-model BIC on the 3' GA channel:
        //   M_bias = ga3_null + ga0_null   (no decay, no spike)
        //   M_DS   = ga3_alt  + ga0_null   (smooth decay at pos 1-10, no pos-0 spike)
        //   M_SS   = ga3_null + ga0_alt    (spike at pos-0, no smooth decay)
        //   M_mix  = ga3_alt  + ga0_alt    (spike + decay; environmental SS + composition bias → SS)
        ChannelDecayFit ga3 = fit_decay_fixed_lambda(
            profile.a_freq_3prime, profile.ag_total_3prime,
            static_cast<float>(baseline_ag), lambda_lib, 1, 10);
        // ga0: single-position test at 3' pos-0 (SS spike signal); min_valid=1 to allow 1-point fit
        ChannelDecayFit ga0 = fit_decay_fixed_lambda(
            profile.a_freq_3prime, profile.ag_total_3prime,
            static_cast<float>(baseline_ag), lambda_lib, 0, 0, 1);

        profile.libtype_fit_amplitude_3prime_ga = ga3.amplitude;
        profile.libtype_fit_amplitude_3prime_ct = ga0.amplitude;  // repurposed: pos-0 spike amplitude
        profile.libtype_delta_bic_3prime_ga     = ga3.delta_bic;
        profile.libtype_delta_bic_3prime_ct     = ga0.delta_bic;  // repurposed: pos-0 spike ΔBIC

        if (ga3.valid && ga0.valid) {
            // True BIC scores (stored for diagnostics/logging).
            profile.library_bic_bias = ga3.bic_null + ga0.bic_null;
            profile.library_bic_ds   = ga3.bic_alt  + ga0.bic_null;  // smooth decay, no spike
            profile.library_bic_ss   = ga3.bic_null + ga0.bic_alt;   // spike only, no decay
            profile.library_bic_mix  = ga3.bic_alt  + ga0.bic_alt;   // spike + decay

            // Guard against BIC hypersensitivity at high n and spurious 3' pos-0 spikes in
            // DS libraries: require a meaningful spike amplitude before ga0 can influence
            // classification. Amplitude is scale-independent unlike ΔBIC/n (which depends on
            // coverage at pos-0, not total reads). Empirical gap on 91 clay-test samples:
            //   DS FP max amplitude = 0.089 (LV7008888557, end-repair artifact)
            //   SS correct minimum  = 0.125 (LV7008890981)
            // Threshold 0.10 sits cleanly in the gap; amplitude < 0.10 collapses SS → bias
            // and mix → DS.
            const bool ga0_informative = (ga0.amplitude >= 0.10f);

            const float eff_bic_bias = profile.library_bic_bias;
            const float eff_bic_ds   = profile.library_bic_ds;
            const float eff_bic_ss   = ga0_informative ? profile.library_bic_ss : profile.library_bic_bias;
            const float eff_bic_mix  = ga0_informative ? profile.library_bic_mix : profile.library_bic_ds;

            double best = eff_bic_bias;
            profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;

            if (eff_bic_ds < best) {
                best = eff_bic_ds;
                profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            }
            if (eff_bic_ss < best) {
                best = eff_bic_ss;
                profile.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            }
            if (eff_bic_mix < best) {
                // M_mix: pos-0 spike (ga0_informative) and smooth 1-10 decay (ga3) are both
                // significant. Two principled rules discriminate DS from SS:
                //
                //   Rule 1 (no 5' damage): d_max_5 < 0.01 → SS.
                //     DS deamination always produces 5' C→T. Absent 5' signal means the library
                //     is SS (complement-orientation reads dominate; C→T is at their 3' end, not 5').
                //
                //   Rule 2 (has 5' damage): ctrl_shift_3prime ≥ 0.05 → SS, else DS.
                //     ctrl_shift = T/(T+C) excess at 3' (C→T channel). SS original-orientation
                //     reads contribute strong 3' C→T (shift = 0.112–0.188 in 91 clay samples).
                //     DS libraries have ctrl_shift ≤ 0.014 — well below the threshold.
                //
                // Validated on 91 clay-test samples (46 DS + 45 SS):
                //   DS M_mix FPs: d5=0.19–0.43, ctrl_shift=0.001–0.014 → DS ✓ (Rule 2)
                //   SS M_mix FNs: d5=0.000                             → SS ✓ (Rule 1)
                //   SS corrects (d5>0): ctrl_shift=0.112–0.188         → SS ✓ (Rule 2)
                const bool no_5p_damage   = (profile.max_damage_5prime < 0.01f);
                const bool ss_ctrl_signal = (profile.ctrl_shift_3prime >= 0.05f);
                profile.library_type = (no_5p_damage || ss_ctrl_signal)
                    ? SampleDamageProfile::LibraryType::SINGLE_STRANDED
                    : SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            }
        } else {
            profile.library_type = SampleDamageProfile::LibraryType::UNKNOWN;
        }
        profile.library_type_auto_detected = true;
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
        float raw_d_max_5prime = std::clamp(profile.damage_rate_5prime[0], 0.0f, 1.0f);
        float raw_d_max_3prime = std::clamp(profile.damage_rate_3prime[0], 0.0f, 1.0f);


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

                        for (int p = 0; p < 15; ++p) {
                            float decay = std::exp(-lambda * p);
                            float delta_p = b.d_max * decay;
                            float pi_undamaged = b.baseline_tc;
                            float pi_damaged = b.baseline_tc + (1.0f - b.baseline_tc) * delta_p;

                            pi_undamaged = std::clamp(pi_undamaged, 0.001f, 0.999f);
                            pi_damaged = std::clamp(pi_damaged, 0.001f, 0.999f);

                            double k = static_cast<double>(b.t_counts[p]);
                            double n = static_cast<double>(b.t_counts[p] + b.c_counts[p]);
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

                double total_c_sites = static_cast<double>(weight_sum);
                std::array<SuperRead, N_GC_BINS> super_reads;
                for (int bin = 0; bin < N_GC_BINS; ++bin) {
                    const auto& b = profile.gc_bins[bin];
                    super_reads[bin].gc_bin = bin;
                    super_reads[bin].c_sites = static_cast<double>(b.c_sites);

                    for (int p = 0; p < N_POSITIONS; ++p) {
                        super_reads[bin].k_tc[p] = static_cast<double>(b.t_counts[p]);
                        super_reads[bin].n_tc[p] = static_cast<double>(b.t_counts[p] + b.c_counts[p]);
                    }

                    // A/(A+G) not tracked per bin; approximate from global proportionally
                    double bin_fraction = b.c_sites > 0 ? static_cast<double>(b.c_sites) / total_c_sites : 0.0;
                    for (int p = 0; p < N_POSITIONS; ++p) {
                        super_reads[bin].k_ag[p] = profile.a_freq_5prime[p] * bin_fraction * base_ag_total;
                        super_reads[bin].n_ag[p] = (profile.a_freq_5prime[p] + profile.g_freq_5prime[p]) * bin_fraction * base_ag_total;
                    }

                    for (int p = 0; p < N_POSITIONS; ++p) {
                        super_reads[bin].k_stop[p] = static_cast<double>(b.stop_counts[p]);
                        super_reads[bin].n_stop[p] = static_cast<double>(b.stop_counts[p] + b.pre_counts[p]);
                    }

                    super_reads[bin].k_tc_int = static_cast<double>(b.t_interior);
                    super_reads[bin].n_tc_int = static_cast<double>(b.t_interior + b.c_interior);
                    super_reads[bin].k_stop_int = static_cast<double>(b.stop_interior);
                    super_reads[bin].n_stop_int = static_cast<double>(b.stop_interior + b.pre_interior);

                    super_reads[bin].k_ag_int = profile.baseline_a_freq * bin_fraction * base_ag_total;
                    super_reads[bin].n_ag_int = (profile.baseline_a_freq + profile.baseline_g_freq) * bin_fraction * base_ag_total;
                }

                auto mixture_result = MixtureDamageModel::fit(super_reads);
                profile.mixture_K = mixture_result.K;
                profile.mixture_d_population = mixture_result.d_population;
                profile.mixture_d_ancient = mixture_result.d_ancient;
                profile.mixture_d_reference = mixture_result.d_reference;
                profile.mixture_pi_ancient = mixture_result.pi_ancient;
                profile.mixture_bic = mixture_result.bic;
                profile.mixture_converged = mixture_result.converged;

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
                // ds libraries: do not fall back to Channel A; that reintroduces the
                // compositional false positives the joint model suppresses.
                profile.d_max_5prime = 0.0f;
                profile.d_max_3prime = 0.0f;
                profile.d_max_combined = 0.0f;
                profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
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
                profile.d_max_5prime = 0.0f;
                profile.d_max_3prime = 0.0f;
                profile.d_max_combined = 0.0f;
                profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
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
}

void FrameSelector::merge_sample_profiles(SampleDamageProfile& dst, const SampleDamageProfile& src) {
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
        dst.g_count_5prime[i] += src.g_count_5prime[i];
    }
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
            db.stop_counts[p] += sb.stop_counts[p];
            db.pre_counts[p] += sb.pre_counts[p];
        }
        db.t_interior += sb.t_interior;
        db.c_interior += sb.c_interior;
        db.stop_interior += sb.stop_interior;
        db.pre_interior += sb.pre_interior;
        db.n_reads += sb.n_reads;
    }

    dst.n_reads += src.n_reads;
}

void FrameSelector::update_sample_profile_weighted(
    SampleDamageProfile& profile,
    const std::string& seq,
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

    size_t mid_start = len / 3;
    size_t mid_end = 2 * len / 3;
    for (size_t i = mid_start; i < mid_end; ++i) {
        char base = fast_upper(seq[i]);
        if (base == 'T') profile.baseline_t_freq += weight;
        else if (base == 'C') profile.baseline_c_freq += weight;
        else if (base == 'A') profile.baseline_a_freq += weight;
        else if (base == 'G') profile.baseline_g_freq += weight;
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
    for (int i = 0; i < 15; ++i) {
        profile.t_freq_5prime[i] = 0.0;
        profile.c_freq_5prime[i] = 0.0;
        profile.a_freq_3prime[i] = 0.0;
        profile.g_freq_3prime[i] = 0.0;
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
    profile.mixture_d_population = 0.0f;
    profile.mixture_d_ancient = 0.0f;
    profile.mixture_d_reference = 0.0f;
    profile.mixture_pi_ancient = 0.0f;
    profile.mixture_bic = 0.0f;
    profile.mixture_converged = false;

    profile.n_reads = 0;
}

SampleDamageProfile FrameSelector::compute_sample_profile(
    const std::vector<std::string>& sequences) {

    SampleDamageProfile profile;

    for (const auto& seq : sequences) {
        update_sample_profile(profile, seq);
    }

    finalize_sample_profile(profile);
    return profile;
}

} // namespace dart
