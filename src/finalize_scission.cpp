// Reference-free scission rate γ (bp⁻¹): left-truncated exponential MLE over the right tail of the
// fine fragment-length histogram. Model: P(L | L ≥ L_mode) = γ·exp(−γ·(L−L_mode)).
//
// L_mode = mean of the peak fine bin (exact via len_sum, not bin midpoint approximation).
// Right tail = peak bin + all bins after it.
// Sufficient statistic: S = Σ_{i≥peak} n_i·(x̄_i − L_mode).
// MLE:  γ̂ = n_tail / S.
// SE:   γ̂ / √n_tail   (Fisher information of left-truncated exponential).
// CI:   Wald 95%, lo = max(0, γ̂ − 1.96·SE), hi = γ̂ + 1.96·SE.
//
// Quality gates: ≥3 non-empty bins in the tail (peak included), ≥50 tail reads, S > 0.

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "taph/sample_damage_profile.hpp"
#include "taph/read_ancient_llr.hpp"

namespace taph {

ScissionEstimate finalize_scission(const SampleDamageProfile& profile) {
    ScissionEstimate out;

    constexpr int NFINE = SampleDamageProfile::N_LEN_FINE;
    const auto& LB = profile.len_bins;

    uint64_t n_total  = 0;
    uint64_t len_total = 0;
    for (int i = 0; i < NFINE; ++i) {
        n_total   += LB[i].n_reads;
        len_total += LB[i].len_sum;
    }
    out.n_total = static_cast<int64_t>(n_total);
    if (n_total == 0) return out;
    out.mean_length = static_cast<double>(len_total) / static_cast<double>(n_total);

    // Peak bin: mode by read count
    int peak_i = 0;
    for (int i = 1; i < NFINE; ++i)
        if (LB[i].n_reads > LB[peak_i].n_reads) peak_i = i;
    if (LB[peak_i].n_reads == 0) return out;

    // L_mode from actual read-length mean of peak bin
    const double L_mode = static_cast<double>(LB[peak_i].len_sum)
                        / static_cast<double>(LB[peak_i].n_reads);
    out.modal_length = L_mode;

    // Right-tail sufficient statistics
    uint64_t n_tail   = 0;
    double   S        = 0.0;  // Σ n_i · (x̄_i − L_mode)
    int      tail_bins = 0;
    for (int i = peak_i; i < NFINE; ++i) {
        const uint64_t ni = LB[i].n_reads;
        if (ni == 0) continue;
        const double xi = static_cast<double>(LB[i].len_sum) / static_cast<double>(ni);
        n_tail += ni;
        S      += static_cast<double>(ni) * (xi - L_mode);
        ++tail_bins;
    }
    out.n_tail = static_cast<int64_t>(n_tail);

    if (tail_bins < 3 || S <= 0.0 || n_tail < 50) return out;

    const double gamma = static_cast<double>(n_tail) / S;
    const double se    = gamma / std::sqrt(static_cast<double>(n_tail));

    out.gamma  = gamma;
    out.lo     = std::max(0.0, gamma - 1.96 * se);
    out.hi     = gamma + 1.96 * se;
    out.fitted = true;
    return out;
}

} // namespace taph
