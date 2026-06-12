// Validated reference-free ancient fraction (profile.pi) + prior-free per-read LLR.
// SOLUTION_pi_delta_dmax.md §6. Additive shadow-mode step 1 — the legacy mixture path is untouched and
// no consumer is wired; finalize_pi only populates the new profile.pi field.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

#include "taph/sample_damage_profile.hpp"
#include "taph/read_ancient_llr.hpp"
#include "damage_estimation_finalize_ctx.hpp"  // finalize_pi declaration

namespace taph {

namespace {

// 95% profile-likelihood interval for δ̂_l from the per-bin curve: invert profile_loglik at Δℓ = −1.92
// (χ²₁(0.95)/2), linearly interpolating the two crossings around the peak. SOLUTION §6.5(2): bound δ,
// then clip+divide — NEVER bootstrap pi_l (that censors the negative tail and reproduces the H0 floor).
// Returns {δ_lo ≥ 0, δ_hi}; rails to the grid edge where the curve does not cross (unidentified bound).
std::pair<double, double> profile_delta_ci(const BulkDamagePerBin& b) {
    const auto& g = b.profile_delta;
    const auto& l = b.profile_loglik;
    if (g.size() < 2 || g.size() != l.size()) return {b.delta, b.delta};
    constexpr double TH = -1.92;

    std::size_t pk = 0;
    for (std::size_t i = 1; i < l.size(); ++i)
        if (l[i] > l[pk]) pk = i;

    double lo = g.front();
    for (std::size_t i = pk; i-- > 0;) {
        if (l[i] <= TH) {
            const double t = (TH - l[i]) / (l[i + 1] - l[i]);
            lo = g[i] + t * (g[i + 1] - g[i]);
            break;
        }
    }
    double hi = g.back();
    for (std::size_t i = pk; i + 1 < l.size(); ++i) {
        if (l[i + 1] <= TH) {
            const double t = (TH - l[i]) / (l[i + 1] - l[i]);
            hi = g[i] + t * (g[i + 1] - g[i]);
            break;
        }
    }
    return {std::max(0.0, lo), hi};
}

}  // namespace

// Populate profile.pi from the validated bulk law: shortest live length bin (bins ascend in length),
// pi_l = clip(δ_0 / D_MAX_CONSERVED), w_length-gated, CI from the profile-likelihood curve. Runs after
// finalize_bulk (needs bulk_damage + w_length). Default state UNDETERMINED — never "modern".
void finalize_pi(SampleDamageProfile& profile) {
    profile.pi = DamageEstimate{};
    if (!profile.bulk_attempted) return;

    const BulkDamageResult& R = profile.bulk_damage;
    // Length-decay gate (replaces the w_length>0.6 gate): τ̂'s 95% CI must sit in the fast-decay regime
    // (state==DETECTED). A pervasive artifact gives τ→∞ (NOT_DETECTED) → pi stays UNDETERMINED.
    if (profile.tau.state != DamageConfidence::DETECTED) return;

    const BulkDamagePerBin* b0 = nullptr;     // shortest live identified bin
    for (const auto& bb : R.bins) {
        if (bb.identified && bb.n_reads > 0) { b0 = &bb; break; }
    }
    if (!b0 || b0->delta <= 0.0) return;      // no positive short-bin signal → UNDETERMINED

    const auto clip01 = [](double x) { return std::clamp(x, 0.0, 1.0); };
    const auto ci = profile_delta_ci(*b0);
    profile.pi.point = clip01(b0->delta  / D_MAX_CONSERVED);
    profile.pi.lo    = clip01(ci.first   / D_MAX_CONSERVED);
    profile.pi.hi    = clip01(ci.second  / D_MAX_CONSERVED);
    profile.pi.state = DamageConfidence::DETECTED;
}

std::optional<double> read_ancient_llr(const ReadDamageObs& obs, const SampleDamageProfile& profile) {
    if (profile.pi.state != DamageConfidence::DETECTED) return std::nullopt;
    // The LLR kernel uses per-position Briggs λ (position axis k). In reference-free mode λ is never
    // fitted (stays at 0.3 default, lambda_5prime_fitted=false). Scoring reads off an unfitted λ produces
    // confident-looking but wrong LLRs — refuse rather than mislead.
    if (!profile.lambda_5prime_fitted && !profile.lambda_3prime_fitted) return std::nullopt;

    const double A = D_MAX_CONSERVED;  // imported amplitude (§6.5(1)); NOT the mixture

    const auto end_llr = [&](const TerminalMismatch* s, std::uint32_t n, double lambda, double b) {
        double llr = 0.0;
        const double qm = std::clamp(b, 1e-6, 1.0 - 1e-6);  // modern: background only
        for (std::uint32_t i = 0; i < n; ++i) {
            const double decay = std::exp(-lambda * static_cast<double>(s[i].pos));
            const double qa = std::clamp(b + (1.0 - b) * A * decay, 1e-6, 1.0 - 1e-6);  // ancient
            llr += s[i].deaminated ? std::log(qa / qm)
                                   : std::log((1.0 - qa) / (1.0 - qm));
        }
        return llr;
    };

    double llr = end_llr(obs.five, obs.n_five, static_cast<double>(profile.lambda_5prime),
                         static_cast<double>(profile.fit_baseline_5prime));
    llr += end_llr(obs.three, obs.n_three, static_cast<double>(profile.lambda_3prime),
                   static_cast<double>(profile.fit_baseline_3prime));
    return llr;
}

}  // namespace taph
