// Reference-free length-decay constant τ of the per-bin terminal-deamination amplitude δ(L)≈A·exp(−L/τ).
//
// The per-position Briggs λ (estimate_briggs_params) fits decay along READ POSITION p — the wrong axis for
// the bulk/length-stratified law, where the discriminating signal is how the per-bin amplitude δ_l falls
// with READ LENGTH L. Genuine terminal deamination occupies a fixed terminal zone, so δ(L) decays fast
// (small τ); a pervasive per-base/whole-molecule artifact keeps δ flat-or-rising in L (τ→∞). finalize_tau
// profiles χ²(τ) over a 1-D grid via closed-form WLS on the live bins, censoring genuine-zero bins at a
// floor, and reads τ̂ + its 95% χ²-profile interval. NO external deps (pure std), NO Eigen.
//
// Live bin    : profile CI at Δℓ=−0.5 excludes zero, i.e. loglik(δ=0) < peak_loglik − 0.5.
//               Uses the profile-peak δ (not b.delta) to avoid point-estimate artefacts.
// Censored bin: boundary-peaked (peak at grid 0) ∧ loglik[1] < −1.0 (genuine zero).
// Uninformative: CI straddles zero but not a genuine zero → dropped silently.
// NOTE: the `identified` flag is intentionally NOT used here — it is an estimator-health boolean,
//   not an observation model, and silently changes τ by excluding flat-tail bins that should
//   instead be treated as genuinely uninformative via the CI criterion.
// Per-live weight w_l = 1/se², se = half-width of the δ profile-likelihood interval at Δℓ = −0.5
//   (se=0 → grid-resolution floor 0.015).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "taph/sample_damage_profile.hpp"
#include "taph/read_ancient_llr.hpp"

namespace taph {

namespace {

// Half-width of the δ profile-likelihood interval at Δℓ = −0.5, around the curve's peak. Linear
// interpolation of the two crossings; rails to the grid edge where the curve does not cross.
double profile_se(const BulkDamagePerBin& b, std::size_t peak) {
    const auto& g = b.profile_delta;
    const auto& l = b.profile_loglik;
    constexpr double TH = -0.5;

    double lo = g.front();
    for (std::size_t i = peak; i-- > 0;) {
        if (l[i] <= TH) {
            const double t = (TH - l[i]) / (l[i + 1] - l[i]);
            lo = g[i] + t * (g[i + 1] - g[i]);
            break;
        }
    }
    double hi = g.back();
    for (std::size_t i = peak; i + 1 < l.size(); ++i) {
        if (l[i + 1] <= TH) {
            const double t = (TH - l[i]) / (l[i + 1] - l[i]);
            hi = g[i] + t * (g[i + 1] - g[i]);
            break;
        }
    }
    return 0.5 * (hi - lo);
}

struct LiveBin {
    double mid;     // (length_lo + length_hi)/2
    double delta;
    double w;       // 1/se²
};

}  // namespace

DamageEstimate finalize_tau(const SampleDamageProfile& profile, const TauConfig& cfg) {
    DamageEstimate out;  // default UNDETERMINED, point/lo/hi = −1
    if (!profile.bulk_attempted) return out;

    const BulkDamageResult& R = profile.bulk_damage;

    std::vector<LiveBin> live;
    std::vector<double>  censored;   // mid lengths of genuine-zero bins
    for (const auto& b : R.bins) {
        const auto& pl = b.profile_loglik;
        const auto& pd = b.profile_delta;
        if (pl.size() < 2 || pl.size() != pd.size()) continue;

        std::size_t peak = 0;
        for (std::size_t i = 1; i < pl.size(); ++i)
            if (pl[i] > pl[peak]) peak = i;

        const double mid = 0.5 * (static_cast<double>(b.length_lo) + static_cast<double>(b.length_hi));

        if (pl.front() < pl[peak] - 0.5) {
            // CI at Δℓ=−0.5 excludes zero → live bin. Use profile-peak δ, not b.delta.
            double se = profile_se(b, peak);
            if (se <= 0.0) se = 0.015;
            live.push_back({mid, pd[peak], 1.0 / (se * se)});
        } else if (peak == 0 && pl[1] < -1.0) {
            // Boundary-peaked with steep loglik drop → genuine zero, censor.
            censored.push_back(mid);
        }
        // else: CI straddles zero, not a genuine zero → uninformative, drop.
    }

    // No live bin carries a positive, non-boundary δ → no terminal-decay signal to fit. This is the
    // amplitude-fail verdict (0 live bins < min_live_bins, Σδ=0 < a_min), so NOT_DETECTED, not the default
    // UNDETERMINED. (FLB10m: position-0-only artifact, every δ_l rails to 0.)
    if (live.empty()) {
        out.state = DamageConfidence::NOT_DETECTED;
        return out;
    }

    const double df = cfg.delta_floor;
    const double w_censor = 1.0 / (df * df);

    double best_tau = 0.0, best_chi2 = 0.0;
    bool   have_best = false;
    std::vector<std::pair<double, double>> curve;  // (τ, χ²)
    for (double tau = 10.0; tau <= 400.0 + 1e-9; tau += 1.0) {
        double sxd = 0.0, sxx = 0.0;
        for (const auto& b : live) {
            const double x = std::exp(-b.mid / tau);
            sxd += b.w * x * b.delta;
            sxx += b.w * x * x;
        }
        const double A = sxx > 0.0 ? sxd / sxx : 0.0;

        double chi2 = 0.0;
        for (const auto& b : live) {
            const double r = b.delta - A * std::exp(-b.mid / tau);
            chi2 += b.w * r * r;
        }
        for (double mid : censored) {
            const double pred = A * std::exp(-mid / tau);
            if (pred > df) {
                const double r = pred - df;
                chi2 += w_censor * r * r;
            }
        }
        curve.emplace_back(tau, chi2);
        if (!have_best || chi2 < best_chi2) { best_chi2 = chi2; best_tau = tau; have_best = true; }
    }

    double tau_lo = best_tau, tau_hi = best_tau;
    for (const auto& tc : curve) {
        if (tc.second - best_chi2 <= 3.84) {
            tau_lo = std::min(tau_lo, tc.first);
            tau_hi = std::max(tau_hi, tc.first);
        }
    }

    double amp = 0.0;
    for (const auto& b : live) amp += b.delta;
    const bool amplitude_ok = amp > cfg.a_min &&
                              static_cast<int>(live.size()) >= cfg.min_live_bins;

    out.point = best_tau;
    out.lo    = tau_lo;
    out.hi    = tau_hi;
    if (!amplitude_ok) {
        out.state = DamageConfidence::NOT_DETECTED;
    } else if (tau_hi < cfg.tau_max_detected) {
        out.state = DamageConfidence::DETECTED;
    } else if (tau_hi < cfg.tau_max_undetermined) {
        out.state = DamageConfidence::UNDETERMINED;
    } else {
        out.state = DamageConfidence::NOT_DETECTED;
    }
    return out;
}

}  // namespace taph
