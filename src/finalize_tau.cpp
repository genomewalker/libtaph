// Reference-free length-decay constant τ and floor-model decomposition for δ(L) = f0 + A·exp(−L/τ).
//
// Model: terminal C→T excess has two additive components:
//   - Overhang (ss-overhang deamination): end-concentrated, length-dependent: A·exp(−L/τ). Small τ → fast
//     decay with read length → genuine terminal deamination mechanism.
//   - Pervasive floor (bulk hydrolysis/oxidation): length-independent baseline f0 ≥ 0. f0>0 indicates
//     in-duplex deamination accumulated over burial time.
//
// At each τ in [10, 400] the optimal {f0, A} is found via closed-form 2-param WLS (non-negativity
// constrained). χ²(τ) = Σ w_l·(δ_l − f0 − A·exp(−mid_l/τ))² + censor terms. Profile CI on τ is the
// set {τ : χ²(τ) − χ²_min ≤ 3.84} (1-DOF; f0 and A are profiled out).
//
// Outputs in DamageEstimate.tau: τ̂ + 95% CI, state, AND {f0, amplitude, overhang_fraction} with
// delta-method CI for overhang_fraction = A/(A+f0).
//
// Live bin    : profile CI at Δℓ=−0.5 excludes zero (pl.front() < pl[peak] − 0.5).
// Censored bin: boundary-peaked ∧ loglik[1] < −1.0 (genuine zero).
// Uninformative: CI straddles zero but not a genuine zero → dropped.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "taph/sample_damage_profile.hpp"
#include "taph/read_ancient_llr.hpp"

namespace taph {

namespace {

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

struct LiveBin { double mid, delta, w; };

// 2-param non-negative-constrained WLS: δ = f0 + A·x, x = exp(-mid/tau).
// Returns {f0, A} — both guaranteed ≥ 0.
struct FloorFit { double f0, A; };

FloorFit wls2(const std::vector<LiveBin>& live, double tau) {
    double sw = 0.0, swx = 0.0, swxx = 0.0, swy = 0.0, swxy = 0.0;
    for (const auto& b : live) {
        const double x = std::exp(-b.mid / tau);
        sw   += b.w;
        swx  += b.w * x;
        swxx += b.w * x * x;
        swy  += b.w * b.delta;
        swxy += b.w * x * b.delta;
    }
    const double det = sw * swxx - swx * swx;
    double f0, A;
    if (det > 0.0) {
        f0 = (swxx * swy  - swx  * swxy) / det;
        A  = (sw   * swxy - swx  * swy ) / det;
    } else {
        // Collinear regressors (1 live bin or τ→∞): fall back to 1-param
        f0 = 0.0;
        A  = swxx > 0.0 ? swxy / swxx : 0.0;
    }
    // Non-negative projection (constrained WLS on boundary of ℝ²₊)
    if (f0 < 0.0 && A >= 0.0) {
        f0 = 0.0;
        A  = swxx > 0.0 ? swxy / swxx : 0.0;
    } else if (f0 >= 0.0 && A < 0.0) {
        A  = 0.0;
        f0 = sw > 0.0 ? swy / sw : 0.0;
    } else if (f0 < 0.0 && A < 0.0) {
        f0 = 0.0; A = 0.0;
    }
    return {f0, A};
}

}  // namespace

DamageEstimate finalize_tau(const SampleDamageProfile& profile, const TauConfig& cfg) {
    DamageEstimate out;
    if (!profile.bulk_attempted) return out;

    const BulkDamageResult& R = profile.bulk_damage;

    std::vector<LiveBin> live;
    std::vector<double>  censored;
    for (const auto& b : R.bins) {
        const auto& pl = b.profile_loglik;
        const auto& pd = b.profile_delta;
        if (pl.size() < 2 || pl.size() != pd.size()) continue;

        std::size_t peak = 0;
        for (std::size_t i = 1; i < pl.size(); ++i)
            if (pl[i] > pl[peak]) peak = i;

        const double mid = 0.5 * (static_cast<double>(b.length_lo) + static_cast<double>(b.length_hi));

        if (pl.front() < pl[peak] - 0.5) {
            double se = profile_se(b, peak);
            if (se <= 0.0) se = 0.015;
            live.push_back({mid, pd[peak], 1.0 / (se * se)});
        } else if (peak == 0 && pl[1] < -1.0) {
            censored.push_back(mid);
        }
    }

    if (live.empty()) {
        out.state = DamageConfidence::NOT_DETECTED;
        return out;
    }

    const double df       = cfg.delta_floor;
    const double w_censor = 1.0 / (df * df);

    double best_tau = 0.0, best_chi2 = 0.0;
    bool   have_best = false;
    std::vector<std::pair<double, double>> curve;

    for (double tau = 10.0; tau <= 400.0 + 1e-9; tau += 1.0) {
        const auto [f0, A] = wls2(live, tau);

        double chi2 = 0.0;
        for (const auto& b : live) {
            const double r = b.delta - f0 - A * std::exp(-b.mid / tau);
            chi2 += b.w * r * r;
        }
        for (double mid : censored) {
            const double pred = f0 + A * std::exp(-mid / tau);
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

    // Amplitude check uses total δ capacity (A+f0 at the best-τ solution)
    const auto [f0_hat, A_hat] = wls2(live, best_tau);
    const double total_amp = A_hat + f0_hat;
    const bool amplitude_ok = total_amp > cfg.a_min &&
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

    // Floor decomposition at best_tau
    out.f0        = f0_hat;
    out.amplitude = A_hat;

    if (total_amp > 0.0) {
        out.overhang_fraction = A_hat / total_amp;

        // Delta-method CI on overhang_fraction when both components are positive (unconstrained interior).
        // Cov(f0,A) = [[S_wxx, -S_wx], [-S_wx, S_w]] / det  (exact WLS covariance, weights = 1/se²).
        // Var(r) = [f0²·S_w + 2·A·f0·S_wx + A²·S_wxx] / (det·(A+f0)⁴).
        if (f0_hat > 0.0 && A_hat > 0.0) {
            double sw = 0.0, swx = 0.0, swxx = 0.0;
            for (const auto& b : live) {
                const double x = std::exp(-b.mid / best_tau);
                sw   += b.w;
                swx  += b.w * x;
                swxx += b.w * x * x;
            }
            const double det4  = sw * swxx - swx * swx;
            const double tot2  = total_amp * total_amp;
            const double var_r = (f0_hat * f0_hat * sw
                                + 2.0 * A_hat * f0_hat * swx
                                + A_hat * A_hat * swxx)
                               / (det4 * tot2 * tot2);
            if (var_r >= 0.0) {
                const double se_r = std::sqrt(var_r);
                out.overhang_lo = std::max(0.0, out.overhang_fraction - 1.96 * se_r);
                out.overhang_hi = std::min(1.0, out.overhang_fraction + 1.96 * se_r);
            }
        }
        // When projected to boundary: f0=0 → overhang_fraction=1 (pure overhang, no CI needed for gate);
        //                             A=0  → overhang_fraction=0 (pure floor, gate will correctly fail).
    }

    // Secondary gate upgrade: promote UNDETERMINED → DETECTED when overhang is dominant (τ was pushed
    // past 35 by the 1-param model's f0 bias; now that f0 is separated, high overhang_fraction is the
    // reliable signal even when τ is moderate).
    if (out.state == DamageConfidence::UNDETERMINED &&
        out.overhang_fraction >= cfg.overhang_fraction_min) {
        out.state = DamageConfidence::DETECTED;
    }

    return out;
}

} // namespace taph
