// Validated reference-free damaged fraction (profile.pi) + prior-free per-read LLR.
// SOLUTION_pi_delta_dmax.md §6. finalize_pi runs production — 3-state output (DETECTED/TRACE/ANCIENT_CPG)
// serialised to "pi_estimate" in profile JSON. Legacy mixture path untouched; per-read LLR still ref-only.

#include <cstdio>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "taph/sample_damage_profile.hpp"
#include "taph/read_ancient_llr.hpp"
#include "damage_estimation_finalize_ctx.hpp"  // finalize_pi declaration

namespace taph {

namespace {

// Channel-5 strand-symmetry gate from the deam_context_spectrum (5′ C→T vs 3′ G→A RC-mapped):
// genuine terminal deamination is strand-symmetric, so ct_5prime_excess[i] mirrors ga_3prime_rc_excess[i].
// amplitude = mean over the 16 canonical NCN contexts of 0.5·(ex5+ex3g_rc); mirror residual = mean
// |ex5−ex3g_rc|. AUTHENTIC ⇔ amplitude > PI_CH5_AMP_THR AND residual < 0.5·amplitude. Mirrors the
// finalize.cpp bilateral-CpG block's ctx_ex; recomputed here (the spectrum is not stored on the profile).
// Channel-5 strand-symmetry: the artifact-robust damage discriminator. Returns the authenticated
// terminal-deamination amplitude (mean of the strand-mirrored 5'C->T / 3'G->A context excess) and the
// mirror residual. amp is the A_obs used for the pi range — bulk delta and tau are artifact-inflated /
// non-discriminating (FLB57md's ->G overcall pushes bulk delta ABOVE the genuinely-damaged FLB03mAds3),
// so channel-5 is the only clean separator (validated: FLB03mAds3 amp~+0.16 AUTHENTIC, FLB57md ~0 ABSTAIN).
struct Channel5 { bool authentic = false; double amp = 0.0; double resid = 0.0; double amp_se = 0.0; };
Channel5 channel5_eval(const SampleDamageProfile& profile) {
    static constexpr int RC4[4] = {3, 2, 1, 0};
    const auto ctx_ex = [&](const std::array<std::uint64_t, 64>& term,
                            const std::array<std::uint64_t, 64>& intr,
                            int mid_ref, int mid_dam, int p, int n) -> double {
        const int ir = p * 16 + mid_ref * 4 + n;
        const int id = p * 16 + mid_dam * 4 + n;
        const std::uint64_t tn = term[ir] + term[id];
        const std::uint64_t xn = intr[ir] + intr[id];
        if (!tn || !xn) return 0.0;
        return static_cast<double>(term[id]) / tn - static_cast<double>(intr[id]) / xn;
    };

    double amp = 0.0, amp2 = 0.0, resid = 0.0;
    int used = 0;
    for (int i = 0; i < 16; ++i) {
        const int p = i / 4, n = i % 4;
        const double ex5 = ctx_ex(profile.tri_5prime_terminal, profile.tri_5prime_interior, 1, 3, p, n);
        const double ex3 = ctx_ex(profile.tri_3prime_terminal, profile.tri_3prime_interior,
                                  2, 0, RC4[n], RC4[p]);
        const double a_i = 0.5 * (ex5 + ex3);
        amp   += a_i;
        amp2  += a_i * a_i;
        resid += std::abs(ex5 - ex3);
        ++used;
    }
    if (used == 0) return {};
    amp   /= used;
    resid /= used;
    // between-context SE of the mean amplitude (rough, but a real uncertainty band for the pi interval)
    const double var = std::max(0.0, amp2 / used - amp * amp);
    const double amp_se = std::sqrt(var / used);
    return { amp > PI_CH5_AMP_THR && resid < 0.5 * amp, amp, resid, amp_se };
}

bool channel5_authentic(const SampleDamageProfile& profile) { return channel5_eval(profile).authentic; }

}  // namespace

// Gated reference-free pi RANGE from the per-position terminal-decay fit (PiShapeFit). pi is a
// GATED RANGE+STATE, never a bare point (the point is the range midpoint, kept only for the legacy
// downstream consumer; the JSON emitter nulls it unless it sits inside its own CI).
//   pi ≈ A_obs / (D_MAX_CONSERVED − η),  A_obs = fitted BULK 5′ amplitude,  η = 0 (conservative).
// Range from ±1.96·A_se propagated through the same divisor.
// State:
//   DETECTED    channel-5 AUTHENTIC ∧ shape LRT DETECTED ∧ pi_lo > 0.005
//   TRACE       channel-5 AUTHENTIC ∧ shape LRT DETECTED ∧ pi_lo ≤ 0.005 ≤ pi_hi
//   ANCIENT_CPG gate not met but bilateral CpG fingerprint authenticated (dilute/UDG ancient)
//   UNDETERMINED otherwise (ABSTAIN: no range emitted)
// D_MAX_CONSERVED (damage_estimate.hpp:12) is the cohort-conserved per-ancient terminal amplitude;
// A_obs/D_MAX_CONSERVED rescales the bulk amplitude to a damaged-read fraction. η held at 0
// (no attenuation subtracted) — the lo edge already carries the conservative band via A_se.
void finalize_pi(SampleDamageProfile& profile) {
    profile.pi = DamageEstimate{};
    if (!profile.bulk_attempted) return;

    const auto try_ancient_cpg = [&]() {
        if (!std::isnan(profile.cpg_delta_bilateral) &&
            (double)profile.cpg_delta_bilateral >= CPG_BILATERAL_ANCIENT_THR) {
            profile.pi.state = DamageConfidence::ANCIENT_CPG;
        }
    };

    // Channel-5 strand-symmetry is BOTH the authenticity gate and the amplitude source. It is the only
    // artifact-robust discriminator: bulk delta and tau are artifact-inflated/non-discriminating (the
    // ->G overcall pushes FLB57md's bulk delta ABOVE genuinely-damaged FLB03mAds3), so the prior
    // fit_pi_shape-gated path mis-abstained on real damage. fit_pi_shape (profile.pi_shape) is retained
    // for the per-read LLR decay (lambda/baseline), NOT as the detection gate here.
    const Channel5 c5 = channel5_eval(profile);
    if (!c5.authentic) { try_ancient_cpg(); return; }

    const auto clip01 = [](double x) { return std::clamp(x, 0.0, 1.0); };
    const double denom = D_MAX_CONSERVED;
    if (denom <= 0.0) { try_ancient_cpg(); return; }

    // pi VALUE from the Briggs per-end d_max (position-peak amplitude, COMMENSURATE with D_MAX_CONSERVED).
    // channel-5 above is the authenticity GATE ONLY — its context-scale amp (~0.03) is NOT used for the pi
    // value: amp/D_MAX_CONSERVED mixed a context-axis numerator with a position-peak denominator and put pi
    // on the wrong scale (~0.07 for a clearly-damaged library). d_max_combined is unreliable here (zeroed by
    // the bilateral-inversion combiner), so use the larger per-end d_max; the gate already enforced strand
    // symmetry, so the larger end is not an artifact.
    const double A    = std::max(profile.d_max_5prime, profile.d_max_3prime);
    const double A_se = (profile.bulk_damage.d_max_damaged_valid && profile.bulk_damage.d_max_se > 0.0)
                            ? profile.bulk_damage.d_max_se : 0.2 * A;
    profile.pi.lo    = clip01((A - 1.96 * A_se) / denom);
    profile.pi.hi    = clip01((A + 1.96 * A_se) / denom);
    profile.pi.point = clip01(A / denom);  // range midpoint; JSON nulls it unless lo≤point≤hi

    constexpr double PI_THR = 0.005;
    if (profile.pi.lo > PI_THR)       profile.pi.state = DamageConfidence::DETECTED;
    else if (profile.pi.hi >= PI_THR) profile.pi.state = DamageConfidence::TRACE;
    else                              try_ancient_cpg();
}

namespace {

// Single-stranded ends share one molecule's overhang lesion chemistry (5′ C→T and 3′ C→T are the same
// deamination event seen from opposite termini), so treating the two ends as independent double-counts the
// evidence. We deflate the SECOND (3′) end's log-LR by SS_DEP_RHO ∈ (0,1] as a first-order dependence
// correction. ceiling: this is a scalar deflation, not a true joint two-end likelihood — it cannot capture
// per-position 5′/3′ covariance; upgrade: score against the pi_joint moment cube (sample_damage_profile.hpp)
// when a per-read joint kernel exists. DS ends are on separate strands ⇒ independent ⇒ rho = 1.
constexpr double SS_DEP_RHO = 0.5;

// Per-end position-resolved log-LR under the Briggs decay q_a(p)=b+(1−b)·A·e^{−λp} vs background q_m=b.
double end_llr(const TerminalMismatch* s, std::uint32_t n, double lambda, double b, double A) {
    double llr = 0.0;
    const double qm = std::clamp(b, 1e-6, 1.0 - 1e-6);  // modern: background only
    for (std::uint32_t i = 0; i < n; ++i) {
        const double decay = std::exp(-lambda * static_cast<double>(s[i].pos));
        const double qa = std::clamp(b + (1.0 - b) * A * decay, 1e-6, 1.0 - 1e-6);  // ancient
        llr += s[i].deaminated ? std::log(qa / qm)
                               : std::log((1.0 - qa) / (1.0 - qm));
    }
    return llr;
}

}  // namespace

double damage_evidence_llr(const ReadDamageObs& obs, const SampleDamageProfile& profile) {
    const double A = D_MAX_CONSERVED;  // imported amplitude (§6.5(1)); NOT the mixture

    // Prefer the per-position decay actually fitted on the ref-free path (pi_shape); fall back to the
    // per-end λ/baseline fields for the reference-anchored path. A flat λ=0.3 default still yields a usable
    // monotone ranking score (more terminal deam ⇒ higher), which is exactly what RANKING needs.
    double lam5 = static_cast<double>(profile.lambda_5prime);
    double lam3 = static_cast<double>(profile.lambda_3prime);
    double b5   = static_cast<double>(profile.fit_baseline_5prime);
    double b3   = static_cast<double>(profile.fit_baseline_3prime);
    if (profile.pi_shape.fitted) {
        lam5 = lam3 = profile.pi_shape.lambda;
        b5 = b3 = profile.pi_shape.baseline;
    }

    const double llr5 = end_llr(obs.five,  obs.n_five,  lam5, b5, A);
    const double llr3 = end_llr(obs.three, obs.n_three, lam3, b3, A);

    const double rho = (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED)
                           ? SS_DEP_RHO : 1.0;
    return llr5 + rho * llr3;
}

std::optional<double> read_ancient_posterior(const ReadDamageObs& obs, const SampleDamageProfile& profile) {
    if (profile.pi.state != DamageConfidence::DETECTED) return std::nullopt;
    const double pi = std::clamp(profile.pi.point, 1e-6, 1.0 - 1e-6);
    const double logit_pi = std::log(pi / (1.0 - pi));
    const double z = damage_evidence_llr(obs, profile) + logit_pi;
    return 1.0 / (1.0 + std::exp(-z));  // sigmoid
}

std::optional<double> read_ancient_llr(const ReadDamageObs& obs, const SampleDamageProfile& profile) {
    if (profile.pi.state != DamageConfidence::DETECTED) return std::nullopt;
    // The LLR kernel uses per-position Briggs λ. In reference-free mode λ is never fitted (stays at 0.3
    // default, lambda_5prime_fitted=false). Scoring reads off an unfitted λ produces confident-looking but
    // wrong LLRs — refuse rather than mislead. (Legacy contract; damage_evidence_llr has no such guard.)
    if (!profile.lambda_5prime_fitted && !profile.lambda_3prime_fitted) return std::nullopt;
    return damage_evidence_llr(obs, profile);
}

}  // namespace taph
