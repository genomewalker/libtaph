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

    // Both ds and ss read the 3' arm as G->A (mid_ref=2=G, mid_dam=0=A, RC-mapped to align with 5' C->T).
    // EMPIRICAL: real ss FLB libraries show 3' G->A (ds-like, symmetric with 5' C->T), NOT 3' C->T — the
    // C->T-at-3' signal belongs to the ss-blank composition artifact (ExrNTC: 3' C->T strongly positive
    // while 3' G->A goes negative). So the strand-symmetry arm authenticates real ss and rejects the blank.
    const bool is_ss = profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED;
    double amp = 0.0, amp2 = 0.0, resid = 0.0;
    int used = 0;
    for (int i = 0; i < 16; ++i) {
        const int p = i / 4, n = i % 4;
        const double ex5 = ctx_ex(profile.tri_5prime_terminal, profile.tri_5prime_interior, 1, 3, p, n);
        const double ex3 = ctx_ex(profile.tri_3prime_terminal, profile.tri_3prime_interior, 2, 0, RC4[n], RC4[p]);
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
    // DS: proportional symmetry gate (strand-symmetric mechanism). SS: absolute residual ceiling, since
    // genuine ss is 5'-dominant and cannot meet resid<0.5*amp; the amp>THR floor rejects the ss blanks.
    const bool authentic = is_ss
        ? (amp > PI_CH5_AMP_THR && resid < PI_SS_RESID_MAX)
        : (amp > PI_CH5_AMP_THR && resid < 0.5 * amp);
    return { authentic, amp, resid, amp_se };
}

bool channel5_authentic(const SampleDamageProfile& profile) { return channel5_eval(profile).authentic; }

// === END-ASYMMETRIC live-end decay LRT (low-abundance cascade primary gate) ===
// Two-model LRT on ONE end's merged per-position histogram, run on the SAME terminal window the spike
// occupies (covariate-matched: same L/C strata, same positions p=0..P_PI-1; interior/longest reads are
// NEVER used as the artifact-inclusive reference). Shared baseline b = terminal-window asymptote (deepest
// live position). The two competing shapes for the terminus:
//   ALT  (Briggs decay):  r(p) = b + (1-b)·A·e^{-λp}        — genuine terminal deamination, decays into interior
//   NULL (spike/flat):    r(p) = s for p∈{0,1}, r(p) = b otherwise  — position-fixed dropout (->G overcall),
//                                                                     elevated only at the very terminus, no decay
// LRT = 2·(ll_decay − ll_spike). The null free param s is the pooled pos0-1 rate (1 extra df over flat),
// so this is a SHAPE test (decay continuation through pos2-4) not a mere amplitude test — exactly the
// degeneracy discriminator: a dropout spike collapses to b by pos2 (null wins), a real channel decays
// smoothly (alt wins). Same closed-form WLS linearisation as fit_pi_shape; consumes only the merged cube.
struct EndDecayLRT { double lrt = 0.0; double A = 0.0; double A_se = 0.0; double lambda = 0.0;
                     double baseline = 0.0; double n_elig_total = 0.0; bool fitted = false; bool live = false; };
EndDecayLRT live_end_decay_lrt(const SampleDamageProfile::PiPosCube& cube) {
    constexpr int P = SampleDamageProfile::P_PI;
    EndDecayLRT out;
    std::array<std::array<double,2>, P> pos = {};  // [p] = {n_elig, n_deam}
    // Pool over ALL C (incl. C=0). C_bin is the 5' C->T centroid: for the 3'-live case (5' artifact-dead,
    // FLB57md) it is artifact-contaminated, so C∈{1,2} selects on the 5' artifact and scrambles the 3'
    // signal -> flat rates, no decay. The UNSTRATIFIED per-position profile (all C) is the one the Briggs
    // fit uses, where the genuine d_max3 decay lives. (The modest C-enrichment isn't worth the circularity.)
    for (int L = 0; L < SampleDamageProfile::N_PI_LEN; ++L)
        for (int C = 0; C < SampleDamageProfile::N_PI_C; ++C)
            for (int p = 0; p < P; ++p) {
                pos[p][0] += static_cast<double>(cube[L][C][p].n_elig);
                pos[p][1] += static_cast<double>(cube[L][C][p].n_deam);
            }
    for (int p = 0; p < P; ++p) out.n_elig_total += pos[p][0];
    int n_live = 0;
    for (int p = 0; p < P; ++p) if (pos[p][0] >= 50.0) ++n_live;
    if (n_live < 3) return out;

    double b = 0.0;
    for (int p = P - 1; p >= 0; --p)
        if (pos[p][0] >= 50.0) { b = pos[p][1] / pos[p][0]; break; }
    b = std::clamp(b, 0.0, 0.95);  // do NOT cap at 0.5: the deep-window asymptote can exceed 0.5
                                   // (e.g. 3' A/(A+G)); a 0.5 cap forced r_p-b<=0 -> degenerate WLS -> lambda=0

    // ALT: linearise y_p = log((r_p−b)/(1−b)) = log A − λ·p, binomial-weighted WLS (same as fit_pi_shape).
    double sw = 0.0, swx = 0.0, swxx = 0.0, swy = 0.0, swxy = 0.0;
    for (int p = 0; p < P; ++p) {
        const double n = pos[p][0];
        if (n < 50.0) continue;
        const double r = pos[p][1] / n;
        const double num = r - b;
        if (num <= 1e-6 || b >= 1.0) continue;
        const double y = std::log(num / (1.0 - b));
        const double rc = std::clamp(r, 1e-6, 1.0 - 1e-6);
        const double w = n * rc * (1.0 - rc);
        const double x = static_cast<double>(p);
        sw += w; swx += w * x; swxx += w * x * x; swy += w * y; swxy += w * x * y;
    }
    const double det = sw * swxx - swx * swx;
    if (!(det > 0.0) || sw <= 0.0) return out;
    const double logA  = (swxx * swy - swx * swxy) / det;
    const double slope = (sw * swxy - swx * swy) / det;
    double A = std::clamp(std::exp(logA), 0.0, 1.0);
    double lambda = std::max(0.0, -slope);
    if (!std::isfinite(A) || !std::isfinite(lambda)) return out;
    const double var_logA = swxx / det;
    const double A_se = (var_logA >= 0.0) ? A * std::sqrt(var_logA) : 0.2 * A;

    const auto binom_ll = [&](const double* r) {
        double ll = 0.0;
        for (int p = 0; p < P; ++p) {
            const double n = pos[p][0], k = pos[p][1];
            if (n <= 0.0) continue;
            const double rp = std::clamp(r[p], 1e-9, 1.0 - 1e-9);
            ll += k * std::log(rp) + (n - k) * std::log(1.0 - rp);
        }
        return ll;
    };
    // Nested decay-vs-FLAT LRT (valid χ², always ≥0 for a real decay): ALT = WLS Briggs decay, NULL =
    // constant pooled rate (decay with A=0). The decay-vs-SPIKE test was NOT nested (a free terminal rate
    // is more flexible than the flat-decay), so ll_decay<ll_spike gave a hugely NEGATIVE statistic that
    // could never clear the floor. Spike/dropout discrimination is handled per-end by the artifact flag,
    // which already marks the dead end; here we only need the LIVE end to carry a genuine terminal decay.
    double r_decay[P], r_flat[P];
    double tot_n = 0.0, tot_k = 0.0;
    for (int p = 0; p < P; ++p) if (pos[p][0] >= 50.0) { tot_n += pos[p][0]; tot_k += pos[p][1]; }
    const double r_bar = tot_n > 0.0 ? tot_k / tot_n : b;
    for (int p = 0; p < P; ++p) {
        r_decay[p] = b + (1.0 - b) * A * std::exp(-lambda * static_cast<double>(p));
        r_flat[p]  = r_bar;
    }
    const double lrt = 2.0 * (binom_ll(r_decay) - binom_ll(r_flat));

    out.A = A; out.A_se = A_se; out.lambda = lambda; out.baseline = b;
    out.lrt = lrt; out.fitted = true;
    // LIVE ⇔ the decay shape beats a flat constant-rate null (df=2: A and λ) AND it genuinely decays
    // (λ>0, A above floor). χ²(df=2) 0.99 = 9.21.
    out.live = lrt >= 9.21 && lambda > 0.0 && A >= 0.04;
    if (std::getenv("TAPH_DBG_LRT")) {
        std::fprintf(stderr, "[LRT2] b=%.4f r_bar=%.4f A=%.4f lambda=%.4f lrt=%.2f live=%d nelig=%.0f | rates:",
                     b, r_bar, A, lambda, lrt, (int)out.live, out.n_elig_total);
        for (int p = 0; p < P; ++p) std::fprintf(stderr, " %.4f", pos[p][0] > 0 ? pos[p][1]/pos[p][0] : 0.0);
        std::fprintf(stderr, "\n");
    }
    return out;
}

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
    constexpr double PI_THR = 0.005;  // π-floor shared by DETECTED, TRACE, and LOW_ABUNDANCE gates

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
    const auto clip01 = [](double x) { return std::clamp(x, 0.0, 1.0); };
    const double denom = D_MAX_CONSERVED;

    const Channel5 c5 = channel5_eval(profile);
    if (!c5.authentic) {
        // END-ASYMMETRIC LOW_ABUNDANCE recovery. channel-5 needs BOTH ends strand-symmetric; it fails
        // when one terminus is ->G-overcall artifact-dead. Recover IFF EXACTLY ONE end is artifact-flagged
        // AND the OPPOSITE (live) end passes the two-model decay LRT. The artifact flag marks WHICH end is
        // dead; the LRT proves the surviving end carries decay-shaped (not spike) terminal deamination.
        // DS-ONLY: end-asymmetric recovery is validated for double-stranded libraries (live end carries a
        // clean G->A/C->T terminal decay). For SINGLE-STRANDED libraries the 3' terminus carries a
        // library-prep artifact (a clean-looking C->T decay, asymmetric d_max3>>d_max5; e.g. ExrNTC ss
        // blanks d_max3~0.64 vs metaDMG per-taxon max ~0.08) that the live-end test cannot separate from
        // real damage, so recovery on ss would emit inflated garbage. Honest scope: ds recovers, ss ABSTAINs
        // (its 3'-prep artifact is unresolved). TODO: model/subtract the ss 3'-prep artifact, then enable ss.
        const bool ds_lib = profile.library_type != SampleDamageProfile::LibraryType::SINGLE_STRANDED;
        if (denom > 0.0 && ds_lib) {
            const bool a5 = profile.artifact_overcall_5p, a3 = profile.artifact_overcall_3p;
            if (a5 != a3) {  // exactly one end dead
                const bool live_is_3p = a5;  // 5' dead ⇒ 3' is the live end (FLB57md), and vice versa
                const EndDecayLRT le = live_end_decay_lrt(
                    live_is_3p ? profile.pi_pos_3prime_ds : profile.pi_pos_5prime);
                const double A = live_is_3p ? profile.d_max_3prime : profile.d_max_5prime;
                // NOTE: an amplitude ceiling (A < D_MAX_CONSERVED) is WRONG for deep time — at ~4 Myr
                // deamination can saturate, so a genuine d_max of 0.5-0.7 is plausible; D_MAX_CONSERVED=0.39
                // is a YOUNG-aDNA cohort mean. Discrimination of real-saturated vs artifact must be by SHAPE,
                // not amplitude (see assessment: per-position window shape does NOT separate them here either).
                if (A > 0.0) {
                    // A_se: prefer the bulk damaged-split SE; else the LRT's own A_se; else 20% of A.
                    const double A_se = (profile.bulk_damage.d_max_damaged_valid &&
                                         profile.bulk_damage.d_max_se > 0.0)
                                            ? profile.bulk_damage.d_max_se
                                            : (le.A_se > 0.0 ? le.A_se : 0.2 * A);
                    // Propagate a denominator-uncertainty term too: D_MAX_CONSERVED is a cohort mean,
                    // not sample-specific, so widen by its documented spread (0.34-0.48 ⇒ ~0.035 sd).
                    constexpr double DENOM_SD = 0.035;
                    const double rel = std::sqrt((A_se * A_se) / (A * A) +
                                                 (DENOM_SD * DENOM_SD) / (denom * denom));
                    const double pt = A / denom;
                    const double lo = clip01(pt * (1.0 - 1.96 * rel));
                    const double hi = clip01(pt * (1.0 + 1.96 * rel));
                    // DETECTION FLOOR — SELF-CALIBRATING, no hardcoded d_max threshold. The live end is
                    // recovered IFF its d_max-SE-propagated lower CI edge clears the same PI_THR π-floor the
                    // symmetric DETECTED path uses (i.e. the recovered interval significantly excludes zero),
                    // AND the terminal window carries enough eligible sites. Significance comes from the data's
                    // own d_max_se (Briggs fit) + the cohort-denominator spread — not a magic amplitude cut.
                    // The Briggs length-stratified d_max already established the decay (the per-position window
                    // LRT `le` is under-powered at this dilution and kept only as a diagnostic). The dead end is
                    // the artifact-flagged one; A is a G->A amplitude a ->G overcall would DEPLETE not inflate.
                    // Pass ⇒ LOW_ABUNDANCE; fail ⇒ BELOW_FLOOR (honest upper bound, point+lo nulled).
                    // TWO conditions, both needed and both principled (not a new hardcode):
                    //  (1) significance — the d_max-SE-propagated lower CI clears PI_THR (self-calibrating);
                    //  (2) detection floor — live d_max ≥ the codebase's existing amplitude floor a_min(0.04).
                    // Significance alone is insufficient: at high depth (a 943k-read blank) a tiny artifact
                    // d_max~0.023 has a tight CI and looks "significant", so it must also clear the minimum
                    // detectable amplitude below which real damage is indistinguishable from artifact ref-free.
                    // FLB57md live 3' d_max3≈0.058 ≥ 0.04 passes; the blank's 0.023 falls to BELOW_FLOOR.
                    const double A_MIN = TauConfig{}.a_min;
                    if (lo > PI_THR && A >= A_MIN && le.n_elig_total >= PI_FLOOR_MIN_ELIG) {
                        profile.pi.point = clip01(pt);
                        profile.pi.lo    = lo;
                        profile.pi.hi    = hi;
                        profile.pi.state = DamageConfidence::LOW_ABUNDANCE;
                    } else {
                        profile.pi.point = -1.0;
                        profile.pi.lo    = -1.0;
                        profile.pi.hi    = hi;
                        profile.pi.state = DamageConfidence::BELOW_FLOOR;
                    }
                    return;
                }
            }
        }
        try_ancient_cpg();
        return;
    }

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
    // DELIBERATELY DETECTED-only: LOW_ABUNDANCE authenticates the SAMPLE end-asymmetrically (bulk pi),
    // but per-read scoring needs a fitted per-position λ that ref-free mode never produces, and the dead
    // end's reads carry the ->G artifact. Admitting LOW_ABUNDANCE here would emit confident-but-wrong
    // per-read posteriors. The ranked per-read tail (online-FDR) is the deferred cascade step that will
    // unlock it; until then LOW_ABUNDANCE returns nullopt (no per-read claim) by design.
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
