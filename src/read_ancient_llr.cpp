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
//   pi ≈ A_obs / (denom − η),  A_obs = fitted BULK 5′ amplitude,  η = 0 (conservative).
// Range from ±1.96·A_se propagated through the same divisor.
// State:
//   DETECTED    channel-5 AUTHENTIC ∧ shape LRT DETECTED ∧ pi_lo > 0.005
//   TRACE       channel-5 AUTHENTIC ∧ shape LRT DETECTED ∧ pi_lo ≤ 0.005 ≤ pi_hi
//   ANCIENT_CPG gate not met but bilateral CpG fingerprint authenticated (dilute/UDG ancient)
//   UNDETERMINED otherwise (ABSTAIN: no range emitted)
// denom is the library's OWN co-occurrence-fished damaged-class terminal amplitude (ev.cooc_dmax_point,
// mirroring finalize_pi's profile.d_max_cooccurrence.point — no hardcoded cohort constant); A_obs/denom
// rescales the bulk amplitude to a damaged-read fraction. η held at 0 (no attenuation subtracted) — the
// lo edge already carries the conservative band via A_se. denom ≤ 0 ⇒ non-identifiable ⇒ abstain.
// Pure byte-identical copy of finalize_pi's body that returns a PiVerdict from a DamageEvidence view
// instead of mutating profile.pi. NOT called yet — added as the consolidation seam so finalize_pi can
// later delegate to it. Control flow / arithmetic mirror finalize_pi line-for-line: the only changes are
// profile.<field> -> ev.<field>, profile.bulk_damage.<f> -> ev.bulk_d_max_<f>, channel5_eval(profile)
// -> the precomputed channel5_* scalars, profile.pi_pos_* -> *ev.pi_pos_*, and profile.pi -> the
// returned PiVerdict v.
PiVerdict combine(const DamageEvidence& ev) {
    PiVerdict v;
    if (!ev.bulk_attempted) return v;
    constexpr double PI_THR = 0.005;  // π-floor shared by DETECTED, TRACE, and LOW_ABUNDANCE gates

    const auto try_ancient_cpg = [&]() {
        if (!std::isnan(ev.cpg_delta_bilateral) &&
            (double)ev.cpg_delta_bilateral >= CPG_BILATERAL_ANCIENT_THR) {
            v.state = DamageConfidence::ANCIENT_CPG;
        }
    };

    const auto clip01 = [](double x) { return std::clamp(x, 0.0, 1.0); };
    const double denom = (ev.cooc_dmax_point > 0.0) ? ev.cooc_dmax_point : -1.0;

    const Channel5 c5 = { ev.channel5_authentic, ev.channel5_amp, ev.channel5_resid, ev.channel5_amp_se };
    if (!c5.authentic) {
        const bool ds_lib = ev.library_type == SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
        if (denom > 0.0 && ds_lib) {
            const bool a5 = ev.artifact_overcall_5p, a3 = ev.artifact_overcall_3p;
            if (a5 != a3) {  // exactly one end dead
                const bool live_is_3p = a5;  // 5' dead ⇒ 3' is the live end (FLB57md), and vice versa
                const EndDecayLRT le = live_end_decay_lrt(
                    live_is_3p ? *ev.pi_pos_3prime_ds : *ev.pi_pos_5prime);
                const double A = live_is_3p ? ev.d_max_3prime : ev.d_max_5prime;
                if (A > 0.0) {
                    const double A_se = (ev.bulk_d_max_damaged_valid &&
                                         ev.bulk_d_max_se > 0.0)
                                            ? ev.bulk_d_max_se
                                            : (le.A_se > 0.0 ? le.A_se : 0.2 * A);
                    constexpr double DENOM_SD = 0.035;
                    const double rel = std::sqrt((A_se * A_se) / (A * A) +
                                                 (DENOM_SD * DENOM_SD) / (denom * denom));
                    const double pt = A / denom;
                    const double lo = clip01(pt * (1.0 - 1.96 * rel));
                    const double hi = clip01(pt * (1.0 + 1.96 * rel));
                    const double A_MIN = TauConfig{}.a_min;
                    if (lo > PI_THR && A >= A_MIN && le.n_elig_total >= PI_FLOOR_MIN_ELIG) {
                        v.point = clip01(pt);
                        v.lo    = lo;
                        v.hi    = hi;
                        v.state = DamageConfidence::LOW_ABUNDANCE;
                    } else {
                        v.point = -1.0;
                        v.lo    = -1.0;
                        v.hi    = hi;
                        v.state = DamageConfidence::BELOW_FLOOR;
                    }
                    return v;
                }
            }
        }
        try_ancient_cpg();
        return v;
    }

    if (denom <= 0.0) { try_ancient_cpg(); return v; }

    const double A    = std::max(ev.d_max_5prime, ev.d_max_3prime);
    const double A_se = (ev.bulk_d_max_damaged_valid && ev.bulk_d_max_se > 0.0)
                            ? ev.bulk_d_max_se : 0.2 * A;
    v.lo    = clip01((A - 1.96 * A_se) / denom);
    v.hi    = clip01((A + 1.96 * A_se) / denom);
    v.point = clip01(A / denom);  // range midpoint; JSON nulls it unless lo≤point≤hi

    if (v.lo > PI_THR)       v.state = DamageConfidence::DETECTED;
    else if (v.hi >= PI_THR) v.state = DamageConfidence::TRACE;
    else                     try_ancient_cpg();
    return v;
}

// ── Within-read co-occurrence π (centered-covariance estimator) ──────────────────────────────────
// Identifies the ancient fraction π SEPARATELY from the per-lineage damage amplitude by exploiting
// that the two read ends co-vary only through the shared latent ancient state. Per length bin L,
// terminal pooled over C vs the interior null:
//   x_L = (μ5T−μ5I)(μ3T−μ3I)  [excess-marginal product ∝ π²δ5δ3]
//   y_L = covT − covI          [excess covariance         ∝ π(1−π)δ5δ3]
// ⇒ y_L/x_L = (1−π)/π = α, π = 1/(1+α). Global zero-intercept inverse-variance fit y_L = α·x_L.
// GATE BEFORE the ratio: a blank has x,y→0 so the raw quotient is 0/0 (noise → π garbage); emit a
// point ONLY when the aggregate excess covariance is significant (z≥2.58) and both excess marginals
// are >0, else ABSTAIN. ds-only; DIAGNOSTIC — never touches profile.pi, d_max, or mixture_model.
static void estimate_pi_cooccurrence(SampleDamageProfile& profile) {
    profile.pi_cooccurrence   = DamageEstimate{};
    profile.d_max_cooccurrence = DamageEstimate{};
    profile.pi_cooccurrence_lift   = -1.0;
    profile.pi_cooccurrence_lift_z = -1.0;
    // PREP PRECEDENCE: the --library-type FORCED call (profile.library_type) is AUTHORITATIVE and already
    // drove the 3' allele counted into pi_cooc at accumulate time (#A ds / #T ss). The data-driven detector
    // below (pooled 3' A-excess vs T-excess over the shuffle null) is a QC CROSS-CHECK only — it can itself
    // misclassify on dilute libraries where the 3' excess is at the noise floor, so it never overrides the
    // arg. A forced/detected mismatch is surfaced as a WARN (wrong metadata or an odd/contaminated library).
    // forced = the user's --library-type arg (forced_library_type; UNKNOWN ⇒ "auto"). The channel that
    // actually drove pi_cooc is profile.library_type (= the forced type when set, else the BIC call).
    const char* forced = profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED ? "ss"
                       : profile.forced_library_type == SampleDamageProfile::LibraryType::DOUBLE_STRANDED ? "ds" : "auto";
    const double qtn = static_cast<double>(profile.cooc_qc_term_n);
    const double qnn = static_cast<double>(profile.cooc_qc_null_n);
    const double a3_excess = (qtn > 0.0 && qnn > 0.0)
        ? static_cast<double>(profile.cooc_qc_term_a3) / qtn - static_cast<double>(profile.cooc_qc_null_a3) / qnn : 0.0;
    const double t3_excess = (qtn > 0.0 && qnn > 0.0)
        ? static_cast<double>(profile.cooc_qc_term_t3) / qtn - static_cast<double>(profile.cooc_qc_null_t3) / qnn : 0.0;
    constexpr double PREP_MIN = 0.01;   // 3' marginal-excess floor to call a detected prep
    const char* detected = (std::max(a3_excess, t3_excess) < PREP_MIN) ? "unknown"
                         : (a3_excess >= t3_excess ? "ds" : "ss");
    // forced ∈ {ds,ss,auto}, detected ∈ {ds,ss,unknown}; auto (no explicit arg) never disagrees.
    const bool agree = profile.forced_library_type == SampleDamageProfile::LibraryType::UNKNOWN
                     || forced[0] == detected[0];   // 'd'/'s' distinguishes ds/ss; 'u'(unknown) matches neither
    if (std::getenv("DEBUG_PICOOC")) {
        std::fprintf(stderr, "PICOOC PREP forced=%s detected=%s a3_excess=%.4f t3_excess=%.4f agree=%s\n",
                     forced, detected, a3_excess, t3_excess, agree ? "yes" : "no");
        if (!agree)
            std::fprintf(stderr, "PICOOC PREP WARN forced=%s but detected=%s (3' geometry disagrees with "
                         "--library-type; forced wins the computation — check library metadata)\n",
                         forced, detected);
    }

    constexpr int P = SampleDamageProfile::P_PI;
    constexpr double PI_COOC_THR = 0.005;
    constexpr uint64_t MIN_N = 200;   // per-bin read floor for a stable second moment
    // ss deaminates C->T at BOTH termini, so the 5'-T and 3'-T co-occurrence counts share ONE base (T):
    // within a fixed-length read they compete for the same finite T budget, making the composition
    // covariance GENUINELY NEGATIVE at both the terminal (covT) and the shuffle null (covI). That negative
    // pedestal must be SUBTRACTED (covT − covI), not clipped: flooring covI→0 deletes the cancellation and
    // drives y = covT − 0 < 0 on real ss damage. ds counts DIFFERENT 3' bases (5'-T vs 3'-A) so its
    // pedestal sign assumption (covI ≥ 0, covT > 0) is a ds-only premise. profile.library_type is the
    // channel that actually drove pi_cooc accumulation (see PREP PRECEDENCE note above).
    const bool is_ss = profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED;

    auto moments = [](const SampleDamageProfile::PiCooc& c,
                      double& mu5, double& mu3, double& cov, double& var_cov) {
        const double n = static_cast<double>(c.n);
        mu5 = static_cast<double>(c.sum_k5) / n;
        mu3 = static_cast<double>(c.sum_k3) / n;
        cov = static_cast<double>(c.sum_k5k3) / n - mu5 * mu3;      // E[k5·k3] − μ5μ3
        double eg2 = 0.0;                                           // E[g²], g = (k5−μ5)(k3−μ3)
        for (int a = 0; a <= P; ++a)
            for (int b = 0; b <= P; ++b) {
                const double h = static_cast<double>(c.hist[a * (P + 1) + b]);
                if (h == 0.0) continue;
                const double g = (static_cast<double>(a) - mu5) * (static_cast<double>(b) - mu3);
                eg2 += h * g * g;
            }
        eg2 /= n;
        var_cov = std::max(0.0, (eg2 - cov * cov) / n);             // Var(cov) = Var(g)/n (centered influence)
    };

    double Sxy = 0.0, Sxx = 0.0;          // inverse-variance accumulators for α, over ADMITTED bins only
    double sumY = 0.0, varSumY = 0.0;     // aggregate excess covariance + variance (ALL bins; legacy lift_z)
    double sumYadm = 0.0, varYadm = 0.0;  // aggregate excess covariance + variance over ADMITTED bins (gate)
    double ex5w = 0.0, ex3w = 0.0;        // read-weighted aggregate excess marginals (admitted bins)
    double jointT = 0.0, indepT = 0.0, nT_all = 0.0;   // raw lift diagnostic
    int used_bins = 0;

    for (int L = 0; L < SampleDamageProfile::N_PI_LEN; ++L) {
        SampleDamageProfile::PiCooc term{};           // terminal pooled over C
        for (int C = 0; C < SampleDamageProfile::N_PI_C; ++C) {
            const auto& s = profile.pi_cooc[L][C];
            term.n += s.n; term.sum_k5 += s.sum_k5; term.sum_k3 += s.sum_k3; term.sum_k5k3 += s.sum_k5k3;
            for (int h = 0; h < (P + 1) * (P + 1); ++h) term.hist[h] += s.hist[h];
        }
        const auto& intr = profile.pi_cooc_interior_ds[L];
        if (term.n < MIN_N || intr.n < MIN_N) continue;

        double mu5T, mu3T, covT, vcT, mu5I, mu3I, covI, vcI;
        moments(term, mu5T, mu3T, covT, vcT);
        moments(intr, mu5I, mu3I, covI, vcI);

        const double d5 = mu5T - mu5I;         // excess terminal 5' T (damage marginal), must be > 0
        const double d3 = mu3T - mu3I;         // excess terminal 3' A (damage marginal), must be > 0
        const double x  = d5 * d3;             // ∝ π²δ5δ3 (SIGN-BLIND product — signs gated below)
        // covI is the shuffle-null compositional pedestal. The WHOLE-read shuffle couples k5,k3 through
        // the read's global n_T/n_A base-competition (negative at fixed length), but the terminal covT
        // only sees a DILUTED version of that coupling through its two 5 bp windows, so the raw shuffle
        // covI is more negative than covT's true pedestal and would over-subtract (→ π biased low). The
        // pedestal is floored at 0 per the fix design ("covI must be ≈0 or ≥0, a clean compositional
        // pedestal"): a negative composition covariance cannot masquerade as extra damage covariance.
        // ss EXCEPTION (is_ss): both termini count T, so the pedestal is legitimately negative and must be
        // subtracted raw — flooring it deletes the base-competition cancellation and sinks y below 0.
        const double covI_eff = is_ss ? covI : std::max(0.0, covI);
        const double y  = covT - covI_eff;     // ∝ π(1−π)δ5δ3
        const double vy = vcT + vcI;
        if (!(vy > 0.0) || !std::isfinite(x) || !std::isfinite(y)) continue;

        // SIGN BACKSTOP admission: both excess marginals positive, terminal covariance positive, and the
        // excess covariance positive. Drops composition-noise bins (e.g. L=2 with d5<0,d3<0 yet x>0 via
        // the sign-blind product, and a hugely negative covI) so they can never enter — and thereby
        // deflate — the fit. A pathological bin SELF-EXCLUDES here. |x| >= X_MIN additionally drops the
        // near-zero-x blowup (bin L=2 had x=1e-6 → alpha_L=89402): a bin whose excess-marginal product is
        // consistent with 0 cannot identify α (y/x explodes) and its inverse-variance weight x²/vy is
        // negligible anyway, so excluding it changes the fit by ~0 while removing the pathology.
        constexpr double X_MIN = 1e-4;
        // covT > 0.0 is a ds sign-assumption (its pedestal is ~0/positive). For ss the terminal covT is
        // itself negative from the shared-T base competition; the damage signal lives in the pedestal-
        // subtracted excess y = covT − covI > 0, which is already required below. So ss admits on y > 0
        // without the covT > 0 term; the other backstops (d5>0, d3>0, |x|>=X_MIN) are prep-independent.
        const bool admit = d5 > 0.0 && d3 > 0.0 && (is_ss || covT > 0.0) && y > 0.0 && std::fabs(x) >= X_MIN;
        const double alpha_L = (x > 0.0) ? y / x : -1.0;   // per-bin (1−π)/π (diagnostic; fit uses Sxy/Sxx)

        // lift_z aggregates ALL usable bins (diagnostic; NOT gated on admission) so its emitted value is
        // comparable to the historical dumps.
        sumY += y; varSumY += vy;
        jointT += static_cast<double>(term.sum_k5k3);
        indepT += static_cast<double>(term.n) * mu5T * mu3T;
        nT_all += static_cast<double>(term.n);

        if (admit) {
            Sxy += x * y / vy;  Sxx += x * x / vy;    // inverse-variance weighted mean of alpha_L = y/x
            sumYadm += y; varYadm += vy;              // admitted-bin excess-cov significance (fit gate)
            ex5w += static_cast<double>(term.n) * d5;
            ex3w += static_cast<double>(term.n) * d3;
            ++used_bins;
        }
        if (std::getenv("DEBUG_PICOOC"))
            std::fprintf(stderr,
                "PICOOC L=%d nT=%llu nI=%llu mu5T=%.4f mu5I=%.4f mu3T=%.4f mu3I=%.4f "
                "covT=%.5f covI=%.5f covIeff=%.5f x=%.6f y=%.6f vy=%.3e alpha_L=%.4f admit=%d\n",
                L, (unsigned long long)term.n, (unsigned long long)intr.n,
                mu5T, mu5I, mu3T, mu3I, covT, covI, covI_eff, x, y, vy, alpha_L, (int)admit);
    }

    if (nT_all > 0.0 && indepT > 0.0)
        profile.pi_cooccurrence_lift = jointT / indepT;    // raw uncentered lift (≈1 on a blank)
    const double lift_z = (varSumY > 0.0) ? sumY / std::sqrt(varSumY) : 0.0;
    profile.pi_cooccurrence_lift_z = lift_z;                // legacy all-bin lift_z: DIAGNOSTIC ONLY, kept
                                                            // for JSON continuity; NOT the gate statistic.

    // SIGNIFICANCE GATE — re-instated on the CORRECTED covI. Fixing the null (shuffle, floored) re-signed
    // the aggregate excess-covariance z so it is discriminating again: blank lift_z≈2.2 vs ancient ≈44.6
    // (the OLD interior null had this backwards — blank 121 > ancient 98 — which is the ONLY reason it was
    // dropped; that premise is void here). Blank ABSTENTION is a hard requirement, so the significance is
    // load-bearing (the covI floor alone is NOT sufficient safety). Gate on the all-bin lift_z ≥ 3.0
    // (stricter than 2.58): rejects the blank (≤2.9) with margin, passes the ancient (44.6) trivially.
    // sig_z (admitted-bin z) is kept as a corroborating diagnostic only.
    const double sig_z = (varYadm > 0.0) ? sumYadm / std::sqrt(varYadm) : 0.0;
    constexpr double LIFT_Z_THR = 3.0;

    const double alpha = (Sxx > 0.0) ? Sxy / Sxx : -1.0;   // (1−π)/π; must be >0 to identify π∈(0,1)
    // GATE: ≥1 admitted (sign-clean) bin, identifiable α>0, positive aggregate excess marginals, AND the
    // aggregate excess covariance significant under the shuffle null (lift_z ≥ 3). Blank fails on lift_z.
    const bool gate = used_bins > 0 && Sxx > 0.0 && ex5w > 0.0 && ex3w > 0.0 && alpha > 0.0
                      && lift_z >= LIFT_Z_THR;

    if (std::getenv("DEBUG_PICOOC"))
        std::fprintf(stderr, "PICOOC AGG used_bins=%d Sxy=%.4e Sxx=%.4e alpha=%.4f ex5w=%.4e ex3w=%.4e "
                     "sig_z=%.2f lift_z=%.2f gate=%d\n",
                     used_bins, Sxy, Sxx, alpha, ex5w, ex3w, sig_z, lift_z, (int)gate);

    if (!gate) { profile.pi_cooccurrence.state = DamageConfidence::UNDETERMINED; return; }

    const double var_alpha = 1.0 / Sxx;
    const double pi_pt     = 1.0 / (1.0 + alpha);
    const double dpi_da    = -1.0 / ((1.0 + alpha) * (1.0 + alpha));     // delta method
    const double se_pi     = std::sqrt(var_alpha) * std::fabs(dpi_da);

    auto& e = profile.pi_cooccurrence;
    e.point = pi_pt;
    e.lo = std::max(0.0, pi_pt - 1.96 * se_pi);
    e.hi = std::min(1.0, pi_pt + 1.96 * se_pi);
    if (e.lo > PI_COOC_THR)       e.state = DamageConfidence::DETECTED;
    else if (e.hi >= PI_COOC_THR) e.state = DamageConfidence::TRACE;
    else { e.point = -1.0; e.lo = -1.0; e.state = DamageConfidence::BELOW_FLOOR; }  // upper bound only

    // ── Intrinsic d_max from the co-occurrence π (f_C-FREE, self-normalized) ─────────────────────────
    // The count-excess route needed the terminal ref-C fraction f_C (unrecoverable reference-free → a
    // ~1.5x bias). We instead use the SELF-NORMALIZED per-position C→T rate rate[p]=ndeam/nelig = T/(C+T)
    // among C-eligible sites, which the two-state mixture makes f_C-free and exact:
    //     rate[p] = baseline + π·d_p·(1−baseline)      (baseline = undamaged composition rate c/(c+t))
    //  ⇒  d_p     = (rate[p] − baseline) / (π·(1−baseline))
    // (1−baseline) is read off the deepest (interior-most) window positions; no f_C, no geometric sum
    // (this is a per-position rate, not a window-count total). d_max = d_p at the peak position.
    if (e.state == DamageConfidence::DETECTED || e.state == DamageConfidence::TRACE) {
        // rate[p] = pooled per-position 5' C→T rate (T/(C+T)); pooled over all L,C (matches the LRT input).
        double rate[SampleDamageProfile::P_PI] = {};
        {
            double nelig[SampleDamageProfile::P_PI] = {}, ndeam[SampleDamageProfile::P_PI] = {};
            for (int L = 0; L < SampleDamageProfile::N_PI_LEN; ++L)
                for (int C = 0; C < SampleDamageProfile::N_PI_C; ++C)
                    for (int p = 0; p < P; ++p) {
                        nelig[p] += static_cast<double>(profile.pi_pos_5prime[L][C][p].n_elig);
                        ndeam[p] += static_cast<double>(profile.pi_pos_5prime[L][C][p].n_deam);
                    }
            for (int p = 0; p < P; ++p) rate[p] = nelig[p] > 0.0 ? ndeam[p] / nelig[p] : 0.0;
        }
        // baseline = interior undamaged C→T rate = T/(C+T) among interior C-eligible sites, READ OFF THE
        // MIDDLE-THIRD census (baseline_t_freq/(baseline_t+baseline_c)). This is the true composition
        // asymptote the identity wants. rate[P-1] is NOT it: the deepest window position (p=P−1) is still
        // inside the oscillating terminal zone and sits ABOVE the interior asymptote (FLB01: rate[4]=0.521
        // vs interior 0.492), which undercounts the excess and biases d_max LOW. Fall back to rate[P-1]
        // only if the interior census is empty.
        const double bt = static_cast<double>(profile.baseline_t_freq);
        const double bc = static_cast<double>(profile.baseline_c_freq);
        const double baseline = (bt + bc > 0.0) ? bt / (bt + bc) : rate[P - 1];
        // rate_peak = MAX over the terminal window: robust to the odd/even oscillation AND to the
        // adapter-suppressed pos-0 (which can sit BELOW baseline).
        double rate_peak = rate[0];
        for (int p = 1; p < P; ++p) rate_peak = std::max(rate_peak, rate[p]);

        // λ retained as a DIAGNOSTIC only (no longer used in the d_max computation).
        double lambda = 0.25;
        if (rate[0] - baseline > 1e-6) {
            const double r = (rate[1] - baseline) / (rate[0] - baseline);
            if (r > 0.0 && r < 1.0) lambda = r;
        }

        // HONEST NUMERATOR NOISE: the excess (rate_peak − baseline) sits on top of a position-to-position
        // COMPOSITIONAL wobble. Estimate its sigma from the scatter of the 5' C→T rate across the FLAT
        // INTERIOR positions (past the terminal damage/fill-in oscillation), read off the 15-position C→T
        // census. This is the systematic floor the peak must clear — the binomial SE understates it ~100×
        // at these depths. The interior baseline is a huge-N census (negligible SE), so σ_excess ≈ σ_peak.
        double sigma_excess = 0.0;
        {
            constexpr int WOBBLE_START = 8;   // positions ≥8 are past the terminal oscillation → flat interior
            double v[SampleDamageProfile::N_POS]; int nw = 0; double m = 0.0;
            for (int p = WOBBLE_START; p < SampleDamageProfile::N_POS; ++p) {
                double t = 0.0, tot = 0.0;
                for (int c = 0; c < SampleDamageProfile::N_CT_CTX; ++c) {
                    t   += profile.ct_ctx_t_5prime[c][p];
                    tot += profile.ct_ctx_total_5prime[c][p];
                }
                if (tot > 0.0) { v[nw] = t / tot; m += v[nw]; ++nw; }
            }
            if (nw >= 2) {
                m /= nw; double s2 = 0.0;
                for (int i = 0; i < nw; ++i) s2 += (v[i] - m) * (v[i] - m);
                sigma_excess = std::sqrt(s2 / (nw - 1));
            }
        }

        auto& dm = profile.d_max_cooccurrence;
        const double norm = pi_pt * (1.0 - baseline);
        const double excess = rate_peak - baseline;
        constexpr double DMAX_Z = 2.0;   // excess must clear Z·σ_wobble to be a recoverable amplitude
        // RELIABILITY GATE (data-driven, NOT a π hardcode): the amplitude is recoverable only when the bulk
        // excess exceeds the composition noise. At low π the excess collapses to the wobble floor (FLB45:
        // excess≈0.0035 vs σ≈wobble) and d_max is NOT trustworthy per-sample → abstain on the point but
        // still report π. FLB01 (excess≈0.045 ≫ σ) passes. This trips on dilution WITHOUT knowing π first.
        const bool reliable = excess > 0.0 && sigma_excess > 0.0 && excess >= DMAX_Z * sigma_excess;
        if (excess > 0.0 && norm > 0.0) {
            const double d_max_raw = excess / norm;
            // Propagate BOTH noise sources: numerator wobble (dominates at low π) AND π's SE. σ_dmax/dmax =
            // sqrt((σ_excess/excess)² + (σ_π/π)²). The old CI propagated ONLY π's band, understating the true
            // uncertainty at low π by ~the amplification factor — this makes the CI honest.
            const double rel_excess = (excess > 0.0)  ? sigma_excess / excess : 0.0;
            const double rel_pi     = (pi_pt  > 0.0)  ? se_pi / pi_pt         : 0.0;
            const double sigma_dmax = d_max_raw * std::sqrt(rel_excess * rel_excess + rel_pi * rel_pi);
            dm.lo = std::clamp(d_max_raw - 1.96 * sigma_dmax, 0.0, 1.0);
            dm.hi = std::clamp(d_max_raw + 1.96 * sigma_dmax, 0.0, 1.0);
            if (reliable) {
                dm.point = std::clamp(d_max_raw, 0.0, 1.0);
                dm.state = e.state;                       // DETECTED / TRACE — a trustworthy point
            } else {
                dm.point = -1.0;                          // numerator at the noise floor → ABSTAIN on point
                dm.state = DamageConfidence::BELOW_FLOOR; // identifiable=false; the (huge) CI still emits
            }
            if (std::getenv("DEBUG_PICOOC"))
                std::fprintf(stderr,
                    "PICOOC DMAX2 rate=[%.4f %.4f %.4f %.4f %.4f] base=%.4f peak=%.4f pi=%.4f lambda=%.4f "
                    "dmax_fcfree=%.4f sig_exc=%.4f exc/sig=%.2f reliable=%d ci=[%.4f,%.4f]\n",
                    rate[0], rate[1], rate[2], rate[3], rate[4], baseline, rate_peak, pi_pt, lambda,
                    d_max_raw, sigma_excess, sigma_excess > 0.0 ? excess / sigma_excess : 0.0,
                    (int)reliable, dm.lo, dm.hi);
        } else {
            // No terminal excess above the interior baseline ⇒ no amplitude signal ⇒ ABSTAIN on d_max.
            dm.state = DamageConfidence::UNDETERMINED;
            if (std::getenv("DEBUG_PICOOC"))
                std::fprintf(stderr,
                    "PICOOC DMAX2 rate=[%.4f %.4f %.4f %.4f %.4f] base=%.4f peak=%.4f pi=%.4f lambda=%.4f "
                    "dmax_fcfree=UNDETERMINED(peak<=base)\n",
                    rate[0], rate[1], rate[2], rate[3], rate[4], baseline, rate_peak, pi_pt, lambda);
        }
    }
}

void finalize_pi(SampleDamageProfile& profile) {
    profile.pi = DamageEstimate{};
    estimate_pi_cooccurrence(profile);   // additive diagnostic; independent of the shape-fit path below
    profile.oxidation_present = false;

    // C2: OXIDATION as an INDEPENDENT damage axis — computed unconditionally here (not nested in the
    // deamination control flow) because it never overrides pi.state. Source: oxo_two_marker δβ PAIRED
    // contrast (delta_beta = beta1 G→T − beta2 C→A), populated by finalize_oxidation_comovement BEFORE
    // finalize_pi. Gate = paired-contrast validity (markers_consistent — a strand-asymmetric prep artifact
    // FAILS it, keeping modern blanks null) AND B1 SE significance (delta_beta − 1.96·delta_beta_se > 0,
    // the contrast significantly excludes the no-oxidation null). Raw beta1/stop-rate thresholds fire on
    // blanks (the modern fixture has delta_beta_z≈24 yet is a clean negative — caught by markers_consistent
    // =false); the paired δβ + consistency gate fires on the oxidised-but-deam-barren FLB45m (consistent,
    // delta_beta−1.96·se > 0) while staying null on modern. Validated against the extended-oracle fixtures.
    {
        const auto& otm = profile.oxidation_comovement;
        profile.oxidation_present =
            otm.valid && otm.markers_consistent &&
            otm.delta_beta_se > 0.0 &&
            (otm.delta_beta - 1.96 * otm.delta_beta_se) > 0.0;
    }

    if (!profile.bulk_attempted) return;
    constexpr double PI_THR = 0.005;  // π-floor shared by DETECTED, TRACE, and LOW_ABUNDANCE gates

    // Composition-immune codon DiD (compute_codon_did): strand-paired effect size min(5',3'). The
    // purine-mirror control cancels TERMINAL positional composition skew (end-repair) that Channel A
    // amplitude and the stop-codon channel both mistake for damage; real one-sided C→T/G→A deamination
    // survives. Calibrated on a 132-lib modern panel + 13-sample FLB cohort: floor +0.010 gives ~2%
    // fresh-modern FP while recovering ROCS (+0.018) and deep FLB (FLB10m/13m/16m/57m). VETO: a clearly
    // NEGATIVE paired DiD means the symmetric terminal amplitude is composition, not damage (coastal
    // water −0.04) → reject. CONFIRM: a positive paired DiD authenticates faint damage channel-5 missed.
    // ds: strand-paired min(5' C→T, 3' G→A) — both ends must agree (composition skew fails this).
    // ss: 5' ONLY — the ss 3' terminus carries a library-prep artifact (a clean-looking C→T decay,
    // did_3prime≈−0.12 on a real ss ancient), so it must NOT gate; the 5' C→T is the reliable end.
    constexpr float DID_FLOOR = 0.010f;
    const bool did_is_ss = (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
    const float did_signal = did_is_ss ? profile.did_5prime
                                       : std::min(profile.did_5prime, profile.did_3prime);
    const bool  did_confirms = profile.did_valid && did_signal >= DID_FLOOR;

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
    // Constant-free denominator (PART A): rescale the fitted amplitude A by the library's OWN
    // fished terminal amplitude (within-read co-occurrence d_max), NOT a hardcoded cohort mean
    // (~0.39). When co-occurrence is non-identifiable (the dilute LOW_ABUNDANCE/DiD
    // regime, where there is no data-derived amplitude) denom<0 falls through to the existing
    // abstention paths below (the `denom > 0.0` guards and the `denom <= 0.0` early return) — the
    // honest behavior is to abstain, never to divide by a magic constant. The DENOM_SD widening
    // term below is now conservative-only; the reviewable follow-up is to source the denominator CI
    // from d_max_cooccurrence.ci_lo/hi instead. (combine() at :212 is the byte-identical per-read
    // mirror; it needs d_max_cooccurrence plumbed into DamageEvidence to match — staged, see report.)
    const bool cooc_dmax_ident =
        (profile.d_max_cooccurrence.state == DamageConfidence::DETECTED ||
         profile.d_max_cooccurrence.state == DamageConfidence::TRACE) &&
        profile.d_max_cooccurrence.point > 0.0;
    const double denom = cooc_dmax_ident
                             ? static_cast<double>(profile.d_max_cooccurrence.point)
                             : -1.0;

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
        // Require CONFIRMED ds: != SINGLE_STRANDED also admitted UNKNOWN (unresolved library type), letting
        // an UNKNOWN-but-actually-ss library run ds recovery over the ds-built 3' G->A accumulator -> garbage.
        const bool ds_lib = profile.library_type == SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
        if (denom > 0.0 && ds_lib) {
            const bool a5 = profile.artifact_overcall_5p, a3 = profile.artifact_overcall_3p;
            if (a5 != a3) {  // exactly one end dead
                const bool live_is_3p = a5;  // 5' dead ⇒ 3' is the live end (FLB57md), and vice versa
                const EndDecayLRT le = live_end_decay_lrt(
                    live_is_3p ? profile.pi_pos_3prime_ds : profile.pi_pos_5prime);
                const double A = live_is_3p ? profile.d_max_3prime : profile.d_max_5prime;
                // NOTE: a hardcoded amplitude ceiling is WRONG for deep time — at ~4 Myr deamination can
                // saturate, so a genuine d_max of 0.5-0.7 is plausible; a young-aDNA cohort mean (~0.39)
                // is not a valid cap. Discrimination of real-saturated vs artifact must be by SHAPE,
                // not amplitude (see assessment: per-position window shape does NOT separate them here either).
                if (A > 0.0) {
                    // A_se: prefer the bulk damaged-split SE; else the LRT's own A_se; else 20% of A.
                    const double A_se = (profile.bulk_damage.d_max_damaged_valid &&
                                         profile.bulk_damage.d_max_se > 0.0)
                                            ? profile.bulk_damage.d_max_se
                                            : (le.A_se > 0.0 ? le.A_se : 0.2 * A);
                    // Propagate a denominator-uncertainty term too: the co-occurrence amplitude is itself
                    // estimated, so widen the pi band by a conservative fixed denom spread (~0.035 sd).
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
                    } else if (did_confirms) {
                        // The amplitude test put this deep faint ancient below floor, but the
                        // composition-immune DiD positively confirms one-sided deamination (FLB10m/FLB57md:
                        // did +0.05/+0.035 with sub-floor live-end amplitude at low depth). DiD is the
                        // stronger, artifact-robust evidence here — promote to LOW_ABUNDANCE.
                        profile.pi.point = clip01(pt);
                        profile.pi.lo    = 0.0;
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
        // DiD RECOVERY (replaces the retired Channel-B amplitude admission): channel-5 strand-symmetry
        // did not authenticate, but the composition-immune codon DiD positively confirms one-sided
        // deamination. This is the faint/low-endogenous regime (ROCS marine, deep FLB). Safe where the
        // old stop-codon path was NOT: the DiD's purine-mirror control nulls the terminal composition
        // skew (coastal water → did_paired<0 → did_confirms=false), so this cannot re-admit the modern
        // composition artifacts that the amplitude-based C1 falsely authenticated. pi value from the faint
        // Channel-A d_max that exists but did not clear the symmetric gate.
        if (denom > 0.0 && did_confirms) {
            const double A = std::max(profile.d_max_5prime, profile.d_max_3prime);
            profile.pi.point = clip01(A / denom);
            profile.pi.lo    = 0.0;
            profile.pi.hi    = clip01((A + 1.96 * 0.2 * A) / denom);
            profile.pi.state = DamageConfidence::LOW_ABUNDANCE;
            return;
        }
        try_ancient_cpg();
        return;
    }

    if (denom <= 0.0) { try_ancient_cpg(); return; }

    // dart A∧B joint table: channel-5 authenticated on SYMMETRIC terminal AMPLITUDE, but the
    // composition-immune codon DiD must CONFIRM real one-sided deamination. did flat (≈0 — e.g.
    // agricultural-soil/gut metagenomes with N-inflated ΔLLR but no net pyrimidine-over-purine excess)
    // OR negative (coastal water, did≈−0.04) → "A fires, B flat/negative → compositional artifact",
    // d_max=0 → ABSTAIN. Only bypass when the DiD is underpowered (did_valid=false, low convertible
    // coverage), where the table's B=n/a falls back to Channel A. Real damage has did ≥ +DID_FLOOR
    // (FLB ancients +0.05, ROCS +0.018, estuary relic +0.01–0.03), so this keeps them.
    if (profile.did_valid && !did_confirms) { profile.pi = DamageEstimate{}; try_ancient_cpg(); return; }

    // pi VALUE from the Briggs per-end d_max (position-peak amplitude, COMMENSURATE with the co-occurrence denom).
    // channel-5 above is the authenticity GATE ONLY — its context-scale amp (~0.03) is NOT used for the pi
    // value: amp/(cohort constant) mixed a context-axis numerator with a position-peak denominator and put pi
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
    // RANKING amplitude, data-derived (no hardcoded cohort mean): prefer the library's
    // own constant-free co-occurrence d_max, else its fitted per-end d_max. This score only RANKS reads
    // (monotone: more terminal deamination ⇒ higher), so the amplitude sets the scale, not a calibrated
    // rate; a small positive floor keeps the LLR monotone when no amplitude was fitted.
    // ceiling: floor 0.1 is a ranking scale only (never a reported damage rate); upgrade: drop when a
    // per-read joint kernel replaces this single-amplitude score.
    double A = std::max({static_cast<double>(profile.d_max_5prime),
                         static_cast<double>(profile.d_max_3prime),
                         profile.d_max_cooccurrence.point});
    if (!(A > 0.0)) A = 0.1;

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
