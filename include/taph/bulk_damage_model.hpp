#pragma once
//
// Bulk reference-free damage law (Phase 1) — count-level binomial GLM.  Design B.
//
//   logit η_{c,e,l,p} = β̂_{l,c} + a_{c,e,p} + s_{l,e}·[c==0 && p==0]   PER-CHANNEL artifact, a_{c,e,14}=0
//   μ_ctl  = sigmoid(η_ctl)
//   μ_dmg  = sigmoid(η_dmg)                        for p==0 (p0 excluded from decay)
//   μ_dmg  = h + (1−h)·δ_l·exp(−λ·p)              for p≥1,  h=sigmoid(η_dmg)
//   y_{c,e,l,p} ~ Binomial(n, μ)
//   r(p) = exp(−λ·p) for p≥1, 0 at p=0  (both ds and ss — no Briggs ss kernel)
//
// β̂_{l,c} = logit(interior C/(C+T)) is FIXED from init() (read-middle interior p∈[L/3,2L/3]);
// step_beta is DISABLED so the terminus models only the excess over read-middle.
//
// δ_l ISOTONIC: generalised isotonic regression (PAVA) non-increasing in median_len; the live bin
// with the LONGEST median length (live.back()) is PINNED to δ=0 (the δ_{Lmax}=0 lower-envelope
// anchor that closes the a_dmg/δ shift-symmetry without projections). Each block's δ is the 1-D
// nonneg binomial MLE over that block's (both-end pooled) damage cells at p≥1.
//
// a_{c,e,p} is a FREE per-(channel,end,position) length-invariant artifact (no orthogonal projection).
// The p0 spike s[l][e] handles p0-only artifacts on the damage channel.
//
// Channels resolved UPSTREAM (ds: damage=TC@5′,AG@3′ / control=AG@5′,TC@3′;
// ss: damage=TC@5′,TC@3′ / control=AG@both).
//
// Solver: coordinate ascent  a → s → δ(isotonic+anchor) → λ(golden),
//   monotone-ascent guarded (β̂ fixed).
//   - a : Fisher scoring (IRLS), step-halved, free per channel (NO projection).
//   - δ : PAVA isotonic non-increasing + δ_{Lmax}=0 anchor, 1-D block MLE (both ends pooled).
//   - λ : golden-section line search [0.05, 0.60].
//
// Threshold-free outputs: per (end,length) δ̂, K_eff, R_damage, Δℓ0/S0, and a FULL profile-likelihood
// curve over δ_l (nuisances re-optimized per grid point, scanned outward from MLE).
//
#include <array>
#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <limits>

namespace taph {

// ----------------------------------------------------------------------------- suff-stats

// Eligibility-stratified per-read damage co-occurrence moments (mixture P2). Conditioning on the
// per-read damage-channel SITE counts (S5,S3) — which are damage-INVARIANT (C→T preserves T+C) —
// removes the GC/composition confound, so Cov(k5,k3 | S5,S3) is pure ancient-mixture covariance and
// the ancient fraction π separates from d_max. JP=5 terminal positions per end; k = #damaged sites.
struct JStrat {
    static constexpr int JP = 5;
    std::array<std::array<uint64_t, JP + 1>, JP + 1> n{};     // [S5][S3] reads
    std::array<std::array<uint64_t, JP + 1>, JP + 1> sk5{};   // Σ k5
    std::array<std::array<uint64_t, JP + 1>, JP + 1> sk3{};   // Σ k3
    std::array<std::array<uint64_t, JP + 1>, JP + 1> sk53{};  // Σ k5·k3
    void add(const JStrat& o) {
        for (int a = 0; a <= JP; ++a)
            for (int b = 0; b <= JP; ++b) {
                n[a][b] += o.n[a][b];   sk5[a][b]  += o.sk5[a][b];
                sk3[a][b] += o.sk3[a][b]; sk53[a][b] += o.sk53[a][b];
            }
    }
};

struct BulkDamageSuffStats {
    static constexpr int N_POS = 15;   // terminal positions 0..14 per end
    static constexpr int N_END = 2;    // 0 = 5′, 1 = 3′
    static constexpr int N_CH  = 2;    // 0 = damage, 1 = control

    struct Cell {
        std::array<uint64_t, N_POS> k{};   // successes (damage-substitution count)
        std::array<uint64_t, N_POS> n{};   // trials
    };
    // bin[l][c][e] terminal counts; interior[c][l] for β init; median_len[l] sets isotonic order.
    std::vector<std::array<std::array<Cell, N_END>, N_CH>> bin;
    std::vector<std::array<uint64_t, N_CH>> k_interior;   // [l][c]
    std::vector<std::array<uint64_t, N_CH>> n_interior;   // [l][c]
    std::vector<double> median_len;                       // [l], ascending

    // Per-read eligibility-stratified damage co-occurrence (mixture P2 suff-stat); jstrat[l] holds the
    // (S5,S3)-keyed moments {n, Σk5, Σk3, Σk5·k3} for bin l (ds/ss channel resolved at aggregation).
    std::vector<JStrat> jstrat;                           // [l]

    bool ss = false;            // diagnostic only (channels already resolved upstream)
    bool skip_3p_pos0 = true;   // ss ligation artifact: drop 3′ p=0 from the damage score

    int L() const { return static_cast<int>(bin.size()); }

    // A length bin contributes if its damage channel carries usable terminal coverage.
    bool bin_valid(int l) const {
        uint64_t min_dmg = std::numeric_limits<uint64_t>::max();
        for (int e = 0; e < N_END; ++e)
            for (int p = 0; p < N_POS; ++p)
                min_dmg = std::min(min_dmg, bin[l][0][e].n[p]);
        return min_dmg >= 50 && n_interior[l][0] >= 500;
    }
};

// ----------------------------------------------------------------------------- per-bin result

struct BulkDamagePerBin {
    int    length_lo = 0, length_hi = 0;
    double median_len = 0.0;
    int64_t n_reads = 0;

    double delta = 0.0;            // δ̂_l (isotonic non-increasing MLE; δ_{Lmax}=0 anchor)
    bool   identified = false;     // δ̂_l estimated independently (singleton PAVA block)
    bool   borrowed = false;       // PAVA pooled with at least one neighbour (block_size > 1)

    // threshold-free diagnostics, per end (0=5′,1=3′)
    std::array<double, 2> k_eff   = {0.0, 0.0};   // effective informative positions (IPR)
    std::array<double, 2> r_damage = {0.0, 0.0};  // (D_ctrl−D_full)/D_ctrl on the damage channel
    double delta_ell0 = 0.0;       // ℓ̂_l − ℓ(δ_l=0), nuisances re-profiled (full profile)
    double s0 = 1.0;               // exp(−delta_ell0) = L(δ_l=0)/L(δ̂_l)

    // FULL profile-likelihood curve over δ_l (nuisances re-optimized per grid point)
    std::vector<double> profile_delta;    // grid
    std::vector<double> profile_loglik;   // ℓ(δ_l=grid) − ℓ̂  (≤ 0, peak 0 at δ̂_l)

    // mixture P2 — eligibility-conditioned ancient/modern split (per-read two-end co-occurrence)
    uint64_t joint_n = 0;                       // reads contributing to the co-occurrence moments
    std::array<double, 2> joint_mean = {0.0, 0.0};  // marginal E[k5], E[k3] (diagnostic)
    double joint_cov = 0.0;                     // marginal Cov(k5,k3) (diagnostic; GC-confounded)
    double pi_ancient = -1.0;                   // ancient read fraction π_l = δ_l/d_max (−1 = undetermined)
    double pi_lo = -1.0, pi_hi = -1.0;          // π_l 95% interval (threshold-free; spans → undetermined)
};

// ----------------------------------------------------------------------------- whole-fit result

struct BulkDamageResult {
    bool   converged = false;
    bool   valid = false;
    int    n_sweeps = 0;
    double log_lik = 0.0;

    double lambda = 0.0;
    std::array<std::array<double, 2>, 2> kappa = {};   // κ[c][e]; control row ≡ 1
    std::array<std::array<double, 15>, 2> artifact = {};  // a[e][p], a[e][14]=0

    std::vector<BulkDamagePerBin> bins;

    // Length-coupling weight (post-fit, threshold-free). Real terminal deamination occupies a fixed
    // terminal zone, so its damaged-base mass per read m(L)=δ_L·L is flat-or-falling in read length L
    // (slope ≤ 0); a per-base/whole-molecule artifact (FLB57m) adds a constant rate, so m(L)∝L
    // (slope ≫ 0). slope_m = OLS slope of m(L) vs L over live bins; the centre slope_m=0 is the
    // mechanistic terminal-vs-pervasive boundary. w_length = σ(−slope_m / (c·SE)) ∈[0,1] with SE from
    // the per-bin profile-likelihood curvature → used as δ_auth = δ·w_length. (Undetermined → w≈0.5
    // when δ_L is too weakly determined to fix the slope sign.)
    double slope_m  = 0.0;
    double w_length = 1.0;

    // Length-coupling DIAGNOSTIC (post-fit, NOT used in the fit). Sign of the slope of the recovered δ_l
    // against median read length over live bins: −1 ⇒ δ decreasing with length (classic length-coupled
    // terminal deamination, e.g. FLB08m), +1 ⇒ increasing, 0 ⇒ flat/undetermined.
    double length_coupling_slope = 0.0;   // OLS slope of δ_l vs median_len over live bins
    int    length_coupling = 0;           // sign(slope): −1 decreasing, +1 increasing, 0 flat

    // mixture P2 — shared per-ancient-read deamination intensity d_max (reference-free analog of
    // metaDMG A_b), estimated from the eligibility-conditioned two-end co-occurrence with δ_l pinned
    // (π_l = δ_l/d_max). d_max = −1 ⇒ undetermined (no length bin carried a usable co-occurrence signal,
    // e.g. δ≈0 everywhere or the artifact gate fired). Threshold-free: d_max_se carries the data's own
    // uncertainty; the per-bin π interval spans → undetermined where the signal can't fix the split.
    double d_max_ancient = -1.0;
    double d_max_se      = 0.0;
    double d_max_raw     = -1.0;   // pre-gate pooled estimate; >1 (non-physical) ⇒ railed ⇒ unidentified
};

// ----------------------------------------------------------------------------- model / solver

class BulkDamageModel {
public:
    static BulkDamageResult fit(const BulkDamageSuffStats& s);

    // grid / bounds
    static constexpr int    N_POS       = BulkDamageSuffStats::N_POS;
    static constexpr int    N_END       = BulkDamageSuffStats::N_END;
    static constexpr double LAMBDA_MIN  = 0.05;
    static constexpr double LAMBDA_MAX  = 0.60;
    static constexpr double DELTA_MIN   = 0.0;
    static constexpr double DELTA_MAX   = 0.60;
    static constexpr double MU_EPS      = 1e-6;
    static constexpr int    MAX_SWEEPS  = 500;
    static constexpr double TOL_PARAM   = 1e-7;   // scale-free: max |Δ param| between sweeps
    // LL-plateau stop: on a weakly-identified (near-zero-δ) ridge the parameters crawl for
    // hundreds of sweeps while the likelihood is already flat. The LL is monotone (every step
    // is ascent-guarded), so we stop once the per-sweep gain is negligible — this terminates on
    // the ridge where the parameter-change test never would.
    static constexpr double TOL_LL_ABS  = 1e-7;   // tiny floor for |ℓ|→0
    static constexpr double TOL_LL_REL  = 1e-8;   // gain threshold = TOL_LL_ABS + TOL_LL_REL·|ℓ|
                                                  // relative ⇒ scale-invariant across N (2M…1e8 reads).
                                                  // δ is the exact block-MLE given nuisances and settles
                                                  // long before this; the bar only bounds nuisance crawl.
    static constexpr int    PROFILE_PTS = 41;   // δ grid resolution for the full profile
    static constexpr int    PROFILE_MAX_SWEEPS = 80;  // profile points are warm-started → few sweeps
    // w_length softness: the centre slope_m=0 is the mechanistic terminal-vs-pervasive boundary; this
    // only sets the transition width (in slope-SE units). High-coverage libraries fix the slope sign at
    // |z|≫1 regardless of c; c=1 (one SE) shapes only the undetermined band of thin libraries.
    static constexpr double W_LENGTH_C = 1.0;

private:
    // free parameters during the fit
    struct Params {
        std::vector<std::array<double, 2>> beta;          // β[l][c]
        std::vector<double>                delta;         // δ[l]  (isotonic non-increasing, δ_{Lmax}=0)
        double                             lambda = 0.3;
        // FREE per-channel, length-INVARIANT terminal artifact (design B — no projection):
        //   logit η_{c,e,l,p} = β̂[l][c] + a[c][e][p] + s[l][e]·[c==0 && p==0]
        // a[c][e][p] is free per (channel,end,position) over live bins; a[c][e][14]=0 anchored.
        std::array<std::array<std::array<double, 15>, N_END>, 2> a = {}; // a[c][e][p], a[c][e][14]=0
        std::vector<std::array<double, N_END>>             s_spike;      // s[l][e]: p0 spike on η_dmg per (e,l)
    };

    static constexpr double SPIKE_RIDGE_N   = 2000.0; // weak shrinkage on the p0 spike s[l][e]

    static double sigmoid(double x) {
        if (x >= 0) { double z = std::exp(-x); return 1.0 / (1.0 + z); }
        double z = std::exp(x); return z / (1.0 + z);
    }
    static double logit(double p) {
        p = std::clamp(p, 1e-9, 1.0 - 1e-9);
        return std::log(p / (1.0 - p));
    }
    static double clamp_mu(double m) { return std::clamp(m, MU_EPS, 1.0 - MU_EPS); }

    // Damage decay kernel r(p;λ).  r(p) = exp(−λ·p) for p≥1, 0 at p=0.
    // Identical for ds and ss (no Briggs ss-overhang term in design B).
    // p=0 is always excluded from the decay: the p0 spike s[l][e] handles terminus-only artifacts there.
    static double r_kernel(const Params& P, bool /*ss*/, int p) {
        return (p >= 1) ? std::exp(-P.lambda * p) : 0.0;
    }

    // η_{c,e,l,p} = σ(β̂[l][c] + a[c][e][p] + s[l][e]·[c==0 && p==0]) — FREE per-channel artifact, p0 spike dmg-only.
    static double eta(const Params& P, int c, int e, int l, int p) {
        double spike = (c == 0 && p == 0) ? P.s_spike[l][e] : 0.0;
        return sigmoid(P.beta[l][c] + P.a[c][e][p] + spike);
    }
    // μ for a cell; damage channel (c==0) carries δ_l·r(p) (r=0 at p0 ⇒ p0 excluded from decay), control plain.
    static double mu(const Params& P, bool ss, int c, int e, int l, int p) {
        double h = eta(P, c, e, l, p);
        if (c == 1) return clamp_mu(h);
        double r = r_kernel(P, ss, p);
        if (r == 0.0) return clamp_mu(h);
        return clamp_mu(h + (1.0 - h) * P.delta[l] * r);
    }

    static bool cell_live(const BulkDamageSuffStats& s, int c, int e, int p) {
        // drop 3′ position 0 from the damage score when requested (ss ligation artifact)
        return !(s.skip_3p_pos0 && c == 0 && e == 1 && p == 0);
    }

    static double log_lik(const BulkDamageSuffStats& s, const Params& P);
    static double penalized_obj(const BulkDamageSuffStats& s, const Params& P);
    static void   init(const BulkDamageSuffStats& s, Params& P, const std::vector<int>& live);
    static void   step_beta(const BulkDamageSuffStats& s, Params& P, const std::vector<int>& live);
    static void   step_spike(const BulkDamageSuffStats& s, Params& P, const std::vector<int>& live);
    static void   step_artifact(const BulkDamageSuffStats& s, Params& P, const std::vector<int>& live);
    // step_delta_isotonic: PAVA isotonic non-increasing in median_len; live.back() pinned to δ=0 (Lmax anchor).
    // Each PAVA block's δ = 1-D nonneg binomial MLE (delta_block_mle) over its pooled damage cells (p≥1).
    static void   step_delta_isotonic(const BulkDamageSuffStats& s, Params& P, const std::vector<int>& live,
                                      std::vector<int>* block_id, int pinned_idx = -1, double pinned_val = 0.0);
    static void   step_lambda(const BulkDamageSuffStats& s, Params& P, const std::vector<int>& live);
    static double max_param_change(const Params& A, const Params& B, const std::vector<int>& live);
    static int    run_sweeps(const BulkDamageSuffStats& s, Params& P, const std::vector<int>& live,
                             int pinned_idx, double pinned_val, std::vector<int>* block_id,
                             int max_sweeps = MAX_SWEEPS);

    // 1-D bounded MLE of a single δ over an arbitrary set of damage cells (nuisances fixed).
    // The binomial score is monotone-decreasing in δ (concave LL) → safeguarded bisection.
    template <class CellSet>
    static double delta_block_mle(const BulkDamageSuffStats& s, const Params& P,
                                  const CellSet& cells);
};

// ------------------------------------------------------------------ log-likelihood

inline double BulkDamageModel::log_lik(const BulkDamageSuffStats& s, const Params& P) {
    double ll = 0.0;
    const int L = s.L();
    for (int l = 0; l < L; ++l) {
        for (int c = 0; c < BulkDamageSuffStats::N_CH; ++c) {
            for (int e = 0; e < N_END; ++e) {
                const auto& cell = s.bin[l][c][e];
                for (int p = 0; p < N_POS; ++p) {
                    if (!cell_live(s, c, e, p)) continue;
                    uint64_t n = cell.n[p];
                    if (n == 0) continue;
                    double m = mu(P, s.ss, c, e, l, p);
                    double k = static_cast<double>(cell.k[p]);
                    ll += k * std::log(m) + (static_cast<double>(n) - k) * std::log(1.0 - m);
                }
            }
        }
    }
    return ll;
}

// Penalized objective = log-likelihood − weak shrink on the p0 spike. The spike shrink keeps s[l][e]
// from absorbing genuine signal; used as the ascent guard for the spike step.
inline double BulkDamageModel::penalized_obj(const BulkDamageSuffStats& s, const Params& P) {
    double obj = log_lik(s, P);
    const int L = s.L();
    for (int l = 0; l < L && l < static_cast<int>(P.s_spike.size()); ++l)
        for (int e = 0; e < N_END; ++e) obj -= SPIKE_RIDGE_N * P.s_spike[l][e] * P.s_spike[l][e];
    return obj;
}

// ------------------------------------------------------------------ init

inline void BulkDamageModel::init(const BulkDamageSuffStats& s, Params& P,
                                  const std::vector<int>& live) {
    const int L = s.L();
    P.beta.assign(L, {0.0, 0.0});
    P.delta.assign(L, 0.0);
    P.lambda = 0.3;
    P.s_spike.assign(L, {0.0, 0.0});                    // s[l][e]: init 0
    for (int c = 0; c < BulkDamageSuffStats::N_CH; ++c)
        for (int e = 0; e < N_END; ++e) P.a[c][e].fill(0.0);   // a[c][e][p]: init 0

    // β_{c,l} = logit(interior rate); interior carries neither artifact nor decay so η≈interior.
    for (int l : live) {
        for (int c = 0; c < BulkDamageSuffStats::N_CH; ++c) {
            uint64_t n = s.n_interior[l][c], k = s.k_interior[l][c];
            P.beta[l][c] = logit(n > 0 ? static_cast<double>(k) / static_cast<double>(n) : 0.5);
        }
    }
}

// ------------------------------------------------------------------ β step (Fisher scoring)

inline void BulkDamageModel::step_beta(const BulkDamageSuffStats& s, Params& P,
                                       const std::vector<int>& live) {
    for (int l : live) {
        for (int c = 0; c < BulkDamageSuffStats::N_CH; ++c) {
            double score = 0.0, info = 0.0;
            for (int e = 0; e < N_END; ++e) {
                const auto& cell = s.bin[l][c][e];
                for (int p = 0; p < N_POS; ++p) {
                    if (!cell_live(s, c, e, p)) continue;
                    uint64_t n = cell.n[p];
                    if (n == 0) continue;
                    double h = eta(P, c, e, l, p);
                    double m = mu(P, s.ss, c, e, l, p);
                    double r = (c == 0) ? r_kernel(P, s.ss, p) : 0.0;
                    double dmu = (1.0 - (c == 0 ? P.delta[l] * r : 0.0)) * h * (1.0 - h); // ∂μ/∂β
                    score += (static_cast<double>(cell.k[p]) - static_cast<double>(n) * m) /
                             (m * (1.0 - m)) * dmu;     // (k − nμ)/(μ(1−μ)) · ∂μ/∂β
                    info  += static_cast<double>(n) / (m * (1.0 - m)) * dmu * dmu;
                }
            }
            if (info > 0) {
                double base = log_lik(s, P);
                double full = score / info, t = 1.0;
                double saved = P.beta[l][c];
                for (int b = 0; b < 12; ++b) {
                    P.beta[l][c] = saved + t * full;
                    if (log_lik(s, P) >= base) break;
                    t *= 0.5;
                    if (b == 11) P.beta[l][c] = saved;
                }
            }
        }
    }
}

// ----------------------------------------------- s[l][e]: p0 spike on η_dmg (per end,length), weakly shrunk
//
// A strong position-0-only damage-channel artifact (FLB10m/Thomaz ligation/end-repair) would otherwise rail
// λ→0 (the decay tail bends to chase p0) or leak into δ. s[l][e] is a sign-free spike on η_dmg at p==0 only,
// per (e,l), weakly shrunk (SPIKE_RIDGE_N). Because the damage decay now starts at p≥1, s absorbs the p0
// artifact while λ and δ are fit from the p≥1 tail. Fisher-scored on the penalized objective. Init 0.
inline void BulkDamageModel::step_spike(const BulkDamageSuffStats& s, Params& P,
                                        const std::vector<int>& live) {
    const int p = 0;
    for (int l : live) {
        for (int e = 0; e < N_END; ++e) {
            if (!cell_live(s, 0, e, p)) continue;        // ss 3′ p0 dropped → no spike there
            const auto& cell = s.bin[l][0][e];
            uint64_t n = cell.n[p];
            if (n == 0) continue;
            double h = eta(P, 0, e, l, p);
            double m = mu(P, s.ss, 0, e, l, p);          // ds p0: μ=η; ss-Briggs p0 carries δ·r(0)
            double r0 = r_kernel(P, s.ss, p);            // ds: 0; ss-Briggs: (1−w)+w (incl p0)
            double dmu = (1.0 - P.delta[l] * r0) * h * (1.0 - h);   // ∂μ/∂s (Briggs factor at p0)
            double score = (static_cast<double>(cell.k[p]) - static_cast<double>(n) * m) /
                           (m * (1.0 - m)) * dmu;
            double info  = static_cast<double>(n) / (m * (1.0 - m)) * dmu * dmu;
            score -= 2.0 * SPIKE_RIDGE_N * P.s_spike[l][e];   // weak shrink gradient
            info  += 2.0 * SPIKE_RIDGE_N;                     // weak shrink curvature
            if (info > 0) {
                double base = log_lik(s, P) - SPIKE_RIDGE_N * P.s_spike[l][e] * P.s_spike[l][e];
                double full = score / info, t = 1.0, saved = P.s_spike[l][e];
                for (int b = 0; b < 12; ++b) {
                    P.s_spike[l][e] = saved + t * full;
                    double cand = log_lik(s, P) - SPIKE_RIDGE_N * P.s_spike[l][e] * P.s_spike[l][e];
                    if (cand >= base) break;
                    t *= 0.5;
                    if (b == 11) P.s_spike[l][e] = saved;
                }
            }
        }
    }
}

// ----------------------------------------------- a[c][e][p]: FREE per-channel length-invariant artifact
//
// a[c][e][p] is a free per-(channel,end,position) log-odds shift on η, SHARED across length bins
// (length-INVARIANT), with a[c][e][14]=0 anchored. Fisher-scored independently per channel over all live
// bins. Design B: NO orthogonal projection — a[c][e][p] is fit freely for both channels. Identification
// of the a_dmg/δ shift symmetry is carried by the δ_{Lmax}=0 anchor (step_delta_isotonic).
inline void BulkDamageModel::step_artifact(const BulkDamageSuffStats& s, Params& P,
                                           const std::vector<int>& live) {
    for (int c = 0; c < BulkDamageSuffStats::N_CH; ++c) {
        for (int e = 0; e < N_END; ++e) {
            for (int p = 0; p < N_POS - 1; ++p) {   // p = 0..13 free; p=14 anchored to 0
                if (!cell_live(s, c, e, p)) continue;
                double score = 0.0, info = 0.0;
                for (int l : live) {
                    const auto& cell = s.bin[l][c][e];
                    uint64_t n = cell.n[p];
                    if (n == 0) continue;
                    double h = eta(P, c, e, l, p);
                    double m = mu(P, s.ss, c, e, l, p);
                    double r = (c == 0) ? r_kernel(P, s.ss, p) : 0.0;
                    double dmu = (1.0 - (c == 0 ? P.delta[l] * r : 0.0)) * h * (1.0 - h);   // ∂μ/∂a
                    score += (static_cast<double>(cell.k[p]) - static_cast<double>(n) * m) /
                             (m * (1.0 - m)) * dmu;
                    info  += static_cast<double>(n) / (m * (1.0 - m)) * dmu * dmu;
                }
                if (info > 0) {
                    double base = log_lik(s, P), full = score / info, t = 1.0, saved = P.a[c][e][p];
                    for (int b = 0; b < 12; ++b) {
                        P.a[c][e][p] = saved + t * full;
                        if (log_lik(s, P) >= base) break;
                        t *= 0.5;
                        if (b == 11) P.a[c][e][p] = saved;
                    }
                }
            }
        }
    }
    // No projection (design B): a[c][e][p] left free; δ_{Lmax}=0 closes the shift symmetry.
}


// ------------------------------------------------------------------ δ block MLE (1-D safeguarded)

template <class CellSet>
inline double BulkDamageModel::delta_block_mle(const BulkDamageSuffStats& s, const Params& P,
                                               const CellSet& cells) {
    // 1-D nonneg binomial MLE for δ over the supplied cell set (nuisances fixed).
    // Identification of the a_dmg/δ shift symmetry is via the δ_{Lmax}=0 anchor in step_delta_isotonic.
    // score(δ) = Σ_{p≥1} (k − nμ)/(μ(1−μ)) · (1−η)·r(p), decreasing in δ. Root in [0, DELTA_MAX].
    auto score = [&](double d) {
        double g = 0.0;
        for (const auto& cl : cells) {
            int e = cl.e, l = cl.l;
            const auto& cell = s.bin[l][0][e];
            for (int p = 0; p < N_POS; ++p) {            // ds: r=0 at p0 (skipped); ss-Briggs: incl p0
                if (!cell_live(s, 0, e, p)) continue;
                double r = r_kernel(P, s.ss, p);
                if (r == 0.0) continue;
                uint64_t n = cell.n[p];
                if (n == 0) continue;
                double h = eta(P, 0, e, l, p);
                double m = clamp_mu(h + (1.0 - h) * d * r);
                g += (static_cast<double>(cell.k[p]) - static_cast<double>(n) * m) /
                     (m * (1.0 - m)) * (1.0 - h) * r;
            }
        }
        return g;
    };
    double s0 = score(DELTA_MIN);
    if (s0 <= 0.0) return DELTA_MIN;                 // boundary at 0
    double s1 = score(DELTA_MAX);
    if (s1 >= 0.0) return DELTA_MAX;                 // boundary at δmax
    double lo = DELTA_MIN, hi = DELTA_MAX;
    for (int it = 0; it < 60; ++it) {
        double mid = 0.5 * (lo + hi);
        if (score(mid) > 0.0) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

// ------------------------------------------------------------------ δ step: isotonic non-increasing + Lmax=0 anchor
//
// Design B: generalised isotonic regression (PAVA, pool-adjacent violators) non-increasing in median_len
// (bins in live[] are already sorted ascending by median_len). The live bin with the LONGEST median length
// (live.back()) is FORCED to δ=0 — the δ_{Lmax}=0 lower-envelope anchor that closes the a_dmg/δ shift
// symmetry without projections. Each PAVA block's δ is the 1-D nonneg binomial MLE (delta_block_mle) over
// all (both-end, p≥1) damage cells in the block. block_id[l] identifies which PAVA block bin l belongs to;
// block_size>1 ⇒ borrowed. pinned_idx≥0 overrides: that bin is held at pinned_val (profile-LL pin).
inline void BulkDamageModel::step_delta_isotonic(const BulkDamageSuffStats& s, Params& P,
                                                  const std::vector<int>& live,
                                                  std::vector<int>* block_id,
                                                  int pinned_idx, double pinned_val) {
    struct LE { int e, l; };
    const int L = s.L();
    if (block_id) block_id->assign(L, -1);
    const int N = static_cast<int>(live.size());
    if (N == 0) return;

    // Index of the longest live bin (Lmax anchor): live is sorted ascending, so it is live.back().
    const int lmax_idx = live.back();

    // PAVA pass: blocks track a list of live-bin indices; δ for a block = block-MLE over all its cells.
    // Invariant after each push: the sequence of block δ values is non-increasing.
    struct Block {
        std::vector<int> bins;   // live bin indices in this block
        double           delta;  // current block δ
    };
    std::vector<Block> blocks;
    blocks.reserve(N);

    auto block_delta = [&](const std::vector<int>& bins) -> double {
        // Forced anchor: if the block contains lmax_idx, return 0.
        for (int bl : bins) if (bl == lmax_idx) return 0.0;
        // Profile pin: if the block contains pinned_idx, return pinned_val.
        for (int bl : bins) if (bl == pinned_idx) return pinned_val;
        // Collect cells for all bins in the block (both ends, pooled).
        std::vector<LE> cells;
        cells.reserve(bins.size() * 2);
        for (int bl : bins) { cells.push_back({0, bl}); cells.push_back({1, bl}); }
        return std::max(DELTA_MIN, delta_block_mle(s, P, cells));
    };

    for (int i = 0; i < N; ++i) {
        int l = live[i];
        Block cur;
        cur.bins = {l};
        cur.delta = block_delta(cur.bins);
        // Pool with previous block while it violates non-increasing order (δ_prev < δ_this).
        while (!blocks.empty() && blocks.back().delta < cur.delta) {
            cur.bins.insert(cur.bins.begin(), blocks.back().bins.begin(), blocks.back().bins.end());
            blocks.pop_back();
            cur.delta = block_delta(cur.bins);
        }
        blocks.push_back(std::move(cur));
    }

    // Write δ values back and assign block ids.
    int bid = 0;
    for (const auto& blk : blocks) {
        for (int bl : blk.bins) {
            P.delta[bl] = blk.delta;
            if (block_id) (*block_id)[bl] = bid;
        }
        ++bid;
    }
}


// ------------------------------------------------------------------ λ step: golden section

inline void BulkDamageModel::step_lambda(const BulkDamageSuffStats& s, Params& P,
                                         const std::vector<int>& /*live*/) {
    const double phi = 0.6180339887498949;
    double a = LAMBDA_MIN, b = LAMBDA_MAX;
    double c = b - phi * (b - a), d = a + phi * (b - a);
    auto f = [&](double lam) { double sv = P.lambda; P.lambda = lam; double v = log_lik(s, P); P.lambda = sv; return v; };
    double fc = f(c), fd = f(d);
    for (int it = 0; it < 40; ++it) {
        if (fc > fd) { b = d; d = c; fd = fc; c = b - phi * (b - a); fc = f(c); }
        else        { a = c; c = d; fc = fd; d = a + phi * (b - a); fd = f(d); }
    }
    P.lambda = 0.5 * (a + b);
}

// ------------------------------------------------------------------ convergence helpers

inline double BulkDamageModel::max_param_change(const Params& A, const Params& B,
                                                const std::vector<int>& live) {
    double m = std::fabs(A.lambda - B.lambda);
    for (int c = 0; c < BulkDamageSuffStats::N_CH; ++c)
        for (int e = 0; e < N_END; ++e)
            for (int p = 0; p < N_POS; ++p)
                m = std::max(m, std::fabs(A.a[c][e][p] - B.a[c][e][p]));
    for (int l : live) {
        m = std::max(m, std::fabs(A.beta[l][0] - B.beta[l][0]));
        m = std::max(m, std::fabs(A.beta[l][1] - B.beta[l][1]));
        m = std::max(m, std::fabs(A.delta[l] - B.delta[l]));
        for (int e = 0; e < N_END; ++e)
            m = std::max(m, std::fabs(A.s_spike[l][e] - B.s_spike[l][e]));
    }
    return m;
}

// Coordinate ascent to convergence (scale-free param-change stop). pinned_idx≥0 holds δ_{pinned}=pinned_val.
inline int BulkDamageModel::run_sweeps(const BulkDamageSuffStats& s, Params& P,
                                       const std::vector<int>& live,
                                       int pinned_idx, double pinned_val,
                                       std::vector<int>* block_id, int max_sweeps) {
    int sweep = 0;
    double ll_prev = log_lik(s, P);
    for (; sweep < max_sweeps; ++sweep) {
        Params old = P;
        // β̂_{c,l} is FIXED at the interior sufficient statistic (logit(interior C/(C+T)), set in init()) —
        // NOT re-fit from terminal data. A pervasive per-base artifact (FLB57m) elevates the INTERIOR
        // baseline; if β were free to the terminus, that elevation would split between β (short bins) and
        // δ, and the δ_{Lmax} anchor would relabel it as short-read excess (false damage). Anchoring β to
        // the interior makes the terminus model only the EXCESS OVER INTERIOR, so a per-base artifact
        // (terminus ≤ interior) gives δ→0 while genuine terminal deamination (terminus > interior) gives
        // δ>0 — the pervasiveness discriminator, baked into the model. (room-cdf15767)
        // step_beta(s, P, live);   // intentionally disabled — β̂ anchored to interior
        step_artifact(s, P, live);                    // a[c][e][p]: free per channel, no projection
        step_spike(s, P, live);                       // s[l][e]: p0 spike on η_dmg
        step_delta_isotonic(s, P, live, block_id, pinned_idx, pinned_val);   // PAVA + δ_{Lmax}=0
        step_lambda(s, P, live);
        double ll = log_lik(s, P);
        // Stop on EITHER a tiny parameter move OR a flat likelihood. The LL is monotone
        // (every step is ascent-guarded) so ll − ll_prev ≥ 0; on the near-zero-δ ridge the
        // parameters keep sliding while the LL has already plateaued — the LL test stops there.
        bool ll_flat   = (ll - ll_prev) <= TOL_LL_ABS + TOL_LL_REL * std::fabs(ll);
        bool par_small = max_param_change(old, P, live) < TOL_PARAM;
        if (par_small || ll_flat) { ++sweep; break; }
        ll_prev = ll;
    }
    return sweep;
}

// ------------------------------------------------------------------ driver

inline BulkDamageResult BulkDamageModel::fit(const BulkDamageSuffStats& s) {
    BulkDamageResult R;
    const int L = s.L();
    if (L == 0) return R;

    std::vector<int> live;
    for (int l = 0; l < L; ++l) if (s.bin_valid(l)) live.push_back(l);
    if (live.empty()) return R;
    R.valid = true;

    Params P;
    init(s, P, live);

    std::vector<int> block_id;
    int sweeps = run_sweeps(s, P, live, -1, 0.0, &block_id);
    R.converged = sweeps < MAX_SWEEPS;
    R.n_sweeps  = sweeps;
    R.log_lik   = log_lik(s, P);
    R.lambda   = P.lambda;
    // R.artifact = the damage-channel free terminal artifact a[0][e][p] (log-odds shift over β̂[l][0]) for
    // JSON. a[0][e][14]≡0 anchored. The p0 spike s[l][e] is length-dependent; fold its live-bin mean into
    // the p0 entry so the reported curve reflects the effective p0 shift used by the fit.
    for (int e = 0; e < N_END; ++e) {
        for (int p = 0; p < N_POS; ++p) {
            double spike_mean = 0.0;
            if (p == 0) {
                int cnt = 0;
                for (int l : live) { spike_mean += P.s_spike[l][e]; ++cnt; }
                spike_mean = cnt > 0 ? spike_mean / cnt : 0.0;
            }
            R.artifact[e][p] = P.a[0][e][p] + spike_mean;
        }
    }
    R.kappa[0] = {1.0, 1.0};             // κ retired (per-channel a); identity rows for JSON compat
    R.kappa[1] = {1.0, 1.0};

    // per-bin outputs + threshold-free diagnostics
    R.bins.resize(L);
    for (int l = 0; l < L; ++l) {
        auto& B = R.bins[l];
        B.median_len = l < static_cast<int>(s.median_len.size()) ? s.median_len[l] : 0.0;
        B.delta = P.delta[l];
        bool is_live = std::find(live.begin(), live.end(), l) != live.end();
        if (!is_live) continue;

        // Isotonic PAVA: block_size>1 ⇒ bin was pooled with a neighbour (borrowed);
        // singleton block ⇒ identified independently. NOT a damage verdict — read δ̂/Δℓ0/K_eff.
        int my_block = block_id[l];
        int block_size = 0;
        for (int q : live) if (block_id[q] == my_block) ++block_size;
        B.borrowed   = block_size > 1;
        B.identified = !B.borrowed;

        // K_eff and R_damage per end (damage channel only)
        for (int e = 0; e < N_END; ++e) {
            const auto& cell = s.bin[l][0][e];
            double sumI = 0.0, sumI2 = 0.0;
            double dev_full = 0.0, dev_ctrl = 0.0;
            for (int p = 0; p < N_POS; ++p) {           // r=0 at p0 (excluded from decay)
                if (!cell_live(s, 0, e, p)) continue;
                uint64_t n = cell.n[p];
                if (n == 0) continue;
                double r = r_kernel(P, s.ss, p);
                if (r == 0.0) continue;
                double h = eta(P, 0, e, l, p);
                double m  = clamp_mu(h + (1.0 - h) * P.delta[l] * r);   // full model
                double m0 = clamp_mu(h);                                // δ≡0 (artifact only)
                double Ip = static_cast<double>(n) * std::pow((1.0 - h) * r, 2.0) / (m * (1.0 - m));
                sumI += Ip; sumI2 += Ip * Ip;
                double k = static_cast<double>(cell.k[p]), nn = static_cast<double>(n);
                double obs = std::clamp(k / nn, MU_EPS, 1.0 - MU_EPS);
                auto dev = [&](double mm) {
                    return 2.0 * (k * std::log(obs / mm) + (nn - k) * std::log((1.0 - obs) / (1.0 - mm)));
                };
                dev_full += dev(m);
                dev_ctrl += dev(m0);
            }
            B.k_eff[e]    = sumI2 > 0 ? (sumI * sumI) / sumI2 : 0.0;
            B.r_damage[e] = dev_ctrl > 1e-12 ? (dev_ctrl - dev_full) / dev_ctrl : 0.0;
        }
    }

    // FULL profile-likelihood per live bin: pin δ_l on a grid, re-optimize all nuisances + other δ,
    // warm-started outward from the MLE. Δℓ0 read off the δ_l=0 grid point.
    for (int l : live) {
        auto& B = R.bins[l];
        B.profile_delta.assign(PROFILE_PTS, 0.0);
        B.profile_loglik.assign(PROFILE_PTS, 0.0);
        for (int g = 0; g < PROFILE_PTS; ++g) B.profile_delta[g] = DELTA_MAX * g / (PROFILE_PTS - 1);
        int g0 = std::clamp((int)std::lround(P.delta[l] / DELTA_MAX * (PROFILE_PTS - 1)), 0, PROFILE_PTS - 1);
        // walk outward from the MLE, each grid point warm-started from the previous converged fit;
        // δ_l is pinned and ALL other params (incl. the isotonic neighbours) are re-optimized → full profile.
        Params Q = P;
        for (int g = g0; g < PROFILE_PTS; ++g) {
            run_sweeps(s, Q, live, l, B.profile_delta[g], nullptr, PROFILE_MAX_SWEEPS);
            B.profile_loglik[g] = log_lik(s, Q) - R.log_lik;
        }
        Q = P;
        for (int g = g0 - 1; g >= 0; --g) {
            run_sweeps(s, Q, live, l, B.profile_delta[g], nullptr, PROFILE_MAX_SWEEPS);
            B.profile_loglik[g] = log_lik(s, Q) - R.log_lik;
        }
        // Anchor the curve to its own maximum so it is a proper relative profile LL with peak
        // exactly 0; also absorbs any finite-precision over-shoot from differencing two large LLs.
        double pk = *std::max_element(B.profile_loglik.begin(), B.profile_loglik.end());
        for (double& v : B.profile_loglik) v -= pk;
        B.delta_ell0 = -B.profile_loglik[0];     // Δℓ0 = ℓ̂ − ℓ(δ_l=0)
        B.s0 = std::exp(-std::max(0.0, B.delta_ell0));
    }

    // ── length-coupling weight w_length (threshold-free) ─────────────────────────────────────────
    // Real terminal deamination occupies a fixed terminal zone → damaged-base mass per read
    // m(L)=δ_L·L is flat-or-falling in read length (slope_m ≤ 0). A per-base/whole-molecule artifact
    // adds a constant rate → m(L)∝L (slope_m ≫ 0). slope_m = OLS slope of m vs L over live bins;
    // w_length = σ(−slope_m / (c·SE_slope)). Var(δ_l) is read from the curvature of each bin's profile
    // likelihood (the data's own Fisher information): a bin whose δ is weakly determined widens SE and
    // pulls w toward 0.5 (undetermined) instead of forcing a verdict.
    {
        std::vector<double> Lv, mv, vd;       // length, mass m=δ·L, Var(δ) per live bin
        for (int l : live) {
            const BulkDamagePerBin& B = R.bins[l];
            int ng = static_cast<int>(B.profile_loglik.size());
            if (ng < 3 || static_cast<int>(B.profile_delta.size()) < 2) continue;
            // curvature at the profile peak: I = −d²ℓ/dδ² ⇒ Var(δ)=1/I (∞ if non-concave/boundary-flat)
            int kpk = 0;
            for (int g = 1; g < ng; ++g) if (B.profile_loglik[g] > B.profile_loglik[kpk]) kpk = g;
            int kc = std::clamp(kpk, 1, ng - 2);          // centre stencil; boundary peak → nearest interior
            double h  = B.profile_delta[1] - B.profile_delta[0];
            double d2 = (B.profile_loglik[kc + 1] - 2.0 * B.profile_loglik[kc] + B.profile_loglik[kc - 1])
                        / (h * h);
            double var = (h > 0.0 && d2 < 0.0) ? 1.0 / (-d2)
                                               : std::numeric_limits<double>::infinity();
            Lv.push_back(B.median_len);
            mv.push_back(B.delta * B.median_len);
            vd.push_back(var);
        }
        int n = static_cast<int>(Lv.size());
        if (n >= 2) {
            double Lbar = 0.0; for (double L : Lv) Lbar += L; Lbar /= n;
            double Sxx = 0.0;  for (double L : Lv) Sxx += (L - Lbar) * (L - Lbar);
            if (Sxx > 0.0) {
                double slope = 0.0, var_slope = 0.0;
                for (int i = 0; i < n; ++i) {
                    double ci = (Lv[i] - Lbar) / Sxx;                 // OLS weight on m_i
                    slope     += ci * mv[i];
                    var_slope += ci * ci * Lv[i] * Lv[i] * vd[i];     // Var(m_i)=L²·Var(δ_i)
                }
                R.slope_m = slope;
                double se = std::sqrt(var_slope);     // +∞ if any bin's δ is undetermined → w→0.5
                if (se > 0.0) R.w_length = sigmoid(-slope / (W_LENGTH_C * se));
            }
        }
    }

    // ── length-coupling DIAGNOSTIC (sign of slope δ_l vs median_len; NOT used in the fit) ────────────
    {
        std::vector<double> Lv, Dv;
        for (int l : live) { Lv.push_back(R.bins[l].median_len); Dv.push_back(R.bins[l].delta); }
        int n = static_cast<int>(Lv.size());
        if (n >= 2) {
            double Lbar = 0.0, Dbar = 0.0;
            for (int i = 0; i < n; ++i) { Lbar += Lv[i]; Dbar += Dv[i]; }
            Lbar /= n; Dbar /= n;
            double Sxx = 0.0, Sxy = 0.0;
            for (int i = 0; i < n; ++i) { Sxx += (Lv[i]-Lbar)*(Lv[i]-Lbar); Sxy += (Lv[i]-Lbar)*(Dv[i]-Dbar); }
            if (Sxx > 0.0) {
                R.length_coupling_slope = Sxy / Sxx;
                R.length_coupling = (R.length_coupling_slope < -1e-6) ? -1
                                  : (R.length_coupling_slope >  1e-6) ?  1 : 0;
            }
        }
    }

    // ── mixture P2: eligibility-conditioned ancient/modern split ─────────────────────────────────
    // Conditioning on the per-read damage-channel SITE counts (S5,S3) — damage-invariant — removes the
    // GC/composition confound, so Cov(k5,k3 | S5,S3) is pure ancient-mixture covariance. With δ_l pinned
    // each (bin,stratum) estimates the shared d_max via
    //     Cov_s = S5·S3·(d_max − δ_l)·δ_l·(1−η_l)²·ē² ,   ē = mean_p exp(−λp) over the JP positions,
    // pooled inverse-variance → d_max ± SE (threshold-free). π_l = δ_l/d_max. A pervasive artifact
    // (w_length ≤ 0.5) has no meaningful split, so the report is gated on the length-coupling weight.
    if (!s.jstrat.empty()) {
        const int JP = JStrat::JP;
        // ē = mean per-position deamination weight over the JP terminal co-occurrence positions. The P2
        // moments count damage at p=0..JP−1 with probability ∝ e^{−λp} at EVERY position (the JP window is
        // the very terminus where the overhang is single-stranded), so ē uses the bare exp law incl p0 —
        // NOT the ds fit-kernel that zeroes p0 (that exclusion is a terminal-artifact guard for the δ fit,
        // not a statement about the true damage probability at the terminal base).
        double ebar = 0.0;
        for (int p = 0; p < JP; ++p) ebar += std::exp(-R.lambda * p);
        ebar /= JP;
        double sum_w = 0.0, sum_wd = 0.0;
        for (int l = 0; l < L && l < static_cast<int>(s.jstrat.size()); ++l) {
            const JStrat& J = s.jstrat[l];
            // marginal diagnostic (pool all strata): joint_n / mean / cov (cov here is GC-confounded)
            uint64_t nt = 0; double sa = 0.0, sb = 0.0, sab = 0.0;
            for (int a = 0; a <= JP; ++a)
                for (int b = 0; b <= JP; ++b) {
                    nt += J.n[a][b]; sa += static_cast<double>(J.sk5[a][b]);
                    sb += static_cast<double>(J.sk3[a][b]); sab += static_cast<double>(J.sk53[a][b]);
                }
            if (nt > 0) {
                double dn = static_cast<double>(nt), m5 = sa / dn, m3 = sb / dn;
                R.bins[l].joint_n = nt;
                R.bins[l].joint_mean = {m5, m3};
                R.bins[l].joint_cov = sab / dn - m5 * m3;
            }
            // d_max contribution — signal bins only (δ→0 bins carry no split information)
            double dl = R.bins[l].delta;
            if (dl <= 1e-3) continue;
            double eta = s.n_interior[l][0] > 0
                ? static_cast<double>(s.k_interior[l][0]) / static_cast<double>(s.n_interior[l][0]) : 0.5;
            double om = 1.0 - eta;
            for (int a = 1; a <= JP; ++a)
                for (int b = 1; b <= JP; ++b) {
                    uint64_t ns = J.n[a][b];
                    if (ns < 200) continue;                          // too few reads for a stable cov
                    double dn = static_cast<double>(ns);
                    double m5 = static_cast<double>(J.sk5[a][b]) / dn;
                    double m3 = static_cast<double>(J.sk3[a][b]) / dn;
                    double cov = static_cast<double>(J.sk53[a][b]) / dn - m5 * m3;
                    double denom = static_cast<double>(a) * b * dl * om * om * ebar * ebar;
                    if (denom <= 1e-12) continue;
                    double dmax_s = dl + cov / denom;                // per-stratum d_max estimate
                    double v5 = m5 * (1.0 - m5 / a), v3 = m3 * (1.0 - m3 / b);   // conditional binom var
                    double var_dmax = (v5 * v3 / dn) / (denom * denom);
                    if (!(var_dmax > 0.0) || !std::isfinite(var_dmax)) continue;
                    double w = 1.0 / var_dmax;
                    sum_w += w; sum_wd += w * dmax_s;
                }
        }
        double dmax_raw = sum_w > 0.0 ? sum_wd / sum_w : -1.0;
        double dmax_floor = 0.0;
        for (int l = 0; l < L; ++l) dmax_floor = std::max(dmax_floor, R.bins[l].delta);
        R.d_max_raw = dmax_raw;
        // Identified only when the shared d_max is INTERIOR: above max δ_l (so π<1) and strictly below
        // the probability bound. A solution railed at d_max→1 means the conditional covariance exceeds
        // what any mixture can produce — residual composition structure that (S5,S3)-conditioning did
        // not remove — so the split is UNIDENTIFIED → undetermined (not a confident d_max=1). Gated by
        // w_length so a pervasive artifact never yields a split.
        if (sum_w > 0.0 && R.w_length > 0.5 && dmax_raw > dmax_floor && dmax_raw < 1.0 - 1e-3) {
            double dmax = dmax_raw;
            double dse  = std::sqrt(1.0 / sum_w);
            R.d_max_ancient = dmax;
            R.d_max_se = dse;
            for (int l = 0; l < L; ++l) {
                double dl = R.bins[l].delta;
                if (dl <= 1e-4) continue;
                R.bins[l].pi_ancient = std::min(1.0, dl / dmax);
                double hi = dmax + 1.96 * dse, lo = std::max(dl, dmax - 1.96 * dse);
                R.bins[l].pi_lo = std::min(1.0, dl / hi);
                R.bins[l].pi_hi = std::min(1.0, dl / lo);
            }
        }
    }

    return R;
}

} // namespace taph
