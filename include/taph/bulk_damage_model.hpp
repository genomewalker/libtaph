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
// removes the GC/composition confound, so Cov(k5,k3 | S5,S3) is pure damaged-mixture covariance and
// the damaged fraction π separates from d_max. JP=5 terminal positions per end; k = #damaged sites.
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

    bool ss = false;            // single-stranded library: enables the 5' single-strand overhang kernel
    bool skip_3p_pos0 = true;   // ss ligation artifact: drop 3′ p=0 from the damage score
    bool ss_p0_overhang = false; // Wave-3: model the 5' terminus (r(0)=1) as genuine ss overhang — set iff
                                 // ss AND p0 is not an adapter/composition artifact (identifiable overhang)

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
    double interior_baseline = 0.0; // g_l: raw interior C→T rate (k_int/n_int); δ_l = terminal_rate − g_l
    bool   identified = false;     // δ̂_l estimated independently (singleton PAVA block)
    bool   borrowed = false;       // PAVA pooled with at least one neighbour (block_size > 1)

    // threshold-free diagnostics, per end (0=5′,1=3′)
    std::array<double, 2> k_eff   = {0.0, 0.0};   // effective informative positions (IPR)
    // (D_ctrl−D_full)/D_ctrl on the damage channel: fraction of the δ≡0 (artifact-only) deviance explained
    // by adding terminal deamination. Normally ∈[0,1] (1 ⇒ δ fully explains the terminal excess; 0 ⇒ no
    // improvement). The 0.0 here is BOTH the genuine "no improvement" value AND the not-computed default
    // (no live cells / D_ctrl≤1e-12) — a consumer cannot distinguish them from this field alone; read it
    // together with k_eff[e]>0 and identified/delta_ell0 to confirm the bin actually carried signal.
    std::array<double, 2> r_damage = {0.0, 0.0};
    double delta_ell0 = 0.0;       // ℓ̂_l − ℓ(δ_l=0), nuisances re-profiled (full profile)
    double s0 = 1.0;               // exp(−delta_ell0) = L(δ_l=0)/L(δ̂_l)

    // FULL profile-likelihood curve over δ_l (nuisances re-optimized per grid point)
    std::vector<double> profile_delta;    // grid
    std::vector<double> profile_loglik;   // ℓ(δ_l=grid) − ℓ̂  (≤ 0, peak 0 at δ̂_l)

    // mixture P2 — eligibility-conditioned damaged/non-damaged split (per-read two-end co-occurrence)
    uint64_t joint_n = 0;                       // reads contributing to the co-occurrence moments
    std::array<double, 2> joint_mean = {0.0, 0.0};  // marginal E[k5], E[k3] (diagnostic)
    double joint_cov = 0.0;                     // marginal Cov(k5,k3) (diagnostic; GC-confounded)
    double pi_damaged = -1.0;                   // damaged read fraction π_l = δ_l/d_max (−1 = undetermined)
    double pi_lo = -1.0, pi_hi = -1.0;          // π_l 95% interval (threshold-free; spans → undetermined)
    // C1 cross-file flag: true ONLY when the mixture split was identified for this bin and pi_damaged/
    // pi_lo/pi_hi hold real probabilities in [0,1]. The −1.0 sentinels are internal-only; profile_json.cpp
    // gates on this flag to emit JSON null instead of −1.0. Default false (reset in fit()).
    bool   pi_identified = false;
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
    // BEHAVIORAL CHANGE (C5): default 0.5 (undetermined), not 1.0. The slope-based update below is
    // gated on n>=2 live bins; with exactly 1 live bin the whole block is skipped and w_length must
    // read as undetermined (w≈0.5), NOT full authenticity (1.0). 1.0 over-certified headline_delta_auth
    // and let the pi_damaged gate (w_length>0.5) pass on a single bin where the slope is undefined.
    double w_length = 0.5;

    // Length-coupling DIAGNOSTIC (post-fit, NOT used in the fit). Sign of the slope of the recovered δ_l
    // against median read length over live bins: −1 ⇒ δ decreasing with length (classic length-coupled
    // terminal deamination, e.g. FLB08m), +1 ⇒ increasing, 0 ⇒ flat/undetermined.
    double length_coupling_slope = 0.0;   // OLS slope of δ_l vs median_len over live bins
    int    length_coupling = 0;           // sign(slope): −1 decreasing, +1 increasing, 0 flat

    // mixture P2 — shared per-damaged-read deamination intensity d_max (reference-free analog of
    // metaDMG A_b), estimated from the eligibility-conditioned two-end co-occurrence with δ_l pinned
    // (π_l = δ_l/d_max). d_max = −1 ⇒ undetermined (no length bin carried a usable co-occurrence signal,
    // e.g. δ≈0 everywhere or the artifact gate fired). Threshold-free: d_max_se carries the data's own
    // uncertainty; the per-bin π interval spans → undetermined where the signal can't fix the split.
    double d_max_damaged = -1.0;
    double d_max_se      = 0.0;
    double d_max_raw     = -1.0;   // pre-gate pooled estimate; >1 (non-physical) ⇒ railed ⇒ unidentified
    // C1 cross-file flags (profile_json.cpp reads these to emit JSON null instead of the sentinels):
    //   d_max_damaged_valid: true ONLY when the split was identified; then d_max_damaged∈(δmax,1) and
    //     d_max_se is a real SE. When false, d_max_damaged=−1.0 and d_max_se=0.0 are NOT-COMPUTED
    //     sentinels (se=0.0 must NOT read as zero uncertainty). Default false (reset in fit()).
    //   d_max_raw_railed: true when d_max_raw is non-physical (<0 unset, or >1 railed). When true the
    //     raw value must not be emitted as a rate. Default false (reset in fit()).
    bool   d_max_damaged_valid = false;
    bool   d_max_raw_railed    = false;
    bool   lambda_at_boundary  = false;  // C4: true when lambda == LAMBDA_MAX (decay rate unidentifiable)
    bool   ss_overhang_modeled = false;  // Wave-3: ss 5' terminal overhang included in the kernel (r(0)=1);
                                         // false ⇒ ds library, or ss whose p0 overhang was not identifiable
};

// ----------------------------------------------------------------------------- model / solver

class BulkDamageModel {
public:
    // n_threads parallelizes ONLY the per-bin profile-likelihood loop (each bin owns its R.bins[l]
    // slot and a private Params Q ⇒ output is index-ordered and bit-identical to the serial run).
    // Default 1 preserves every existing caller; the length-strat per-bin fits MUST keep 1 (their
    // bins are already distributed across a worker pool — no nested oversubscription).
    static BulkDamageResult fit(const BulkDamageSuffStats& s, int n_threads = 1);

    // Terminal (p=0) kernel weight — the one place the ds and ss decay kernels DIFFER.
    //   ds: the 5' terminus is excluded from the decay (0); the p0 spike s[l][e] absorbs any
    //       terminus-only adapter/composition artifact there.
    //   ss: the 5' terminus is genuine single-strand overhang — the maximally exposed, fully
    //       deaminated base — so it carries the full kernel weight r(0)=exp(0)=1, BUT only when
    //       the overhang is identifiable (ss library AND p0 not flagged as an artifact). When not
    //       identifiable the ss kernel falls back to the ds form (0) and ss_overhang_degenerate is
    //       logged. Public + branch-free so it is directly unit-testable (ds≠ss).
    static double terminal_kernel_weight(bool ss, bool ss_p0_overhang) {
        return (ss && ss_p0_overhang) ? 1.0 : 0.0;
    }

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
        bool                               ss_p0_overhang = false;  // Wave-3: ss 5' overhang at p0 (r(0)=1)
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

    // Damage decay kernel r(p;λ). r(p)=exp(−λ·p) for p≥1; at p=0 the ds/ss terminus split applies
    // (terminal_kernel_weight). The interior decay (p≥1) is identical for ds and ss.
    static double r_kernel(const Params& P, bool ss, int p) {
        return (p >= 1) ? std::exp(-P.lambda * p)
                        : terminal_kernel_weight(ss, P.ss_p0_overhang);
    }

    // FIX B: λ is fixed within a single log_lik / golden-search f() evaluation, so the 15 decay
    // weights r(p) are constant across all cells of that evaluation. Precompute them ONCE per λ and
    // reuse — each entry is the EXACT bits std::exp(-P.lambda*p) (or terminal_kernel_weight at p0)
    // that r_kernel would return, so every consumer reads identical values (no math reordering).
    using DecayTable = std::array<double, N_POS>;
    static DecayTable make_decay(const Params& P, bool ss) {
        DecayTable d;
        d[0] = terminal_kernel_weight(ss, P.ss_p0_overhang);
        for (int p = 1; p < N_POS; ++p) d[p] = std::exp(-P.lambda * p);
        return d;
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
    // μ using a precomputed decay table (FIX B). r comes from decay[p] — identical bits to mu() above.
    static double mu_d(const Params& P, const DecayTable& decay, int c, int e, int l, int p) {
        double h = eta(P, c, e, l, p);
        if (c == 1) return clamp_mu(h);
        double r = decay[p];
        if (r == 0.0) return clamp_mu(h);
        return clamp_mu(h + (1.0 - h) * P.delta[l] * r);
    }

    static bool cell_live(const BulkDamageSuffStats& s, int c, int e, int p) {
        // drop 3′ position 0 from the damage score when requested (ss ligation artifact)
        return !(s.skip_3p_pos0 && c == 0 && e == 1 && p == 0);
    }

    static double log_lik(const BulkDamageSuffStats& s, const Params& P);

    // FIX B (golden search): eta=σ(β̂+a+spike) is LAMBDA-INDEPENDENT. step_lambda's golden search
    // varies ONLY λ across ~40 log_lik calls, so cache eta per live cell once and reuse — only the
    // decay table changes per λ. EtaCache[l][c][e][p] holds the exact sigmoid bits eta() returns.
    using EtaCache = std::vector<std::array<std::array<std::array<double, N_POS>, N_END>, 2>>;
    static void build_eta_cache(const BulkDamageSuffStats& s, const Params& P,
                                const std::vector<int>& live, EtaCache& ec);
    // log_lik with eta read from the cache and r from the decay table — same per-cell arithmetic,
    // same summation order, identical bits to log_lik().
    static double log_lik_cached(const BulkDamageSuffStats& s, const Params& P,
                                 const EtaCache& ec, const DecayTable& decay);

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

} // namespace taph
