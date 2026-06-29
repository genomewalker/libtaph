// BulkDamageModel implementation (moved from bulk_damage_model.hpp).
#include "taph/bulk_damage_model.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <thread>
#include <vector>
namespace taph {
// ------------------------------------------------------------------ log-likelihood

double BulkDamageModel::log_lik(const BulkDamageSuffStats& s, const Params& P) {
    double ll = 0.0;
    const int L = s.L();
    const DecayTable decay = make_decay(P, s.ss);   // FIX B: 15 exp(-λp) once per evaluation
    for (int l = 0; l < L; ++l) {
        for (int c = 0; c < BulkDamageSuffStats::N_CH; ++c) {
            for (int e = 0; e < N_END; ++e) {
                const auto& cell = s.bin[l][c][e];
                for (int p = 0; p < N_POS; ++p) {
                    if (!cell_live(s, c, e, p)) continue;
                    uint64_t n = cell.n[p];
                    if (n == 0) continue;
                    double m = mu_d(P, decay, c, e, l, p);
                    double k = static_cast<double>(cell.k[p]);
                    ll += k * std::log(m) + (static_cast<double>(n) - k) * std::log(1.0 - m);
                }
            }
        }
    }
    return ll;
}

// FIX B: build the λ-independent eta cache for every cell of every bin (control + damage).
// eta() depends only on β̂,a,spike — none of which the golden search touches — so this is computed
// once before the search and reused across all ~40 evaluations. Values are the exact eta() bits.
void BulkDamageModel::build_eta_cache(const BulkDamageSuffStats& s, const Params& P,
                                      const std::vector<int>& /*live*/, EtaCache& ec) {
    const int L = s.L();
    ec.resize(L);
    for (int l = 0; l < L; ++l)
        for (int c = 0; c < BulkDamageSuffStats::N_CH; ++c)
            for (int e = 0; e < N_END; ++e)
                for (int p = 0; p < N_POS; ++p)
                    ec[l][c][e][p] = eta(P, c, e, l, p);
}

// FIX B: log_lik reading eta from the cache and r from the decay table. Same loop nesting, same
// per-cell formula, same accumulation order ⇒ bit-identical to log_lik().
double BulkDamageModel::log_lik_cached(const BulkDamageSuffStats& s, const Params& P,
                                       const EtaCache& ec, const DecayTable& decay) {
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
                    double h = ec[l][c][e][p];
                    double m;
                    if (c == 1) {
                        m = clamp_mu(h);
                    } else {
                        double r = decay[p];
                        m = (r == 0.0) ? clamp_mu(h)
                                       : clamp_mu(h + (1.0 - h) * P.delta[l] * r);
                    }
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
double BulkDamageModel::penalized_obj(const BulkDamageSuffStats& s, const Params& P) {
    double obj = log_lik(s, P);
    const int L = s.L();
    for (int l = 0; l < L && l < static_cast<int>(P.s_spike.size()); ++l)
        for (int e = 0; e < N_END; ++e) obj -= SPIKE_RIDGE_N * P.s_spike[l][e] * P.s_spike[l][e];
    return obj;
}

// ------------------------------------------------------------------ init

void BulkDamageModel::init(const BulkDamageSuffStats& s, Params& P,
                                  const std::vector<int>& live) {
    const int L = s.L();
    P.beta.assign(L, {0.0, 0.0});
    P.delta.assign(L, 0.0);
    P.lambda = 0.3;
    P.ss_p0_overhang = s.ss_p0_overhang;                // Wave-3: carry the ss-overhang gate into r_kernel
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

void BulkDamageModel::step_beta(const BulkDamageSuffStats& s, Params& P,
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
void BulkDamageModel::step_spike(const BulkDamageSuffStats& s, Params& P,
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
void BulkDamageModel::step_artifact(const BulkDamageSuffStats& s, Params& P,
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
double BulkDamageModel::delta_block_mle(const BulkDamageSuffStats& s, const Params& P,
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
void BulkDamageModel::step_delta_isotonic(const BulkDamageSuffStats& s, Params& P,
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

    // FIX C: reuse one cell-scratch buffer across all block_delta calls (clear+reuse, reserved once)
    // instead of allocating a std::vector<LE> per call. Cells are pushed in the SAME order, so the
    // delta_block_mle sum over them is bit-identical.
    std::vector<LE> cells;
    cells.reserve(static_cast<std::size_t>(N) * 2);
    auto block_delta = [&](const std::vector<int>& bins) -> double {
        // Forced anchor: if the block contains lmax_idx, return 0.
        for (int bl : bins) if (bl == lmax_idx) return 0.0;
        // Profile pin: if the block contains pinned_idx, return pinned_val.
        for (int bl : bins) if (bl == pinned_idx) return pinned_val;
        // Collect cells for all bins in the block (both ends, pooled).
        cells.clear();
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

void BulkDamageModel::step_lambda(const BulkDamageSuffStats& s, Params& P,
                                         const std::vector<int>& live) {
    const double phi = 0.6180339887498949;
    double a = LAMBDA_MIN, b = LAMBDA_MAX;
    double c = b - phi * (b - a), d = a + phi * (b - a);
    // FIX B: eta is λ-independent — cache it once; each f(lam) only rebuilds the 15-entry decay table.
    EtaCache ec;
    build_eta_cache(s, P, live, ec);
    auto f = [&](double lam) {
        double sv = P.lambda; P.lambda = lam;
        DecayTable decay = make_decay(P, s.ss);
        double v = log_lik_cached(s, P, ec, decay);
        P.lambda = sv; return v;
    };
    double fc = f(c), fd = f(d);
    for (int it = 0; it < 40; ++it) {
        if (fc > fd) { b = d; d = c; fd = fc; c = b - phi * (b - a); fc = f(c); }
        else        { a = c; c = d; fc = fd; d = a + phi * (b - a); fd = f(d); }
    }
    P.lambda = 0.5 * (a + b);
}

// ------------------------------------------------------------------ convergence helpers

double BulkDamageModel::max_param_change(const Params& A, const Params& B,
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
int BulkDamageModel::run_sweeps(const BulkDamageSuffStats& s, Params& P,
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
        // Stop only when both the parameters and likelihood are flat.  Individual
        // coordinate updates can backtrack to the saved value, so a tiny net parameter
        // move alone is not reliable convergence evidence on the near-zero-delta ridge.
        bool ll_flat   = (ll - ll_prev) <= TOL_LL_ABS + TOL_LL_REL * std::fabs(ll);
        bool par_small = max_param_change(old, P, live) < TOL_PARAM;
        if (par_small && ll_flat) { ++sweep; break; }
        ll_prev = ll;
    }
    return sweep;
}

// ------------------------------------------------------------------ driver

BulkDamageResult BulkDamageModel::fit(const BulkDamageSuffStats& s, int n_threads) {
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
    R.lambda_at_boundary = (P.lambda >= LAMBDA_MAX - 1e-9);
    R.ss_overhang_modeled = s.ss && s.ss_p0_overhang;   // Wave-3: did r_kernel model the 5' ss overhang?
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
        {
            const uint64_t ni = s.n_interior[l][0], ki = s.k_interior[l][0];
            B.interior_baseline = ni > 0 ? static_cast<double>(ki) / static_cast<double>(ni) : 0.0;
        }
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
    // The body is read-only on (s, P, R.log_lik) and writes ONLY R.bins[l] for its own l with a
    // PRIVATE Params Q — zero cross-bin dependency. Distribute over a fixed worker pool by index into
    // `live`; each result lands in its fixed R.bins[l] slot, so output is completion-order-independent
    // ⇒ byte-identical to the serial loop. Same pattern as LengthBinStats::finalize_all.
    auto profile_bin = [&](int l) {
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
    };
    const int workers = std::min<int>(n_threads < 1 ? 1 : n_threads,
                                      static_cast<int>(live.size()));
    if (workers <= 1) {
        for (int l : live) profile_bin(l);
    } else {
        std::atomic<size_t> next{0};
        auto run = [&] {
            for (;;) {
                size_t i = next.fetch_add(1, std::memory_order_relaxed);
                if (i >= live.size()) break;
                profile_bin(live[i]);
            }
        };
        std::vector<std::thread> pool;
        pool.reserve(static_cast<size_t>(workers) - 1);
        for (int t = 1; t < workers; ++t) pool.emplace_back(run);
        run();
        for (auto& th : pool) th.join();
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

    // ── mixture P2: eligibility-conditioned damaged/non-damaged split ─────────────────────────────────
    // Conditioning on the per-read damage-channel SITE counts (S5,S3) — damage-invariant — removes the
    // GC/composition confound, so Cov(k5,k3 | S5,S3) is pure damaged-mixture covariance. With δ_l pinned
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
        // C1: railed when non-physical (unset <0, or >1 ⇒ conditional cov exceeds any mixture). Gates the
        // JSON emitter to null instead of letting a >1 value be read as a rate.
        R.d_max_raw_railed = !(dmax_raw >= 0.0 && dmax_raw <= 1.0);
        // Identified only when the shared d_max is INTERIOR: above max δ_l (so π<1) and strictly below
        // the probability bound. A solution railed at d_max→1 means the conditional covariance exceeds
        // what any mixture can produce — residual composition structure that (S5,S3)-conditioning did
        // not remove — so the split is UNIDENTIFIED → undetermined (not a confident d_max=1). Gated by
        // w_length so a pervasive artifact never yields a split.
        if (sum_w > 0.0 && R.w_length > 0.5 && dmax_raw > dmax_floor && dmax_raw < 1.0 - 1e-3) {
            double dmax = dmax_raw;
            double dse  = std::sqrt(1.0 / sum_w);
            R.d_max_damaged = dmax;
            R.d_max_se = dse;
            R.d_max_damaged_valid = true;   // C1: split identified ⇒ d_max_damaged/d_max_se are real
            for (int l = 0; l < L; ++l) {
                double dl = R.bins[l].delta;
                if (dl <= 1e-4) continue;
                R.bins[l].pi_damaged = std::min(1.0, dl / dmax);
                double hi = dmax + 1.96 * dse, lo = std::max(dl, dmax - 1.96 * dse);
                R.bins[l].pi_lo = std::min(1.0, dl / hi);
                R.bins[l].pi_hi = std::min(1.0, dl / lo);
                R.bins[l].pi_identified = true;   // C1: this bin's π_l interval holds real probabilities
            }
        }
    }

    return R;
}

} // namespace taph
