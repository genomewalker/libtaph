#pragma once
#include <cstdint>
#include <utility>

namespace taph {

// Interior background estimate from count arrays t[]/tc[].
// mean = unweighted mean of T/TC over positions [start_pos, n_pos) with tc[p] >= 50.
// var  = Var(mean) = (1/m^2) * sum_i rate_i*(1-rate_i)/tc_i  (propagated binomial variance).
// Falls back to {fallback_bg, 0.0} when fewer than 2 qualifying positions.
//
// NOTE: bg is subtracted from every fitting position in fit_exp_decay_irls, so the
// shared bg uncertainty induces a rank-1 covariance diag(var_rate) + var_bg*J across
// fitted positions. fit_exp_decay_irls uses that covariance in its second pass via
// a Sherman-Morrison rank-1 GLS update rather than treating bg.var as independent
// per-position noise.
struct BgEstimate { double mean; double var; };
BgEstimate pool_interior_bg(const int64_t* t, const int64_t* tc,
                             int n_pos, double fallback_bg,
                             int start_pos = 8);

// Free-lambda exponential decay fit via 2-pass IRLS.
// Fits excess(p) = d_max * exp(-lambda*p) where excess = T/TC(p) - bg.mean,
// over positions p=1..6 (pos 0 always skipped: adapter composition artifact).
// Returns {d_max, lambda}, or {NaN, NaN} when the fit fails.
// skip_pos0: when true, return d(pos1)=d_max*exp(-lambda) rather than
//   the extrapolated d(pos0). Use for SS or blunted-end libraries.
//
// The first pass uses Var(rate[p]) + Var(bg) as a conservative detection gate.
// The second pass fits log-excess with the full rank-1 covariance induced by the
// shared background estimate.
std::pair<double,double> fit_exp_decay_irls(
    const int64_t* t, const int64_t* tc, int n_pos,
    BgEstimate bg, bool skip_pos0 = false);

} // namespace taph
