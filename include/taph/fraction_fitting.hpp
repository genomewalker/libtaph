#pragma once
#include <cstdint>
#include <utility>

namespace taph {

// Estimate interior background rate from count arrays t[]/tc[].
// Returns mean T/TC over positions [start_pos, n_pos) with tc[p] >= 50.
// Falls back to fallback_bg when fewer than 2 qualifying positions.
double pool_interior_bg(const int64_t* t, const int64_t* tc,
                        int n_pos, double fallback_bg,
                        int start_pos = 8);

// Free-lambda exponential decay fit via 2-pass IRLS.
// Fits excess(p) = d_max * exp(-lambda*p) where excess = T/TC(p) - bg,
// over positions p=1..6 (pos 0 always skipped: adapter composition artifact).
// Returns {d_max, lambda}, or {NaN, NaN} when the fit fails.
// skip_pos0: when true, return d(pos1)=d_max*exp(-lambda) rather than
//   the extrapolated d(pos0). Use for SS or blunted-end libraries.
std::pair<double,double> fit_exp_decay_irls(
    const int64_t* t, const int64_t* tc, int n_pos,
    double bg, bool skip_pos0 = false);

} // namespace taph
