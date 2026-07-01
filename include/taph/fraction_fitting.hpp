#pragma once
#include <array>
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

// CT-specific decay fit: fits the decay model directly on the CT-vs-(CA,CG)
// difference series, replacing the external "interior baseline" estimate
// (pool_interior_bg) entirely.
//
// Motivation: pool_interior_bg assumes positions past a fixed window (default
// 8) are damage-free. They are not, in proportion to true damage — ssDNA
// overhang length scales with deamination extent, so more-damaged samples
// leak real signal deep into the "interior," inflating the externally
// estimated baseline and mechanically deflating fit_exp_decay_irls's d_max in
// proportion to the very quantity being measured (validated on a 31-sample
// metaDMG benchmark: corr(true damage, baseline) = +0.64; corr(fitted d_max,
// baseline) = -0.58; magnitude correlation with truth was NEGATIVE, r=-0.49,
// for genuinely damaged samples). Two prior attempts to patch this — (1) an
// adaptive/deeper baseline window driven by a first-pass lambda estimate, and
// (2) an unconstrained joint fit of baseline+amplitude+lambda — both failed
// validation: (1) is circular (the correcting lambda is itself biased by the
// contamination it's meant to correct), (2) is numerically unstable (rejects
// genuine strong positives) and not CT-specific (refits the same generic
// terminal-composition artifact the CT-vs-CA/CG specificity gate exists to
// reject).
//
// This function sidesteps the whole problem: instead of estimating a
// baseline from CT alone, it fits on
//   g(p) = rate_CT(p) - 0.5*(rate_CA(p) + rate_CG(p))
// CA and CG are non-C->T substitutions at the same terminal positions. Any
// position-dependent effect that is NOT deamination-specific — including the
// overhang-driven composition drift that contaminates pool_interior_bg's
// window — hits CT, CA, and CG equally and cancels in the difference. True
// deamination is CT-specific and survives as a clean decaying signal in g(p),
// so no external baseline estimate (and nothing for damage to contaminate) is
// needed at all. amplitude = fitted intercept (terminal excess at position 1,
// where decay curves of any lambda coincide), so a strong true signal is read
// off directly rather than compressed against an external reference.
//
// Validated (hand-check, 9 documented benchmark samples): 4 false-positive/
// low-truth samples (metaDMG truth 0.005-0.111) fail to fit (fewer than 3
// positions clear the positivity filter) or fit near-zero; 5 true-positive
// samples (truth 0.451-0.558) fit cleanly with Pearson r=+0.65 against truth
// (prior methods: r=-0.34 to -0.49).
//
// tri_pos[p] holds trinucleotide-context substitution counts at position p:
// tri_pos[p][prev*16 + base*4 + next], base in {A=0,C=1,G=2,T=3}. Marginal
// rate for substitution FROM->TO at position p is aggregated over all
// prev/next flanking bases.
struct CtSpecificDecayFit { double A; double lambda; int n_points; };
CtSpecificDecayFit fit_ct_specific_decay(
    const std::array<uint64_t, 64>* tri_pos, int n_pos);

} // namespace taph
