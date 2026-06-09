#pragma once
// Oxidation co-movement score: reference-free G→T + C→A ancient-vs-modern contrast.
//
// Usage pattern (mirrors ancient_fraction.hpp):
//   1. Allocate OxBinAcc bins[N_GC_BINS] = {} per worker thread.
//   2. Call update_ox_bins(bins[gc_bin], q, T, G, Cv, Av) for each read's middle-third.
//   3. Merge across threads with merge_ox_bins(dst, src).
//   4. Call compute_ox_scores(merged_bins, N_GC_BINS) for S_oxog, s_ca, D_oriented.
//
// Interpretation:
//   s_oxog > 0 AND s_ca > 0 → oxidation-consistent co-movement (G→T and its complement
//   C→A both enriched in ancient-weighted reads vs modern-weighted reads).
//   Neither positive → no reference-free oxidation signal detected.
#include <cmath>
#include <cstdint>
#include <type_traits>

namespace taph {

// Per-GC-bin accumulator.  Trivially copyable; holds q-weighted moments for G→T
// (s_oxog) and C→A (s_ca) channels plus unweighted totals for D_oriented.
struct OxBinAcc {
    double A=0,B=0,C=0,D=0;        // q-weighted G→T: A=Σq·T, B=Σq·(T+G)
    double AA=0,AB=0,BB=0;          // second moments (ancient G→T)
    double CC=0,CD=0,DD=0;          // second moments (modern G→T)
    double AC=0,AD=0,BC=0,BD=0;     // cross moments
    double M_total=0;
    uint64_t reads=0;
    double T_all=0,TG_all=0;        // unweighted T, T+G (for D_oriented)
    double C_all=0,CA_all=0;        // unweighted C, C+A (for D_oriented)
    double A_ca=0,B_ca=0;           // q-weighted C→A ancient: A_ca=Σq·Cv, B_ca=Σq·(Cv+Av)
    double C_ca=0,D_ca=0;           // q-weighted C→A modern
};
static_assert(std::is_trivially_copyable_v<OxBinAcc>);

struct OxScoreResult {
    double s_oxog     = 0.0;  // ancient G→T rate minus modern (co-movement numerator)
    double se_s_oxog  = 0.0;  // standard error of s_oxog
    double s_ca       = 0.0;  // ancient C→A rate minus modern (co-movement complement)
    double d_oriented = 0.0;  // log(T·A/G·C) ilr contrast (log-odds D_oriented)
    bool   has_score  = false;
};

// Hot per-read update.  q = ancient weight in [0,1].
// T,G,Cv,Av = orientation-corrected middle-third base counts (T+G must be > 0).
inline void update_ox_bins(OxBinAcc& x, double q, int T, int G, int Cv, int Av) {
    double M   = static_cast<double>(T + G);
    double a   = q*T,         b_ = q*M;
    double c   = (1.0-q)*T,   d  = (1.0-q)*M;
    x.A += a;  x.B += b_;  x.C += c;  x.D += d;
    x.AA += a*a; x.AB += a*b_; x.BB += b_*b_;
    x.CC += c*c; x.CD += c*d;  x.DD += d*d;
    x.AC += a*c; x.AD += a*d;  x.BC += b_*c; x.BD += b_*d;
    x.M_total += M; ++x.reads;
    x.T_all += T; x.TG_all += M;
    x.C_all += Cv; x.CA_all += Cv + Av;
    double a_ca=q*Cv, b_ca=q*(Cv+Av), c_ca=(1.0-q)*Cv, d_ca=(1.0-q)*(Cv+Av);
    x.A_ca += a_ca; x.B_ca += b_ca; x.C_ca += c_ca; x.D_ca += d_ca;
}

// Additive merge (for cross-thread reduction).
inline void merge_ox_bins(OxBinAcc& dst, const OxBinAcc& s) {
    dst.A+=s.A; dst.B+=s.B; dst.C+=s.C; dst.D+=s.D;
    dst.AA+=s.AA; dst.AB+=s.AB; dst.BB+=s.BB;
    dst.CC+=s.CC; dst.CD+=s.CD; dst.DD+=s.DD;
    dst.AC+=s.AC; dst.AD+=s.AD; dst.BC+=s.BC; dst.BD+=s.BD;
    dst.M_total+=s.M_total; dst.reads+=s.reads;
    dst.T_all+=s.T_all; dst.TG_all+=s.TG_all;
    dst.C_all+=s.C_all; dst.CA_all+=s.CA_all;
    dst.A_ca+=s.A_ca; dst.B_ca+=s.B_ca; dst.C_ca+=s.C_ca; dst.D_ca+=s.D_ca;
}

// Pure derivation over merged per-GC-bin accumulators.
// n_bins is typically N_GC_BINS (10).
OxScoreResult compute_ox_scores(const OxBinAcc* bins, int n_bins);

} // namespace taph
