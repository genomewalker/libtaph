// Damage scoring: CpG, oxoG interior, trinuc oxoG, depurination, damage mask.
#include "taph/library_interpretation.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>

namespace taph {

// Display cap for exploratory binomial z-scores. Detection gates use z>3, well
// inside ±12, so behaviour is unchanged; consumers no longer see absurd
// magnitudes from sqrt(N)-inflated, correlated-read Bernoulli statistics.
static constexpr double kExploratoryZCap = 12.0;
static constexpr double kDepurZCap       = kExploratoryZCap;
HexStats compute_hex_stats(const SampleDamageProfile& dp) {
    HexStats r;
    double tot_t = static_cast<double>(dp.n_hexamers_5prime);
    double tot_i = static_cast<double>(dp.n_hexamers_interior);

    // Entropy + JSD
    if (tot_t > 0 && tot_i > 0) {
        double h_mix = 0.0;
        for (int i = 0; i < 4096; ++i) {
            double p = dp.hexamer_count_5prime[i]  / tot_t;
            double q = dp.hexamer_count_interior[i] / tot_i;
            double m = 0.5 * (p + q);
            if (p > 0) r.entropy_terminal -= p * std::log2(p);
            if (q > 0) r.entropy_interior -= q * std::log2(q);
            if (m > 0) h_mix              -= m * std::log2(m);
        }
        r.jsd = h_mix - 0.5 * (r.entropy_terminal + r.entropy_interior);
    }

    // Multinomial G-test (terminal vs interior)
    if (tot_t > 0 && tot_i > 0) {
        double N_hex = tot_t + tot_i;
        double G_stat = 0.0;
        int k_eff = 0;
        for (int i = 0; i < 4096; ++i) {
            double c = static_cast<double>(dp.hexamer_count_5prime[i]);
            double d = static_cast<double>(dp.hexamer_count_interior[i]);
            if (c + d < 5.0) continue;
            double Ec = tot_t * (c + d) / N_hex;
            double Ed = tot_i * (c + d) / N_hex;
            if (c > 0) G_stat += 2.0 * c * std::log(c / Ec);
            if (d > 0) G_stat += 2.0 * d * std::log(d / Ed);
            ++k_eff;
        }
        r.shift_g = G_stat;
        if (k_eff > 1) {
            double df = static_cast<double>(k_eff - 1);
            r.shift_z = (G_stat - df) / std::sqrt(2.0 * df);
            r.shift_p = 0.5 * std::erfc(r.shift_z / std::sqrt(2.0));
        }
    }
    return r;
}

// ── Scores ────────────────────────────────────────────────────────────────────

CpgScore compute_cpg_score(const SampleDamageProfile& dp) {
    CpgScore r;
    if (std::isnan(dp.log2_cpg_ratio)) return r;
    double ecpg  = static_cast<double>(dp.effcov_ct5_cpg_like_terminal);
    double encpg = static_cast<double>(dp.effcov_ct5_noncpg_like_terminal);
    double se = std::sqrt(1.0 / (ecpg + 1.0) + 1.0 / (encpg + 1.0));
    if (se > 1e-9) {
        r.z = static_cast<double>(dp.log2_cpg_ratio) / se;
        r.p = std::erfc(std::abs(r.z) / std::sqrt(2.0));
        r.computed = true;
    }
    return r;
}

OxogInteriorScore compute_oxog_interior_score(const SampleDamageProfile& dp) {
    OxogInteriorScore r;
    static constexpr int RC4[4] = {3, 2, 1, 0};
    double sc = 0.0, vr = 0.0;
    for (const auto* tri : {&dp.tri_5prime_interior, &dp.tri_3prime_interior}) {
        for (int p = 0; p < 4; ++p) {
            for (int n = 0; n < 4; ++n) {
                double k  = (*tri)[p*16 + 3*4 + n];
                double g  = (*tri)[p*16 + 2*4 + n];
                double nt = k + g;
                if (nt < 10.0) continue;
                int rp = RC4[n], rn = RC4[p];
                double a_rc = (*tri)[rp*16 + 0*4 + rn];
                double c_rc = (*tri)[rp*16 + 1*4 + rn];
                double ca   = a_rc + c_rc;
                if (ca < 10.0) continue;
                double theta = a_rc / ca;
                sc += k - nt * theta;
                vr += nt * theta * (1.0 - theta);
            }
        }
    }
    if (vr > 0.0) {
        // exploratory; clamped, not a calibrated p-value (correlated reads).
        // Raw Binomial z scales ~sqrt(read count); gt_rate (compute_oxog_trinuc)
        // is the primary effect-size measure. Clamp magnitude and floor p so the
        // erfc underflow-to-0 case stays distinguishable from a genuine zero.
        r.z = std::max(-kExploratoryZCap,
                       std::min(kExploratoryZCap, sc / std::sqrt(vr)));
        r.p = std::max(0.5 * std::erfc(r.z / std::sqrt(2.0)),
                       std::numeric_limits<double>::min());
    }
    return r;
}

OxogTrinucResult compute_oxog_trinuc(const SampleDamageProfile& dp) {
    static constexpr double OXOG_REF[16] = {
        0.0512, 0.0388, 0.0257, 0.0601,
        0.0788, 0.0621, 0.0302, 0.1124,
        0.0389, 0.0296, 0.0198, 0.0453,
        0.0863, 0.0682, 0.0388, 0.1138,
    };
    static constexpr int RC4[4] = {3, 2, 1, 0};
    OxogTrinucResult r;
    double v[16] = {};
    // F1: per-class strand-symmetric excess. exc[c] accumulates (k - nt*theta) — the
    // G→T-minus-C→A-complement residual already corrected by theta (the same oxo_two_marker
    // correction used for v[]); opp[c] accumulates nt (eligible opportunities). Classes:
    // 5'-flanker signature {oh_radical=C/T (•OH), m1g=A, hole=G} and G-run {GGG,GG,G} keyed on whether
    // the G has a 5'-G and/or 3'-G neighbour (single-trinuc approximation of run length).
    enum { CL_OH, CL_M1G, CL_HOLE, CL_GGG, CL_GG, CL_G, N_CL };
    double exc[N_CL] = {}, opp[N_CL] = {};
    for (const auto* tri : {&dp.tri_5prime_interior, &dp.tri_3prime_interior}) {
        for (int p = 0; p < 4; ++p) {
            for (int n = 0; n < 4; ++n) {
                double k  = static_cast<double>((*tri)[p*16 + 3*4 + n]);
                double g  = static_cast<double>((*tri)[p*16 + 2*4 + n]);
                double nt = k + g;
                if (nt < 10.0) continue;
                int rp = RC4[n], rn = RC4[p];
                double a_rc = static_cast<double>((*tri)[rp*16 + 0*4 + rn]);
                double c_rc = static_cast<double>((*tri)[rp*16 + 1*4 + rn]);
                double ca = a_rc + c_rc;
                if (ca < 10.0) continue;
                double theta = a_rc / ca;
                double resid = k - nt * theta;   // signed strand-symmetric excess count
                v[p*4 + n] += std::max(0.0, resid);
                ++r.n_ctx;
                // 5'-flanker mechanism (p is the 5' neighbour of the central G).
                int mech = (p == 2) ? CL_HOLE : (p == 0) ? CL_M1G : CL_OH;  // C/T -> oh_radical (•OH signature)
                exc[mech] += resid; opp[mech] += nt;
                // G-run subgradient: count G neighbours flanking the central G.
                int gruns = (p == 2 ? 1 : 0) + (n == 2 ? 1 : 0);
                int run = (gruns == 2) ? CL_GGG : (gruns == 1) ? CL_GG : CL_G;
                exc[run] += resid; opp[run] += nt;
            }
        }
    }
    // Finalize each class as excess = sum_resid/sum_opp with binomial SE; gate on a minimum
    // eligible-opportunity count so sparse classes stay NaN (-> JSON null), never 0.
    static constexpr double kMinOpp = 50.0;
    auto finalize_cl = [&](int c, OxogContextExcess& o) {
        o.n = opp[c];
        if (opp[c] < kMinOpp) return;
        double rate = exc[c] / opp[c];
        o.excess = rate;
        o.se = std::sqrt(std::max(0.0, rate * (1.0 - rate)) / opp[c]);
    };
    finalize_cl(CL_OH, r.oh_radical);
    finalize_cl(CL_M1G, r.m1g);
    finalize_cl(CL_HOLE, r.hole);
    finalize_cl(CL_GGG, r.grun_ggg);
    finalize_cl(CL_GG,  r.grun_gg);
    finalize_cl(CL_G,   r.grun_g);
    if (std::isfinite(r.grun_ggg.excess) && std::isfinite(r.grun_g.excess)) {
        r.grun_gradient = r.grun_ggg.excess - r.grun_g.excess;
        double sg = std::sqrt(r.grun_ggg.se * r.grun_ggg.se + r.grun_g.se * r.grun_g.se);
        if (sg > 1e-30) r.grun_gradient_z = r.grun_gradient / sg;
    }
    double dot = 0.0, nv = 0.0, nr = 0.0;
    for (int i = 0; i < 16; ++i) {
        dot += v[i] * OXOG_REF[i];
        nv  += v[i] * v[i];
        nr  += OXOG_REF[i] * OXOG_REF[i];
    }
    if (nv > 1e-30 && nr > 0)
        r.cosine = dot / (std::sqrt(nv) * std::sqrt(nr));

    // RC log-ratio estimator: IVW mean of log(rc(XGY)/XGY) across 16 XGY contexts.
    // Uses interior trinucs only (5'+3' summed) to avoid 5' deamination confound.
    // gt_asymmetry > 0 means XGY depleted relative to complement → net G→T oxidative damage.
    {
        static constexpr double ALPHA = 0.5;
        double sum_w = 0.0, sum_we = 0.0;
        for (int p = 0; p < 4; ++p) {
            for (int n = 0; n < 4; ++n) {
                double g_i = ALPHA, c_i = ALPHA;
                for (const auto* tri : {&dp.tri_5prime_interior, &dp.tri_3prime_interior}) {
                    g_i += static_cast<double>((*tri)[p*16 + 2*4 + n]);            // XGn
                    c_i += static_cast<double>((*tri)[RC4[n]*16 + 1*4 + RC4[p]]); // rc(XGn) = rc(n)Crc(p)
                }
                double w_i  = (g_i * c_i) / (g_i + c_i);
                double eta_i = std::log(c_i / g_i);
                sum_w  += w_i;
                sum_we += w_i * eta_i;
            }
        }
        if (sum_w > 0.0) {
            r.gt_asymmetry = sum_we / sum_w;
            r.gt_rate   = std::max(0.0, 1.0 - std::exp(-r.gt_asymmetry));
        }
    }

    return r;
}

// F3: resolve the GG-breakpoint scission contrast by purine dinucleotide context.
// Trinuc idx = prev*16 + mid*4 + next (A=0,C=1,G=2,T=3). The breakpoint context is the
// (mid,next) dinucleotide; the depurination-neutral comparator is mid=A (matches the legacy
// compute_oxscission_delta A-context subtraction). Each class is the composition-corrected
// terminal-vs-interior double difference, gated on a minimum terminal context count.
ScissionContextResult compute_scission_by_context(const std::array<uint64_t, 64>& term,
                                                  const std::array<uint64_t, 64>& inter) {
    auto den = [](const std::array<uint64_t, 64>& s) {
        uint64_t d = 0;
        for (int i = 0; i < 64; ++i) d += s[i];
        return static_cast<double>(d);
    };
    auto ctx_count = [](const std::array<uint64_t, 64>& s, int mid, int next) {
        uint64_t c = 0;
        for (int prev = 0; prev < 4; ++prev) c += s[prev*16 + mid*4 + next];
        return static_cast<double>(c);
    };
    auto mid_count = [](const std::array<uint64_t, 64>& s, int mid) {
        uint64_t c = 0;
        for (int prev = 0; prev < 4; ++prev)
            for (int next = 0; next < 4; ++next) c += s[prev*16 + mid*4 + next];
        return static_cast<double>(c);
    };

    const double dt = den(term), di = den(inter);
    ScissionContextResult r;
    if (dt < 1.0 || di < 1.0) return r;   // no eligible termini -> all classes stay null

    // Depurination-neutral baseline contrast: mid=A terminal-vs-interior fraction shift.
    const double base = mid_count(term, /*A*/0) / dt - mid_count(inter, 0) / di;

    static constexpr double kMinTermCtx = 50.0;
    auto finalize_ctx = [&](int mid, int next, OxogContextExcess& o) {
        const double kt = ctx_count(term, mid, next);
        o.n = kt;
        if (kt < kMinTermCtx) return;     // sparse class -> NaN -> JSON null
        const double pt = kt / dt;
        const double pi = ctx_count(inter, mid, next) / di;
        o.excess = (pt - pi) - base;
        o.se = std::sqrt(std::max(0.0, pt * (1.0 - pt)) / dt);  // dominant terminal-fraction term
    };
    finalize_ctx(0, 2, r.ag);   // AG
    finalize_ctx(2, 2, r.gg);   // GG (legacy GG-breakpoint contrast)
    finalize_ctx(0, 0, r.aa);   // AA
    finalize_ctx(2, 0, r.ga);   // GA
    return r;
}

// F4: provisional oxidative-vs-hydrolytic C->T pathway split. Reuses the oxo_two_marker bins
// (no new accumulator). `total` is the composition-corrected interior C->T excess: the
// harmonic-mean-weighted (effective-N) IVW mean of the per-cell strand contrast
// D = T/(T+G) - A/(A+C), the same per-cell quantity compute_oxo_two_marker regresses; its
// aggregate is the C->T-driven top-strand T enrichment net of the complement baseline. The
// oxidative component is the part of `total` that co-moves with the reference-free interior
// G->T oxidation rate (otr.gt_rate), gated on the oxidation coupling being marker-consistent
// (otm.markers_consistent) and significant (|beta1_z| via otm). Because both routes give the
// identical C->T substitution, this attribution is NOT identifiable reference-free: emitted
// provisional (provisional=true, separable=false), never as a hard decomposition.
CtPathwaySplit compute_ct_pathway_split(const SampleDamageProfile& dp,
                                        const OxoTwoMarkerResult& otm,
                                        const OxogTrinucResult& otr,
                                        bool is_ss) {
    using Bins = SampleDamageProfile::OxoTwoMarkerBins;
    CtPathwaySplit r;
    const auto& panel = is_ss ? dp.oxo_two_marker_ss : dp.oxo_two_marker;

    // IVW aggregate of the per-cell strand contrast D, weighted by harmonic-mean effective N
    // (matching compute_oxo_two_marker's WLS weight). var(D) per cell from the two binomial
    // terms; SE of the weighted mean = sqrt(1/Σw_var) with w_var = 1/var(D).
    double sum_wD = 0.0, sum_w = 0.0, sum_wvar = 0.0;
    int n_cells = 0;
    for (int i = 0; i < Bins::TOTAL; ++i) {
        const auto& c = panel.cells[i];
        if (c.sum_nGT < 200 || c.sum_nAC < 200) continue;  // sparse cell -> excluded, not 0
        const double pT = static_cast<double>(c.sum_T) / c.sum_nGT;
        const double pA = static_cast<double>(c.sum_A) / c.sum_nAC;
        const double D  = pT - pA;
        const double w  = 2.0 * c.sum_nGT * c.sum_nAC
                        / static_cast<double>(c.sum_nGT + c.sum_nAC);
        const double varD = pT * (1.0 - pT) / c.sum_nGT + pA * (1.0 - pA) / c.sum_nAC;
        sum_wD += w * D;
        sum_w  += w;
        if (varD > 1e-30) sum_wvar += 1.0 / varD;
        ++n_cells;
    }
    r.n_cells = n_cells;
    if (n_cells < 6 || sum_w <= 0.0) return r;  // under-determined -> all NaN (JSON null)

    r.total    = sum_wD / sum_w;
    r.total_se = (sum_wvar > 1e-30) ? std::sqrt(1.0 / sum_wvar)
                                    : std::numeric_limits<double>::quiet_NaN();
    r.valid    = true;

    // Oxidative attribution: the C->T component co-moving with the interior G->T oxidation
    // rate, gated on a marker-consistent, significant oxidation coupling. Withheld (NaN) when
    // the coupling is not validated -> JSON null, so hydrolytic also stays null (we will not
    // claim "all hydrolytic" merely because oxidation coupling was unmeasurable).
    const double gt = otr.gt_rate;  // reference-free interior 8-oxoG rate
    const bool coupling_ok = otm.valid && otm.markers_consistent
                             && std::isfinite(otm.beta1) && std::isfinite(otm.beta1_se)
                             && otm.beta1_se > 0.0;
    if (coupling_ok) r.coupling_z = otm.beta1 / otm.beta1_se;
    if (coupling_ok && std::isfinite(gt) && gt > 0.0 && std::isfinite(r.total)) {
        // Oxidation-coupled C->T excess = coupling slope x oxidation rate, clamped to [0,total].
        double ox = otm.beta1 * gt;
        ox = std::max(0.0, std::min(ox, std::max(0.0, r.total)));
        r.oxidative    = ox;
        r.oxidative_se = std::sqrt(std::max(0.0, otm.beta1_se * otm.beta1_se * gt * gt));
        r.hydrolytic   = std::max(0.0, r.total - ox);
    }
    return r;
}

DepurScore compute_depur_score(const SampleDamageProfile& dp, bool is_ss) {
    DepurScore r;
    // exploratory; clamped, not a calibrated p-value (correlated reads).
    // Bernoulli z from raw read counts scales ~sqrt(N); clamp the reported
    // magnitude and floor p so the underflow-to-0 case is distinguishable from
    // a genuine zero. shift5/shift3 are the primary interpretable effect sizes.
    auto clamp_z = [](double z) { return std::max(-kDepurZCap, std::min(kDepurZCap, z)); };
    auto erfc_p  = [](double z) {
        double p = 0.5 * std::erfc(z / std::sqrt(2.0));
        return std::max(p, std::numeric_limits<double>::min());
    };
    r.z5     = clamp_z(static_cast<double>(dp.purine_z_5prime));
    r.shift5 = static_cast<double>(dp.purine_enrichment_5prime);
    r.shift3 = static_cast<double>(dp.purine_enrichment_3prime);

    double p5 = erfc_p(r.z5);
    if (!std::isfinite(r.z5)) {
        r.z = std::numeric_limits<double>::quiet_NaN();
        r.p = std::numeric_limits<double>::quiet_NaN();
    } else if (!is_ss && std::isfinite(dp.purine_z_3prime)) {
        // DS libraries can expose both molecular ends. The 3' term is the same
        // A+G-over-all-bases endpoint contrast, not the A/(A+G) deamination axis
        // and not the T/(T+C) negative control.
        r.z3      = clamp_z(static_cast<double>(dp.purine_z_3prime));
        double p3 = erfc_p(r.z3);
        r.z = std::min(r.z5, r.z3);
        r.p = std::max(p5, p3);
    } else {
        r.z = r.z5;
        r.p = p5;
    }
    return r;
}

// ── Damage mask ───────────────────────────────────────────────────────────────

DamageMask compute_damage_mask(const SampleDamageProfile& dp,
                                bool is_ss, double threshold, int min_cov) {
    DamageMask r;
    const double bg_5  = dp.fit_baseline_5prime;
    const double bg_3  = dp.fit_baseline_3prime;
    const double bg_tc = dp.baseline_t_freq /
                         (dp.baseline_t_freq + dp.baseline_c_freq + 1e-9);

    for (int p = 0; p < INTERP_N_POS; ++p) {
        double excess_5 = 0.0, excess_3 = 0.0;
        if (dp.tc_total_5prime[p] >= static_cast<double>(min_cov))
            excess_5 = dp.t_freq_5prime[p] - bg_5;
        if (is_ss) {
            if (dp.tc_total_3prime[p] >= static_cast<double>(min_cov))
                excess_3 = dp.t_freq_3prime[p] / dp.tc_total_3prime[p] - bg_tc;
        } else {
            if (dp.ag_total_3prime[p] >= static_cast<double>(min_cov))
                excess_3 = dp.a_freq_3prime[p] - bg_3;
        }
        r.pos[p] = (excess_5 > threshold) || (excess_3 > threshold);
    }
    for (int p = 0; p < INTERP_N_POS; ++p) {
        if (!r.pos[p]) continue;
        if (r.n_masked) r.masked_str += ',';
        r.masked_str += std::to_string(p);
        ++r.n_masked;
    }
    return r;
}
} // namespace taph
