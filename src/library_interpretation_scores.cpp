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
                v[p*4 + n] += std::max(0.0, k - nt * theta);
                ++r.n_ctx;
            }
        }
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
    r.z5     = clamp_z(static_cast<double>(dp.ctrl_z_5prime));
    r.shift5 = static_cast<double>(dp.purine_enrichment_5prime);
    r.shift3 = static_cast<double>(dp.purine_enrichment_3prime);

    double p5 = erfc_p(r.z5);
    if (!is_ss) {
        // BEHAVIORAL CHANGE (C5): use terminal_z_3prime (A/(A+G) purine channel)
        // not ctrl_z_3prime (T/(T+C) pyrimidine negative-control) for the DS
        // 3' depurination signal; the conjunction must mix two purine ratios.
        r.z3      = clamp_z(static_cast<double>(dp.terminal_z_3prime));
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
