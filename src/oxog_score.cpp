#include "taph/oxog_score.hpp"
#include <cmath>

namespace taph {

OxScoreResult compute_ox_scores(const OxBinAcc* bins, int n_bins) {
    OxScoreResult r;

    // Kish ESS gate: require min effective sample size in both ancient and modern
    // partitions so a single low-weight read cannot drive S_oxog.
    static constexpr double MIN_ESS = 5.0;

    double M_all = 0.0;
    for (int b = 0; b < n_bins; ++b) {
        const auto& x = bins[b];
        bool anc_ok = x.B > 0 && x.BB > 0 && (x.B * x.B / x.BB) >= MIN_ESS;
        bool bg_ok  = x.D > 0 && x.DD > 0 && (x.D * x.D / x.DD) >= MIN_ESS;
        if (anc_ok && bg_ok) M_all += x.M_total;
    }

    if (M_all > 0) {
        double var_S = 0.0;
        for (int b = 0; b < n_bins; ++b) {
            const auto& x = bins[b];
            bool anc_ok = x.B > 0 && x.BB > 0 && (x.B * x.B / x.BB) >= MIN_ESS;
            bool bg_ok  = x.D > 0 && x.DD > 0 && (x.D * x.D / x.DD) >= MIN_ESS;
            if (!anc_ok || !bg_ok) continue;
            double alpha = x.M_total / M_all;
            double p_anc = x.A / x.B;
            double p_bg  = x.C / x.D;
            r.s_oxog += alpha * (p_anc - p_bg);

            // C→A co-movement: same q-weighting on the C/(C+A) channel.
            // Oxidation predicts s_oxog > 0 AND s_ca > 0 (G→T and its complement
            // both enriched in ancient reads vs modern).
            if (x.B_ca > 0 && x.D_ca > 0)
                r.s_ca += alpha * (x.A_ca / x.B_ca - x.C_ca / x.D_ca);

            // Delta-method variance for S_oxog.
            double v_anc = (x.AA - 2*p_anc*x.AB + p_anc*p_anc*x.BB) / (x.B*x.B);
            double v_bg  = (x.CC - 2*p_bg *x.CD + p_bg *p_bg *x.DD) / (x.D*x.D);
            double cov   = (x.AC - p_bg*x.AD - p_anc*x.BC + p_anc*p_bg*x.BD) / (x.B*x.D);
            double var_b = v_anc + v_bg - 2*cov;
            if (var_b < 0) var_b = 0;
            var_S += alpha * alpha * var_b;
        }
        r.se_s_oxog = std::sqrt(var_S);
        r.has_score = true;
    }

    // D_oriented: log-odds Chargaff contrast log(T·A / G·C) across GC bins.
    // Inverse-variance weight 1/(1/Tp + 1/Gp + 1/Cp + 1/Ap) (Haldane +0.5).
    double w_dor = 0.0;
    for (int b = 0; b < n_bins; ++b) {
        const auto& x = bins[b];
        double G  = x.TG_all - x.T_all;
        double A  = x.CA_all - x.C_all;
        double Tp = x.T_all  + 0.5, Gp = G + 0.5;
        double Cp = x.C_all  + 0.5, Ap = A + 0.5;
        double log_contrast = std::log(Tp * Ap / (Gp * Cp));
        double w = 1.0 / (1.0/Tp + 1.0/Gp + 1.0/Cp + 1.0/Ap);
        r.d_oriented += w * log_contrast;
        w_dor += w;
    }
    if (w_dor > 0) r.d_oriented /= w_dor;

    return r;
}

} // namespace taph
