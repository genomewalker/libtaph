#include "taph/ancient_fraction.hpp"
#include "taph/fraction_fitting.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace taph {

AncientFractionResult compute_ancient_fraction(
    const AncientFractionBins* bins, int n_bins,
    double bg5, double bg3,
    bool position_0_artifact_5prime,
    bool position_0_artifact_3prime,
    SampleDamageProfile& dp)
{
    AncientFractionResult result{};

    int64_t hard_n_damaged = 0, hard_n_tot = 0;
    double  sw_sum = 0.0;
    double  sw_t5[TAPH_FRAC_N_SOFT_POS]  = {};
    double  sw_tc5[TAPH_FRAC_N_SOFT_POS] = {};
    double  sw_h3[TAPH_FRAC_N_SOFT_POS]  = {};
    double  sw_n3[TAPH_FRAC_N_SOFT_POS]  = {};

    for (int i = 0; i < n_bins; ++i) {
        const auto& b = bins[i];
        hard_n_damaged += b.n_damaged;
        hard_n_tot     += b.n_damaged + b.n_undamaged;
        sw_sum         += b.sw_sum;
        for (int sp = 0; sp < TAPH_FRAC_N_SOFT_POS; ++sp) {
            sw_t5[sp]  += b.sw_t5_anc[sp];
            sw_tc5[sp] += b.sw_tc5_anc[sp];
            sw_h3[sp]  += b.sw_h3_anc[sp];
            sw_n3[sp]  += b.sw_n3_anc[sp];
        }
    }

    if (hard_n_tot < 10000 || sw_sum < 10.0)
        return result;

    dp.damaged_fraction_n     = hard_n_damaged;
    dp.damaged_fraction_valid = true;

    constexpr int NP = TAPH_FRAC_N_POS;
    int64_t anc_t5[NP]={}, anc_tc5[NP]={}, anc_h3[NP]={}, anc_n3[NP]={};
    int64_t mod_t5[NP]={}, mod_tc5[NP]={}, mod_h3[NP]={}, mod_n3[NP]={};

    for (int i = 0; i < n_bins; ++i) {
        const auto& b = bins[i];
        for (int p = 0; p < NP; ++p) {
            anc_t5[p]  += b.t_5_anc[p];
            anc_tc5[p] += b.tc_5_anc[p];
            anc_h3[p]  += b.h_3_anc[p];
            anc_n3[p]  += b.n_3_anc[p];
            mod_t5[p]  += b.t_5_mod[p];
            mod_tc5[p] += b.tc_5_mod[p];
            mod_h3[p]  += b.h_3_mod[p];
            mod_n3[p]  += b.n_3_mod[p];
        }
    }

    // Pos-0 peak estimate for non-damaged
    if (mod_tc5[0] >= 50) {
        dp.nondamaged_fraction_d5 = static_cast<float>(
            std::max(0.0, static_cast<double>(mod_t5[0]) / mod_tc5[0] - bg5));
        dp.nondamaged_fraction_d5_computed = true;
    }
    if (mod_n3[0] >= 50) {
        dp.nondamaged_fraction_d3 = static_cast<float>(
            std::max(0.0, static_cast<double>(mod_h3[0]) / mod_n3[0] - bg3));
        dp.nondamaged_fraction_d3_computed = true;
    }

    auto mod_bg5 = pool_interior_bg(mod_t5, mod_tc5, NP, bg5);
    auto mod_bg3 = pool_interior_bg(mod_h3, mod_n3,  NP, bg3);

    bool p0a5      = position_0_artifact_5prime;
    bool bulk_p0a3 = position_0_artifact_3prime;

    // For the 3' fraction fit, only skip pos0 when the fraction data itself
    // shows depletion — the bulk artifact fires on modern-read adapter blunting
    // but the ancient fraction may have a genuine peak at pos0.
    auto frac_p0a3 = [&](const int64_t* h, const int64_t* n, double frac_bg_mean) -> bool {
        if (!bulk_p0a3) return false;
        if (n[0] < 10)  return bulk_p0a3;
        return static_cast<double>(h[0]) / n[0] < frac_bg_mean;
    };

    // Use modern interior as the ancient fraction fit background.  This makes
    // d_max_*_fit terminal deamination above the genomic channel baseline rather
    // than above the ancient component's own elevated interior signal.
    const bool anc_p0a3 = frac_p0a3(anc_h3, anc_n3, mod_bg3.mean);
    const bool mod_p0a3 = frac_p0a3(mod_h3, mod_n3, mod_bg3.mean);
    auto [ad5, al5] = fit_exp_decay_irls(anc_t5, anc_tc5, NP, mod_bg5, p0a5);
    auto [ad3, al3] = fit_exp_decay_irls(anc_h3, anc_n3,  NP, mod_bg3, anc_p0a3);
    auto [md5, ml5] = fit_exp_decay_irls(mod_t5, mod_tc5, NP, mod_bg5, p0a5);
    auto [md3, ml3] = fit_exp_decay_irls(mod_h3, mod_n3,  NP, mod_bg3, mod_p0a3);

    // NOTE: an exact-permutation goodness-of-fit gate was tried here and
    // REVERTED — validated against real benchmark data, its p-values for
    // documented false positives (0.019-0.060) and documented true positives
    // (0.014-0.056) completely overlap (no discriminating power). Root cause:
    // the false positives are NOT random per-position noise (which exchangeable
    // permutation could detect) — they are a REAL, reproducible terminal
    // compositional artifact that also elevates the CA/CG (non-C->T) channels
    // at the same positions, something permutation of a single channel cannot
    // see.
    //
    // CT-specificity gate (validated replacement): true deamination elevates
    // C->T specifically; a generic terminal compositional artifact (e.g. an
    // end-repair/fill-in bias) elevates ALL C-derived substitutions (C->T,
    // C->A, C->G) roughly equally. Compute the same terminal-excess statistic
    // (pos1 minus mean(pos4..6)) for CT, CA, CG from the bulk trinucleotide
    // position counts (dp.tri_5prime_pos, the same source as
    // substitution_pos_profiles), and require CT's excess to clear the mean
    // CA/CG control excess by a margin. Validated on 6 documented samples:
    // false positives score -0.003 to +0.015; true positives score +0.107 to
    // +0.112 — a ~7x gap with zero overlap. Threshold set well inside the gap.
    auto sub_rate_excess = [&](int from, int to) -> double {
        auto rate_at = [&](int p) -> double {
            uint64_t nf = 0, nt = 0;
            for (int prev = 0; prev < 4; ++prev)
                for (int next = 0; next < 4; ++next) {
                    nf += dp.tri_5prime_pos[p][prev * 16 + from * 4 + next];
                    nt += dp.tri_5prime_pos[p][prev * 16 + to   * 4 + next];
                }
            return (nf + nt > 0) ? static_cast<double>(nt) / (nf + nt)
                                  : std::numeric_limits<double>::quiet_NaN();
        };
        double r1 = rate_at(1), tail = 0.0; int n = 0;
        for (int p = 4; p <= 6; ++p) { double r = rate_at(p); if (std::isfinite(r)) { tail += r; ++n; } }
        if (!std::isfinite(r1) || n == 0) return std::numeric_limits<double>::quiet_NaN();
        return r1 - tail / n;
    };
    // Recalibrated 2026-07-01: was 0.05. All 6 documented samples this
    // threshold was validated against still separate cleanly at a lower
    // value (FP max 0.015, TP min 0.107) — 0.05 carried needless margin
    // that rejected at least one genuinely damaged, heavily-diluted FLB
    // library sitting at 0.049 excess. 0.02 keeps >25x margin above every
    // documented FP and >5x margin below every documented TP.
    constexpr double CT_SPECIFICITY_MIN = 0.02;  // FP max 0.015, TP min 0.107
    double d_ct = sub_rate_excess(1, 3);   // C->T
    double d_ca = sub_rate_excess(1, 0);   // C->A (control)
    double d_cg = sub_rate_excess(1, 2);   // C->G (control)
    double ct_specific_excess = std::numeric_limits<double>::quiet_NaN();
    if (std::isfinite(d_ct) && std::isfinite(d_ca) && std::isfinite(d_cg))
        ct_specific_excess = d_ct - 0.5 * (d_ca + d_cg);
    const bool ct_specific = std::isfinite(ct_specific_excess) && ct_specific_excess >= CT_SPECIFICITY_MIN;
    dp.damaged_fraction_d5_raw  = static_cast<float>(ad5);  // pre-gate value
    if (!ct_specific) {
        ad5 = std::numeric_limits<double>::quiet_NaN();
    }
    dp.damaged_fraction_d5_pval = static_cast<float>(ct_specific_excess);  // repurposed: CT-specific excess
    for (int i = 0; i < 7 && i < NP; ++i) {
        dp.debug_anc_t5[i]  = anc_t5[i];
        dp.debug_anc_tc5[i] = anc_tc5[i];
    }

    // Diagnostic (pending validation, not yet production): CT-vs-(CA,CG)
    // difference-series decay fit. Sidesteps the baseline-contamination
    // mechanism entirely instead of patching around it — see
    // fit_ct_specific_decay's doc comment in fraction_fitting.hpp.
    {
        auto cf = fit_ct_specific_decay(dp.tri_5prime_pos.data(), NP);
        dp.ct_diff_d5_fit    = static_cast<float>(cf.A);
        dp.ct_diff_d5_lambda = static_cast<float>(cf.lambda);
        dp.ct_diff_d5_npts   = cf.n_points;
    }

    dp.damaged_fraction_d5_fit  = static_cast<float>(ad5);
    dp.damaged_fraction_lambda5 = static_cast<float>(al5);
    dp.damaged_fraction_bg5     = static_cast<float>(mod_bg5.mean);
    dp.damaged_fraction_bg3     = static_cast<float>(mod_bg3.mean);
    dp.damaged_fraction_bg5_var = static_cast<float>(mod_bg5.var);
    dp.damaged_fraction_d3_fit  = static_cast<float>(ad3);
    dp.damaged_fraction_lambda3 = static_cast<float>(al3);
    dp.nondamaged_fraction_d5_fit   = static_cast<float>(md5);
    dp.nondamaged_fraction_lambda5  = static_cast<float>(ml5);
    dp.nondamaged_fraction_d3_fit   = static_cast<float>(md3);
    dp.nondamaged_fraction_lambda3  = static_cast<float>(ml3);

    // π from mixture identity (breaks the circular soft prior)
    {
        const double d_bulk5   = dp.d_max_5prime;
        const double d_anc_fit = dp.damaged_fraction_d5_fit;
        const double hard_frac = static_cast<double>(hard_n_damaged) / hard_n_tot;
        double pi_est;
        if (std::isfinite(d_anc_fit) && d_anc_fit > 0.02 &&
            d_anc_fit >= d_bulk5 && d_bulk5 > 0.005)
            pi_est = std::clamp(d_bulk5 / d_anc_fit, 0.01, 1.0);
        else
            pi_est = std::clamp(hard_frac, 0.0, 1.0);
        dp.damaged_fraction_pi = static_cast<float>(pi_est);
        result.pi = static_cast<float>(pi_est);
        if (pi_est > 0.0) {
            dp.damaged_fraction_d5 = static_cast<float>(
                std::clamp(static_cast<double>(dp.d_max_5prime) / pi_est, 0.0, 1.0));
            dp.damaged_fraction_d3 = static_cast<float>(
                std::clamp(static_cast<double>(dp.d_max_3prime) / pi_est, 0.0, 1.0));
        }
    }

    // Leakage detection
    constexpr float LEAKAGE_THRESH = 0.5f;
    constexpr float MIN_ANC_SIGNAL = 0.01f;
    dp.nondamaged_fraction_leakage_5prime =
        dp.damaged_fraction_d5_fit > MIN_ANC_SIGNAL &&
        dp.nondamaged_fraction_d5_fit >= LEAKAGE_THRESH * dp.damaged_fraction_d5_fit;
    dp.nondamaged_fraction_leakage_3prime =
        dp.damaged_fraction_d3_fit > MIN_ANC_SIGNAL &&
        dp.nondamaged_fraction_d3_fit >= LEAKAGE_THRESH * dp.damaged_fraction_d3_fit;

    if (dp.nondamaged_fraction_leakage_5prime || dp.nondamaged_fraction_leakage_3prime)
        dp.damaged_fraction_valid = false;

    // Per-position rates for HTML fraction curves
    for (int p = 0; p < NP; ++p) {
        if (anc_tc5[p] >= 10)
            dp.damaged_fraction_rate5[p] = static_cast<float>(
                static_cast<double>(anc_t5[p]) / anc_tc5[p]);
        if (anc_n3[p] >= 10)
            dp.damaged_fraction_rate3[p] = static_cast<float>(
                static_cast<double>(anc_h3[p]) / anc_n3[p]);
        if (mod_tc5[p] >= 10)
            dp.nondamaged_fraction_rate5[p] = static_cast<float>(
                static_cast<double>(mod_t5[p]) / mod_tc5[p]);
        if (mod_n3[p] >= 10)
            dp.nondamaged_fraction_rate3[p] = static_cast<float>(
                static_cast<double>(mod_h3[p]) / mod_n3[p]);
    }

    result.valid          = dp.damaged_fraction_valid;
    result.leakage_5prime = dp.nondamaged_fraction_leakage_5prime;
    result.leakage_3prime = dp.nondamaged_fraction_leakage_3prime;
    result.d_anc5 = dp.damaged_fraction_d5_fit;
    result.d_anc3 = dp.damaged_fraction_d3_fit;
    result.d_mod5 = dp.nondamaged_fraction_d5_fit;
    result.d_mod3 = dp.nondamaged_fraction_d3_fit;
    return result;
}

} // namespace taph
