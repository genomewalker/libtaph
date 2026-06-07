#include "taph/ancient_fraction.hpp"
#include "taph/fraction_fitting.hpp"
#include <algorithm>
#include <cmath>

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

    // Pos-0 peak estimate for modern
    if (mod_tc5[0] >= 50) {
        dp.modern_fraction_d5 = static_cast<float>(
            std::max(0.0, static_cast<double>(mod_t5[0]) / mod_tc5[0] - bg5));
        dp.modern_fraction_d5_computed = true;
    }
    if (mod_n3[0] >= 50) {
        dp.modern_fraction_d3 = static_cast<float>(
            std::max(0.0, static_cast<double>(mod_h3[0]) / mod_n3[0] - bg3));
        dp.modern_fraction_d3_computed = true;
    }

    double anc_bg5 = pool_interior_bg(anc_t5, anc_tc5, NP, bg5);
    double anc_bg3 = pool_interior_bg(anc_h3, anc_n3,  NP, bg3);
    double mod_bg5 = pool_interior_bg(mod_t5, mod_tc5, NP, bg5);
    double mod_bg3 = pool_interior_bg(mod_h3, mod_n3,  NP, bg3);

    bool p0a5      = position_0_artifact_5prime;
    bool bulk_p0a3 = position_0_artifact_3prime;

    // For the 3' fraction fit, only skip pos0 when the fraction data itself
    // shows depletion — the bulk artifact fires on modern-read adapter blunting
    // but the ancient fraction may have a genuine peak at pos0.
    auto frac_p0a3 = [&](const int64_t* h, const int64_t* n, double frac_bg) -> bool {
        if (!bulk_p0a3) return false;
        if (n[0] < 10)  return bulk_p0a3;
        return static_cast<double>(h[0]) / n[0] < frac_bg;
    };

    // Use modern interior (genomic T proxy) as the ancient fraction fit background.
    // This makes d_max_5prime_fit = terminal deamination above genomic T rather than
    // above the ancient component's own (elevated) interior, which previously
    // confounded damage with interior deamination signal.
    auto [ad5, al5] = fit_exp_decay_irls(anc_t5, anc_tc5, NP, mod_bg5, p0a5);
    auto [ad3, al3] = fit_exp_decay_irls(anc_h3, anc_n3,  NP, anc_bg3,
                                          frac_p0a3(anc_h3, anc_n3, anc_bg3));
    auto [md5, ml5] = fit_exp_decay_irls(mod_t5, mod_tc5, NP, mod_bg5, p0a5);
    auto [md3, ml3] = fit_exp_decay_irls(mod_h3, mod_n3,  NP, mod_bg3,
                                          frac_p0a3(mod_h3, mod_n3, mod_bg3));

    dp.damaged_fraction_d5_fit  = static_cast<float>(ad5);
    dp.damaged_fraction_lambda5 = static_cast<float>(al5);
    dp.damaged_fraction_d3_fit  = static_cast<float>(ad3);
    dp.damaged_fraction_lambda3 = static_cast<float>(al3);
    dp.modern_fraction_d5_fit   = static_cast<float>(md5);
    dp.modern_fraction_lambda5  = static_cast<float>(ml5);
    dp.modern_fraction_d3_fit   = static_cast<float>(md3);
    dp.modern_fraction_lambda3  = static_cast<float>(ml3);

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
    dp.modern_fraction_leakage_5prime =
        dp.damaged_fraction_d5_fit > MIN_ANC_SIGNAL &&
        dp.modern_fraction_d5_fit >= LEAKAGE_THRESH * dp.damaged_fraction_d5_fit;
    dp.modern_fraction_leakage_3prime =
        dp.damaged_fraction_d3_fit > MIN_ANC_SIGNAL &&
        dp.modern_fraction_d3_fit >= LEAKAGE_THRESH * dp.damaged_fraction_d3_fit;

    if (dp.modern_fraction_leakage_5prime || dp.modern_fraction_leakage_3prime)
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
            dp.modern_fraction_rate5[p] = static_cast<float>(
                static_cast<double>(mod_t5[p]) / mod_tc5[p]);
        if (mod_n3[p] >= 10)
            dp.modern_fraction_rate3[p] = static_cast<float>(
                static_cast<double>(mod_h3[p]) / mod_n3[p]);
    }

    result.valid          = dp.damaged_fraction_valid;
    result.leakage_5prime = dp.modern_fraction_leakage_5prime;
    result.leakage_3prime = dp.modern_fraction_leakage_3prime;
    result.d_anc5 = dp.damaged_fraction_d5_fit;
    result.d_anc3 = dp.damaged_fraction_d3_fit;
    result.d_mod5 = dp.modern_fraction_d5_fit;
    result.d_mod3 = dp.modern_fraction_d3_fit;
    return result;
}

} // namespace taph
