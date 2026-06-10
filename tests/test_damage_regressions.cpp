#include "taph/ancient_fraction.hpp"
#include "taph/oxog_score.hpp"
#include "taph/sample_damage_profile.hpp"
#include "../src/damage_estimation_finalize_ctx.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace {

int64_t count_for_rate(double rate, int64_t n) {
    return static_cast<int64_t>(std::llround(rate * static_cast<double>(n)));
}

void require(bool cond, const char* msg) {
    if (!cond) {
        std::fprintf(stderr, "FAIL: %s\n", msg);
        std::abort();
    }
}

void test_oxog_ca_uses_product_base() {
    taph::OxBinAcc bin{};

    for (int i = 0; i < 20; ++i) {
        taph::update_ox_bins(bin, 0.9, 30, 70, 70, 30);
        taph::update_ox_bins(bin, 0.1,  5, 95, 95,  5);
    }

    auto r = taph::compute_ox_scores(&bin, 1);
    require(r.has_score, "OxoG synthetic case should pass ESS gate");
    require(r.s_oxog > 0.0, "G->T score should be ancient-enriched");
    require(r.s_ca > 0.0, "C->A score should use A/(C+A), not C/(C+A)");
}

void test_ancient_fraction_3prime_uses_modern_background() {
    taph::AncientFractionBins bin{};
    bin.n_damaged = 8000;
    bin.n_undamaged = 4000;
    bin.sw_sum = 100.0;

    constexpr int64_t n = 10000;
    for (int p = 0; p < taph::TAPH_FRAC_N_POS; ++p) {
        const double decay = std::exp(-0.30 * static_cast<double>(p));
        const double anc5_rate = (p >= 1 && p <= 6) ? 0.10 + 0.15 * decay : 0.16;
        const double anc3_rate = (p >= 1 && p <= 6) ? 0.10 + 0.20 * decay : 0.23;

        bin.tc_5_anc[p] = n;
        bin.t_5_anc[p]  = count_for_rate(anc5_rate, n);
        bin.tc_5_mod[p] = n;
        bin.t_5_mod[p]  = count_for_rate(0.10, n);

        bin.n_3_anc[p] = n;
        bin.h_3_anc[p] = count_for_rate(anc3_rate, n);
        bin.n_3_mod[p] = n;
        bin.h_3_mod[p] = count_for_rate(0.10, n);
    }

    taph::SampleDamageProfile dp;
    auto r = taph::compute_ancient_fraction(
        &bin, 1, 0.10, 0.10, false, false, dp);

    require(std::isfinite(r.d_anc3), "ancient 3prime fit should be finite");
    require(r.d_anc3 > 0.15f, "ancient 3prime fit should use modern background");
}

void test_gc_bin_llr_keeps_bin_specific_baseline() {
    taph::SampleDamageProfile dp;
    taph::FinalCtx ctx;
    ctx.baseline_tc = 0.50;
    ctx.baseline_ag = 0.50;
    dp.lambda_5prime = 0.30f;
    dp.lambda_3prime = 0.30f;

    auto& b = dp.gc_bins[2];
    b.n_reads = 2000;
    b.t_interior = 2000;
    b.c_interior = 20000;
    b.a_interior = 10000;
    b.g_interior = 10000;

    for (int p = 0; p < 15; ++p) {
        const double decay = std::exp(-0.30 * static_cast<double>(p));
        const double rate = 0.09 + 0.20 * decay;
        const uint64_t n_tc = 10000;
        b.t_counts[p] = static_cast<uint64_t>(std::llround(rate * n_tc));
        b.c_counts[p] = n_tc - b.t_counts[p];
        b.a_counts[p] = 5000;
        b.g_counts[p] = 5000;
    }

    taph::finalize_dmax(dp, ctx);

    require(dp.gc_bins[2].valid, "GC bin should be valid");
    require(dp.gc_bins[2].baseline_tc < 0.11f, "GC-bin baseline should not be overwritten by global baseline");
    require(dp.gc_bins[2].baseline_tc > 0.08f, "GC-bin baseline should retain bin-specific value");
}

} // namespace

int main() {
    test_oxog_ca_uses_product_base();
    test_ancient_fraction_3prime_uses_modern_background();
    test_gc_bin_llr_keeps_bin_specific_baseline();
    std::printf("damage_regressions: PASS\n");
    return 0;
}
