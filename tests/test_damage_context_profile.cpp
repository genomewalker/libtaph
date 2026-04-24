// Unit tests for taph::compute_damage_context_profile.
//
// Self-contained: no gtest / catch dependency; pure <cassert>. Each test
// builds a minimal SampleDamageProfile, calls compute_damage_context_profile,
// and asserts on score finiteness and dominant_process.
//
// Invoked via `ctest` after cmake -DTAPH_BUILD_TESTS=ON.

#include "taph/library_interpretation.hpp"
#include "taph/sample_damage_profile.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

using taph::DamageContextProfile;
using taph::SampleDamageProfile;
using taph::compute_damage_context_profile;
using D = DamageContextProfile::DominantProcess;

namespace {

constexpr float NaNf = std::numeric_limits<float>::quiet_NaN();
constexpr double NaNd = std::numeric_limits<double>::quiet_NaN();
constexpr double Infd = std::numeric_limits<double>::infinity();

SampleDamageProfile make_base(size_t n_reads = 100'000) {
    SampleDamageProfile dp;
    dp.n_reads = n_reads;
    dp.d_max_5prime = 0.0f;
    dp.d_max_3prime = 0.0f;
    dp.lambda_5prime = 0.3f;
    dp.lambda_3prime = 0.3f;
    dp.log2_cpg_ratio = 0.0f;
    dp.dipyr_contrast = 0.0f;
    dp.ox_gt_asymmetry = 0.0f;
    dp.purine_enrichment_5prime = 0.0f;
    dp.position_0_artifact_5prime = false;
    dp.position_0_artifact_3prime = false;
    dp.fit_offset_5prime = 1;
    dp.fit_offset_3prime = 1;
    dp.s_oxog_16ctx.fill(0.0f);
    return dp;
}

int failures = 0;

#define EXPECT(cond) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); \
        ++failures; \
    } \
} while (0)

// ── Rule coverage ─────────────────────────────────────────────────────────────

void test_insufficient_coverage() {
    auto dp = make_base(500);
    dp.d_max_5prime = 0.2f;
    auto r = compute_damage_context_profile(dp, 0.0, 0.0, false, false, false);
    EXPECT(r.dominant_process == D::None);
    EXPECT(std::isfinite(r.terminal_deamination_score));  // score still reported
}

void test_terminal_nan_forces_none() {
    auto dp = make_base();
    dp.d_max_5prime = NaNf;
    dp.d_max_3prime = NaNf;
    auto r = compute_damage_context_profile(dp, 5.0, 0.0, false, false, false);
    EXPECT(r.dominant_process == D::None);
    EXPECT(std::isnan(r.terminal_deamination_score));
}

void test_low_damage() {
    auto dp = make_base();
    dp.d_max_5prime = 0.005f;
    dp.d_max_3prime = 0.005f;
    auto r = compute_damage_context_profile(dp, 0.0, 0.0, false, false, false);
    EXPECT(r.dominant_process == D::LowDamage);
}

void test_cytosine_deamination_fallback() {
    auto dp = make_base();
    dp.d_max_5prime = 0.20f;
    dp.d_max_3prime = 0.18f;
    auto r = compute_damage_context_profile(dp, 0.0, 0.0, false, false, false);
    EXPECT(r.dominant_process == D::CytosineDeamination);
}

void test_cpg_enriched() {
    auto dp = make_base();
    dp.d_max_5prime = 0.20f;
    dp.log2_cpg_ratio = 0.40f;
    auto r = compute_damage_context_profile(dp, 6.0, 0.0, false, false, false);
    EXPECT(r.dominant_process == D::CpgEnrichedDeamination);
}

void test_cpg_z_high_but_effect_too_small() {
    // Large-N sample: z high, log2 ratio below 0.15 floor → falls through.
    auto dp = make_base(10'000'000);
    dp.d_max_5prime = 0.20f;
    dp.log2_cpg_ratio = 0.05f;
    auto r = compute_damage_context_profile(dp, 8.0, 0.0, false, false, false);
    EXPECT(r.dominant_process == D::CytosineDeamination);
    EXPECT(std::isfinite(r.cpg_context_score) && r.cpg_context_score > 0.7f);
}

void test_oxidative_like() {
    auto dp = make_base();
    dp.d_max_5prime = 0.04f;  // td ≈ 0.33, passes lt(td, 0.5)
    dp.d_max_3prime = 0.02f;
    dp.ox_gt_asymmetry = 0.06f;
    auto r = compute_damage_context_profile(dp, 0.0, 0.0, false, false, false);
    EXPECT(r.dominant_process == D::OxidativeLike);
}

void test_dipyrimidine_biased() {
    auto dp = make_base();
    dp.d_max_5prime = 0.20f;
    dp.d_max_3prime = 0.05f;
    dp.dipyr_contrast = 0.04f;  // > 0.4 * kDipyrNorm(0.05)
    auto r = compute_damage_context_profile(dp, 0.0, 0.0, false, false, false);
    EXPECT(r.dominant_process == D::DipyrimidineBiased);
}

void test_fragmentation_bias() {
    auto dp = make_base();
    dp.d_max_5prime = 0.01f;
    dp.d_max_3prime = 0.01f;
    dp.purine_enrichment_5prime = 0.10f;  // > 0.5 * kFragNorm(0.15)
    auto r = compute_damage_context_profile(dp, 0.0, 0.0, false, false, false);
    EXPECT(r.dominant_process == D::FragmentationBias);
}

// ── library_artifact gating ───────────────────────────────────────────────────

void test_high_hex_z_alone_does_not_fire_artifact() {
    // Synthetic-clean scenario: z=20, no flags, genuine damage present.
    // Must NOT label library_artifact_likely; score is reported but label
    // requires flag corroboration.
    auto dp = make_base();
    dp.d_max_5prime = 0.20f;
    dp.d_max_3prime = 0.18f;
    auto r = compute_damage_context_profile(dp, 0.0, 20.0, false, false, false);
    EXPECT(r.dominant_process != D::LibraryArtifactLikely);
    EXPECT(r.library_artifact_score > 0.9f);
}

void test_artifact_flag_plus_weak_damage_fires() {
    auto dp = make_base();
    dp.d_max_5prime = 0.03f;
    dp.d_max_3prime = 0.03f;
    auto r = compute_damage_context_profile(dp, 0.0, 0.0, true, false, false);
    EXPECT(r.dominant_process == D::LibraryArtifactLikely);
}

void test_artifact_flag_plus_strong_damage_keeps_damage_label() {
    auto dp = make_base();
    dp.d_max_5prime = 0.25f;
    dp.d_max_3prime = 0.22f;
    auto r = compute_damage_context_profile(dp, 0.0, 0.0, true, false, false);
    EXPECT(r.dominant_process != D::LibraryArtifactLikely);
    EXPECT(r.library_artifact_score >= 1.0f);  // score still surfaces flag
    EXPECT(r.evidence.adapter_clipped == true);
}

void test_fit_offset_drives_artifact() {
    auto dp = make_base();
    dp.d_max_5prime = 0.02f;
    dp.fit_offset_5prime = 3;
    auto r = compute_damage_context_profile(dp, 0.0, 0.0, false, false, false);
    EXPECT(r.dominant_process == D::LibraryArtifactLikely);
    EXPECT(r.evidence.fit_offset_5prime == 3);
}

// ── NaN/inf score behaviour ───────────────────────────────────────────────────

void test_inf_inputs_yield_nan_scores() {
    auto dp = make_base();
    dp.d_max_5prime = 0.15f;
    dp.log2_cpg_ratio = static_cast<float>(Infd);
    dp.dipyr_contrast = NaNf;
    dp.purine_enrichment_5prime = static_cast<float>(Infd);
    auto r = compute_damage_context_profile(dp, Infd, 0.0, false, false, false);
    EXPECT(std::isnan(r.cpg_context_score));
    EXPECT(std::isnan(r.dipyrimidine_context_score));
    EXPECT(std::isnan(r.fragmentation_context_score));
    // Terminal score is finite because d_max_5prime is finite.
    EXPECT(std::isfinite(r.terminal_deamination_score));
}

void test_one_dmax_nan_other_finite() {
    auto dp = make_base();
    dp.d_max_5prime = 0.20f;
    dp.d_max_3prime = NaNf;
    auto r = compute_damage_context_profile(dp, 0.0, 0.0, false, false, false);
    EXPECT(std::isfinite(r.terminal_deamination_score));
    EXPECT(r.dominant_process == D::CytosineDeamination);
}

// ── Boundary semantics (strict > / <) ─────────────────────────────────────────

void test_exactly_threshold_does_not_fire_gt() {
    // terminal score = 1 - exp(-d_max / 0.10). Pick d_max so score == 0.10
    // exactly (td < 0.10 is strict, so this must NOT collapse to low_damage).
    auto dp = make_base();
    // 1 - exp(-d_max / 0.10) = 0.10  =>  d_max = -0.10 * ln(0.90) ≈ 0.01053605
    dp.d_max_5prime = 0.01053605f;
    dp.d_max_3prime = 0.01053605f;
    auto r = compute_damage_context_profile(dp, 0.0, 0.0, false, false, false);
    // score should be ~0.10; strict lt(td, 0.10) must fail.
    // With float rounding the branch can go either way; guard only the
    // semantics: if td >= 0.10f, result must NOT be LowDamage.
    if (r.terminal_deamination_score >= 0.10f) {
        EXPECT(r.dominant_process != D::LowDamage);
    }
}

// ── Evidence plumbing ─────────────────────────────────────────────────────────

void test_evidence_is_populated() {
    auto dp = make_base();
    dp.d_max_5prime = 0.15f;
    dp.d_max_3prime = 0.04f;
    dp.lambda_5prime = 0.25f;
    dp.lambda_3prime = 0.07f;
    dp.log2_cpg_ratio = 0.38f;
    dp.dipyr_contrast = 0.01f;
    dp.ox_gt_asymmetry = 0.02f;
    dp.purine_enrichment_5prime = 0.05f;
    dp.position_0_artifact_5prime = true;
    dp.fit_offset_3prime = 2;
    auto r = compute_damage_context_profile(dp, 3.5, 1.2, false, false, false);
    EXPECT(r.evidence.d_max_5 == 0.15f);
    EXPECT(r.evidence.d_max_3 == 0.04f);
    EXPECT(r.evidence.lambda_5 == 0.25f);
    EXPECT(r.evidence.lambda_3 == 0.07f);
    EXPECT(r.evidence.log2_cpg_ratio == 0.38f);
    EXPECT(std::abs(r.evidence.cpg_z - 3.5f) < 1e-5f);
    EXPECT(std::abs(r.evidence.hex_shift_z - 1.2f) < 1e-5f);
    EXPECT(r.evidence.position_0_artifact_5prime == true);
    EXPECT(r.evidence.fit_offset_3prime == 2);
    EXPECT(r.evidence.n_reads == 100'000);
    EXPECT(!r.dominant_process_str.empty());
}

void test_to_string_covers_all_enumerators() {
    EXPECT(std::string(taph::to_string(D::None)) == "none");
    EXPECT(std::string(taph::to_string(D::LowDamage)) == "low_damage");
    EXPECT(std::string(taph::to_string(D::CytosineDeamination)) == "cytosine_deamination");
    EXPECT(std::string(taph::to_string(D::CpgEnrichedDeamination)) == "cpg_enriched_deamination");
    EXPECT(std::string(taph::to_string(D::DipyrimidineBiased)) == "dipyrimidine_biased");
    EXPECT(std::string(taph::to_string(D::OxidativeLike)) == "oxidative_like");
    EXPECT(std::string(taph::to_string(D::FragmentationBias)) == "fragmentation_bias");
    EXPECT(std::string(taph::to_string(D::LibraryArtifactLikely)) == "library_artifact_likely");
}

} // namespace

int main() {
    test_insufficient_coverage();
    test_terminal_nan_forces_none();
    test_low_damage();
    test_cytosine_deamination_fallback();
    test_cpg_enriched();
    test_cpg_z_high_but_effect_too_small();
    test_oxidative_like();
    test_dipyrimidine_biased();
    test_fragmentation_bias();
    test_high_hex_z_alone_does_not_fire_artifact();
    test_artifact_flag_plus_weak_damage_fires();
    test_artifact_flag_plus_strong_damage_keeps_damage_label();
    test_fit_offset_drives_artifact();
    test_inf_inputs_yield_nan_scores();
    test_one_dmax_nan_other_finite();
    test_exactly_threshold_does_not_fire_gt();
    test_evidence_is_populated();
    test_to_string_covers_all_enumerators();

    if (failures == 0) {
        std::printf("all DamageContextProfile tests passed\n");
        return 0;
    }
    std::fprintf(stderr, "%d assertion(s) failed\n", failures);
    return 1;
}
