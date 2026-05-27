// Unit tests for taph::compute_damage_context_profile.
//
// Self-contained: no gtest / catch dependency; pure <cassert>. Each test
// builds a minimal SampleDamageProfile, calls compute_damage_context_profile,
// and asserts on score finiteness and dominant_process.
//
// Invoked via `ctest` after cmake -DTAPH_BUILD_TESTS=ON.

#include "taph/library_interpretation.hpp"
#include "taph/profile_json.hpp"
#include "taph/sample_damage_profile.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
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

// ── Hexamer confound certification ───────────────────────────────────────────
// Tests the d3_confounded logic in profile_to_json:
//   rc_overlap_topk==0  AND  dmg_frac_3 < 0.5  AND  fit_offset_3prime >= 1
//   => confounded_outputs: ["d_max_3"]

static taph::HexEnrichment make_hex(int idx, double lfc, bool dc) {
    taph::HexEnrichment h;
    h.idx              = idx;
    h.log2fc           = lfc;
    h.damage_consistent = dc;
    return h;
}

void test_circligase_confound_certified() {
    auto dp = make_base();
    dp.d_max_5prime  = 0.072f;
    dp.d_max_3prime  = 0.197f;
    dp.fit_offset_3prime = 1;

    taph::ProfileJsonInput in;
    // RC-asymmetric ends (CircLigase signature)
    in.hex_end_asymmetry.rc_overlap_topk = 0;
    in.hex_end_asymmetry.rc_excess_jsd   = 0.596;
    in.hex_end_asymmetry.status          = "ok";
    // Top 3' hexamers: AAACCC, ATACCC, AGACCC — last base C, not damage-consistent
    in.top_hex_enriched_3prime = {
        make_hex(0, 2.45, false),  // AAACCC
        make_hex(1, 2.41, false),  // ATACCC
        make_hex(2, 2.28, false),  // AGACCC
    };

    std::ostringstream out;
    taph::profile_to_json(dp, out, in);
    const std::string json = out.str();

    EXPECT(json.find("\"confounded_outputs\": [\"d_max_3\"]") != std::string::npos);
    EXPECT(json.find("\"rc_overlap_topk\": 0") != std::string::npos);
    EXPECT(json.find("\"top_damage_consistent_fraction_3prime\": 0") != std::string::npos);
}

void test_t4_ligase_no_confound() {
    auto dp = make_base();
    dp.d_max_5prime  = 0.043f;
    dp.d_max_3prime  = 0.045f;
    dp.fit_offset_3prime = 1;

    taph::ProfileJsonInput in;
    // RC-symmetric ends (T4 ligase signature)
    in.hex_end_asymmetry.rc_overlap_topk = 4;
    in.hex_end_asymmetry.rc_excess_jsd   = 0.011;
    in.hex_end_asymmetry.status          = "ok";
    // Top 3' hexamers: mixed damage-consistent/not — fraction >= 0.5
    in.top_hex_enriched_3prime = {
        make_hex(10, 1.8, true),   // ends in A — damage-consistent
        make_hex(11, 1.7, true),
        make_hex(12, 1.6, false),
        make_hex(13, 1.5, true),
    };

    std::ostringstream out;
    taph::profile_to_json(dp, out, in);
    const std::string json = out.str();

    EXPECT(json.find("\"confounded_outputs\": []") != std::string::npos);
    EXPECT(json.find("\"rc_overlap_topk\": 4") != std::string::npos);
}

// ── Fix 1: per-position flatness ─────────────────────────────────────────────

void test_d5_flatness_blunting_suspected() {
    auto dp = make_base();
    dp.library_type   = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
    dp.d_max_5prime   = 0.003f;
    dp.d_max_3prime   = 0.12f;
    // Set flatness fields directly (normally populated by finalize_sample_profile)
    dp.d5_profile_flat         = true;
    dp.d3_profile_flat         = false;
    dp.d5_max_rate_pos0_4      = 0.002f;
    dp.d3_max_rate_pos0_4      = 0.11f;
    dp.d5_profile_var_pos0_9   = 1e-6f;
    dp.d3_profile_var_pos0_9   = 5e-4f;
    dp.d5_blunting_suspected   = true;

    taph::ProfileJsonInput in;
    in.hex_end_asymmetry.rc_overlap_topk = 0;
    in.hex_end_asymmetry.rc_excess_jsd   = 0.82;
    in.hex_end_asymmetry.status          = "high";

    std::ostringstream out;
    taph::profile_to_json(dp, out, in);
    const std::string json = out.str();

    EXPECT(json.find("\"d5_profile_flat\": true")       != std::string::npos);
    EXPECT(json.find("\"d3_profile_flat\": false")      != std::string::npos);
    EXPECT(json.find("\"d5_blunting_suspected\": true") != std::string::npos);
    // ss_end_asymmetry: d5 suppressed → recommended = d_max_3prime
    EXPECT(json.find("\"d5_suppressed\": true")                   != std::string::npos);
    EXPECT(json.find("\"recommended_estimate\": \"d_max_3prime\"") != std::string::npos);
}

// ── Fix 3: CircLigase selection bias corrects d_max_combined ─────────────────

void test_ss_selection_bias_corrects_combined() {
    // Typical Ellesmere-style SS sample: rc_overlap_topk=0, d3 >> d5, d5 not flat.
    // Expected: d_max_combined in JSON = d_max_5prime, source = "5prime_conservative_ss".
    auto dp = make_base();
    dp.library_type   = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
    dp.d_max_5prime   = 0.010f;
    dp.d_max_3prime   = 0.085f;   // d3 > d5 * 1.5  AND  d5 < 0.02
    dp.d_max_combined = 0.085f;   // as set by MAX_SS_ASYMMETRY branch
    dp.d_max_source   = SampleDamageProfile::DmaxSource::MAX_SS_ASYMMETRY;
    dp.d5_profile_flat = false;

    taph::ProfileJsonInput in;
    in.hex_end_asymmetry.rc_overlap_topk = 0;
    in.hex_end_asymmetry.rc_excess_jsd   = 0.79;
    in.hex_end_asymmetry.status          = "high";
    // damage-consistent 3' hexamers (real deamination, not pure CCC)
    in.top_hex_enriched_3prime = {
        make_hex(20, 2.1, true),
        make_hex(21, 1.9, true),
        make_hex(22, 1.7, true),
    };

    std::ostringstream out;
    taph::profile_to_json(dp, out, in);
    const std::string json = out.str();

    // Combined must be overridden to d_max_5prime (0.010), not 0.085
    EXPECT(json.find("\"d_max_combined\": 0.01") != std::string::npos);
    EXPECT(json.find("\"source\": \"5prime_conservative_ss\"") != std::string::npos);
    EXPECT(json.find("\"d3_selection_biased\": true")          != std::string::npos);
    EXPECT(json.find("\"selection_biased_outputs\": [\"d_max_3\"]") != std::string::npos);
    EXPECT(json.find("\"confounded_outputs\": []")             != std::string::npos);
    // output_effects.d_max_3 status must be "selection_biased"
    EXPECT(json.find("\"status\": \"selection_biased\"") != std::string::npos);
}

void test_ss_symmetric_no_selection_bias() {
    // Symmetric SS library: d5 ≈ d3, rc_overlap_topk=0 by construction but
    // the ratio condition (d3 > d5 * 1.5) does not hold → no override.
    auto dp = make_base();
    dp.library_type   = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
    dp.d_max_5prime   = 0.068f;
    dp.d_max_3prime   = 0.075f;
    dp.d_max_combined = 0.075f;
    dp.d_max_source   = SampleDamageProfile::DmaxSource::MAX_SS_ASYMMETRY;
    dp.d5_profile_flat = false;

    taph::ProfileJsonInput in;
    in.hex_end_asymmetry.rc_overlap_topk = 0;
    in.hex_end_asymmetry.rc_excess_jsd   = 0.49;
    in.hex_end_asymmetry.status          = "ok";

    std::ostringstream out;
    taph::profile_to_json(dp, out, in);
    const std::string json = out.str();

    // No override — combined stays at 0.075
    EXPECT(json.find("\"d_max_combined\": 0.075") != std::string::npos);
    EXPECT(json.find("\"source\": \"max_ss_asymmetry\"") != std::string::npos);
    EXPECT(json.find("\"d3_selection_biased\": false") != std::string::npos);
    EXPECT(json.find("\"selection_biased_outputs\": []") != std::string::npos);
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
    test_circligase_confound_certified();
    test_t4_ligase_no_confound();
    test_d5_flatness_blunting_suspected();
    test_ss_selection_bias_corrects_combined();
    test_ss_symmetric_no_selection_bias();

    if (failures == 0) {
        std::printf("all DamageContextProfile tests passed\n");
        return 0;
    }
    std::fprintf(stderr, "%d assertion(s) failed\n", failures);
    return 1;
}
