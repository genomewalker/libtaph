#pragma once
#include "taph/sample_damage_profile.hpp"
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace taph {

// ── Hexamer utilities ─────────────────────────────────────────────────────────

// Decode a 12-bit hexamer code to a null-terminated 7-char array (ACGT alphabet).
std::array<char,7> decode_hex(int code);

// Encode the 6-mer at position pos in s to a 12-bit code, or -1 on ambiguous base.
int encode_hex_at(const std::string& s, int pos);

struct HexEnrichment {
    int    idx;
    double log2fc;
    bool   damage_consistent;  // first base == 'T' (C→T deamination product)
};

// Compute terminal-vs-interior log2FC for all 6-mers; return sorted descending.
// lfc_threshold: minimum log2fc to include (default 1.5 = 3× enrichment).
std::vector<HexEnrichment> compute_hex_enriched_5prime(
    const SampleDamageProfile& dp,
    float lfc_threshold = 1.5f);

struct AdapterStubs {
    std::vector<std::string> stubs5;           // up to 5 enriched 5' stubs
    std::vector<std::string> stubs3;           // up to 5 enriched 3' stubs
    bool adapter_clipped   = false;            // 5' stubs present
    bool adapter3_clipped  = false;            // 3' stubs present
    bool flag_hex_artifact = false;            // non-T-leading 5' enrichment
    std::vector<HexEnrichment> top_enriched;   // top hexamers from dp (post-clip)
};

// Detect 5' adapter stubs from dp.hexamer_count_5prime vs interior,
// and 3' stubs from hex3_terminal (4096 counts from pre-scan last-6bp)
// vs dp.hexamer_count_interior. Both capped at 5 stubs.
// 5' threshold: lfc > 3.0, first base != 'T'.
// 3' threshold: lfc > 3.0, last base != 'A' (avoids genuine G→A termini).
// 3' detection gated on adapter_clipped (5' stubs must exist first).
AdapterStubs detect_adapter_stubs(
    const SampleDamageProfile& dp,
    const uint32_t hex3_terminal[4096],
    uint64_t n_hex3);

struct HexStats {
    double entropy_terminal = 0.0;   // Shannon entropy of 5' hexamer dist (bits)
    double entropy_interior = 0.0;
    double jsd              = 0.0;   // Jensen-Shannon divergence (terminal vs interior)
    double shift_g          = 0.0;   // G-test statistic
    double shift_z          = 0.0;   // z-score of G-test vs chi2(df) null
    double shift_p          = 1.0;
};

// Compute hexamer entropy, JSD, and multinomial G-test (terminal vs interior).
HexStats compute_hex_stats(const SampleDamageProfile& dp);

// ── Score functions (pure functions of finalized SampleDamageProfile) ─────────

struct CpgScore { double z = 0.0; double p = 1.0; };

// log2(CpG-like / non-CpG-like d_max) / SE(log2 ratio).
CpgScore compute_cpg_score(const SampleDamageProfile& dp);

struct OxogInteriorScore { double z = 0.0; double p = 1.0; };

// Chargaff-symmetric G→T excess test over 16 trinucleotide contexts.
OxogInteriorScore compute_oxog_interior_score(const SampleDamageProfile& dp);

struct OxogTrinucResult {
    double cosine = std::numeric_limits<double>::quiet_NaN();
    int    n_ctx  = 0;
};

// Cosine similarity of per-context G→T residuals to empirical 8-oxoG reference.
OxogTrinucResult compute_oxog_trinuc(const SampleDamageProfile& dp);

struct DepurScore {
    double z      = 0.0;
    double p      = 1.0;
    double z5     = 0.0;
    double z3     = std::numeric_limits<double>::quiet_NaN();
    double shift5 = 0.0;
    double shift3 = 0.0;
};

// Conjunction test on 5' A/(A+G) and (DS only) 3' T/(T+C) terminal enrichment.
DepurScore compute_depur_score(const SampleDamageProfile& dp, bool is_ss);

// ── Damage mask ───────────────────────────────────────────────────────────────

static constexpr int INTERP_N_POS = 15;

struct DamageMask {
    bool        pos[INTERP_N_POS] = {};
    int         n_masked = 0;
    std::string masked_str;   // comma-separated 0-based position indices
};

// threshold: excess above background to call a position damaged (default 0.05).
// min_cov: minimum count to use a position (default 100).
DamageMask compute_damage_mask(const SampleDamageProfile& dp,
                                bool is_ss,
                                double threshold = 0.05,
                                int min_cov = 100);

// ── Library QC flags ──────────────────────────────────────────────────────────

struct LibraryQcFlags {
    bool adapter_remnant_5prime    = false;
    bool adapter_remnant_3prime    = false;
    bool hexamer_composition_bias  = false;
    bool hexamer_terminal_shift    = false;
    bool short_read_spike          = false;
    bool depurination              = false;
    bool ds_3prime_signal_absent   = false;
    bool hexamer_artifact_bias     = false;
    bool ga3_inward_displaced      = false;
};

LibraryQcFlags compute_library_qc_flags(
    const SampleDamageProfile& dp,
    bool is_ss,
    bool flag_hex_artifact,
    double jsd,
    double h_term,
    double short_read_frac);

// ── Preservation / authenticity ───────────────────────────────────────────────

struct PreservationSummary {
    double      authenticity_eff      = 0.0;
    double      authenticity_evidence = 0.0;
    double      d5_raw                = 0.0;
    double      d5_hexamer_corrected  = 0.0;
    bool        d5_was_corrected      = false;
    double      oxidation_eff         = 0.0;
    double      oxidation_evidence    = 0.0;
    double      qc_risk_eff           = 0.0;
    double      qc_evidence           = 0.0;
    const char* label                 = "modern-like";
};

PreservationSummary compute_preservation_summary(
    const SampleDamageProfile& dp,
    bool   is_ss,
    bool   adapter_clipped,
    bool   flag_hex_artifact,
    double cpg_score_z,
    double oxog_score_z,
    double oxog_trinuc_cosine,
    double hex_shift_p);

} // namespace taph
