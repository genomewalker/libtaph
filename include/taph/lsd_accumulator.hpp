#pragma once

// Shared per-read length-stratified-damage (LSD) accumulation and ancient/modern
// LLR classification. Moved here from fqdup so DART (and any other libtaph
// consumer) gets IDENTICAL per-read math instead of reimplementing it —
// previously fqdup was the only tool that could produce these numbers, which
// caused cross-tool mismatches (a session investigating DART's own damage
// estimates found different results than fqdup produced from the same reads,
// because DART had no access to this accumulation logic at all).
//
// This header intentionally contains ONLY pure per-read math (no FASTQ I/O,
// no threading, no fqdup-specific masking/QC concerns) so it has zero
// dependency on any calling tool. Each consumer (fqdup, dart) keeps its own
// file-reading/threading driver and calls these functions once per read.

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "taph/length_bin_damage_profile.hpp"

namespace taph {

inline constexpr int LSD_L_MIN             = 30;
inline constexpr int LSD_L_MAX             = 500;
inline constexpr int LSD_HIST_BINS         = 128;
inline constexpr int LSD_MIN_READS_PER_BIN = 100;

// Per-bin ancient/modern LLR-classified accumulator. Field layout matches the
// legacy fqdup LlrBinAccum so prebuilt data can be moved without copying.
struct LsdLlrBinAccum {
    int64_t n_damaged   = 0;
    int64_t n_undamaged = 0;
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> tc_5_anc{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> h_3_anc{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> n_3_anc{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc_cpg{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> tc_5_anc_cpg{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc_noncpg{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> tc_5_anc_noncpg{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc_g{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> tg_5_anc{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_mod_g{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> tg_5_mod{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_mod{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> tc_5_mod{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> h_3_mod{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> n_3_mod{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> a_5_anc_all{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> c_5_anc_all{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> g_5_anc_all{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc_all{};
    // Soft-EM posterior-weighted accumulators for ancient-fraction d_max.
    static constexpr int N_SOFT_POS = 4;
    double sw_t5_anc[N_SOFT_POS]  = {};
    double sw_tc5_anc[N_SOFT_POS] = {};
    double sw_h3_anc[N_SOFT_POS]  = {};
    double sw_n3_anc[N_SOFT_POS]  = {};
    double sw_sum = 0.0;
};

// Parameters for per-read ancient/modern classification (LLR scorer).
struct LsdClassifyParams {
    double d_anc        = 0.0;
    double lam_5        = 0.3;
    double lam_3        = 0.3;
    double bg_5         = 0.5;
    double bg_3         = 0.5;
    bool   is_ss           = false;
    bool   skip_pos0_5prime = false;  // skip pos 0 in 5' LLR (hexamer artifact)
    double cpg_scale    = 1.0;
    double noncpg_scale = 1.0;
    double d_anc_contract     = -1.0;
    bool   contract_gated_off = false;
};

// Minimal, tool-independent view of a finalized bulk damage fit — exactly the
// scalar fields make_lsd_classify_params needs, so this header has no
// dependency on any consumer's own bulk-profile type (fqdup's DamageProfile,
// or DART's equivalent).
struct BulkDamageForClassify {
    bool   mixture_converged   = false;
    double mixture_d_damaged   = 0.0;
    double d_max_5prime        = 0.0;
    double d_max_3prime        = 0.0;
    double lambda_5prime       = 0.5;
    double lambda_3prime       = 0.5;
    double bg_5_tc             = 0.5;
    double bg_3_channel        = 0.5;
    bool   ss_mode             = false;
    double d_cpg_5prime        = -1.0;
    double d_noncpg_5prime     = -1.0;
    bool   pi_detected         = false;
};

// Build classify params from a finalized bulk damage fit.
LsdClassifyParams make_lsd_classify_params(const BulkDamageForClassify& bulk);

// Return raw LLR score for one read (positive = ancient-consistent).
double lsd_llr_score(const std::string& seq, const LsdClassifyParams& p);

// Classify one read as ancient (true) or modern (false) using the LLR.
bool lsd_classify_read(const std::string& seq, const LsdClassifyParams& p);

// Accumulate terminal stats for one read into an LsdLlrBinAccum.
void lsd_accumulate(const std::string& seq, LsdLlrBinAccum& acc,
                    bool ancient, bool is_ss);

// Soft-EM accumulation: add posterior-weighted pos-0 counts to sw_* fields.
void lsd_accumulate_soft(const std::string& seq, LsdLlrBinAccum& acc,
                         double w, bool is_ss);

// Length-bin edge selection mode (auto GMM / quantile / explicit edges).
struct LengthBinOptions {
    enum class Mode { DISABLED, AUTO, QUANTILE, EXPLICIT };
    Mode mode = Mode::DISABLED;
    int quantile_bins = 1;
    std::vector<int> explicit_edges;
    bool enabled() const { return mode != Mode::DISABLED; }
};

std::vector<int> compute_lsd_edges(const std::vector<uint64_t>& hist,
                                   const LengthBinOptions& options);

} // namespace taph
