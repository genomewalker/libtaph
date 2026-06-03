#pragma once
#include "taph/sample_damage_profile.hpp"
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace taph {

// Per-length-bin damage characterization.  All fields default to 0 / NaN to
// allow safe reads before finalization.  Populated by estimate_damage_by_length().
struct LengthBinDamageProfile {
    static constexpr int N_POS    = 15;
    static constexpr int N_GC_BINS = 10;

    int     length_lo = 0;
    int     length_hi = 0;
    int64_t n_reads   = 0;

    double d_max_5prime  = 0.0;
    double d_max_3prime  = 0.0;
    double lambda_5prime = 0.0;
    double lambda_3prime = 0.0;
    double bg_5prime     = 0.0;
    double bg_3prime     = 0.0;
    // DEPRECATED (removal in next major release): was upstream-GC (GpC) contrast, not CpG.
    // Use log2_cpg_ratio / cpg_z. JSON emits "cpg_contrast": null as migration shim.
    double cpg_contrast  = std::numeric_limits<double>::quiet_NaN();
    bool   validated     = false;
    bool   ss_mode       = false;
    std::string source   = "none";

    std::array<double, N_POS> per_pos_5prime_ct{};
    std::array<double, N_POS> per_pos_3prime{};

    // GC mixture model
    double mixture_d_damaged    = 0.0;
    double mixture_d_reference  = 0.0;
    double mixture_d_population = 0.0;
    double mixture_pi_damaged   = 0.0;
    int    mixture_n_components = 0;
    bool   mixture_converged    = false;
    bool   mixture_identifiable = false;

    std::array<double,  N_GC_BINS> gc_d_max{};
    std::array<int64_t, N_GC_BINS> gc_n_reads{};
    std::array<double,  N_GC_BINS> gc_p_damaged{};

    // Per-read LLR unmixing
    int64_t n_damaged   = 0;
    int64_t n_undamaged = 0;
    double  d_max_5_damaged          = -1.0;
    double  d_max_3_damaged          = -1.0;
    double  d_max_5_cpg_damaged      = -1.0;
    double  d_max_5_noncpg_damaged   = -1.0;
    double  s_gt_5_damaged_vs_undamaged = std::numeric_limits<double>::quiet_NaN();
    double  g_to_t_5_damaged         = std::numeric_limits<double>::quiet_NaN();
    double  pG_terminal_5_damaged    = std::numeric_limits<double>::quiet_NaN();
    double  pG_interior_5_damaged    = std::numeric_limits<double>::quiet_NaN();

    std::array<double, N_POS> per_pos_5prime_ct_damaged{};
    std::array<double, N_POS> per_pos_3prime_damaged{};
    std::array<double, N_POS> per_pos_5prime_ct_cpg_damaged{};
    std::array<double, N_POS> per_pos_5prime_ct_noncpg_damaged{};
    std::array<double, N_POS> per_pos_5prime_gt_damaged{};
    std::array<double, N_POS> per_pos_5prime_gt_undamaged{};

    // Trinucleotide spectra
    std::array<int64_t, 64> tri_5prime_terminal{};
    std::array<int64_t, 64> tri_5prime_interior{};
    std::array<int64_t, 64> tri_3prime_terminal{};
    std::array<int64_t, 64> tri_3prime_interior{};

    // OxoG channels per bin
    float  ox_stop_rate_baseline  = 0.0f;
    float  ox_stop_rate_terminal  = 0.0f;
    float  ox_uniformity_ratio    = 0.0f;
    bool   channel_c_valid        = false;
    float  ca_stop_rate_baseline  = 0.0f;
    float  ca_stop_rate_terminal  = 0.0f;
    float  ca_uniformity_ratio    = 0.0f;
    float  channel_f_z            = 0.0f;
    float  channel_f_mh_z         = 0.0f;
    float  channel_f_common_or    = 0.0f;
    bool   channel_f_valid        = false;
    float  cg_stop_rate_baseline  = 0.0f;
    float  cg_stop_rate_terminal  = 0.0f;
    float  cg_uniformity_ratio    = 0.0f;
    float  channel_g_z            = 0.0f;
    bool   channel_g_valid        = false;
    float  at_stop_rate_baseline  = 0.0f;
    float  at_stop_rate_terminal  = 0.0f;
    float  at_uniformity_ratio    = 0.0f;
    float  channel_h_z            = 0.0f;
    float  channel_h_z_p2plus     = 0.0f;
    bool   channel_h_valid        = false;
};

// Length-stratified result: a vector of bins + joint mixture model.
struct LengthStratifiedDamageProfile {
    std::string         method = "single";
    std::vector<int>    edges;
    std::vector<LengthBinDamageProfile> bins;
    int64_t reads_scanned = 0;
    int     min_length    = 0;
    int     max_length    = 0;

    // Joint length × GC 2-component mixture
    double  d_joint_ancient    = 0.0;
    double  pi_joint_ancient   = 0.0;
    double  d_joint_population = 0.0;
    bool    joint_converged    = false;
    bool    joint_separated    = false;
    std::vector<std::array<double, LengthBinDamageProfile::N_GC_BINS>> cell_w_ancient;
};

} // namespace taph
