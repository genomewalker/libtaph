#pragma once

#include "types.hpp"
#include "joint_damage_model.hpp"
#include "mixture_damage_model.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <utility>

namespace dart {

// Forward declarations
class DamageModel;
struct UnifiedDamageContext;

struct SampleDamageProfile {
    // Position-specific base counts at 5' end (positions 0-14)
    // Using double to avoid float precision loss at >16M reads
    std::array<double, 15> t_freq_5prime = {};  // T count at each position
    std::array<double, 15> c_freq_5prime = {};  // C count (baseline)
    // Negative control counts at 5' end (for A/(A+G) ratio)
    std::array<double, 15> a_freq_5prime = {};  // A count at each position (control)
    std::array<double, 15> g_freq_5prime = {};  // G count at each position (control)

    // Position-specific base counts at 3' end (positions 0-14 from end)
    std::array<double, 15> a_freq_3prime = {};  // A count at each position
    std::array<double, 15> g_freq_3prime = {};  // G count (baseline)
    // Negative control counts at 3' end (for T/(T+C) ratio)
    std::array<double, 15> t_freq_3prime = {};  // T count at each position (control)
    std::array<double, 15> c_freq_3prime = {};  // C count at each position (control)

    // Middle-of-read baseline counts (undamaged)
    double baseline_t_freq = 0.0;
    double baseline_c_freq = 0.0;
    double baseline_a_freq = 0.0;
    double baseline_g_freq = 0.0;

    // Computed damage rates (excess over baseline)
    std::array<float, 15> damage_rate_5prime = {};  // C→T rate at each 5' position
    std::array<float, 15> damage_rate_3prime = {};  // G→A rate at each 3' position

    // Codon-position-aware damage tracking (positions 1,2,3 in codon)
    // At 5' end: T/(T+C) ratio by codon position
    std::array<float, 3> codon_pos_t_rate_5prime = {0.5f, 0.5f, 0.5f};
    // At 3' end: A/(A+G) ratio by codon position
    std::array<float, 3> codon_pos_a_rate_3prime = {0.5f, 0.5f, 0.5f};
    // Raw counts for aggregation
    std::array<size_t, 3> codon_pos_t_count_5prime = {};
    std::array<size_t, 3> codon_pos_c_count_5prime = {};
    std::array<size_t, 3> codon_pos_a_count_3prime = {};
    std::array<size_t, 3> codon_pos_g_count_3prime = {};
    // Raw totals for significance testing
    std::array<double, 15> tc_total_5prime = {};  // T+C counts at 5'
    std::array<double, 15> ag_total_3prime = {};  // A+G counts at 3'

    // CpG context damage tracking
    float cpg_damage_rate = 0.0f;      // C→T rate in CpG context
    float non_cpg_damage_rate = 0.0f;  // C→T rate outside CpG
    size_t cpg_c_count = 0;            // C's in CpG context
    size_t cpg_t_count = 0;            // T's where C expected in CpG
    size_t non_cpg_c_count = 0;
    size_t non_cpg_t_count = 0;

    // Summary statistics
    float max_damage_5prime = 0.0f;  // Maximum C→T rate at position 0
    float max_damage_3prime = 0.0f;  // Maximum G→A rate at last position
    float sample_damage_prob = 0.0f; // Overall probability sample is ancient

    // Estimated decay constants (sample-specific)
    float lambda_5prime = 0.3f;  // Decay constant for 5' end (estimated from data)
    float lambda_3prime = 0.3f;  // Decay constant for 3' end (estimated from data)

    // Briggs damage model parameters: δ(pos) = δ_s·P_overhang(pos) + δ_d·(1-P_overhang(pos))
    float delta_s_5prime = 0.0f;  // Single-stranded deamination rate at 5' end
    float delta_d_5prime = 0.0f;  // Double-stranded (background) deamination at 5'
    float delta_s_3prime = 0.0f;  // Single-stranded deamination rate at 3' end
    float delta_d_3prime = 0.0f;  // Double-stranded (background) deamination at 3'
    float r_squared_5prime = 0.0f;  // Goodness of fit for 5' model
    float r_squared_3prime = 0.0f;  // Goodness of fit for 3' model

    // Codon-position-specific damage analysis
    float codon_pos1_damage = 0.0f;  // C→T rate at codon position 1
    float codon_pos2_damage = 0.0f;  // C→T rate at codon position 2
    float codon_pos3_damage = 0.0f;  // C→T rate at codon position 3 (wobble)
    float wobble_ratio = 1.0f;       // pos3 / ((pos1 + pos2) / 2), >1 indicates ancient
    float hexamer_damage_llr = 0.0f; // Hexamer-based damage log-likelihood ratio
    float terminal_shift_5prime = 0.0f;  // terminal T/(T+C) - interior baseline
    float terminal_shift_3prime = 0.0f;  // terminal A/(A+G) - interior baseline
    float terminal_z_5prime = 0.0f;      // z-score for 5' terminal enrichment
    float terminal_z_3prime = 0.0f;      // z-score for 3' terminal enrichment
    bool terminal_inversion = false;     // true if terminals show significant depletion
    bool position_0_artifact_5prime = false;  // pos0 depleted but pos1 enriched (adapter bias)
    bool position_0_artifact_3prime = false;  // pos0 depleted but pos1 enriched (adapter bias)

    // Negative control statistics (should NOT show enrichment if damage is real)
    // 5' control: A/(A+G) at 5' end - real C→T damage shouldn't affect this
    // 3' control: T/(T+C) at 3' end - real G→A damage shouldn't affect this
    float ctrl_shift_5prime = 0.0f;   // 5' A/(A+G) terminal - interior
    float ctrl_shift_3prime = 0.0f;   // 3' T/(T+C) terminal - interior
    float ctrl_z_5prime = 0.0f;       // z-score for 5' control enrichment
    float ctrl_z_3prime = 0.0f;       // z-score for 3' control enrichment
    bool composition_bias_5prime = false;  // control comparable to damage → bias
    bool composition_bias_3prime = false;  // control comparable to damage → bias

    // 3' control channel per-position data (T/(T+C) at each position from 3' end).
    // For DS libraries this is the negative control; for SS it is the damage signal.
    // tc_total_3prime[p] = T+C count at position p from 3' (coverage gate for fqdup masking).
    std::array<double, 15> tc_total_3prime = {};

    // Library-type classifier: 4-channel joint BIC.
    // Channels: ct5 = 5' C→T, ga3 = 3' G→A smooth, ga0 = 3' G→A pos-0 spike, ct3 = 3' C→T.
    // BIC stored as double: values reach ~1e9 at high coverage, exceeding float precision.
    float  libtype_amp_ct5  = 0.0f;   // 5' C→T fitted amplitude
    float  libtype_amp_ga3  = 0.0f;   // 3' G→A smooth decay amplitude (pos 1-10)
    float  libtype_amp_ga0  = 0.0f;   // 3' G→A pos-0 spike amplitude
    float  libtype_amp_ct3  = 0.0f;   // 3' C→T fitted amplitude
    double libtype_dbic_ct5 = 0.0;    // ΔBIC ct5  (positive = decay favoured over flat)
    double libtype_dbic_ga3 = 0.0;
    double libtype_dbic_ga0 = 0.0;
    double libtype_dbic_ct3 = 0.0;
    double library_bic_bias = 0.0;  // M_bias:   ct5=null, ga3=null, ga0=null, ct3=null
    double library_bic_ds   = 0.0;  // M_DS:     ct5=alt,  ga3=alt,  ga0=null, ct3=null
    double library_bic_ss   = 0.0;  // best SS model BIC (min of SS_comp, SS_orig, SS_full)
    double library_bic_mix  = 0.0;  // M_SS_full: ct5=alt, ga3=alt, ga0=alt,  ct3=alt

    // Library type detection
    enum class LibraryType { UNKNOWN, DOUBLE_STRANDED, SINGLE_STRANDED };
    LibraryType library_type = LibraryType::DOUBLE_STRANDED;  // Default to double-stranded
    LibraryType forced_library_type = LibraryType::UNKNOWN;  // User override (UNKNOWN = auto-detect)
    bool library_type_auto_detected = false;  // true when set by auto-detection, not user override

    // Inverted pattern: terminal T/(T+C) < interior (reference-free detection unreliable)
    bool inverted_pattern_5prime = false;  // 5' terminal T/(T+C) < interior
    bool inverted_pattern_3prime = false;  // 3' terminal A/(A+G) < interior
    float terminal_gradient_5prime = 0.0f;  // pos0 - pos10-14 average (negative = inverted)
    float terminal_gradient_3prime = 0.0f;  // pos0 - pos10-14 average (negative = inverted)

    // Exponential fit parameters: p(pos) = b + A * exp(-lambda * pos)
    // Fitted baseline (asymptotic T/(T+C) or A/(A+G))
    float fit_baseline_5prime = 0.0f;   // b parameter for 5' end
    float fit_baseline_3prime = 0.0f;   // b parameter for 3' end
    // Fitted amplitude (damage signal above baseline)
    float fit_amplitude_5prime = 0.0f;  // A parameter for 5' end
    float fit_amplitude_3prime = 0.0f;  // A parameter for 3' end
    // Fit quality (RMSE of residuals)
    float fit_rmse_5prime = 0.0f;
    float fit_rmse_3prime = 0.0f;

    // Detected adapter offset: start_pos that gave the best BIC channel fit.
    // 1 = no adapter (default); 2 = 1-bp adapter remnant shifted signal to pos 1;
    // 3 = 2-bp adapter remnant shifted signal to pos 2.
    // Used to correct d_max estimation when terminal positions carry adapter sequence
    // rather than biological damage signal.
    int fit_offset_5prime = 1;
    int fit_offset_3prime = 1;

    // Calibrated D_max values (comparable to metaDMG)
    // D = A / (1 - b), the fraction of C that became T
    float d_max_5prime = 0.0f;  // Calibrated D_max for 5' end
    float d_max_3prime = 0.0f;  // Calibrated D_max for 3' end
    float d_max_combined = 0.0f;  // Final D_max using asymmetry-aware combination
    float asymmetry = 0.0f;  // |D_5p - D_3p| / ((D_5p + D_3p) / 2)
    bool high_asymmetry = false;  // True if asymmetry > 0.5 (possible artifact)

    // Track source of d_max_combined estimate
    enum class DmaxSource { AVERAGE, MIN_ASYMMETRY, MAX_SS_ASYMMETRY, FIVE_PRIME_ONLY, THREE_PRIME_ONLY, CHANNEL_B_STRUCTURAL, CHANNEL_B3_STRUCTURAL, NONE };
    DmaxSource d_max_source = DmaxSource::AVERAGE;

    const char* d_max_source_str() const {
        switch (d_max_source) {
            case DmaxSource::AVERAGE: return "average";
            case DmaxSource::MIN_ASYMMETRY: return "min_asymmetry";
            case DmaxSource::MAX_SS_ASYMMETRY: return "max_ss_asymmetry";
            case DmaxSource::FIVE_PRIME_ONLY: return "5prime_only";
            case DmaxSource::THREE_PRIME_ONLY: return "3prime_only";
            case DmaxSource::CHANNEL_B_STRUCTURAL: return "channel_b_structural";
            case DmaxSource::CHANNEL_B3_STRUCTURAL: return "channel_b3_structural";
            case DmaxSource::NONE: return "none";
            default: return "unknown";
        }
    }

    size_t n_reads = 0;  // Number of reads used in computation
    size_t n_reads_gc_filtered = 0;  // Reads skipped due to low GC content
    size_t n_reads_sampled = 0;  // Total reads sampled for GC histogram

    // GC content histogram (100 bins: 0-1%, 1-2%, ..., 99-100%)
    // Used to compute adaptive GC threshold for damage detection
    std::array<size_t, 100> gc_histogram = {};
    float adaptive_gc_threshold = 0.0f;  // Computed from histogram (70th percentile)
    bool gc_threshold_computed = false;

    // Hexamer-based damage detection
    // Track hexamer counts at 5' terminal positions (first 6 bases)
    // For each hexamer, we count occurrences at terminal vs interior positions
    // C→T damage should show excess T-hexamers at terminals relative to expected
    std::array<double, 4096> hexamer_count_5prime = {};  // Hexamer counts at 5' (pos 0-5)
    std::array<double, 4096> hexamer_count_interior = {}; // Hexamer counts at interior
    size_t n_hexamers_5prime = 0;    // Total hexamers counted at 5' terminal
    size_t n_hexamers_interior = 0;  // Total hexamers counted at interior

    // Hexamer-based T/(T+C) ratios (more reliable than position 0 or 1 alone)
    // These average positions 1-6 and are less affected by first-base artifacts
    float hexamer_terminal_tc = 0.0f;   // T/(T+C) at terminal from hexamer analysis
    float hexamer_interior_tc = 0.0f;   // T/(T+C) at interior from hexamer analysis
    float hexamer_excess_tc = 0.0f;     // Terminal - interior (negative = inverted)

    // Likelihood-based model comparison (exponential decay vs constant)
    // Positive LLR = exponential fits better (real decay pattern)
    // Negative/zero LLR = constant fits better (no decay, likely composition bias)
    float decay_llr_5prime = 0.0f;      // Log-likelihood ratio for 5' damage channel
    float decay_llr_3prime = 0.0f;      // Log-likelihood ratio for 3' damage channel
    float ctrl_decay_llr_5prime = 0.0f; // Log-likelihood ratio for 5' control channel
    float ctrl_decay_llr_3prime = 0.0f; // Log-likelihood ratio for 3' control channel

    float delta_llr_5prime = 0.0f;      // decay_llr - ctrl_decay_llr at 5'
    float delta_llr_3prime = 0.0f;      // decay_llr - ctrl_decay_llr at 3'

    float channel_divergence_5prime = 0.0f;  // |damage_shift - control_shift| at 5'
    float channel_divergence_3prime = 0.0f;  // |damage_shift - control_shift| at 3'

    // Channel B: convertible stop codon counts at 5' end by nucleotide position (0-14)
    // Position = nucleotide position of the C/T in the codon (from read start)
    // For CAA/TAA: position of the first base (C or T)
    // Exposure = CAA + TAA, Stops = TAA
    std::array<double, 15> convertible_caa_5prime = {};  // CAA (Gln) codons
    std::array<double, 15> convertible_taa_5prime = {};  // TAA (Stop) codons
    std::array<double, 15> convertible_cag_5prime = {};  // CAG (Gln) codons
    std::array<double, 15> convertible_tag_5prime = {};  // TAG (Stop) codons
    std::array<double, 15> convertible_cga_5prime = {};  // CGA (Arg) codons
    std::array<double, 15> convertible_tga_5prime = {};  // TGA (Stop) codons

    // Total codons observed at each position (denominator for exposure)
    std::array<double, 15> total_codons_5prime = {};

    // Interior reference counts (for baseline estimation, positions 30+)
    double convertible_caa_interior = 0.0;
    double convertible_taa_interior = 0.0;
    double convertible_cag_interior = 0.0;
    double convertible_tag_interior = 0.0;
    double convertible_cga_interior = 0.0;
    double convertible_tga_interior = 0.0;
    double total_codons_interior = 0.0;

    // Computed statistics for Channel B
    float stop_conversion_rate_baseline = 0.0f;  // Interior stop/(pre+stop) ratio
    float stop_decay_llr_5prime = 0.0f;  // LLR for stop position decay (Channel B)
    float stop_amplitude_5prime = 0.0f;  // Fitted amplitude of stop excess
    bool channel_b_valid = false;  // True if sufficient data for Channel B

    // Channel B structural d_max from multi-position stop codon conversion
    // WLS model: r_p = b0 + (1-b0) * d_max * exp(-λp)
    float d_max_from_channel_b = 0.0f;   // Structural d_max estimate from stop codons
    float channel_b_weight = 0.0f;       // Exposure weight W_B for joint likelihood
    float channel_b_slope = 0.0f;        // Raw WLS slope (positive = damage, negative = inverted)
    bool channel_b_quantifiable = false; // True if Channel B can provide d_max estimate
    bool channel_b_inverted = false;     // True if slope <= 0 (terminal stops LOWER than baseline)

    // Channel B₃': G→A stop codon conversion at 3' end (validates SS library damage)
    // TGG (Trp) is the only non-stop codon convertible to a stop via single G→A:
    //   TGG + b1 G→A → TAG (amber stop)
    //   TGG + b2 G→A → TGA (opal stop)
    // Position p = nucleotide distance of codon's last base from the 3' terminus
    std::array<double, 15> convertible_tgg_3prime = {};      // TGG (Trp) codons
    std::array<double, 15> convertible_tag_ga_3prime = {};   // TAG from TGG b1 G→A
    std::array<double, 15> convertible_tga_ga_3prime = {};   // TGA from TGG b2 G→A

    double convertible_tgg_interior = 0.0;
    double convertible_tag_ga_interior = 0.0;
    double convertible_tga_ga_interior = 0.0;

    float stop_conversion_rate_baseline_3prime = 0.0f;
    float stop_decay_llr_3prime = 0.0f;
    float stop_amplitude_3prime = 0.0f;
    bool  channel_b3_valid = false;

    float d_max_from_channel_b3 = 0.0f;
    float channel_b3_weight = 0.0f;
    float channel_b3_slope = 0.0f;
    bool  channel_b3_quantifiable = false;
    bool  channel_b3_inverted = false;

    // Channel C: oxidative stop codon tracking (G→T transversions, uniform across reads)
    std::array<double, 15> convertible_gag_5prime = {};      // GAG (Glu) codons at 5'
    std::array<double, 15> convertible_tag_ox_5prime = {};   // TAG (Stop) from G→T at 5'
    std::array<double, 15> convertible_gaa_5prime = {};      // GAA (Glu) codons at 5'
    std::array<double, 15> convertible_taa_ox_5prime = {};   // TAA (Stop) from G→T at 5'
    std::array<double, 15> convertible_gga_5prime = {};      // GGA (Gly) codons at 5'
    std::array<double, 15> convertible_tga_ox_5prime = {};   // TGA (Stop) from G→T at 5'

    double convertible_gag_interior = 0.0;
    double convertible_tag_ox_interior = 0.0;
    double convertible_gaa_interior = 0.0;
    double convertible_taa_ox_interior = 0.0;
    double convertible_gga_interior = 0.0;
    double convertible_tga_ox_interior = 0.0;

    float ox_stop_conversion_rate_baseline = 0.0f;
    float ox_stop_rate_terminal = 0.0f;
    float ox_stop_rate_interior = 0.0f;
    float ox_uniformity_ratio = 0.0f;   // terminal/interior (≈1 = uniform = real oxidation)
    bool channel_c_valid = false;

    // Channel D: G→T / C→A transversion tracking (8-oxoG, uniform across read)
    std::array<double, 15> g_count_5prime = {};     // G count at each 5' position
    std::array<double, 15> t_from_g_5prime = {};    // T where G expected (from G→T)
    std::array<double, 15> c_count_ox_5prime = {};  // C count for oxidation tracking
    std::array<double, 15> a_from_c_5prime = {};    // A where C expected (from C→A)

    double baseline_g_to_t_count = 0.0;
    double baseline_g_total = 0.0;
    double baseline_c_to_a_count = 0.0;
    double baseline_c_ox_total = 0.0;

    float ox_gt_rate_terminal = 0.0f;
    float ox_gt_rate_interior = 0.0f;
    float ox_gt_baseline = 0.0f;
    float ox_gt_uniformity = 0.0f;
    float ox_gt_asymmetry = 0.0f;

    float ox_ca_rate_terminal = 0.0f;
    float ox_ca_rate_interior = 0.0f;
    float ox_ca_baseline = 0.0f;
    float ox_ca_uniformity = 0.0f;

    float ox_d_max = 0.0f;
    bool ox_damage_detected = false;
    bool ox_is_artifact = false;

    // GT exponential-background fit: GT(p) = A*exp(-mu*p) + B
    // B = uniform background (8-oxoG estimate); A = terminal excess (artifact)
    float g_bg_fitted  = 0.0f;   // B: uniform G→T background
    float g_term_fitted = 0.0f;  // A: terminal excess
    float g_decay_fitted = 0.0f; // mu: terminal decay rate
    float s_gt = 0.0f;           // B - ox_ca_baseline: Chargaff contrast (signal for SS; ~0 for DS)

    // Channel E: depurination (purine loss at strand breaks)
    float purine_rate_terminal_5prime = 0.0f;
    float purine_rate_interior = 0.0f;
    float purine_enrichment_5prime = 0.0f;
    float purine_enrichment_3prime = 0.0f;
    bool depurination_detected = false;

    // GC-stratified damage: separate estimation per GC bin (interior GC to avoid bias)
    static constexpr int N_GC_BINS = 10;  // 0-10%, 10-20%, ..., 90-100%

    struct GCBinStats {
        // Channel A: T and C counts at terminal positions (0-14)
        std::array<uint64_t, 15> t_counts = {};
        std::array<uint64_t, 15> c_counts = {};
        // Channel A: interior baseline counts
        uint64_t t_interior = 0;
        uint64_t c_interior = 0;

        // Channel B: stop codon counts at terminal positions
        std::array<uint64_t, 15> stop_counts = {};  // TAA+TAG+TGA
        std::array<uint64_t, 15> pre_counts = {};   // CAA+CAG+CGA
        // Channel B: interior baseline
        uint64_t stop_interior = 0;
        uint64_t pre_interior = 0;

        // Computed results
        float d_max = 0.0f;           // Estimated damage for this bin
        float d_max_channel_b = 0.0f; // Channel B estimate for this bin
        uint64_t n_reads = 0;         // Number of reads in this bin
        uint64_t c_sites = 0;         // Total C sites (for weighting)
        bool valid = false;           // Sufficient data for estimation

        // GC-conditional damage parameters (for per-read inference)
        float p_damaged = 0.0f;       // P(damaged) for this bin from LLR classification
        float baseline_tc = 0.5f;     // Interior T/(T+C) baseline for this bin
        float llr = 0.0f;             // Log-likelihood ratio (positive = damaged)
        bool classified_damaged = false;  // Hard classification from LLR threshold

        // Terminal observation counts (for per-read LLR update)
        uint64_t n_terminal_obs() const {
            uint64_t n = 0;
            for (int p = 0; p < 15; ++p) n += t_counts[p] + c_counts[p];
            return n;
        }
    };

    std::array<GCBinStats, N_GC_BINS> gc_bins = {};

    // Aggregated GC-stratified results
    float gc_stratified_d_max_weighted = 0.0f;  // Weighted average across bins
    float gc_stratified_d_max_peak = 0.0f;      // Max d_max across valid bins
    int gc_peak_bin = -1;                        // Which bin has peak damage
    bool gc_stratified_valid = false;            // At least one bin has valid estimate

    float pi_damaged = 0.0f;          // Fraction of terminal obs from damaged bins
    float d_ancient = 0.0f;           // E[δ | damaged bins] - severity among damaged
    float d_population = 0.0f;        // E[δ] over all bins - average across all DNA
    int n_damaged_bins = 0;           // Number of bins classified as damaged

    // Joint evidence decision (legacy two-channel)
    bool damage_validated = false;  // True if both channels agree on damage
    bool damage_artifact = false;   // True if Channel A fires but Channel B doesn't

    // Joint probabilistic model results (BIC comparison of damage vs no-damage)
    float joint_delta_max = 0.0f;      // MLE estimate of damage rate
    float joint_lambda = 0.0f;         // Decay constant
    float joint_a_max = 0.0f;          // Artifact amplitude (signed)
    float joint_log_lik_m1 = 0.0f;     // Log-likelihood for M1 (damage)
    float joint_log_lik_m0 = 0.0f;     // Log-likelihood for M0 (no damage)
    float joint_delta_bic = 0.0f;      // BIC_M0 - BIC_M1 (positive = damage)
    float joint_bayes_factor = 0.0f;   // BF_10 ≈ exp(ΔBIC/2)
    float joint_p_damage = 0.0f;       // P(damage | data)
    bool joint_model_valid = false;    // Sufficient data for joint model

    // Mixture model results (K-component EM over GC-stratified bins)
    int mixture_K = 0;                 // Number of classes selected by BIC
    float mixture_d_population = 0.0f; // E[δ] over all C-sites
    float mixture_d_ancient = 0.0f;    // E[δ | δ > 5%] (ancient tail)
    float mixture_d_reference = 0.0f;  // E[δ | GC > 50%] (metaDMG proxy)
    float mixture_pi_ancient = 0.0f;   // Fraction of C-sites in high-damage classes
    float mixture_bic = 0.0f;          // BIC for model selection
    bool mixture_converged = false;    // Did EM converge?

    // Alignability-weighted damage estimate (proxy for reference-based tools)
    float d_metamatch = 0.0f;              // Calibrated metaDMG-comparable estimate
    float d_alignability_weighted = 0.0f;  // Raw alignability-weighted d_max
    float metamatch_gamma = 0.0f;          // Blending coefficient (0 = use d_global, 1 = use weighted)
    float mean_alignability = 0.0f;        // Mean alignability score across reads
    float alignability_damage_corr = 0.0f; // Correlation between alignability and per-read damage

    // Alignability-weighted accumulators (for incremental computation)
    double alignability_weighted_t_sum = 0.0;  // Σ(alignability × T_terminal)
    double alignability_weighted_c_sum = 0.0;  // Σ(alignability × C_terminal)
    double alignability_sum = 0.0;             // Σ(alignability)
    double alignability_sq_sum = 0.0;          // Σ(alignability²) for variance

    // Compute adaptive GC threshold from histogram
    float compute_adaptive_gc_threshold(float target_percentile = 0.70f) {
        size_t total = 0;
        for (auto c : gc_histogram) total += c;
        if (total < 1000) return 0.40f;  // Default if insufficient data

        size_t target_count = static_cast<size_t>(total * target_percentile);
        size_t cumulative = 0;
        for (int i = 0; i < 100; ++i) {
            cumulative += gc_histogram[i];
            if (cumulative >= target_count) {
                return static_cast<float>(i) / 100.0f;
            }
        }
        return 0.50f;  // Fallback
    }

    bool is_valid() const { return n_reads >= 1000; }

    bool is_detection_unreliable() const {
        return inverted_pattern_5prime || inverted_pattern_3prime ||
               composition_bias_5prime || composition_bias_3prime;
    }

    // Library type string for output
    const char* library_type_str() const {
        switch (library_type) {
            case LibraryType::DOUBLE_STRANDED: return "double-stranded";
            case LibraryType::SINGLE_STRANDED: return "single-stranded";
            default: return "unknown";
        }
    }

    // Get GC bin index (0-9) for a sequence based on interior GC content
    static int get_gc_bin(const std::string& seq) {
        if (seq.length() < 20) return 4;
        size_t gc = 0, total = 0;
        size_t start = std::min(size_t(5), seq.length() / 4);
        size_t end = seq.length() - start;

        for (size_t i = start; i < end; ++i) {
            char c = seq[i];
            if (c == 'G' || c == 'g' || c == 'C' || c == 'c') ++gc;
            if (c == 'A' || c == 'a' || c == 'T' || c == 't' ||
                c == 'G' || c == 'g' || c == 'C' || c == 'c') ++total;
        }

        if (total == 0) return 4;
        float gc_frac = static_cast<float>(gc) / static_cast<float>(total);
        int bin = static_cast<int>(gc_frac * 10.0f);
        return std::clamp(bin, 0, N_GC_BINS - 1);
    }

    // Get GC-conditional damage parameters for a read
    struct GCDamageParams {
        float p_damaged;     // Prior P(damaged) from bin classification
        float delta_s;       // Damage rate δ_s for this bin
        float baseline_tc;   // Interior T/(T+C) baseline
        float lambda;        // Decay constant (shared across bins)
        bool bin_valid;      // Whether this bin has valid estimates
    };

    GCDamageParams get_gc_params(const std::string& seq) const {
        int bin = get_gc_bin(seq);
        const auto& stats = gc_bins[bin];

        GCDamageParams params;
        params.p_damaged = stats.p_damaged;
        params.delta_s = stats.d_max;
        params.baseline_tc = stats.baseline_tc;
        params.lambda = lambda_5prime;  // Use fitted lambda
        params.bin_valid = stats.valid;

        return params;
    }

    // Get effective damage rate for a read (posterior-weighted)
    // δ_eff = P(damaged|read) * δ_s
    float get_effective_damage(const std::string& seq, float read_ancient_prob) const {
        int bin = get_gc_bin(seq);
        return read_ancient_prob * gc_bins[bin].d_max;
    }
};

// Tri-state validation: VALIDATED (both channels), CONTRADICTED (A only), UNVALIDATED (no B data)
enum class DamageValidationState {
    VALIDATED,      // Both channels fired → full damage
    CONTRADICTED,   // Channel A fired but Channel B negative → artifact
    UNVALIDATED     // Insufficient data for Channel B → soft suppression
};

inline DamageValidationState get_damage_validation_state(const SampleDamageProfile& profile) {
    if (profile.damage_validated) {
        return DamageValidationState::VALIDATED;
    }
    if (profile.damage_artifact) {
        return DamageValidationState::CONTRADICTED;
    }
    // Channel B has sufficient data but contradicts Channel A
    if (profile.channel_b_valid && profile.stop_decay_llr_5prime < -100.0f) {
        return DamageValidationState::CONTRADICTED;
    }
    return DamageValidationState::UNVALIDATED;
}

inline float get_damage_suppression_factor(DamageValidationState state) {
    switch (state) {
        case DamageValidationState::VALIDATED:    return 1.0f;
        case DamageValidationState::CONTRADICTED: return 0.0f;
        case DamageValidationState::UNVALIDATED:  return 0.5f;
    }
    return 1.0f;  // Unreachable, but satisfies compiler
}

inline float get_damage_suppression_factor(const SampleDamageProfile& profile) {
    return get_damage_suppression_factor(get_damage_validation_state(profile));
}

} // namespace dart
