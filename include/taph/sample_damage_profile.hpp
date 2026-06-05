#pragma once
#include <cstdint>

#include "types.hpp"
#include "joint_damage_model.hpp"
#include "mixture_damage_model.hpp"
#include "bulk_damage_model.hpp"
#include "channel_count_table.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <limits>
#include <utility>

namespace taph {

// Forward declarations
class DamageModel;
struct UnifiedDamageContext;

struct SampleDamageProfile {
    static constexpr int N_POS = 15;
    static constexpr int N_CT_CTX = 2;
    static constexpr int N_OXOG16 = 16;
    static constexpr int N_UPSTREAM_CTX = 4;  // AC, CC, GC, TC (upstream base before C)
    static constexpr int N_OX_LEN_BINS = 4;   // <=50, <=75, <=100, >100
    static constexpr int N_OX_GC_BINS = 10;
    static constexpr int N_OX_DEAM_STRATA = 5;
    static constexpr int N_OX_BINS = N_OX_LEN_BINS * N_OX_GC_BINS;

    enum CtContext : int { CPG_LIKE = 0, NONCPG_LIKE = 1 };
    enum UpstreamContext : int { CTX_AC = 0, CTX_CC = 1, CTX_GC = 2, CTX_TC = 3 };

    // C2: emitted clamp for explosive/uncalibrated test statistics (binom_z,
    // MH z, etc.). Detection gates use z>3, well inside ±12, so clamping the
    // EMITTED magnitude does not change any detection decision; it only prevents
    // absurd sqrt(N)-scaled magnitudes (hundreds–thousands) from reaching
    // consumers as if calibrated. Exploratory, not a calibrated p-value.
    static constexpr float kZCap = 12.0f;
    // C2: raw binomial-LLR sums scale linearly with N (|LLR|≈1e5–1e6 on the panel)
    // and are NOT chi-squared(1) calibrated. Clamp the emitted LLR magnitude so
    // consumers do not read it as a p-value; the sign/effect-size remains informative.
    static constexpr float kLlrCap = 50.0f;

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

    // === Tail-anchored background tracking (chemistry-aware Briggs fit) ===
    // Positions BG_TAIL_LO..BG_TAIL_HI from each terminus give a stable
    // C->T (G->A) baseline far from terminal damage. Filled only when read
    // length permits (read_length > pos+1).
    static constexpr int BG_TAIL_LO = 20;
    static constexpr int BG_TAIL_HI = 49;
    static constexpr int BG_TAIL_N  = BG_TAIL_HI - BG_TAIL_LO + 1;  // 30
    std::array<double, BG_TAIL_N> tail_t_5prime  = {};
    std::array<double, BG_TAIL_N> tail_tc_5prime = {};
    std::array<double, BG_TAIL_N> tail_a_3prime  = {};
    std::array<double, BG_TAIL_N> tail_ag_3prime = {};

    // CpG context damage tracking
    float cpg_damage_rate = 0.0f;      // C→T rate in CpG context
    float non_cpg_damage_rate = 0.0f;  // C→T rate outside CpG
    size_t cpg_c_count = 0;            // C's in CpG context
    size_t cpg_t_count = 0;            // T's where C expected in CpG
    size_t non_cpg_c_count = 0;
    size_t non_cpg_t_count = 0;

    // 5' context-split C/T accumulators (CpG-like vs non-CpG-like)
    // double to avoid float precision loss at >16M observations per position
    std::array<std::array<double, N_POS>, N_CT_CTX> ct_ctx_t_5prime = {};
    std::array<std::array<double, N_POS>, N_CT_CTX> ct_ctx_total_5prime = {};

    // Interior baseline accumulators for context split
    std::array<double, N_CT_CTX> ct_ctx_t_interior = {};
    std::array<double, N_CT_CTX> ct_ctx_total_interior = {};

    // Fitted baselines and amplitudes for context split
    float fit_baseline_ct5_cpg_like    = std::numeric_limits<float>::quiet_NaN();
    float fit_baseline_ct5_noncpg_like = std::numeric_limits<float>::quiet_NaN();
    float dmax_ct5_cpg_like    = std::numeric_limits<float>::quiet_NaN();
    float dmax_ct5_noncpg_like = std::numeric_limits<float>::quiet_NaN();
    float cpg_ratio     = std::numeric_limits<float>::quiet_NaN();
    float log2_cpg_ratio = std::numeric_limits<float>::quiet_NaN();

    // Coverage diagnostics for context split
    float effcov_ct5_cpg_like_terminal    = 0.0f;
    float effcov_ct5_noncpg_like_terminal = 0.0f;
    float effcov_ct5_cpg_like_interior    = 0.0f;
    float effcov_ct5_noncpg_like_interior = 0.0f;
    float cov_ct5_cpg_like_terminal       = 0.0f;
    float cov_ct5_noncpg_like_terminal    = 0.0f;
    float cov_ct5_cpg_like_interior       = 0.0f;
    float cov_ct5_noncpg_like_interior    = 0.0f;
    int   fit_positions_ct5_cpg_like    = 0;
    int   fit_positions_ct5_noncpg_like = 0;

    // --- Experimental (no stability guarantee) --------------------------------
    // Fields below may be renamed, retyped, or removed in any release.
    // Do not rely on them in production consumers without pinning the library version.

    // Upstream-context-aware C→T tracking
    // 4 disjoint bins by upstream base: AC, CC, GC, TC
    std::array<std::array<double, N_POS>, N_UPSTREAM_CTX> ct5_t_by_upstream = {};
    std::array<std::array<double, N_POS>, N_UPSTREAM_CTX> ct5_total_by_upstream = {};
    std::array<double, N_UPSTREAM_CTX> ct5_t_interior_by_upstream = {};
    std::array<double, N_UPSTREAM_CTX> ct5_total_interior_by_upstream = {};

    // Fitted amplitudes per upstream context (shrinkage-estimated)
    std::array<float, N_UPSTREAM_CTX> dmax_ct5_by_upstream = {
        std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::quiet_NaN()
    };
    std::array<float, N_UPSTREAM_CTX> baseline_ct5_by_upstream = {
        std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::quiet_NaN()
    };
    std::array<float, N_UPSTREAM_CTX> cov_ct5_terminal_by_upstream = {};
    std::array<float, N_UPSTREAM_CTX> cov_ct5_interior_by_upstream = {};

    // Derived contrasts
    float dipyr_contrast = std::numeric_limits<float>::quiet_NaN();  // mean(CC,TC) - mean(AC,GC); upstream (5') base; correct CPD geometry

    float context_heterogeneity_chi2 = 0.0f;  // C2: coverage-weighted index, clamped to kZCap (NOT a calibrated chi2(3))
    float context_heterogeneity_chi2_raw = 0.0f;  // C2: unclamped coverage-scaled sum (∝ reads; non-distribution-free)
    float context_heterogeneity_p    = 1.0f;  // nominal p-value only; not calibrated (prefer dipyr_contrast effect size)
    bool  context_heterogeneity_detected = false;  // true if p < 0.05
    bool  context_heterogeneity_computed = false;  // true only when valid_ctx_count==4 and mean_d>0.001

    // oxoG 16-context interior panel
    std::array<float, N_OXOG16> oxog16_t        = {};
    std::array<float, N_OXOG16> oxog16_a_rc     = {};
    std::array<float, N_OXOG16> s_oxog_16ctx    = {};
    std::array<float, N_OXOG16> cov_oxog_16ctx  = {};

    // Reference-free trinucleotide spectrum (64 contexts, prev*16 + mid*4 + next;
    // A=0,C=1,G=2,T=3). Counts observed trinucleotides at the 5' terminal zone
    // (read positions 1..4) and the interior null-distribution zone (pos 10..14).
    // Non-ACGT bases skip the sample. Mirror counters at 3' end use read-offset
    // positions 1..4 and 10..14 from the 3' terminus (orientation as observed —
    // downstream analysis reverse-complements contexts for strand collapsing).
    // Downstream post-processing contrasts terminal vs interior context rates:
    //   d_ct_ctx[XCY] = T_rate_term(XCY ∪ XTY) - T_rate_int(XCY ∪ XTY)
    //   d_gt_ctx[XGY] = T_rate_term(XGY ∪ XTY) - T_rate_int(XGY ∪ XTY)
    // These feed into the damage-context scores in library_interpretation.
    static constexpr int N_TRINUC = 64;
    std::array<uint64_t, N_TRINUC> tri_5prime_terminal = {};
    std::array<uint64_t, N_TRINUC> tri_5prime_interior = {};
    std::array<uint64_t, N_TRINUC> tri_3prime_terminal = {};
    std::array<uint64_t, N_TRINUC> tri_3prime_interior = {};

    // Per-position trinucleotide counts (positions 1..N_POS_TRI-1 from each end).
    // p=1 = second base from end (first with a valid left flank).
    // Encoding: prev*16 + mid*4 + next  (A=0,C=1,G=2,T=3).
    // 3' counts use read orientation (not reference); reverse-complement
    // contexts downstream for strand-collapsing.
    static constexpr int N_POS_TRI = 15;
    std::array<std::array<uint64_t, N_TRINUC>, N_POS_TRI> tri_5prime_pos = {};
    std::array<std::array<uint64_t, N_TRINUC>, N_POS_TRI> tri_3prime_pos = {};

    // Reference-free oxidation-like composition contrast.
    // Reads are first stratified by read-local terminal deamination excess.
    // Finalization then calibrates those strata against bin-matched interior
    // baselines and compares high-deamination vs low-deamination strata within
    // length x GC bins. The counters are composition-only: no reference-derived
    // substitution calls are implied.
    struct OxidationLikeBin {
        struct Stratum {
            double term_t5 = 0.0, term_tc5 = 0.0;
            double term_a3 = 0.0, term_ag3 = 0.0;
            double int_t = 0.0, int_tc = 0.0;
            double int_a = 0.0, int_ag = 0.0;
            double sig_t = 0.0, sig_tg = 0.0;
            double sig_a = 0.0, sig_ac = 0.0;
            double ctrl_a = 0.0, ctrl_at = 0.0;
            double ctrl_c = 0.0, ctrl_cg = 0.0;
            uint64_t reads = 0;
        };
        std::array<Stratum, N_OX_DEAM_STRATA> strata = {};
    };
    std::array<OxidationLikeBin, N_OX_BINS> oxidation_like_bins = {};
    float oxidation_like_signal = 0.0f;
    float oxidation_like_signal_se = 0.0f;
    float oxidation_like_control = 0.0f;
    float oxidation_like_control_se = 0.0f;
    float oxidation_like_adjusted = 0.0f;
    float oxidation_like_excess = 0.0f;
    float oxidation_like_se = 0.0f;
    float oxidation_like_z = 0.0f;
    float oxidation_like_reliability = 0.0f;
    int oxidation_like_bins_used = 0;
    float oxidation_like_effective_bins = 0.0f;
    float oxidation_like_heterogeneity = 0.0f;
    bool oxidation_like_artifact_suspect = false;

    // Interior clustered C→T: excess short-range co-occurrence of T at non-CpG {C,T} sites
    struct InteriorCtClusterAccumulator {
        uint64_t short_reads_skipped = 0;
        uint64_t reads_used_ct = 0;
        uint64_t reads_used_ag = 0;
        uint64_t obs_ct[11]   = {};
        uint64_t pairs_ct[11] = {};
        double   exp_ct[11]   = {};
        double   var_ct[11]   = {};
        uint64_t obs_ag[11]   = {};
        uint64_t pairs_ag[11] = {};
        double   exp_ag[11]   = {};
        double   var_ag[11]   = {};
    };
    InteriorCtClusterAccumulator interior_ct_cluster = {};

    float    interior_ct_cluster_short_log2oe       = 0.0f;  // CT log2(obs/exp), d=1..5
    float    interior_ct_cluster_short_asym_log2oe  = 0.0f;  // CT minus AG control
    float    interior_ct_cluster_short_z            = 0.0f;  // Heuristic z-score
    float    interior_ct_cluster_sep_log2oe[10]     = {};    // Per-separation d=1..10
    uint64_t interior_ct_cluster_short_obs          = 0;
    uint64_t interior_ct_cluster_short_pairs        = 0;
    double   interior_ct_cluster_short_exp          = 0.0;
    uint64_t interior_ct_cluster_reads_used         = 0;
    uint64_t interior_ct_cluster_reads_used_control = 0;
    uint64_t interior_ct_cluster_reads_skipped      = 0;

    // Summary statistics
    float max_damage_5prime = 0.0f;  // Maximum C→T rate at position 0
    float max_damage_3prime = 0.0f;  // Maximum G→A rate at last position
    float sample_damage_prob = 0.0f; // Overall probability sample is ancient

    // Estimated decay constants (sample-specific)
    float lambda_5prime = 0.3f;  // Decay constant for 5' end (estimated from data)
    float lambda_3prime = 0.3f;  // Decay constant for 3' end (estimated from data)
    bool  lambda_5prime_fitted = false;  // D22: true only when the fit updated lambda_5prime (distinguishes converged fit from 0.3 init)
    bool  lambda_3prime_fitted = false;  // D22: true only when the fit updated lambda_3prime

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
    bool library_type_rescued = false;  // true when DS call recovered from BIC failure via rescue rule

    // Posterior probabilities derived from the 7-model BIC tournament via
    // softmax(-BIC/2). Sum to 1.0 within numerical precision when evaluable.
    // p_winner = max of the three; library_type_confidence_threshold below
    // gates the UNKNOWN fallback. When library_type_evaluable=false (e.g. an
    // input channel was invalid), all four are 0 and library_type=UNKNOWN.
    float library_p_ds      = 0.0f;
    float library_p_ss      = 0.0f;
    float library_p_bias    = 0.0f;
    float library_p_winner  = 0.0f;
    // F3: post-hoc-vetoed final probabilities. When override fires, these are
    // one-hot to the surviving library_type; otherwise they mirror raw probs.
    float library_p_ds_final     = 0.0f;
    float library_p_ss_final     = 0.0f;
    float library_p_bias_final   = 0.0f;
    float library_p_winner_final = 0.0f;
    bool  library_type_evaluable = false;
    static constexpr float kLibraryTypeConfidenceThreshold = 0.60f;

    // === ADDITIVE: P0-1 forced-type diagnostics ===
    LibraryType library_auto_type = LibraryType::UNKNOWN;       // BIC tournament winner regardless of override
    bool        library_auto_evaluable = false;                 // true if the BIC tournament produced a valid call
    LibraryType library_forced_type = LibraryType::UNKNOWN;     // mirror of forced_library_type for downstream JSON

    // === ADDITIVE: P1-A calibrated call quality ===
    std::string library_bic_winner_model;     // e.g. "M_DS_symm_art"
    std::string library_bic_second_model;
    // NOTE: library_bic_margin scales linearly with read depth (ΔBIC ≈ n·2·ΔLLR/read);
    // it is NOT comparable across libraries of different coverage. Use the per-read
    // companion below for cross-library confidence comparison.
    double      library_bic_margin = 0.0;     // second_BIC - winner_BIC (>=0; bigger = more confident, same-library only)
    double      library_bic_margin_per_read = 0.0;  // margin / bg_denominator_5prime (coverage-normalised; comparable across libraries)
    // C5: top-level winner-model can contradict library_type after a post-hoc veto.
    // These co-locate the resolving signals so the contradiction is machine-visible
    // without reading the nested bic section.
    std::string library_bic_winner_model_class;   // "ds" | "ss" | "bias" of winner_model
    std::string library_bic_override_reason;      // mirror of final_library_bic_override_reason ("" if none)
    // C5: library_p_ds/p_ss/p_bias are PRE-override diagnostic posteriors; after a
    // veto they can disagree with library_type. True when an override fired, telling
    // consumers these raw probs must not be read as a type-confidence score (use *_final).
    bool        library_p_is_pre_override = false;
    float       library_p_ds_class_min   = 0.0f;  // softmax over best-BIC-per-class only
    float       library_p_ss_class_min   = 0.0f;
    float       library_p_bias_class_min = 0.0f;
    bool        library_artifact_contaminated = false;
    std::vector<std::string> library_artifact_reasons;

    // Protocol-tag classifier evidence: top 5' hexamer matches a known
    // library-prep oligo fingerprint (e.g. TGTAGA = SS Santa-Cruz / Gansauge,
    // CGATCT = DS TruSeq). When the BIC verdict disagrees with the chemistry
    // fingerprint, the chemistry is the stronger evidence and the cascade
    // overrides library_type. Fields populated whenever the top hex matches a
    // table entry; protocol_tag_applied=true only when the override fired.
    std::string protocol_tag_5prime;       // observed top 5' hexamer (e.g. "TGTAGA")
    std::string protocol_tag_protocol;     // matched protocol name
    LibraryType protocol_tag_class    = LibraryType::UNKNOWN;
    float       protocol_tag_log2fc   = 0.0f;
    float       protocol_tag_log_lr   = 0.0f;
    bool        protocol_tag_applied  = false;  // true when override flipped library_type

    // === Chemistry-aware Briggs fit + headline area-excess statistics ===
    // bg_*_anchored: trimmed-mean C->T (G->A) rate over tail positions
    // [BG_TAIL_LO..BG_TAIL_HI] using only positions with denom >= 100.
    // briggs_pos0_masked_*: pos 0 excluded from Briggs fit because protocol
    // tag dominates that position (5' ligation footprint suppresses real signal).
    // damage_*_area_excess: sum_{i=k}^{14} max(0, rate[i] - bg) where k=1 if
    // pos0 masked, else k=0. Robust headline statistic for cross-library
    // comparison; preferred over d_max when chemistry tag is set.
    // damage_*_lr: log-likelihood ratio of binomial(rate, bg) per-position
    // model vs bg-only null over positions k..14. Coverage-scaling companion.
    float  bg_5prime_anchored        = 0.0f;
    float  bg_3prime_anchored        = 0.0f;
    int    bg_n_positions_5prime     = 0;
    int    bg_n_positions_3prime     = 0;
    double bg_denominator_5prime     = 0.0;
    double bg_denominator_3prime     = 0.0;
    bool   briggs_pos0_masked_5prime = false;
    bool   briggs_pos0_masked_3prime = false;
    float  damage_5prime_area_excess = 0.0f;
    float  damage_3prime_area_excess = 0.0f;
    float  damage_5prime_lr          = 0.0f;
    float  damage_3prime_lr          = 0.0f;

    // Input mode: single-end (-i) or paired-end (-1/-2). Paired routes R1
    // 5' to per_pos_5prime, R2 5' to per_pos_3prime; read 3' ends ignored.
    enum class InputMode { SINGLE, PAIRED };
    InputMode input_mode = InputMode::SINGLE;
    uint64_t  pe_short_insert_skipped = 0;  // pairs skipped for short-insert read-through

    // === ADDITIVE: P1-B all 7 submodel BICs ===
    double library_bic_M_bias        = 0.0;
    double library_bic_M_DS_symm     = 0.0;
    double library_bic_M_DS_spike    = 0.0;
    double library_bic_M_DS_symm_art = 0.0;
    double library_bic_M_SS_comp     = 0.0;
    double library_bic_M_SS_orig     = 0.0;
    double library_bic_M_SS_asym     = 0.0;
    double library_bic_M_SS_full     = 0.0;
    bool   library_M_SS_orig_active  = false;  // ct3.delta_bic > 0
    bool   library_M_SS_asym_active  = false;  // spike_is_ss

    // === ADDITIVE: P2 spike-gate diagnostics ===
    bool  library_spike_is_ss = false;
    float library_spike_gate_ga0_amp = 0.0f;
    bool  library_spike_gate_structural_bilateral = false;
    bool  library_joint_lambda_restricted = false;

    // === ADDITIVE: S1 telemetry (frozen schema, never affects classification) ===
    // Counterfactual M_DS_asym (independent ct5/ga3 amps, no shared-amp constraint).
    // Telemetry only: never enters cascade or softmax.
    double library_bic_M_DS_asym_art = 0.0;
    // M_DS_symm rebuilt from a no-offset joint fit (start_pos=1 forced) — invariant
    // probe to catch joint best-offset regressions.
    double library_bic_M_DS_symm_art_no_offset = 0.0;

    // Raw 9-candidate winner ranking (8 cascade contenders + M_DS_asym_art),
    // ignores cascade gating and post-hoc rescue.
    std::string library_bic_raw_winner_model;
    std::string library_bic_raw_winner_class;   // "ds" | "ss" | "bias"
    std::string library_bic_raw_second_model;
    double      library_bic_raw_margin = 0.0;   // raw_2nd - raw_winner (>=0)
    bool        library_bic_raw_winner_in_cascade = false;

    // Cascade exclusion booleans (multi-valued: each gate evaluated independently).
    bool library_bic_excl_structural_bilateral = false;  // !structural_bilateral required for some gates
    bool library_bic_excl_spike_is_ss = false;           // spike_is_ss flipped routing
    bool library_bic_excl_ct3_zero = false;              // ct3.delta_bic <= 0 -> M_SS_orig inactive
    bool library_bic_excl_M_SS_full_hardcoded = false;   // M_SS_full structurally excluded
    bool library_bic_excl_in_cascade = false;            // raw winner WAS in cascade (no exclusion)

    // Post-hoc final winner (after veto/rescue logic mutates library_type).
    std::string final_library_bic_winner_model;
    std::string final_library_bic_override_reason;  // "" if final == cascade winner

    // Gate inputs (already computed; expose for audit).
    bool  library_gate_artifact_5 = false;
    bool  library_gate_artifact_3 = false;
    bool  library_gate_position_0_artifact_5prime = false;
    bool  library_gate_position_0_artifact_3prime = false;
    bool  library_gate_inverted_pattern_5prime = false;
    bool  library_gate_inverted_pattern_3prime = false;
    float library_gate_max_damage_5prime = 0.0f;
    bool  library_gate_ss_orientation_evidence = false;
    bool  library_gate_ga_spike_dominant = false;
    bool  library_gate_ds_spike_won = false;
    bool  library_gate_ga0_dominates_ct5 = false;
    bool  library_gate_structural_bilateral = false;

    // Per-channel offsets (already computed; expose).
    int library_ct3_offset = 1;
    int library_ds_symm_offset = 1;

    // ds_symm joint-fit diagnostics.
    float library_ds_symm_lambda_used = 0.0f;
    float library_ds_symm_amp = 0.0f;
    float library_ds_symm_ct5_resid = 0.0f;  // ct5.amplitude - ds_symm.amplitude
    float library_ds_symm_ga3_resid = 0.0f;  // ga3.amplitude - ds_symm.amplitude

    // Damage status: independent of library type, based on effect size + confidence interval
    enum class DamageStatus { ABSENT, WEAK, PRESENT };
    DamageStatus damage_status = DamageStatus::ABSENT;

    // C5: which evidence path set damage_status. Lets consumers tell an empirical
    // excess-CI confirmation from a fit-amplitude-only WEAK (inverted/inward 5'
    // pattern where both CI lower bounds return -1 but the fit amplitude is real).
    enum class DamageStatusBasis { NONE, CI_5PRIME, CI_3PRIME, FIT_AMPLITUDE_ONLY };
    DamageStatusBasis damage_status_basis = DamageStatusBasis::NONE;

    const char* damage_status_basis_str() const {
        switch (damage_status_basis) {
            case DamageStatusBasis::CI_5PRIME:          return "ci_5prime";
            case DamageStatusBasis::CI_3PRIME:          return "ci_3prime";
            case DamageStatusBasis::FIT_AMPLITUDE_ONLY: return "fit_amplitude_only";
            default:                                    return "none";
        }
    }

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

    // --- Stable outputs -------------------------------------------------------
    // Fields below are part of the public API and will not be renamed or removed
    // without a major version bump and changelog entry.

    // Calibrated D_max values (comparable to metaDMG)
    // D = A / (1 - b), the fraction of C that became T
    float d_max_5prime = 0.0f;  // Calibrated D_max for 5' end
    float d_max_3prime = 0.0f;  // Calibrated D_max for 3' end
    float d_max_combined = 0.0f;  // Final D_max using asymmetry-aware combination
    float asymmetry = 0.0f;  // |D_5p - D_3p| / ((D_5p + D_3p) / 2)
    bool high_asymmetry = false;  // True if asymmetry > 0.5 (possible artifact)

    // Per-position decay flatness metrics (populated in finalize_sample_profile).
    // A flat profile (max rate < 0.005 AND variance < 2.5e-5 over positions 0-9)
    // indicates the end signal is indistinguishable from noise — suggesting library
    // construction erasure (5' blunting) rather than genuinely zero ancient damage.
    bool  d5_profile_flat         = false;
    bool  d3_profile_flat         = false;
    float d5_max_rate_pos0_4      = 0.0f;  // max excess C→T over positions 0-4
    float d3_max_rate_pos0_4      = 0.0f;  // max excess rate over positions 0-4
    float d5_profile_var_pos0_9   = 0.0f;  // variance of rate over positions 0-9
    float d3_profile_var_pos0_9   = 0.0f;
    // d5 flat while d3 has real decay → 5' overhang blunted before adapter ligation
    bool  d5_blunting_suspected   = false;

    // Ancient-fraction d_max: damage estimated only from reads classified as ancient
    // by the per-read LLR scorer (fused into the oxoG pass inside fqdup profile).
    // Comparable to metaDMG per-reference "damaged fraction" values.
    bool    damaged_fraction_valid = false;
    float   damaged_fraction_d5    = 0.0f;
    float   damaged_fraction_d3    = 0.0f;
    float   damaged_fraction_pi    = 0.0f;   // fraction of reads classified as ancient
    int64_t damaged_fraction_n     = 0;      // count of reads classified as ancient
    float   modern_fraction_d5     = 0.0f;   // d_max_5prime for reads classified as modern
    float   modern_fraction_d3     = 0.0f;
    // Independently-fitted deamination profiles for each fraction (log-linear regression
    // over all N_POS positions). _fit fields are the fitted d_max; _lambda are decay rates.
    float   damaged_fraction_d5_fit  = 0.0f;
    float   damaged_fraction_lambda5 = 0.0f;
    float   damaged_fraction_d3_fit  = 0.0f;
    float   damaged_fraction_lambda3 = 0.0f;
    float   modern_fraction_d5_fit   = 0.0f;
    float   modern_fraction_lambda5  = 0.0f;
    float   modern_fraction_d3_fit   = 0.0f;
    float   modern_fraction_lambda3  = 0.0f;
    // Per-position raw rates for fraction damage curves (0 = no data / below coverage gate)
    std::array<float, 15> damaged_fraction_rate5{};  // T/TC(p) for ancient 5'
    std::array<float, 15> damaged_fraction_rate3{};  // signal/denom for ancient 3'
    std::array<float, 15> modern_fraction_rate5{};   // T/TC(p) for modern 5'
    std::array<float, 15> modern_fraction_rate3{};   // signal/denom for modern 3'
    bool    modern_fraction_leakage_5prime = false;
    bool    modern_fraction_leakage_3prime = false;
    bool    modern_fraction_d5_computed = false;  // true only when mod_tc5[0] >= 50
    bool    modern_fraction_d3_computed = false;  // true only when mod_n3[0]  >= 50

    // Track source of d_max_combined estimate
    enum class DmaxSource { AVERAGE, MIN_ASYMMETRY, MAX_SS_ASYMMETRY,
                            FIVE_PRIME_ONLY, THREE_PRIME_ONLY,
                            FIVE_PRIME_CONSERVATIVE_SS,
                            CHANNEL_B_STRUCTURAL, CHANNEL_B3_STRUCTURAL,
                            JOINT_BILATERAL, NONE };
    DmaxSource d_max_source = DmaxSource::AVERAGE;

    const char* d_max_source_str() const {
        switch (d_max_source) {
            case DmaxSource::AVERAGE:                  return "average";
            case DmaxSource::MIN_ASYMMETRY:            return "min_asymmetry";
            case DmaxSource::MAX_SS_ASYMMETRY:         return "max_ss_asymmetry";
            case DmaxSource::FIVE_PRIME_ONLY:          return "5prime_only";
            case DmaxSource::THREE_PRIME_ONLY:         return "3prime_only";
            case DmaxSource::FIVE_PRIME_CONSERVATIVE_SS: return "5prime_conservative_ss";
            case DmaxSource::CHANNEL_B_STRUCTURAL:     return "channel_b_structural";
            case DmaxSource::CHANNEL_B3_STRUCTURAL:    return "channel_b3_structural";
            case DmaxSource::JOINT_BILATERAL:          return "joint_bilateral_inversion";
            case DmaxSource::NONE:                     return "none";
            default:                                   return "unknown";
        }
    }

    size_t n_reads = 0;  // Number of reads used in computation
    size_t n_reads_gc_filtered = 0;  // Reads skipped due to low GC content
    size_t n_reads_sampled = 0;  // Total reads sampled for GC histogram

    // Lifecycle guard: finalize_sample_profile() converts raw count arrays
    // (t_freq_5prime, etc.) into rates IN PLACE. Calling finalize twice, or
    // merging a finalized profile, silently corrupts the math. This flag is
    // checked at the entry of finalize_sample_profile() and merge_sample_profiles()
    // so misuse becomes a hard error instead of a silent garbage result.
    bool finalized = false;

    // GC content histogram (100 bins: 0-1%, 1-2%, ..., 99-100%)
    // Used to compute adaptive GC threshold for damage detection
    std::array<size_t, 100> gc_histogram = {};
    float adaptive_gc_threshold = 0.0f;  // Computed from histogram (70th percentile)
    bool gc_threshold_computed = false;

    // Hexamer-based damage detection
    // Track hexamer counts at 5' terminal positions (first 6 bases)
    // For each hexamer, we count occurrences at terminal vs interior positions
    // C→T damage should show excess T-hexamers at terminals relative to expected
    std::array<double, 4096> hexamer_count_5prime = {};   // Hexamer counts at 5' (pos 0-5)
    std::array<double, 4096> hexamer_count_interior = {};  // Hexamer counts at interior
    std::array<double, 4096> hexamer_count_3prime = {};    // Hexamer counts at 3' (last 6 bases)
    size_t n_hexamers_5prime = 0;    // Total hexamers counted at 5' terminal
    size_t n_hexamers_interior = 0;  // Total hexamers counted at interior
    size_t n_hexamers_3prime = 0;    // Total hexamers counted at 3' terminal

    // Hexamer-based T/(T+C) ratios (more reliable than position 0 or 1 alone)
    // These average positions 1-6 and are less affected by first-base artifacts
    float hexamer_terminal_tc = 0.0f;   // T/(T+C) at terminal from hexamer analysis
    float hexamer_interior_tc = 0.0f;   // T/(T+C) at interior from hexamer analysis
    float hexamer_excess_tc = 0.0f;     // Terminal - interior (negative = inverted)

    // Prefix-conditioned F/G/H terminal sums (positions 0–4 aggregated, keyed by
    // first-hexamer code at 5' end and last-hexamer code at 3' end).
    // Populated during update_sample_profile; consumed by
    // recompute_fgh_excluding_adapter_prefixes() to exclude adapter-matching reads.
    std::array<double, 4096> ca_pre_terminal_by_pfx   = {};  // Channel F 5' pre
    std::array<double, 4096> ca_stop_terminal_by_pfx  = {};  // Channel F 5' stop
    std::array<double, 4096> ca_deam_shadow_terminal_by_pfx = {};  // Channel F 5' C→T shadow (TTA/TTG/TAT/TGT)
    std::array<double, 4096> cg_pre_terminal_by_pfx   = {};  // Channel G 5' pre
    std::array<double, 4096> cg_stop_terminal_by_pfx  = {};  // Channel G 5' stop
    std::array<double, 4096> at_pre_terminal_by_pfx        = {};  // Channel H 5' pre  (pos 0-4)
    std::array<double, 4096> at_stop_terminal_by_pfx       = {};  // Channel H 5' stop (pos 0-4)
    std::array<double, 4096> at_pre_terminal_p2plus_by_pfx  = {};  // Channel H 5' pre  (pos 2-4)
    std::array<double, 4096> at_stop_terminal_p2plus_by_pfx = {};  // Channel H 5' stop (pos 2-4)
    std::array<double, 4096> ca_pre_terminal_3p_by_pfx  = {};  // Channel F 3' pre
    std::array<double, 4096> ca_stop_terminal_3p_by_pfx = {};  // Channel F 3' stop
    std::array<double, 4096> cg_pre_terminal_3p_by_pfx  = {};  // Channel G 3' pre
    std::array<double, 4096> cg_stop_terminal_3p_by_pfx = {};  // Channel G 3' stop
    std::array<double, 4096> at_pre_terminal_3p_by_pfx  = {};  // Channel H 3' pre
    std::array<double, 4096> at_stop_terminal_3p_by_pfx = {};  // Channel H 3' stop
    uint32_t fgh_adapter_prefixes_excluded = 0; // set by recompute_fgh_excluding_adapter_prefixes

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

    // --- Diagnostics ----------------------------------------------------------
    // Informational counters and intermediate statistics. Stable names but
    // values may change when algorithms are updated; not recommended as primary
    // analysis inputs.

    // Computed statistics for Channel B
    float stop_conversion_rate_baseline = 0.0f;  // Interior stop/(pre+stop) ratio
    float stop_decay_llr_5prime = 0.0f;  // LLR for stop position decay (Channel B, pooled)
    float stop_amplitude_5prime = 0.0f;  // Fitted amplitude of stop excess
    bool channel_b_valid = false;  // True if sufficient data for Channel B (pooled)

    // Per-stop-type Channel B LLRs (split by deamination parent codon)
    float stop_decay_llr_taa_5prime = 0.0f;  // CAA→TAA specific LLR
    float stop_decay_llr_tag_5prime = 0.0f;  // CAG→TAG specific LLR
    float stop_decay_llr_tga_5prime = 0.0f;  // CGA→TGA specific LLR
    // Exposure counts for each stop type (interior, used to gate per-type validity)
    double n_convertible_caa = 0.0;  // CAA+TAA interior total
    double n_convertible_cag = 0.0;  // CAG+TAG interior total
    double n_convertible_cga = 0.0;  // CGA+TGA interior total
    bool channel_b_valid_tga = false;  // True if CGA exposure sufficient for TGA-specific LLR
    // D24: symmetrical per-type validity flags (only valid_tga existed). True when
    // the type-specific interior exposure (n_convertible_caa / n_convertible_cag)
    // reaches 50; emitter nulls lrt_taa / lrt_tag otherwise.
    bool channel_b_valid_taa = false;  // True if CAA exposure sufficient for TAA-specific LLR
    bool channel_b_valid_tag = false;  // True if CAG exposure sufficient for TAG-specific LLR

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

    uint64_t convertible_gag_interior = 0;
    uint64_t convertible_tag_ox_interior = 0;
    uint64_t convertible_gaa_interior = 0;
    uint64_t convertible_taa_ox_interior = 0;
    uint64_t convertible_gga_interior = 0;
    uint64_t convertible_tga_ox_interior = 0;

    float ox_stop_conversion_rate_baseline = 0.0f;
    float ox_stop_rate_terminal = 0.0f;
    float ox_stop_rate_interior = 0.0f;
    float ox_uniformity_ratio = 0.0f;   // terminal/interior (≈1 = uniform = real oxidation)
    bool channel_c_valid = false;
    // Channel C positional-coverage validity (C1/C3). channel_c_valid gates only the
    // interior baseline; these gate the terminal/interior positional rates and the
    // uniformity ratio so the emitter can null an uncomputed value vs a measured one.
    bool ox_stop_rate_positional_computed = false;  // terminal+mid positional rates computed
    bool ox_uniformity_ratio_computed     = false;  // ox_uniformity_ratio actually measured

    // Channel C3': G→T stop codon conversion at 3' end (validates oxidation at 3' overhang)
    // Symmetric to Channel B3' (G→A stops at 3') but for G→T transversions.
    // GAG/GAA/GGA near 3' end → TAG/TAA/TGA if G→T oxidation occurred.
    // Position p = nucleotide distance of codon start from the 3' terminus.
    std::array<double, 15> convertible_gag_3prime = {};
    std::array<double, 15> convertible_tag_ox_3prime = {};
    std::array<double, 15> convertible_gaa_3prime = {};
    std::array<double, 15> convertible_taa_ox_3prime = {};
    std::array<double, 15> convertible_gga_3prime = {};
    std::array<double, 15> convertible_tga_ox_3prime = {};

    float ox_stop_rate_terminal_3prime = 0.0f;
    float ox_stop_rate_interior_3prime = 0.0f;
    float ox_stop_baseline_3prime = 0.0f;
    float ox_uniformity_ratio_3prime = 0.0f;
    bool  channel_c3_valid = false;
    // C1: ox_uniformity_ratio_3prime default 0.0f is an ambiguous sentinel; emitter
    // nulls it when this flag is false (mirrors the 5' ox_uniformity_ratio_computed).
    bool  ox_uniformity_ratio_3prime_computed = false;


    // Channel F: C→A oxidation stop codon tracking (complement-strand G→T).
    // Precursors where C→A at codon pos 1 or pos 2 creates a stop:
    //   TCA→TAA (Ser→Ochre, pos1), TCG→TAG (Ser→Amber, pos1)
    //   TAC→TAA (Tyr→Ochre, pos2), TGC→TGA (Cys→Opal,  pos2)
    // The Channel F/C ratio is a reference-free symmetric oxidation test.
    std::array<double, 15> convertible_tca_5prime = {};     // TCA (Ser) precursor
    std::array<double, 15> convertible_tcg_5prime = {};     // TCG (Ser) precursor
    std::array<double, 15> convertible_tac_5prime = {};     // TAC (Tyr) precursor
    std::array<double, 15> convertible_tgc_5prime = {};     // TGC (Cys) precursor
    std::array<double, 15> convertible_taa_ca_5prime = {};  // TAA stop from C→A
    std::array<double, 15> convertible_tag_ca_5prime = {};  // TAG stop from C→A
    std::array<double, 15> convertible_tga_ca_5prime = {};  // TGA stop from C→A
    std::array<double, 15> ca_deam_shadow_5prime = {};      // C→T shadows of TCA/TCG/TAC/TGC (TTA/TTG/TAT/TGT)
    std::array<double, 15> ca_shadow_5prime_ctx0 = {};  // TCA+TAC shadows (TTA+TAT), per position
    std::array<double, 15> ca_shadow_5prime_ctx1 = {};  // TCG shadows (TTG), per position
    std::array<double, 15> ca_shadow_5prime_ctx2 = {};  // TGC shadows (TGT), per position

    uint64_t convertible_tca_interior = 0;
    uint64_t convertible_tcg_interior = 0;
    uint64_t convertible_tac_interior = 0;
    uint64_t convertible_tgc_interior = 0;
    uint64_t convertible_taa_ca_interior = 0;
    uint64_t convertible_tag_ca_interior = 0;
    uint64_t convertible_tga_ca_interior = 0;

    // Channel F 3' end
    std::array<double, 15> convertible_tca_3prime = {};
    std::array<double, 15> convertible_tcg_3prime = {};
    std::array<double, 15> convertible_tac_3prime = {};
    std::array<double, 15> convertible_tgc_3prime = {};
    std::array<double, 15> convertible_taa_ca_3prime = {};
    std::array<double, 15> convertible_tag_ca_3prime = {};
    std::array<double, 15> convertible_tga_ca_3prime = {};
    std::array<double, 15> ca_deam_shadow_3prime = {};      // C→T shadows at 3' end

    float ca_stop_rate_baseline     = 0.0f;
    float ca_stop_rate_terminal     = 0.0f;
    float ca_stop_rate_interior     = 0.0f;
    // C1/C2: binom_z fields default to NaN (not 0.0f) so 'not computed' is
    // distinguishable from a genuine z==0; producer clamps the computed value to
    // [-kZCap, kZCap] (exploratory, not a calibrated p-value — correlated reads).
    float channel_f_z               = std::numeric_limits<float>::quiet_NaN();
    float channel_f_mh_z            = std::numeric_limits<float>::quiet_NaN();
    float channel_f_common_or       = 0.0f;
    float ca_uniformity_ratio       = 0.0f;
    float ca_stop_rate_baseline_3prime = 0.0f;
    float ca_stop_rate_terminal_3prime = 0.0f;
    float ca_stop_rate_interior_3prime = 0.0f;
    float ca_uniformity_ratio_3prime   = 0.0f;
    bool  channel_f_valid           = false;
    bool  channel_f3_valid          = false;

    // Layer-0 stop-channel count tables (F/G/H, 5' and 3') — the pre-clamp source of truth the
    // golden regression gate diffs. Additive diagnostics; does not affect any emitted estimator.
    std::vector<StopChannelCountTable> count_tables;

    // Channel G: C→G stop codon conversion (hydantoin-class advanced oxidation)
    // TCA (Ser) → TGA via C→G at codon pos 2
    // TAC (Tyr) → TAG via C→G at codon pos 3
    std::array<double, 15> convertible_tca_cg_5prime = {};
    std::array<double, 15> convertible_tac_cg_5prime = {};
    std::array<double, 15> convertible_tga_cg_5prime = {};
    std::array<double, 15> convertible_tag_cg_5prime = {};
    std::array<double, 15> convertible_tca_cg_3prime = {};
    std::array<double, 15> convertible_tac_cg_3prime = {};
    std::array<double, 15> convertible_tga_cg_3prime = {};
    std::array<double, 15> convertible_tag_cg_3prime = {};

    float cg_stop_rate_terminal    = 0.0f;
    float cg_stop_rate_interior    = 0.0f;
    float cg_stop_rate_baseline    = 0.0f;
    float channel_g_z              = std::numeric_limits<float>::quiet_NaN();  // C1/C2: NaN = not computed; clamped when computed
    float cg_uniformity_ratio      = 0.0f;
    float cg_stop_rate_terminal_3prime = 0.0f;
    float cg_stop_rate_interior_3prime = 0.0f;
    float cg_stop_rate_baseline_3prime = 0.0f;
    float cg_uniformity_ratio_3prime   = 0.0f;
    bool  channel_g_valid          = false;
    bool  channel_g3_valid         = false;

    // Channel H: A→T stop codon conversion (adenine oxidation / trans-lesion artifact)
    // AAA (Lys) → TAA via A→T at codon pos 1
    // AAG (Lys) → TAG via A→T at codon pos 1
    // AGA (Arg) → TGA via A→T at codon pos 1
    std::array<double, 15> convertible_aaa_h_5prime = {};
    std::array<double, 15> convertible_aag_h_5prime = {};
    std::array<double, 15> convertible_aga_h_5prime = {};
    std::array<double, 15> convertible_taa_at_5prime = {};
    std::array<double, 15> convertible_tag_at_5prime = {};
    std::array<double, 15> convertible_tga_at_5prime = {};
    std::array<double, 15> convertible_aaa_h_3prime = {};
    std::array<double, 15> convertible_aag_h_3prime = {};
    std::array<double, 15> convertible_aga_h_3prime = {};
    std::array<double, 15> convertible_taa_at_3prime = {};
    std::array<double, 15> convertible_tag_at_3prime = {};
    std::array<double, 15> convertible_tga_at_3prime = {};

    float at_stop_rate_terminal    = 0.0f;
    float at_stop_rate_interior    = 0.0f;
    float at_stop_rate_baseline    = 0.0f;
    float channel_h_z              = std::numeric_limits<float>::quiet_NaN();  // C1/C2: NaN = not computed; clamped when computed
    float channel_h_z_p2plus       = std::numeric_limits<float>::quiet_NaN();  // C1/C2: NaN = not computed; clamped when computed
    float at_uniformity_ratio      = 0.0f;
    // C5: h_z and h_z_p2plus can have opposite signs; this flag (same-sign) lets
    // the emitter/consumer see the contradiction. AND-gate detection lives in the emitter.
    bool  channel_h_z_consistent   = false;
    float at_stop_rate_terminal_3prime = 0.0f;
    float at_stop_rate_interior_3prime = 0.0f;
    float at_stop_rate_baseline_3prime = 0.0f;
    float at_uniformity_ratio_3prime   = 0.0f;
    bool  channel_h_valid          = false;
    bool  channel_h3_valid         = false;

    uint64_t ca_pre_interior  = 0;  // Channel F far-interior (pos 30+)
    uint64_t ca_stop_interior = 0;
    uint64_t ca_deam_shadow_interior = 0;  // far-interior C→T shadows (TTA/TTG/TAT/TGT)
    std::array<double, 3> ca_pre_interior_by_ctx    = {};  // [0]=TCA+TAC, [1]=TCG, [2]=TGC
    std::array<double, 3> ca_stop_interior_by_ctx   = {};
    std::array<double, 3> ca_shadow_interior_by_ctx = {};
    uint64_t cg_pre_interior  = 0;  // Channel G far-interior
    uint64_t cg_stop_interior = 0;
    uint64_t at_pre_interior  = 0;  // Channel H far-interior
    uint64_t at_stop_interior = 0;



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
    // C1: ox_gt_uniformity default 0.0f is ambiguous (no companion validity flag).
    bool  ox_gt_uniformity_computed = false;

    float ox_ca_rate_terminal = 0.0f;
    float ox_ca_rate_interior = 0.0f;
    float ox_ca_baseline = 0.0f;
    float ox_ca_uniformity = 0.0f;

    // Channel D Chargaff-contrast validity (C1). True only when both interior
    // baselines reach coverage (gt_base_total>=500 AND ca_base_total>=500); gates
    // D / tg_interior / ac_interior / s_gt to null when false.
    bool  d_computed = false;

    float ox_d_max = 0.0f;
    bool ox_damage_detected = false;
    bool ox_is_artifact = false;
    // C5: ox_damage_detected is written by two disagreeing paths (codon block sets
    // it WITH ox_d_max; GT model-fit overwrites it WITHOUT ox_d_max). Keep both
    // verdicts explicit so channel_c_detected and ox_d_max can stay coherent.
    // ox_d_max remains coupled to the codon path (feeds the preservation gate via
    // ox_is_artifact). ox_damage_detected continues to hold the model-fit verdict.
    bool ox_damage_detected_codon = false;  // codon-block (G→T stop) verdict; coupled to ox_d_max
    bool ox_damage_detected_model = false;  // GT model-fit (s_gt / best_B) verdict

    // Two-marker oxidation bins: s1 (C→T 5' pos 1-3) × s2 (G→A 3' pos 1-3) × GC × length.
    // Used for the deamination-coupled G→T regression (β₁ and β₂ consistency check).
    struct OxoTwoMarkerBins {
        static constexpr int N_S  = 4;   // 0,1,2,3 bases matching marker at pos 1-3
        static constexpr int N_GC = 4;   // GC bins: 0-25%, 25-50%, 50-75%, 75-100%
        static constexpr int N_L  = 4;   // length bins: <45, 45-70, 70-110, 110+
        static constexpr int TOTAL = N_S * N_S * N_GC * N_L;  // 256

        struct Cell {
            uint32_t n_reads = 0;
            uint32_t sum_nGT = 0;  // interior T+G count
            uint32_t sum_T   = 0;  // interior T count  (T→T+G numerator)
            uint32_t sum_nAC = 0;  // interior A+C count
            uint32_t sum_A   = 0;  // interior A count  (A→A+C numerator)
        };

        std::array<Cell, TOTAL> cells = {};

        static int idx(int s1, int s2, int gc, int l) {
            return s1 * (N_S * N_GC * N_L) + s2 * (N_GC * N_L) + gc * N_L + l;
        }
    };
    OxoTwoMarkerBins oxo_two_marker;

    // GT exponential-background fit: GT(p) = A*exp(-mu*p) + B
    // B = uniform background (8-oxoG estimate); A = terminal excess (artifact)
    float g_bg_fitted        = 0.0f;  // B: uniform G→T background
    float g_term_fitted      = 0.0f;  // A: terminal excess
    float g_decay_fitted     = 0.0f;  // mu: terminal decay rate
    float g_bg_fitted_ci_lo  = 0.0f;  // WLS CI95 lower bound on B
    float g_bg_fitted_ci_hi  = 0.0f;  // WLS CI95 upper bound on B
    float g_bg_interior_mean = 0.0f;  // interior-mean G→T (positions 5-14, fallback when degenerate)
    bool  g_fit_degenerate   = false; // true when best_mu <= 0.10 (flat model preferred)
    // Boundary/degeneracy flags for the GT exponential-background fit (C4).
    // Emitter gates g_decay_fitted/g_bg_fitted/ox_theta to null when these fire.
    bool  gt_decay_at_upper_boundary = false;  // best_mu == grid max (true mu unidentifiable)
    bool  gt_term_zero_clamped       = false;  // best_A clamped to 0 -> decay meaningless
    bool  gt_bg_at_upper_boundary    = false;  // B_raw > 0.5 -> background clamped at 0.5
    float g_bg_fitted_unclamped      = 0.0f;   // B_raw point estimate (auditable when clamped)
    bool  ox_theta_at_clamp          = false;  // g_bg_fitted at 0.5 clamp while degenerate
    float s_gt = 0.0f;               // B - ox_ca_baseline: Chargaff contrast (signal for SS; ~0 for DS)

    // Channel E: depurination (purine loss at strand breaks)
    float purine_rate_terminal_5prime = 0.0f;
    float purine_rate_interior = 0.0f;
    float purine_enrichment_5prime = 0.0f;
    float purine_enrichment_3prime = 0.0f;
    bool depurination_detected = false;
    // Channel E validity (C1): true only when total_terminal > 100 AND purine_baseline > 0.01.
    // Emitter nulls rate_interior/enrichment_* and three-states `detected` when false.
    bool channel_e_valid = false;

    // GC-stratified damage: separate estimation per GC bin (interior GC to avoid bias)
    static constexpr int N_GC_BINS = 10;  // 0-10%, 10-20%, ..., 90-100%

    struct GCBinStats {
        // Channel A: T and C counts at 5' terminal positions (0-14)
        std::array<uint64_t, 15> t_counts = {};
        std::array<uint64_t, 15> c_counts = {};
        // Control channel: A and G counts at 5' terminal positions (0-14)
        std::array<uint64_t, 15> a_counts = {};
        std::array<uint64_t, 15> g_counts = {};
        // 3' terminal counts (position 0 = last base, 1 = second-to-last, ...).
        // For SS libraries damage shows as C→T on both ends; for DS it's G→A at 3'.
        std::array<uint64_t, 15> t_counts_3prime = {};
        std::array<uint64_t, 15> c_counts_3prime = {};
        std::array<uint64_t, 15> a_counts_3prime = {};
        std::array<uint64_t, 15> g_counts_3prime = {};
        // Channel A: interior baseline counts
        uint64_t t_interior = 0;
        uint64_t c_interior = 0;
        // Control channel: interior baseline counts
        uint64_t a_interior = 0;
        uint64_t g_interior = 0;

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

    // ── Bulk damage law (Phase 1): per-length-bin terminal counts ─────────────
    // The bulk count-level law needs per-read-length terminal substitution
    // counts at both ends. The pooled *_freq_* arrays and the GC bins carry no
    // length structure, so accumulate fine fixed length bins here during the
    // streaming pass; finalize_sample_profile aggregates them into data-driven
    // adaptive (quantile) length bins for BulkDamageModel::fit.
    static constexpr int LEN_FINE_MIN = 30;   // reads < 30 rejected upstream
    static constexpr int LEN_FINE_W   = 5;    // fine bin width (bp)
    static constexpr int N_LEN_FINE   = 40;   // covers [30, 230); >= 230 -> top bin

    struct LenBinStats {
        std::array<uint64_t, 15> t_counts = {};          // 5' terminal
        std::array<uint64_t, 15> c_counts = {};
        std::array<uint64_t, 15> a_counts = {};
        std::array<uint64_t, 15> g_counts = {};
        std::array<uint64_t, 15> t_counts_3prime = {};   // 3' terminal (i=0 last base)
        std::array<uint64_t, 15> c_counts_3prime = {};
        std::array<uint64_t, 15> a_counts_3prime = {};
        std::array<uint64_t, 15> g_counts_3prime = {};
        uint64_t t_interior = 0, c_interior = 0, a_interior = 0, g_interior = 0;
        uint64_t n_reads = 0;
        uint64_t len_sum = 0;   // sum of read lengths (mean-length ordering key)

        // Per-read eligibility-stratified damage co-occurrence (mixture P2). JStrat keys {n, Σk5, Σk3,
        // Σk5·k3} by the per-read damage-channel site counts (S5,S3) over JP=5 terminal positions
        // (damage-invariant ⇒ conditioning removes the GC confound). ds/ss resolved at aggregation.
        JStrat jstrat_ds, jstrat_ss;
    };
    std::array<LenBinStats, N_LEN_FINE> len_bins = {};

    static int len_fine_bin(size_t len) {
        if (len < static_cast<size_t>(LEN_FINE_MIN)) return 0;
        int b = static_cast<int>((len - static_cast<size_t>(LEN_FINE_MIN)) / LEN_FINE_W);
        return b < N_LEN_FINE ? b : N_LEN_FINE - 1;
    }

    // Aggregated GC-stratified results
    float gc_stratified_d_max_weighted = 0.0f;  // C-site-weighted mean of CT-only b.d_max (Channel A rate; source="average")
    float gc_stratified_d_max_joint = 0.0f;     // C5: C-site-weighted mean of max(Channel A, Channel B) per bin (can exceed d_max_5prime)
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
    float joint_z_delta = 0.0f;        // Wald z for δ̂ (observed Fisher information)
    bool  joint_z_delta_capped = false; // |raw Wald z| exceeded Z_CAP before clamping
    int   joint_n_informative = 0;     // informative 5' positions (cov≥100 & excess > 2·noise)
    bool joint_model_valid = false;    // Sufficient data for joint model

    // ── Bulk damage law (Phase 1) result: threshold-free δ(L) ─────────────────
    // Count-level binomial GLM over data-driven length bins (BulkDamageModel).
    // Emitted as the `bulk_damage` JSON block; not yet wired into d_max_combined.
    BulkDamageResult bulk_damage;
    double bulk_headline_delta = 0.0;  // read-weighted δ̂ over identified length bins
    bool   bulk_attempted = false;     // fit ran (>= 1 valid length bin)

    // Mixture model results (K-component EM over GC-stratified bins)
    int mixture_K = 0;                 // Number of classes selected by BIC
    int mixture_n_components = 0;      // Number of classes selected by BIC
    float mixture_d_population = 0.0f; // E[δ] over all C-sites
    float mixture_d_ancient = 0.0f;    // E[δ | δ > 5%] (ancient tail)
    float mixture_d_population_highgc = 0.0f;  // E[δ | GC >= 50%] — high-GC-weighted mean (diagnostic only)
    float mixture_pi_ancient = 0.0f;   // Fraction of C-sites in high-damage classes
    float mixture_bic = 0.0f;          // BIC for model selection
    bool mixture_converged = false;    // Did EM converge?
    bool mixture_identifiable = false; // Are non-trivial classes well-separated?

    // Preservation index
    enum class PreservationLabel {
        INSUFFICIENT,       // too few reads for reliable estimate
        ARTIFACT_SUSPECTED, // ox_is_artifact=true
        MODERN_LIKE,        // score < 0.15
        WEAK,               // 0.15–0.35
        MODERATE,           // 0.35–0.60
        STRONG,             // 0.60–0.80
        EXCEPTIONAL         // ≥ 0.80
    };
    float preservation_f5          = 0.0f;  // 5' terminal damage factor
    float preservation_f3          = 0.0f;  // 3' terminal damage factor
    float preservation_f_coh       = 0.0f;  // mixture coherence factor
    float preservation_f_cpg       = std::numeric_limits<float>::quiet_NaN();  // D21: NaN default so an early exit does not feed a wrong 0.0 into the geo-mean
    bool  cpg_ratio_evaluable      = false;  // D21: true when f_cpg came from the computed CpG-ratio branch, false when the 0.3 uninformative prior was used
    float preservation_evidence    = 0.0f;  // geometric mean of factors
    float preservation_reliability = 0.0f;  // g_N * g_fit * g_ox
    float preservation_score       = 0.0f;  // evidence * reliability [0,1]
    PreservationLabel preservation_label = PreservationLabel::INSUFFICIENT;

    const char* preservation_label_str() const {
        switch (preservation_label) {
            case PreservationLabel::INSUFFICIENT:       return "insufficient";
            case PreservationLabel::ARTIFACT_SUSPECTED: return "artifact-suspected";
            case PreservationLabel::MODERN_LIKE:        return "modern-like";
            case PreservationLabel::WEAK:               return "weak";
            case PreservationLabel::MODERATE:           return "moderate";
            case PreservationLabel::STRONG:             return "strong";
            case PreservationLabel::EXCEPTIONAL:        return "exceptional";
            default:                                    return "unknown";
        }
    }

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

    // Sentinel returned by try_get_gc_bin() when the sequence is too short or
    // entirely ambiguous (no ACGT bases in the interior). Callers that ignore
    // validity can use the bin-4 fallback via get_gc_bin().
    static constexpr int GC_BIN_UNKNOWN = -1;

    // Compute the interior GC bin. Returns GC_BIN_UNKNOWN when the sequence is
    // shorter than 20 bp or has no ACGT bases in the interior window — letting
    // callers distinguish "GC ≈ middle" from "no usable signal".
    static int try_get_gc_bin(const std::string& seq) {
        if (seq.length() < 20) return GC_BIN_UNKNOWN;
        size_t gc = 0, total = 0;
        size_t start = std::min(size_t(5), seq.length() / 4);
        size_t end = seq.length() - start;

        for (size_t i = start; i < end; ++i) {
            char c = seq[i];
            if (c == 'G' || c == 'g' || c == 'C' || c == 'c') ++gc;
            if (c == 'A' || c == 'a' || c == 'T' || c == 't' ||
                c == 'G' || c == 'g' || c == 'C' || c == 'c') ++total;
        }

        if (total == 0) return GC_BIN_UNKNOWN;
        float gc_frac = static_cast<float>(gc) / static_cast<float>(total);
        int bin = static_cast<int>(gc_frac * 10.0f);
        return std::clamp(bin, 0, N_GC_BINS - 1);
    }

    // Backwards-compatible wrapper: collapses the unknown case to bin 4 (the
    // historical silent fallback). Prefer try_get_gc_bin() in new code.
    static int get_gc_bin(const std::string& seq) {
        int b = try_get_gc_bin(seq);
        return b == GC_BIN_UNKNOWN ? 4 : b;
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
        int bin = try_get_gc_bin(seq);
        GCDamageParams params{};
        if (bin == GC_BIN_UNKNOWN) {
            params.bin_valid = false;
            return params;
        }
        const auto& stats = gc_bins[bin];
        params.p_damaged   = stats.p_damaged;
        params.delta_s     = stats.d_max;
        params.baseline_tc = stats.baseline_tc;
        params.lambda      = lambda_5prime;
        params.bin_valid   = stats.valid;
        return params;
    }

    // Get effective damage rate for a read (posterior-weighted)
    // δ_eff = P(damaged|read) * δ_s. Returns 0 when the GC bin is unknown.
    float get_effective_damage(const std::string& seq, float read_ancient_prob) const {
        int bin = try_get_gc_bin(seq);
        if (bin == GC_BIN_UNKNOWN) return 0.0f;
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

} // namespace taph
