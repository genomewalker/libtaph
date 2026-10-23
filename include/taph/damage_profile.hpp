#pragma once

// BulkDamageProfile (fqdup: DamageProfile): bulk exponential damage model plus a per-position empirical mask,
// as fitted by the library profiler. Moved from fqdup so the profiler can live in
// libtaph; fqdup's derep/extend/split use it through fqdup/damage_profile.hpp.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "taph/library_interpretation.hpp"
#include "taph/length_bin_damage_profile.hpp"
#include "taph/length_stratified_profile.hpp"
#include "taph/sample_damage_profile.hpp"
#include "taph/lsd_accumulator.hpp"
#include "taph/seq_reader.hpp"

namespace taph {

// Oxidative channel detection threshold (binomial z-score)
inline constexpr float kOxChannelZDetect = 3.0f;


inline int lsd_hist_bin(int L) {
    static const double kLogMin = std::log(static_cast<double>(LSD_L_MIN));
    static const double kStep   = (std::log(static_cast<double>(LSD_L_MAX + 1)) - kLogMin)
                                  / LSD_HIST_BINS;
    int Lcap = L < LSD_L_MAX ? L : LSD_L_MAX;
    int hb   = static_cast<int>((std::log(static_cast<double>(Lcap)) - kLogMin) / kStep);
    if (hb < 0)              hb = 0;
    if (hb >= LSD_HIST_BINS) hb = LSD_HIST_BINS - 1;
    return hb;
}

// 3' damage prior mode. ds_ga = use 3' G→A signal (DS libraries),
// ss_ct = use 3' C→T signal (SS libraries), neutral = no 3' damage prior
// (UNKNOWN library type — drop the 3' channel from masking entirely).
enum class Damage3PrimeMode : uint8_t { ds_ga, ss_ct, neutral };

struct BulkDamageProfile {
    // 3' damage prior mode — set from taph::SampleDamageProfile::library_type. UNKNOWN maps to
    // neutral so we don't silently apply DS or SS priors. ss_mode (below) is
    // kept for compatibility and derived from this field.
    Damage3PrimeMode mode_3prime = Damage3PrimeMode::ds_ga;

    // Fitted exponential model parameters (for logging and expected_mismatches only).
    double d_max_5prime   = 0.0;
    double d_max_3prime   = 0.0;
    double lambda_5prime  = 0.5;
    double lambda_3prime  = 0.5;
    double background     = 0.02;
    // Interior baseline fractions used by per-read LLR (null-distribution models).
    // bg_5_tc = 5' interior T/(T+C); bg_3_channel = 3' interior T/(T+C) for SS
    // libraries, A/(A+G) for DS libraries.
    double bg_5_tc        = 0.5;
    double bg_3_channel   = 0.5;
    // Bulk CpG-split 5' C→T d_max (from taph::dmax_ct5_cpg_like / _noncpg_like).
    // Negative (NaN) = not available. Used to weight the per-read 5' LLR: CpG
    // C/T sites are stronger evidence of ancient than non-CpG sites.
    double d_cpg_5prime    = -1.0;
    double d_noncpg_5prime = -1.0;
    double mask_threshold = 0.05;
    double pcr_error_rate = 0.0;
    double mixture_d_damaged = 0.0;
    double mixture_pi_damaged = 0.0;
    double mixture_d_population_highgc = 0.0;
    int    mixture_K = 0;
    int    mixture_n_components = 0;
    bool   mixture_converged = false;
    bool   mixture_identifiable = false;
    // Validated reference-free ancient fraction (libtaph SampleDamageProfile::pi, SOLUTION §6.3).
    // Shadow-mode step 2: carried for the contract-vs-mixture divergence shadow; not yet a verdict source.
    double pi_point    = -1.0;
    double pi_lo       = -1.0;
    double pi_hi       = -1.0;
    taph::DamageConfidence pi_state = taph::DamageConfidence::UNDETERMINED;
    bool   pi_detected = false;
    bool   enabled        = false;
    bool   ss_mode        = false;  // true = single-stranded library (C→T at both ends)

    // Empirical per-position mask (primary masking authority).
    //
    // mask_pos[p] = true means position p (measured from either terminus) sits
    // in the damage zone and will be replaced with a neutral byte before hashing.
    //
    // Symmetry invariant: canonical_hash(seq) == canonical_hash(revcomp(seq)).
    // apply_damage_mask() uses max(mask_pos[i], mask_pos[j]) semantics so that
    // position p from the 5' end and position p from the 3' end are treated
    // identically regardless of which strand is processed.
    // Upper bound on terminal positions tracked for damage masking.
    // 32 covers extreme aDNA scenarios (very heavy deamination > 20 bp deep)
    // while staying tiny in memory (~32 B per profile). All loops use this
    // constant directly — bump and recompile to extend further.
    static constexpr int MASK_POSITIONS = 32;
    bool mask_pos[MASK_POSITIONS] = {};

    // Populate mask_pos from the fitted exponential model.
    // Used when damage parameters are given manually rather than estimated.
    void populate_mask_from_model() {
        for (int p = 0; p < MASK_POSITIONS; ++p) {
            double exc5 = d_max_5prime * std::exp(-lambda_5prime * p);
            double exc3 = d_max_3prime * std::exp(-lambda_3prime * p);
            mask_pos[p] = std::max(exc5, exc3) > mask_threshold;
        }
    }

    // Helpers retained for expected_mismatches() and informational logging only.
    double dmg_5(int p) const { return d_max_5prime * std::exp(-lambda_5prime * p); }
    double dmg_3(int p) const { return d_max_3prime * std::exp(-lambda_3prime * p); }

    // P(damage-driven mismatch at pos) for a read of length L. Used by Phase 3
    // LR scoring (T5.2) to weigh damage-channel mismatches. ss_mode applies
    // C→T at both ends; ds_mode applies C→T at 5' and G→A at 3'. Caller is
    // responsible for checking that the mismatch alt is in the damage channel.
    double p_damage_at(int pos, int L) const {
        if (!enabled || L <= 0) return 0.0;
        double d5 = dmg_5(pos);
        double d3 = dmg_3(L - 1 - pos);
        return std::max(d5, d3);  // either end can drive damage at this position
    }

    double expected_mismatches(int L) const {
        // d_max and the exponential model are defined on the conditional T/(T+C)
        // scale: dmg_5(p) is the excess deamination rate among C/T bases only.
        // The per-position mismatch probability between two reads is therefore
        // dmg × (1-dmg) × f_C, where f_C = 1-background is the C fraction within
        // the T+C pool (background = T/(T+C) at undamaged positions). Without this
        // factor the expectation is inflated by ~1/(1-background) ≈ 2× for typical
        // aDNA base composition. The same factor applies to the 3' G→A channel by
        // complement symmetry (G fraction in A+G pool ≈ 1-background).
        const double f_suscep = 1.0 - background;
        double e = 0.0;
        for (int p = 0; p < L; ++p) {
            double d5 = dmg_5(p);
            double d3 = dmg_3(L - 1 - p);
            e += 2.0 * d5 * (1.0 - d5) * f_suscep;
            e += 2.0 * d3 * (1.0 - d3) * f_suscep;
        }
        e += 2.0 * L * pcr_error_rate;
        return e;
    }

    int mismatch_tolerance(int L) const {
        double lam = expected_mismatches(L);
        if (lam < 1e-9) return 0;
        double cumP = 0.0, pois = std::exp(-lam);
        int k = 0;
        while (cumP + pois < 0.99 && k < 100) {
            cumP += pois;
            ++k;
            pois *= lam / k;
        }
        return k;
    }

    void print_info(int typical_read_length,
                    const std::function<void(const std::string&)>& log_info) const {
        log_info("--- Damage-Aware Deduplication ---");
        log_info("  5'-end d_max:  " + std::to_string(d_max_5prime));
        log_info("  3'-end d_max:  " + std::to_string(d_max_3prime));
        log_info("  5'-end lambda: " + std::to_string(lambda_5prime));
        log_info("  3'-end lambda: " + std::to_string(lambda_3prime));
        log_info("  Background:    " + std::to_string(background));
        log_info("  Mask threshold:" + std::to_string(mask_threshold));
        if (pcr_error_rate > 0.0)
            log_info("  PCR error rate:" + std::to_string(pcr_error_rate) + " per bp");
        int n_masked = 0;
        std::string pos_str;
        for (int p = 0; p < MASK_POSITIONS; ++p) {
            if (mask_pos[p]) {
                if (n_masked > 0) pos_str += ',';
                pos_str += std::to_string(p);
                ++n_masked;
            }
        }
        if (n_masked > 0)
            log_info("  Masked positions: " + pos_str + " (" +
                     std::to_string(n_masked) + " bp each end)");
        else
            log_info("  Masked positions: none");
        if (typical_read_length > 0) {
            double e   = expected_mismatches(typical_read_length);
            int    tol = mismatch_tolerance(typical_read_length);
            log_info("  Expected mismatches (L=" + std::to_string(typical_read_length) +
                     "): " + std::to_string(e) +
                     ", 99th-pct tolerance: " + std::to_string(tol));
        }
    }
};


// Build classify params from a finalized bulk BulkDamageProfile. Thin wrapper:
// extracts the scalar fields taph::make_lsd_classify_params needs into the
// tool-independent taph::BulkDamageForClassify view, so DART (or any other
// consumer) can call the shared libtaph function directly with its own
// equivalent bulk-fit type instead of needing fqdup's BulkDamageProfile.
inline LsdClassifyParams make_lsd_classify_params(const BulkDamageProfile& bulk) {
    taph::BulkDamageForClassify b;
    b.mixture_converged = bulk.mixture_converged;
    b.mixture_d_damaged  = bulk.mixture_d_damaged;
    b.d_max_5prime       = bulk.d_max_5prime;
    b.d_max_3prime       = bulk.d_max_3prime;
    b.lambda_5prime      = bulk.lambda_5prime;
    b.lambda_3prime      = bulk.lambda_3prime;
    b.bg_5_tc            = bulk.bg_5_tc;
    b.bg_3_channel       = bulk.bg_3_channel;
    b.ss_mode            = bulk.ss_mode;
    b.d_cpg_5prime       = bulk.d_cpg_5prime;
    b.d_noncpg_5prime    = bulk.d_noncpg_5prime;
    b.pi_detected        = bulk.pi_detected;
    return taph::make_lsd_classify_params(b);
}

// Prebuilt LSD data produced by a fused oxoG+LSD pass.  Pass via the
// prebuilt parameter of estimate_damage_by_length to skip the FASTQ reader.
struct LsdPrebuilt {
    std::vector<int>            edges;         // from compute_lsd_edges
    taph::LengthBinStats        merged_stats;  // accumulated but NOT finalized
    std::vector<LsdLlrBinAccum> llr_bins;      // per-bin, merged across threads
};

// Length-stratified damage. Pass 1 builds a log-length histogram and picks bin edges
// (explicit | quantile | auto GMM+BIC); pass 2 routes each read to its per-bin profile.
// With `prebuilt` set (the fused pass of the library profiler) neither pass reads `path`.
LengthStratifiedDamageProfile estimate_damage_by_length(
    const std::string& path,
    const SeqReaderOpener& open_reader,
    SampleDamageProfile::LibraryType forced_lib,
    const LengthBinOptions& options,
    const std::vector<uint64_t>* prebuilt_hist = nullptr,
    size_t reader_threads = 0,
    int64_t max_reads = 0,
    const BulkDamageProfile* bulk_for_llr = nullptr,
    const LsdPrebuilt* prebuilt = nullptr,
    const std::function<void(const std::string&)>& log_info = {},
    const std::function<void(const std::string&)>& log_warn = {});

} // namespace taph
