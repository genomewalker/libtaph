#pragma once
#include "taph/sample_damage_profile.hpp"
#include "taph/library_interpretation.hpp"
#include "taph/length_bin_damage_profile.hpp"
#include <iosfwd>
#include <string>
#include <vector>

namespace taph {

// External/programmatic context that cannot be derived from SampleDamageProfile alone.
// All fields are optional with sensible defaults; leave unset for minimal output.
struct ProfileJsonInput {
    std::string sample_name;          // input path / sample identifier
    std::string version;              // program version string
    uint64_t    n_reads      = 0;

    // Pre-rendered JSON fragment (e.g. `"cooccur53": {...}`, no leading/trailing
    // comma) spliced verbatim before the root close. Empty = omitted.
    std::string cooccur53_json;

    // Pre-rendered `"library": {...}` fragment from fqdup merge's .bsubst v3 trailer
    // (library_type, biological_termini, linker, adapters, merge counts, balance).
    // No leading/trailing comma; spliced before the root close. Empty = omitted.
    std::string library_json;

    // Adapter detection results (optional; leave default when not performed)
    std::vector<std::string>  adapter_stubs_5prime;
    std::vector<std::string>  adapter_stubs_3prime;
    std::vector<HexEnrichment> top_hex_enriched;
    std::vector<HexEnrichment> top_hex_enriched_3prime;
    HexEndAsymmetry            hex_end_asymmetry;
    bool adapter_clipped    = false;
    bool adapter3_clipped   = false;
    bool flag_hex_artifact  = false;

    // Adapter stub read fraction (optional; -1 = not computed)
    // Fraction of pre-scan reads whose terminal 6 bp matched a detected stub.
    double  adapter_stub5_read_fraction = -1.0;
    double  adapter_stub3_read_fraction = -1.0;
    int64_t adapter_stub_reads_checked  = 0;

    // Fraction of reads shorter than 50 bp (optional)
    double short_read_frac = 0.0;

    // OxoG EM-based score (optional; fqdup computes this from joint EM state)
    double s_oxog       = 0.0;
    double se_s_oxog    = 0.0;
    double s_ca         = 0.0;   // C→A co-movement score (complement of G→T oxidation)
    double d_oriented   = 0.0;
    bool   has_oxog_score = false;

    // Length-stratified profiles (optional; full per-bin summary from fqdup)
    const LengthStratifiedDamageProfile* lsd = nullptr;

    // Per-position overlap substitution rates from fqdup merge --subst-out.
    // Computed from R1 vs RC(R2) mismatches before consensus — strand-resolved.
    // Empty = not provided.
    std::vector<double> paired_ct_decay;   // C→T rate at pos 0..N from 5' of insert
    std::vector<double> paired_ga_decay;   // G→A rate at pos 0..N from 3' of insert
    // OxoG: G→T (fwd[*][T][G]) vs Chargaff control C→A (fwd[*][C][A]).
    // PCR fixation makes the ct/ga decay channels unreliable for authentic aDNA deamination;
    // OxoG is detectable because 8-oxoG persists into sequencing and causes misreading at
    // the base-caller level (not fixed concordantly by PCR).
    int64_t paired_tg_count  = 0;   // sum_pos fwd[pos][T][G]
    int64_t paired_tg_denom  = 0;   // sum_pos sum_r1 fwd[pos][r1][G]
    int64_t paired_ca_count  = 0;   // sum_pos fwd[pos][C][A]  (Chargaff control)
    int64_t paired_ca_denom  = 0;   // sum_pos sum_r1 fwd[pos][r1][A]
    double paired_oxog_rate  = -1.0; // = paired_tg_count / paired_tg_denom (-1 = not computed)
    int64_t paired_n_pairs   = 0;
    int64_t paired_n_bases   = 0;
    std::vector<double> paired_tg_pos;
    std::vector<double> paired_ca_pos;

    // Fitted bsubst decay: bg + d_max * exp(-lambda * p), per-C-site scale (bg≈0.005).
    // -1.0 = not computed (e.g. no bsubst, or <3 non-zero positions).
    double paired_ct_d_max_5prime  = -1.0;
    double paired_ct_lambda_5prime = -1.0;
    double paired_ct_bg_5prime     = -1.0;
    double paired_ga_d_max_3prime  = -1.0;
    double paired_ga_lambda_3prime = -1.0;
    double paired_ga_bg_3prime     = -1.0;
};

// Serialize a finalized SampleDamageProfile to JSON on `out`.
// Computes all derived scores (CpG, hex, OxoG trinuc, depurination, preservation)
// internally. Callers supply only ProfileJsonInput for programmatic context.
// Output format is compatible with fqdup's JSON schema.
void profile_to_json(const SampleDamageProfile& dp,
                     std::ostream& out,
                     const ProfileJsonInput& in = {});

// Serialize the Layer-0 stop-channel count tables to JSONL (one object per line) — the pre-clamp
// source of truth the golden regression gate diffs. Deterministic field order, round-trippable
// precision (setprecision(17)); non-finite -> null. Emits nothing when there are no count tables.
void count_tables_to_jsonl(const SampleDamageProfile& dp, std::ostream& out);

} // namespace taph
