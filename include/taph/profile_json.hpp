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
    double d_oriented   = 0.0;
    bool   has_oxog_score = false;

    // Length-stratified profiles (optional; full per-bin summary from fqdup)
    const LengthStratifiedDamageProfile* lsd = nullptr;
};

// Serialize a finalized SampleDamageProfile to JSON on `out`.
// Computes all derived scores (CpG, hex, OxoG trinuc, depurination, preservation)
// internally. Callers supply only ProfileJsonInput for programmatic context.
// Output format is compatible with fqdup's JSON schema.
void profile_to_json(const SampleDamageProfile& dp,
                     std::ostream& out,
                     const ProfileJsonInput& in = {});

} // namespace taph
