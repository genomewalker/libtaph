#pragma once

// Library damage profiler: the pipeline behind `fqdup profile` and `dart profile`.
// Both tools only parse arguments and supply a SeqReaderOpener; everything that
// decides a damage number lives here, so the two produce byte-identical output.

#include "taph/damage_profile.hpp"
#include "taph/seq_reader.hpp"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace taph {

// --length-bins spec: auto | off | N | e1,e2,... ; ok is false for a bad spec.
LengthBinOptions parse_length_bins(const std::string& spec, bool& ok);

struct LibraryProfileOptions {
    std::string in_path;                       // single-end / merged FASTQ
    std::string r1_path, r2_path;              // paired raw reads (R1 5', R2 3')
    int         n_threads = 1;
    SampleDamageProfile::LibraryType forced_lib = SampleDamageProfile::LibraryType::UNKNOWN;
    double      mask_threshold = 0.05;
    std::string tsv_path, json_path, damage_summary_path, count_table_path, html_path;
    int         json_level = 2;                // 1 summary, 2 standard, 3 full
    std::vector<std::string> subst_in_paths;   // fqdup merge .bsubst files
    bool        run_oxog = true;
    int         min_length = 0;                // <30 lowers the short-fragment floor to max(16, N)
    int         max_length = 0;                // 0 = off
    LengthBinOptions lb_opts = [] { bool ok; return parse_length_bins("auto", ok); }();
    int64_t     adapter_scan_reads = 1'000'000;  // 0 = scan all reads
    std::string producer;                      // JSON version field, e.g. "fqdup/<sha>"
    std::function<void(const std::string&)> log_info, log_warn;
};

// Handed to the caller at the point `fqdup profile` writes its split-model export.
struct ProfileModelExportInput {
    bool                         paired_mode;
    const std::string&           in_path;
    const BulkDamageProfile&         bulk_dp;
    int64_t                      reads_scanned;
    const std::vector<uint64_t>& merged_lsd_hist;
};
// Returns nonzero to abort the run.
using ProfileModelExport = std::function<int(const ProfileModelExportInput&)>;

// Adapter pre-scan, one pass over the reads, finalize, length-stratified and OxoG
// channels, then the requested outputs (JSON, .tdmg, TSV, count tables, HTML) and the
// text report on stdout. Returns 0 on success, 1 on an input or output error.
int run_library_profile(const LibraryProfileOptions& opt, const SeqReaderOpener& open_reader,
                        const ProfileModelExport& export_model = {});

} // namespace taph
