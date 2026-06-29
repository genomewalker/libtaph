#pragma once

#include "sample_damage_profile.hpp"
#include <string>
#include <string_view>
#include <vector>

namespace taph {

class FrameSelector {
public:
    static SampleDamageProfile compute_sample_profile(
        const std::vector<std::string>& sequences);

    static void update_sample_profile(
        SampleDamageProfile& profile,
        std::string_view seq);

    // Paired-end variant: r1 contributes to molecule 5'-end counters,
    // r2 contributes to molecule 3'-end counters via complement-mapping
    // (R2 reads the bottom strand from the molecule 3' inward, so each
    // R2[i] is complement(top_strand[len-1-i])). Read 3' ends ignored
    // (adapter contamination zone for inserts < read length).
    // Returns true if pair was scored, false if skipped (short or filtered).
    static bool update_sample_profile_paired(
        SampleDamageProfile& profile,
        std::string_view r1,
        std::string_view r2);

    static void update_sample_profile_weighted(
        SampleDamageProfile& profile,
        std::string_view seq,
        float weight);

    // bulk_fit_threads parallelizes ONLY the main bulk MLE's per-bin profile loop (BulkDamageModel::fit).
    // Default 1 keeps every caller serial; the length-strat per-bin path leaves it at 1 (its bins are
    // already distributed across a worker pool — passing >1 here would nest-oversubscribe).
    static void finalize_sample_profile(SampleDamageProfile& profile, int bulk_fit_threads = 1);

    static void merge_sample_profiles(SampleDamageProfile& dst, const SampleDamageProfile& src);

    static void reset_sample_profile(SampleDamageProfile& profile);

    // After finalize_sample_profile and detect_adapter_stubs, recompute channel
    // F/G/H terminal sums and z-scores excluding reads whose 5' (or 3') first
    // hexamer is in the supplied adapter-stub code sets.  Updates channel_f/g/h_z
    // and associated rate fields in-place; sets fgh_adapter_prefixes_excluded.
    static void recompute_fgh_excluding_adapter_prefixes(
        SampleDamageProfile& profile,
        const std::vector<uint32_t>& excl_5p,
        const std::vector<uint32_t>& excl_3p);

    static int get_gc_bin(const std::string& seq) {
        return SampleDamageProfile::get_gc_bin(seq);
    }
};

} // namespace taph
