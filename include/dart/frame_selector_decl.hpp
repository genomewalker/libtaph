#pragma once

#include "sample_damage_profile.hpp"
#include <string>
#include <vector>

namespace dart {

class FrameSelector {
public:
    /**
     * Compute sample-level damage profile from a collection of reads
     * This should be called once on the first pass through the data
     *
     * @param sequences Vector of DNA sequences
     * @return Computed SampleDamageProfile with aggregate statistics
     */
    static SampleDamageProfile compute_sample_profile(
        const std::vector<std::string>& sequences);

    /**
     * Update sample profile incrementally with a new read
     * Thread-safe for parallel processing
     *
     * @param profile Profile to update (should be protected by mutex)
     * @param seq New sequence to add
     */
    static void update_sample_profile(
        SampleDamageProfile& profile,
        const std::string& seq);

    /**
     * Update sample profile with a weighted contribution
     * Used for iterative damage refinement - weights by coding probability
     *
     * @param profile Profile to update
     * @param seq DNA sequence to add
     * @param weight Weight for this sequence (0.0-1.0, typically coding_prob)
     */
    static void update_sample_profile_weighted(
        SampleDamageProfile& profile,
        const std::string& seq,
        float weight);

    /**
     * Finalize sample profile after all reads have been added
     * Computes damage rates from raw frequencies
     *
     * @param profile Profile to finalize
     */
    static void finalize_sample_profile(SampleDamageProfile& profile);

    /**
     * Merge two sample profiles (for parallel aggregation)
     * Adds counts from src into dst
     *
     * @param dst Destination profile (modified in place)
     * @param src Source profile to merge from
     */
    static void merge_sample_profiles(SampleDamageProfile& dst, const SampleDamageProfile& src);

    /**
     * Reset sample profile for a fresh damage estimation pass
     * Clears all counts while preserving structure
     *
     * @param profile Profile to reset
     */
    static void reset_sample_profile(SampleDamageProfile& profile);

    /**
     * Get GC bin index (0-9) for a sequence based on interior GC content
     * Delegates to SampleDamageProfile::get_gc_bin
     */
    static int get_gc_bin(const std::string& seq) {
        return SampleDamageProfile::get_gc_bin(seq);
    }
};

} // namespace dart
