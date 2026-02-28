#pragma once

#include "sample_damage_profile.hpp"
#include <string>
#include <vector>

namespace dart {

class FrameSelector {
public:
    static SampleDamageProfile compute_sample_profile(
        const std::vector<std::string>& sequences);

    static void update_sample_profile(
        SampleDamageProfile& profile,
        const std::string& seq);

    static void update_sample_profile_weighted(
        SampleDamageProfile& profile,
        const std::string& seq,
        float weight);

    static void finalize_sample_profile(SampleDamageProfile& profile);

    static void merge_sample_profiles(SampleDamageProfile& dst, const SampleDamageProfile& src);

    static void reset_sample_profile(SampleDamageProfile& profile);

    static int get_gc_bin(const std::string& seq) {
        return SampleDamageProfile::get_gc_bin(seq);
    }
};

} // namespace dart
