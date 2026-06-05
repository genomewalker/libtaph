#pragma once
// Layer-0 count-table contract — the canonical PRE-CLAMP source of truth for a stop channel.
//
// The emitted-JSON regression gate is cap-masked: a stop channel's z rails at +/-kZCap, so a
// refactor that corrupts the underlying z still renders the identical clamped string and the gate
// sees nothing. The count table is therefore the real source of truth: every emitted stop-channel
// metric must be reproducible from these fields alone, and the golden gate diffs THESE (pre-clamp,
// pre-format) so denominator-membership swaps and cap-masked z regressions cannot hide.
#include <array>
#include <string>

namespace taph {

// One Mantel-Haenszel context stratum (Channel F only): the 2x2 terminal-vs-interior table.
struct StopChannelStratum {
    int    stratum_id  = -1;
    double a_term_stop  = 0.0;  // terminal damage (stop) count
    double b_term_other = 0.0;  // terminal non-stop count (pre + shadow)
    double c_int_stop   = 0.0;  // interior damage (stop) count
    double d_int_other  = 0.0;  // interior non-stop count (pre + shadow)
};

// Canonical per-stop-channel count table. Shadow-free 2x2 building blocks plus the exact
// numerator/denominator that fed the pooled binomial z, the pre-clamp z, and the cap decision.
struct StopChannelCountTable {
    std::string channel_id;          // "F","F3","G","G3","H","H3"
    char        channel_type = '?';  // 'F' | 'G' | 'H'

    // shadow-free terminal/interior building blocks
    double term_pre  = 0.0;
    double term_stop = 0.0;
    double int_pre   = 0.0;
    double int_stop  = 0.0;

    // Channel F deamination shadow: in the z denominator, NOT the rate denominator
    bool   has_shadow     = false;
    double term_shadow    = 0.0;
    double int_shadow     = 0.0;
    bool   shadow_in_z    = false;
    bool   shadow_in_rate = false;

    // resolved numerator/denominator actually used for the pooled z (shadow folded in for F)
    double z_num_term = 0.0;   // = term_stop
    double z_den_term = 0.0;   // = term_pre + term_stop (+ term_shadow if shadow_in_z)
    double z_num_int  = 0.0;   // = int_stop
    double z_den_int  = 0.0;   // = int_pre + int_stop (+ int_shadow if shadow_in_z)

    // shadow-free rates (what the emitted ca_stop_rate_* compare)
    double raw_rate_term = 0.0;
    double raw_rate_int  = 0.0;

    // the statistic: pre_clamp_z is the gate's source of truth; post_clamp_z is what gets emitted
    double pre_clamp_z  = 0.0;
    double z_cap        = 0.0;
    bool   cap_applied  = false;
    double post_clamp_z = 0.0;

    // Channel F context strata (lets the gate recompute the Mantel-Haenszel OR from counts)
    bool                              has_strata = false;
    std::array<StopChannelStratum, 3> strata{};
};

}  // namespace taph
