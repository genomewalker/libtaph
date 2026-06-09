// FrameSelector::recompute_fgh_excluding_adapter_prefixes.
#include "taph/frame_selector_decl.hpp"
#include "taph/channel_count_table.hpp"
#include "taph/channel_registry.hpp"
#include "damage_estimation_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
namespace taph {
void FrameSelector::recompute_fgh_excluding_adapter_prefixes(
    SampleDamageProfile& profile,
    const std::vector<uint32_t>& excl_5p,
    const std::vector<uint32_t>& excl_3p)
{
    if (excl_5p.empty() && excl_3p.empty()) return;

    bool skip5[4096] = {}, skip3[4096] = {};
    for (auto c : excl_5p) if (c < 4096) skip5[c] = true;
    for (auto c : excl_3p) if (c < 4096) skip3[c] = true;

    // Helper: resum prefix-binned arrays excluding flagged codes, then recompute
    // binom_z and fin_oxog-style rate fields.
    auto resum = [](const std::array<double,4096>& pre_arr,
                    const std::array<double,4096>& stop_arr,
                    const bool skip[4096]) -> std::pair<double,double> {
        double pre = 0, stop = 0;
        for (uint32_t i = 0; i < 4096; ++i) {
            if (!skip[i]) { pre += pre_arr[i]; stop += stop_arr[i]; }
        }
        return {pre, stop};
    };

    // Interior baselines are not prefix-split; reuse existing values.
    // Channel F/G/H z is recomputed from the adapter-excluded terminal counts vs unchanged interior
    // via the single binom_z_clamped (the duplicated binom_z_f lambda was deleted in P3).

    uint32_t n_excl = 0;

    // Single-producer invariant across the exclusion path: the emitted z/rate are overwritten below
    // from a Layer-0 row rebuilt on the adapter-excluded counts, and that row REPLACES the primary
    // row in profile.count_tables, so the golden gate sees exactly what is emitted. The
    // Mantel-Haenszel strata are not adapter-excluded (channel_f_mh_z is a separate stratified
    // statistic computed on the full window and not recomputed here), so preserve them from the
    // primary row.
    auto replace_ct_row = [&profile](char type, StopChannelCountTable nct) {
        for (auto& row : profile.count_tables) {
            if (row.channel_type == type) {
                nct.strata = row.strata;
                nct.has_strata = row.has_strata;
                row = std::move(nct);
                return;
            }
        }
    };

    if (!excl_5p.empty()) {
        // Channel F 5' — rebuild the Layer-0 row from the adapter-excluded terminal counts (interior
        // is not prefix-split), then take z and shadow-free rates straight from that row.
        auto [pf, sf] = resum(profile.ca_pre_terminal_by_pfx,
                               profile.ca_stop_terminal_by_pfx, skip5);
        double shadow_f = 0;
        for (uint32_t i = 0; i < 4096; ++i)
            if (!skip5[i]) shadow_f += profile.ca_deam_shadow_terminal_by_pfx[i];
        StopChannelCountTable f_ct = make_stop_count_table(stop_channel_spec('F'),
            pf, sf, profile.ca_pre_interior, profile.ca_stop_interior,
            shadow_f, profile.ca_deam_shadow_interior);
        replace_ct_row('F', f_ct);
        profile.channel_f_z = static_cast<float>(f_ct.post_clamp_z);
        profile.ca_stop_rate_terminal = static_cast<float>(f_ct.raw_rate_term);  // shadow-free (registry shadow_in_rate=false), matches primary path
        profile.ca_stop_rate_interior = static_cast<float>(f_ct.raw_rate_int);
        if (profile.ca_stop_rate_interior > 1e-6f && profile.ca_stop_rate_terminal > 0.0f)
            profile.ca_uniformity_ratio = profile.ca_stop_rate_terminal / profile.ca_stop_rate_interior;
        profile.channel_f_valid = (f_ct.z_den_term >= 10 && f_ct.z_den_int >= 10);
        // Channel G 5'
        auto [pg, sg] = resum(profile.cg_pre_terminal_by_pfx,
                               profile.cg_stop_terminal_by_pfx, skip5);
        StopChannelCountTable g_ct = make_stop_count_table(stop_channel_spec('G'),
            pg, sg, profile.cg_pre_interior, profile.cg_stop_interior, 0.0, 0.0);
        replace_ct_row('G', g_ct);
        profile.channel_g_z = static_cast<float>(g_ct.post_clamp_z);
        profile.channel_g_or = static_cast<float>(stop_channel_or_haldane(g_ct));
        profile.cg_stop_rate_terminal = static_cast<float>(g_ct.raw_rate_term);
        profile.cg_stop_rate_interior = static_cast<float>(g_ct.raw_rate_int);
        profile.channel_g_valid = (g_ct.z_den_term >= 10 && g_ct.z_den_int >= 10);
        // Channel H 5'
        auto [ph, sh] = resum(profile.at_pre_terminal_by_pfx,
                               profile.at_stop_terminal_by_pfx, skip5);
        double ni_h = profile.at_pre_interior + profile.at_stop_interior;
        StopChannelCountTable h_ct = make_stop_count_table(stop_channel_spec('H'),
            ph, sh, profile.at_pre_interior, profile.at_stop_interior, 0.0, 0.0);
        replace_ct_row('H', h_ct);
        profile.channel_h_z = static_cast<float>(h_ct.post_clamp_z);
        profile.channel_h_or = static_cast<float>(stop_channel_or_haldane(h_ct));
        auto [ph2, sh2] = resum(profile.at_pre_terminal_p2plus_by_pfx,
                                 profile.at_stop_terminal_p2plus_by_pfx, skip5);
        profile.channel_h_z_p2plus = binom_z_clamped(sh2, ph2+sh2,
                                                 profile.at_stop_interior, ni_h);
        profile.at_stop_rate_terminal = static_cast<float>(h_ct.raw_rate_term);
        profile.at_stop_rate_interior = static_cast<float>(h_ct.raw_rate_int);
        profile.channel_h_valid = (h_ct.z_den_term >= 10 && h_ct.z_den_int >= 10);
        // C5: recompute h_z sign-consistency after adapter-prefix exclusion.
        profile.channel_h_z_consistent =
            std::isfinite(profile.channel_h_z) && std::isfinite(profile.channel_h_z_p2plus) &&
            ((profile.channel_h_z >= 0.0f) == (profile.channel_h_z_p2plus >= 0.0f));
        for (auto c : excl_5p) if (c < 4096) ++n_excl;
    }

    if (!excl_3p.empty()) {
        // 3' side: channels F/G/H have no separate z-score fields; update rate
        // fields so JSON output reflects the adapter-excluded terminal window.
        auto [pf3, sf3] = resum(profile.ca_pre_terminal_3p_by_pfx,
                                 profile.ca_stop_terminal_3p_by_pfx, skip3);
        profile.ca_stop_rate_terminal_3prime = (pf3+sf3 > 0) ? sf3/(pf3+sf3) : 0.0;
        profile.channel_f3_valid = (pf3+sf3 >= 10);
        auto [pg3, sg3] = resum(profile.cg_pre_terminal_3p_by_pfx,
                                 profile.cg_stop_terminal_3p_by_pfx, skip3);
        profile.cg_stop_rate_terminal_3prime = (pg3+sg3 > 0) ? sg3/(pg3+sg3) : 0.0;
        profile.channel_g3_valid = (pg3+sg3 >= 10);
        auto [ph3, sh3] = resum(profile.at_pre_terminal_3p_by_pfx,
                                 profile.at_stop_terminal_3p_by_pfx, skip3);
        profile.at_stop_rate_terminal_3prime = (ph3+sh3 > 0) ? sh3/(ph3+sh3) : 0.0;
        profile.channel_h3_valid = (ph3+sh3 >= 10);
        for (auto c : excl_3p) if (c < 4096) ++n_excl;
    }

    profile.fgh_adapter_prefixes_excluded = n_excl;
}

} // namespace taph
