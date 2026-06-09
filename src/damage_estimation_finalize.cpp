// FrameSelector::finalize_sample_profile and supporting decay/ctx helpers.
#include "taph/frame_selector_decl.hpp"
#include "taph/codon_tables.hpp"
#include "taph/hexamer_tables.hpp"
#include "taph/library_interpretation.hpp"
#include "taph/channel_count_table.hpp"
#include "taph/channel_registry.hpp"
#include "damage_estimation_detail.hpp"
#include <algorithm>
#include <cmath>
#include <array>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>
#include "damage_estimation_finalize_ctx.hpp"
#include "taph/frame_selector_decl.hpp"

namespace taph {

void FrameSelector::finalize_sample_profile(SampleDamageProfile& profile) {
    if (profile.n_reads == 0) return;
    FinalCtx ctx;
    finalize_init(profile, ctx);
    finalize_decay(profile, ctx);
    finalize_oxidation(profile, ctx);
    finalize_context(profile, ctx);
    finalize_libtype(profile, ctx);
    finalize_dmax(profile, ctx);
    finalize_preservation(profile);
}

} // namespace taph
