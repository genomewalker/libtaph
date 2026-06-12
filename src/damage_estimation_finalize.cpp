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
#include "taph/read_ancient_llr.hpp"  // finalize_tau

namespace taph {

void FrameSelector::finalize_sample_profile(SampleDamageProfile& profile) {
    if (profile.n_reads == 0) return;
    FinalCtx ctx;
    finalize_init(profile, ctx);
    finalize_decay(profile, ctx);
    finalize_oxidation(profile, ctx);
    finalize_context(profile, ctx);
    finalize_libtype(profile, ctx);
    finalize_bulk(profile);
    profile.tau                   = finalize_tau(profile);
    profile.scission              = finalize_scission(profile);
    profile.epsilon               = finalize_epsilon(profile);
    profile.gc_depletion          = finalize_gc_depletion(profile);
    profile.oxidation_comovement  = finalize_oxidation_comovement(profile);

    // Layer-1 burial fingerprint: Θ = ln(γ/f0), φ_share = σ₀/(σ₀+f0)
    {
        const auto& sc  = profile.scission;
        const auto& tau = profile.tau;
        const auto& ox  = profile.gc_depletion;
        BurialFingerprint bf;
        if (sc.fitted && sc.gamma > 0.0 && tau.f0 > 0.0) {
            bf.theta = std::log(sc.gamma / tau.f0);
            if (sc.hi > sc.lo) {
                const double se_gamma = 0.5 * (sc.hi - sc.lo) / 1.96;
                const double se_f0    = 0.05 * tau.f0;
                const double se_theta = std::sqrt((se_gamma / sc.gamma) * (se_gamma / sc.gamma)
                                                + (se_f0 / tau.f0)     * (se_f0 / tau.f0));
                bf.theta_lo = bf.theta - 1.96 * se_theta;
                bf.theta_hi = bf.theta + 1.96 * se_theta;
            }
            bf.fitted = true;
        }
        bf.overhang_fraction = tau.overhang_fraction;
        bf.tau_hat           = tau.point;
        // phi_share = σ₀/(σ₀+f0): relative GC-depletion proxy (composition-confounded upper bound).
        // Requires f0 > 0.005: when f0≈0 the ratio is uninformative. NOT an absolute oxidation fraction.
        if (ox.fitted && ox.sigma0 > 0.0 && tau.f0 > 0.005) {
            const double denom = ox.sigma0 + tau.f0;
            bf.phi_share = denom > 0.0 ? ox.sigma0 / denom : -1.0;
        }
        profile.burial = bf;
    }

    finalize_pi(profile);
    finalize_dmax(profile, ctx);
    finalize_preservation(profile);
}

} // namespace taph
