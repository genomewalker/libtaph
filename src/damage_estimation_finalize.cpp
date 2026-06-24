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
    // finalize_oxidation must run AFTER finalize_libtype: its s_gt threshold
    // (is_ss_lib) and oxidation_like 3'-deam channel (use_3prime_deam) read
    // profile.library_type, which defaults to DOUBLE_STRANDED until finalize_libtype
    // resolves it. Running oxidation first gave auto-detected SS samples DS logic.
    // finalize_libtype itself depends only on finalize_init + finalize_decay +
    // finalize_context (damage_rate_5/3prime), none of which oxidation produces.
    finalize_context(profile, ctx);
    finalize_libtype(profile, ctx);
    finalize_oxidation(profile, ctx);
    finalize_bulk(profile);
    profile.tau                   = finalize_tau(profile);
    profile.pi_shape              = fit_pi_shape(profile);
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

    // Bilateral CpG Δ: same logic as profile_json.cpp deam_context_spectrum.
    // Uses tri arrays (not tetra) to include the d=1 position's CpG signal.
    {
        static constexpr int RC4[4] = {3,2,1,0};
        auto ctx_ex = [&](const std::array<uint64_t,64>& term,
                          const std::array<uint64_t,64>& intr,
                          int mid_ref, int mid_dam, int p, int n) -> double {
            int ir = p*16 + mid_ref*4 + n;
            int id = p*16 + mid_dam*4 + n;
            uint64_t tn = term[ir] + term[id];
            uint64_t xn = intr[ir] + intr[id];
            if (!tn || !xn) return 0.0;
            return (double)term[id]/tn - (double)intr[id]/xn;
        };
        double ex5[16], ex3g_rc[16];
        for (int i = 0; i < 16; ++i) {
            int p = i/4, n = i%4;
            ex5[i]     = ctx_ex(profile.tri_5prime_terminal, profile.tri_5prime_interior, 1,3, p,n);
            ex3g_rc[i] = ctx_ex(profile.tri_3prime_terminal, profile.tri_3prime_interior, 2,0, RC4[n],RC4[p]);
        }
        auto cpgd_of = [](const double* ex) {
            double cpg = 0.0, non = 0.0;
            for (int i = 0; i < 16; ++i) {
                if (i % 4 == 2) cpg += ex[i]; else non += ex[i];
            }
            return cpg / 4.0 - non / 12.0;
        };
        profile.cpg_delta_bilateral = (float)std::min(cpgd_of(ex5), cpgd_of(ex3g_rc));
    }
    finalize_dmax(profile, ctx);   // must precede finalize_pi: pi reads the calibrated per-end d_max
    finalize_pi(profile);
    finalize_preservation(profile);
}

} // namespace taph
