// SQUAREM on the bulk damage fit:
//  - alpha = 1 is the plain two-sweep point;
//  - on a synthetic ds library the accelerated fit ends at least as high as plain;
//  - on a synthetic pervasive-artifact library (plain fit at the sweep cap) the penalised objective
//    never decreases across the accelerated loop. Convergence is printed, not asserted: over 12
//    variants of this library SQUAREM cleared the cap in only 4, so it is not a cap-clearing check;
//  - with alpha forced far outside its bounds every extrapolation falls back to the plain point.
#include "taph/bulk_damage_model.hpp"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

using namespace taph;

// A synthetic library shaped like a pervasive-artifact ds library (no deamination): both channels sit at
// an interior rate of 0.47 and share a terminal artifact that deepens with read length, so delta stays
// below 0.04 and plain coordinate ascent crawls along the a/delta ridge to the sweep cap. Counts carry
// binomial-scale noise from a fixed LCG; without it the ridge is exactly flat and SQUAREM stalls too.
static BulkDamageSuffStats artifact_like() {
    BulkDamageSuffStats s;
    uint64_t rng = 2;
    auto z = [&] {
        double u = 0.0;
        for (int i = 0; i < 4; ++i) {
            rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
            u += (rng >> 11) * 0x1.0p-53;
        }
        return (u - 2.0) * std::sqrt(3.0);
    };
    const double lens[6] = {40, 60, 85, 110, 140, 165}, scale[6] = {0.65, 0.82, 0.93, 1.0, 1.02, 1.03};
    const double n_int[6] = {6e7, 1.5e8, 2.1e8, 2.7e8, 3.5e8, 2.1e8};
    const double shape[2][15] = {
        {-0.10, 0.00, -0.08, -0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.00, 0.01, 0.01, -0.01, 0.01, 0.00},
        {-0.19, -0.14, -0.04, 0.01, 0.02, 0.03, 0.01, 0.02, 0.01, -0.01, 0.00, 0.01, -0.01, 0.00, 0.01}};
    const int L = 6;
    s.bin.resize(L); s.k_interior.resize(L); s.n_interior.resize(L); s.jstrat.resize(L);
    s.median_len.assign(lens, lens + L);
    for (int l = 0; l < L; ++l) {
        const double interior = 0.47;
        for (int c = 0; c < 2; ++c) {
            s.n_interior[l][c] = static_cast<uint64_t>(n_int[l]);
            s.k_interior[l][c] = static_cast<uint64_t>(std::llround(n_int[l] * interior));
            for (int e = 0; e < 2; ++e)
                for (int p = 0; p < BulkDamageSuffStats::N_POS; ++p) {
                    const double n = c == 0 ? 6e6 : 8.5e6;
                    s.bin[l][c][e].n[p] = static_cast<uint64_t>(n);
                    const double r = interior + scale[l] * shape[c][p];
                    s.bin[l][c][e].k[p] = static_cast<uint64_t>(std::llround(n * r + z() * std::sqrt(n * r * (1 - r))));
                }
        }
    }
    return s;
}

// Plain coordinate ascent itself dips by ~1e-14 relative (summation rounding over ~1e8 counts), so the
// check allows 1e-12 of |obj| per step.
static bool monotone(const BulkFitTrace& t, double& worst) {
    worst = 0.0;
    for (size_t i = 1; i < t.obj.size(); ++i) worst = std::min(worst, t.obj[i] - t.obj[i - 1]);
    return worst >= -1e-12 * std::fabs(t.obj.front());
}

int main() {
    const double x0[3] = {0.1, -2.5, 7.0}, x1[3] = {0.3, -2.0, 6.5}, x2[3] = {0.45, -1.7, 6.2};
    double r[3], v[3], sr2, sv2, srv, x[3];
    squarem_build_directions(x0, x1, x2, r, v, 3, sr2, sv2, srv);
    std::memcpy(x, x0, sizeof x);
    squarem_extrapolate(x, r, v, 3, 1.0);
    double plain_err = 0.0;
    for (int i = 0; i < 3; ++i) plain_err = std::max(plain_err, std::fabs(x[i] - x2[i]));
    const bool plain_ok = plain_err <= 1e-15;

    BulkDamageSuffStats s;
    const double lens[4] = {40, 60, 90, 140}, dl[4] = {0.30, 0.20, 0.10, 0.0};
    const int L = 4;
    s.bin.resize(L); s.k_interior.resize(L); s.n_interior.resize(L); s.median_len.assign(lens, lens + L);
    for (int l = 0; l < L; ++l) {
        for (int c = 0; c < 2; ++c) {
            s.n_interior[l][c] = 2000000;
            s.k_interior[l][c] = 20000;   // baseline 1%
            for (int e = 0; e < 2; ++e)
                for (int p = 0; p < BulkDamageSuffStats::N_POS; ++p) {
                    const double art = (p == 0 ? 0.02 : 0.0) + 0.004 * std::exp(-0.5 * p);
                    const double rate = 0.01 + art + (c == 0 ? dl[l] * std::exp(-0.35 * p) : 0.0);
                    s.bin[l][c][e].n[p] = 200000;
                    s.bin[l][c][e].k[p] = static_cast<uint64_t>(std::llround(200000 * rate));
                }
        }
    }
    const BulkDamageResult plain = BulkDamageModel::fit(s, 1, false);
    const BulkDamageResult accel = BulkDamageModel::fit(s, 1, true);
    std::printf("synthetic: plain_err=%.3g plain sweeps=%d ll=%.6f | squarem sweeps=%d ll=%.6f\n",
                plain_err, plain.n_sweeps, plain.log_lik, accel.n_sweeps, accel.log_lik);
    const bool ll_ok = accel.log_lik >= plain.log_lik - 1e-6 * std::fabs(plain.log_lik);

    const BulkDamageSuffStats art = artifact_like();
    const BulkDamageResult art_plain = BulkDamageModel::fit(art, 1, false);
    BulkFitTrace tr;
    const BulkDamageResult ra = BulkDamageModel::fit(art, 1, true, &tr);
    double worst;
    const bool mono_ok = monotone(tr, worst) && tr.accelerated > 0;
    std::printf("artifact-like: plain sweeps=%d ll=%.6f | obj=%.6f sweeps=%d converged=%d iterations=%zu accelerated=%d fallbacks=%d worst_step=%.3g\n",
                art_plain.n_sweeps, art_plain.log_lik, tr.obj.back(), ra.n_sweeps, ra.converged, tr.obj.size() - 1, tr.accelerated, tr.fallbacks, worst);

    BulkFitTrace forced;
    forced.force_alpha = 1e3;
    BulkDamageModel::fit(art, 1, true, &forced);
    double worst_f;
    const bool fb_ok = forced.accelerated > 0 && forced.fallbacks == forced.accelerated && monotone(forced, worst_f);
    std::printf("forced alpha=1e3: accelerated=%d fallbacks=%d worst_step=%.3g\n",
                forced.accelerated, forced.fallbacks, worst_f);

    return (plain_ok && ll_ok && mono_ok && fb_ok) ? 0 : 1;
}
