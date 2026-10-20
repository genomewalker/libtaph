// SQUAREM on the bulk damage fit: alpha = 1 is the plain two-sweep point, and the accelerated fit
// ends at least as high as plain coordinate ascent on a synthetic ds library.
#include "taph/bulk_damage_model.hpp"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

int main() {
    using namespace taph;
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
    std::printf("plain_err=%.3g plain: sweeps=%d ll=%.6f lambda=%.6f | squarem: sweeps=%d ll=%.6f lambda=%.6f\n",
                plain_err, plain.n_sweeps, plain.log_lik, plain.lambda, accel.n_sweeps, accel.log_lik, accel.lambda);
    const bool ll_ok = accel.log_lik >= plain.log_lik - 1e-6 * std::fabs(plain.log_lik);
    return (plain_ok && ll_ok) ? 0 : 1;
}
