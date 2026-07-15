// Closed-form proof, as a test: k× PCR duplication scales every cube count by k, so the fitted
// terminal-decay shape {A, λ, baseline} is EXACTLY invariant, the amplitude SE shrinks as A_se/√k,
// and the shape LRT (and n_elig) grow ~linearly in k. This is the concrete statement of "duplicates
// buy CI width and LRT significance, not rate" — the premise behind profiling dedup vs merged reads,
// and the reason pi_shape.detected (an LRT threshold) is duplicate-purchasable while {A,λ,b} are not.
#include "taph/read_ancient_llr.hpp"        // fit_pi_shape_cube, TauConfig
#include "taph/sample_damage_profile.hpp"
#include <cmath>
#include <cstdint>
#include <cstdio>

using taph::SampleDamageProfile;

// Deterministic decay cube scaled by an integer factor. Base counts are computed once at scale=1 and
// multiplied by `scale`, so the k× cube is EXACTLY k× the base (no per-scale rounding drift).
static SampleDamageProfile::PiPosCube make_cube(std::int64_t scale) {
    SampleDamageProfile::PiPosCube cube{};
    constexpr int P = SampleDamageProfile::P_PI;
    const double b = 0.05, A = 0.30, lam = 0.25;
    for (int p = 0; p < P; ++p) {
        const double r = b + (1.0 - b) * A * std::exp(-lam * p);
        const std::int64_t n0 = 1000;
        const std::int64_t k0 = std::llround(r * n0);
        cube[0][1][p].n_elig = static_cast<std::uint64_t>(n0 * scale);
        cube[0][1][p].n_deam = static_cast<std::uint64_t>(k0 * scale);
    }
    return cube;
}

static int fails = 0;
static void check(bool ok, const char* what) {
    if (!ok) { std::printf("FAIL: %s\n", what); ++fails; }
}

int main() {
    const auto f1 = taph::fit_pi_shape_cube(make_cube(1), /*exclude_c0=*/true);
    check(f1.fitted, "base fit produced a finite shape");

    for (std::int64_t k : {std::int64_t(2), std::int64_t(4), std::int64_t(8)}) {
        const auto fk = taph::fit_pi_shape_cube(make_cube(k), /*exclude_c0=*/true);
        check(fk.fitted, "k-scaled fit produced a finite shape");
        // Rate params EXACTLY invariant (scale cancels in every WLS moment).
        check(std::fabs(fk.A - f1.A) <= 1e-9,               "A invariant under k× duplication");
        check(std::fabs(fk.lambda - f1.lambda) <= 1e-9,     "lambda invariant under k× duplication");
        check(std::fabs(fk.baseline - f1.baseline) <= 1e-12,"baseline invariant under k× duplication");
        // CI width shrinks as 1/√k; LRT and eligible-site count grow linearly.
        check(std::fabs(fk.A_se * std::sqrt(double(k)) - f1.A_se) <= 1e-6 * f1.A_se + 1e-9,
              "A_se scales as A_se/sqrt(k)");
        check(std::fabs(fk.lrt / double(k) - f1.lrt) <= 1e-6 * std::fabs(f1.lrt) + 1e-9,
              "lrt scales linearly in k (duplicate-purchasable)");
        check(fk.n_elig == k * f1.n_elig, "n_elig scales linearly in k");
    }

    if (fails == 0) std::printf("pi_shape k-invariance OK\n");
    return fails == 0 ? 0 : 1;
}
