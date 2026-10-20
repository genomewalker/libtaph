// The d_max split needs the lower 95% edge of d_max above max δ, not just the point.
#include "taph/bulk_damage_model.hpp"
#include <cstdio>

int main() {
    const double sum_w = 1e4;   // SE 0.01
    const bool point_only = taph::dmax_split_identified(sum_w, 0.30, 0.29);   // 0.30 > 0.29, 0.2804 < 0.29
    const bool separated  = taph::dmax_split_identified(sum_w, 0.30, 0.25);
    const bool no_weight  = taph::dmax_split_identified(0.0, 0.30, 0.25);
    std::printf("point_only=%d separated=%d no_weight=%d\n", point_only, separated, no_weight);
    return (!point_only && separated && !no_weight) ? 0 : 1;
}
