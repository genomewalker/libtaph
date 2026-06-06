// Unit test for the Wave-3 single-strand overhang kernel: the bulk damage decay kernel must
// behave DIFFERENTLY for ds vs ss at the 5' terminus (p=0).
//
//   ds: terminus excluded from the decay (0) -- the p0 spike absorbs terminal artifacts.
//   ss: terminus is genuine single-strand overhang, included with full weight r(0)=1, but ONLY
//       when the overhang is identifiable; otherwise it falls back to the ds form (0).
//
// Self-contained: no gtest / catch, pure <cassert>. Invoked via ctest after -DTAPH_BUILD_TESTS=ON.

#include "taph/bulk_damage_model.hpp"

#include <cassert>
#include <cstdio>

int main() {
    using taph::BulkDamageModel;
    const double W = BulkDamageModel::terminal_kernel_weight(true,  true);   // ss, identifiable
    const double Wss_degen = BulkDamageModel::terminal_kernel_weight(true,  false);  // ss, not identifiable
    const double Wds = BulkDamageModel::terminal_kernel_weight(false, true);   // ds (gate irrelevant)
    const double Wds0 = BulkDamageModel::terminal_kernel_weight(false, false);  // ds

    // The headline taphonomy A-bar: ds != ss at the terminus.
    assert(Wds != W && "ds and ss terminal kernels must differ");
    assert(W   == 1.0 && "identifiable ss overhang carries full terminal weight r(0)=1");
    assert(Wds == 0.0 && "ds terminus is excluded from the decay");
    assert(Wds0 == 0.0 && "ds terminus is excluded regardless of the gate flag");

    // The degenerate fallback: a non-identifiable ss overhang behaves like ds (no terminal weight).
    assert(Wss_degen == 0.0 && "non-identifiable ss overhang falls back to the ds form");
    assert(Wss_degen == Wds0 && "ss fallback is byte-identical to the ds terminal kernel");

    std::printf("ss_overhang_kernel: PASS (ds=%.1f != ss=%.1f; degenerate ss=%.1f == ds)\n",
                Wds, W, Wss_degen);
    return 0;
}
