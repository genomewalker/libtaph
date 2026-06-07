#pragma once
#include "taph/sample_damage_profile.hpp"

namespace taph {

// Pre-computed per-GC-bin LLR scoring table.
// Populated once after pass-1 SampleDamageProfile finalization.
// Eliminates exp/log/clamp from the per-read worker inner loop.
struct LLRTable {
    float addT[15];   // log(pd/bg)         at each position
    float addC[15];   // log((1-pd)/(1-bg)) at each position
    int   max_pos;    // positions 0..max_pos-1 have excess >= 0.01
    float logit_pi;
};

// Populate n_tables entries from dp.gc_bins[].
// n_tables must equal SampleDamageProfile::N_GC_BINS.
void llr_table_build(LLRTable* tables, int n_tables,
                     const SampleDamageProfile& dp);

// Per-read posterior P(ancient|read) — forward orientation (5'→3').
float llr_table_weight_fwd(const char* seq, int len, const LLRTable& t);

// Per-read posterior — reverse-complement (3'→5' scan, A↔T G↔C).
float llr_table_weight_rev(const char* seq, int len, const LLRTable& t);

} // namespace taph
