#pragma once
#include "taph/sample_damage_profile.hpp"
#include <array>
namespace taph {
struct FinalCtx {
    double baseline_tc   = 0.0;
    double baseline_ag   = 0.0;
    float  fit_lambda_5p = 0.0f;
    float  fit_lambda_3p = 0.0f;
    std::array<double,15> ctrl_freq_3p  = {};
    std::array<double,15> ctrl_total_3p = {};
};
void finalize_init        (SampleDamageProfile&, FinalCtx&);
void finalize_decay       (SampleDamageProfile&, FinalCtx&);
void finalize_oxidation   (SampleDamageProfile&, const FinalCtx&);
void finalize_context     (SampleDamageProfile&, FinalCtx&);
void finalize_libtype     (SampleDamageProfile&, const FinalCtx&);
void finalize_bulk        (SampleDamageProfile&);
void compute_codon_did    (SampleDamageProfile&);  // composition-immune codon DiD validator (finalize_decay.cpp)
void finalize_pi          (SampleDamageProfile&);  // validated reference-free pi (read_ancient_llr.cpp)
void finalize_dmax        (SampleDamageProfile&, const FinalCtx&);
void finalize_preservation(SampleDamageProfile&);
} // namespace taph
