// Primary reference-free oxidation estimator: deamination-coupled G→T regression.
//
// Wraps compute_oxo_two_marker from library_interpretation.hpp so that the result is stored
// in SampleDamageProfile::oxidation_comovement during the standard finalize pass, rather than
// being computed inline at JSON-serialization time (profile_json.cpp) and discarded.
//
// Why this is the correct level: the OxoTwoMarkerBins accumulator (s1×s2×GC×length, 256 cells)
// is populated by the fqdup damage pass via damage_estimation_update. compute_oxo_two_marker
// regresses the interior Chargaff D over these cells: beta1 = deam-marker-coupled G→T slope.
// beta1_z≈12 across FLB samples confirms that the signal is strong and already in the data.
// Neither finalize_sigma nor finalize_epsilon can recover this: marginal base counts cannot
// separate oxidation from composition (algebraically identical in a pooled DS library).

#include "taph/library_interpretation.hpp"
#include "taph/read_ancient_llr.hpp"

namespace taph {

OxoTwoMarkerResult finalize_oxidation_comovement(const SampleDamageProfile& profile) {
    const bool is_ss = (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
    return compute_oxo_two_marker(profile, is_ss);
}

} // namespace taph
