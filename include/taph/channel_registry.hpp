#pragma once
// Layer-1 of the stop-channel registry: a declarative ChannelSpec table that states, ONCE, each
// stop channel's mechanism, observable-vs-inferred semantics, denominator/shadow policy, context
// stratification, variance family and cap policy. This is what ends the "four comparable stop
// channels" fiction: F's deamination-shadow denominator and Mantel-Haenszel stratification are
// TYPED capability fields here (has_deam_shadow / has_mh_stratification), not parallelism implied by
// six hand-rolled code blocks. compute_stop_channel(spec, counts) (P3) becomes the single producer
// driven by this table; until then the live blocks and this table are kept in lock-step by the
// Layer-0 count-table golden gate.
#include <array>
#include <cstddef>

namespace taph {

enum class ChannelMechanism {
    DEAMINATION_OX_SHADOW,     // F: C->A complement-strand 8-oxoG, deamination-confounded
    STOP_ENRICHMENT_EMPIRICAL, // G/H: terminal stop-enrichment with no established named lesion
};

enum class MechanismStatus { ESTABLISHED, EMPIRICAL };

enum class StrandFrame { TOP_5P, TOP_3P, COMPLEMENT };

// Primary effect-size / variance family for the channel's adjusted statistic.
//   MANTEL_HAENSZEL_OR  - F: stratified over context; common OR + RBG variance is interpretable.
//   PLAIN_2x2_OR        - G/H: single-stratum 2x2 odds ratio (Haldane-Anscombe +0.5) is the primary
//                         effect size; the pooled-Bernoulli z stays exploratory (correlated reads).
//   POOLED_BERNOULLI    - retained for channels whose primary statistic is still only the pooled z.
enum class VarianceFamily { POOLED_BERNOULLI, MANTEL_HAENSZEL_OR, PLAIN_2x2_OR };

enum class CapPolicy { CLAMP_ZCAP, NONE };

struct ChannelSpec {
    const char*      channel_id;            // "F","G","H" (+ "F3"/"G3"/"H3" when 3' z is migrated)
    char             channel_type;          // 'F' | 'G' | 'H'
    const char*      json_name;             // emitted legend "name"
    const char*      json_description;      // emitted legend "description"
    const char*      json_mechanism;        // emitted legend "mechanism" string
    ChannelMechanism mechanism;
    MechanismStatus  mechanism_status;
    const char*      observable_name;       // measured substitution, e.g. "C_to_A_complement_asymmetry"
    const char*      inferred_lesion;       // earned lesion or nullptr (G/H assert none)
    StrandFrame      strand;
    bool             has_deam_shadow;       // F only: deamination shadow exists for this channel
    bool             shadow_in_z;           // F: shadow folded into the z denominator
    bool             shadow_in_rate;        // F: NOT folded into the rate denominator (rates are shadow-free)
    bool             has_mh_stratification; // F only: 3 context strata feed a Mantel-Haenszel test
    int              n_strata;              // F: 3, G/H: 1
    bool             applicable_in_ss;      // F/G/H: false (applicable == !is_ss)
    VarianceFamily   variance;
    CapPolicy        cap;
};

// The 5' stop channels as they behave TODAY (F: shadow + MH; G/H: bare single-stratum z). The
// asymmetry is declared, not implied. 3' (F3/G3/H3) join when their z is migrated in P3.
inline constexpr std::array<ChannelSpec, 3> kStopChannels5p = {{
    { "F", 'F', "8_oxog_complement",
      "C to A oxidative stop codons (TCA/TCG/TAC/TGC); bottom-strand 8-oxoguanine",
      "oxidative_guanine_8_oxog",
      ChannelMechanism::DEAMINATION_OX_SHADOW, MechanismStatus::ESTABLISHED,
      "C_to_A_complement_asymmetry", "bottom_strand_8oxoG", StrandFrame::COMPLEMENT,
      /*has_deam_shadow*/ true, /*shadow_in_z*/ true, /*shadow_in_rate*/ false,
      /*has_mh_stratification*/ true, /*n_strata*/ 3, /*applicable_in_ss*/ false,
      VarianceFamily::MANTEL_HAENSZEL_OR, CapPolicy::CLAMP_ZCAP },

    { "G", 'G', "cg_stop_enrichment",
      "C to G stop codons (TCA/TAC to TGA/TAG); empirical terminal stop-enrichment. NOT an earned hydantoin assignment (hydantoins are guanine over-oxidation products; this is a complement-strand C to G observation)",
      "empirical_stop_enrichment",
      ChannelMechanism::STOP_ENRICHMENT_EMPIRICAL, MechanismStatus::EMPIRICAL,
      "C_to_G_stop_enrichment", nullptr, StrandFrame::COMPLEMENT,
      /*has_deam_shadow*/ false, /*shadow_in_z*/ false, /*shadow_in_rate*/ false,
      /*has_mh_stratification*/ false, /*n_strata*/ 1, /*applicable_in_ss*/ false,
      VarianceFamily::PLAIN_2x2_OR, CapPolicy::CLAMP_ZCAP },

    { "H", 'H', "at_stop_enrichment",
      "A to T stop codons (AAA/AAG/AGA to TAA/TAG/TGA); empirical terminal stop-enrichment. No established direct adenine-oxidation to A to T pathway",
      "empirical_stop_enrichment",
      ChannelMechanism::STOP_ENRICHMENT_EMPIRICAL, MechanismStatus::EMPIRICAL,
      "A_to_T_stop_enrichment", nullptr, StrandFrame::TOP_5P,
      /*has_deam_shadow*/ false, /*shadow_in_z*/ false, /*shadow_in_rate*/ false,
      /*has_mh_stratification*/ false, /*n_strata*/ 1, /*applicable_in_ss*/ false,
      VarianceFamily::PLAIN_2x2_OR, CapPolicy::CLAMP_ZCAP },
}};

inline constexpr const ChannelSpec& stop_channel_spec(char type) {
    for (const auto& s : kStopChannels5p)
        if (s.channel_type == type) return s;
    return kStopChannels5p[0];
}

}  // namespace taph
