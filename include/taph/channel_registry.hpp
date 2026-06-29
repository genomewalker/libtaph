#pragma once
// Layer-1 of the stop-channel registry: a declarative ChannelSpec table that states, ONCE, each
// channel's mechanism, observable-vs-inferred semantics, denominator/shadow policy, context
// stratification, variance family, cap policy and inference class. This is what ends the "comparable
// channels" fiction: F's deamination-shadow denominator and Mantel-Haenszel stratification are TYPED
// capability fields here (has_deam_shadow / has_mh_stratification), not parallelism implied by
// hand-rolled code blocks. compute_stop_channel(spec, counts) is the single producer driven by the
// 5' stop channels; the 3' channels are declared here too with inference_class=RATE_ONLY so their
// non-membership in the stop-channel z/OR inference is EXPLICIT, not implicit.
#include <array>
#include <cstddef>

namespace taph {

enum class MechanismStatus { ESTABLISHED, EMPIRICAL };

enum class StrandFrame { TOP_5P, TOP_3P, COMPLEMENT };

// Primary effect-size / variance family for the channel's adjusted statistic.
//   MANTEL_HAENSZEL_OR  - F: stratified over context; common OR + RBG variance is interpretable.
//   PLAIN_2x2_OR        - G/H: single-stratum 2x2 odds ratio (Haldane-Anscombe +0.5) is the primary
//                         effect size; the pooled-Bernoulli z stays exploratory (correlated reads).
//   POOLED_BERNOULLI    - retained for channels whose primary statistic is still only the pooled z.
enum class VarianceFamily { POOLED_BERNOULLI, MANTEL_HAENSZEL_OR, PLAIN_2x2_OR };

enum class CapPolicy { CLAMP_ZCAP, NONE };

// How the channel participates in inference. STOP_CHANNEL: full Layer-0 count table + pooled z +
// effect-size OR (the 5' F/G/H). RATE_ONLY: terminal/baseline rates only, no z/OR/count-table
// (the 3' F3/G3/H3) — declared so "no z here" is a stated contract, not a silent omission.
enum class InferenceClass { STOP_CHANNEL, RATE_ONLY };

struct ChannelSpec {
    const char*      channel_id;            // "F","G","H","F3","G3","H3"
    char             channel_type;          // 'F' | 'G' | 'H'
    const char*      json_name;             // emitted legend "name"
    const char*      json_description;      // emitted legend "description"
    const char*      json_mechanism;        // emitted legend "mechanism" string
    MechanismStatus  mechanism_status;
    const char*      observable_name;       // measured substitution, e.g. "C_to_A_complement_asymmetry"
    const char*      inferred_lesion;       // earned lesion or nullptr (G/H assert none)
    StrandFrame      strand;
    bool             has_deam_shadow;       // F/F3: deamination shadow exists for this channel
    bool             shadow_in_z;           // F: shadow folded into the z denominator
    bool             shadow_in_rate;        // F: NOT folded into the rate denominator (rates are shadow-free)
    bool             has_mh_stratification; // F only: 3 context strata feed a Mantel-Haenszel test
    int              n_strata;              // F: 3, others: 1
    bool             applicable_in_ss;      // F/G/H: false (applicable == !is_ss)
    VarianceFamily   variance;
    CapPolicy        cap;
    InferenceClass   inference_class;       // STOP_CHANNEL (5' F/G/H) | RATE_ONLY (3' F3/G3/H3)
    // Correction 5: estimand metadata as single source of truth. The emitter READS these from the
    // registry rather than hardcoding the legend, so the registry is the one authority for what each
    // channel's z/OR actually estimates and under what assumptions it holds.
    const char*      estimand;              // what the channel's primary statistic estimates
    const char*      assumptions;           // semicolon-separated conditions under which the estimand holds
};

// The 5' stop channels as they behave TODAY (F: shadow + MH; G/H: bare single-stratum z). The
// asymmetry is declared, not implied.
inline constexpr std::array<ChannelSpec, 3> kStopChannels5p = {{
    { "F", 'F', "8_oxog_complement",
      "C to A oxidative stop codons (TCA/TCG/TAC/TGC); bottom-strand 8-oxoguanine. Scope: terminal stop-enrichment proxy for an established lesion -- a null means no TERMINAL oxidation enrichment and is insensitive to interior/non-terminal oxidation (the primary reference-free oxidation readout is oxo_two_marker)",
      "oxidative_guanine_8_oxog",
      MechanismStatus::ESTABLISHED,
      "C_to_A_complement_asymmetry", "bottom_strand_8oxoG", StrandFrame::COMPLEMENT,
      /*has_deam_shadow*/ true, /*shadow_in_z*/ true, /*shadow_in_rate*/ false,
      /*has_mh_stratification*/ true, /*n_strata*/ 3, /*applicable_in_ss*/ false,
      VarianceFamily::MANTEL_HAENSZEL_OR, CapPolicy::CLAMP_ZCAP, InferenceClass::STOP_CHANNEL,
      "terminal_vs_interior_C_to_A_stop_enrichment_common_odds_ratio",
      "stable_terminal_vs_interior_context_composition;deamination_shadow_separable;interior_reference_is_oxidation_free" },

    { "G", 'G', "cg_stop_enrichment",
      "C to G stop codons (TCA/TAC to TGA/TAG); empirical terminal stop-enrichment. NOT an earned hydantoin assignment (hydantoins are guanine over-oxidation products; this is a complement-strand C to G observation)",
      "empirical_stop_enrichment",
      MechanismStatus::EMPIRICAL,
      "C_to_G_stop_enrichment", nullptr, StrandFrame::COMPLEMENT,
      /*has_deam_shadow*/ false, /*shadow_in_z*/ false, /*shadow_in_rate*/ false,
      /*has_mh_stratification*/ false, /*n_strata*/ 1, /*applicable_in_ss*/ false,
      VarianceFamily::PLAIN_2x2_OR, CapPolicy::CLAMP_ZCAP, InferenceClass::STOP_CHANNEL,
      "terminal_vs_interior_C_to_G_stop_enrichment_odds_ratio",
      "stable_terminal_vs_interior_context_composition;empirical_not_an_earned_hydantoin_lesion" },

    { "H", 'H', "at_stop_enrichment",
      "A to T stop codons (AAA/AAG/AGA to TAA/TAG/TGA); empirical terminal stop-enrichment. No established direct adenine-oxidation to A to T pathway",
      "empirical_stop_enrichment",
      MechanismStatus::EMPIRICAL,
      "A_to_T_stop_enrichment", nullptr, StrandFrame::TOP_5P,
      /*has_deam_shadow*/ false, /*shadow_in_z*/ false, /*shadow_in_rate*/ false,
      /*has_mh_stratification*/ false, /*n_strata*/ 1, /*applicable_in_ss*/ false,
      VarianceFamily::PLAIN_2x2_OR, CapPolicy::CLAMP_ZCAP, InferenceClass::STOP_CHANNEL,
      "terminal_vs_interior_A_to_T_stop_enrichment_odds_ratio",
      "stable_terminal_vs_interior_context_composition;empirical_no_established_adenine_oxidation_pathway" },
}};

// The 3' channels: terminal/baseline RATES only. They have no terminal-vs-interior stop-enrichment
// z/OR (the interior reference at the 3' end is not separated the way it is at 5'), so they are
// declared RATE_ONLY rather than left as implicit non-members of the stop-channel inference.
inline constexpr std::array<ChannelSpec, 3> kStopChannels3p = {{
    { "F3", 'F', "8_oxog_complement_3prime",
      "C to A oxidative stop rate at the 3' end (rate-only; no terminal-enrichment z/OR inference)",
      "oxidative_guanine_8_oxog",
      MechanismStatus::ESTABLISHED,
      "C_to_A_complement_asymmetry", "bottom_strand_8oxoG", StrandFrame::COMPLEMENT,
      /*has_deam_shadow*/ true, /*shadow_in_z*/ false, /*shadow_in_rate*/ false,
      /*has_mh_stratification*/ false, /*n_strata*/ 1, /*applicable_in_ss*/ false,
      VarianceFamily::POOLED_BERNOULLI, CapPolicy::NONE, InferenceClass::RATE_ONLY,
      "terminal_and_interior_C_to_A_stop_rate_3prime",
      "rate_only;no_terminal_enrichment_inference_at_3prime" },

    { "G3", 'G', "cg_stop_enrichment_3prime",
      "C to G stop rate at the 3' end (rate-only; empirical, no z/OR inference)",
      "empirical_stop_enrichment",
      MechanismStatus::EMPIRICAL,
      "C_to_G_stop_enrichment", nullptr, StrandFrame::COMPLEMENT,
      /*has_deam_shadow*/ false, /*shadow_in_z*/ false, /*shadow_in_rate*/ false,
      /*has_mh_stratification*/ false, /*n_strata*/ 1, /*applicable_in_ss*/ false,
      VarianceFamily::POOLED_BERNOULLI, CapPolicy::NONE, InferenceClass::RATE_ONLY,
      "terminal_and_interior_C_to_G_stop_rate_3prime",
      "rate_only;empirical;no_terminal_enrichment_inference_at_3prime" },

    { "H3", 'H', "at_stop_enrichment_3prime",
      "A to T stop rate at the 3' end (rate-only; empirical, no z/OR inference)",
      "empirical_stop_enrichment",
      MechanismStatus::EMPIRICAL,
      "A_to_T_stop_enrichment", nullptr, StrandFrame::TOP_3P,
      /*has_deam_shadow*/ false, /*shadow_in_z*/ false, /*shadow_in_rate*/ false,
      /*has_mh_stratification*/ false, /*n_strata*/ 1, /*applicable_in_ss*/ false,
      VarianceFamily::POOLED_BERNOULLI, CapPolicy::NONE, InferenceClass::RATE_ONLY,
      "terminal_and_interior_A_to_T_stop_rate_3prime",
      "rate_only;empirical;no_terminal_enrichment_inference_at_3prime" },
}};

// F4: the 5-hydroxycytosine oxidative-deamination channel. Registered declaratively here as the
// single authority for its semantics, but kept OUT of kStopChannels5p/3p: it is NOT a terminal
// stop channel and does NOT participate in the F/G/H pooled-z / OR stop-channel inference
// (inference_class = RATE_ONLY). Its readout is the interior C->T pathway split
// (context_primitives.deamination.ct_pathway_split), which is PROVISIONAL: hydrolytic and
// oxidative deamination give the identical C->T substitution, so the route is not separable
// reference-free (mechanism_status = EMPIRICAL, inferred_lesion asserts the candidate lesion is
// not earned from the substitution alone).
inline constexpr ChannelSpec kOxidativeDeaminationChannel = {
    "O", 'O', "oxidative_deamination_5ohc",
    "Interior C to T via 5-hydroxycytosine (oxidative deamination), as opposed to hydrolytic C to T. "
    "PROVISIONAL: both routes yield the same C to T substitution, so they are not separable "
    "reference-free; the split is an oxidation-coupling attribution, not an identified decomposition.",
    "oxidative_deamination_5hydroxycytosine",
    MechanismStatus::EMPIRICAL,
    "interior_C_to_T_oxidation_coupling", nullptr, StrandFrame::TOP_5P,
    /*has_deam_shadow*/ false, /*shadow_in_z*/ false, /*shadow_in_rate*/ false,
    /*has_mh_stratification*/ false, /*n_strata*/ 1, /*applicable_in_ss*/ true,
    VarianceFamily::POOLED_BERNOULLI, CapPolicy::NONE, InferenceClass::RATE_ONLY,
    "oxidation_coupled_fraction_of_interior_C_to_T_excess",
    "provisional;routes_not_separable_reference_free;requires_external_validation"
};

inline constexpr const ChannelSpec& stop_channel_spec(char type) {
    for (const auto& s : kStopChannels5p)
        if (s.channel_type == type) return s;
    return kStopChannels5p[0];
}

}  // namespace taph
