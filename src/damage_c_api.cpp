// damage_c_api.cpp — extern "C" wrapper for libdart-damage (C++17)

#include "dart/damage_c_api.h"
#include "dart/frame_selector_decl.hpp"
#include "dart/sample_damage_profile.hpp"

#include <cassert>
#include <cmath>
#include <cstring>
#include <new>

// ---------------------------------------------------------------------------
// Internal state bundled behind the opaque handle

struct dart_profile_t {
    dart::SampleDamageProfile profile;
    bool finalized = false;
};

// ---------------------------------------------------------------------------
// Pass 1 – estimation

dart_profile_t *dart_profile_create(void) {
    return new (std::nothrow) dart_profile_t{};
}

void dart_profile_destroy(dart_profile_t *p) {
    delete p;
}

void dart_profile_add_read(dart_profile_t *p, const char *seq, size_t len) {
    if (!p || !seq || len == 0) { return; }
    dart::FrameSelector::update_sample_profile(p->profile,
                                               std::string(seq, len));
}

void dart_profile_finalize(dart_profile_t *p) {
    if (!p || p->finalized) { return; }
    dart::FrameSelector::finalize_sample_profile(p->profile);
    p->finalized = true;
}

// ---------------------------------------------------------------------------
// Accessors

float dart_profile_dmax(const dart_profile_t *p) {
    return p ? p->profile.d_max_combined : 0.0f;
}

float dart_profile_lambda5(const dart_profile_t *p) {
    return p ? p->profile.lambda_5prime : 0.0f;
}

float dart_profile_lambda3(const dart_profile_t *p) {
    return p ? p->profile.lambda_3prime : 0.0f;
}

int dart_profile_library_type(const dart_profile_t *p) {
    if (!p) { return 0; }
    using LT = dart::SampleDamageProfile::LibraryType;
    switch (p->profile.library_type) {
        case LT::DOUBLE_STRANDED: return 1;
        case LT::SINGLE_STRANDED: return 2;
        default:                  return 0;
    }
}

int dart_profile_damage_validated(const dart_profile_t *p) {
    return (p && p->profile.damage_validated) ? 1 : 0;
}

int dart_profile_damage_artifact(const dart_profile_t *p) {
    return (p && p->profile.damage_artifact) ? 1 : 0;
}

int dart_profile_is_reliable(const dart_profile_t *p) {
    return (p && !p->profile.is_detection_unreliable()) ? 1 : 0;
}

// ---------------------------------------------------------------------------
// Pass 2 – per-read correction

// Return position-dependent C→T damage probability at dist bases from 5' end.
static inline float ct_prob(const dart::SampleDamageProfile &sp,
                             size_t dist_from_5prime) {
    if (dist_from_5prime < dart::SampleDamageProfile::N_POS) {
        return sp.damage_rate_5prime[dist_from_5prime];
    }
    // Exponential extrapolation beyond the measured window
    float d = static_cast<float>(dist_from_5prime);
    return sp.d_max_5prime * std::exp(-sp.lambda_5prime * d);
}

// Return position-dependent G→A damage probability at dist bases from 3' end.
static inline float ga_prob(const dart::SampleDamageProfile &sp,
                             size_t dist_from_3prime) {
    if (dist_from_3prime < dart::SampleDamageProfile::N_POS) {
        return sp.damage_rate_3prime[dist_from_3prime];
    }
    float d = static_cast<float>(dist_from_3prime);
    return sp.d_max_3prime * std::exp(-sp.lambda_3prime * d);
}

size_t dart_correct_read(const dart_profile_t *p,
                         const char           *seq,
                         size_t                len,
                         char                 *out_buf,
                         float                 confidence_threshold) {
    if (!p || !seq || !out_buf || len == 0) { return 0; }

    std::memcpy(out_buf, seq, len);
    out_buf[len] = '\0';

    size_t corrections = 0;
    const dart::SampleDamageProfile &sp = p->profile;

    for (size_t i = 0; i < len; ++i) {
        char c = seq[i];
        size_t dist3 = len - 1 - i;

        if ((c == 'T' || c == 't') && ct_prob(sp, i) >= confidence_threshold) {
            out_buf[i] = (c == 'T') ? 'C' : 'c';
            ++corrections;
        } else if ((c == 'A' || c == 'a') && ga_prob(sp, dist3) >= confidence_threshold) {
            out_buf[i] = (c == 'A') ? 'G' : 'g';
            ++corrections;
        }
    }
    return corrections;
}
