// damage_c_api.cpp — extern "C" wrapper for libtaph (C++17)

#include "taph/damage_c_api.h"
#include "taph/frame_selector_decl.hpp"
#include "taph/sample_damage_profile.hpp"

#include <cmath>
#include <cstring>
#include <new>
#include <string_view>

// ---------------------------------------------------------------------------
// Internal state bundled behind the opaque handle

struct taph_profile_t {
    taph::SampleDamageProfile profile;
    bool finalized = false;
};

// ---------------------------------------------------------------------------
// Pass 1 – estimation

taph_profile_t *taph_profile_create(void) {
    return new (std::nothrow) taph_profile_t{};
    // Caller must check for NULL (OOM).
}

void taph_profile_destroy(taph_profile_t *p) {
    delete p;
}

void taph_profile_add_read(taph_profile_t *p, const char *seq, size_t len) {
    if (!p || p->finalized || !seq || len == 0) { return; }
    taph::FrameSelector::update_sample_profile(p->profile,
                                               std::string_view(seq, len));
}

void taph_profile_finalize(taph_profile_t *p) {
    if (!p || p->finalized) { return; }
    taph::FrameSelector::finalize_sample_profile(p->profile);
    p->finalized = true;
}

// ---------------------------------------------------------------------------
// Accessors — all require a finalized profile; return safe defaults otherwise.

float taph_profile_dmax(const taph_profile_t *p) {
    return (p && p->finalized) ? p->profile.d_max_combined : 0.0f;
}

float taph_profile_lambda5(const taph_profile_t *p) {
    return (p && p->finalized) ? p->profile.lambda_5prime : 0.0f;
}

float taph_profile_lambda3(const taph_profile_t *p) {
    return (p && p->finalized) ? p->profile.lambda_3prime : 0.0f;
}

int taph_profile_library_type(const taph_profile_t *p) {
    if (!p || !p->finalized) { return 0; }
    using LT = taph::SampleDamageProfile::LibraryType;
    switch (p->profile.library_type) {
        case LT::DOUBLE_STRANDED: return 1;
        case LT::SINGLE_STRANDED: return 2;
        default:                  return 0;
    }
}

int taph_profile_damage_validated(const taph_profile_t *p) {
    return (p && p->finalized && p->profile.damage_validated) ? 1 : 0;
}

int taph_profile_damage_artifact(const taph_profile_t *p) {
    return (p && p->finalized && p->profile.damage_artifact) ? 1 : 0;
}

int taph_profile_is_reliable(const taph_profile_t *p) {
    return (p && p->finalized && !p->profile.is_detection_unreliable()) ? 1 : 0;
}

// ---------------------------------------------------------------------------
// Pass 2 – per-read correction

// Position-dependent C→T probability at dist bases from 5' end.
static inline float ct_prob(const taph::SampleDamageProfile &sp,
                             size_t dist_from_5prime) {
    if (dist_from_5prime < static_cast<size_t>(taph::SampleDamageProfile::N_POS)) {
        return sp.damage_rate_5prime[dist_from_5prime];
    }
    float d = static_cast<float>(dist_from_5prime);
    return sp.d_max_5prime * std::exp(-sp.lambda_5prime * d);
}

// Position-dependent G→A probability at dist bases from 3' end.
static inline float ga_prob(const taph::SampleDamageProfile &sp,
                             size_t dist_from_3prime) {
    if (dist_from_3prime < static_cast<size_t>(taph::SampleDamageProfile::N_POS)) {
        return sp.damage_rate_3prime[dist_from_3prime];
    }
    float d = static_cast<float>(dist_from_3prime);
    return sp.d_max_3prime * std::exp(-sp.lambda_3prime * d);
}

size_t taph_mask_read(const taph_profile_t *p,
                      const char           *seq,
                      size_t                len,
                      char                 *out_buf,
                      float                 confidence_threshold,
                      char                  mask_char) {
    if (!p || !p->finalized || !seq || !out_buf) { return 0; }
    if (len == 0) { out_buf[0] = '\0'; return 0; }

    std::memcpy(out_buf, seq, len);
    out_buf[len] = '\0';

    size_t masked = 0;
    const taph::SampleDamageProfile &sp = p->profile;

    for (size_t i = 0; i < len; ++i) {
        char c = seq[i];
        size_t dist3 = len - 1 - i;

        bool is_ct = (c == 'T' || c == 't') && ct_prob(sp, i)    >= confidence_threshold;
        bool is_ga = (c == 'A' || c == 'a') && ga_prob(sp, dist3) >= confidence_threshold;

        if (is_ct || is_ga) {
            out_buf[i] = mask_char;
            ++masked;
        }
    }
    return masked;
}

size_t taph_correct_read(const taph_profile_t *p,
                         const char           *seq,
                         size_t                len,
                         char                 *out_buf,
                         float                 confidence_threshold) {
    if (!p || !p->finalized || !seq || !out_buf) { return 0; }
    if (len == 0) { out_buf[0] = '\0'; return 0; }

    std::memcpy(out_buf, seq, len);
    out_buf[len] = '\0';

    size_t corrections = 0;
    const taph::SampleDamageProfile &sp = p->profile;

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
