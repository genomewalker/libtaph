#pragma once

// DamageSummary — the library-level damage profile as a fixed-layout binary sidecar.
//
// This is the single-engine serialization seam: fqdup fits a SampleDamageProfile at
// merge/profile time and WRITES this summary; DART READS it via --library-profile and
// applies it instead of re-deriving the profile from the reads it annotates. The scalar
// set is exactly the AGD override contract that DART's DamageIndexWriter stamps
// (damage_index_writer.cpp lines 17-54): whatever DART reads off a SampleDamageProfile to
// build the AGD header, this summary carries. libtaph owns both the struct and its I/O so
// neither consumer re-implements the mapping.
//
// Layout is explicit little-endian, field-by-field (no struct memcpy) so padding/ABI never
// enter the on-disk format.

#include <cstdint>
#include <cstring>
#include <fstream>
#include <optional>
#include <string>
#include <utility>

#include "taph/sample_damage_profile.hpp"

namespace taph {

struct DamageSummary {
    static constexpr uint32_t MAGIC   = 0x314d4454;  // 'TDM1'
    // v2 adds the 3' terminal-decay shape. Without it a consumer restores pi_shape (5')
    // but scores the 3' end off the lambda_3prime/fit_baseline_3prime FALLBACK, while the
    // producer used the fitted 3' shape -- i.e. the sidecar silently reintroduced the very
    // per-end divergence the b3=b5 fix removed. No migration cost: nothing had ever called
    // write_damage_summary, so no v1 sidecar exists anywhere.
    // v3 adds pi_damaged_point: the mixture-based damaged-MOLECULE prevalence
    // LengthStratifiedDamageProfile::pi_joint_damaged, not a terminal-rate ratio. Unlike
    // pi_point it is NOT gated by the DamageConfidence verdict, so it survives LOW_ABUNDANCE
    // and UNDETERMINED, and it is the empirical-Bayes prior for the consumer's per-protein
    // authentication probability.
    static constexpr uint32_t VERSION = 3;

    // continuous scalars the AGD header carries
    float  d_max_combined       = 0.0f;
    float  d_max_5prime         = 0.0f;
    float  d_max_3prime         = 0.0f;
    float  lambda_5prime        = 0.3f;
    float  lambda_3prime        = 0.3f;
    float  fit_baseline_5prime  = 0.0f;
    float  fit_baseline_3prime  = 0.0f;
    float  terminal_shift_5prime = 0.0f;
    float  stop_decay_llr_5prime = 0.0f;
    double d_max_cooc_point     = -1.0;
    double pi_point             = -1.0;
    double pi_shape_lambda      = -1.0;
    double pi_shape_baseline    = -1.0;
    // v2: the 3' end's OWN fitted decay (DS G->A / SS C->T). The per-read kernel takes
    // lam3/b3 from here when fitted; carrying only the 5' shape made the consumer diverge.
    double pi_shape_3p_lambda   = -1.0;
    double pi_shape_3p_baseline = -1.0;
    // v3: prevalence of damaged MOLECULES from the length-stratified joint mixture
    // (LengthStratifiedDamageProfile::pi_joint_damaged) -- two components fitted over the
    // length x GC cells, this being the c_sites-weighted fraction of cells assigned to the
    // high-d_max component. Empirical-Bayes PRIOR for the consumer's per-protein probability.
    // Differs from SampleDamageProfile::pi_damaged (damaged_obs/total_obs) two ways that
    // matter: it is a mixture assignment over cells rather than one pooled terminal ratio, and
    // it is ungated by DamageConfidence so it is populated where pi_point is withheld.
    // NOT independent evidence. The mixture E-step likelihood reads only bins[i].d_max and
    // bins[i].c_sites (joint_damage_model.hpp) -- no length or GC term enters it -- so this is
    // a different FUNCTIONAL of the same C->T/G->A channel the per-read LLR scores, not a
    // second channel. Measured rho(pi_damaged_point, d_max) ~ 0.91 over 38 KapK libraries.
    // Clamped to [1e-3, 1-1e-3] by PI_FLOOR (joint_damage_model.hpp:494): 0.999 is a pinned
    // bound, not a fit (8/38 KapK libraries pin there).
    // Left at -1.0 (undetermined) when the LSD joint is null or <= 0 so the consumer's
    // apply_to `>= 0` guard skips it (no-authentication tier).
    float  pi_damaged_point     = -1.0f;

    // discrete/flag scalars
    uint8_t library_type      = 2;  // 0 unknown, 1 ss, 2 ds (matches AGD encoding)
    uint8_t damage_validated  = 0;
    uint8_t damage_artifact   = 0;
    uint8_t channel_b_valid   = 0;
    uint8_t pi_shape_fitted   = 0;
    uint8_t pi_shape_3p_fitted = 0;  // v2
    uint8_t pi_state          = static_cast<uint8_t>(DamageConfidence::UNDETERMINED);

    // ── build from a finalized profile (fqdup side) ─────────────────────────────
    // pi_joint_damaged is LengthStratifiedDamageProfile::pi_joint_damaged (the latent
    // damaged-molecule prevalence). Pass <= 0 (default) when the LSD is null/unavailable or
    // the joint mixture did not resolve -- pi_damaged_point then stays at its -1.0 default so
    // the consumer's apply_to guard skips the prior. Never pass p.pi_damaged here: it is the
    // pooled terminal ratio (damaged_obs/total_obs), a different and DamageConfidence-gated
    // estimand -- consumers key off this field being the mixture value.
    static DamageSummary from_profile(const SampleDamageProfile& p,
                                      double pi_joint_damaged = -1.0) {
        DamageSummary s;
        s.d_max_combined       = p.d_max_combined;
        s.d_max_5prime         = p.d_max_5prime;
        s.d_max_3prime         = p.d_max_3prime;
        s.lambda_5prime        = p.lambda_5prime;
        s.lambda_3prime        = p.lambda_3prime;
        s.fit_baseline_5prime  = p.fit_baseline_5prime;
        s.fit_baseline_3prime  = p.fit_baseline_3prime;
        s.terminal_shift_5prime = p.terminal_shift_5prime;
        s.stop_decay_llr_5prime = p.stop_decay_llr_5prime;
        s.d_max_cooc_point     = p.d_max_cooccurrence.point;
        s.pi_point             = p.pi.point;
        s.pi_shape_lambda      = p.pi_shape.lambda;
        s.pi_shape_baseline    = p.pi_shape.baseline;
        s.pi_shape_3p_lambda   = p.pi_shape_3prime.lambda;
        s.pi_shape_3p_baseline = p.pi_shape_3prime.baseline;
        s.pi_damaged_point     = (pi_joint_damaged > 0.0)   // v3: LSD joint prior; -1.0 = undetermined
                                     ? static_cast<float>(pi_joint_damaged) : -1.0f;
        switch (p.library_type) {
            case SampleDamageProfile::LibraryType::SINGLE_STRANDED: s.library_type = 1; break;
            case SampleDamageProfile::LibraryType::DOUBLE_STRANDED: s.library_type = 2; break;
            default:                                                s.library_type = 0; break;
        }
        s.damage_validated = p.damage_validated ? 1 : 0;
        s.damage_artifact  = p.damage_artifact ? 1 : 0;
        s.channel_b_valid  = p.channel_b_valid ? 1 : 0;
        s.pi_shape_fitted  = p.pi_shape.fitted ? 1 : 0;
        s.pi_shape_3p_fitted = p.pi_shape_3prime.fitted ? 1 : 0;
        s.pi_state         = static_cast<uint8_t>(p.pi.state);
        return s;
    }

    // ── apply onto a profile the consumer will score/stamp with (DART side) ─────
    // Sets exactly the fields DamageIndexWriter reads AND the fields read_ancient_posterior
    // needs (pi.point + pi.state), so the per-read prior fires before the annotate loop.
    void apply_to(SampleDamageProfile& p) const {
        p.d_max_combined       = d_max_combined;
        p.d_max_5prime         = d_max_5prime;
        p.d_max_3prime         = d_max_3prime;
        p.lambda_5prime        = lambda_5prime;
        p.lambda_3prime        = lambda_3prime;
        p.fit_baseline_5prime  = fit_baseline_5prime;
        p.fit_baseline_3prime  = fit_baseline_3prime;
        p.terminal_shift_5prime = terminal_shift_5prime;
        p.stop_decay_llr_5prime = stop_decay_llr_5prime;
        p.d_max_cooccurrence.point = d_max_cooc_point;
        p.pi.point             = pi_point;
        p.pi.state             = static_cast<DamageConfidence>(pi_state);
        p.pi_shape.lambda      = pi_shape_lambda;
        p.pi_shape.baseline    = pi_shape_baseline;
        p.pi_shape.fitted      = pi_shape_fitted != 0;
        // v2: restore the 3' end's own shape, so the consumer's per-read kernel takes the
        // SAME lam3/b3 the producer fitted instead of silently falling back to
        // lambda_3prime/fit_baseline_3prime (a 5'-derived background on a 3' channel).
        p.pi_shape_3prime.lambda   = pi_shape_3p_lambda;
        p.pi_shape_3prime.baseline = pi_shape_3p_baseline;
        p.pi_shape_3prime.fitted   = pi_shape_3p_fitted != 0;
        switch (library_type) {
            case 1:  p.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED; break;
            case 2:  p.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED; break;
            default: p.library_type = SampleDamageProfile::LibraryType::UNKNOWN; break;
        }
        p.damage_validated = damage_validated != 0;
        p.damage_artifact  = damage_artifact != 0;
        p.channel_b_valid  = channel_b_valid != 0;
        // v3: restore the ungated mixture-based damaged-molecule prevalence (LSD joint) so the
        // consumer's AGD writer stamps it as the per-protein authentication prior. It re-expresses
        // the same terminal-deamination signal the per-read LLR scores, so it is a prior, not
        // independent corroboration. Negative ⇒ absent ⇒ leave the profile default.
        if (pi_damaged_point >= 0.0f) p.pi_damaged = pi_damaged_point;
    }
};

namespace detail {
// Little-endian on the wire regardless of host byte order. On x86 (LE) the
// reverse is skipped, so bytes are byte-for-byte identical to a native memcpy
// and existing sidecars stay valid; on a big-endian host the field is swapped
// to/from LE so the on-disk format matches the documented contract.
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
constexpr bool host_is_little_endian = false;
#else
constexpr bool host_is_little_endian = true;
#endif
template <typename T>
inline void put(std::string& b, T v) {
    char buf[sizeof(T)];
    std::memcpy(buf, &v, sizeof(T));
    if (!host_is_little_endian)
        for (std::size_t i = 0; i < sizeof(T) / 2; ++i)
            std::swap(buf[i], buf[sizeof(T) - 1 - i]);
    b.append(buf, sizeof(T));
}
template <typename T>
inline bool get(const char*& p, const char* end, T& v) {
    if (p + sizeof(T) > end) return false;
    char buf[sizeof(T)];
    std::memcpy(buf, p, sizeof(T));
    if (!host_is_little_endian)
        for (std::size_t i = 0; i < sizeof(T) / 2; ++i)
            std::swap(buf[i], buf[sizeof(T) - 1 - i]);
    std::memcpy(&v, buf, sizeof(T));
    p += sizeof(T);
    return true;
}
}  // namespace detail

inline bool write_damage_summary(const std::string& path, const DamageSummary& s) {
    std::string b;
    detail::put(b, DamageSummary::MAGIC);
    detail::put(b, DamageSummary::VERSION);
    detail::put(b, s.d_max_combined);
    detail::put(b, s.d_max_5prime);
    detail::put(b, s.d_max_3prime);
    detail::put(b, s.lambda_5prime);
    detail::put(b, s.lambda_3prime);
    detail::put(b, s.fit_baseline_5prime);
    detail::put(b, s.fit_baseline_3prime);
    detail::put(b, s.terminal_shift_5prime);
    detail::put(b, s.stop_decay_llr_5prime);
    detail::put(b, s.d_max_cooc_point);
    detail::put(b, s.pi_point);
    detail::put(b, s.pi_shape_lambda);
    detail::put(b, s.pi_shape_baseline);
    detail::put(b, s.pi_shape_3p_lambda);    // v2
    detail::put(b, s.pi_shape_3p_baseline);  // v2
    detail::put(b, s.pi_damaged_point);      // v3
    detail::put(b, s.library_type);
    detail::put(b, s.damage_validated);
    detail::put(b, s.damage_artifact);
    detail::put(b, s.channel_b_valid);
    detail::put(b, s.pi_shape_fitted);
    detail::put(b, s.pi_shape_3p_fitted);    // v2
    detail::put(b, s.pi_state);
    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    if (!f) return false;
    f.write(b.data(), static_cast<std::streamsize>(b.size()));
    return static_cast<bool>(f);
}

inline std::optional<DamageSummary> read_damage_summary(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return std::nullopt;
    std::string b((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    const char* p = b.data();
    const char* end = p + b.size();
    uint32_t magic = 0, version = 0;
    if (!detail::get(p, end, magic) || magic != DamageSummary::MAGIC) return std::nullopt;
    // Version-tolerant: the layout is append-only, so a v2 file is a v3 file WITHOUT the
    // trailing pi_damaged_point float. Accept v2 or v3 (reject only <2 or >VERSION); read the
    // shared fields for both, and pull pi_damaged_point ONLY when version >= 3. A v2 file leaves
    // pi_damaged_point at its -1.0 default, so apply_to's `>= 0` guard skips it (existing .tdmg
    // keep loading; those libraries fall through to annotate's no_library/fallback tier).
    if (!detail::get(p, end, version) || version < 2 || version > DamageSummary::VERSION)
        return std::nullopt;
    DamageSummary s;
    bool ok = detail::get(p, end, s.d_max_combined)
           && detail::get(p, end, s.d_max_5prime)
           && detail::get(p, end, s.d_max_3prime)
           && detail::get(p, end, s.lambda_5prime)
           && detail::get(p, end, s.lambda_3prime)
           && detail::get(p, end, s.fit_baseline_5prime)
           && detail::get(p, end, s.fit_baseline_3prime)
           && detail::get(p, end, s.terminal_shift_5prime)
           && detail::get(p, end, s.stop_decay_llr_5prime)
           && detail::get(p, end, s.d_max_cooc_point)
           && detail::get(p, end, s.pi_point)
           && detail::get(p, end, s.pi_shape_lambda)
           && detail::get(p, end, s.pi_shape_baseline)
           && detail::get(p, end, s.pi_shape_3p_lambda)    // v2
           && detail::get(p, end, s.pi_shape_3p_baseline); // v2
    if (ok && version >= 3)
        ok = detail::get(p, end, s.pi_damaged_point);       // v3 (append-only trailing float)
    ok = ok
           && detail::get(p, end, s.library_type)
           && detail::get(p, end, s.damage_validated)
           && detail::get(p, end, s.damage_artifact)
           && detail::get(p, end, s.channel_b_valid)
           && detail::get(p, end, s.pi_shape_fitted)
           && detail::get(p, end, s.pi_shape_3p_fitted)    // v2
           && detail::get(p, end, s.pi_state);
    if (!ok) return std::nullopt;
    return s;
}

}  // namespace taph
