// library_interpretation.cpp — all damage interpretation logic for taph.
// Pure functions of taph::SampleDamageProfile; no FASTQ I/O.

#include "taph/library_interpretation.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>

namespace taph {

// ── Hexamer utilities ─────────────────────────────────────────────────────────

std::array<char,7> decode_hex(int code) {
    const char bases[] = "ACGT";
    std::array<char,7> s{}; s[6] = '\0';
    if (code < 0 || code >= 4096) { s[0] = '\0'; return s; }
    for (int i = 5; i >= 0; --i) { s[i] = bases[code & 3]; code >>= 2; }
    return s;
}

int encode_hex_at(const std::string& s, int pos) {
    if (pos < 0 || static_cast<size_t>(pos) + 6 > s.size()) return -1;
    int code = 0;
    for (int i = pos; i < pos + 6; ++i) {
        int b;
        switch (static_cast<unsigned char>(s[i]) | 0x20) {  // tolower via bit
            case 'a': b = 0; break;
            case 'c': b = 1; break;
            case 'g': b = 2; break;
            case 't': b = 3; break;
            default: return -1;
        }
        code = (code << 2) | b;
    }
    return code;
}

std::vector<HexEnrichment> compute_hex_enriched_5prime(
        const SampleDamageProfile& dp, float lfc_threshold) {
    std::vector<HexEnrichment> out;
    double tot_t = static_cast<double>(dp.n_hexamers_5prime);
    double tot_i = static_cast<double>(dp.n_hexamers_interior);
    if (tot_t <= 0 || tot_i <= 0) return out;
    for (int i = 0; i < 4096; ++i) {
        double t = dp.hexamer_count_5prime[i];
        double q = dp.hexamer_count_interior[i];
        if (t < 20.0 || q < 20.0) continue;
        double lfc = std::log2((t / tot_t + 1e-12) / (q / tot_i + 1e-12));
        if (lfc > static_cast<double>(lfc_threshold)) {
            auto seq = decode_hex(i);
            out.push_back({i, lfc, seq[0] == 'T'});
        }
    }
    std::sort(out.begin(), out.end(),
              [](const HexEnrichment& a, const HexEnrichment& b) {
                  if (a.log2fc != b.log2fc) return a.log2fc > b.log2fc;
                  return a.idx < b.idx;
              });
    return out;
}

AdapterStubs detect_adapter_stubs(
        const SampleDamageProfile& dp,
        const uint32_t hex3_terminal[4096],
        uint64_t n_hex3) {

    AdapterStubs r;

    // 5' stubs: non-T leading, lfc>3.0 vs interior, cap at 5
    auto enriched = compute_hex_enriched_5prime(dp, 1.5f);

    if (!enriched.empty()) {
        auto top = decode_hex(enriched[0].idx);
        if (top[0] != 'T' && enriched[0].log2fc > 1.5)
            r.flag_hex_artifact = true;
    }

    if (r.flag_hex_artifact) {
        for (const auto& hr : enriched) {
            if ((int)r.stubs5.size() >= 5) break;
            auto s = decode_hex(hr.idx);
            if (s[0] != 'T' && hr.log2fc > 3.0)
                r.stubs5.push_back(std::string(s.data(), 6));
        }
        if (!r.stubs5.empty()) r.adapter_clipped = true;
    }

    r.top_enriched = std::move(enriched);

    // 3' stubs: only when 5' adapter stubs exist (gating avoids false positives
    // on SS libraries where 3' enrichment is genuine G→A signal)
    if (r.adapter_clipped && n_hex3 > 0) {
        double tot_t3 = static_cast<double>(n_hex3);
        double tot_i  = static_cast<double>(dp.n_hexamers_interior);
        if (tot_i > 0) {
            struct H3 { double lfc; int idx; };
            std::vector<H3> hex3_enriched;
            for (int i = 0; i < 4096; ++i) {
                double t = static_cast<double>(hex3_terminal[i]);
                double q = static_cast<double>(dp.hexamer_count_interior[i]);
                if (t < 20.0 || q < 20.0) continue;
                double lfc = std::log2((t / tot_t3 + 1e-12) / (q / tot_i + 1e-12));
                if (lfc > 3.0) hex3_enriched.push_back({lfc, i});
            }
            std::sort(hex3_enriched.begin(), hex3_enriched.end(),
                      [](const H3& a, const H3& b) {
                          if (a.lfc != b.lfc) return a.lfc > b.lfc;
                          return a.idx < b.idx;
                      });
            for (const auto& hr : hex3_enriched) {
                if ((int)r.stubs3.size() >= 5) break;
                auto s = decode_hex(hr.idx);
                // last base 'A' excluded: G→A deamination enriches A at genuine 3' termini
                if (s[5] != 'A')
                    r.stubs3.push_back(std::string(s.data(), 6));
            }
            if (!r.stubs3.empty()) r.adapter3_clipped = true;
        }
    }

    return r;
}

HexStats compute_hex_stats(const SampleDamageProfile& dp) {
    HexStats r;
    double tot_t = static_cast<double>(dp.n_hexamers_5prime);
    double tot_i = static_cast<double>(dp.n_hexamers_interior);

    // Entropy + JSD
    if (tot_t > 0 && tot_i > 0) {
        double h_mix = 0.0;
        for (int i = 0; i < 4096; ++i) {
            double p = dp.hexamer_count_5prime[i]  / tot_t;
            double q = dp.hexamer_count_interior[i] / tot_i;
            double m = 0.5 * (p + q);
            if (p > 0) r.entropy_terminal -= p * std::log2(p);
            if (q > 0) r.entropy_interior -= q * std::log2(q);
            if (m > 0) h_mix              -= m * std::log2(m);
        }
        r.jsd = h_mix - 0.5 * (r.entropy_terminal + r.entropy_interior);
    }

    // Multinomial G-test (terminal vs interior)
    if (tot_t > 0 && tot_i > 0) {
        double N_hex = tot_t + tot_i;
        double G_stat = 0.0;
        int k_eff = 0;
        for (int i = 0; i < 4096; ++i) {
            double c = static_cast<double>(dp.hexamer_count_5prime[i]);
            double d = static_cast<double>(dp.hexamer_count_interior[i]);
            if (c + d < 5.0) continue;
            double Ec = tot_t * (c + d) / N_hex;
            double Ed = tot_i * (c + d) / N_hex;
            if (c > 0) G_stat += 2.0 * c * std::log(c / Ec);
            if (d > 0) G_stat += 2.0 * d * std::log(d / Ed);
            ++k_eff;
        }
        r.shift_g = G_stat;
        if (k_eff > 1) {
            double df = static_cast<double>(k_eff - 1);
            r.shift_z = (G_stat - df) / std::sqrt(2.0 * df);
            r.shift_p = 0.5 * std::erfc(r.shift_z / std::sqrt(2.0));
        }
    }
    return r;
}

// ── Scores ────────────────────────────────────────────────────────────────────

CpgScore compute_cpg_score(const SampleDamageProfile& dp) {
    CpgScore r;
    if (std::isnan(dp.log2_cpg_ratio)) return r;
    double ecpg  = static_cast<double>(dp.effcov_ct5_cpg_like_terminal);
    double encpg = static_cast<double>(dp.effcov_ct5_noncpg_like_terminal);
    double se = std::sqrt(1.0 / (ecpg + 1.0) + 1.0 / (encpg + 1.0));
    if (se > 1e-9) {
        r.z = static_cast<double>(dp.log2_cpg_ratio) / se;
        r.p = std::erfc(std::abs(r.z) / std::sqrt(2.0));
    }
    return r;
}

OxogInteriorScore compute_oxog_interior_score(const SampleDamageProfile& dp) {
    OxogInteriorScore r;
    static constexpr int RC4[4] = {3, 2, 1, 0};
    double sc = 0.0, vr = 0.0;
    for (const auto* tri : {&dp.tri_5prime_interior, &dp.tri_3prime_interior}) {
        for (int p = 0; p < 4; ++p) {
            for (int n = 0; n < 4; ++n) {
                double k  = (*tri)[p*16 + 3*4 + n];
                double g  = (*tri)[p*16 + 2*4 + n];
                double nt = k + g;
                if (nt < 10.0) continue;
                int rp = RC4[n], rn = RC4[p];
                double a_rc = (*tri)[rp*16 + 0*4 + rn];
                double c_rc = (*tri)[rp*16 + 1*4 + rn];
                double ca   = a_rc + c_rc;
                if (ca < 10.0) continue;
                double theta = a_rc / ca;
                sc += k - nt * theta;
                vr += nt * theta * (1.0 - theta);
            }
        }
    }
    if (vr > 0.0) {
        r.z = sc / std::sqrt(vr);
        r.p = 0.5 * std::erfc(r.z / std::sqrt(2.0));
    }
    return r;
}

OxogTrinucResult compute_oxog_trinuc(const SampleDamageProfile& dp) {
    static constexpr double OXOG_REF[16] = {
        0.0512, 0.0388, 0.0257, 0.0601,
        0.0788, 0.0621, 0.0302, 0.1124,
        0.0389, 0.0296, 0.0198, 0.0453,
        0.0863, 0.0682, 0.0388, 0.1138,
    };
    static constexpr int RC4[4] = {3, 2, 1, 0};
    OxogTrinucResult r;
    double v[16] = {};
    for (const auto* tri : {&dp.tri_5prime_interior, &dp.tri_3prime_interior}) {
        for (int p = 0; p < 4; ++p) {
            for (int n = 0; n < 4; ++n) {
                double k  = static_cast<double>((*tri)[p*16 + 3*4 + n]);
                double g  = static_cast<double>((*tri)[p*16 + 2*4 + n]);
                double nt = k + g;
                if (nt < 10.0) continue;
                int rp = RC4[n], rn = RC4[p];
                double a_rc = static_cast<double>((*tri)[rp*16 + 0*4 + rn]);
                double c_rc = static_cast<double>((*tri)[rp*16 + 1*4 + rn]);
                double ca = a_rc + c_rc;
                if (ca < 10.0) continue;
                double theta = a_rc / ca;
                v[p*4 + n] += std::max(0.0, k - nt * theta);
                ++r.n_ctx;
            }
        }
    }
    double dot = 0.0, nv = 0.0, nr = 0.0;
    for (int i = 0; i < 16; ++i) {
        dot += v[i] * OXOG_REF[i];
        nv  += v[i] * v[i];
        nr  += OXOG_REF[i] * OXOG_REF[i];
    }
    if (nv > 1e-30 && nr > 0)
        r.cosine = dot / (std::sqrt(nv) * std::sqrt(nr));
    return r;
}

DepurScore compute_depur_score(const SampleDamageProfile& dp, bool is_ss) {
    DepurScore r;
    r.z5     = static_cast<double>(dp.ctrl_z_5prime);
    r.shift5 = static_cast<double>(dp.purine_enrichment_5prime);
    r.shift3 = static_cast<double>(dp.purine_enrichment_3prime);

    double p5 = 0.5 * std::erfc(r.z5 / std::sqrt(2.0));
    if (!is_ss) {
        r.z3      = static_cast<double>(dp.ctrl_z_3prime);
        double p3 = 0.5 * std::erfc(r.z3 / std::sqrt(2.0));
        r.z = std::min(r.z5, r.z3);
        r.p = std::max(p5, p3);
    } else {
        r.z = r.z5;
        r.p = p5;
    }
    return r;
}

// ── Damage mask ───────────────────────────────────────────────────────────────

DamageMask compute_damage_mask(const SampleDamageProfile& dp,
                                bool is_ss, double threshold, int min_cov) {
    DamageMask r;
    const double bg_5  = dp.fit_baseline_5prime;
    const double bg_3  = dp.fit_baseline_3prime;
    const double bg_tc = dp.baseline_t_freq /
                         (dp.baseline_t_freq + dp.baseline_c_freq + 1e-9);

    for (int p = 0; p < INTERP_N_POS; ++p) {
        double excess_5 = 0.0, excess_3 = 0.0;
        if (dp.tc_total_5prime[p] >= static_cast<double>(min_cov))
            excess_5 = dp.t_freq_5prime[p] - bg_5;
        if (is_ss) {
            if (dp.tc_total_3prime[p] >= static_cast<double>(min_cov))
                excess_3 = dp.t_freq_3prime[p] / dp.tc_total_3prime[p] - bg_tc;
        } else {
            if (dp.ag_total_3prime[p] >= static_cast<double>(min_cov))
                excess_3 = dp.a_freq_3prime[p] - bg_3;
        }
        r.pos[p] = (excess_5 > threshold) || (excess_3 > threshold);
    }
    for (int p = 0; p < INTERP_N_POS; ++p) {
        if (!r.pos[p]) continue;
        if (r.n_masked) r.masked_str += ',';
        r.masked_str += std::to_string(p);
        ++r.n_masked;
    }
    return r;
}

// ── Library QC flags ──────────────────────────────────────────────────────────

LibraryQcFlags compute_library_qc_flags(
        const SampleDamageProfile& dp,
        bool is_ss,
        bool flag_hex_artifact,
        double jsd,
        double h_term,
        double short_read_frac) {
    LibraryQcFlags r;

    double tot_term = static_cast<double>(dp.n_hexamers_5prime);

    r.adapter_remnant_5prime   = dp.fit_offset_5prime > 1 || dp.position_0_artifact_5prime;
    r.adapter_remnant_3prime   = dp.fit_offset_3prime > 1 || dp.position_0_artifact_3prime;
    r.hexamer_composition_bias = tot_term > 0 && h_term < 4.5;
    r.hexamer_terminal_shift   = tot_term > 0 && dp.n_hexamers_interior > 0 && jsd > 0.05;
    r.short_read_spike         = short_read_frac > 0.20;
    r.depurination             = dp.depurination_detected;
    r.ds_3prime_signal_absent  = !is_ss
                                  && dp.d_max_5prime > 0.05f
                                  && dp.d_max_3prime < 0.03f;
    r.hexamer_artifact_bias    = flag_hex_artifact;

    // Inward-displaced G→A: pos0 depleted while pos1-3 elevated (trimmer artifact)
    if (!is_ss && dp.d_max_5prime > 0.05f) {
        constexpr double MIN_GA = 100.0;
        double bg_n = 0.0, bg_d = 0.0;
        for (int p = 8; p < 15; ++p) {
            bg_n += dp.a_freq_3prime[p];
            bg_d += dp.ag_total_3prime[p];
        }
        double bg = (bg_d >= 500.0) ? bg_n / bg_d : 0.47;
        double ga0 = (dp.ag_total_3prime[0] >= MIN_GA)
            ? dp.a_freq_3prime[0] / dp.ag_total_3prime[0] : -1.0;
        double ga_peak = 0.0;
        for (int p = 1; p <= 3; ++p)
            if (dp.ag_total_3prime[p] >= MIN_GA)
                ga_peak = std::max(ga_peak,
                                   dp.a_freq_3prime[p] / dp.ag_total_3prime[p]);
        if (ga0 >= 0.0)
            r.ga3_inward_displaced = ((ga0 - bg) < -0.02) && ((ga_peak - bg) > 0.02);
    }
    return r;
}

// ── Preservation / authenticity ───────────────────────────────────────────────

PreservationSummary compute_preservation_summary(
        const SampleDamageProfile& dp,
        bool   is_ss,
        bool   adapter_clipped,
        bool   flag_hex_artifact,
        double cpg_score_z,
        double oxog_score_z,
        double oxog_trinuc_cosine,
        double hex_shift_p) {

    PreservationSummary r;

    auto clamp01 = [](double x) -> double { return std::clamp(x, 0.0, 1.0); };
    auto sat = [](double x, double half_sat) -> double {
        x = std::max(0.0, x);
        return x / (x + half_sat);
    };
    auto z_evidence = [](double z, double scale = 8.0) -> double {
        return 1.0 - std::exp(-std::max(0.0, z) / scale);
    };
    auto p_evidence = [](double p) -> double {
        if (std::isnan(p)) return 0.0;
        p = std::clamp(p, 1e-300, 1.0);
        double x = -std::log10(p);
        return x / (x + 3.0);
    };

    // d5_hexamer_corrected: when adapter artifact suppresses pos0 C→T,
    // extrapolate the decay curve back: d5_corr = d5 * exp(λ × fit_offset).
    double d5_raw  = static_cast<double>(dp.d_max_5prime);
    double d5_corr = d5_raw;
    if (!is_ss && !adapter_clipped && flag_hex_artifact
            && dp.position_0_artifact_5prime
            && dp.lambda_5prime > 0.0f && dp.fit_offset_5prime >= 1) {
        d5_corr = std::min(d5_raw * std::exp(static_cast<double>(dp.lambda_5prime)
                                              * dp.fit_offset_5prime), 0.50);
    }
    r.d5_raw               = d5_raw;
    r.d5_hexamer_corrected = d5_corr;
    r.d5_was_corrected     = (d5_corr != d5_raw);

    double a5 = sat(d5_corr, 0.05);
    double a3 = sat(static_cast<double>(dp.d_max_3prime), 0.05);
    double w5 = is_ss ? 1.0 : 2.0, w3 = is_ss ? 1.0 : 0.5;

    double authenticity_eff;
    if (dp.mixture_converged) {
        double pi = clamp01(static_cast<double>(dp.mixture_pi_ancient));
        authenticity_eff = clamp01(0.20 * a5 + 0.80 * pi);
    } else {
        authenticity_eff = clamp01((w5*a5 + w3*a3) / (w5 + w3));
    }
    r.authenticity_eff = authenticity_eff;

    // authenticity_evidence
    double mix_term = 0.0, wmix = 0.0;
    if (dp.mixture_converged && dp.mixture_identifiable) {
        mix_term = sat(static_cast<double>(dp.mixture_d_ancient), 0.05) *
                   clamp01(static_cast<double>(dp.mixture_pi_ancient));
        wmix = 1.0;
    }
    double cpg_cov = static_cast<double>(dp.effcov_ct5_cpg_like_terminal) +
                     static_cast<double>(dp.effcov_ct5_noncpg_like_terminal);
    double wcpg = cpg_cov / (cpg_cov + 2000.0);

    double z5e = z_evidence(static_cast<double>(dp.terminal_z_5prime));
    double z3e = is_ss ? 0.0 : z_evidence(static_cast<double>(dp.terminal_z_3prime));
    double mix_e = 0.0;
    if (dp.mixture_converged && dp.mixture_identifiable) {
        mix_e = mix_term;
    } else if (dp.mixture_converged) {
        double mix_mag = sat(static_cast<double>(dp.mixture_d_ancient), 0.05) *
                         clamp01(static_cast<double>(dp.mixture_pi_ancient));
        mix_e = 0.5 * mix_mag;
    }
    double cpg_e = wcpg * z_evidence(std::max(0.0, cpg_score_z), 5.0);
    double wmix_e = dp.mixture_converged ? 1.0 : 0.0;
    double auth_n_terms = 1.0 + (is_ss ? 0.0 : 1.0) + wmix_e + wcpg;
    r.authenticity_evidence = (auth_n_terms > 1e-9)
        ? clamp01((z5e + z3e + mix_e + cpg_e) / auth_n_terms)
        : 0.0;

    // oxidation_eff / oxidation_evidence
    bool has_trinuc = !std::isnan(oxog_trinuc_cosine);
    double ox_shape = has_trinuc
        ? clamp01((oxog_trinuc_cosine - 0.65) / 0.20)
        : 0.0;
    r.oxidation_eff = clamp01(
        sat(static_cast<double>(dp.ox_d_max), 0.02) *
        (has_trinuc ? ox_shape : 1.0));
    double oz = (!is_ss && static_cast<double>(dp.ox_d_max) > 0.01) ? oxog_score_z : 0.0;
    r.oxidation_evidence = has_trinuc
        ? clamp01(ox_shape * 0.5 * (z_evidence(oz) + 1.0))
        : clamp01(z_evidence(oz));

    // qc_risk_eff / qc_evidence
    double d5_eff = static_cast<double>(dp.d_max_5prime);
    double d3_eff = static_cast<double>(dp.d_max_3prime);
    bool adapter_5 = (dp.fit_offset_5prime > 1) ||
                     (dp.position_0_artifact_5prime && d5_eff < 0.05 && !is_ss);
    bool adapter_3 = (dp.fit_offset_3prime > 1) ||
                     (dp.position_0_artifact_3prime && d3_eff < 0.05);
    double sig_eff     = is_ss ? d3_eff : d5_eff;
    double hex_dam_disc = sat(sig_eff, 0.05);
    double qhex  = sat(std::abs(static_cast<double>(dp.hexamer_excess_tc)), 0.05)
                   * (1.0 - hex_dam_disc);
    double qart  = dp.damage_artifact ? 1.0 : 0.0;
    double qadapt = (adapter_5 || adapter_3) ? 1.0 : 0.0;
    r.qc_risk_eff = clamp01((2.0*qhex + qart + qadapt) / 4.0);

    double hex_ev = qhex * p_evidence(hex_shift_p);
    r.qc_evidence = clamp01(std::max({hex_ev, qart, qadapt * 0.75}));

    // label
    if      (authenticity_eff < 0.10) r.label = "modern-like";
    else if (authenticity_eff < 0.30) r.label = "weak";
    else if (authenticity_eff < 0.55) r.label = "moderate";
    else                               r.label = "ancient";

    return r;
}

// ── Damage-context profile ────────────────────────────────────────────────────

namespace {

constexpr float kDeamNorm      = 0.10f;  // d_max at which deamination score ≈ 0.63
constexpr float kDipyrNorm     = 0.05f;  // dipyr contrast at which score saturates
                                          // (dp.dipyr_contrast uses 0.5 factors)
constexpr float kOxidativeNorm = 0.05f;  // oxidative signal at which score saturates
constexpr float kFragNorm      = 0.15f;  // purine enrichment at which score saturates
constexpr uint64_t kMinReads   = 1000;   // SampleDamageProfile::is_valid() threshold

inline float clamp01f(float x) {
    if (std::isnan(x)) return x;
    return x < 0.0f ? 0.0f : (x > 1.0f ? 1.0f : x);
}
inline float sigmoidf(float x) { return 1.0f / (1.0f + std::exp(-x)); }

} // namespace

const char* to_string(DamageContextProfile::DominantProcess p) {
    using D = DamageContextProfile::DominantProcess;
    switch (p) {
        case D::None:                   return "none";
        case D::LowDamage:              return "low_damage";
        case D::CytosineDeamination:    return "cytosine_deamination";
        case D::CpgEnrichedDeamination: return "cpg_enriched_deamination";
        case D::DipyrimidineBiased:     return "dipyrimidine_biased";
        case D::OxidativeLike:          return "oxidative_like";
        case D::FragmentationBias:      return "fragmentation_bias";
        case D::LibraryArtifactLikely:  return "library_artifact_likely";
    }
    return "none";
}

DamageContextProfile compute_damage_context_profile(
        const SampleDamageProfile& dp,
        double cpg_z,
        double hex_shift_z,
        bool   adapter_clipped,
        bool   adapter3_clipped,
        bool   flag_hex_artifact) {

    DamageContextProfile r;

    // Evidence block: raw underlying numbers for downstream auditing.
    r.evidence.d_max_5                    = dp.d_max_5prime;
    r.evidence.d_max_3                    = dp.d_max_3prime;
    r.evidence.lambda_5                   = dp.lambda_5prime;
    r.evidence.lambda_3                   = dp.lambda_3prime;
    r.evidence.log2_cpg_ratio             = dp.log2_cpg_ratio;
    r.evidence.cpg_z                      = static_cast<float>(cpg_z);
    r.evidence.ox_gt_asymmetry            = dp.ox_gt_asymmetry;
    r.evidence.purine_enrichment_5prime   = dp.purine_enrichment_5prime;
    r.evidence.hex_shift_z                = static_cast<float>(hex_shift_z);
    r.evidence.adapter_clipped            = adapter_clipped;
    r.evidence.adapter3_clipped           = adapter3_clipped;
    r.evidence.flag_hex_artifact          = flag_hex_artifact;
    r.evidence.position_0_artifact_5prime = dp.position_0_artifact_5prime;
    r.evidence.position_0_artifact_3prime = dp.position_0_artifact_3prime;
    r.evidence.fit_offset_5prime          = dp.fit_offset_5prime;
    r.evidence.fit_offset_3prime          = dp.fit_offset_3prime;
    r.evidence.n_reads                    = dp.n_reads;

    // 8-oxoG 16-context panel summary (mean/max of the per-context signal).
    double s_sum = 0.0; float s_max = 0.0f; int n_ctx = 0;
    for (int i = 0; i < SampleDamageProfile::N_OXOG16; ++i) {
        float v = dp.s_oxog_16ctx[i];
        if (!std::isfinite(v)) continue;
        s_sum += v;
        if (v > s_max) s_max = v;
        ++n_ctx;
    }
    r.evidence.s_oxog_mean = n_ctx > 0 ? static_cast<float>(s_sum / n_ctx) : 0.0f;
    r.evidence.s_oxog_max  = s_max;

    // Use the already-finalized contrast from SampleDamageProfile:
    //   dipyr_contrast = 0.5*(CC + TC) - 0.5*(AC + GC).
    r.evidence.dipyr_contrast = dp.dipyr_contrast;

    const bool any_art_flag = flag_hex_artifact || adapter_clipped ||
                              adapter3_clipped ||
                              dp.position_0_artifact_5prime ||
                              dp.position_0_artifact_3prime ||
                              dp.fit_offset_5prime > 1 ||
                              dp.fit_offset_3prime > 1;
    const auto NaN = std::numeric_limits<float>::quiet_NaN();
    auto finite_max2 = [](float a, float b) {
        bool fa = std::isfinite(a), fb = std::isfinite(b);
        if (fa && fb) return std::max(a, b);
        if (fa)       return a;
        if (fb)       return b;
        return std::numeric_limits<float>::quiet_NaN();
    };

    // Six scores. Computed unconditionally so that low-coverage samples still
    // surface any evaluable raw signal; dominant_process is set to None below.

    // terminal deamination: 1 - exp(-max(d5, d3) / kDeamNorm). NaN-safe.
    float dmax_term = finite_max2(dp.d_max_5prime, dp.d_max_3prime);
    r.terminal_deamination_score = std::isnan(dmax_term)
        ? NaN
        : clamp01f(1.0f - std::exp(-dmax_term / kDeamNorm));

    // CpG context: sigmoid on the CpG/non-CpG z-score from compute_cpg_score.
    r.cpg_context_score = (!std::isfinite(dp.log2_cpg_ratio) || !std::isfinite(cpg_z))
        ? NaN
        : sigmoidf(static_cast<float>(cpg_z));

    // Dipyrimidine context: normalized upstream-context excess.
    r.dipyrimidine_context_score = !std::isfinite(r.evidence.dipyr_contrast)
        ? NaN
        : clamp01f(r.evidence.dipyr_contrast / kDipyrNorm);

    // Oxidative context: max of strand-asymmetric G->T signal and mean NGN panel.
    float ox_signed = finite_max2(std::fabs(dp.ox_gt_asymmetry),
                                  r.evidence.s_oxog_mean);
    r.oxidative_context_score = std::isnan(ox_signed)
        ? NaN
        : clamp01f(ox_signed / kOxidativeNorm);

    // Fragmentation context: purine enrichment at read starts vs interior.
    r.fragmentation_context_score = !std::isfinite(dp.purine_enrichment_5prime)
        ? NaN
        : clamp01f(dp.purine_enrichment_5prime / kFragNorm);

    // Library artifact: max of flags and a sigmoid on the composition shift z.
    // If hex_shift_z is NaN and no flags are set the signal is unknown; the
    // score stays NaN rather than silently reading as 0.
    float art_flag = any_art_flag ? 1.0f : 0.0f;
    float art_sig  = std::isfinite(hex_shift_z)
        ? sigmoidf(static_cast<float>(hex_shift_z) - 4.0f)
        : NaN;
    if (any_art_flag) {
        r.library_artifact_score = clamp01f(std::max(art_flag,
                                     std::isnan(art_sig) ? 0.0f : art_sig));
    } else if (std::isnan(art_sig)) {
        r.library_artifact_score = NaN;
    } else {
        r.library_artifact_score = clamp01f(art_sig);
    }

    // Deterministic dominant-process rule. Unevaluable signals (NaN) are not
    // allowed to trigger a branch; a NaN terminal-deamination score in
    // particular must not collapse to low_damage.
    using D = DamageContextProfile::DominantProcess;
    float td  = r.terminal_deamination_score;
    float cpg = r.cpg_context_score;
    float dip = r.dipyrimidine_context_score;
    float ox  = r.oxidative_context_score;
    float fr  = r.fragmentation_context_score;
    float art = r.library_artifact_score;
    auto gt = [](float x, float t){ return std::isfinite(x) && x >  t; };
    auto lt = [](float x, float t){ return std::isfinite(x) && x <  t; };

    if (dp.n_reads < kMinReads) {
        r.dominant_process = D::None;
        r.interpretation   = "insufficient coverage for damage-context scoring";
    } else if (std::isnan(td)) {
        r.dominant_process = D::None;
        r.interpretation   = "terminal deamination signal not evaluable";
    } else if (any_art_flag && lt(td, 0.5f)) {
        // Rule fires when a boolean flag is set (hex-artifact detector, adapter
        // stub, position-0 artifact, or fit-offset adapter remnant) AND the
        // terminal deamination signal is not dominating. The score itself
        // (library_artifact_score) saturates at 1.0 whenever a flag is present,
        // so `gt(art, 0.7)` would be tautological; the meaningful second gate
        // is "damage is not the primary driver", i.e. td < 0.5. When strong
        // genuine damage is present alongside adapter evidence, the damage
        // label wins (artifact booleans remain visible in evidence).
        r.dominant_process = D::LibraryArtifactLikely;
        r.interpretation = "composition or adapter-stub evidence dominates over damage signal";
    } else if (gt(fr, 0.5f) && lt(td, 0.3f)) {
        r.dominant_process = D::FragmentationBias;
        r.interpretation = "purine enrichment at fragment starts without matching terminal deamination";
    } else if (lt(td, 0.10f)) {
        r.dominant_process = D::LowDamage;
        r.interpretation = "terminal deamination signal below detection threshold";
    } else if (gt(cpg, 0.7f) && gt(td, 0.3f)
               && gt(r.evidence.log2_cpg_ratio, 0.15f)) {
        // CpG z-scores can exceed the sigmoid threshold at modest effect sizes
        // on large samples; require a minimum log2 ratio so the label reflects
        // a real CpG-vs-non-CpG enrichment, not just statistical significance.
        r.dominant_process = D::CpgEnrichedDeamination;
        r.interpretation = "terminal deamination with elevated CpG-context contribution "
                           "consistent with methylated-cytosine deamination";
    } else if (gt(ox, 0.5f) && lt(td, 0.5f)) {
        r.dominant_process = D::OxidativeLike;
        r.interpretation = "strand-asymmetric G->T / C->A excess consistent with oxidative damage";
    } else if (gt(dip, 0.4f)) {
        r.dominant_process = D::DipyrimidineBiased;
        r.interpretation = "dipyrimidine upstream-context excess in terminal C->T rates";
    } else {
        r.dominant_process = D::CytosineDeamination;
        r.interpretation = "terminal C->T / G->A enrichment consistent with post-mortem cytosine deamination";
    }
    r.dominant_process_str = to_string(r.dominant_process);
    return r;
}

} // namespace taph
