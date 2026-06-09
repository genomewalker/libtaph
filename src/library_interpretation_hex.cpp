// Hexamer encoding, protocol tag lookup, end asymmetry.
#include "taph/library_interpretation.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>

namespace taph {
static constexpr double kExploratoryZCap = 12.0;
static constexpr double kDepurZCap       = kExploratoryZCap;

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

// ── Protocol-tag table ────────────────────────────────────────────────────────
// Empirical mapping derived from observed 5' top-hexamer fingerprints across
// labelled libraries. Each entry is a deterministic chemistry artefact, not a
// statistical pattern. log_lr = +4.0 (≈55:1 odds shift) reflects strong but
// not absolute evidence — leaves room for residual sample-mislabel risk.
//
// SS protocols (Santa-Cruz / Gansauge SS, SRSLY, splint-ligation variants):
//   TGTAGA  Santa-Cruz/Gansauge SS      (CL78/CL128 oligo footprint)
//   ACACTC  SRSLY                       (splint-ligation tag)
//   GTCACT  SS-variant                  (observed in SS-labelled libraries)
//   AACCCC  SS-variant                  (observed; tentative)
//
// DS protocols (TruSeq Y-adapter / NEBNext / similar):
//   CGATCT  TruSeq                      (Y-adapter remnant)
//   TGCTCT  NEBNext                     (Y-adapter variant)
//   TCAGTG  DS-UMI                      (end-repair + A-tail signature)
//   ACAATC  DS-variant
//   TCCGAT  DS-variant
//   CGCTCT  DS-variant
static const ProtocolTag kProtocolTags[] = {
    {"TGTAGA", SampleDamageProfile::LibraryType::SINGLE_STRANDED, 4.0f, "SantaCruz/Gansauge_SS"},
    {"ACACTC", SampleDamageProfile::LibraryType::SINGLE_STRANDED, 4.0f, "SRSLY"},
    {"GTCACT", SampleDamageProfile::LibraryType::SINGLE_STRANDED, 4.0f, "SS_variant"},
    {"AACCCC", SampleDamageProfile::LibraryType::SINGLE_STRANDED, 3.0f, "SS_variant_tentative"},
    {"CGATCT", SampleDamageProfile::LibraryType::DOUBLE_STRANDED, 4.0f, "TruSeq"},
    {"TGCTCT", SampleDamageProfile::LibraryType::DOUBLE_STRANDED, 4.0f, "NEBNext"},
    {"TCAGTG", SampleDamageProfile::LibraryType::DOUBLE_STRANDED, 4.0f, "DS_UMI"},
    {"ACAATC", SampleDamageProfile::LibraryType::DOUBLE_STRANDED, 3.0f, "DS_variant"},
    {"TCCGAT", SampleDamageProfile::LibraryType::DOUBLE_STRANDED, 3.0f, "DS_variant"},
    {"CGCTCT", SampleDamageProfile::LibraryType::DOUBLE_STRANDED, 3.0f, "DS_variant"},
};
static constexpr size_t kProtocolTagsCount = sizeof(kProtocolTags) / sizeof(kProtocolTags[0]);

const ProtocolTag* lookup_protocol_tag(const char* hex) {
    if (!hex || !hex[0]) return nullptr;
    for (size_t i = 0; i < kProtocolTagsCount; ++i)
        if (std::strcmp(kProtocolTags[i].hex, hex) == 0)
            return &kProtocolTags[i];
    return nullptr;
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

std::vector<HexEnrichment> compute_hex_enriched_3prime(
        const SampleDamageProfile& dp, float lfc_threshold) {
    std::vector<HexEnrichment> out;
    double tot_t = static_cast<double>(dp.n_hexamers_3prime);
    double tot_i = static_cast<double>(dp.n_hexamers_interior);
    if (tot_t <= 0 || tot_i <= 0) return out;
    for (int i = 0; i < 4096; ++i) {
        double t = dp.hexamer_count_3prime[i];
        double q = dp.hexamer_count_interior[i];
        if (t < 20.0 || q < 20.0) continue;
        double lfc = std::log2((t / tot_t + 1e-12) / (q / tot_i + 1e-12));
        if (lfc > static_cast<double>(lfc_threshold)) {
            auto seq = decode_hex(i);
            out.push_back({i, lfc, seq[5] == 'A'});  // damage_consistent = last base 'A' (G→A)
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
    r.top_enriched_3prime = compute_hex_enriched_3prime(dp, 1.5f);

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


// ── End hexamer asymmetry ─────────────────────────────────────────────────────

static int rc_hex_code(int h) {
    int rc = 0;
    for (int i = 0; i < 6; ++i) {
        int base = h & 3;
        h >>= 2;
        rc = (rc << 2) | (3 - base);
    }
    return rc;
}

static double jsd_pair(const double* p, const double* q, int n) {
    double jsd = 0.0;
    for (int i = 0; i < n; ++i) {
        double m = 0.5 * (p[i] + q[i]);
        if (p[i] > 0.0) jsd += 0.5 * p[i] * std::log2(p[i] / m);
        if (q[i] > 0.0) jsd += 0.5 * q[i] * std::log2(q[i] / m);
    }
    return std::max(0.0, std::min(1.0, jsd));
}

HexEndAsymmetry compute_hex_end_asymmetry(
        const SampleDamageProfile&        dp,
        const std::vector<HexEnrichment>& top5,
        const std::vector<HexEnrichment>& top3,
        int topk) {
    HexEndAsymmetry r;
    r.n_5prime = dp.n_hexamers_5prime;
    r.n_3prime = dp.n_hexamers_3prime;
    r.topk     = topk;

    double tot5 = static_cast<double>(dp.n_hexamers_5prime);
    double tot3 = static_cast<double>(dp.n_hexamers_3prime);
    double toti = static_cast<double>(dp.n_hexamers_interior);
    if (tot5 <= 0.0 || tot3 <= 0.0 || toti <= 0.0) {
        r.status = "insufficient_excess";
        return r;
    }

    static constexpr int N = 4096;
    double E5[N], E3[N];
    double sum_E5 = 0.0, sum_E3 = 0.0;
    for (int i = 0; i < N; ++i) {
        double p5 = dp.hexamer_count_5prime[i]   / tot5;
        double p3 = dp.hexamer_count_3prime[i]   / tot3;
        double pi = dp.hexamer_count_interior[i] / toti;
        E5[i] = std::max(0.0, p5 - pi);
        E3[i] = std::max(0.0, p3 - pi);
        sum_E5 += E5[i];
        sum_E3 += E3[i];
    }

    r.excess_mass_5prime = sum_E5;
    r.excess_mass_3prime = sum_E3;

    static constexpr double kMinExcess = 1e-4;
    if (sum_E5 < kMinExcess || sum_E3 < kMinExcess) {
        r.status = "insufficient_excess";
        return r;
    }

    double En5[N], En3[N];
    for (int i = 0; i < N; ++i) { En5[i] = E5[i] / sum_E5; En3[i] = E3[i] / sum_E3; }

    r.fwd_excess_jsd = jsd_pair(En5, En3, N);

    double Enc3[N] = {};
    for (int i = 0; i < N; ++i) Enc3[rc_hex_code(i)] = En3[i];

    r.rc_excess_jsd = jsd_pair(En5, Enc3, N);

    int k = topk;
    if ((int)top5.size() < k) k = (int)top5.size();
    if ((int)top3.size() < k) k = (int)top3.size();
    r.topk = k;
    if (k > 0) {
        int overlap = 0;
        for (int i = 0; i < k; ++i) {
            int c5 = top5[i].idx;
            for (int j = 0; j < k; ++j)
                if (rc_hex_code(top3[j].idx) == c5) { ++overlap; break; }
        }
        r.rc_overlap_topk = overlap;
    }

    r.status = "ok";
    return r;
}

} // namespace taph
