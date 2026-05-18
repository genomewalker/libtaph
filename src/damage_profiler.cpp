#include "taph/damage_profiler.hpp"
#include "taph/frame_selector_decl.hpp"
#include <algorithm>
#include <array>

namespace taph {

PrescanResult run_prescan(BatchFn          batch_fn,
                          SampleDamageProfile::LibraryType forced_lib,
                          int    min_read_len,
                          size_t prescan_reads) {
    SampleDamageProfile profile;
    profile.forced_library_type = forced_lib;

    std::array<uint32_t, 4096> hex3{};
    uint64_t n_hex3 = 0;
    size_t   done   = 0;

    std::vector<std::string> batch;
    bool more = true;
    while (more && done < prescan_reads) {
        batch.clear();
        size_t want = std::min(size_t(50000), prescan_reads - done);
        more = batch_fn(batch, want);
        for (const auto& seq : batch) {
            if ((int)seq.size() < min_read_len) continue;
            FrameSelector::update_sample_profile(profile, seq);
            int L = (int)seq.size();
            if (L >= 6) {
                int code = encode_hex_at(seq, L - 6);
                if (code >= 0) { ++hex3[code]; ++n_hex3; }
            }
            if (++done >= prescan_reads) break;
        }
    }

    FrameSelector::finalize_sample_profile(profile);
    AdapterStubs stubs = detect_adapter_stubs(profile, hex3.data(), n_hex3);
    return {std::move(profile), std::move(stubs)};
}

std::vector<uint32_t> encode_stub_codes(const std::vector<std::string>& stubs) {
    std::vector<uint32_t> codes;
    codes.reserve(stubs.size());
    for (const auto& s : stubs) {
        int code = encode_hex_at(s, 0);
        if (code >= 0) codes.push_back(static_cast<uint32_t>(code));
    }
    return codes;
}

std::string clip_adapters(const std::string&           seq,
                          const std::vector<uint32_t>& stub5_codes,
                          const std::vector<uint32_t>& stub3_codes,
                          int                          min_len) {
    std::string s = seq;

    if (!stub5_codes.empty() && (int)s.size() >= 6) {
        int code = encode_hex_at(s, 0);
        for (uint32_t sc : stub5_codes) {
            if (code == (int)sc) { s = s.substr(6); break; }
        }
    }
    if (!stub3_codes.empty() && (int)s.size() >= 6) {
        int L = (int)s.size();
        int code = encode_hex_at(s, L - 6);
        for (uint32_t sc : stub3_codes) {
            if (code == (int)sc) { s = s.substr(0, L - 6); break; }
        }
    }

    if ((int)s.size() < min_len) return {};
    return s;
}

} // namespace taph
