#include "taph/llr_table.hpp"
#include <algorithm>
#include <cmath>

namespace taph {

void llr_table_build(LLRTable* tables, int n_tables, const SampleDamageProfile& dp)
{
    for (int b = 0; b < n_tables; ++b) {
        auto& t = tables[b];
        const auto& gc = dp.gc_bins[b];
        float bg   = (gc.valid && gc.baseline_tc > 0.0f) ? gc.baseline_tc : dp.fit_baseline_5prime;
        float dmax = (gc.valid && gc.d_max > 0.005f)     ? gc.d_max       : dp.d_max_5prime;
        float lam  = dp.lambda_5prime;
        float pi   = (gc.valid && gc.p_damaged > 0.0f)   ? gc.p_damaged   : dp.pi_damaged;
        pi = std::clamp(pi, 0.001f, 0.999f);
        t.logit_pi = std::log(pi / (1.0f - pi));
        t.max_pos  = 0;
        for (int p = 0; p < 15; ++p) {
            float excess = dmax * std::exp(-lam * static_cast<float>(p));
            if (excess < 0.01f) { t.max_pos = p; break; }
            float pd = std::clamp(bg + (1.0f - bg) * excess, 1e-6f, 1.0f - 1e-6f);
            t.addT[p] = std::log(pd / bg);
            t.addC[p] = std::log((1.0f - pd) / (1.0f - bg));
            t.max_pos  = p + 1;
        }
    }
}

float llr_table_weight_fwd(const char* seq, int len, const LLRTable& t)
{
    int mp = std::min(t.max_pos, len);
    float llr = 0.0f;
    for (int p = 0; p < mp; ++p) {
        char b = static_cast<char>(seq[p] & ~0x20u);
        if      (b == 'T') llr += t.addT[p];
        else if (b == 'C') llr += t.addC[p];
    }
    return 1.0f / (1.0f + std::exp(-(t.logit_pi + llr)));
}

float llr_table_weight_rev(const char* seq, int len, const LLRTable& t)
{
    int mp = std::min(t.max_pos, len);
    float llr = 0.0f;
    for (int p = 0; p < mp; ++p) {
        char b = static_cast<char>(seq[len - 1 - p] & ~0x20u);
        if      (b == 'A') llr += t.addT[p];   // A in original → T in RC
        else if (b == 'G') llr += t.addC[p];   // G in original → C in RC
    }
    return 1.0f / (1.0f + std::exp(-(t.logit_pi + llr)));
}

} // namespace taph
