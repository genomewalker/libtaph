// A paired read whose R1 is shorter than the first provisional edge must put its LSD
// signature in the top stratum, the same bin LengthBinStats::update_paired uses for its
// aggregate terminal statistics.
#include "library_profile/worker.hpp"
#include <cstdio>
#include <numeric>

int main() {
    using namespace taph;
    const std::vector<int> edges = {40, 60};
    SigTables tables;
    tables.lsd_cnt.assign((edges.size() + 1) * N_LSD_SIG5, 0);
    WorkerState st;
    st.tables = &tables;
    st.lsd_on = true;
    st.lsd_prov_edges = edges;
    st.lbs.configure(edges);

    PairedWorkQueue queue(4);
    PairedBatch batch;
    batch.r1.push_back(std::string("TTAC") + std::string(31, 'G'));   // 35 bp, below the first edge
    batch.r2.push_back(std::string("AAGT") + std::string(76, 'C'));   // 80 bp
    queue.push(std::move(batch));
    queue.set_done();
    paired_worker_fn(queue, st);

    std::vector<long> per_bin(edges.size() + 1, 0);
    for (size_t b = 0; b < per_bin.size(); ++b)
        per_bin[b] = std::accumulate(tables.lsd_cnt.begin() + b * N_LSD_SIG5,
                                     tables.lsd_cnt.begin() + (b + 1) * N_LSD_SIG5, 0L);
    const bool ok = st.reads_scanned == 1 && per_bin[0] == 0 && per_bin[1] == 0 && per_bin[2] == 1;
    std::printf("reads_scanned=%lld per-bin LSD signatures: %ld %ld %ld -> %s\n",
                static_cast<long long>(st.reads_scanned), per_bin[0], per_bin[1], per_bin[2],
                ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
