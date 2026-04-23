# libtaph

**Reference-free ancient DNA damage estimation and library-type classification from raw FASTQ reads.**

---

## What it does

libtaph estimates ancient-DNA terminal damage directly from raw reads, without a reference genome or read alignment. It measures six damage and fragmentation channels, fits exponential terminal-decay models, cross-validates apparent C→T damage against composition-robust stop-codon evidence, and reports an asymmetry-aware `d_max_combined`.

It also classifies each library as **double-stranded (DS)**, **single-stranded (SS)**, or **UNKNOWN** using a four-channel BIC classifier. `UNKNOWN` is the expected result when no damage model clearly beats the null.

---

## Quick start

```cpp
#include <taph/frame_selector_decl.hpp>

// One-shot: compute profile from a vector of sequences
std::vector<std::string> reads = { "ACGTCTAGCT...", ... };
taph::SampleDamageProfile profile =
    taph::FrameSelector::compute_sample_profile(reads);

// Streaming: accumulate reads incrementally
taph::SampleDamageProfile profile{};
for (const auto& seq : reads)
    taph::FrameSelector::update_sample_profile(profile, seq);
taph::FrameSelector::finalize_sample_profile(profile);

// Results
std::cout << "D_max (5'): " << profile.d_max_5prime << "\n";
std::cout << "D_max (3'): " << profile.d_max_3prime << "\n";
std::cout << "Library:    " << profile.library_type_str() << "\n";
// "double-stranded" | "single-stranded" | "unknown"
```

---

## Key features

| Feature | Description |
|---------|-------------|
| Reference-free | Works directly on FASTQ, no BAM, no alignment |
| D_max estimation | Calibrated, metaDMG-comparable damage rates |
| Library-type detection | BIC classifier: DS / SS / UNKNOWN |
| Multi-channel validation | 6 damage channels (A, B, B₃′, C, D, E) cross-validate signal |
| GC-stratified estimation | Separates ancient from modern DNA in mixed samples |
| Length-stratified estimation | Up to four length bins with BIC-selected edges; shared-component length × GC joint mixture |
| Context-aware C→T | Per-upstream-base amplitudes with dipyrimidine / CpG contrasts and chi-squared heterogeneity test |
| Streaming API | Incremental updates for memory-efficient processing |

---

## Validated performance

Tested on 315 Mediterranean sediment aDNA libraries (two independent datasets):

| Dataset | Correct | UNKNOWN | Wrong | Accuracy (determined) |
|---------|---------|---------|-------|-----------------------|
| Dataset 1 (91 samples) | 88 | 3 | 0 | **100%** |
| Dataset 2 (224 samples) | 193 | 28 | 3 | **98.5%** |

UNKNOWN: no detectable signal above the null model (zero-damage libraries where no library type can be inferred from sequence alone).

---

## Build

```cmake
# CMakeLists.txt
find_package(taph REQUIRED)
target_link_libraries(your_target PRIVATE taph)
```

Or as a CMake FetchContent dependency:

```cmake
include(FetchContent)
FetchContent_Declare(libtaph
    GIT_REPOSITORY https://github.com/genomewalker/libtaph.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(libtaph)
target_link_libraries(your_target PRIVATE taph)
```
