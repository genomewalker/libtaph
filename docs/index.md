# libdart-damage

**Reference-free ancient DNA damage estimation and library-type classification from raw FASTQ reads.**

---

## What it does

libdart-damage estimates cytosine deamination damage in ancient DNA without a reference genome or read alignment. It scans terminal positions of raw FASTQ sequences for C→T (5' end) and G→A (3' end) substitution patterns and fits an exponential decay model to each end.

It also classifies each library as **double-stranded (DS)** or **single-stranded (SS)** using a 4-channel joint BIC model, so downstream tools can apply library-appropriate deduplication and damage masking.

---

## Quick start

```cpp
#include <dart/frame_selector_decl.hpp>

// One-shot: compute profile from a vector of sequences
std::vector<std::string> reads = { "ACGTCTAGCT...", ... };
dart::SampleDamageProfile profile =
    dart::FrameSelector::compute_sample_profile(reads);

// Streaming: accumulate reads incrementally
dart::SampleDamageProfile profile;
dart::FrameSelector::reset_sample_profile(profile);
for (const auto& seq : reads)
    dart::FrameSelector::update_sample_profile(profile, seq);
dart::FrameSelector::finalize_sample_profile(profile);

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
find_package(dart-damage REQUIRED)
target_link_libraries(your_target PRIVATE dart-damage)
```

Or as a CMake FetchContent dependency:

```cmake
include(FetchContent)
FetchContent_Declare(libdart_damage
    GIT_REPOSITORY https://github.com/genomewalker/libdart-damage.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(libdart_damage)
target_link_libraries(your_target PRIVATE dart-damage)
```
