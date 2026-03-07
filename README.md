# libdart-damage

**Reference-free ancient DNA damage estimation and library-type classification from raw FASTQ reads.**

Ancient DNA carries a chemical fingerprint of time: cytosine deamination converts C→T at 5' read ends and G→A at 3' ends, with rates decaying exponentially away from each terminus. libdart-damage quantifies this damage directly from FASTQ sequences — no reference genome, no alignment — and automatically classifies each library as double-stranded, single-stranded, or undetermined using a 4-channel joint BIC model.

---

## Why reference-free?

Reference-based damage estimation (mapDamage, metaDMG) requires a mapped BAM file. In many aDNA workflows — particularly metagenomics, sediment cores, or pre-screening pipelines — you want damage information before alignment, either to guide deduplication parameters or to triage samples. libdart-damage fills this gap.

---

## What it classifies

Ancient DNA library preparation matters for deduplication and damage masking:

| Library type | Damage pattern |
|---|---|
| Double-stranded (DS) | Symmetric C→T at 5' and G→A at 3', exponential decay |
| DS + end-repair artifact | As above plus a G→A spike at 3' position 0 |
| Single-stranded (SS) | Asymmetric: spike-only, C→T at both ends, or mixed |
| UNKNOWN | No signal above null model — zero-damage or heavily modern |

The BIC classifier uses four channels (ct5, ga3, ga0, ct3) and seven composite models to distinguish these patterns without user input.

---

## Quick start

```cpp
#include <dart/frame_selector_decl.hpp>

// One-shot
std::vector<std::string> reads = { "ACGTCTAGCT...", ... };
dart::SampleDamageProfile profile =
    dart::FrameSelector::compute_sample_profile(reads);

// Streaming
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

## Validated performance

Tested on 315 Mediterranean sediment aDNA libraries (two independent datasets):

| Dataset | Correct | UNKNOWN | Wrong | Accuracy (determined) |
|---|---|---|---|---|
| Dataset 1 (91 samples) | 88 | 3 | 0 | **100%** |
| Dataset 2 (224 samples) | 193 | 28 | 3 | **98.5%** |

UNKNOWN = no model beats the null (zero-damage libraries where library type cannot be inferred from sequence alone).

---

## Build

As a CMake dependency:

```cmake
find_package(dart-damage REQUIRED)
target_link_libraries(your_target PRIVATE dart-damage)
```

Or via FetchContent:

```cmake
include(FetchContent)
FetchContent_Declare(libdart_damage
    GIT_REPOSITORY https://github.com/genomewalker/libdart-damage.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(libdart_damage)
target_link_libraries(your_target PRIVATE dart-damage)
```

---

## Documentation

Full methods, API reference, and changelog: **https://genomewalker.github.io/libdart-damage**

- [Methods](https://genomewalker.github.io/libdart-damage/methods/) — damage model equations, BIC classifier design, rescue rules
- [API Reference](https://genomewalker.github.io/libdart-damage/api/) — FrameSelector methods, SampleDamageProfile fields
- [Changelog](https://genomewalker.github.io/libdart-damage/changelog/) — recent changes and milestones
