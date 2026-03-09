# libdart-damage

Reference-free ancient DNA damage estimation and library-type classification from raw FASTQ reads.

libdart-damage scans raw FASTQ reads for five independent biochemical damage signals, quantifies each, and classifies each library as double-stranded (DS), single-stranded (SS), or undetermined — all without a reference genome or alignment.

---

## Why reference-free?

Reference-based damage estimation (mapDamage, metaDMG) requires a mapped BAM file. In aDNA workflows where damage information is needed before alignment (metagenomics, sediment cores, pre-screening), libdart-damage operates on raw sequences to guide deduplication parameters and sample triage.

---

## What it measures

Five biochemical damage channels, each tracking a distinct degradation process:

| Channel | Chemistry | Signal | Read position |
|---------|-----------|--------|---------------|
| A | Cytosine deamination | C→T at 5' (ct5), G→A at 3' (ga3/ga0), C→T at 3' (ct3) | Terminal, exponential decay |
| B | Stop codon conversion via C→T | CAA/CAG/CGA → TAA/TAG/TGA excess at 5' terminus | Terminal, exponential decay |
| C | 8-oxoG oxidation (stop codon) | GAG/GAA/GGA → TAG/TAA/TGA, uniform across read | Uniform |
| D | 8-oxoG oxidation (direct) | G→T and C→A transversions, uniformity test | Uniform |
| E | Depurination / AP-site fragmentation | Purine enrichment at 5' strand-break junctions | 5' terminal |

Channels B–E cross-validate Channel A: if Channel A detects deamination but Channel B
contradicts it, the signal is flagged as a composition artifact rather than genuine ancient
damage.

## What it classifies

The library preparation protocol determines which damage channels are active,
because only the sequenced strand contributes to the observed substitution
pattern.

**Double-stranded (DS) libraries** (NEBNext, TruSeq) ligate adapters to both
strands. Reads from the original strand show C→T at the 5' end; reads from
the complementary strand map the same deaminated cytosines onto the 3' end as
G→A. Both ends show smooth exponential decay with similar amplitude, making
the signal symmetric.

**Single-stranded (SS) libraries** (Gansauge & Meyer 2013, SRSLY) circularize
and sequence only one strand per molecule. Which damage channels are active
depends on which strand is sequenced:

- **Complement-orientation**: the ligation junction leaves a position-0 G→A
  spike at the 3' end without a smooth ga3 decay. No 5' C→T signal because
  the original strand is not sequenced.
- **Original-orientation**: the original damaged strand is sequenced, so C→T
  appears at both the 5' end (from 5'-overhang deamination) and the 3' end
  (ct3; from 3'-overhang deamination on the same strand). No G→A signal.
- **Mixed**: libraries where both orientations are present produce a 5' C→T
  decay (from original-orientation reads) plus a 3' position-0 spike (from
  complement-orientation reads), but no smooth ga3 decay.

The BIC classifier evaluates seven models jointly against all four sub-channels
and selects the best-fitting description:

| Library type | Active channels | Damage plot appearance |
|---|---|---|
| DS | ct5 + ga3 | Symmetric C→T (5') and G→A (3'), both smooth decay |
| DS + end-repair artifact | ct5 + ga3 + ga0 | As above + isolated spike at 3' position 0 |
| SS complement-orientation | ga0 | Spike at 3' position 0 only; 5' flat |
| SS original-orientation | ct5 + ct3 | C→T at both 5' and 3' ends; no G→A |
| SS mixed orientations | ct5 + ga0 | C→T decay at 5' + spike at 3' pos 0; no smooth ga3 |
| UNKNOWN | — | No channel above null; standard exact-match deduplication |

UNKNOWN is the correct call for zero-damage or near-zero-damage libraries where
the library type cannot be inferred from sequence alone — it is not an error.

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

UNKNOWN: no model beats the null (zero-damage libraries where library type cannot be inferred from sequence alone).

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

- [Methods](https://genomewalker.github.io/libdart-damage/methods/): damage model equations, BIC classifier design, rescue rules
- [API Reference](https://genomewalker.github.io/libdart-damage/api/): FrameSelector methods, SampleDamageProfile fields
- [Changelog](https://genomewalker.github.io/libdart-damage/changelog/): recent changes and milestones
