# libdart-damage

Reference-free ancient DNA damage estimation and library-type classification from raw FASTQ reads.

libdart-damage scans raw FASTQ reads for six independent damage and fragmentation channels, estimates terminal deamination, cross-validates C→T damage with composition-robust stop-codon signals, and classifies each library as double-stranded (DS), single-stranded (SS), or UNKNOWN — all without a reference genome or read alignment.

---

## Why reference-free?

Reference-based damage estimation (mapDamage, metaDMG) requires a mapped BAM file. In aDNA workflows where damage information is needed before alignment (metagenomics, sediment cores, pre-screening), libdart-damage operates on raw sequences to guide deduplication parameters and sample triage.

---

## What it measures

Six biochemical damage channels plus three context-sensitive metrics, together
covering all reference-free FASTQ-detectable aDNA damage types:

### Channel A: Cytosine deamination (primary damage signal)

Cytosine in single-stranded DNA overhangs spontaneously deaminates to uracil,
which the polymerase reads as thymine. The rate is highest at the fragment
terminus and decays exponentially inward. This produces:

- **ct5**: C→T excess at 5' terminal positions (positions 0–14)
- **ga3**: G→A excess at 3' terminal positions, the complementary strand's
  deaminated cytosines appear as G→A when read 5'→3' from that end
- **ga0**: isolated G→A spike at 3' position 0. Common in single-stranded
  complement-orientation libraries, but can also arise as a double-stranded
  end-repair artifact; the classifier distinguishes these cases using the
  full four-channel pattern.
- **ct3**: C→T excess at 3' terminal positions (single-stranded original-strand
  libraries, where the same strand is damaged at both ends)

Channel A drives position masking in fqdup. The decay model is
`d_max × exp(−λ × pos) + bg`; positions where the excess above background
exceeds a threshold are masked before hashing.

### Channel B: Stop codon conversion at 5' (composition-independent C→T validation)

The same C→T deamination converts specific sense codons to stop codons:
CAA/CAG/CGA → TAA/TAG/TGA. Because this test uses triplet context rather
than raw base frequency, it is insensitive to GC-composition bias that can
confound Channel A. If Channel A fires but Channel B contradicts it, the
signal is flagged as `damage_artifact` (likely a composition artifact, not
genuine ancient damage).

### Channel B₃′: Stop codon conversion at 3' (SS G→A validation)

The complement of Channel B at the 3' end: TGG → TAG or TGA via G→A
deamination. Validates the 3' damage signal, particularly in SS libraries
where smooth ga3 decay may be absent and the G→A signal is limited to
position 0. `d_max` from Channel B₃′ can serve as an alternative structural
estimate when Channel A is unreliable at the 3' end.

### Channel C: 8-oxoG oxidative damage (stop codon uniformity)

8-Oxoguanine forms when guanine is oxidised. The polymerase misreads it as
adenine, producing G→T transversions. Unlike deamination, oxidation
accumulates throughout the molecule rather than preferentially at the ends.
Channel C tracks G→T stop codon conversions (GAG/GAA/GGA → TAG/TAA/TGA)
and tests whether they are uniformly distributed (genuine oxidation) or
terminal-enriched (likely co-occurring with deamination rather than
independent 8-oxoG).

### Channel D: 8-oxoG direct transversion rate

Measures G→T and its complement C→A as raw transversion frequencies at
terminal vs. interior positions. High terminal/interior ratio with correlated
G→T and C→A signals on opposite strands confirms genuine 8-oxoG oxidation.
Flags `ox_damage_detected` and `ox_is_artifact` independently of Channel C.

### Channel E: Depurination / AP-site fragmentation

Purines (A and G) are lost by hydrolysis under acidic or warm burial
conditions, leaving apurinic (AP) sites that become strand-break points.
Fragmentation at AP sites enriches purines at the newly exposed 5' ends.
Channel E measures this purine enrichment at 5' terminal positions versus
the interior baseline. High enrichment confirms that fragmentation occurred
at AP sites, independent evidence of genuine ancient origin even when
deamination is low.

---

Channels B–E cross-validate Channel A: if Channel A detects deamination but
Channel B contradicts it, the signal is flagged as a composition artifact
rather than genuine ancient damage. Channels C–E are reported for sample
characterisation but do not directly affect position masking.

### CpG-like context split

The same cytosine deamination that drives Channel A acts faster at CpG
dinucleotides when the cytosine was methylated (5mC → T). libdart-damage fits
the 5' C→T amplitude separately for positions where the next base is G
(CpG-like context) and positions where it is not (non-CpG). A reference-free
interior baseline (middle third of each read) is used to estimate the
background T fraction in each context. The ratio `cpg_ratio = dmax_cpg /
dmax_noncpg` quantifies the methylation-enhanced component. Values near 1
indicate uniform deamination independent of methylation status; values
substantially above 1 suggest enriched CpG methylation in the source organism.
Reported in the JSON block `deamination.cpg_like`.

### Interior C→T clustering

Beyond terminal overhangs, deamination can cluster in read interiors when
closely spaced cytosines in single-stranded micro-domains are co-deaminated in
the same hydrolytic event. libdart-damage measures excess co-occurrence of T at
adjacent non-CpG `{C,T}` sites within the read interior (middle third) at
pair separations d = 1–10 bp, using the within-read T fraction as the null.
An AG-track control (analogous G→A co-occurrence) corrects for strand-composition
biases. The summary statistic `short_asym_log2oe` is the background-corrected
log₂ observed/expected ratio; `short_z` is the normalised contrast. Significant
positive values indicate genuine clustered interior deamination not captured by
the terminal exponential model. Reported in the JSON block
`interior_ct_cluster`.

### 8-oxoG 16-context panel

The overall `s_oxog` statistic (G→T strand asymmetry) is split across all 16
flanking-dinucleotide bins (N**G**N contexts). The resulting 16-element vector
`s_oxog_16ctx` characterises the trinucleotide specificity of oxidation: genuine
ancient 8-oxoG shows broad context enrichment, while modern oxidation from
sample preparation tends to be context-biased. Reported in the JSON block
`complement_asymmetry.s_oxog_16ctx`.

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

| Library type | Dominant channels | Typical damage plot appearance |
|---|---|---|
| DS | ct5 + ga3 | Symmetric 5' C→T and 3' G→A, both smooth decay |
| DS + end-repair artifact | ct5 + ga3 + ga0 | DS pattern plus isolated 3' pos-0 G→A spike |
| SS complement-orientation | ga0-dominated | Strong 3' pos-0 G→A spike; little or no smooth ga3; 5' flat |
| SS original-orientation | ct5 + ct3 | C→T at both 5' and 3' ends; no G→A |
| SS mixed orientations | ct5 + ga0 | 5' C→T decay plus 3' pos-0 G→A spike; weak residual ga3 may occur |
| UNKNOWN | none | No channel clearly beats the null model |

UNKNOWN is the correct call for zero-damage or near-zero-damage libraries where
the library type cannot be inferred from sequence alone, it is not an error.

---

## fqdup integration

[fqdup](https://github.com/genomewalker/fqdup) uses libdart-damage for
reference-free damage estimation before deduplication. The standalone
`fqdup damage` subcommand runs the full profiler and prints a human-readable
report:

```
$ fqdup damage -i merged.fq.gz

=== fqdup damage ===
Input:   merged.fq.gz
Threads: 16
Reads:   5,582,073 scanned
Length:  min=30  mean=91.2  max=150

Library: DS (auto-detected)
  BIC  bias=0.0  DS=125432.1  SS=42.0  SS_full=89.3
  fit  CT5_amp=0.1928  ΔBIC=125432  GA3_amp=0.0403  ΔBIC=1249  GA0_amp=0.0201  ΔBIC=22  CT3_amp=0.0012  ΔBIC=0.1
  5' terminal shift: +0.0193  (z=18.4)
  3' terminal shift: +0.0040  (z=3.9)

5'-end   d_max=0.1928  lambda=0.246  bg=0.4872
3'-end   d_max=0.0403  lambda=0.069  bg=0.5091
combined d_max=0.1928  (source=5prime) [validated]

Mask threshold: 0.05 → 1 position masked (pos 0)

pos  5'_CT   5'_GA   3'_GA
---  ------  ------  ------
  0  0.2194  0.4895  0.0521  *
  1  0.1571  0.4888  0.0458
  2  0.1048  0.4891  0.0432
  ...
```

The report shows library type with BIC evidence, decay parameters per end,
per-position frequencies, and which positions exceed the mask threshold.
Run it before `fqdup extend` or `fqdup derep --damage-auto` to verify
auto-detection or supply parameters manually.

---

## Quick start

```cpp
#include <dart/frame_selector_decl.hpp>

// One-shot
std::vector<std::string> reads = { "ACGTCTAGCT...", ... };
dart::SampleDamageProfile profile =
    dart::FrameSelector::compute_sample_profile(reads);

// Streaming
dart::SampleDamageProfile profile{};
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
