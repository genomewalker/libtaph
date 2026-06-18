---
type: libtaph-json-field
title: depurination_deconvolution
tier: standard
estimand: AP-site C→T pathway signal deconvolved from bulk deamination
stability: stable
emitted_by: profile_json.cpp
---

### `depurination_deconvolution` block — cut-site base preference via trinucleotide deconvolution (Channel E extended)

The hydrolytic cleavage of the glycosidic bond at a purine leaves an apurinic (AP) site that
undergoes β-elimination to produce a strand break. The base lost by this process — adenine or
guanine — sits at genomic position −1 relative to the read start: it is not the first sequenced
base, but the base immediately upstream of it on the cleaved strand. Because that base is outside
the read, its identity cannot be measured directly. However, the read-start trinucleotide
composition is a mixture of all genomic trinucleotides (x, b₀, b₁) weighted by the cut-site
preference w(x) for the hidden −1 base x. Interior trinucleotide frequencies, which are
accumulated from positions well away from the terminus, provide an unbiased estimate of the
genomic background P_gen3(x, b₀, b₁). Inverting this relationship via non-negative least squares
(NNLS) recovers w(x) — the relative preference for each nucleotide at the hidden cut-site —
without requiring a reference genome.

The NNLS system is strongly overdetermined: 16 observed read-start dinucleotide counts (b₀, b₁)
constrain 4 unknowns w(A), w(G), w(C), w(T), subject to non-negativity and Σw = 1. The resulting
`depurination_index` is the purine-to-pyrimidine weight ratio (w(A) + w(G)) / (w(C) + w(T)).
A value substantially above 1 indicates that fragment ends are enriched for sites where a purine
occupied the −1 position, which is the expected signature of depurination-driven scission.
A value near or below 1 is consistent with sequence-neutral cutting (e.g., nuclease activity or
mechanical shearing) rather than base-loss chemistry.

| Field | Type | Description |
|-------|------|-------------|
| `w_A` | float | NNLS mixture weight for adenine at the hidden −1 cut-site; range [0, 1], constrained non-negative |
| `w_C` | float | NNLS mixture weight for cytosine at the hidden −1 cut-site; range [0, 1], constrained non-negative |
| `w_G` | float | NNLS mixture weight for guanine at the hidden −1 cut-site; range [0, 1], constrained non-negative |
| `w_T` | float | NNLS mixture weight for thymine at the hidden −1 cut-site; range [0, 1], constrained non-negative |
| `depurination_index` | float | (w_A + w_G) / (w_C + w_T); values >1 indicate purine-preferential fragmentation; 0.0 when the purine weight is zero |

The four weights sum to 1.0 by construction. A `depurination_index` of 0.0 arises when the NNLS
solution places no weight on either purine, which can reflect genuinely pyrimidine-biased cuts or
collapse of the deconvolution due to collinear background columns — cross-check with the parent
`depurination` block's `valid` and `detected` fields before concluding absence of signal. High
index values (>2) with roughly equal `w_A` and `w_G` are broadly consistent with hydrolytic AP
scission; asymmetric solutions (e.g., dominant `w_G`) may reflect guanine-specific oxidative
fragmentation rather than classical acid-catalysed depurination and should be interpreted
alongside the `oxog_estimate` block. Because the deconvolution uses read-start composition, any
residual C→T deamination inflation at position 0 can inflate the apparent pyrimidine weight and
suppress the index; the index is most reliable when terminal deamination rates are low or have
been accounted for.


---