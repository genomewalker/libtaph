---
type: libtaph-json-field
title: damage_channel_stats
tier: standard
estimand: 12-channel substitution terminal excess, interior fraction, decay lambda
stability: stable
emitted_by: profile_json.cpp
---

### `damage_channel_stats` block — per-channel terminal vs interior statistics

A flat summary of terminal excess, interior enrichment, and rate-decay parameters
for all 12 substitution channels (`CT CA CG GA GT GC TC TA TG AG AT AC`) at both
the 5′ and 3′ ends. This provides a channel-level view complementary to the
context-specific blocks above.

The block has two sub-objects, `5prime` and `3prime`, each containing:

| Field | Description |
|-------|-------------|
| `interior_start_pos` | Read position at which interior counting begins (default 3) |
| `ct_interior_fraction` | Interior C→T fraction — the composition baseline used for excess computation |
| `ct_decay_lambda` | Exponential decay constant fitted to the CT terminal profile; `null` when fit is not reliable |
| `channels` | Object keyed by two-letter substitution code; see sub-fields below |

Per-channel sub-fields:

| Sub-field | Description |
|-----------|-------------|
| `interior_fraction` | Fraction of reads showing this substitution at interior positions |
| `interior_log_ratio` | log(`interior_fraction` / 0.5) — signed distance from the uniform expectation; negative = interior-depleted (terminal-enriched) |
| `terminal_excess` | Terminal substitution rate minus interior rate; positive = enriched at the terminus |
| `decay_lambda` | Exponential decay constant for this channel's terminal profile; `null` when fit is unreliable or excess is near zero |
| `lambda_ratio_vs_ct` | `decay_lambda` of this channel divided by `ct_decay_lambda` — normalises for read-length differences; 1.0 = same spatial decay as C→T |
| `coupled_balance` | For paired channels (GT↔TG, CT↔TC): ratio of terminal excess on the 3′ complement to expected from the 5′ excess under Chargaff symmetry; `null` when the expected excess is near zero. Values near 1 indicate strand-balanced damage; values ≫ 1 indicate strand-asymmetric terminal enrichment |

**Interpreting `coupled_balance` for oxidative damage.** For genuine in-situ 8-oxoG
oxidation of double-stranded DNA, G→T damage is distributed across the molecule and
is not strongly terminal-enriched (`terminal_excess` near zero or small negative for GT).
A `coupled_balance` value ≫ 1 for the GT channel at 5′ would indicate strand-asymmetric
terminal G→T enrichment consistent with an OxoG library-preparation artefact (Costello
2013). Values near 1 or negative `terminal_excess` support genuine in-situ oxidation.

---