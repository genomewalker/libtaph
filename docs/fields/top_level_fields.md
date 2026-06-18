---
type: libtaph-json-field
title: top_level_fields
tier: summary
estimand: scalar top-level fields: d_max, library_type, n_reads, schema_version
stability: stable
emitted_by: profile_json.cpp
---

### Top-level fields

| Field | Type | Description |
|-------|------|-------------|
| `input` | string | Input FASTQ path |
| `n_reads` | int | Total reads scanned |
| `library_type` | string | `"double-stranded"`, `"single-stranded"`, or `"unknown"` |
| `library_type_auto` | bool | `true` if type was auto-detected (not forced by `--library-type`) |
| `library_type_rescued` | bool | `true` if a rescue rule overrode the primary BIC winner |
| `damage_status` | string | `"present"`, `"weak"`, or `"absent"` |

---