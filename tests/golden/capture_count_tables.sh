#!/usr/bin/env bash
# Capture Layer-0 count-table golden fixtures from the validated panel at -p8 (count tables are
# thread-order-independent, verified). Fixtures are named by pathology.
# Inputs live in the ellesmere project (external to the repo); see manifest.tsv.
set -uo pipefail
FQ="${FQDUP:-}"
ROOT="${ELLESMERE_ROOT:-}"
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/count_tables"
mkdir -p "$OUT"

profile() {  # pathology_name  input_relpath
  local name="$1" in="$ROOT/$2"
  echo "[$(date +%H:%M:%S)] capturing $name"
  "$FQ" profile -i "$in" -p 8 --library-type auto --length-bins auto \
      --count-table-json "$OUT/$name.jsonl" >/dev/null 2>&1
  echo "[$(date +%H:%M:%S)]   $name: $(wc -l < "$OUT/$name.jsonl" 2>/dev/null || echo MISSING) rows"
}

profile cap_masked_f_reversal_ds  zihao_kap_capture/fastq/merged.fq.gz
profile ds_artifact_cap_masked    derep_split/flb/FLB57mdBds3LV7005225899PI2GK/deduped.fq.gz
profile ss_ancient                retrim/FLB01mAss4_merged.fastq.gz
profile marine_lowcov             derep_split/rocs/LV7001898061/deduped.fq.gz
echo "[$(date +%H:%M:%S)] === CAPTURE DONE ==="
