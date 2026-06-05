#!/usr/bin/env bash
# Layer-0 count-table golden gate.
#
# Re-profile each fixtured sample at -p8 and diff its Layer-0 stop-channel count tables against the
# checked-in fixture. The count tables are thread-order-independent (integer-count sums, verified) and
# emitted at setprecision(17), so a byte-exact diff IS the pre-rounded numeric check: it catches
# cap-masked z regressions (e.g. F's pre_clamp_z drifting while the emitted clamped z stays at -12),
# denominator-membership swaps, and shadow-policy changes that the emitted-JSON gate cannot see.
#
# Self-skips (exit 0) when fqdup or the external inputs are absent so CI without the data is green;
# fails (exit 1) on any count-table drift or missing fixture.
set -uo pipefail
FQ="${FQDUP:-/path/to/fqdup/build-icelake/fqdup}"
ROOT="${ELLESMERE_ROOT:-/path/to/ellesmere}"
HERE="$(cd "$(dirname "$0")" && pwd)"
MAN="$HERE/manifest.tsv"
FIX="$HERE/count_tables"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

if [[ ! -x "$FQ" ]]; then echo "SKIP: fqdup binary not found ($FQ)"; exit 0; fi

fail=0; n=0; skipped=0
while IFS=$'\t' read -r name lib in note; do
  [[ -z "${name:-}" || "$name" == \#* ]] && continue
  fix="$FIX/$name.jsonl"
  if [[ ! -f "$fix" ]]; then echo "FAIL $name: fixture missing ($fix)"; fail=1; continue; fi
  if [[ ! -f "$ROOT/$in" ]]; then echo "SKIP $name: input absent ($in)"; skipped=$((skipped+1)); continue; fi
  n=$((n+1))
  "$FQ" profile -i "$ROOT/$in" -p 8 --library-type auto --length-bins auto \
      --count-table-json "$TMP/$name.jsonl" >/dev/null 2>&1
  if diff -q "$fix" "$TMP/$name.jsonl" >/dev/null 2>&1; then
    echo "PASS $name"
  else
    echo "FAIL $name: count-table drift"
    diff "$fix" "$TMP/$name.jsonl" 2>&1 | head -12
    fail=1
  fi
done < "$MAN"
echo "--- count-table golden gate: $n checked, $skipped skipped, result=$([[ $fail -eq 0 ]] && echo PASS || echo FAIL) ---"
exit $fail
