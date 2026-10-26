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
# fails (exit 1) on any count-table drift, missing fixture, or non-zero fqdup exit.
#
# NOT HERMETIC -- read before trusting a failure:
#  * The default FQ is build-icelake, an AVX-512 (-march=icelake) build, while this gate runs on
#    whatever node the scheduler hands out. On a CPU without that ISA it SIGILLs, and only on the
#    fixtures whose code paths reach the vectorised code -- so it fails per-fixture and looks
#    exactly like drift. Seen 2026-07-15: ss_ancient SIGILL on dandycmpn01fl (96c, no avx512)
#    while marine_lowcov passed on the same node, and all fixtures passed on dandycmpn17fl (192c).
#    Pin the job (--nodelist/--constraint) or set FQDUP to a portable build.
#  * FQ is a PRE-BUILT binary, NOT built from the tree under test. A stale FQ silently gates the
#    wrong code (the default was 12 days old when this note was written). Rebuild it, or set FQDUP.
set -uo pipefail
FQ="${FQDUP:-}"
ROOT="${ELLESMERE_ROOT:-}"
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
  # Check the exit status. Ignoring it made ANY failure -- crash, OOM, bad flag -- come out as
  # "count-table drift", because the diff then compared against a file the dead process never
  # wrote. A binary that died verified nothing; say so instead of blaming the counts.
  "$FQ" profile -i "$ROOT/$in" -p 8 --library-type auto --length-bins auto \
      --count-table-json "$TMP/$name.jsonl" >/dev/null 2>"$TMP/$name.err"
  rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "FAIL $name: fqdup EXITED rc=$rc (NOT count-table drift -- it produced no output)"
    if [[ $rc -eq 132 ]]; then
      # 128+4 = SIGILL. The default FQ is a -march=icelake (AVX-512) build, and this gate runs
      # on whatever node the scheduler picks; a CPU without that ISA dies only on the code paths
      # that actually reach it (so some fixtures still pass, which is what makes it look like
      # per-fixture drift). Pin the job to a capable node or point FQDUP at a portable build.
      echo "     SIGILL: '$FQ' targets an ISA this CPU lacks ($(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //'))."
      echo "     avx512f present: $(grep -qo avx512f /proc/cpuinfo && echo yes || echo NO)"
      echo "     Fix: run on a node with the ISA (sbatch --nodelist/--constraint), or FQDUP=/path/to/portable/fqdup"
    fi
    sed 's/^/     /' "$TMP/$name.err" | head -6
    fail=1
    continue
  fi
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
