#!/usr/bin/env bash
# Full-schema regression diff, companion to oracle.sh. oracle.sh deliberately checks
# only 8 verdict-critical fields (see its header) -- a narrow, high-signal set for
# "did the DETECTED/LOW_ABUNDANCE/ABSTAIN call change". It does NOT cover oxidation,
# the joint model, hexamer diagnostics, or (as of the Channel B consolidation into
# context_primitives) most of the schema at all -- a bug there passes the oracle
# silently. This script recursively diffs EVERY leaf value in the JSON output against
# the same per-fixture baselines oracle.sh uses, so ongoing schema-consolidation work
# (renames, moves into context_primitives, new deprecation aliases) has a standing,
# comprehensive safety net instead of a one-off diff script per change.
#
# This is NOT a replacement for oracle.sh (different purpose: schema completeness vs.
# verdict correctness) -- run both. Uses the SAME fixtures/baselines as oracle.sh, so
# CAPTURE=1 in oracle.sh keeps both scripts pointed at a consistent, current baseline.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FQDUP_REPO="${FQDUP_REPO:-}"
FQDUP_GIVEN="${FQDUP:+1}"
FQDUP="${FQDUP:-${FQDUP_REPO:+$FQDUP_REPO/build/fqdup}}"
TMPDIR_ORACLE="${TMPDIR_ORACLE:-${TMPDIR:-/tmp}}"

FIXTURES=(
  "det_ds  ds"
  "lowab   ds"
  "abstain ds"
  "ss      ss"
  "rocs    ds"
  "modern  ds"
)

if [ -z "$FQDUP" ]; then echo "SKIP: set FQDUP or FQDUP_REPO"; exit 0; fi
if [ -z "$FQDUP_GIVEN" ]; then
    ( cd "$FQDUP_REPO" && cmake --build build --parallel ) >/dev/null 2>&1 \
        || { echo "FULL_DIFF BUILD FAIL"; exit 2; }
fi

fail=0
for row in "${FIXTURES[@]}"; do
    read -r tag libtype <<<"$row"
    fixture="$HERE/fixture_${tag}.fq.gz"
    baseline="$HERE/baseline_${tag}.json"
    [ -f "$fixture" ] || { echo "MISSING FIXTURE $fixture"; fail=1; continue; }
    [ -f "$baseline" ] || { echo "MISSING BASELINE $baseline"; fail=1; continue; }

    candidate="$(mktemp "$TMPDIR_ORACLE/full_diff_candidate.${tag}.XXXXXX.json")"
    "$FQDUP" profile -i "$fixture" --library-type "$libtype" -p 8 --json "$candidate" >/dev/null 2>&1 \
        || { echo "PROFILE FAIL $tag"; fail=1; rm -f "$candidate"; continue; }

    python3 - "$tag" "$baseline" "$candidate" <<'PY' || fail=1
import json, sys

tag, base_p, cand_p = sys.argv[1], sys.argv[2], sys.argv[3]
base = json.load(open(base_p))
cand = json.load(open(cand_p))
TOL = 1e-6
mismatches = []

def diff(path, b, c):
    if isinstance(b, dict) and isinstance(c, dict):
        for k in sorted(set(b) | set(c)):
            if k not in b:
                mismatches.append(f"{path}.{k}: MISSING in baseline, present in candidate={c[k]!r}")
            elif k not in c:
                mismatches.append(f"{path}.{k}: present in baseline={b[k]!r}, MISSING in candidate")
            else:
                diff(f"{path}.{k}", b[k], c[k])
    elif isinstance(b, list) and isinstance(c, list):
        if len(b) != len(c):
            mismatches.append(f"{path}: length baseline={len(b)} candidate={len(c)}")
            return
        for i, (bi, ci) in enumerate(zip(b, c)):
            diff(f"{path}[{i}]", bi, ci)
    elif isinstance(b, (int, float)) and isinstance(c, (int, float)) and not isinstance(b, bool) and not isinstance(c, bool):
        if abs(float(b) - float(c)) > TOL:
            mismatches.append(f"{path}: baseline={b} candidate={c} (|d|={abs(float(b)-float(c)):.3e})")
    else:
        if b != c:
            mismatches.append(f"{path}: baseline={b!r} candidate={c!r}")

diff("$", base, cand)

if mismatches:
    print(f"[{tag}] FULL_DIFF FAIL ({len(mismatches)} mismatches)")
    for m in mismatches[:50]:
        print(f"  {m}")
    if len(mismatches) > 50:
        print(f"  ... and {len(mismatches)-50} more")
    sys.exit(1)
print(f"[{tag}] FULL_DIFF PASS (0 mismatches)")
PY
    rm -f "$candidate"
done

[ "$fail" = 0 ] && echo "FULL_DIFF PASS" || echo "FULL_DIFF FAIL"
exit "$fail"
