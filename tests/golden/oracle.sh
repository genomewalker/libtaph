#!/usr/bin/env bash
# Behavior-preservation oracle for the reference-free damage verdict.
# Re-profiles the deterministic fixture and compares VERDICT-CRITICAL JSON fields
# against baseline.json (captured from the pre-refactor HEAD). Exit 0 = identical.
#
# Verdict-critical fields:
#   pi_estimate.state          (exact string: DETECTED/LOW_ABUNDANCE/ABSTAIN/...)
#   pi_estimate.point          (float, abs tol 1e-6)
#   deamination.d_max_5prime   (float)
#   deamination.d_max_3prime   (float)
#   deamination.lambda_5prime  (float)
#   channel_b.d_max            (float)
#   channel_b.valid            (bool, exact)
#   channel_b.inverted         (bool, exact)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FQDUP="${FQDUP:-/path/to/fqdup/build/fqdup}"
FIXTURE="$HERE/fixture_ds.fq.gz"
BASELINE="$HERE/baseline.json"
CANDIDATE="$(mktemp /path/to/tmp/oracle_candidate.XXXXXX.json)"
trap 'rm -f "$CANDIDATE"' EXIT

# Rebuild fqdup from the working-tree libtaph (FETCHCONTENT_SOURCE_DIR_LIBTAPH) so the oracle ALWAYS
# reflects the current edits — incremental (~2s), and without it the oracle silently tests a stale binary.
( cd /path/to/fqdup && cmake --build build --parallel ) >/dev/null 2>&1 \
    || { echo "ORACLE BUILD FAIL"; exit 2; }

"$FQDUP" profile -i "$FIXTURE" --library-type ds -p 8 --json "$CANDIDATE" >/dev/null 2>&1

python3 - "$BASELINE" "$CANDIDATE" <<'PY'
import json, sys
base = json.load(open(sys.argv[1]))
cand = json.load(open(sys.argv[2]))
TOL = 1e-6
EXACT = ["pi_estimate.state", "channel_b.valid", "channel_b.inverted"]
FLOAT = ["pi_estimate.point", "deamination.d_max_5prime", "deamination.d_max_3prime",
         "deamination.lambda_5prime", "channel_b.d_max"]

def get(d, path):
    for k in path.split("."):
        d = d[k]
    return d

fail = 0
for p in EXACT:
    b, c = get(base, p), get(cand, p)
    if b != c:
        print(f"MISMATCH {p}: baseline={b!r} candidate={c!r}")
        fail = 1
for p in FLOAT:
    b, c = get(base, p), get(cand, p)
    if b is None or c is None:
        if b != c:
            print(f"MISMATCH {p}: baseline={b!r} candidate={c!r}")
            fail = 1
        continue
    if abs(float(b) - float(c)) > TOL:
        print(f"MISMATCH {p}: baseline={b} candidate={c} (|d|={abs(float(b)-float(c)):.3e} > {TOL})")
        fail = 1

if fail:
    print("ORACLE FAIL")
    sys.exit(1)
print("ORACLE PASS")
PY
