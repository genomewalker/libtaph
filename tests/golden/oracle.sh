#!/usr/bin/env bash
# Behavior-preservation oracle for the reference-free damage verdict.
# Re-profiles deterministic fixtures across ALL regimes (DETECTED ds/ss, LOW_ABUNDANCE,
# ABSTAIN, ROCS, modern negative) and compares VERDICT-CRITICAL JSON fields against
# per-fixture baselines captured from current behavior. Exit 0 = identical for every fixture.
#
# Verdict-critical fields (8):
#   pi_estimate.state          (exact string: DETECTED/LOW_ABUNDANCE/ABSTAIN/...)
#   pi_estimate.point          (float, abs tol 1e-6)
#   deamination.d_max_5prime   (float)
#   deamination.d_max_3prime   (float)
#   deamination.lambda_5prime  (float)
#   channel_b.d_max            (float)
#   channel_b.valid            (bool, exact)
#   channel_b.inverted         (bool, exact)
#
# Fixtures (fixture, library_type, baseline) — regime in comments:
#   fixture_det_ds  ds  baseline_det_ds.json   # DETECTED ds (FLB03mAds3)
#   fixture_lowab   ds  baseline_lowab.json    # LOW_ABUNDANCE (FLB10mAds3)
#   fixture_abstain ds  baseline_abstain.json  # ABSTAIN (FLB45mAds1)
#   fixture_ss      ss  baseline_ss.json       # ss (FLB03mAss4)
#   fixture_rocs    ds  baseline_rocs.json     # ABSTAIN now (LV3003046936)
#   fixture_modern  ds  baseline_modern.json   # ABSTAIN, clean negative (ERR5526852 head)
#
# Run with CAPTURE=1 to (re)capture every baseline from the current binary instead of diffing.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FQDUP="${FQDUP:-/path/to/fqdup/build/fqdup}"
TMPDIR_ORACLE="${TMPDIR_ORACLE:-/path/to/tmp}"
CAPTURE="${CAPTURE:-0}"

# fixture_tag  library_type
FIXTURES=(
  "det_ds  ds"
  "lowab   ds"
  "abstain ds"
  "ss      ss"
  "rocs    ds"
  "modern  ds"
)

# Rebuild fqdup from the working-tree libtaph (FETCHCONTENT_SOURCE_DIR_LIBTAPH) so the oracle ALWAYS
# reflects the current edits — incremental (~2s), and without it the oracle silently tests a stale binary.
( cd /path/to/fqdup && cmake --build build --parallel ) >/dev/null 2>&1 \
    || { echo "ORACLE BUILD FAIL"; exit 2; }

fail=0
for row in "${FIXTURES[@]}"; do
    read -r tag libtype <<<"$row"
    fixture="$HERE/fixture_${tag}.fq.gz"
    baseline="$HERE/baseline_${tag}.json"
    [ -f "$fixture" ] || { echo "MISSING FIXTURE $fixture"; fail=1; continue; }

    if [ "$CAPTURE" = "1" ]; then
        "$FQDUP" profile -i "$fixture" --library-type "$libtype" -p 8 --json "$baseline" >/dev/null 2>&1 \
            || { echo "CAPTURE FAIL $tag"; fail=1; continue; }
        echo "CAPTURED $tag -> $(basename "$baseline")"
        continue
    fi

    [ -f "$baseline" ] || { echo "MISSING BASELINE $baseline"; fail=1; continue; }
    candidate="$(mktemp "$TMPDIR_ORACLE/oracle_candidate.${tag}.XXXXXX.json")"
    "$FQDUP" profile -i "$fixture" --library-type "$libtype" -p 8 --json "$candidate" >/dev/null 2>&1 \
        || { echo "PROFILE FAIL $tag"; fail=1; rm -f "$candidate"; continue; }

    python3 - "$tag" "$baseline" "$candidate" <<'PY' || fail=1
import json, sys
tag, base_p, cand_p = sys.argv[1], sys.argv[2], sys.argv[3]
base = json.load(open(base_p))
cand = json.load(open(cand_p))
TOL = 1e-6
EXACT = ["pi_estimate.state", "channel_b.valid", "channel_b.inverted"]
FLOAT = ["pi_estimate.point", "deamination.d_max_5prime", "deamination.d_max_3prime",
         "deamination.lambda_5prime", "channel_b.d_max"]

def get(d, path):
    for k in path.split("."):
        d = d[k]
    return d

bad = 0
for p in EXACT:
    b, c = get(base, p), get(cand, p)
    if b != c:
        print(f"[{tag}] MISMATCH {p}: baseline={b!r} candidate={c!r}")
        bad = 1
for p in FLOAT:
    b, c = get(base, p), get(cand, p)
    if b is None or c is None:
        if b != c:
            print(f"[{tag}] MISMATCH {p}: baseline={b!r} candidate={c!r}")
            bad = 1
        continue
    if abs(float(b) - float(c)) > TOL:
        print(f"[{tag}] MISMATCH {p}: baseline={b} candidate={c} (|d|={abs(float(b)-float(c)):.3e} > {TOL})")
        bad = 1

state = get(cand, "pi_estimate.state")
if bad:
    print(f"[{tag}] FAIL (state={state})")
    sys.exit(1)
print(f"[{tag}] PASS (state={state})")
PY
    rm -f "$candidate"
done

if [ "$CAPTURE" = "1" ]; then
    [ "$fail" = 0 ] && echo "ORACLE CAPTURE COMPLETE" || echo "ORACLE CAPTURE INCOMPLETE"
    exit "$fail"
fi
[ "$fail" = 0 ] && echo "ORACLE PASS" || echo "ORACLE FAIL"
exit "$fail"
