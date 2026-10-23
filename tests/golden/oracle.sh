#!/usr/bin/env bash
# Behavior-preservation oracle for the reference-free damage verdict.
# Re-profiles deterministic fixtures across ALL regimes (DETECTED ds/ss, LOW_ABUNDANCE,
# ABSTAIN, ROCS, modern negative) and compares VERDICT-CRITICAL JSON fields against
# per-fixture baselines captured from current behavior. Exit 0 = identical for every fixture.
#
# Verdict-critical fields (7). Each is read at its v10 path, then its v11 path, so the pre-v11
# baselines keep their full assertion power and a future schema move cannot silently blind this
# oracle again (the v11 six-section refactor did exactly that: it left every lookup throwing
# KeyError, which reads as "broken test" rather than "the caller changed its mind", and a negative
# control flipped to "damaged" behind it):
#   verdict.deamination.call   (exact: present/absent)  <- THE call; catches a negative control flipping
#   deamination.d_max_5prime_raw_ungated  (float, abs tol 1e-6)
#   deamination.d_max_3prime_raw_ungated  (float)
#   deamination.lambda_5prime             (float)
#   channel_b.d_max                       (float)
#   channel_b.valid                       (bool, exact)
#   channel_b.inverted                    (bool, exact)
# pi_estimate.{state,point} are GONE, not moved: v11 replaced the estimator (shape-fit pi demoted to
# a diagnostic, co-occurrence promoted), so the old values have no successor to assert against.
#
# authenticated_by IS asserted, against the explicit EXPECTED_ROUTE table below -- not against the
# baselines, which predate the current routes. Asserting only `call` lets a fixture pass for the
# WRONG REASON: the verdict combiner is an OR over authentication routes, so a route that breaks and
# a route that over-fires cancel in `call` and the oracle sees nothing. That is not hypothetical --
# it is exactly how the modern negative control went "present" (its d5 fit authenticated on a
# context-marginal control) while every call-level assertion stayed green. The route is the claim;
# pin it. A deliberate route change must be re-pinned here in the same commit that causes it.
#
# Fixtures (fixture, library_type, expected_route, baseline) — regime in comments:
#   fixture_det_ds  ds  pi_estimate        baseline_det_ds.json   # DETECTED ds (FLB03mAds3)
#   fixture_lowab   ds  did5_length_decay  baseline_lowab.json    # LOW_ABUNDANCE (FLB10mAds3)
#   fixture_abstain ds  none               baseline_abstain.json  # ABSTAIN (FLB45mAds1)
#   fixture_ss      ss  d5_decay_fit       baseline_ss.json       # ss (FLB03mAss4)
#   fixture_rocs    ds  channel_b          baseline_rocs.json     # ABSTAIN now (LV3003046936)
#   fixture_modern  ds  none               baseline_modern.json   # ABSTAIN, clean negative (ERR5526852 head)
#
# lowab authenticates via did5_length_decay because pi co-occurrence is SECOND-order in damage
# (x ~ pi^2 * d5 * d3) and abstains by construction at its ~4% damage, while the DiD is first-order.
# rocs rests on channel_b, which finalize_dmax.cpp:361-363 states "cannot authenticate" -- it is the
# last uncontrolled route in the OR. That is a filed, open defect; this table pins the status quo so
# that fixing it shows up here as a loud, deliberate change rather than a silent drift.
#
# Run with CAPTURE=1 to (re)capture every baseline from the current binary instead of diffing.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FQDUP_GIVEN="${FQDUP:+1}"
FQDUP="${FQDUP:-/path/to/fqdup/build/fqdup}"
TMPDIR_ORACLE="${TMPDIR_ORACLE:-/path/to/tmp}"
CAPTURE="${CAPTURE:-0}"

# fixture_tag  library_type  expected authenticated_by route
FIXTURES=(
  "det_ds  ds  pi_estimate"
  "lowab   ds  did5_length_decay"
  "abstain ds  none"
  "ss      ss  d5_decay_fit"
  "rocs    ds  channel_b"
  "modern  ds  none"
)

# With the default fqdup, rebuild it from the working-tree libtaph (FETCHCONTENT_SOURCE_DIR_LIBTAPH)
# so the oracle reflects the current edits (incremental, ~2s); otherwise it silently tests a stale
# binary. An explicit FQDUP is used as given: it is never rebuilt, so the oracle cannot touch another
# build directory.
if [ -z "$FQDUP_GIVEN" ]; then
    ( cd /path/to/fqdup && cmake --build build --parallel ) >/dev/null 2>&1 \
        || { echo "ORACLE BUILD FAIL"; exit 2; }
fi

fail=0
for row in "${FIXTURES[@]}"; do
    read -r tag libtype want_route <<<"$row"
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

    python3 - "$tag" "$baseline" "$candidate" "$want_route" <<'PY' || fail=1
import json, sys
tag, base_p, cand_p, want_route = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
base = json.load(open(base_p))
cand = json.load(open(cand_p))
TOL = 1e-6

# (v10 path, v11 path, exact?) -- the v11 six-section refactor MOVED these; every value is
# identical at the new path across all 6 fixtures, so the pre-v11 baselines stay valid and keep
# their full assertion power. Each key is looked up at its v10 path, then its v11 path, so this
# oracle reads BOTH baseline generations and cannot rot the way it did when v11 landed without a
# consumer migration and left it throwing KeyError instead of comparing anything.
FIELDS = [
    ("deamination.d_max_5prime_raw_ungated", "channels.deamination.d_max_5prime_raw_ungated", False),
    ("deamination.d_max_3prime_raw_ungated", "channels.deamination.d_max_3prime_raw_ungated", False),
    ("deamination.lambda_5prime",            "channels.deamination.lambda_5prime",            False),
    ("channel_b.d_max",                      "diagnostics.channel_b.d_max",                   False),
    ("channel_b.valid",                      "diagnostics.channel_b.valid",                   True),
    ("channel_b.inverted",                   "diagnostics.channel_b.inverted",                True),
    # THE verdict-critical field. pi_estimate.{state,point} used to sit here, but v11 did not move
    # pi -- it REPLACED the estimator (shape-fit demoted to a diagnostic, co-occurrence promoted),
    # so its old values have no successor to assert against. The deamination CALL is what the science
    # turns on, it survived the refactor intact, and asserting it is what catches a negative control
    # flipping to "present". The ROUTE that produced the call is asserted separately, below.
    ("verdict.deamination.call", "characterization.verdict.deamination.call", True),
]

MISSING = object()

def get(d, path):
    for k in path.split("."):
        if not isinstance(d, dict) or k not in d:
            return MISSING
        d = d[k]
    return d

def get2(d, v10, v11):
    v = get(d, v10)
    return get(d, v11) if v is MISSING else v

bad = 0
for v10, v11, exact in FIELDS:
    b, c = get2(base, v10, v11), get2(cand, v10, v11)
    if b is MISSING or c is MISSING:
        print(f"[{tag}] ABSENT {v10} (v11: {v11}): baseline={b is not MISSING} candidate={c is not MISSING}")
        bad = 1
    elif exact or b is None or c is None:
        if b != c:
            print(f"[{tag}] MISMATCH {v10}: baseline={b!r} candidate={c!r}")
            bad = 1
    elif abs(float(b) - float(c)) > TOL:
        print(f"[{tag}] MISMATCH {v10}: baseline={b} candidate={c} (|d|={abs(float(b)-float(c)):.3e} > {TOL})")
        bad = 1

call = get2(cand, "verdict.deamination.call", "characterization.verdict.deamination.call")
via  = get2(cand, "verdict.deamination.authenticated_by", "characterization.verdict.deamination.authenticated_by")

# Assert the ROUTE, not just the call: an OR-combiner can reach the right call by the wrong path
# (a broken route and an over-firing one cancel in `call`). Compared against the script's table,
# not the baselines -- the baselines predate these routes.
if via is MISSING:
    print(f"[{tag}] ABSENT authenticated_by (expected {want_route!r})")
    bad = 1
elif via != want_route:
    print(f"[{tag}] ROUTE DRIFT authenticated_by: expected={want_route!r} got={via!r} (call={call!r})")
    bad = 1

if bad:
    print(f"[{tag}] FAIL (call={call}, authenticated_by={via})")
    sys.exit(1)
print(f"[{tag}] PASS (call={call}, authenticated_by={via})")
PY
    rm -f "$candidate"
done

if [ "$CAPTURE" = "1" ]; then
    [ "$fail" = 0 ] && echo "ORACLE CAPTURE COMPLETE" || echo "ORACLE CAPTURE INCOMPLETE"
    exit "$fail"
fi
[ "$fail" = 0 ] && echo "ORACLE PASS" || echo "ORACLE FAIL"
exit "$fail"
