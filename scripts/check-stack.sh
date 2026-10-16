#!/usr/bin/env bash
# Check that every libtaph consumer is on the same libtaph commit. Read-only.
#
# Compares: the LIBTAPH_PIN of DART and fqdup, the shared pinned worktree's HEAD (and that it
# is clean), that the pin is reachable from libtaph's origin (so CI can fetch it), and the
# libtaph SHA each built binary reports in `--version`. Exit 0 when all agree, 1 otherwise.
#
# Environment (defaults in brackets):
#   REPOS           [/maps/projects/fernandezguerra/apps/repos]
#   LIBTAPH_DIR     [$REPOS/libtaph]
#   LIBTAPH_PINNED  [/projects/caeg/scratch/kbd606/tmp/libtaph_pinned]
#   DART_DIR        [$REPOS/dart]            DART_BIN   [$DART_DIR/build/dart]
#   FQDUP_DIR       [$REPOS/fqdup]           FQDUP_BIN  [$FQDUP_DIR/build/fqdup]
#   NO_FETCH=1      skip `git fetch` of libtaph origin
set -uo pipefail
REPOS=${REPOS:-/maps/projects/fernandezguerra/apps/repos}
LIBTAPH_DIR=${LIBTAPH_DIR:-$REPOS/libtaph}
LIBTAPH_PINNED=${LIBTAPH_PINNED:-/projects/caeg/scratch/kbd606/tmp/libtaph_pinned}
DART_DIR=${DART_DIR:-$REPOS/dart}
DART_BIN=${DART_BIN:-$DART_DIR/build/dart}
FQDUP_DIR=${FQDUP_DIR:-$REPOS/fqdup}
FQDUP_BIN=${FQDUP_BIN:-$FQDUP_DIR/build/fqdup}

fail=0
row() { printf '%-28s %-42s %s\n' "$1" "$2" "$3"; }
bad() { row "$1" "$2" "MISMATCH${3:+: $3}"; fail=1; }

pin_of() { [ -f "$1/LIBTAPH_PIN" ] && tr -d '[:space:]' < "$1/LIBTAPH_PIN"; }
bin_sha() { [ -x "$1" ] && "$1" --version 2>/dev/null | awk '$1=="libtaph"{print $2; exit}'; }

dart_pin=$(pin_of "$DART_DIR")
fqdup_pin=$(pin_of "$FQDUP_DIR")
ref=${dart_pin:-$fqdup_pin}
[ -n "$ref" ] || { echo "no LIBTAPH_PIN in $DART_DIR or $FQDUP_DIR" >&2; exit 1; }

row component value status
if [ "$dart_pin" = "$ref" ]; then
    row "DART LIBTAPH_PIN" "$dart_pin" "ok ($(git -C "$DART_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null))"
else
    bad "DART LIBTAPH_PIN" "${dart_pin:-missing}" "branch $(git -C "$DART_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null)"
fi
if [ "$fqdup_pin" = "$ref" ]; then
    row "fqdup LIBTAPH_PIN" "$fqdup_pin" "ok ($(git -C "$FQDUP_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null))"
else
    bad "fqdup LIBTAPH_PIN" "${fqdup_pin:-missing}" "branch $(git -C "$FQDUP_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null)"
fi

pinned_head=$(git -C "$LIBTAPH_PINNED" rev-parse HEAD 2>/dev/null)
pinned_dirty=$(git -C "$LIBTAPH_PINNED" status --porcelain 2>/dev/null | wc -l)
if [ "$pinned_head" = "$ref" ] && [ "$pinned_dirty" -eq 0 ]; then
    row "libtaph_pinned HEAD" "$pinned_head" "ok, clean"
else
    bad "libtaph_pinned HEAD" "${pinned_head:-missing}" "dirty files: $pinned_dirty"
fi

[ "${NO_FETCH:-0}" = 1 ] || git -C "$LIBTAPH_DIR" fetch -q origin 2>/dev/null
if git -C "$LIBTAPH_DIR" merge-base --is-ancestor "$ref" origin/main 2>/dev/null; then
    row "pin on origin/main" "$ref" "ok"
elif [ -n "$(git -C "$LIBTAPH_DIR" branch -r --contains "$ref" 2>/dev/null; git -C "$LIBTAPH_DIR" tag --contains "$ref" 2>/dev/null)" ]; then
    row "pin on origin" "$ref" "ok (not on main)"
else
    bad "pin on origin" "$ref" "not reachable from any origin ref or tag; CI cannot fetch it"
fi

for pair in "dart:$DART_BIN" "fqdup:$FQDUP_BIN"; do
    name=${pair%%:*}; bin=${pair#*:}
    sha=$(bin_sha "$bin")
    if [ -z "$sha" ]; then bad "$name binary" "-" "no libtaph line from $bin --version"
    elif [ "$sha" = "$ref" ]; then row "$name binary" "$sha" "ok"
    else bad "$name binary" "$sha" "rebuild $bin"; fi
done

[ $fail -eq 0 ] && echo "stack in sync on libtaph $ref" || echo "stack OUT OF SYNC" >&2
exit $fail
