#!/usr/bin/env bash
# Move every libtaph consumer (DART, fqdup) to one libtaph commit, test them together, commit.
#
#   scripts/bump-libtaph.sh <libtaph-sha>
#
# 1. The SHA must be reachable from libtaph's origin (CI fetches it by SHA).
# 2. LIBTAPH_PINNED (the shared clean worktree local builds use) is checked out at it.
# 3. Each consumer gets a new worktree on branch libtaph-bump-<short> from its current HEAD
#    (the user's checkouts are not touched), LIBTAPH_PIN is rewritten, and the consumer's own
#    scripts/libtaph-consumer-test.sh <build_dir> <libtaph_dir> configures, builds and runs ctest.
# 4. fqdup's scripts/make-tdm1-fixture.sh regenerates the TDM1 fixture with the new fqdup build;
#    it is copied into DART's tests/data/tdm1 and DART's library_profile_roundtrip is rerun.
# 5. Each branch commits LIBTAPH_PIN (and the DART fixture). Nothing is pushed or merged.
# 6. check-stack.sh runs against the new branches and binaries.
#
# Stops at the first failure and leaves worktrees/builds in BUILD_ROOT for inspection.
# Environment: REPOS, LIBTAPH_DIR, LIBTAPH_PINNED, DART_DIR, FQDUP_DIR as in check-stack.sh;
# DART_BASE / FQDUP_BASE [HEAD]: ref each bump branch starts from (it must carry LIBTAPH_PIN and
# scripts/libtaph-consumer-test.sh); BUILD_ROOT [/path/to/tmp/libtaph-bump-<short>];
# DART_CMAKE_ARGS, FQDUP_* etc. pass through to the consumer scripts.
set -euo pipefail
[ $# -eq 1 ] || { echo "usage: $0 <libtaph-sha>" >&2; exit 2; }
HERE=$(cd "$(dirname "$0")" && pwd)
REPOS=${REPOS:-/path/to/repos}
LIBTAPH_DIR=${LIBTAPH_DIR:-$REPOS/libtaph}
LIBTAPH_PINNED=${LIBTAPH_PINNED:-/path/to/tmp/libtaph_pinned}
DART_DIR=${DART_DIR:-$REPOS/dart}
FQDUP_DIR=${FQDUP_DIR:-$REPOS/fqdup}

git -C "$LIBTAPH_DIR" fetch -q origin
SHA=$(git -C "$LIBTAPH_DIR" rev-parse --verify "$1^{commit}")
SHORT=${SHA:0:7}
BRANCH=libtaph-bump-$SHORT
BUILD_ROOT=${BUILD_ROOT:-/path/to/tmp/libtaph-bump-$SHORT}
log() { printf '== %s\n' "$*"; }

if ! git -C "$LIBTAPH_DIR" merge-base --is-ancestor "$SHA" origin/main &&
   [ -z "$(git -C "$LIBTAPH_DIR" branch -r --contains "$SHA"; git -C "$LIBTAPH_DIR" tag --contains "$SHA")" ]; then
    echo "$SHA is not reachable from libtaph's origin; push it before pinning consumers to it" >&2
    exit 1
fi
[ -z "$(git -C "$LIBTAPH_PINNED" status --porcelain)" ] || { echo "$LIBTAPH_PINNED is dirty" >&2; exit 1; }
mkdir -p "$BUILD_ROOT"

log "libtaph_pinned -> $SHA (was $(git -C "$LIBTAPH_PINNED" rev-parse --short HEAD))"
git -C "$LIBTAPH_PINNED" checkout -q --detach "$SHA"

declare -A SRC
for name in dart fqdup; do
    repo=$([ $name = dart ] && echo "$DART_DIR" || echo "$FQDUP_DIR")
    base=$([ $name = dart ] && echo "${DART_BASE:-HEAD}" || echo "${FQDUP_BASE:-HEAD}")
    wt=$BUILD_ROOT/$name-src
    git -C "$repo" show-ref --quiet "refs/heads/$BRANCH" && { echo "$name: branch $BRANCH exists" >&2; exit 1; }
    log "$name: worktree $wt on $BRANCH from $base ($(git -C "$repo" rev-parse --short "$base"))"
    git -C "$repo" worktree add -q -b "$BRANCH" "$wt" "$base"
    [ -x "$wt/scripts/libtaph-consumer-test.sh" ] || { echo "$name: no scripts/libtaph-consumer-test.sh on $BRANCH" >&2; exit 1; }
    echo "$SHA" > "$wt/LIBTAPH_PIN"
    SRC[$name]=$wt
done

for name in fqdup dart; do
    log "$name: build + ctest against $SHORT"
    "${SRC[$name]}/scripts/libtaph-consumer-test.sh" "$BUILD_ROOT/$name-build" "$LIBTAPH_PINNED"
done

log "fqdup: regenerate TDM1 fixture"
"${SRC[fqdup]}/scripts/make-tdm1-fixture.sh" "$BUILD_ROOT/fqdup-build/fqdup" "$BUILD_ROOT/fixture"
mkdir -p "${SRC[dart]}/tests/data/tdm1"
cp "$BUILD_ROOT/fixture/fqdup_sample.tdm1" "$BUILD_ROOT/fixture/fqdup_sample.expected.tsv" "${SRC[dart]}/tests/data/tdm1/"
log "dart: library_profile_roundtrip with the new fixture"
(cd "$BUILD_ROOT/dart-build" && ctest -R library_profile_roundtrip --output-on-failure)

git -C "${SRC[fqdup]}" add LIBTAPH_PIN
git -C "${SRC[fqdup]}" diff --cached --quiet || git -C "${SRC[fqdup]}" commit -q -m "Bump libtaph to $SHORT"
git -C "${SRC[dart]}" add LIBTAPH_PIN tests/data/tdm1
git -C "${SRC[dart]}" diff --cached --quiet ||
    git -C "${SRC[dart]}" commit -q -m "Bump libtaph to $SHORT; regenerate the fqdup TDM1 fixture"
log "committed: dart $(git -C "${SRC[dart]}" rev-parse --short HEAD), fqdup $(git -C "${SRC[fqdup]}" rev-parse --short HEAD) on $BRANCH (not pushed)"

NO_FETCH=1 DART_DIR=${SRC[dart]} FQDUP_DIR=${SRC[fqdup]} \
    DART_BIN=$BUILD_ROOT/dart-build/dart FQDUP_BIN=$BUILD_ROOT/fqdup-build/fqdup \
    "$HERE/check-stack.sh"
cat <<EOF
Next: review and merge $BRANCH in each repo, then rebuild the usual build dirs.
Roll back: git -C <repo> worktree remove $BUILD_ROOT/<name>-src && git -C <repo> branch -D $BRANCH;
           git -C $LIBTAPH_PINNED checkout --detach <previous sha>
EOF
