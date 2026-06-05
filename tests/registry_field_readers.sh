#!/usr/bin/env bash
# Architecture A-bar (wire-iff-consumer): every field declared in ChannelSpec must have at least one
# READER outside the registry header. A field with zero readers is a decorative constexpr promise the
# count-table gate cannot witness (the same bug class as cap-masking) -> wire it or delete it. This
# test fails CI otherwise, so the rule is self-enforcing rather than a judgement call.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
REG="$ROOT/include/taph/channel_registry.hpp"
[[ -f "$REG" ]] || { echo "SKIP: registry header not found ($REG)"; exit 0; }

# Field names declared inside `struct ChannelSpec { ... };`
fields=$(awk '/struct ChannelSpec[[:space:]]*\{/{f=1;next} f&&/^\};/{f=0} f' "$REG" \
         | sed -n 's/^[[:space:]]*[A-Za-z_][A-Za-z0-9_:<>* ]*[[:space:]*&]\([A-Za-z_][A-Za-z0-9_]*\)[[:space:]]*;.*/\1/p')

fail=0; n=0
for fld in $fields; do
  n=$((n+1))
  # readers = member accesses `.<field>` anywhere in src/ + include/, EXCLUDING the registry header itself
  readers=$(grep -rn --include='*.cpp' --include='*.hpp' "\.${fld}\b" "$ROOT/src" "$ROOT/include" 2>/dev/null \
            | grep -v "channel_registry.hpp" | wc -l)
  if [[ "$readers" -ge 1 ]]; then
    printf '  OK   %-22s %s reader(s)\n' "$fld" "$readers"
  else
    printf '  FAIL %-22s 0 readers (decorative -> wire or delete)\n' "$fld"
    fail=1
  fi
done
echo "--- ChannelSpec field-reader gate: $n fields checked, result=$([[ $fail -eq 0 ]] && echo PASS || echo FAIL) ---"
exit $fail
