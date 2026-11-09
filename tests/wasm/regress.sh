#!/bin/bash
# The regression decks through the WebAssembly build under Node. Wasm's maths
# library (musl) differs from glibc's in the last bit of some functions, so
# decks are either byte-identical to the reference or compared statistically
# (tests/stat/statcmp, measured errors, 0.25 % floor).
#   tests/wasm/regress.sh [deck names...]     default: all decks
# Needs emcc on the PATH (source emsdk_env.sh) and the reference outputs of
# tests/regress/run.sh (fast, NSCALE 4·10⁶) in $REF.
# Env: REF (default /tmp/nbregress56/ref-fast), WORK (default /tmp/nbwasm),
#      JOBS, TIMEOUT (s, default 2400).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
ref=${REF:-/tmp/nbregress56/ref-fast}
work=${WORK:-/tmp/nbwasm}
mkdir -p "$work/out"
echo "== build"
NSCALE=4000000 NIMCACHE="$work/nimcache" "$root/web/build.sh" cli "$work"
nim c --hints:off --warnings:off -d:release --nimcache:"$work/nimcache-statcmp" \
  -o:"$work/statcmp" "$root/tests/stat/statcmp.nim"

one() {   # deck file -> one result line
  local d=$(basename "$1" .in) rc=0 res
  timeout "${TIMEOUT:-2400}" node --stack-size=65500 "$work/magboltz.js" < "$1" \
    > "$work/out/$d.out" 2> "$work/out/$d.err" || rc=$?
  if [ ! -s "$ref/$d.out" ]; then echo "noref $d rc=$rc"
  elif cmp -s "$work/out/$d.out" "$ref/$d.out"; then echo "same  $d"
  elif [ $rc != 0 ]; then echo "FAIL  $d rc=$rc: $(tail -1 "$work/out/$d.err")"
  elif res=$("$work/statcmp" --floor=0.25 --errors="$root/tests/stat/spread-fast.tsv" --deck="$d" \
               "$ref/$d.out" "$work/out/$d.out" | head -1); then echo "stat  $d $res"
  else echo "FAIL  $d $res"; fi
}
export -f one; export work ref root
if [ $# -gt 0 ]; then decks=$(printf "$root/tests/regress/decks/%s.in\n" "$@")
else decks=$(ls "$root"/tests/regress/decks/*.in); fi
echo "$decks" | xargs -P "${JOBS:-$(nproc)}" -I{} bash -c 'one {}' | sort -k2 | tee "$work/result.txt"
echo "== $(grep -c '^same' "$work/result.txt") byte-identical, $(grep -c '^stat' "$work/result.txt")" \
     "statistically equal, $(grep -c '^FAIL' "$work/result.txt") failed"
! grep -q '^FAIL' "$work/result.txt"
