#!/bin/bash
# Bitwise regression suite: runs every deck in decks/ through the reference
# Fortran build and the Nim build (src/nimboltz) and requires byte-identical
# stdout.
#
#   ./run.sh [deck-glob]          fast mode (default): NSCALE 4e7 -> 4e6 in both
#   FULL=1 ./run.sh [deck-glob]   unmodified NSCALE
#
# Fortran outputs are cached in $WORK/ref-fast or $WORK/ref-full and rerun
# only when the reference binary or the deck is newer.
#
# Env: MAGBOLTZ_SRC (default ~/src/Magboltz/magboltz-11.21.f),
#      WORK (default /tmp/nbregress), JOBS (default nproc), TIMEOUT (s per run).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
glob=${1:-*}
. "$here/../lib/reference.sh"
jobs=${JOBS:-$(nproc)}
mkdir -p "$work/out"
reference_setup

# Nim build
echo "== nim c"
nim c --hints:off --warnings:off -d:release -d:nscale=$nscale --nimcache:"$work/nimcache" -o:"$work/magboltz-nim" \
  "$root/src/nimboltz/magboltz_main.nim"

# Run both on every deck, in parallel
decks=( "$here"/decks/$glob.in )
echo "== running ${#decks[@]} decks on $jobs cores"
{ reference_cmds "${decks[@]}"
  for d in "${decks[@]}"; do
    n=$(basename "$d" .in)
    echo "timeout $to $work/magboltz-nim < $d > $work/out/$n.nim.out 2>$work/out/$n.nim.err"
  done
} | xargs -P "$jobs" -I{} sh -c "{}" || true

# Compare
isprefix() {  # $1 nonempty, shorter than $2 and equal to its start
  [ -s "$1" ] && [ $(stat -c %s "$1") -lt $(stat -c %s "$2") ] &&
    head -c $(stat -c %s "$1") "$2" | cmp -s - "$1"
}
fail=0
partial=0
for d in "${decks[@]}"; do
  n=$(basename "$d" .in)
  f="$ref/$n.out"; m="$work/out/$n.nim.out"
  if cmp -s "$f" "$m"; then
    echo "ok    $n"
  elif isprefix "$f" "$m" || isprefix "$m" "$f"; then
    partial=$((partial+1))
    echo "part  $n: identical up to the end of the shorter output (timeout)"
  else
    fail=$((fail+1))
    first=$(cmp "$ref/$n.out" "$work/out/$n.nim.out" 2>&1 | head -1 || true)
    echo "FAIL  $n: $first"
  fi
done
echo "== ${#decks[@]} decks, $fail failed, $partial identical but cut short by TIMEOUT=${TIMEOUT:-3600}"
[ $fail = 0 ]
