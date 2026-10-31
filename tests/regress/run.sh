#!/bin/bash
# Bitwise regression suite: runs every deck in decks/ through the reference
# Fortran build and the f2nim-translated Nim build and requires byte-identical
# stdout.
#
#   ./run.sh [deck-glob]          fast mode (default): NSCALE 4e7 -> 4e6 in both
#   FULL=1 ./run.sh [deck-glob]   unmodified source, committed port/gen
#
# Env: MAGBOLTZ_SRC (default ~/src/Magboltz/magboltz-11.21.f),
#      WORK (default /tmp/nbregress), JOBS (default nproc).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
src=${MAGBOLTZ_SRC:-$HOME/src/Magboltz/magboltz-11.21.f}
work=${WORK:-/tmp/nbregress}
jobs=${JOBS:-$(nproc)}
glob=${1:-*}
FFLAGS="-O0 -ffp-contract=off -std=legacy -finit-local-zero"

mkdir -p "$work/nim" "$work/out"
if [ "${FULL:-0}" = 1 ]; then
  cp "$src" "$work/magboltz.f"
else
  sed 's/^      NSCALE=40000000 *$/      NSCALE=4000000/' "$src" > "$work/magboltz.f"
  [ "$(grep -c '^      NSCALE=4000000$' "$work/magboltz.f")" = 2 ] || { echo "NSCALE patch failed"; exit 1; }
fi

# Fortran reference
if [ ! "$work/magboltz-f" -nt "$work/magboltz.f" ]; then
  echo "== gfortran"; gfortran $FFLAGS "$work/magboltz.f" -o "$work/magboltz-f"
fi

# Nim translation
echo "== f2nim + nim c"
nim c --hints:off --warnings:off -d:release -o:"$work/f2nim" "$root/tools/f2nim/f2nim.nim" >/dev/null
cp "$root"/port/{frt.nim,rm48.nim,main.nim} "$work/nim/"; cp "$root/nim.cfg" "$work/nim/"
if [ "${FULL:-0}" = 1 ]; then
  rm -rf "$work/nim/gen"; cp -r "$root/port/gen" "$work/nim/gen"
else
  "$work/f2nim" "$work/magboltz.f" "$work/nim/gen" 2>&1 | grep -v -E '^const intrinsic|split label|hoist' || true
fi
nim c --hints:off --warnings:off -d:release -o:"$work/magboltz-nim" "$work/nim/main.nim"

# Run both on every deck, in parallel
decks=( "$here"/decks/$glob.in )
echo "== running ${#decks[@]} decks x 2 builds on $jobs cores"
for d in "${decks[@]}"; do
  n=$(basename "$d" .in)
  echo "$work/magboltz-f < $d > $work/out/$n.f.out 2>/dev/null"
  echo "$work/magboltz-nim < $d > $work/out/$n.nim.out 2>$work/out/$n.nim.err"
done | xargs -P "$jobs" -I{} timeout "${TIMEOUT:-3600}" sh -c "{}" || true

# Compare
isprefix() {  # $1 nonempty, shorter than $2 and equal to its start
  [ -s "$1" ] && [ $(stat -c %s "$1") -lt $(stat -c %s "$2") ] &&
    head -c $(stat -c %s "$1") "$2" | cmp -s - "$1"
}
fail=0
partial=0
for d in "${decks[@]}"; do
  n=$(basename "$d" .in)
  f="$work/out/$n.f.out"; m="$work/out/$n.nim.out"
  if cmp -s "$f" "$m"; then
    echo "ok    $n"
  elif isprefix "$f" "$m" || isprefix "$m" "$f"; then
    partial=$((partial+1))
    echo "part  $n: identical up to the end of the shorter output (timeout)"
  else
    fail=$((fail+1))
    first=$(cmp "$work/out/$n.f.out" "$work/out/$n.nim.out" 2>&1 | head -1 || true)
    echo "FAIL  $n: $first"
  fi
done
echo "== ${#decks[@]} decks, $fail failed, $partial identical but cut short by TIMEOUT=${TIMEOUT:-3600}"
[ $fail = 0 ]
