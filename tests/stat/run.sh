#!/bin/bash
# Statistical comparison with the reference Fortran build: the Nim build
# runs every deck with its own random numbers (NIMBOLTZ_SEED) and each
# result must agree with the cached Fortran result within ZMAX combined
# standard errors (statcmp.nim, with an error floor of FLOOR percent). The
# errors are the printed ones or, where larger, the spread measured over
# several random streams (spread-$mode.tsv, written by spread.sh): Magboltz
# underestimates some errors by up to a factor 16 (docs/DESIGN.md).
#
# Both builds run a copy of each deck with a fixed energy range: twice the
# EFINAL that Fortran chose for the original deck (fixefinal.nim). The
# automatic choice depends on the random numbers and is often too low in
# fast mode, which stops the run with ENERGY OUT OF RANGE ("abort" below).
# Decks where either run timed out are reported as "part" and not compared.
#
#   ./run.sh [deck-glob]          decks of tests/regress
#
# Env: as tests/regress (MAGBOLTZ_SRC, WORK, JOBS, FULL, TIMEOUT), plus
#      SEED (default 1), ZMAX (default 5), FLOOR (default 0.25 with a spread
#      table, else 1).
#      Calibration (docs/DESIGN.md): unchanged port, seed 1, fast mode:
#      90 of 96 decks pass; expected failures are SST fits whose printed
#      errors are too small (penning_neco2_t0; gas15 with other seeds).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
glob=${1:-*}
. "$here/../lib/reference.sh"
jobs=${JOBS:-$(nproc)}
seed=${SEED:-1}
zmax=${ZMAX:-5}
errors=$here/spread-$mode.tsv
if [ -f "$errors" ]; then floor=${FLOOR:-0.25}; else floor=${FLOOR:-1}; fi
out=$work/stat-$mode
mkdir -p "$out"
reference_setup

echo "== nim c"
nim c --hints:off --warnings:off -d:release -d:nscale=$nscale --nimcache:"$out/nimcache" -o:"$out/magboltz-nim" \
  "$root/src/nimboltz/magboltz_main.nim"
nim c --hints:off --warnings:off -d:release --nimcache:"$out/nimcache-statcmp" -o:"$out/statcmp" "$here/statcmp.nim"
nim c --hints:off --warnings:off -d:release --nimcache:"$out/nimcache-fixefinal" -o:"$out/fixefinal" "$here/fixefinal.nim"

orig=( "$root"/tests/regress/decks/$glob.in )
echo "== Fortran on the original decks (for their energy range)"
reference_cmds "${orig[@]}" | xargs -P "$jobs" -I{} sh -c "{}" || true
mkdir -p "$out/decks"
decks=()
for d in "${orig[@]}"; do
  n=$(basename "$d" .in)
  "$out/fixefinal" "$d" "$ref/$n.out" > "$out/decks/$n.in.new"
  cmp -s "$out/decks/$n.in.new" "$out/decks/$n.in" || mv "$out/decks/$n.in.new" "$out/decks/$n.in"
  rm -f "$out/decks/$n.in.new"
  decks+=( "$out/decks/$n.in" )
done
refout=$ref/stat
echo "== running ${#decks[@]} decks on $jobs cores, seed $seed"
{ reference_cmds "${decks[@]}"
  for d in "${decks[@]}"; do
    n=$(basename "$d" .in)
    echo "NIMBOLTZ_SEED=$seed timeout $to $out/magboltz-nim < $d > $out/$n.s$seed.out 2>$out/$n.s$seed.err; echo \$? > $out/$n.s$seed.rc"
  done
} | xargs -P "$jobs" -I{} sh -c "{}" || true

fail=0
partial=0
aborted=0
for d in "${decks[@]}"; do
  n=$(basename "$d" .in)
  if [ "$(cat "$refout/$n.rc" 2>/dev/null)" = 124 ] || [ "$(cat "$out/$n.s$seed.rc")" = 124 ]; then
    partial=$((partial+1)); echo "part  $n (timeout)"; continue
  fi
  opts=( --z=$zmax --floor=$floor )
  [ -f "$errors" ] && opts+=( --errors="$errors" --deck=$n )
  rc=0; res=$("$out/statcmp" "${opts[@]}" "$refout/$n.out" "$out/$n.s$seed.out") || rc=$?
  case $rc in
    0) echo "ok    $n: $res" ;;
    2) aborted=$((aborted+1)); echo "abort $n: $res" ;;
    *) fail=$((fail+1)); echo "FAIL  $n: $res" ;;
  esac
done
echo "== ${#decks[@]} decks, $fail failed, $aborted aborted, $partial not compared (timeout ${to}s), |z| <= $zmax, floor $floor%"
[ $fail = 0 ]
