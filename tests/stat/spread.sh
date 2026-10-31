#!/bin/bash
# Measures the real statistical error of every printed result: runs the
# fixed-EFINAL decks of run.sh with SEEDS more random streams of the Nim
# build and writes, per deck and result, the spread over those runs and the
# Fortran reference (spread.nim) to spread-$mode.tsv, which run.sh uses
# instead of the printed errors where they are smaller.
#
#   ./spread.sh [deck-glob]      (run ./run.sh first: it writes the decks)
#
# Env: as run.sh, plus SEEDS (default "2 3 4 5 6 7 8 9"). Only meaningful
# for a build that is statistically equal to Fortran (the unchanged port).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
glob=${1:-*}
. "$here/../lib/reference.sh"
jobs=${JOBS:-$(nproc)}
seeds=${SEEDS:-2 3 4 5 6 7 8 9}
out=$work/stat-$mode
spread=$out/spread
mkdir -p "$spread"
reference_setup

echo "== nim c"
nim c --hints:off --warnings:off -d:release -d:nscale=$nscale --nimcache:"$out/nimcache" -o:"$out/magboltz-nim" \
  "$root/src/nimboltz/magboltz_main.nim"
nim c --hints:off --warnings:off -d:release --nimcache:"$out/nimcache-spread" -o:"$out/spread-bin" "$here/spread.nim"

decks=( "$out"/decks/$glob.in )
echo "== running ${#decks[@]} decks x seeds $seeds on $jobs cores"
for s in $seeds; do
  for d in "${decks[@]}"; do
    n=$(basename "$d" .in)
    echo "NIMBOLTZ_SEED=$s timeout $to $out/magboltz-nim < $d > $spread/$n.s$s.out 2>/dev/null; echo \$? > $spread/$n.s$s.rc"
  done
done | xargs -P "$jobs" -I{} sh -c "{}" || true

tsv=$here/spread-$mode.tsv
: > "$tsv.new"
skipped=0
for d in "${decks[@]}"; do
  n=$(basename "$d" .in)
  runs=()
  [ "$(cat "$ref/stat/$n.rc" 2>/dev/null)" = 0 ] && runs+=( "$ref/stat/$n.out" )
  for s in $seeds; do
    [ "$(cat "$spread/$n.s$s.rc")" = 0 ] && runs+=( "$spread/$n.s$s.out" )
  done
  if [ ${#runs[@]} -lt 4 ]; then skipped=$((skipped+1)); echo "skip  $n (${#runs[@]} complete runs)"; continue; fi
  "$out/spread-bin" "$n" "${runs[@]}" >> "$tsv.new"
done
# rows of decks not run this time are kept
if [ -f "$tsv" ]; then
  for d in "${decks[@]}"; do basename "$d" .in; done > "$tsv.done"
  awk -F'\t' 'NR == FNR { done[$1]; next } !($1 in done)' "$tsv.done" "$tsv" >> "$tsv.new"
  rm -f "$tsv.done"
fi
mv "$tsv.new" "$tsv"
echo "== $tsv: $(cut -f1 "$tsv" | sort -u | wc -l) decks, $skipped skipped"
"$out/spread-bin" --summary "$tsv"
