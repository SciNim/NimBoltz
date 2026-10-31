#!/bin/bash
# Bitwise check of every GAS routine: gasdmp.f (reference Fortran) and
# gasdmp.nim (translated port) are fed the same cases and their raw dumps
# are compared case by case (one dump is ~19 MB).
#   ./run.sh [gas numbers...]     default: 1..61
# Env: MAGBOLTZ_SRC, WORK (default /tmp/gasdump), JOBS.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
src=${MAGBOLTZ_SRC:-$HOME/src/Magboltz/magboltz-11.21.f}
work=${WORK:-/tmp/gasdump}
jobs=${JOBS:-$(nproc)}
mkdir -p "$work/cases"

# library = everything but PROGRAM MAGBOLTZ (lines 472..625 of 11.21)
start=$(grep -n '^      PROGRAM MAGBOLTZ' "$src" | cut -d: -f1)
end=$(awk -v s="$start" 'NR>s && /^      SUBROUTINE ANGCUT/{print NR-1; exit}' "$src")
awk -v s="$start" -v e="$end" 'NR<s || NR>e' "$src" > "$work/mblib.f"
gfortran -O0 -ffp-contract=off -std=legacy -finit-local-zero "$here/gasdmp.f" "$work/mblib.f" -o "$work/gasdmp-f" 2>&1 | grep -v -i warning || true
nim c --hints:off --warnings:off -d:release --path:"$root/src" --nimcache:"$work/nimcache" -o:"$work/gasdmp-nim" "$here/gasdmp.nim"
nim c --hints:off --warnings:off -d:release --nimcache:"$work/nimcache-dumpdiff" -o:"$work/dumpdiff" "$here/dumpdiff.nim"

gases=${*:-$(seq 1 61)}
for g in $gases; do
  : > "$work/cases/$g.in"
  for ef in 0.3 4.0 40.0 400.0 8000.0 200000.0; do
    for an in 0 1 2; do
      for t in 20.0 -50.0; do
        for fr in 100.0 10.0; do
          echo "$g $an $ef $t 760.0 $fr" >> "$work/cases/$g.in"
        done
      done
    done
  done
done

check() {   # gas number -> one result line; details in $work/gas<g>.log
  local g=$1 d=$work/tmp.$1 log=$work/gas$1.log i=0 nfail=0 l r
  mkdir -p $d; : > $log
  while read -r l; do
    i=$((i+1))
    echo "$l" > $d/c.in
    "$work/gasdmp-f" $d/f.bin < $d/c.in > $d/f.out 2>/dev/null
    if ! "$work/gasdmp-nim" $d/n.bin < $d/c.in > $d/n.out 2>$d/n.err; then
      nfail=$((nfail+1)); echo "case $i [$l]: nim failed: $(tail -1 $d/n.err)" >> $log; continue
    fi
    if ! cmp -s $d/f.bin $d/n.bin; then
      nfail=$((nfail+1)); echo "case $i [$l]:" >> $log
      "$work/dumpdiff" $d/f.bin $d/n.bin 5 >> $log
    elif ! cmp -s $d/f.out $d/n.out; then
      nfail=$((nfail+1)); echo "case $i [$l]: stdout differs" >> $log
    fi
  done < "$work/cases/$g.in"
  rm -rf $d
  if [ $nfail = 0 ]; then echo "ok    gas$g ($i cases)"; else echo "FAIL  gas$g: $nfail of $i cases, see $log"; fi
}
export -f check; export work
printf '%s\n' $gases | xargs -P "$jobs" -I{} bash -c 'check {}' | sort -V | tee "$work/result.txt"
n=$(grep -c FAIL "$work/result.txt" || true)
echo "== $(wc -l < "$work/result.txt") gases x $(wc -l < "$work/cases/1.in" 2>/dev/null || echo ?) cases, $n failed"
[ "$n" = 0 ]
