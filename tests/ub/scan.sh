#!/bin/bash
# Undefined-behaviour scan of a Magboltz source (default: the fixes branch
# checkout). Two checks, both on the reference Fortran only:
#  1. out-of-bounds: build with -fcheck=bounds, run every regression deck
#     and a grid of per-gas cases (gasdmp), collect "At line ..." errors;
#  2. uninitialised locals: build with zero-initialised and with poisoned
#     locals (NaN reals, huge integers), run every deck, compare outputs.
# Decks run with NSCALE 4e7 -> 4e6 (as in the fast regression mode).
#   ./scan.sh            Env: MAGBOLTZ_SRC, WORK (default /tmp/ubscan), JOBS
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
src=${MAGBOLTZ_SRC:-$HOME/src/Magboltz/magboltz-11.21.f}
work=${WORK:-/tmp/ubscan}
jobs=${JOBS:-$(nproc)}
rm -rf "$work"; mkdir -p "$work/out"
sed 's/NSCALE=40000000/NSCALE=4000000/' "$src" > "$work/magboltz.f"
start=$(grep -n '^      PROGRAM MAGBOLTZ' "$work/magboltz.f" | cut -d: -f1)
end=$(awk -v s="$start" 'NR>s && /^      SUBROUTINE ANGCUT/{print NR-1; exit}' "$work/magboltz.f")
awk -v s="$start" -v e="$end" 'NR<s || NR>e' "$work/magboltz.f" > "$work/mblib.f"
F="gfortran -O0 -ffp-contract=off -std=legacy -w"
$F -fcheck=bounds -g "$work/magboltz.f" -o "$work/mb-bounds"
$F -fcheck=bounds -g "$root/tests/gasdump/gasdmp.f" "$work/mblib.f" -o "$work/gasdmp-bounds"
$F -finit-local-zero "$work/magboltz.f" -o "$work/mb-zero"
$F -finit-real=snan -finit-integer=-1234567 -finit-logical=true -finit-character=63 \
   "$work/magboltz.f" -o "$work/mb-poison"
cd "$work"
{
for d in "$root"/tests/regress/decks/*.in; do
  n=$(basename "$d" .in)
  for b in bounds zero poison; do echo "./mb-$b < $d > out/$n.$b 2> out/$n.$b.err"; done
done
for g in $(seq 1 61); do
  for ef in 0.3 4.0 40.0 400.0 8000.0 200000.0; do for an in 0 1 2; do for fr in 100.0 10.0; do
    echo "echo '$g $an $ef 20.0 760.0 $fr' | ./gasdmp-bounds /dev/null > /dev/null 2> out/g${g}_${an}_${ef}_$fr.err"
  done; done; done
done
} | xargs -P "$jobs" -I{} timeout "${TIMEOUT:-7200}" sh -c "{}" || true
{
echo "== out-of-bounds (count, first error per run)"
grep -h -A1 'At line' out/*.err | grep -v '^--' | paste - - | sed 's/Fortran runtime error: //' | sort | uniq -c | sort -rn
echo "== zero vs poison initialisation"
for f in out/*.zero; do n=${f%.zero}; cmp -s "$f" "$n.poison" || echo "differs: $(basename $n)"; done
echo "== runs that did not finish normally"
for f in out/*.bounds out/*.zero out/*.poison; do grep -q 'RM48' "$f" 2>/dev/null || echo "no output: $f"; done | head -20
} > report.txt
cat report.txt
