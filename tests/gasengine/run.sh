#!/bin/sh
# Compares the generic gas engine (src/nimboltz/gases) with the per-gas
# modules of commit $REF, which match the Magboltz fork bit for bit
# (tests/gasdump), on 96 cases per gas (the 72 of tests/gasdump plus
# EFINAL 1 and 15 eV).
#   tests/gasengine/run.sh [-t tol] [gas...]    (default: all gases, 1e-12)
# env: REF (default f336489), WORK (/tmp/gasengine)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
ref=${REF:-f336489}
work=${WORK:-/tmp/gasengine}
rm -rf "$work/reference"
mkdir -p "$work/reference"
git -C "$root" archive "$ref" src/nimboltz/gases src/nimboltz/gases.nim |
  tar -x -C "$work/reference"
mv "$work/reference/src/nimboltz/gases" "$work/reference/gases"
mv "$work/reference/src/nimboltz/gases.nim" "$work/reference/refgases.nim"
rm -r "$work/reference/src"
sed -i 's#^import ../frt, ../mixture, ../tables#import nimboltz/frt, refxsec#; s#^export math, mixture, tables#export math, refxsec#' \
  "$work/reference/gases/common.nim"
sed -i 's#\.\./frt#nimboltz/frt#' "$work/reference/gases/"*.nim
nim c --hints:off --warnings:off -d:release --path:"$root/src" --path:"$work/reference" --path:"$here" \
  --nimcache:"$work/nimcache" -o:"$work/compare" "$here/compare.nim"
tol=""
if [ "$1" = "-t" ]; then tol="-t $2"; shift 2; fi
"$work/compare" $tol ${*:-$(seq 1 80)}
