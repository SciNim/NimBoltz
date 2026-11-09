#!/bin/bash
# Builds NimBoltz for WebAssembly with Emscripten (emcc on the PATH, e.g.
# `source emsdk_env.sh`):
#   web/build.sh          the web page's module: web/nimboltz.{js,wasm}
#   web/build.sh cli DIR  the program, to run decks under Node:
#                         node DIR/magboltz.js < deck.in
# NSCALE (default 4·10⁷) is Magboltz's collisions per unit of NMAX.
set -euo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
nscale=${NSCALE:-40000000}
nimcache=${NIMCACHE:-/tmp/nimboltz-wasm}
# One thread: the avalanches' threads need SharedArrayBuffer in a browser.
common=(--hints:off --warnings:off -d:release -d:nscale=$nscale --threads:off
        --os:linux --cpu:wasm32 --cc:clang --clang.exe:emcc --clang.linkerexe:emcc
        --clang.options.speed:-O3 -d:emscripten)
mem="-sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=256MB -sMAXIMUM_MEMORY=4GB -sSTACK_SIZE=16MB"
if [ "${1:-}" = cli ]; then
  out=${2:?output directory}
  mkdir -p "$out"
  nim c "${common[@]}" --nimcache:"$nimcache/cli" \
    --passL:"$mem -sNODERAWFS=1 -sEXIT_RUNTIME=1" \
    -o:"$out/magboltz.js" "$here/src/nimboltz/magboltz_main.nim"
else
  nim c "${common[@]}" --nimcache:"$nimcache/web" \
    --passL:"$mem -sENVIRONMENT=web,worker,node -sEXIT_RUNTIME=0 \
-sEXPORTED_FUNCTIONS=_main,_nimboltzRun,_nimboltzGases,_nimboltzCollisionScale \
-sEXPORTED_RUNTIME_METHODS=ccall,UTF8ToString,HEAP8 -sMODULARIZE=1 -sEXPORT_NAME=createNimBoltz" \
    -o:"$here/web/nimboltz.js" "$here/src/nimboltz/wasm.nim"
fi
