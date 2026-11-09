#!/bin/bash
# WebGPU walkers against CPU walkers (tests/wasm/twebgpu.nim) in headless
# Chromium: TATTLETALE=dir (tattletale root) tests/wasm/webgpu.sh [OUT]
# Needs emcc on the PATH, and playwright's Chromium (twebgpu.mjs).
set -euo pipefail
here=$(cd "$(dirname "$0")/../.." && pwd)
out=${1:-/tmp/nimboltz-twebgpu}
mkdir -p "$out"
nim c --hints:off --warnings:off -d:release -d:nscale=4000000 --threads:off \
  --os:linux --cpu:wasm32 --cc:clang --clang.exe:emcc --clang.linkerexe:emcc \
  --clang.options.speed:-O3 -d:emscripten -d:nimboltzWebGpu -d:wgslPointerParams -d:wgslPreciseMath ${NIMFLAGS:-} \
  --path:"$here/src" --path:"${TATTLETALE:?the tattletale root}" --nimcache:"$out/nc" \
  --passL:"-sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=512MB -sMAXIMUM_MEMORY=4GB -sSTACK_SIZE=16MB \
--js-library=$here/web/webgpu.js -sASYNCIFY=1 -sASYNCIFY_STACK_SIZE=65536 -sASYNCIFY_IGNORE_INDIRECT=1 -sENVIRONMENT=web \
-sEXPORTED_FUNCTIONS=_main,_twebgpuRun -sEXPORTED_RUNTIME_METHODS=ccall -sMODULARIZE=1 \
-sEXPORT_NAME=createTest" -o:"$out/twebgpu.js" "$here/tests/wasm/twebgpu.nim"
cat > "$out/twebgpu.html" <<'HTML'
<!doctype html><meta charset="utf-8"><script src="twebgpu.js"></script><script>
createTest({print: (t) => console.log(t), printErr: (t) => console.log(t)}).then(async (m) => {
  try {
    const a = await navigator.gpu.requestAdapter();
    m.gpuDevice = await a.requestDevice({requiredLimits: {maxStorageBufferBindingSize: a.limits.maxStorageBufferBindingSize,
                                                          maxBufferSize: a.limits.maxBufferSize}});
    m.gpuDevice.addEventListener('uncapturederror', (e) => { m.gpuError = e.error.message; console.log('WebGPU: ' + e.error.message); });
    await m.ccall('twebgpuRun', null, [], [], {async: true});
  } catch (e) { console.log('FAILED ' + e); }
});
</script>
HTML
port=$((8800 + RANDOM % 100))
python3 -m http.server -d "$out" $port > /dev/null 2>&1 &
server=$!
trap 'kill $server' EXIT
sleep 1
node "$here/tests/wasm/twebgpu.mjs" "http://127.0.0.1:$port/twebgpu.html"
