// Runs NimBoltz (nimboltz.js/.wasm) off the page's thread. Messages:
//   in:  {setup}                        one run (see src/nimboltz/wasm.nim)
//   out: {type: 'ready', gases, scale}  once the module is loaded
//        {type: 'line', text}           a line of Magboltz's printout
//        {type: 'result', result, ms}   the run's results (or {error})
// A run can't be interrupted from inside: the page terminates the worker.
importScripts('nimboltz.js');

// A WebGPU device for the walkers, if the browser has one (null if not).
// The kernel needs WGSL's pointer parameters (Chrome/Edge 123 and later).
async function openGpu(m) {
  if (!navigator.gpu) return null;
  const f = navigator.gpu.wgslLanguageFeatures;
  if (!f || !f.has('unrestricted_pointer_parameters') || !f.has('pointer_composite_access')) return null;
  const adapter = await navigator.gpu.requestAdapter({powerPreference: 'high-performance'});
  if (!adapter) return null;
  const l = adapter.limits;
  const device = await adapter.requestDevice({requiredLimits: {
    maxStorageBufferBindingSize: l.maxStorageBufferBindingSize, maxBufferSize: l.maxBufferSize}});
  device.addEventListener('uncapturederror', (e) => {
    m.gpuError = e.error.message;  // the run stops (webgpu.js)
    postMessage({type: 'line', text: 'WebGPU: ' + e.error.message});
  });
  device.lost.then((i) => { m.gpuError = 'device lost: ' + i.message; });
  const i = adapter.info || {};
  device.label = [i.vendor, i.architecture, i.device || i.description].filter(Boolean).join(' ') || 'WebGPU';
  // A software adapter (SwiftShader) runs the kernels on the CPU: the
  // browser falls back to it when it can't use the real GPU.
  device.software = !!(i.isFallbackAdapter || adapter.isFallbackAdapter ||
                       /swiftshader|llvmpipe|lavapipe/i.test(device.label) || /^0xc0de$/i.test(i.device || ''));
  return device;
}

const ready = createNimBoltz({
  print: (text) => postMessage({type: 'line', text}),
  printErr: (text) => postMessage({type: 'line', text}),
}).then(async (m) => {
  const gases = JSON.parse(m.UTF8ToString(m.ccall('nimboltzGases', 'number', [], [])));
  m.webgpu = m.ccall('nimboltzWebGpu', 'number', [], []) !== 0;  // built with it
  if (m.webgpu) {
    try { m.gpuDevice = await openGpu(m); } catch (err) { m.gpuDevice = null; }
  }
  postMessage({type: 'ready', gases, scale: m.ccall('nimboltzCollisionScale', 'number', [], []),
               gpu: m.gpuDevice ? m.gpuDevice.label : null,
               software: m.gpuDevice ? m.gpuDevice.software : false});
  return m;
});

onmessage = async (e) => {
  const m = await ready;
  const t = performance.now();
  let result;
  try {
    // asynchronous with WebGPU: the module waits for the GPU (Asyncify)
    result = JSON.parse(await m.ccall('nimboltzRun', 'string', ['string'], [JSON.stringify(e.data.setup)],
                                      m.webgpu ? {async: true} : undefined));
  } catch (err) {  // a Nim defect or a trap: the module is unusable now
    result = {error: String(err), fatal: true};
  }
  postMessage({type: 'result', result, ms: performance.now() - t});
};
