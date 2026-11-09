// WebGPU for the walkers of src/nimboltz/webgpuwalk.nim (an Emscripten
// JS library: build.sh links it with --js-library). The worker opens the
// device before a run (`Module.gpuDevice`); buffers and pipelines live in
// tables here, Nim holds their numbers. Compiling and reading back are
// asynchronous, the module waits for them with Asyncify.
addToLibrary({
  $nbgpu: {buffers: [null], pipelines: [null], bySource: new Map()},

  // 1 if the worker opened a device; clears the last error (Module.gpuError,
  // set by the device's uncapturederror handler or here)
  nbgpuReady__deps: ['$nbgpu'],
  nbgpuReady: () => { Module.gpuError = null; return Module.gpuDevice ? 1 : 0; },

  // the compute pipeline of a WGSL source (one per source), 0 if it fails
  nbgpuPipeline__deps: ['$nbgpu'],
  nbgpuPipeline__async: 'auto',
  nbgpuPipeline: async (src, len) => {
    const code = UTF8ToString(src, len);
    let id = nbgpu.bySource.get(code);
    if (id === undefined) {
      const dev = Module.gpuDevice;
      try {
        const module = dev.createShaderModule({code});
        const pipeline = await dev.createComputePipelineAsync({layout: 'auto', compute: {module, entryPoint: 'walkKernel'}});
        id = nbgpu.pipelines.push(pipeline) - 1;
        nbgpu.bySource.set(code, id);
      } catch (e) {
        Module.gpuError = String(e.message || e);
        err('WebGPU: ' + Module.gpuError);
        return 0;
      }
    }
    return id;
  },

  nbgpuBuffer__deps: ['$nbgpu'],
  nbgpuBuffer: (bytes) => {
    const size = Math.ceil(bytes / 4) * 4;
    const buffer = Module.gpuDevice.createBuffer({size,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC});
    const free = nbgpu.buffers.indexOf(null, 1);
    if (free > 0) { nbgpu.buffers[free] = buffer; return free; }
    return nbgpu.buffers.push(buffer) - 1;
  },

  nbgpuWrite__deps: ['$nbgpu'],
  nbgpuWrite: (id, offset, ptr, len) => {
    Module.gpuDevice.queue.writeBuffer(nbgpu.buffers[id], offset, HEAPU8, ptr, Math.ceil(len / 4) * 4);
  },

  nbgpuRun__deps: ['$nbgpu'],
  nbgpuRun: (pipe, ws, tab, fm, args, tally, groups) => {
    const dev = Module.gpuDevice;
    const pipeline = nbgpu.pipelines[pipe];
    const entries = [ws, tab, fm, args, tally].map((b, binding) => ({binding, resource: {buffer: nbgpu.buffers[b]}}));
    const bindGroup = dev.createBindGroup({layout: pipeline.getBindGroupLayout(0), entries});
    const enc = dev.createCommandEncoder();
    const pass = enc.beginComputePass();
    pass.setPipeline(pipeline);
    pass.setBindGroup(0, bindGroup);
    pass.dispatchWorkgroups(groups);
    pass.end();
    dev.queue.submit([enc.finish()]);
  },

  // copies a buffer to the module's memory; 0 if the device reported an error
  nbgpuRead__deps: ['$nbgpu'],
  nbgpuRead__async: 'auto',
  nbgpuRead: async (id, ptr, len) => {
    const dev = Module.gpuDevice;
    const size = Math.ceil(len / 4) * 4;
    const staging = dev.createBuffer({size, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST});
    const enc = dev.createCommandEncoder();
    enc.copyBufferToBuffer(nbgpu.buffers[id], 0, staging, 0, size);
    dev.queue.submit([enc.finish()]);
    await staging.mapAsync(GPUMapMode.READ);
    HEAPU8.set(new Uint8Array(staging.getMappedRange(), 0, len), ptr);
    staging.unmap();
    staging.destroy();
    return Module.gpuError ? 0 : 1;
  },

  nbgpuFree__deps: ['$nbgpu'],
  nbgpuFree: (id) => {
    if (nbgpu.buffers[id]) nbgpu.buffers[id].destroy();
    nbgpu.buffers[id] = null;
  },
});
