## The avalanches' waves (`cascade.driveWaves`) on an NVIDIA GPU: each GPU
## thread follows electrons of the wave one after the other (`cascade.chase`,
## compiled to CUDA by Crucible as `gpuwalk.nim`'s walkers), with two chases
## of its own (counted, burn-in); the electrons they create go into the next
## wave through an atomic counter, the collision counts and the spectrum
## into 32-bit counters, fetched after each launch. Double precision only.

import std/[monotimes, times, strutils]
import workspace/crucible
from workspace/crucible/src/abis/nvidia_abi import CUdeviceptr
import units, rm48, xoshiro, tables, vec3, fields, freeflight, walkers, drift, cascade
from gpuwalk import releaseEngine

template TabsOf(T: typedesc): typedesc =
  when T is ColdBound: MixtureTable else: GasTables

const
  kindSlots = 30  ## process kinds (IARRY) per gas, then its processes

template slotsPerGas(T: typedesc): int =
  kindSlots + (when T is ColdBound: maxMixProcesses else: maxGasProcesses)

template gasesOf(T: typedesc): int =
  when T is ColdBound: 1 else: maxGases

type
  GpuSink = object
    ## the next wave and the counters, in GPU memory; `counting` 0: burn-in
    next: ptr UncheckedArray[Pending]
    nextCount: ptr UncheckedArray[uint32]
    cap: int32
    times, spectrum, counts: ptr UncheckedArray[uint32]
    stride, counting: int32

proc push(s: GpuSink, p: Pending): int {.inline.} =
  ## index of the new electron, -1 if there is no room
  let i = int32(atomicFetchAdd(addr s.nextCount[0], 1'u32))
  result = -1
  if i < s.cap:
    s.next[i] = p
    result = int(i)

proc isFull(s: GpuSink): bool {.inline.} = int32(s.nextCount[0]) >= s.cap
proc setDir(s: GpuSink, i: int, dir: Vec3) {.inline.} = s.next[i].dir = dir

proc countFlight(s: GpuSink, it: int) {.inline.} =
  if s.counting != 0: discard atomicFetchAdd(addr s.times[it - 1], 1'u32)

proc countBin(s: GpuSink, ie: int) {.inline.} =
  if s.counting != 0: discard atomicFetchAdd(addr s.spectrum[ie - 1], 1'u32)

proc countProcess(s: GpuSink, gas, kind, i: int) {.inline.} =
  if s.counting != 0:
    let base = (gas - 1) * int(s.stride)
    discard atomicFetchAdd(addr s.counts[base + kind - 1], 1'u32)
    discard atomicFetchAdd(addr s.counts[base + kindSlots + i - 1], 1'u32)

proc followOne[G: static FieldGeometry; P: static PlaneKind; T, R](
    ch: var Chase[G, P, T], p: Pending, tabs: TabsOf(T), fm: ptr FlightMaxima,
    sink: GpuSink): Outcome {.inline.} =
  var rng = initRng(R, p.seed)
  let id0 = ch.id
  if startChase(ch, p):
    chase[G, P, T, R, GpuSink](ch, tabs, fm, gpuWindows, rng, sink)
  outcomeOf(ch, p, id0)

template kernel(G: static FieldGeometry, P: static PlaneKind, T, R: typedesc): string =
  cuda:
    proc waveKernel(wave: ptr UncheckedArray[Pending], n, nthreads: int32,
                    init: ptr Chase[G, P, T], chases: ptr UncheckedArray[Chase[G, P, T]],
                    outcomes: ptr UncheckedArray[Outcome],
                    next: ptr UncheckedArray[Pending], nextCount: ptr UncheckedArray[uint32],
                    cap: int32, tabs: ptr TabsOf(T), fm: ptr FlightMaxima,
                    times, spectrum, counts: ptr UncheckedArray[uint32],
                    stride, counting: int32) {.global.} =
      ## the electrons of one kind (`counting` 0: the chains' burn-in, not
      ## counted), one chase per thread, kept from wave to wave
      let tid = int32(blockIdx.x * blockDim.x + threadIdx.x)
      if tid < nthreads:
        var ch = chases[tid]
        let sink = GpuSink(next: next, nextCount: nextCount, cap: cap, times: times,
                           spectrum: spectrum, counts: counts, stride: stride,
                           counting: counting)
        var k = tid
        while k < n:
          outcomes[k] = followOne[G, P, T, R](ch, wave[k], tabs[], fm, sink)
          k = k + nthreads
        chases[tid] = ch

const
  timeKernels = [
    [kernel(fgElectric, pkTime, ColdBound, Xoshiro128StarStar),
     kernel(fgParallel, pkTime, ColdBound, Xoshiro128StarStar),
     kernel(fgCrossed, pkTime, ColdBound, Xoshiro128StarStar),
     kernel(fgOblique, pkTime, ColdBound, Xoshiro128StarStar)],
    [kernel(fgElectric, pkTime, WarmGas, Xoshiro128StarStar),
     kernel(fgParallel, pkTime, WarmGas, Xoshiro128StarStar),
     kernel(fgCrossed, pkTime, WarmGas, Xoshiro128StarStar),
     kernel(fgOblique, pkTime, WarmGas, Xoshiro128StarStar)]]
  spaceKernels = [kernel(fgElectric, pkSpace, ColdBound, Xoshiro128StarStar),
                  kernel(fgElectric, pkSpace, WarmGas, Xoshiro128StarStar)]

func waveKernelSource*(G: static FieldGeometry, P: static PlaneKind, T: typedesc): string =
  ## the CUDA source of the waves' kernel
  when P == pkTime: timeKernels[ord(T is WarmGas)][ord(G)]
  else: spaceKernels[ord(T is WarmGas)]

type
  GpuWaves*[G: static FieldGeometry; P: static PlaneKind; T] = object
    ## waves on the GPU; `c`: the run (and the chases' start)
    engine: CudaEngine
    c: Cascade[G, P, T]
    tab, fm, init: DeviceMem[CUdeviceptr, byte]
    waveBuf, nextBuf, outBuf, chaseBuf: DeviceMem[CUdeviceptr, byte]
    burnBuf: DeviceMem[CUdeviceptr, byte]   ## the burn-ins' chases (dropped)
    burnReady: bool
    waveCap, nextCap: int
    nextCount, times, spectrum, counts: DeviceMem[CUdeviceptr, uint32]
    id: int                ## counted collisions (from the outcomes)
    t0: MonoTime
    kernelTime: Duration
    launches, electrons: int
    work, slots: int       ## debugGpuWaves: collisions, and threads × the longest
    unfetched: int         ## collisions since the counters were fetched

const
  threadsPerBlock = 128
  maxThreads* = 32768  ## GPU threads per wave, at most

proc upload(e: CudaEngine, p: pointer, bytes: int): DeviceMem[CUdeviceptr, byte] =
  result = e.allocDevice(byte, max(bytes, 1))
  if bytes > 0: e.copyToDeviceFromHost(result, cast[ptr byte](p), bytes)

proc zero(e: CudaEngine, m: DeviceMem[CUdeviceptr, uint32], n: int) =
  var z = newSeq[uint32](max(n, 1))
  e.copyToDeviceFromHost(m, addr z[0], z.len)

proc zeroed(e: CudaEngine, n: int): DeviceMem[CUdeviceptr, uint32] =
  result = e.allocDevice(uint32, max(n, 1))
  e.zero(result, n)

proc initGpuWaves*[G: static FieldGeometry; P: static PlaneKind; T](
    c: sink Cascade[G, P, T], tabs: TabsOf(T), fm: FlightMaxima): GpuWaves[G, P, T] =
  ## for the run c (its electrons on the GPU with the exact bounds `fm`)
  releaseEngine()  # the drift's kernel: room on the GPU
  let t0 = getMonoTime()
  result.engine = newEngine(runCuda)
  result.engine.ingest(waveKernelSource(G, P, T))
  when defined(debugGpuWaves):
    stderr.writeLine "gpuwaves: NVRTC ", (getMonoTime() - t0).inMilliseconds, " ms"
  result.t0 = getMonoTime()
  result.c = c
  let ch = initChase(result.c)
  result.tab = result.engine.upload(unsafeAddr tabs, sizeof(tabs))
  result.fm = result.engine.upload(unsafeAddr fm, sizeof(fm))
  result.init = result.engine.upload(unsafeAddr ch, sizeof(ch))
  result.nextCount = result.engine.zeroed(1)
  result.times = result.engine.zeroed(nTimeBins)
  result.spectrum = result.engine.zeroed(nEnergyBins)
  result.counts = result.engine.zeroed(gasesOf(T) * slotsPerGas(T))
  # the threads' chases, kept on the GPU from wave to wave
  var chases = newSeq[Chase[G, P, T]](maxThreads)
  for x in chases.mitems: x = ch
  result.chaseBuf = result.engine.upload(addr chases[0], chases.len * sizeof(ch))

proc gpuCollisions[G: static FieldGeometry; P: static PlaneKind; T](g: GpuWaves[G, P, T]): int =
  g.id
proc gpuStopped[G: static FieldGeometry; P: static PlaneKind; T](g: GpuWaves[G, P, T]): bool =
  g.c.stop != stNone

proc fetch[X](e: CudaEngine, m: DeviceMem[CUdeviceptr, X], n: int): seq[X] =
  result = newSeq[X](n)
  if n > 0: e.copyToHostFromDevice(addr result[0], m, n)

proc fetchCounters[G: static FieldGeometry; P: static PlaneKind; T](g: var GpuWaves[G, P, T]) =
  ## the 32-bit counters into the run, then zeroed (before they could
  ## overflow, and at the end)
  let e = g.engine
  let times = e.fetch(g.times, nTimeBins)
  let spectrum = e.fetch(g.spectrum, nEnergyBins)
  let counts = e.fetch(g.counts, gasesOf(T) * slotsPerGas(T))
  for k in 1 .. nTimeBins: g.c.times[k] = g.c.times[k] + float64(times[k - 1])
  for k in 1 .. nEnergyBins: g.c.spectrum[k] = g.c.spectrum[k] + float64(spectrum[k - 1])
  for gas in 1 .. gasesOf(T):
    let base = (gas - 1) * slotsPerGas(T)
    template tally: untyped =
      when T is ColdBound: g.c.tally else: g.c.tally[gas]
    for k in 1 .. kindSlots:
      tally.kinds[k] = tally.kinds[k] + int(counts[base + k - 1])
    for i in 1 .. slotsPerGas(T) - kindSlots:
      tally.processes[i] = tally.processes[i] + int(counts[base + kindSlots + i - 1])
  e.zero(g.times, nTimeBins)
  e.zero(g.spectrum, nEnergyBins)
  e.zero(g.counts, gasesOf(T) * slotsPerGas(T))
  g.unfetched = 0

proc launch[G: static FieldGeometry; P: static PlaneKind; T](
    g: var GpuWaves[G, P, T], part: seq[Pending], counting: bool,
    next: var seq[Pending], outcomes: var seq[Outcome]) =
  ## one kind of electrons: those they create appended to `next`, their
  ## outcomes in part's order
  let n = part.len
  if n == 0: return
  let e = g.engine
  if n > g.waveCap:
    g.waveCap = max(2 * n, 1024)
    g.waveBuf = e.allocDevice(byte, g.waveCap * sizeof(Pending))
    g.outBuf = e.allocDevice(byte, g.waveCap * sizeof(Outcome))
  e.copyToDeviceFromHost(g.waveBuf, cast[ptr byte](unsafeAddr part[0]), n * sizeof(Pending))
  let nthreads = min(n, maxThreads)
  if g.nextCap < 4 * n:
    g.nextCap = max(4 * n, 65536)
    g.nextBuf = e.allocDevice(byte, g.nextCap * sizeof(Pending))
  if counting and g.unfetched > 1_000_000_000: g.fetchCounters()
  if not counting and not g.burnReady:
    # (only for small avalanches: room on the GPU otherwise)
    let ch = initChase(g.c)
    var chases = newSeq[Chase[G, P, T]](maxThreads)
    for x in chases.mitems: x = ch
    g.burnBuf = e.upload(addr chases[0], chases.len * sizeof(ch))
    g.burnReady = true
  var made: uint32
  while true:
    var zero = 0'u32
    e.copyToDeviceFromHost(g.nextCount, addr zero, 1)
    let blocks = (nthreads + threadsPerBlock - 1) div threadsPerBlock
    let tk = getMonoTime()
    inc g.launches
    g.electrons = g.electrons + n
    defer: g.kernelTime = g.kernelTime + (getMonoTime() - tk)
    if counting:
      e.launchKernel<<(blocks, threadsPerBlock)>>("waveKernel", g.waveBuf,
        (int32(n), int32(nthreads), g.init, g.chaseBuf, g.outBuf, g.nextBuf, g.nextCount,
         int32(g.nextCap), g.tab, g.fm, g.times, g.spectrum, g.counts, int32(slotsPerGas(T)),
         1'i32))
    else:
      e.launchKernel<<(blocks, threadsPerBlock)>>("waveKernel", g.waveBuf,
        (int32(n), int32(nthreads), g.init, g.burnBuf, g.outBuf, g.nextBuf, g.nextCount,
         int32(g.nextCap), g.tab, g.fm, g.times, g.spectrum, g.counts, int32(slotsPerGas(T)),
         0'i32))
    e.copyToHostFromDevice(addr made, g.nextCount, 1)
    if int(made) <= g.nextCap: break
    # no room for all the electrons created
    g.c.stop = stStackFull
    made = uint32(g.nextCap)
    break
  let first = next.len
  next.setLen first + int(made)
  if made > 0:
    e.copyToHostFromDevice(cast[ptr byte](addr next[first]), g.nextBuf,
                           int(made) * sizeof(Pending))
  outcomes.setLen n
  e.copyToHostFromDevice(cast[ptr byte](addr outcomes[0]), g.outBuf, n * sizeof(Outcome))
  if counting:
    for o in outcomes:
      g.id = g.id + int(o.collisions)
      g.unfetched = g.unfetched + int(o.collisions)
  when defined(debugGpuWaves):
    var mx = 0
    var thr = newSeq[int](nthreads)
    for k, o in outcomes:
      thr[k mod nthreads] += int(o.collisions)
      g.work = g.work + int(o.collisions)
    for x in thr: mx = max(mx, x)
    g.slots = g.slots + mx * nthreads

proc gpuProcess[G: static FieldGeometry; P: static PlaneKind; T](
    g: var GpuWaves[G, P, T], wave, next: var seq[Pending], outcomes: var seq[Outcome]) =
  ## follows the wave's electrons on the GPU (see `CpuWaves.process`):
  ## the counted ones, then those of burn-ins
  var counted, burning: seq[Pending]
  var at: seq[(bool, int)]  # where each electron of the wave went
  for p in wave:
    if p.burn != 0:
      at.add (true, burning.len)
      burning.add p
    else:
      at.add (false, counted.len)
      counted.add p
  next.setLen 0
  var oc, ob: seq[Outcome]
  g.launch(counted, true, next, oc)
  g.launch(burning, false, next, ob)
  outcomes.setLen wave.len
  for k, (burn, i) in at: outcomes[k] = (if burn: ob[i] else: oc[i])

proc finish*[G: static FieldGeometry; P: static PlaneKind; T](
    g: var GpuWaves[G, P, T]): Cascade[G, P, T] =
  ## the run: the threads' chases and the counters added up
  when defined(debugGpuWaves):
    stderr.writeLine "gpuwaves: ", g.launches, " launches, ", g.electrons, " electrons, kernels ",
                     g.kernelTime.inMilliseconds, " ms of ", (getMonoTime() - g.t0).inMilliseconds,
                     " ms, collisions ", g.id, ", busy ",
                     formatFloat(100 * g.work / max(g.slots, 1), ffDecimal, 1), " %"
  g.fetchCounters()
  var chases = newSeq[Chase[G, P, T]](maxThreads)
  g.engine.copyToHostFromDevice(cast[ptr byte](addr chases[0]), g.chaseBuf,
                                chases.len * sizeof(Chase[G, P, T]))
  for ch in chases: g.c.absorb(ch)
  move g.c

proc runGpuWaves*[G: static FieldGeometry; P: static PlaneKind; T](
    c: sink Cascade[G, P, T], tabs: TabsOf(T), fm: FlightMaxima, start: Electron,
    seed: var MagRng): Cascade[G, P, T] =
  ## `cascade.driveWaves` on the GPU, for the run c
  let nmax = c.nmax
  var b = initGpuWaves[G, P, T](c, tabs, fm)
  let (_, stop) = driveWaves(b, nmax, start, seed, gpuProcess[G, P, T], gpuCollisions[G, P, T],
                             gpuStopped[G, P, T], maxThreads)
  result = b.finish()
  if stop == stTooManyPrimaries: result.stop = stop
  result.nclus = result.nclus + result.iprim - 1  # primaries after the first, as `run`
  result.j1 = 1 + result.id div max(nmax div 10, 1)
