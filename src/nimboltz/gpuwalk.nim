## Walkers (`walkers.nim`) on an NVIDIA GPU, one per thread.
##
## The kernel is NimBoltz's own `walk`, compiled to CUDA at Nim compile
## time by Crucible, the Nim-to-GPU compiler of tattletale
## (github.com/mratsim/tattletale); NVRTC turns it into PTX at run time.
## Build with `--path:` to the tattletale root.

import workspace/crucible
from workspace/crucible/src/abis/nvidia_abi import CUdeviceptr
import units, rm48, xoshiro, tables, fields, freeflight, walkers, single, soa

type
  GpuTally* = object
    ## real collisions per process and energy bin, in GPU memory shared by
    ## all walkers (atomic adds; counted on the host after each launch, so
    ## 32 bits are plenty)
    processes*, spectrum*: ptr UncheckedArray[uint32]

proc count*(t: GpuTally, i, ie: int) {.inline.} =
  discard atomicFetchAdd(addr t.processes[i - 1], 1'u32)
  discard atomicFetchAdd(addr t.spectrum[ie - 1], 1'u32)

# One kernel per field geometry, target and generator (Crucible compiles
# one kernel per block), in double precision (`cuda`) or demoted to single
# precision (`cudaSingle`: every float64 a float32, see single.nim). The
# walkers lie in slots (soa.nim): a thread copies its walker into a
# variable of its own (registers, as far as they go), walks, and copies it
# back.
template kernel(G: static FieldGeometry, T, R: typedesc, compiler: untyped): string =
  compiler:
    proc walkKernel(fs: ptr UncheckedArray[float64], us: ptr UncheckedArray[uint32],
                    n, stride: int32, f: Field,
                    tab: ptr TablesOf(T), fm: ptr FlightMaxima, estep: eV,
                    collisions, penning, counting: int32,
                    processes, spectrum: ptr UncheckedArray[uint32]) {.global.} =
      let k = int32(blockIdx.x * blockDim.x + threadIdx.x)
      if k < n:
        var w: Walker[T, R]
        load(w, fs, us, k, stride)
        if counting != 0:
          walkWith[G, T, GpuTally, R](w, f, tab[], fm, gpuWindows, estep, collisions, penning,
                                      GpuTally(processes: processes, spectrum: spectrum))
        else:
          walkWith[G, T, NoTally, R](w, f, tab[], fm, gpuWindows, estep, collisions, penning,
                                     NoTally())
        store(w, fs, us, k, stride)

template kernels(R: typedesc, compiler: untyped): untyped =
  [[kernel(fgElectric, ColdGas, R, compiler), kernel(fgParallel, ColdGas, R, compiler),
    kernel(fgCrossed, ColdGas, R, compiler), kernel(fgOblique, ColdGas, R, compiler)],
   [kernel(fgElectric, WarmGas, R, compiler), kernel(fgParallel, WarmGas, R, compiler),
    kernel(fgCrossed, WarmGas, R, compiler), kernel(fgOblique, WarmGas, R, compiler)]]

const
  doubleKernels = kernels(Xoshiro128StarStar, cuda)
  singleKernels = kernels(Xoshiro128Plus, cudaSingle)

when defined(gpuRanmar):
  # Ranmar on the GPU too (1 KB of state per walker: slower), for comparisons
  const
    ranmarDoubleKernels = kernels(Ranmar, cuda)
    ranmarSingleKernels = kernels(Ranmar, cudaSingle)

func kernelSource*(G: static FieldGeometry, T, R: typedesc, single = false): string =
  ## the CUDA source of the kernel for geometry G, target T and generator
  ## R: xoshiro128** in double precision, xoshiro128+ in single precision
  ## (23-bit numbers, exact in a float32); Ranmar in both (`toSingle`) if
  ## built with -d:gpuRanmar
  const t = ord(T is WarmGas)
  when R is Xoshiro128StarStar:
    if single: raise newException(ValueError, "walkers on a GPU in single precision: xoshiro128+")
    doubleKernels[t][ord(G)]
  elif R is Xoshiro128Plus:
    if not single: raise newException(ValueError, "walkers on a GPU in double precision: xoshiro128**")
    singleKernels[t][ord(G)]
  elif R is Ranmar and defined(gpuRanmar):
    if single: ranmarSingleKernels[t][ord(G)] else: ranmarDoubleKernels[t][ord(G)]
  else:
    {.error: "walkers on a GPU: xoshiro128** or xoshiro128+ (Ranmar with -d:gpuRanmar)".}

type
  GpuWalkers*[G: static FieldGeometry; T, R] = object
    ## walkers, their collision tables and the compiled kernel on the GPU;
    ## `single`: all of it in single precision
    engine: CudaEngine
    single: bool
    stride: int            ## walkers room is made for (each slot's length)
    fs, us, tab, fm: DeviceMem[CUdeviceptr, byte]
    n: int
    f: Field
    estep: eV
    penning: int32
    processes, spectrum: DeviceMem[CUdeviceptr, uint32]
    tallies: ref Tallies   ## the counts so far, on the host

const threadsPerBlock = 128

var
  lastEngine: CudaEngine
  lastKernel: (FieldGeometry, bool, string, bool, bool)
    ## the last compiled kernel (geometry, thermal motion, generator, exact, single),
    ## kept for the next run (NVRTC takes a while); only one, as each
    ## engine has a CUDA context of its own, of some 100 MB

proc releaseEngine*() =
  ## frees the last compiled kernel's context (and its GPU memory)
  lastEngine = nil

proc engineFor(G: static FieldGeometry, T, R: typedesc, exact, single: bool): CudaEngine =
  let kernel = (G, T is WarmGas, $R, exact, single)
  if lastEngine == nil or lastKernel != kernel:
    let src = kernelSource(G, T, R, single)
    lastEngine = nil  # the old context goes once its walkers are gone
    lastEngine = newEngine(runCuda)
    if exact: lastEngine.ingest(src, ["--fmad=false"])
    else: lastEngine.ingest(src)
    lastKernel = kernel
  lastEngine

proc bytesOf[X](x: openArray[X], inSingle: bool): seq[byte] =
  ## as the GPU has them
  if inSingle: return single(x)
  result = newSeq[byte](x.len * sizeof(X))
  if x.len > 0: copyMem(addr result[0], unsafeAddr x[0], result.len)

template one[X](x: X): untyped =
  ## x as an openArray of one, in place (tables are large for the stack)
  cast[ptr UncheckedArray[X]](unsafeAddr x).toOpenArray(0, 0)

proc upload(e: CudaEngine, bytes: seq[byte]): DeviceMem[CUdeviceptr, byte] =
  result = e.allocDevice(byte, max(bytes.len, 1))
  if bytes.len > 0: e.copyToDeviceFromHost(result, unsafeAddr bytes[0], bytes.len)

proc clearTally*(g: var GpuWalkers) =
  ## zeroes the collision counts
  g.tallies[] = Tallies()

proc readTally*(g: GpuWalkers, t: var Tallies) =
  ## the collision counts so far
  t = g.tallies[]

proc fetchTally(g: var GpuWalkers) =
  ## adds the GPU's counts to the host's, and zeroes them
  var p = newSeq[uint32](tallyProcesses)
  var s = newSeq[uint32](nEnergyBins)
  g.engine.copyToHostFromDevice(addr p[0], g.processes, tallyProcesses)
  g.engine.copyToHostFromDevice(addr s[0], g.spectrum, nEnergyBins)
  for i in 1 .. tallyProcesses:
    g.tallies.processes[i] = g.tallies.processes[i] + float64(p[i - 1])
  for i in 1 .. nEnergyBins:
    g.tallies.spectrum[i] = g.tallies.spectrum[i] + float64(s[i - 1])
  zeroMem(addr p[0], p.len * sizeof(uint32))
  zeroMem(addr s[0], s.len * sizeof(uint32))
  g.engine.copyToDeviceFromHost(g.processes, addr p[0], tallyProcesses)
  g.engine.copyToDeviceFromHost(g.spectrum, addr s[0], nEnergyBins)

proc initGpuWalkers*[G: static FieldGeometry; T, R](ws: openArray[Walker[T, R]], f: Field,
                                                    tab: TablesOf(T), fm: FlightMaxima,
                                                    estep: eV, penning = false,
                                                    exact = false, single = false): GpuWalkers[G, T, R] =
  ## Copies the walkers, the table and its range maxima (`flightMaxima`,
  ## for the exact bounds at B = 0 and B ∥ E) to the GPU. `exact`: no fused
  ## multiply-adds, so the arithmetic is that of the CPU (but `log`, `sin`,
  ## ... of CUDA, which may differ from the CPU's in the last bit).
  ## `single`: the kernel and all data in single precision (and the
  ## walkers' random numbers xoshiro128+, or the 24-bit RANMAR:
  ## `rm48.toSingle`).
  result.engine = engineFor(G, T, R, exact, single)
  result.single = single
  result.stride = ws.len
  result.penning = int32(ord(penning))
  result.n = ws.len
  result.f = f
  result.estep = estep
  let (fb, ub) = slotBytes(ws, result.stride, single)
  result.fs = result.engine.upload(fb)
  result.us = result.engine.upload(ub)
  result.tab = result.engine.upload(bytesOf(one(tab), single))
  result.fm = result.engine.upload(bytesOf(one(fm), single))
  result.processes = result.engine.allocDevice(uint32, tallyProcesses)
  result.spectrum = result.engine.allocDevice(uint32, nEnergyBins)
  new result.tallies
  var p = newSeq[uint32](tallyProcesses)
  var s = newSeq[uint32](nEnergyBins)
  result.engine.copyToDeviceFromHost(result.processes, addr p[0], tallyProcesses)
  result.engine.copyToDeviceFromHost(result.spectrum, addr s[0], nEnergyBins)

func singleField(f: Field): array[sizeof(Field) div 8, float32] =
  ## the field in single precision (all its members are floats)
  var k = 0
  for _, x in f.fieldPairs:
    result[k] = float32(x.float)
    k = k + 1
  doAssert k == result.len

proc walk*(g: var GpuWalkers, collisions: int32, counting = true) =
  ## `collisions` real collisions of every walker, counted in the tally
  ## unless not `counting`
  let blocks = (g.n + threadsPerBlock - 1) div threadsPerBlock
  if g.single:
    g.engine.launchKernel<<(blocks, threadsPerBlock)>>("walkKernel", g.fs,
      (g.us, int32(g.n), int32(g.stride), singleField(g.f), g.tab, g.fm, float32(g.estep.float),
       collisions, g.penning, int32(ord(counting)), g.processes, g.spectrum))
  else:
    g.engine.launchKernel<<(blocks, threadsPerBlock)>>("walkKernel", g.fs,
      (g.us, int32(g.n), int32(g.stride), g.f, g.tab, g.fm, g.estep, collisions, g.penning,
       int32(ord(counting)), g.processes, g.spectrum))
  if counting: g.fetchTally()

func dropsWalkers*(g: GpuWalkers): bool = false
  ## time goes with the collisions per walker: fewer walkers, no less time

proc shrink*(g: var GpuWalkers, n: int) =
  ## keeps the first n walkers
  doAssert n <= g.n
  g.n = n

proc read*[T, R](g: GpuWalkers, ws: var openArray[Walker[T, R]]) =
  ## the walkers' state, back from the GPU
  doAssert ws.len == g.n
  const c = slotCounts(Walker[T, R])
  let size = if g.single: 4 else: 8
  var fb = newSeq[byte](max(c[0] * g.stride * size, 1))
  var ub = newSeq[byte](max(c[1] * g.stride * 4, 1))
  g.engine.copyToHostFromDevice(addr fb[0], g.fs, fb.len)
  g.engine.copyToHostFromDevice(addr ub[0], g.us, ub.len)
  fromSlotBytes(ws, fb, ub, g.stride, g.single)

proc write*[T, R](g: var GpuWalkers, ws: openArray[Walker[T, R]]) =
  ## replaces the walkers' state on the GPU
  doAssert ws.len == g.n
  let (fb, ub) = slotBytes(ws, g.stride, g.single)
  g.engine.copyToDeviceFromHost(g.fs, unsafeAddr fb[0], fb.len)
  g.engine.copyToDeviceFromHost(g.us, unsafeAddr ub[0], ub.len)
