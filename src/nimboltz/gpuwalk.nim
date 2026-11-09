## Walkers (`walkers.nim`) on an NVIDIA GPU, one per thread.
##
## The kernel is NimBoltz's own `walk`, compiled to CUDA at Nim compile
## time by Crucible, the Nim-to-GPU compiler of tattletale
## (github.com/mratsim/tattletale); NVRTC turns it into PTX at run time.
## Build with `--path:` to the tattletale root.

import workspace/crucible
from workspace/crucible/src/abis/nvidia_abi import CUdeviceptr
import units, tables, fields, freeflight, walkers, single

type
  GpuTally* = object
    ## real collisions per process and energy bin, in GPU memory shared by
    ## all walkers (atomic adds; counted on the host after each launch, so
    ## 32 bits are plenty)
    processes*, spectrum*: ptr UncheckedArray[uint32]

proc count*(t: GpuTally, i, ie: int) {.inline.} =
  discard atomicFetchAdd(addr t.processes[i - 1], 1'u32)
  discard atomicFetchAdd(addr t.spectrum[ie - 1], 1'u32)

# One kernel per field geometry and target (Crucible compiles one kernel
# per block), in double precision (`cuda`) or demoted to single precision
# (`cudaSingle`: every float64 a float32, see single.nim).
template kernel(G: static FieldGeometry, T: typedesc, compiler: untyped): string =
  compiler:
    proc walkKernel(ws: ptr UncheckedArray[Walker[T]], n: int32, f: Field,
                    tab: ptr TablesOf(T), fm: ptr FlightMaxima, estep: eV,
                    collisions, penning, counting: int32,
                    processes, spectrum: ptr UncheckedArray[uint32]) {.global.} =
      let k = int32(blockIdx.x * blockDim.x + threadIdx.x)
      if k < n:
        if counting != 0:
          walkWith[G, T, GpuTally](ws[k], f, tab[], fm, true, estep, collisions, penning,
                               GpuTally(processes: processes, spectrum: spectrum))
        else:
          walkWith[G, T, NoTally](ws[k], f, tab[], fm, true, estep, collisions, penning,
                                  NoTally())

template kernels(compiler: untyped): untyped =
  [[kernel(fgElectric, ColdGas, compiler), kernel(fgParallel, ColdGas, compiler),
    kernel(fgCrossed, ColdGas, compiler), kernel(fgOblique, ColdGas, compiler)],
   [kernel(fgElectric, WarmGas, compiler), kernel(fgParallel, WarmGas, compiler),
    kernel(fgCrossed, WarmGas, compiler), kernel(fgOblique, WarmGas, compiler)]]

const
  doubleKernels = kernels(cuda)
  singleKernels = kernels(cudaSingle)

func kernelSource*(G: static FieldGeometry, T: typedesc, single = false): string =
  ## the CUDA source of the kernel for geometry G and target T
  let k = if single: singleKernels else: doubleKernels
  k[ord(T is WarmGas)][ord(G)]

type
  GpuWalkers*[G: static FieldGeometry; T] = object
    ## walkers, their collision tables and the compiled kernel on the GPU;
    ## `single`: all of it in single precision
    engine: CudaEngine
    single: bool
    stride: int            ## bytes per walker on the GPU
    ws, tab, fm: DeviceMem[CUdeviceptr, byte]
    n: int
    f: Field
    estep: eV
    penning: int32
    processes, spectrum: DeviceMem[CUdeviceptr, uint32]
    tallies: ref Tallies   ## the counts so far, on the host

const threadsPerBlock = 128

var
  lastEngine: CudaEngine
  lastKernel: (FieldGeometry, bool, bool, bool)
    ## the last compiled kernel (geometry, thermal motion, exact, single),
    ## kept for the next run (NVRTC takes a while); only one, as each
    ## engine has a CUDA context of its own, of some 100 MB

proc engineFor(G: static FieldGeometry, T: typedesc, exact, single: bool): CudaEngine =
  let kernel = (G, T is WarmGas, exact, single)
  if lastEngine == nil or lastKernel != kernel:
    lastEngine = nil  # the old context goes once its walkers are gone
    lastEngine = newEngine(runCuda)
    if exact: lastEngine.ingest(kernelSource(G, T, single), ["--fmad=false"])
    else: lastEngine.ingest(kernelSource(G, T, single))
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

proc initGpuWalkers*[G: static FieldGeometry; T](ws: openArray[Walker[T]], f: Field,
                                                 tab: TablesOf(T), fm: FlightMaxima,
                                                 estep: eV, penning = false,
                                                 exact = false, single = false): GpuWalkers[G, T] =
  ## Copies the walkers, the table and its range maxima (`flightMaxima`,
  ## for the exact bounds at B = 0 and B ∥ E) to the GPU. `exact`: no fused
  ## multiply-adds, so the arithmetic is that of the CPU (but `log`, `sin`,
  ## ... of CUDA, which may differ from the CPU's in the last bit).
  ## `single`: the kernel and all data in single precision (and the
  ## walkers' random numbers the 24-bit RANMAR: `rm48.toSingle`).
  result.engine = engineFor(G, T, exact, single)
  result.single = single
  result.stride = if single: singleSize(Walker[T]) else: sizeof(Walker[T])
  result.penning = int32(ord(penning))
  result.n = ws.len
  result.f = f
  result.estep = estep
  result.ws = result.engine.upload(bytesOf(ws, single))
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
    g.engine.launchKernel<<(blocks, threadsPerBlock)>>("walkKernel", g.ws,
      (int32(g.n), singleField(g.f), g.tab, g.fm, float32(g.estep.float), collisions, g.penning,
       int32(ord(counting)), g.processes, g.spectrum))
  else:
    g.engine.launchKernel<<(blocks, threadsPerBlock)>>("walkKernel", g.ws,
      (int32(g.n), g.f, g.tab, g.fm, g.estep, collisions, g.penning, int32(ord(counting)),
       g.processes, g.spectrum))
  if counting: g.fetchTally()

func dropsWalkers*(g: GpuWalkers): bool = false
  ## time goes with the collisions per walker: fewer walkers, no less time

proc shrink*(g: var GpuWalkers, n: int) =
  ## keeps the first n walkers
  doAssert n <= g.n
  g.n = n

proc read*[T](g: GpuWalkers, ws: var openArray[Walker[T]]) =
  ## the walkers' state, back from the GPU
  doAssert ws.len == g.n
  var bytes = newSeq[byte](g.n * g.stride)
  g.engine.copyToHostFromDevice(addr bytes[0], g.ws, bytes.len)
  if g.single:
    var pos = 0
    for w in ws.mitems: unpack32(w, bytes, pos)
  else:
    copyMem(addr ws[0], addr bytes[0], bytes.len)

proc write*[T](g: var GpuWalkers, ws: openArray[Walker[T]]) =
  ## replaces the walkers' state on the GPU
  doAssert ws.len == g.n
  let bytes = bytesOf(ws, g.single)
  g.engine.copyToDeviceFromHost(g.ws, unsafeAddr bytes[0], bytes.len)
