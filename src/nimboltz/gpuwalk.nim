## Walkers (`walkers.nim`) on an NVIDIA GPU, one per thread.
##
## The kernel is NimBoltz's own `walk`, compiled to CUDA at Nim compile
## time by Crucible, the Nim-to-GPU compiler of tattletale
## (github.com/mratsim/tattletale); NVRTC turns it into PTX at run time.
## Build with `--path:` to the tattletale root.

import workspace/crucible
from workspace/crucible/src/abis/nvidia_abi import CUdeviceptr
import units, tables, fields, freeflight, walkers

type
  GpuTally* = object
    ## real collisions per process and energy bin, in GPU memory shared by
    ## all walkers (atomic adds)
    processes*, spectrum*: ptr UncheckedArray[float64]

proc count*(t: GpuTally, i, ie: int) {.inline.} =
  discard atomicFetchAdd(addr t.processes[i - 1], 1.0)
  discard atomicFetchAdd(addr t.spectrum[ie - 1], 1.0)

# One kernel per field geometry and target (Crucible compiles one kernel
# per block).
template kernel(G: static FieldGeometry, T: typedesc): string =
  cuda:
    proc walkKernel(ws: ptr UncheckedArray[Walker[T]], n: int32, f: Field,
                    tab: ptr TablesOf(T), fm: ptr FlightMaxima, estep: eV,
                    collisions, penning, counting: int32,
                    processes, spectrum: ptr UncheckedArray[float64]) {.global.} =
      let k = int32(blockIdx.x * blockDim.x + threadIdx.x)
      if k < n:
        if counting != 0:
          walk[G, T, GpuTally](ws[k], f, tab[], fm, estep, collisions, penning,
                               GpuTally(processes: processes, spectrum: spectrum))
        else:
          walk[G, T, NoTally](ws[k], f, tab[], fm, estep, collisions, penning, NoTally())

const
  coldElectric = kernel(fgElectric, ColdGas)
  coldParallel = kernel(fgParallel, ColdGas)
  coldCrossed = kernel(fgCrossed, ColdGas)
  coldOblique = kernel(fgOblique, ColdGas)
  warmElectric = kernel(fgElectric, WarmGas)
  warmParallel = kernel(fgParallel, WarmGas)
  warmCrossed = kernel(fgCrossed, WarmGas)
  warmOblique = kernel(fgOblique, WarmGas)

func kernelSource*(G: static FieldGeometry, T: typedesc): string =
  ## the CUDA source of the kernel for geometry G and target T
  when T is ColdGas:
    when G == fgElectric: coldElectric
    elif G == fgParallel: coldParallel
    elif G == fgCrossed: coldCrossed
    else: coldOblique
  else:
    when G == fgElectric: warmElectric
    elif G == fgParallel: warmParallel
    elif G == fgCrossed: warmCrossed
    else: warmOblique

type
  GpuWalkers*[G: static FieldGeometry; T] = object
    ## walkers, their collision tables and the compiled kernel on the GPU
    engine: CudaEngine
    ws: DeviceMem[CUdeviceptr, Walker[T]]
    tab: DeviceMem[CUdeviceptr, TablesOf(T)]
    fm: DeviceMem[CUdeviceptr, FlightMaxima]
    n: int
    f: Field
    estep: eV
    penning: int32
    processes, spectrum: DeviceMem[CUdeviceptr, float64]

const threadsPerBlock = 128

var
  lastEngine: CudaEngine
  lastKernel: (FieldGeometry, bool, bool)
    ## the last compiled kernel (geometry, thermal motion, exact), kept for
    ## the next run (NVRTC takes a while); only one, as each engine has a
    ## CUDA context of its own, of some 100 MB

proc engineFor(G: static FieldGeometry, T: typedesc, exact: bool): CudaEngine =
  let kernel = (G, T is WarmGas, exact)
  if lastEngine == nil or lastKernel != kernel:
    lastEngine = nil  # the old context goes once its walkers are gone
    lastEngine = newEngine(runCuda)
    if exact: lastEngine.ingest(kernelSource(G, T), ["--fmad=false"])
    else: lastEngine.ingest(kernelSource(G, T))
    lastKernel = kernel
  lastEngine

proc clearTally*(g: var GpuWalkers) =
  ## zeroes the collision counts
  var t: Tallies
  g.engine.copyToDeviceFromHost(g.processes, addr t.processes[1], tallyProcesses)
  g.engine.copyToDeviceFromHost(g.spectrum, addr t.spectrum[1], nEnergyBins)

proc readTally*(g: GpuWalkers, t: var Tallies) =
  ## the collision counts so far
  g.engine.copyToHostFromDevice(addr t.processes[1], g.processes, tallyProcesses)
  g.engine.copyToHostFromDevice(addr t.spectrum[1], g.spectrum, nEnergyBins)

proc initGpuWalkers*[G: static FieldGeometry; T](ws: openArray[Walker[T]], f: Field,
                                                 tab: TablesOf(T), fm: FlightMaxima,
                                                 estep: eV, penning = false,
                                                 exact = false): GpuWalkers[G, T] =
  ## Copies the walkers, the table and its range maxima (`flightMaxima`,
  ## for the exact bounds at B = 0 and B ∥ E) to the GPU. `exact`: no fused
  ## multiply-adds, so the arithmetic is that of the CPU (but `log`, `sin`,
  ## ... of CUDA, which may differ from the CPU's in the last bit).
  result.engine = engineFor(G, T, exact)
  result.penning = int32(ord(penning))
  result.n = ws.len
  result.f = f
  result.estep = estep
  result.ws = result.engine.allocDevice(Walker[T], ws.len)
  result.engine.copyToDeviceFromHost(result.ws, unsafeAddr ws[0], ws.len)
  result.tab = result.engine.allocDevice(TablesOf(T), 1)
  result.engine.copyToDeviceFromHost(result.tab, unsafeAddr tab, 1)
  result.fm = result.engine.allocDevice(FlightMaxima, 1)
  result.engine.copyToDeviceFromHost(result.fm, unsafeAddr fm, 1)
  result.processes = result.engine.allocDevice(float64, tallyProcesses)
  result.spectrum = result.engine.allocDevice(float64, nEnergyBins)
  result.clearTally()

proc walk*(g: var GpuWalkers, collisions: int32, counting = true) =
  ## `collisions` real collisions of every walker, counted in the tally
  ## unless not `counting`
  let blocks = (g.n + threadsPerBlock - 1) div threadsPerBlock
  g.engine.launchKernel<<(blocks, threadsPerBlock)>>("walkKernel", g.ws,
    (int32(g.n), g.f, g.tab, g.fm, g.estep, collisions, g.penning, int32(ord(counting)),
     g.processes, g.spectrum))

func dropsWalkers*(g: GpuWalkers): bool = false
  ## time goes with the collisions per walker: fewer walkers, no less time

proc shrink*(g: var GpuWalkers, n: int) =
  ## keeps the first n walkers
  doAssert n <= g.n
  g.n = n

proc read*[T](g: GpuWalkers, ws: var openArray[Walker[T]]) =
  ## the walkers' state, back from the GPU
  doAssert ws.len == g.n
  g.engine.copyToHostFromDevice(addr ws[0], g.ws, g.n)

proc write*[T](g: var GpuWalkers, ws: openArray[Walker[T]]) =
  ## replaces the walkers' state on the GPU
  doAssert ws.len == g.n
  g.engine.copyToDeviceFromHost(g.ws, unsafeAddr ws[0], g.n)
