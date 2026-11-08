## Walkers (`walkers.nim`) on an NVIDIA GPU, one per thread.
##
## The kernel is NimBoltz's own `walk`, compiled to CUDA at Nim compile
## time by Crucible, the Nim-to-GPU compiler of tattletale
## (github.com/mratsim/tattletale); NVRTC turns it into PTX at run time.
## Build with `--path:` to the tattletale root.

import workspace/crucible
from workspace/crucible/src/abis/nvidia_abi import CUdeviceptr
import units, tables, fields, freeflight, walkers

const walkKernel = cuda:
  proc walkKernel(ws: ptr UncheckedArray[Walker], n: int32, f: Field,
                  tab: ptr MixtureTable, estep: eV, collisions: int32) {.global.} =
    let k = int32(blockIdx.x * blockDim.x + threadIdx.x)
    if k < n:
      walk(ws[k], f, tab[], estep, collisions)

type
  GpuWalkers* = object
    ## walkers, a collision table and the compiled kernel on the GPU
    engine: CudaEngine
    ws: DeviceMem[CUdeviceptr, Walker]
    tab: DeviceMem[CUdeviceptr, MixtureTable]
    n: int
    f: Field
    estep: eV

const threadsPerBlock = 128

proc kernelSource*(): string = walkKernel
  ## the CUDA source of the kernel

proc initGpuWalkers*(ws: openArray[Walker], f: Field, tab: MixtureTable, estep: eV,
                     exact = false): GpuWalkers =
  ## Copies the walkers and the table to the GPU. `exact`: no fused
  ## multiply-adds, so the arithmetic is that of the CPU (but `log`, `sin`,
  ## ... of CUDA, which may differ from the CPU's in the last bit).
  result.engine = newEngine(runCuda)
  if exact: result.engine.ingest(walkKernel, ["--fmad=false"])
  else: result.engine.ingest(walkKernel)
  result.n = ws.len
  result.f = f
  result.estep = estep
  result.ws = result.engine.allocDevice(Walker, ws.len)
  result.engine.copyToDeviceFromHost(result.ws, unsafeAddr ws[0], ws.len)
  result.tab = result.engine.allocDevice(MixtureTable, 1)
  result.engine.copyToDeviceFromHost(result.tab, unsafeAddr tab, 1)

proc walk*(g: var GpuWalkers, collisions: int32) =
  ## `collisions` real collisions of every walker
  let blocks = (g.n + threadsPerBlock - 1) div threadsPerBlock
  g.engine.launchKernel<<(blocks, threadsPerBlock)>>("walkKernel", g.ws,
    (int32(g.n), g.f, g.tab, g.estep, collisions))

proc read*(g: GpuWalkers, ws: var openArray[Walker]) =
  ## the walkers' state, back from the GPU
  doAssert ws.len == g.n
  g.engine.copyToHostFromDevice(addr ws[0], g.ws, g.n)

proc write*(g: var GpuWalkers, ws: openArray[Walker]) =
  ## replaces the walkers' state on the GPU
  doAssert ws.len == g.n
  g.engine.copyToDeviceFromHost(g.ws, unsafeAddr ws[0], g.n)
