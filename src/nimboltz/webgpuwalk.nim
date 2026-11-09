## Walkers (`walkers.nim`) on a GPU in a web page, through WebGPU: the
## WebAssembly build's counterpart of `gpuwalk.nim`, with its interface.
##
## The kernel is NimBoltz's own `walkWith`, compiled to WGSL at Nim compile
## time by Crucible (`webgpuSingle`, with -d:wgslPointerParams: the tables
## pass by pointer into storage, which WGSL allows with `requires
## unrestricted_pointer_parameters`). WGSL has no double precision, so the
## walkers walk in single precision (single.nim), as `gpu32` on CUDA.
##
## The host side is web/webgpu.js (an Emscripten JS library): the page's
## worker opens a WebGPU device before a run; reading back from the GPU is
## asynchronous, so the module is built with Asyncify. Build with
## -d:nimboltzWebGpu -d:wgslPointerParams -d:wgslPreciseMath and `--path:` to
## the tattletale root (web/build.sh does). WGSL's own acos, asin, sin, cos
## and log are coarse (acos and asin to 4096 ULP): with them the walkers
## drifted 0.8 % fast with thermal motion; with Crucible's precise ones a
## walker's first collision agrees with the CPU's to rounding.

import workspace/crucible
import units, tables, fields, freeflight, walkers, single
from simulation import MagboltzError

type
  KernelArgs = object
    ## the kernel's scalars and the field, in one binding (WebGPU has 8
    ## storage buffers per shader stage)
    f: Field
    estep: eV
    n, collisions, penning, counting: int32

const spectrumAt = tallyProcesses
  ## the tally: real collisions per process, then per energy bin

proc count(t: ptr UncheckedArray[uint32], i, ie: int) {.inline.} =
  discard atomicFetchAdd(addr t[i - 1], 1'u32)
  discard atomicFetchAdd(addr t[spectrumAt + ie - 1], 1'u32)

const threadsPerGroup = 64

template kernel(G: static FieldGeometry, T: typedesc): string =
  webgpuSingle:
    proc walkKernel(ws: ptr UncheckedArray[Walker[T]], tab: ptr UncheckedArray[TablesOf(T)],
                    fm: ptr UncheckedArray[FlightMaxima], a: KernelArgs,
                    tally: ptr UncheckedArray[uint32]) {.global, workgroup: (64, 1, 1).} =
      let k = int32(thread_position_in_grid.x)
      if k < a.n:
        if a.counting != 0:
          walkWith[G, T, ptr UncheckedArray[uint32]](ws[k], a.f, tab[0], addr fm[0], true,
                                                    a.estep, a.collisions, a.penning, tally)
        else:
          walkWith[G, T, NoTally](ws[k], a.f, tab[0], addr fm[0], true,
                                  a.estep, a.collisions, a.penning, NoTally())

const kernels = [
  [kernel(fgElectric, ColdGas), kernel(fgParallel, ColdGas),
   kernel(fgCrossed, ColdGas), kernel(fgOblique, ColdGas)],
  [kernel(fgElectric, WarmGas), kernel(fgParallel, WarmGas),
   kernel(fgCrossed, WarmGas), kernel(fgOblique, WarmGas)]]

func kernelSource*(G: static FieldGeometry, T: typedesc): string =
  ## the WGSL source of the kernel for geometry G and target T
  kernels[ord(T is WarmGas)][ord(G)]

# web/webgpu.js
proc nbgpuReady(): cint {.importc.}
proc nbgpuPipeline(src: cstring, len: cint): cint {.importc.}
proc nbgpuBuffer(bytes: cdouble): cint {.importc.}
proc nbgpuWrite(id: cint, offset: cdouble, p: pointer, len: cdouble) {.importc.}
proc nbgpuRun(pipe, ws, tab, fm, args, tally, groups: cint) {.importc.}
proc nbgpuRead(id: cint, p: pointer, len: cdouble): cint {.importc.}
proc nbgpuFree(id: cint) {.importc.}

proc webGpuReady*(): bool = nbgpuReady() != 0
  ## whether the page opened a WebGPU device

proc failed() =
  raise newException(MagboltzError, "walkers on a GPU: WebGPU failed (see the printout)")

proc readBack(id: cint, p: pointer, len: int) =
  if nbgpuRead(id, p, cdouble(len)) == 0: failed()

type
  Buffers = object
    ws, tab, fm, args, tally: cint

  GpuWalkers*[G: static FieldGeometry; T] = object
    ## walkers, their collision tables and the kernel on the GPU, all in
    ## single precision
    pipe: cint
    bufs: ref Buffers
    stride: int            ## bytes per walker
    n: int
    f: Field
    estep: eV
    penning: int32
    tallies: ref Tallies   ## the counts so far, on the host

proc `=destroy`(b: Buffers) =
  for id in [b.ws, b.tab, b.fm, b.args, b.tally]:
    if id > 0: nbgpuFree(id)

proc upload(bytes: openArray[byte]): cint =
  result = nbgpuBuffer(cdouble(max(bytes.len, 4)))
  if bytes.len > 0: nbgpuWrite(result, 0, unsafeAddr bytes[0], cdouble(bytes.len))

proc packed[T](x: T): seq[byte] = single(x)
proc packed[T](xs: openArray[T]): seq[byte] = single(xs)
  ## (`single` is also a parameter of `initGpuWalkers`)

template one[X](x: X): untyped =
  ## x as an openArray of one, in place (tables are large for the stack)
  cast[ptr UncheckedArray[X]](unsafeAddr x).toOpenArray(0, 0)

proc clearTally*(g: var GpuWalkers) =
  ## zeroes the collision counts
  g.tallies[] = Tallies()

proc readTally*(g: GpuWalkers, t: var Tallies) =
  ## the collision counts so far
  t = g.tallies[]

proc fetchTally(g: var GpuWalkers) =
  ## adds the GPU's counts to the host's, and zeroes them
  var c = newSeq[uint32](tallyProcesses + nEnergyBins)
  readBack(g.bufs.tally, addr c[0], c.len * 4)
  for i in 1 .. tallyProcesses:
    g.tallies.processes[i] = g.tallies.processes[i] + float64(c[i - 1])
  for i in 1 .. nEnergyBins:
    g.tallies.spectrum[i] = g.tallies.spectrum[i] + float64(c[spectrumAt + i - 1])
  zeroMem(addr c[0], c.len * 4)
  nbgpuWrite(g.bufs.tally, 0, addr c[0], cdouble(c.len * 4))

proc initGpuWalkers*[G: static FieldGeometry; T](ws: openArray[Walker[T]], f: Field,
                                                 tab: TablesOf(T), fm: FlightMaxima,
                                                 estep: eV, penning = false,
                                                 exact = false, single = true): GpuWalkers[G, T] =
  ## Copies the walkers, the table and its range maxima to the GPU, in
  ## single precision (`exact`, `single`: for the interface of gpuwalk.nim)
  if nbgpuReady() == 0:
    raise newException(MagboltzError, "walkers on a GPU: no WebGPU device")
  const src = kernelSource(G, T)
  result.pipe = nbgpuPipeline(cstring(src), cint(src.len))
  if result.pipe == 0: failed()
  result.stride = singleSize(Walker[T])
  result.penning = int32(ord(penning))
  result.n = ws.len
  result.f = f
  result.estep = estep
  new result.bufs
  result.bufs.ws = upload(packed(ws))
  result.bufs.tab = upload(packed(one(tab)))
  result.bufs.fm = upload(packed(one(fm)))
  result.bufs.args = upload(packed(KernelArgs()))
  result.bufs.tally = upload(newSeq[byte](4 * (tallyProcesses + nEnergyBins)))
  new result.tallies

const perLaunch = 1000'i32
  ## collisions per walker in one dispatch: a long one may trip the
  ## browser's GPU watchdog

proc walk*(g: var GpuWalkers, collisions: int32, counting = true) =
  ## `collisions` real collisions of every walker, counted in the tally
  ## unless not `counting`
  let groups = cint((g.n + threadsPerGroup - 1) div threadsPerGroup)
  var left = collisions
  while left > 0:
    let c = min(left, perLaunch)
    let a = single(KernelArgs(f: g.f, estep: g.estep, n: int32(g.n), collisions: c,
                              penning: g.penning, counting: int32(ord(counting))))
    nbgpuWrite(g.bufs.args, 0, unsafeAddr a[0], cdouble(a.len))
    nbgpuRun(g.pipe, g.bufs.ws, g.bufs.tab, g.bufs.fm, g.bufs.args, g.bufs.tally, groups)
    left = left - c
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
  readBack(g.bufs.ws, addr bytes[0], bytes.len)
  var pos = 0
  for w in ws.mitems: unpack32(w, bytes, pos)

proc write*[T](g: var GpuWalkers, ws: openArray[Walker[T]]) =
  ## replaces the walkers' state on the GPU
  doAssert ws.len == g.n
  let bytes = single(ws)
  nbgpuWrite(g.bufs.ws, 0, unsafeAddr bytes[0], cdouble(bytes.len))
