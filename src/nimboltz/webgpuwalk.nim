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
import units, rm48, xoshiro, tables, fields, freeflight, walkers, single, soa
from simulation import MagboltzError

type
  KernelArgs = object
    ## the kernel's scalars and the field, in one binding (WebGPU has 8
    ## storage buffers per shader stage)
    f: Field
    estep: eV
    n, stride, collisions, penning, counting: int32

const spectrumAt = tallyProcesses
  ## the tally: real collisions per process, then per energy bin

proc count(t: ptr UncheckedArray[uint32], i, ie: int) {.inline.} =
  discard atomicFetchAdd(addr t[i - 1], 1'u32)
  discard atomicFetchAdd(addr t[spectrumAt + ie - 1], 1'u32)

const threadsPerGroup = 64

# The walkers lie in slots (soa.nim), as in gpuwalk.nim: a thread copies
# its walker into a variable of its own, walks, and copies it back.
template kernel(G: static FieldGeometry, T, R: typedesc): string =
  webgpuSingle:
    proc walkKernel(fs: ptr UncheckedArray[float64], us: ptr UncheckedArray[uint32],
                    tab: ptr UncheckedArray[TablesOf(T)],
                    fm: ptr UncheckedArray[FlightMaxima], a: KernelArgs,
                    tally: ptr UncheckedArray[uint32]) {.global, workgroup: (64, 1, 1).} =
      let k = int32(thread_position_in_grid.x)
      if k < a.n:
        var w: Walker[T, R]
        load(w, fs, us, k, a.stride)
        if a.counting != 0:
          walkWith[G, T, ptr UncheckedArray[uint32], R](w, a.f, tab[0], addr fm[0], gpuWindows,
                                                        a.estep, a.collisions, a.penning, tally)
        else:
          walkWith[G, T, NoTally, R](w, a.f, tab[0], addr fm[0], gpuWindows,
                                     a.estep, a.collisions, a.penning, NoTally())
        store(w, fs, us, k, a.stride)

template kernels(R: typedesc): untyped =
  [[kernel(fgElectric, ColdGas, R), kernel(fgParallel, ColdGas, R),
    kernel(fgCrossed, ColdGas, R), kernel(fgOblique, ColdGas, R)],
   [kernel(fgElectric, WarmGas, R), kernel(fgParallel, WarmGas, R),
    kernel(fgCrossed, WarmGas, R), kernel(fgOblique, WarmGas, R)]]

const xoshiroKernels = kernels(Xoshiro128Plus)
when defined(gpuRanmar):
  const ranmarKernels = kernels(Ranmar)

func kernelSource*(G: static FieldGeometry, T, R: typedesc): string =
  ## the WGSL source of the kernel for geometry G, target T and generator R
  ## (xoshiro128+; Ranmar `toSingle` if built with -d:gpuRanmar)
  when R is Xoshiro128Plus: xoshiroKernels[ord(T is WarmGas)][ord(G)]
  elif R is Ranmar and defined(gpuRanmar): ranmarKernels[ord(T is WarmGas)][ord(G)]
  else:
    {.error: "walkers with WebGPU: xoshiro128+ (Ranmar with -d:gpuRanmar)".}

# web/webgpu.js
proc nbgpuReady(): cint {.importc.}
proc nbgpuPipeline(src: cstring, len: cint): cint {.importc.}
proc nbgpuBuffer(bytes: cdouble): cint {.importc.}
proc nbgpuWrite(id: cint, offset: cdouble, p: pointer, len: cdouble) {.importc.}
proc nbgpuRun(pipe, fs, us, tab, fm, args, tally, groups: cint) {.importc.}
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
    fs, us, tab, fm, args, tally: cint

  GpuWalkers*[G: static FieldGeometry; T, R] = object
    ## walkers, their collision tables and the kernel on the GPU, all in
    ## single precision
    pipe: cint
    bufs: ref Buffers
    stride: int            ## walkers room is made for (each slot's length)
    n: int
    f: Field
    estep: eV
    penning: int32
    tallies: ref Tallies   ## the counts so far, on the host

proc `=destroy`(b: Buffers) =
  for id in [b.fs, b.us, b.tab, b.fm, b.args, b.tally]:
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

proc initGpuWalkers*[G: static FieldGeometry; T, R](ws: openArray[Walker[T, R]], f: Field,
                                                    tab: TablesOf(T), fm: FlightMaxima,
                                                    estep: eV, penning = false,
                                                    exact = false, single = true): GpuWalkers[G, T, R] =
  ## Copies the walkers, the table and its range maxima to the GPU, in
  ## single precision (`exact`, `single`: for the interface of gpuwalk.nim)
  if nbgpuReady() == 0:
    raise newException(MagboltzError, "walkers on a GPU: no WebGPU device")
  const src = kernelSource(G, T, R)
  result.pipe = nbgpuPipeline(cstring(src), cint(src.len))
  if result.pipe == 0: failed()
  result.stride = ws.len
  result.penning = int32(ord(penning))
  result.n = ws.len
  result.f = f
  result.estep = estep
  new result.bufs
  let (fb, ub) = slotBytes(ws, result.stride, single = true)
  result.bufs.fs = upload(fb)
  result.bufs.us = upload(ub)
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
    let a = single(KernelArgs(f: g.f, estep: g.estep, n: int32(g.n), stride: int32(g.stride),
                              collisions: c,
                              penning: g.penning, counting: int32(ord(counting))))
    nbgpuWrite(g.bufs.args, 0, unsafeAddr a[0], cdouble(a.len))
    nbgpuRun(g.pipe, g.bufs.fs, g.bufs.us, g.bufs.tab, g.bufs.fm, g.bufs.args, g.bufs.tally, groups)
    left = left - c
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
  var fb = newSeq[byte](max(c[0] * g.stride * 4, 4))
  var ub = newSeq[byte](max(c[1] * g.stride * 4, 4))
  readBack(g.bufs.fs, addr fb[0], fb.len)
  readBack(g.bufs.us, addr ub[0], ub.len)
  fromSlotBytes(ws, fb, ub, g.stride, single = true)

proc write*[T, R](g: var GpuWalkers, ws: openArray[Walker[T, R]]) =
  ## replaces the walkers' state on the GPU
  doAssert ws.len == g.n
  let (fb, ub) = slotBytes(ws, g.stride, single = true)
  nbgpuWrite(g.bufs.fs, 0, unsafeAddr fb[0], cdouble(fb.len))
  nbgpuWrite(g.bufs.us, 0, unsafeAddr ub[0], cdouble(ub.len))
