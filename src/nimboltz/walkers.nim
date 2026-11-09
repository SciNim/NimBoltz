## Independent electrons drifting in a gas, in any field geometry, at 0 K
## (`Walker[ColdGas]`) or with thermal motion of the gas (`Walker[WarmGas]`).
##
## MONTE's collision loop (`drift.runSegment`) without its
## bookkeeping: each walker is one electron with its own random numbers
## (`Ranmar`), and keeps only the sums that give the transport coefficients
## of the ensemble (`ensemble.nim`). Many walkers give the same estimates
## as one long trajectory, not the same numbers; they can run side by side,
## on CPU threads or one per GPU thread.
##
## Plain Nim the GPU compiler takes too: no `break`, `continue` or `case`,
## no references, no `bool` arguments, `int32` counters.

import std/math
import units, rm48, tables, vec3, fields, collisions, freeflight, thermal

const tallyProcesses* = max(maxMixProcesses, maxGases * maxGasProcesses)
  ## processes counted: of the mixture's table, or of the per-gas tables
  ## (gas k's process i at (k - 1) · maxGasProcesses + i)

type
  Walker*[T] = object
    el*: Electron        ## state after the last collision
    target*: T           ## ColdGas: null-collision bounds, raised as needed;
                         ## WarmGas: also the walker's Maxwell numbers
    rng*: Ranmar
    pos*: Vector[m]      ## position
    st*: ps              ## time
    e2*: eV•ps           ## ∫ e dt (B ⊥ E and oblique B: Σ e/ν over the
                         ## collisions, as MONTEB and MONTEC estimate it)
    ions*: int32         ## ionisations (with Auger electrons and Penning transfers)
    attachments*: int32
    collisions*: int32   ## real collisions
    nulls*: int32        ## null collisions
    tmax*: ps            ## longest flight
    pos1*: Vector[m]     ## position and time at the mark (`mark`), for
    st1*: ps             ## the diffusion estimate of `ensemble.nim`

  Tally* = object
    ## real collisions per process of the table and per energy bin, kept
    ## outside the walkers (a GPU adds to one tally atomically: see
    ## gpuwalk.nim); 0-based: process i at [i - 1], bin ie at [ie - 1]
    processes*: ptr UncheckedArray[float64]
    spectrum*: ptr UncheckedArray[float64]

  NoTally* = object
    ## counts nothing

proc count*(t: Tally, i, ie: int) {.inline.} =
  t.processes[i - 1] = t.processes[i - 1] + 1.0
  t.spectrum[ie - 1] = t.spectrum[ie - 1] + 1.0

proc count*(t: NoTally, i, ie: int) {.inline.} = discard

template TablesOf*(T: typedesc): typedesc =
  ## the collision tables a walker on target T uses
  when T is ColdGas: MixtureTable else: GasTables

proc initWalker*(start: Electron, target: ColdGas, seed: int): Walker[ColdGas] =
  Walker[ColdGas](el: start, target: target, rng: initRanmar(seed))

proc initWalker*(start: Electron, gas: ThermalGas, seed: int): Walker[WarmGas] =
  result = Walker[WarmGas](el: start, rng: initRanmar(seed))
  result.target = initWarmGas(gas, result.rng)

func collisionFrequency(tab: MixtureTable, target: ColdGas, ie: int): ps⁻¹ {.inline.} =
  tab.tcf[ie]

func collisionFrequency(tabs: GasTables, target: WarmGas, ie: int): ps⁻¹ {.inline.} =
  ## of all gases
  result = 0.0.ps⁻¹
  for k in 1 .. int(target.gas.count): result = result + tabs[k].tcf[ie]

# Free flights of walkers: the null-collision method of `freeflight.nextHit`,
# the same distributions, other numbers. What depends only on the start of a
# flight is computed once, divisions become multiplications, and the
# logarithm of the flight time is single precision: double precision is
# slow on GPUs (1/64 of single on gaming ones), and a flight time needs no
# more digits.

func lnSingle(r: float64): float64 {.inline.} = float64(ln(float32(r)))

func binOf(e: eV, rstep: eV⁻¹): int {.inline.} =
  ## energy bin of e, rstep = 1 / estep
  clamp(int(trunc(float64(e * rstep)) + 1.0), 1, nEnergyBins)

proc nextHit[G: static FieldGeometry](f: Field, tab: MixtureTable, target: var ColdGas,
                                      rstep: eV⁻¹, el: Electron, rng: var Ranmar,
                                      nnull: var int32): Hit =
  ## gas at 0 K: next real collision of an electron leaving `el`
  let e0 = el  # a copy: not reloaded after each random number
  let band = clamp(int(trunc(e0.e / target.bandWidth) + 1.0), 1, target.bounds.len)
  var tlim = target.bounds[band]
  var itlim = 1.0 / tlim
  when G in {fgElectric, fgParallel}:
    let ap = e0.dir.z * f.f2 * sqrt(e0.e)
  var tdash = 0.0.ps
  var real = false
  while not real:
    let lnr = lnSingle(drand48(rng))
    tdash = tdash - lnr * itlim
    when G in {fgElectric, fgParallel}:
      result.fl.t = tdash
      result.fl.e = e0.e + (ap + f.bp * tdash) * tdash
    else:
      result.fl = fly[G](f, e0, tdash)
    result.ie = binOf(result.fl.e, rstep)
    let tcf = tab.tcf[result.ie]
    if tcf > tlim:
      # bound too low: undo this flight, raise the bound
      tdash = tdash + lnr * itlim
      target.bounds[band] = 1.05 * target.bounds[band]
      tlim = target.bounds[band]
      itlim = 1.0 / tlim
    else:
      real = drand48(rng) * tlim <= tcf
      if not real: inc nnull
  finishFlight[G](f, result.fl)
  result.e = result.fl.e
  result.din = directionAt[G](f, e0, result.fl)
  result.gas = 1

proc nextHit[G: static FieldGeometry](f: Field, tabs: GasTables, target: var WarmGas,
                                      rstep: eV⁻¹, el: Electron, rng: var Ranmar,
                                      nnull: var int32): Hit =
  ## thermal motion: next real collision; energy and direction in the rest
  ## frame of the molecule hit
  let e0 = el
  let itcfmx = 1.0 / target.gas.tcfmx
  let ivfac2 = 1.0 / f.vfac2
  when G == fgElectric:
    # v(t) = v0 + a t along E, no square roots in the loop
    let v0 = e0.dir * (f.vfac * sqrt(e0.e))
    let az = f.efield * f.const5 * f.vfac
    let ap = e0.dir.z * f.f2 * sqrt(e0.e)
  var tdash = 0.0.ps
  var vrel: Vector[m•ps⁻¹]
  var real = false
  while not real:
    tdash = tdash - lnSingle(drand48(rng)) * itcfmx
    when G == fgElectric:
      result.fl.t = tdash
      result.fl.e = e0.e + (ap + f.bp * tdash) * tdash
      var vel = v0
      vel.z = vel.z + az * tdash
    else:
      result.fl = fly[G](f, e0, tdash)
      finishFlight[G](f, result.fl)
      let vel = velocityAt[G](f, e0, result.fl)
    result.gas = target.gas.pickGas(drand48(rng))
    result.vgas = target.maxwell.velocity(target.gas.vtmb[result.gas], rng)
    vrel = vel - result.vgas
    result.e = norm2(vrel) * ivfac2
    result.ie = binOf(result.e, rstep)
    real = not (drand48(rng) * target.gas.tcfmax[result.gas] > tabs[result.gas].tcf[result.ie])
    if not real: inc nnull
  result.din = vrel * (1.0 / (f.vfac * sqrt(result.e)))

func reset*(w: var Walker) =
  ## clears the sums (after the walker has forgotten its start)
  w.pos = vec3(0.0.m, 0.0.m, 0.0.m)
  w.st = 0.0.ps
  w.e2 = 0.0.eV•ps
  w.ions = 0
  w.attachments = 0
  w.collisions = 0
  w.nulls = 0
  w.tmax = 0.0.ps
  w.pos1 = w.pos
  w.st1 = w.st

func mark*(w: var Walker) =
  ## remembers the position and time (halfway through a run): the spread
  ## of the positions grows as 2 D t + const, and the constant cancels in
  ## the growth from the mark on
  w.pos1 = w.pos
  w.st1 = w.st

proc walk*[G: static FieldGeometry; T, C](w: var Walker[T], f: Field, tabs: TablesOf(T),
                                          estep: eV, collisions, penning: int32, tally: C) =
  ## `collisions` real collisions, with the random numbers drawn in the
  ## order of MONTE (MONTET); `penning` ≠ 0: Penning transfers; each
  ## collision is counted in `tally` (`Tally`, `NoTally` or a GPU's)
  mixin count
  let rstep = 1.0 / estep
  for _ in 0'i32 ..< collisions:
    let el = w.el  # a copy, so its square roots are computed once
    let h = nextHit[G](f, tabs, w.target, rstep, el, w.rng, w.nulls)
    if h.fl.t > w.tmax: w.tmax = h.fl.t
    when G in {fgElectric, fgParallel}:
      w.e2 = w.e2 + energyIntegral[G](f, el, h.fl)
    else:
      w.e2 = w.e2 + (estep / 2.0 + estep * float64(h.ie - 1)) /
                    collisionFrequency(tabs, w.target, h.ie)
    discard advance[G](f, el, h.fl, w.pos)
    w.st = w.st + h.fl.t
    template tab: untyped = tableOf(tabs, h)
    let r2 = drand48(w.rng)
    let i = findProcess(tab, h.ie, r2)
    let ei = energyLoss(tab, i, h.e, w.rng, moderation = true)
    when T is ColdGas: tally.count(i, h.ie)
    else: tally.count((h.gas - 1) * maxGasProcesses + i, h.ie)
    let kind = (tab.iarry[i] - 1'i32) mod 5'i32 + 1'i32  # 2: ionisation, 3: attachment
    if kind == 2'i32: w.ions = w.ions + 1'i32
    if kind == 3'i32: w.attachments = w.attachments + 1'i32
    if tab.ipn[i] > 0'i32: w.ions = w.ions + tab.augerCount[i]
    if penning != 0'i32 and tab.penningFraction[i] != 0.0:
      let ran = drand48(w.rng)
      if not (ran > tab.penningFraction[i]): w.ions = w.ions + 1'i32
    let s = scatter(tab, i, h.ie, h.e, ei, h.din, w.rng,
                    r31First = (G == fgElectric and T is ColdGas))
    w.el = leave[G](w.target, f, s, h)
  w.collisions = w.collisions + collisions

func driftVelocity*(ws: openArray[Walker]): μm•ns⁻¹ =
  ## Σ z / Σ t over all walkers
  var z = 0.0.m
  var t = 0.0.ps
  for w in ws:
    z = z + w.pos.z
    t = t + w.st
  (z / t).to(μm•ns⁻¹)

func meanEnergy*(ws: openArray[Walker]): eV =
  ## Σ ∫ e dt / Σ t over all walkers
  var e2 = 0.0.eV•ps
  var t = 0.0.ps
  for w in ws:
    e2 = e2 + w.e2
    t = t + w.st
  e2 / t

type
  Tallies* = object
    ## a CPU tally: real collisions per process and per energy bin
    processes*: array[1..tallyProcesses, float64]
    spectrum*: array[1..nEnergyBins, float64]

func add*(a: var Tallies, b: Tallies) =
  for i in 1 .. tallyProcesses: a.processes[i] = a.processes[i] + b.processes[i]
  for i in 1 .. nEnergyBins: a.spectrum[i] = a.spectrum[i] + b.spectrum[i]

func tally(t: var Tallies): Tally =
  Tally(processes: cast[ptr UncheckedArray[float64]](addr t.processes[1]),
        spectrum: cast[ptr UncheckedArray[float64]](addr t.spectrum[1]))

when compileOption("threads"):
  type
    Job[T] = object
      ws: ptr UncheckedArray[Walker[T]]
      n: int
      f: ptr Field
      tab: ptr TablesOf(T)
      estep: eV
      collisions, penning: int32
      tallies: ptr Tallies  # nil: count nothing

  proc walkJob[G: static FieldGeometry; T](j: ptr Job[T]) {.thread.} =
    for k in 0 ..< j.n:
      if j.tallies == nil:
        walk[G, T, NoTally](j.ws[k], j.f[], j.tab[], j.estep, j.collisions, j.penning,
                            NoTally())
      else:
        walk[G, T, Tally](j.ws[k], j.f[], j.tab[], j.estep, j.collisions, j.penning,
                          tally(j.tallies[]))

proc walkAll*[G: static FieldGeometry; T](ws: var openArray[Walker[T]], f: Field,
                                          tab: TablesOf(T), estep: eV, collisions,
                                          penning: int32, threads = 1,
                                          tallies: ptr Tallies = nil) =
  ## `collisions` real collisions of every walker, on `threads` CPU threads;
  ## counted in `tallies` unless nil
  when compileOption("threads"):
    if threads > 1 and ws.len > 1:
      let n = min(threads, ws.len)
      var jobs = newSeq[Job[T]](n)
      var ts = newSeq[Thread[ptr Job[T]]](n)
      var parts = newSeq[Tallies](if tallies == nil: 0 else: n)
      for k in 0 ..< n:
        let a = k * ws.len div n
        let b = (k + 1) * ws.len div n
        jobs[k] = Job[T](ws: cast[ptr UncheckedArray[Walker[T]]](addr ws[a]), n: b - a,
                      f: unsafeAddr f, tab: unsafeAddr tab, estep: estep,
                      collisions: collisions, penning: penning,
                      tallies: if tallies == nil: nil else: addr parts[k])
        createThread(ts[k], walkJob[G, T], addr jobs[k])
      joinThreads(ts)
      for p in parts: tallies[].add p
      return
  for w in ws.mitems:
    if tallies == nil: walk[G, T, NoTally](w, f, tab, estep, collisions, penning, NoTally())
    else: walk[G, T, Tally](w, f, tab, estep, collisions, penning, tally(tallies[]))

type
  CpuWalkers*[G: static FieldGeometry; T] = object
    ## walkers on CPU threads, with the interface of `gpuwalk.GpuWalkers`
    ws: seq[Walker[T]]
    f: Field
    tab: ptr TablesOf(T)
    estep: eV
    penning: int32
    threads: int
    tallies: ref Tallies

proc initCpuWalkers*[G: static FieldGeometry; T](ws: openArray[Walker[T]], f: Field,
                                                 tab: ptr TablesOf(T), estep: eV,
                                                 penning = false,
                                                 threads = 1): CpuWalkers[G, T] =
  ## `tab` must outlive the walkers
  result = CpuWalkers[G, T](ws: @ws, f: f, tab: tab, estep: estep,
                            penning: int32(ord(penning)), threads: threads)
  new result.tallies

proc walk*(c: var CpuWalkers, collisions: int32, counting = true) =
  ## `collisions` real collisions of every walker, counted in the tally
  ## unless not `counting`
  walkAll[c.G, c.T](c.ws, c.f, c.tab[], c.estep, collisions, c.penning, c.threads,
                    if counting: addr c.tallies[] else: nil)

proc read*[T](c: CpuWalkers, ws: var openArray[Walker[T]]) =
  for k in 0 ..< ws.len: ws[k] = c.ws[k]

proc write*[T](c: var CpuWalkers, ws: openArray[Walker[T]]) =
  for k in 0 ..< ws.len: c.ws[k] = ws[k]

proc shrink*(c: var CpuWalkers, n: int) =
  ## keeps the first n walkers
  c.ws.setLen(n)

proc clearTally*(c: var CpuWalkers) = c.tallies[] = Tallies()

proc readTally*(c: CpuWalkers, t: var Tallies) = t = c.tallies[]
