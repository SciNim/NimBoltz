## Independent electrons drifting in a gas, in any field geometry, at 0 K
## (`Walker[ColdGas]`) or with thermal motion of the gas (`Walker[WarmGas]`).
##
## MONTE's collision loop (`drift.runSegment`) without its
## bookkeeping: each walker is one electron with its own random numbers
## (`Ranmar`, or another `WalkerRng`), and keeps only the sums that give the transport coefficients
## of the ensemble (`ensemble.nim`). Many walkers give the same estimates
## as one long trajectory, not the same numbers; they can run side by side,
## on CPU threads or one per GPU thread.
##
## Plain Nim the GPU compiler takes too: no `break`, `continue` or `case`,
## no references, no `bool` arguments, `int32` counters.

import std/math
import units, rm48, xoshiro, tables, vec3, fields, collisions, freeflight, thermal
export xoshiro

const tallyProcesses* = max(maxMixProcesses, maxGases * maxGasProcesses)
  ## processes counted: of the mixture's table, or of the per-gas tables
  ## (gas k's process i at (k - 1) · maxGasProcesses + i)

type
  Walker*[T; R = Ranmar] = object
    el*: Electron        ## state after the last collision
    target*: T           ## ColdGas: null-collision bounds, raised as needed;
                         ## WarmGas: also the walker's Maxwell numbers
    rng*: R              ## a `WalkerRng`
    pos*: Vector[m]      ## position
    st*: ps              ## time
    e2*: eV•ps           ## ∫ e dt (B ⊥ E and oblique B: Σ e/ν over the
                         ## collisions, as MONTEB and MONTEC estimate it)
    ions*: int32         ## ionisations (with Auger electrons and Penning transfers)
    attachments*: int32
    collisions*: int32   ## real collisions
    nulls*: int32        ## null collisions
    beyond*: int32       ## thermal motion: collisions above the bound (`windowHit`)
    tmax*: ps            ## longest flight

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

proc initWalker*[R: WalkerRng](start: Electron, target: ColdGas, seed: int,
                               rng: typedesc[R]): Walker[ColdGas, R] =
  ## with random numbers of type R
  Walker[ColdGas, R](el: start, target: target, rng: initRng(R, seed))

proc initWalker*(start: Electron, target: ColdGas, seed: int): Walker[ColdGas] =
  initWalker(start, target, seed, Ranmar)

proc initWalker*[R: WalkerRng](start: Electron, gas: ThermalGas,
                               bands: array[1..maxGases, NullBounds], efinal: eV, seed: int,
                               rng: typedesc[R]): Walker[WarmGas, R] =
  ## `bands`: per-gas null-collision bounds in 8 bands of the starting
  ## energy (`mixture.nullBounds` of each gas's table)
  result = Walker[WarmGas, R](el: start, rng: initRng(R, seed))
  result.target = initWarmGas(gas, result.rng)
  result.target.bands = bands
  result.target.bandWidth = efinal / float64(bands[1].len)

proc initWalker*(start: Electron, gas: ThermalGas, bands: array[1..maxGases, NullBounds],
                 efinal: eV, seed: int): Walker[WarmGas] =
  initWalker(start, gas, bands, efinal, seed, Ranmar)

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
                                      rstep: eV⁻¹, el: Electron, rng: var WalkerRng,
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
                                      rstep: eV⁻¹, el: Electron, rng: var WalkerRng,
                                      nnull: var int32): Hit =
  ## thermal motion: next real collision; energy and direction in the rest
  ## frame of the molecule hit. Unlike Magboltz (one bound per gas over all
  ## energies), the bounds of each gas depend on the starting energy, as at
  ## 0 K, and are raised when too low: fewer null collisions.
  let e0 = el
  let band = clamp(int(trunc(e0.e / target.bandWidth) + 1.0), 1, target.bands[1].len)
  let n = int(target.gas.count)
  var total = 0.0.ps⁻¹
  for k in 1 .. n: total = total + target.bands[k][band]
  var itotal = 1.0 / total
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
    let lnr = lnSingle(drand48(rng))
    tdash = tdash - lnr * itotal
    when G == fgElectric:
      result.fl.t = tdash
      result.fl.e = e0.e + (ap + f.bp * tdash) * tdash
      var vel = v0
      vel.z = vel.z + az * tdash
    else:
      result.fl = fly[G](f, e0, tdash)
      finishFlight[G](f, result.fl)
      let vel = velocityAt[G](f, e0, result.fl)
    # the gas, by its share of the bound
    let rg = drand48(rng) * total
    var k = 1
    var below = target.bands[1][band]
    while below < rg and k < n:
      k = k + 1
      below = below + target.bands[k][band]
    result.gas = k
    result.vgas = target.maxwell.velocity(target.gas.vtmb[k], rng)
    vrel = vel - result.vgas
    result.e = norm2(vrel) * ivfac2
    result.ie = binOf(result.e, rstep)
    let tcf = tabs[k].tcf[result.ie]
    let bound = target.bands[k][band]
    if tcf > bound:
      # bound too low: undo this flight, raise the bound
      tdash = tdash + lnr * itotal
      target.bands[k][band] = 1.05 * bound
      total = total + 0.05 * bound
      itotal = 1.0 / total
    else:
      real = drand48(rng) * bound <= tcf
      if not real: inc nnull
  result.din = vrel * (1.0 / (f.vfac * sqrt(result.e)))

# Exact bounds per time window: from the electron's state, the energies of
# the next `windowCollisions` mean free times are bounded, and the largest
# collision frequency over them (two lookups in a `RangeMax`) bounds the
# frequency strictly. B = 0 and B ∥ E: the energy along a flight is
# e0 + ap t + bp t², its range over a window exact. Any other B: the
# magnetic force does no work, and the electric one changes √e at most at
# the rate `Field.slope`, so √e stays within slope · τ of its value at the
# start of the window (a few percent over a window). Where the bound is
# well above the frequency at the start (a Ramsauer minimum: more than
# `shrinkAbove` times as many trials), the window shrinks to
# `windowCollisions` times 1 / bound and gets its own, lower bound. Longer
# windows give more null collisions but fewer bounds. On a CPU, where a
# bound costs more than a trial, 3 and 2 are about the fastest (1 and 1:
# 25 % slower, 2 and 1: 15 %; six decks, 2026-11-12); on a GPU 2 and 1 (3
# and 2: Ar/CF₄ in double precision 17 % slower, in single 1–5 %). The flight time to the next trial is one exponential "depth"
# (−ln r), used up window by window at each window's bound (a piecewise
# constant bound of the frequency): an empty window costs neither a
# logarithm nor a random number. Nothing to raise, nothing missed, few null
# collisions. Thermal motion: per gas, in the molecule's frame, with a
# margin of 6 thermal speeds (beyond it: once in 10¹⁵ collisions, taken as
# real and counted in `Walker.beyond`).

type
  Bounds* = enum
    ## of the walkers' flights: Magboltz's (`nextHit`), or exact per time
    ## window (`windowHit`) with the window lengths for CPUs or for GPUs
    magboltzBounds, cpuWindows, gpuWindows

const
  rangeLevels = 12         ## 2¹² ≥ nEnergyBins
  windowCollisions: array[cpuWindows..gpuWindows, float64] = [3.0, 2.0]
    ## window length in mean free times (at most)
  shrinkAbove: array[cpuWindows..gpuWindows, float64] = [2.0, 1.0]
    ## a window shrinks above this many times its trials
  speedMargin = 6.0        ## molecule speeds in thermal speeds (`vtmb`)

type
  RangeMax* = object
    ## largest collision frequency over a range of energy bins in two
    ## lookups: level k has the maxima of 2^k bins from each bin on
    level*: array[0..rangeLevels - 1, array[1..nEnergyBins, ps⁻¹]]

  FlightMaxima* = object
    ## per gas (at 0 K: the mixture's, as gas 1)
    gas*: array[1..maxGases, RangeMax]
    levelOf*: array[1..nEnergyBins, int32]  ## ⌊log₂ n⌋: the level for n bins
    widthOf*: array[1..nEnergyBins, int32]  ## 2^⌊log₂ n⌋ (WGSL has no signed shifts)

proc fill(r: var RangeMax, tcf: array[1..nEnergyBins, ps⁻¹]) =
  r.level[0] = tcf
  var w = 1
  for k in 1 ..< rangeLevels:
    for i in 1 .. nEnergyBins:
      r.level[k][i] = if i + w <= nEnergyBins: max(r.level[k - 1][i], r.level[k - 1][i + w])
                      else: r.level[k - 1][i]
    w = 2 * w

proc fillLevels(fm: var FlightMaxima) =
  var k = 0'i32
  for n in 1 .. nEnergyBins:
    if 2 shl k <= n: k = k + 1
    fm.levelOf[n] = k
    fm.widthOf[n] = 1'i32 shl k

proc flightMaxima*(tab: MixtureTable): ref FlightMaxima =
  new result
  result.gas[1].fill(tab.tcf)
  result[].fillLevels()

proc flightMaxima*(tabs: GasTables): ref FlightMaxima =
  new result
  for k in 1 .. maxGases: result.gas[k].fill(tabs[k].tcf)
  result[].fillLevels()

func rangeMax(fm: FlightMaxima, gas, lo, hi: int): ps⁻¹ {.inline.} =
  ## of gas over the bins lo .. hi
  let n = hi - lo + 1
  let k = int(fm.levelOf[n])
  max(fm.gas[gas].level[k][lo], fm.gas[gas].level[k][hi - int(fm.widthOf[n]) + 1])

func energyRange(e0: eV, ap: eV•ps⁻¹, bp: eV•ps⁻², t1, t2: ps): (eV, eV) {.inline.} =
  ## lowest and highest of e0 + ap t + bp t² for t in t1 .. t2
  let a = e0 + (ap + bp * t1) * t1
  let b = e0 + (ap + bp * t2) * t2
  var lo = min(a, b)
  if bp > 0.0.eV•ps⁻²:
    let tv = -ap / (2.0 * bp)  # the vertex, a minimum
    if tv > t1 and tv < t2: lo = min(lo, e0 + (ap + bp * tv) * tv)
  (lo, max(a, b))

func energyRange(ew: eV, slope: eV¹ᐟ²•ps⁻¹, tau: ps): (eV, eV) {.inline.} =
  ## bounds of the energy over a time tau from energy ew, the change of √e
  ## at most `slope`
  let s = sqrt(max(ew, 0.0.eV))
  let d = slope * tau
  let lo = max(s - d, 0.0.eV¹ᐟ²)
  (lo * lo, (s + d) * (s + d))

func flightEnergy[G: static FieldGeometry](f: Field, e0: Electron, ap: eV•ps⁻¹,
                                           t: ps): eV {.inline.} =
  ## energy after a flight of time t (B = 0, B ∥ E: ap from e0)
  when G in {fgElectric, fgParallel}: e0.e + (ap + f.bp * t) * t
  else: fly[G](f, e0, t).e

proc windowHit*[G: static FieldGeometry; W: static Bounds](f: Field, tab: MixtureTable, fm: FlightMaxima,
                                           rstep: eV⁻¹, el: Electron, rng: var (WalkerRng | MagRng),
                                           nnull: var int32): Hit =
  ## gas at 0 K: next real collision, exact bounds per window
  let e0 = el
  let ap = (when G in {fgElectric, fgParallel}: e0.dir.z * f.f2 * sqrt(e0.e)
            else: 0.0.eV•ps⁻¹)
  var tw = 0.0.ps    # start of the window
  var ew = e0.e      # energy there
  var depth = -lnSingle(drand48(rng))
  var real = false
  while not real:
    const wc = windowCollisions[W]
    var tau = wc / tab.tcf[binOf(ew, rstep)]
    template boundOver(tau: ps): ps⁻¹ =
      when G in {fgElectric, fgParallel}:
        let (elo, ehi) = energyRange(e0.e, ap, f.bp, tw, tw + tau)
      else:
        let (elo, ehi) = energyRange(ew, f.slope, tau)
      rangeMax(fm, 1, binOf(elo, rstep), binOf(ehi, rstep))
    var bound = boundOver(tau)
    if bound * tau > shrinkAbove[W] * wc:
      # shorter, by the bound: as many trials, less of the range
      tau = wc / bound
      bound = boundOver(tau)
    let reach = depth / bound
    if reach > tau:
      depth = depth - bound * tau
      tw = tw + tau
      ew = flightEnergy[G](f, e0, ap, tw)
    else:
      tw = tw + reach
      when G in {fgElectric, fgParallel}:
        result.fl.t = tw
        result.fl.e = e0.e + (ap + f.bp * tw) * tw
      else:
        result.fl = fly[G](f, e0, tw)
      ew = result.fl.e
      result.ie = binOf(result.fl.e, rstep)
      real = drand48(rng) * bound <= tab.tcf[result.ie]
      if not real:
        inc nnull
        depth = -lnSingle(drand48(rng))
  finishFlight[G](f, result.fl)
  result.e = result.fl.e
  result.din = directionAt[G](f, e0, result.fl)
  result.gas = 1

proc windowHit*[G: static FieldGeometry; W: static Bounds](f: Field, tabs: GasTables, fm: FlightMaxima,
                                           rstep: eV⁻¹, el: Electron, target: var WarmGas,
                                           rng: var (WalkerRng | MagRng), nnull, beyond: var int32): Hit =
  ## thermal motion: next real collision, exact bounds per window; energy
  ## and direction in the rest frame of the molecule hit
  let e0 = el
  let n = int(target.gas.count)
  let ap = (when G in {fgElectric, fgParallel}: e0.dir.z * f.f2 * sqrt(e0.e)
            else: 0.0.eV•ps⁻¹)
  let ivfac2 = 1.0 / f.vfac2
  when G == fgElectric:
    let v0 = e0.dir * (f.vfac * sqrt(e0.e))
    let az = f.efield * f.const5 * f.vfac
  var bounds: array[1..maxGases, ps⁻¹]
  var tw = 0.0.ps
  var ew = e0.e
  var depth = -lnSingle(drand48(rng))
  var vrel: Vector[m•ps⁻¹]
  var real = false
  while not real:
    var nu = 0.0.ps⁻¹
    for k in 1 .. n: nu = nu + tabs[k].tcf[binOf(ew, rstep)]
    const wc = windowCollisions[W]
    var tau = wc / nu
    var total = 0.0.ps⁻¹
    template boundsOver(tau: ps) =
      when G in {fgElectric, fgParallel}:
        let (elo, ehi) = energyRange(e0.e, ap, f.bp, tw, tw + tau)
      else:
        let (elo, ehi) = energyRange(ew, f.slope, tau)
      total = 0.0.ps⁻¹
      for k in 1 .. n:
        # in the frame of a molecule slower than the margin
        let u = speedMargin * target.gas.vtmb[k] / f.vfac
        let slo = max(sqrt(elo) - u, 0.0.eV¹ᐟ²)
        let shi = sqrt(ehi) + u
        bounds[k] = rangeMax(fm, k, binOf(slo * slo, rstep), binOf(shi * shi, rstep))
        total = total + bounds[k]
    boundsOver(tau)
    if total * tau > shrinkAbove[W] * wc:
      # shorter, by the bound: as many trials, less of the range
      tau = wc / total
      boundsOver(tau)
    let reach = depth / total
    if reach > tau:
      depth = depth - total * tau
      tw = tw + tau
      ew = flightEnergy[G](f, e0, ap, tw)
    else:
      tw = tw + reach
      when G in {fgElectric, fgParallel}:
        result.fl.t = tw
        result.fl.e = e0.e + (ap + f.bp * tw) * tw
      else:
        result.fl = fly[G](f, e0, tw)
      ew = result.fl.e
      when G == fgElectric:
        var vel = v0
        vel.z = vel.z + az * tw
      else:
        finishFlight[G](f, result.fl)
        let vel = velocityAt[G](f, e0, result.fl)
      let rg = drand48(rng) * total
      var k = 1
      var below = bounds[1]
      while below < rg and k < n:
        k = k + 1
        below = below + bounds[k]
      result.gas = k
      result.vgas = target.maxwell.velocity(target.gas.vtmb[k], rng)
      vrel = vel - result.vgas
      result.e = norm2(vrel) * ivfac2
      result.ie = binOf(result.e, rstep)
      let tcf = tabs[k].tcf[result.ie]
      if tcf > bounds[k]:
        inc beyond
        real = true
      else:
        real = drand48(rng) * bounds[k] <= tcf
        if not real:
          inc nnull
          depth = -lnSingle(drand48(rng))
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
  w.beyond = 0
  w.tmax = 0.0.ps

proc walkWith*[G: static FieldGeometry; T, C, R](w: var Walker[T, R], f: Field,
                                                          tabs: TablesOf(T),
                                              fm: ptr FlightMaxima, bounds: static Bounds,
                                              estep: eV, collisions, penning: int32,
                                              tally: C) =
  ## `walk` with the exact bounds per window (from `fm`, window lengths of
  ## `bounds`) or Magboltz's bounds (`fm` unused). GPU kernels call it
  ## directly (`gpuWindows`): no null pointers in WGSL.
  mixin count
  let rstep = 1.0 / estep
  for _ in 0'i32 ..< collisions:
    let el = w.el  # a copy, so its square roots are computed once
    var h: Hit
    when bounds != magboltzBounds:
      when T is ColdGas:
        h = windowHit[G, bounds](f, tabs, fm[], rstep, el, w.rng, w.nulls)
      else:
        h = windowHit[G, bounds](f, tabs, fm[], rstep, el, w.target, w.rng, w.nulls, w.beyond)
    else:
      h = nextHit[G](f, tabs, w.target, rstep, el, w.rng, w.nulls)
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
    let s = scatterFast(tab, i, h.ie, h.e, ei, h.din, w.rng)
    w.el = leave[G](w.target, f, s, h)
  w.collisions = w.collisions + collisions

proc walk*[G: static FieldGeometry; T, C, R](w: var Walker[T, R], f: Field,
                                                      tabs: TablesOf(T),
                                          fm: ptr FlightMaxima, estep: eV,
                                          collisions, penning: int32, tally: C) =
  ## `collisions` real collisions (angles by `scatterFast`); `penning` ≠ 0:
  ## Penning transfers; each collision is counted in `tally` (`Tally`,
  ## `NoTally` or a GPU's).
  ## `fm`: exact bounds per window; nil: Magboltz's bounds.
  if fm != nil:
    walkWith[G, T, C, R](w, f, tabs, fm, cpuWindows, estep, collisions, penning, tally)
  else:
    walkWith[G, T, C, R](w, f, tabs, fm, magboltzBounds, estep, collisions, penning, tally)

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
    Job[T, R] = object
      ws: ptr UncheckedArray[Walker[T, R]]
      n: int
      f: ptr Field
      tab: ptr TablesOf(T)
      fm: ptr FlightMaxima
      estep: eV
      collisions, penning: int32
      tallies: ptr Tallies  # nil: count nothing

  proc walkJob[G: static FieldGeometry; T, R](j: ptr Job[T, R]) {.thread.} =
    for k in 0 ..< j.n:
      if j.tallies == nil:
        walk[G, T, NoTally, R](j.ws[k], j.f[], j.tab[], j.fm, j.estep, j.collisions, j.penning,
                            NoTally())
      else:
        walk[G, T, Tally, R](j.ws[k], j.f[], j.tab[], j.fm, j.estep, j.collisions, j.penning,
                          tally(j.tallies[]))

proc walkAll*[G: static FieldGeometry; T, R](ws: var openArray[Walker[T, R]], f: Field,
                                          tab: TablesOf(T), estep: eV, collisions,
                                          penning: int32, threads = 1,
                                          tallies: ptr Tallies = nil,
                                          fm: ptr FlightMaxima = nil) =
  ## `collisions` real collisions of every walker, on `threads` CPU threads;
  ## counted in `tallies` unless nil
  when compileOption("threads"):
    if threads > 1 and ws.len > 1:
      let n = min(threads, ws.len)
      var jobs = newSeq[Job[T, R]](n)
      var ts = newSeq[Thread[ptr Job[T, R]]](n)
      var parts = newSeq[Tallies](if tallies == nil: 0 else: n)
      for k in 0 ..< n:
        let a = k * ws.len div n
        let b = (k + 1) * ws.len div n
        jobs[k] = Job[T, R](ws: cast[ptr UncheckedArray[Walker[T, R]]](addr ws[a]), n: b - a,
                      f: unsafeAddr f, tab: unsafeAddr tab, fm: fm, estep: estep,
                      collisions: collisions, penning: penning,
                      tallies: if tallies == nil: nil else: addr parts[k])
        createThread(ts[k], walkJob[G, T, R], addr jobs[k])
      joinThreads(ts)
      for p in parts: tallies[].add p
      return
  for w in ws.mitems:
    if tallies == nil:
      walk[G, T, NoTally, R](w, f, tab, fm, estep, collisions, penning, NoTally())
    else: walk[G, T, Tally, R](w, f, tab, fm, estep, collisions, penning, tally(tallies[]))

type
  CpuWalkers*[G: static FieldGeometry; T; R = Ranmar] = object
    ## walkers on CPU threads, with the interface of `gpuwalk.GpuWalkers`
    ws: seq[Walker[T, R]]
    f: Field
    tab: ptr TablesOf(T)
    fm: ptr FlightMaxima
    estep: eV
    penning: int32
    threads: int
    tallies: ref Tallies

proc initCpuWalkers*[G: static FieldGeometry; T, R](ws: openArray[Walker[T, R]],
                                                            f: Field,
                                                 tab: ptr TablesOf(T), estep: eV,
                                                 penning = false, threads = 1,
                                                 fm: ptr FlightMaxima = nil): CpuWalkers[G, T, R] =
  ## `tab` and `fm` must outlive the walkers
  result = CpuWalkers[G, T, R](ws: @ws, f: f, tab: tab, fm: fm, estep: estep,
                            penning: int32(ord(penning)), threads: threads)
  new result.tallies

proc walk*(c: var CpuWalkers, collisions: int32, counting = true) =
  ## `collisions` real collisions of every walker, counted in the tally
  ## unless not `counting`
  walkAll[c.G, c.T, c.R](c.ws, c.f, c.tab[], c.estep, collisions, c.penning, c.threads,
                    if counting: addr c.tallies[] else: nil, c.fm)

proc read*[T, R](c: CpuWalkers, ws: var openArray[Walker[T, R]]) =
  for k in 0 ..< ws.len: ws[k] = c.ws[k]

proc write*[T, R](c: var CpuWalkers, ws: openArray[Walker[T, R]]) =
  for k in 0 ..< ws.len: c.ws[k] = ws[k]

func dropsWalkers*(c: CpuWalkers): bool = true
  ## time goes with all collisions: fewer walkers, less time

proc shrink*(c: var CpuWalkers, n: int) =
  ## keeps the first n walkers
  c.ws.setLen(n)

proc clearTally*(c: var CpuWalkers) = c.tallies[] = Tallies()

proc readTally*(c: CpuWalkers, t: var Tallies) = t = c.tallies[]
