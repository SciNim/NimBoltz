## Electron avalanches sampled at planes in time or in space: the runs behind
## Magboltz's pulsed-Townsend, time-of-flight and steady-state Townsend
## (SST) results (MONTEFT, MONTEFTA, MONTEFTG, MONTEFTH, MONTEFD and their
## thermal versions MONTEFTT, MONTEFTAT, MONTEFTGT, MONTEFTHT, MONTEFDT).
##
## Primaries start one after another at the origin, each from the state of
## the previous avalanche after a multiple of 200 collisions. Ionisation
## adds the secondary electron (and Auger electrons) to a stack, Penning
## transfer adds an electron of 1 eV; after an electron ends, the next one
## is taken from the stack. Electrons are sampled where they cross
##
## - time planes (`pkTime`): at the times k·TSTEP after the primary started;
##   an electron ends after the last plane (TFINAL),
## - space planes (`pkSpace`, B = 0 only): at the distances k·ZSTEP along E,
##   weighted by the time spent there (1/|vz|); an electron ends after the
##   last plane or when it can no longer reach the next one.
##
## Electrons also end when attached. The run stops at the first primary
## that starts after NMAX real collisions.
##
## Unlike MONTE, the null-collision method uses one bound for all energies
## (`ColdBound`, gas at 0 K) or the bounds of the thermal gas (`WarmGas`).
##
## Formulas and their evaluation order are those of the Fortran, so results
## agree bit for bit, including where the ten routines differ from each
## other (marked "as MONTEF..."). Units as in Magboltz: positions in m,
## times in ps, energies in eV.
##
## With `FlightMaxima` (`initCascade`'s `fm`), the flights take the exact
## bounds per time window of the walkers (`walkers.windowHit`), and the
## angles `scatterFast`: the same distributions, far fewer null collisions,
## not Magboltz's numbers.

import std/math
import units, rm48, xoshiro, tables, vec3, fields, collisions, thermal, freeflight, drift
from walkers import FlightMaxima, windowHit, cpuWindows

const
  nPlanes* = 8               ## planes of either kind
  stackSize* = 8000          ## electrons waiting to be followed (ISPT)
  maxPrimaries* = 10000000

type
  PlaneKind* = enum
    pkTime   ## planes at fixed times (TPLANE)
    pkSpace  ## planes at fixed distances along E (SPLANE, TCALC); B = 0 only

  ColdBound* = object
    ## gas at 0 K: one upper bound of the collision frequency for all
    ## energies, the largest of the 8 bands of MONTE
    tlim*: ps⁻¹

  Planes* = object
    ## where the electrons are sampled (/CION/, /CTCALC/ from ALPCALC)
    tstep*, tfinal*: ps          ## pkTime: time between planes, last plane
    itfinal*: int                ## pkTime: number of planes
    zstep*, zfinal*: m           ## pkSpace: distance between planes, last plane
    zplane*: array[1..nPlanes, m]
    izfinal*: int                ## pkSpace: number of planes

  TimeSums* = object
    ## sums over the electrons crossing each time plane (TPLOUT, TPLOUTG,
    ## TPLOUTH); fgOblique: positions and velocities in the frame of E
    x*, y*, z*: array[1..nPlanes, m]
    xx*, yy*, zz*, yz*, xz*, xy*: array[1..nPlanes, m²]
    e*: array[1..nPlanes, eV]
    t*: array[1..nPlanes, ps]
    vx*, vy*, vz*: array[1..nPlanes, m•ps⁻¹]
    n*: array[1..nPlanes, int]

  SpaceSums* = object
    ## sums over the electrons crossing each space plane, weighted by
    ## w = 1/|vz| (SPLOUT, SPL1)
    x*, y*, r*, z*: array[1..nPlanes, ps]       ## w·x, ...
    xx*, yy*, rr*, zz*: array[1..nPlanes, m•ps]  ## w·x², ...
    t*: array[1..nPlanes, ps²•m⁻¹]       ## w·t
    tt*: array[1..nPlanes, ps³•m⁻¹]      ## w·t²
    rate*: array[1..nPlanes, m⁻¹]        ## w/t
    e*: array[1..nPlanes, eV•ps•m⁻¹]
    vz*: array[1..nPlanes, float64]      ## w·vz
    w*: array[1..nPlanes, ps•m⁻¹]
    w2*: array[1..nPlanes, ps²•m⁻²]
    born*: array[1..nPlanes + 1, int]
      ## per ZSTEP slice: electrons created minus electrons attached (NESST)

  Pending = object
    ## an electron waiting on the stack
    pos: Vector[m]
    st: ps
    e: eV
    dir: Vec3
    plane: int  ## pkTime: time planes passed; pkSpace: slice of the creation point

  NegativeEnergy* = object
    ## a flight that ended with a negative energy (round-off); Magboltz
    ## prints it and continues with 0.001 eV
    iter*: int
    el*: Electron
    fl*: Flight

  Stop* = enum
    stNone
    stTooManyPrimaries  ## more than maxPrimaries
    stStackFull         ## more than stackSize electrons waiting

  Cascade*[G: static FieldGeometry; P: static PlaneKind; T] = object
    ## The state of a run; T is `ColdBound` or `WarmGas`.
    f: Field
    const2: m•cm•V⁻¹•ps⁻²    ## CONST2, CONST3 of SETUP (plane formulas)
    const3: cm•ps⁻¹•eV⁻¹ᐟ²
    target*: T
    estep: eV
    nmax: int
    penning: bool
    planes: Planes
    rcs, rsn: float64        ## fgOblique: rotation to the frame of E
    # the electron followed
    el: Electron
    pos*: Vector[m]
    st*: ps                  ## time since its primary started
    zstrt: m                 ## where and when it started
    tsstrt: ps
    plane: int               ## pkTime: planes passed; pkSpace: next plane (IZPLANE)
    tstop: ps                ## pkTime: time of the next plane
    tzstop, tzstop1: ps      ## pkSpace: flight times to the next plane (TCALC)
    isol: int                ## pkSpace: 2 if it crosses the plane below first
    stack: seq[Pending]
    next: Electron           ## start of the next primary
    i100: int
    # counts and results
    iprim*: int              ## primaries
    nelec*: int              ## electrons followed
    neion*: int              ## attachments
    nclus*: int              ## electrons created (incl. primaries after the first)
    id*: int                 ## real collisions (ITER)
    iprint, jprint: int
    j1*: int                 ## 1 + number of times NMAX/10 collisions passed
    nnull*: int64
    tmax*: ps                ## longest flight
    times*: array[1..nTimeBins, float64]       ## flights per 1 ps bin
    spectrum*: array[1..nEnergyBins, float64]  ## real collisions per energy bin
    when T is ColdBound:
      tally*: Tally[maxMixProcesses]
    else:
      tally*: GasTallies
    ztot*: m                 ## sums of the end position and time of the electrons
    ttot*: ps
    ztots*: m                ## same, from where each electron started (TTRM)
    ttots*: ps
    primaryE*: eV            ## sums of the start energies of the primaries
    primaryE2*: eV²
    theta*: float64          ## polar angle of the last scattering (THETA)
    when P == pkTime:
      sums*: TimeSums
    else:
      sums*: SpaceSums
    warnings*: seq[NegativeEnergy]
    stop*: Stop
    primaryLimit: int        ## stTooManyPrimaries above (maxPrimaries, or its
                             ## share of a walker)
    fm: ptr FlightMaxima     ## nil: Magboltz's flights and angles; else
                             ## exact bounds per window, `scatterFast`
    cosTheta: float64        ## of the last scattering (with `fm`)

func initCascade*[G: static FieldGeometry; P: static PlaneKind; T](
    f: Field, const2: m•cm•V⁻¹•ps⁻², const3: cm•ps⁻¹•eV⁻¹ᐟ², start: Electron, target: T,
    estep: eV,
    nmax: int, penning: bool, planes: Planes, btheta: float64,
    primaryLimit = maxPrimaries, fm: ptr FlightMaxima = nil): Cascade[G, P, T] =
  ## `start`: the first primary; `primaryLimit`: of this run (a walker gets
  ## its share of maxPrimaries); `fm`: the exact bounds of the walkers (it
  ## must outlive the run), nil: Magboltz's
  static: doAssert P == pkTime or G == fgElectric
  result.f = f
  result.const2 = const2
  result.const3 = const3
  result.target = target
  result.estep = estep
  result.nmax = nmax
  result.penning = penning
  result.planes = planes
  when G == fgOblique:
    let api = arccos(-1.0)
    result.rcs = cos((btheta - 90.0) * api / 180.0)
    result.rsn = sin((btheta - 90.0) * api / 180.0)
  result.el = start
  result.next = start
  result.j1 = 1
  result.jprint = nmax div 10
  result.primaryLimit = primaryLimit
  result.fm = fm

func stacked*[G: static FieldGeometry; P: static PlaneKind; T](c: Cascade[G, P, T]): int =
  c.stack.len

func lastTheta*[G: static FieldGeometry; P: static PlaneKind; T](c: Cascade[G, P, T]): float64 =
  ## polar angle of the last scattering (Magboltz's THETA; `scatterFast`
  ## leaves `theta` out)
  if c.fm == nil: c.theta else: arccos(clamp(c.cosTheta, -1.0, 1.0))

# --- sampling at the planes -------------------------------------------------

proc sampleTime[G: static FieldGeometry; P: static PlaneKind; T](c: var Cascade[G, P, T], k: int) =
  ## TPLANE: the electron crosses time plane k
  let el = c.el
  let tau = float64(k) * c.planes.tstep - c.st
  let t2 = tau * tau
  var p: Vector[m]
  var v: Vector[m•ps⁻¹]
  var e: eV
  when G in {fgElectric, fgParallel}:
    let ap = el.dir.z * c.f.f2 * sqrt(el.e)
    let a = ap * tau
    let b = c.f.bp * t2
    e = el.e + a + b
    let const6 = sqrt(el.e / e)
    let dcz2 = el.dir.z * const6 + c.f.efield * tau * c.f.const5 / sqrt(e)
    when G == fgElectric:
      p.x = c.pos.x + el.dir.x * tau * sqrt(el.e) * c.const3 * unity(m•cm⁻¹)
      p.y = c.pos.y + el.dir.y * tau * sqrt(el.e) * c.const3 * unity(m•cm⁻¹)
    else:
      let wbt = c.f.wb * tau
      let cwt = cos(wbt)
      let swt = sin(wbt)
      p.x = c.pos.x + (el.vel.x * swt - el.vel.y * (1.0 - cwt)) / c.f.wb
      p.y = c.pos.y + (el.vel.y * swt + el.vel.x * (1.0 - cwt)) / c.f.wb
    p.z = c.pos.z + el.dir.z * tau * sqrt(el.e) * c.const3 * unity(m•cm⁻¹) +
          t2 * c.f.efield * c.const2
    v.z = dcz2 * sqrt(e) * c.const3 * unity(m•cm⁻¹)
  else:
    let wbt = c.f.wb * tau
    let cwt = cos(wbt)
    let swt = sin(wbt)
    let eovb = c.f.eovb
    var cx2 = el.vel.x
    when G == fgOblique:
      cx2 = el.vel.x + 2.0 * c.f.f1 * tau
    let cy2 = (el.vel.y - eovb) * cwt + el.vel.z * swt + eovb
    let cz2 = el.vel.z * cwt - (el.vel.y - eovb) * swt
    let vtot = sqrt(cx2 * cx2 + cy2 * cy2 + cz2 * cz2)
    let dcz2 = float64(cz2 / vtot)
    let dcy2 = float64(cy2 / vtot)
    when G == fgCrossed:
      p.x = c.pos.x + el.vel.x * tau
    else:
      let dx = el.vel.x * tau + c.f.f1 * tau * tau
      p.x = c.pos.x + dx
    p.y = c.pos.y + eovb * tau + ((el.vel.y - eovb) * swt + el.vel.z * (1.0 - cwt)) / c.f.wb
    let dz = (el.vel.z * swt + (eovb - el.vel.y) * (1.0 - cwt)) / c.f.wb
    p.z = c.pos.z + dz
    when G == fgCrossed:
      e = el.e + dz * c.f.efield * perElectron  # as Magboltz: (dz·E)·100
    else:
      e = el.e + dz * c.f.ez100 + dx * c.f.ex100
    v.z = dcz2 * sqrt(e) * c.const3 * unity(m•cm⁻¹)
    v.y = dcy2 * sqrt(e) * c.const3 * unity(m•cm⁻¹)
    when G == fgOblique:
      let dcx2 = float64(cx2 / vtot)
      v.x = dcx2 * sqrt(e) * c.const3 * unity(m•cm⁻¹)
      # to the frame of E
      p = vec3(p.z * c.rsn + p.x * c.rcs, p.y, p.z * c.rcs - p.x * c.rsn)
      v = vec3(v.z * c.rsn + v.x * c.rcs, v.y, v.z * c.rcs - v.x * c.rsn)
  template s: untyped = c.sums
  s.x[k] = s.x[k] + p.x
  s.y[k] = s.y[k] + p.y
  s.z[k] = s.z[k] + p.z
  s.xx[k] = s.xx[k] + p.x * p.x
  s.yy[k] = s.yy[k] + p.y * p.y
  s.zz[k] = s.zz[k] + p.z * p.z
  when G in {fgCrossed, fgOblique}:
    s.yz[k] = s.yz[k] + p.y * p.z
  when G == fgOblique:
    s.xz[k] = s.xz[k] + p.x * p.z
    s.xy[k] = s.xy[k] + p.x * p.y
  s.e[k] = s.e[k] + e
  s.t[k] = s.t[k] + c.st + tau
  s.vz[k] = s.vz[k] + v.z
  when G in {fgCrossed, fgOblique}:
    s.vy[k] = s.vy[k] + v.y
  when G == fgOblique:
    s.vx[k] = s.vx[k] + v.x
  inc s.n[k]

proc sampleSpace[G: static FieldGeometry; P: static PlaneKind; T](c: var Cascade[G, P, T],
                                                                  tau: ps, k: int) =
  ## SPLANE: the electron crosses space plane k at the flight time tau
  if k > nPlanes: return
  let el = c.el
  let t2 = tau * tau
  let ap = el.dir.z * c.f.f2 * sqrt(el.e)
  let a = ap * tau
  let b = c.f.bp * t2
  let e = el.e + a + b
  let const6 = sqrt(el.e / e)
  let dcz2 = el.dir.z * const6 + c.f.efield * tau * c.f.const5 / sqrt(e)
  let x = c.pos.x + el.dir.x * tau * sqrt(el.e) * c.const3 * unity(m•cm⁻¹)
  let y = c.pos.y + el.dir.y * tau * sqrt(el.e) * c.const3 * unity(m•cm⁻¹)
  let z = c.pos.z + el.dir.z * tau * sqrt(el.e) * c.const3 * unity(m•cm⁻¹) +
          t2 * c.f.efield * c.const2
  let vz = dcz2 * sqrt(e) * c.const3 * unity(m•cm⁻¹)
  let w = abs(1.0 / vz)
  let r = sqrt(x*x + y*y)
  let t = c.st + tau
  template s: untyped = c.sums
  s.x[k] = s.x[k] + x * w
  s.y[k] = s.y[k] + y * w
  s.r[k] = s.r[k] + r * w
  s.z[k] = s.z[k] + z * w
  s.t[k] = s.t[k] + t * w
  s.tt[k] = s.tt[k] + t * t * w
  s.xx[k] = s.xx[k] + x * x * w
  s.yy[k] = s.yy[k] + y * y * w
  s.rr[k] = s.rr[k] + r * r * w
  s.zz[k] = s.zz[k] + z * z * w
  s.e[k] = s.e[k] + e * w
  s.rate[k] = s.rate[k] + w / t
  s.vz[k] = s.vz[k] + vz * w
  s.w[k] = s.w[k] + w
  s.w2[k] = s.w2[k] + w * w

proc reachNextPlane[G: static FieldGeometry; P: static PlaneKind; T](
    c: var Cascade[G, P, T]): bool =
  ## TCALC: the next space plane and the flight time to it (assuming no
  ## collision); false if the electron cannot reach a plane. With isol = 2
  ## the electron first falls back through the plane below (tzstop) and
  ## crosses it again on the way up (tzstop1).
  c.isol = 1
  let a = c.f.efield * c.const2
  let b = sqrt(c.el.e) * c.const3 * unity(m•cm⁻¹) * c.el.dir.z
  let b2 = b * b
  let z = c.pos.z
  let zp = c.planes.zplane
  var k = 1
  while k <= nPlanes and not (z < zp[k]): inc k
  c.plane = k
  var c1, c2: m
  if k == 1:
    c1 = z - zp[1]
  elif k <= nPlanes:
    c1 = z - zp[k]
    c2 = z - zp[k - 1]
  else:
    c1 = z - zp[nPlanes] - 10.0 * c.planes.zstep
    c2 = z - zp[nPlanes]
  var fac = b2 - 4.0 * a * c1
  if fac < 0.0.m²•ps⁻²:
    c.tzstop = -99.0.ps
    return false
  var tstop1 = (-b + sqrt(b2 - 4.0 * a * c1)) / (2.0 * a)
  var tstop2 = (-b - sqrt(b2 - 4.0 * a * c1)) / (2.0 * a)
  if tstop1 < tstop2:
    c.tzstop = (if tstop1 >= 0.0.ps: tstop1 else: tstop2)
  else:
    c.tzstop = (if tstop2 >= 0.0.ps: tstop2 else: tstop1)
  if k == 1: return true
  fac = b2 - 4.0 * a * c2
  if fac < 0.0.m²•ps⁻²: return true
  tstop1 = (-b + sqrt(fac)) / (2.0 * a)
  tstop2 = (-b - sqrt(fac)) / (2.0 * a)
  if tstop1 < 0.0.ps: return true
  c.isol = 2
  dec c.plane
  if tstop1 < tstop2:
    c.tzstop = tstop1
    c.tzstop1 = tstop2
  else:
    c.tzstop = tstop2
    c.tzstop1 = tstop1
  true

proc stranded[G: static FieldGeometry; P: static PlaneKind; T](c: var Cascade[G, P, T]): bool =
  ## pkSpace: the electron cannot get back to the last plane, or cannot reach
  ## the next one; otherwise prepares the next plane crossing
  if c.pos.z > c.planes.zfinal:
    let epot = c.f.efield * (c.pos.z - c.planes.zfinal) * perElectron  # (E·Δz)·100
    if c.el.e < epot: return true
  not c.reachNextPlane()

func slice(planes: Planes, z: m): int =
  ## ZSTEP slice of z (1..9)
  int(clamp(1.0 + trunc(float64(z / planes.zstep)), 1.0, float64(nPlanes + 1)))

# --- electrons ---------------------------------------------------------------

proc finishElectron[G: static FieldGeometry; P: static PlaneKind; T](
    c: var Cascade[G, P, T]) {.inline.} =
  c.ztot = c.ztot + c.pos.z
  c.ttot = c.ttot + c.st
  c.ztots = c.ztots + c.pos.z - c.zstrt
  c.ttots = c.ttots + c.st - c.tsstrt

proc push[G: static FieldGeometry; P: static PlaneKind; T](c: var Cascade[G, P, T], pos: Vector[m],
                                                           st: ps, e: eV, dir: Vec3,
                                                           plane: int): bool {.inline.} =
  ## adds an electron to the stack; false (and stop) if the stack is full
  if c.stack.len >= stackSize:
    c.stop = stStackFull
    return false
  c.stack.add Pending(pos: pos, st: st, e: e, dir: dir, plane: plane)
  true

proc count[G: static FieldGeometry; P: static PlaneKind; T](c: var Cascade[G, P, T],
                                                           gas, kind, i: int) {.inline.} =
  when T is ColdBound:
    inc c.tally.kinds[kind]
    inc c.tally.processes[i]
  else:
    inc c.tally[gas].kinds[kind]
    inc c.tally[gas].processes[i]

func hitVelocity[G: static FieldGeometry; P: static PlaneKind](
    f: Field, el: Electron, fl: Flight): Vector[m•ps⁻¹] {.inline.} =
  ## lab velocity at the end of the flight (thermal motion)
  when G == fgElectric and P == pkSpace:
    # as MONTEFDT: the direction is renormalised through its z component
    var d = directionAt[G](f, el, fl)
    let ais = (if d.z < 0.0: -1.0 else: 1.0)
    # NIMBOLTZ FIX MB-43: near a turning point rounding makes the argument
    # -1e-16; d.z was NaN and the avalanche could not end
    d.z = ais * sqrt(max(0.0, 1.0 - d.x * d.x - d.y * d.y))
    let vtot = f.vfac * sqrt(fl.e)
    d * vtot
  else:
    velocityAt[G](f, el, fl)

func lateLossClamp[G: static FieldGeometry; P: static PlaneKind; T](ei, e: eV): eV {.inline.} =
  ## energy loss used in the kinematics when the loss (with the secondary's
  ## energy) exceeds the energy; the ten routines differ here
  result = ei
  when T is ColdBound:
    if e < ei:
      result = (when G == fgParallel: e else: 0.0.eV)  # as MONTEFTA / the others
  elif G == fgElectric:
    if e < ei:
      result = (when P == pkTime: e - 0.0001.eV else: 0.0.eV)  # as MONTEFTT / MONTEFDT
  # MONTEFTAT, MONTEFTGT, MONTEFTHT: no clamp

proc follow[G: static FieldGeometry; P: static PlaneKind; T](c: var Cascade[G, P, T],
                                                             tabs: MixtureTable | GasTables,
                                                             rng: var (MagRng | WalkerRng)) =
  ## follows the current electron until it ends (or the run stops)
  let fast = c.fm != nil
  let rstep = 1.0 / c.estep
  var tdash = 0.0.ps
  while true:
    var t, told: ps
    var h: Hit
    if fast:
      # the next real collision at once; the planes before it as below
      var nulls = 0'i32
      when T is ColdBound:
        h = windowHit[G, cpuWindows](c.f, tabs, c.fm[], rstep, c.el, rng, nulls)
      else:
        var beyond = 0'i32
        h = windowHit[G, cpuWindows](c.f, tabs, c.fm[], rstep, c.el, c.target, rng, nulls, beyond)
      c.nnull = c.nnull + nulls
      t = h.fl.t
    else:
      let r1 = drand48(rng)
      t = (when T is ColdBound: -(ln(r1) / c.target.tlim) + tdash
           else: -(ln(r1) / c.target.gas.tcfmx) + tdash)
      told = tdash
    tdash = t
    # planes crossed during the flight
    when P == pkTime:
      while t + c.st >= c.tstop:
        inc c.plane
        c.tstop = c.tstop + c.planes.tstep
        c.sampleTime(c.plane)
        if t + c.st >= c.tstop and c.tstop <= c.planes.tfinal: continue
        if t + c.st >= c.planes.tfinal:
          c.finishElectron()
          return
        break
    else:
      while t >= c.tzstop and told < c.tzstop:
        c.sampleSpace(c.tzstop, c.plane)
        if c.plane >= c.planes.izfinal + 1:
          c.finishElectron()
          return
        if c.isol == 2:
          c.tzstop = c.tzstop1
          c.isol = 1
          continue
        break
    var fl: Flight
    if fast:
      fl = h.fl
    else:
      fl = fly[G](c.f, c.el, t)
      when G == fgCrossed:
        fl.e = c.el.e + fl.dz * c.f.efield * perElectron  # as MONTEFTG: (dz·E)·100
      if fl.e < 0.0.eV:
        c.warnings.add NegativeEnergy(iter: c.id, el: c.el, fl: fl)
        fl.e = 0.001.eV
      # real or null collision
      when T is ColdBound:
        let ie = energyBin(fl.e, c.estep)
        let r5 = drand48(rng)
        if not (r5 <= tabs.tcf[ie] / c.target.tlim):
          inc c.nnull
          continue
        finishFlight[G](c.f, fl)
        h = Hit(fl: fl, ie: ie, e: fl.e, din: directionAt[G](c.f, c.el, fl), gas: 1)
      else:
        finishFlight[G](c.f, fl)
        let vel = hitVelocity[G, P](c.f, c.el, fl)
        h.gas = c.target.gas.pickGas(drand48(rng))
        h.vgas = c.target.maxwell.velocity(c.target.gas.vtmb[h.gas], rng)
        let vrel = vel - h.vgas
        h.e = norm2(vrel) / c.f.vfac2
        h.ie = energyBin(h.e, c.estep)
        let r4 = drand48(rng)
        if not (r4 <= tabs[h.gas].tcf[h.ie] / c.target.gas.tcfmax[h.gas]):
          inc c.nnull
          continue
        h.fl = fl
        h.din = vrel * (1.0 / (c.f.vfac * sqrt(h.e)))
    # a real collision
    template tab: untyped = tableOf(tabs, h)
    if t >= c.tmax: c.tmax = t
    tdash = 0.0.ps
    discard advance[G](c.f, c.el, fl, c.pos)
    c.st = c.st + t
    let it = int(min(trunc(t.float + 1.0), float64(nTimeBins)))  # 1 ps bins
    c.times[it] = c.times[it] + 1.0
    c.spectrum[h.ie] = c.spectrum[h.ie] + 1.0
    let r2 = drand48(rng)
    let i = findProcess(tab, h.ie, r2)
    var ei = tab.ein[i]
    if h.e < ei: ei = h.e - 0.0001.eV
    if tab.moderation[i] == 1'i32: ei = h.e - vibrationMark  # MB-58
    var secondary = -1
    if tab.ipn[i] != 0:
      if tab.ipn[i] == -1:
        # attachment: the electron ends
        inc c.neion
        c.count(h.gas, tab.iarry[i], i)
        inc c.id
        inc c.iprint
        c.times[it] = c.times[it] + 1.0
        c.finishElectron()
        when P == pkSpace:
          dec c.sums.born[c.planes.slice(c.pos.z)]
        return
      # ionisation: the secondary takes part of the energy above threshold
      let eistr = ei
      let r9 = drand48(rng)
      var esec = tab.wpl[i] * tan(r9 * arctan((h.e - ei) / (2.0 * tab.wpl[i])))
      esec = tab.wpl[i] * pow(esec / tab.wpl[i], 0.9524000287055969)
      ei = esec + ei
      let plane = (when P == pkTime: c.plane else: c.planes.slice(c.pos.z))
      inc c.nclus
      if not c.push(c.pos, c.st, esec, Vec3(), plane): return  # direction set below
      secondary = c.stack.high
      when P == pkSpace: inc c.sums.born[plane]
      if eistr > 30.0.eV:
        # Auger electrons, isotropic
        let naug = tab.augerCount[i]
        let eavaug = tab.augerEnergy[i] / float64(naug)
        for _ in 1 .. naug:
          inc c.nclus
          let r3 = drand48(rng)
          let f3 = 1.0 - 2.0 * r3
          let theta0 = arccos(f3)
          let f6 = cos(theta0)
          let f5 = sin(theta0)
          let r4 = drand48(rng)
          let phi0 = 2.0 * PI * r4
          if not c.push(c.pos, c.st, eavaug, vec3(cos(phi0) * f5, sin(phi0) * f5, f6), plane):
            return
          when P == pkSpace: inc c.sums.born[plane]
    c.count(h.gas, tab.iarry[i], i)
    inc c.id
    inc c.iprint
    # Penning transfer: an electron of 1 eV, displaced in space and time
    if c.penning and tab.penningFraction[i] != 0.0:
      let ran = drand48(rng)
      if not (ran > tab.penningFraction[i]):
        inc c.nclus
        if c.stack.len >= stackSize:
          c.stop = stStackFull
          return
        var p = c.pos
        var kept = true
        if tab.penningRange[i].float == 0.0:
          when P == pkSpace:
            if p.z > c.planes.zfinal or p.z < 0.0.m: kept = false
        else:
          var asign = 1.0
          var ran = drand48(rng)
          var ran1 = drand48(rng)
          if ran1 < 0.5: asign = -asign
          p.x = c.pos.x - ln(ran) * tab.penningRange[i] * asign
          ran = drand48(rng)
          ran1 = drand48(rng)
          if ran1 < 0.5: asign = -asign
          p.y = c.pos.y - ln(ran) * tab.penningRange[i] * asign
          ran = drand48(rng)
          ran1 = drand48(rng)
          if ran1 < 0.5: asign = -asign
          p.z = c.pos.z - ln(ran) * tab.penningRange[i] * asign
          when P == pkSpace:
            if p.z < 0.0.m or p.z > c.planes.zfinal: kept = false
        if kept:
          var tpen = c.st
          if tab.penningTime[i].float != 0.0:
            let ran = drand48(rng)
            tpen = c.st - ln(ran) * tab.penningTime[i]
          var plane = 0
          when P == pkTime:
            # planes before tpen; dropped after the last
            var tstop1 = 0.0.ps
            kept = false
            for _ in 1 .. c.planes.itfinal:
              tstop1 = c.planes.tstep + tstop1
              if tpen < tstop1:
                kept = true
                break
              inc plane
          else:
            plane = c.planes.slice(p.z)
          if kept:
            discard c.push(p, tpen, 1.0.eV, c.el.dir, plane)
            when P == pkSpace: inc c.sums.born[plane]
        if not kept: dec c.nclus
    # scattering
    ei = lateLossClamp[G, P, T](ei, h.e)
    let s = (if fast: scatterFast(tab, i, h.ie, h.e, ei, h.din, rng)
             else: scatter(tab, i, h.ie, h.e, ei, h.din, rng))
    c.theta = s.theta
    c.cosTheta = s.cosTheta
    if secondary >= 0:
      # the ionisation secondary: transverse momentum balance, opposite azimuth
      var f5s = s.sinTheta * sqrt(float64(s.e / c.stack[secondary].e))
      if f5s > 1.0: f5s = 1.0
      if fast:
        var f6s = sqrt(1.0 - f5s * f5s)
        if s.cosTheta < 0.0: f6s = -f6s
        c.stack[secondary].dir = rotate(h.din, f6s, f5s, -s.cosPhi, -s.sinPhi)
      else:
        let thsec = arcsin(f5s)
        f5s = sin(thsec)
        var f6s = cos(thsec)
        if s.cosTheta < 0.0: f6s = -f6s
        var phis = s.phi + PI
        if phis > 2.0 * PI: phis = phis - 2.0 * PI  # MB-39
        c.stack[secondary].dir = rotate(h.din, f6s, f5s, cos(phis), sin(phis))
    when T is ColdBound:
      c.el = afterCollision[G](c.f, s.e, s.dir)
    else:
      c.el = leave[G](c.target, c.f, s, h)
    inc c.i100
    if c.i100 == 200:
      c.next = c.el
      c.i100 = 0
    when P == pkSpace:
      if c.stranded():
        c.finishElectron()
        return
    if c.iprint > c.jprint:
      c.iprint = 0
      inc c.j1

proc pop[G: static FieldGeometry; P: static PlaneKind; T](c: var Cascade[G, P, T]): bool =
  ## Takes the next electron from the stack. pkSpace: electrons that cannot
  ## reach a plane end at once; false if that ends the primary's avalanche.
  while true:
    let p = c.stack.pop()
    c.pos = p.pos
    c.st = p.st
    c.el = afterCollision[G](c.f, p.e, p.dir)
    c.plane = p.plane
    c.zstrt = c.pos.z
    c.tsstrt = c.st
    when P == pkSpace:
      if c.stranded():
        inc c.nelec
        c.isol = 1
        c.finishElectron()
        if c.nelec == c.nclus + 1: return false
        continue
    return true

proc run*[G: static FieldGeometry; P: static PlaneKind; T](c: var Cascade[G, P, T],
                                                           tabs: MixtureTable | GasTables,
                                                           rng: var (MagRng | WalkerRng)) =
  ## Primaries until NMAX real collisions. `tabs`: `MixtureTable` for
  ## `ColdBound`, `GasTables` for `WarmGas`.
  while true:
    inc c.iprim
    when P == pkSpace:
      c.plane = 0
      c.tzstop = 1000.0.ps
    if c.iprim > 1:
      if c.id > c.nmax:
        dec c.iprim
        return
      c.pos = Vector[m]()
      c.el = afterCollision[G](c.f, c.next.e, c.next.dir)
      inc c.nclus
      c.st = 0.0.ps
      c.tsstrt = 0.0.ps
      c.zstrt = 0.0.m
      when P == pkTime: c.plane = 0
    if c.iprim > c.primaryLimit:
      c.stop = stTooManyPrimaries
      return
    c.primaryE2 = c.primaryE2 + c.el.e * c.el.e
    c.primaryE = c.primaryE + c.el.e
    while true:
      inc c.nelec
      when P == pkTime:
        c.tstop = c.planes.tstep + float64(c.plane) * c.planes.tstep
      c.follow(tabs, rng)
      if c.stop != stNone: return
      if c.nelec == c.nclus + 1: break
      if not c.pop(): break

# --- walkers ------------------------------------------------------------------

proc add[T: TimeSums | SpaceSums](a: var T, b: T) =
  for x, y in fields(a, b):
    for k in low(x) .. high(x): x[k] = x[k] + y[k]

proc add[N: static int](a: var Tally[N], b: Tally[N]) =
  for k in low(a.kinds) .. high(a.kinds): a.kinds[k] = a.kinds[k] + b.kinds[k]
  for k in low(a.processes) .. high(a.processes):
    a.processes[k] = a.processes[k] + b.processes[k]

proc merge*[G: static FieldGeometry; P: static PlaneKind; T](a: var Cascade[G, P, T],
                                                             b: Cascade[G, P, T]) =
  ## adds the run of another walker (independent random numbers, its share
  ## of NMAX) to a: counts, histograms and plane sums add up; the last
  ## electron (THETA, ST) is b's
  a.iprim = a.iprim + b.iprim
  a.nelec = a.nelec + b.nelec
  a.neion = a.neion + b.neion
  a.nclus = a.nclus + b.nclus
  a.id = a.id + b.id
  a.j1 = a.j1 + b.j1 - 1
  a.nnull = a.nnull + b.nnull
  a.tmax = max(a.tmax, b.tmax)
  for k in 1 .. nTimeBins: a.times[k] = a.times[k] + b.times[k]
  for k in 1 .. nEnergyBins: a.spectrum[k] = a.spectrum[k] + b.spectrum[k]
  when T is ColdBound:
    a.tally.add b.tally
  else:
    for g in 1 .. maxGases: a.tally[g].add b.tally[g]
  a.ztot = a.ztot + b.ztot
  a.ttot = a.ttot + b.ttot
  a.ztots = a.ztots + b.ztots
  a.ttots = a.ttots + b.ttots
  a.primaryE = a.primaryE + b.primaryE
  a.primaryE2 = a.primaryE2 + b.primaryE2
  a.sums.add b.sums
  a.theta = b.theta
  a.cosTheta = b.cosTheta
  a.st = b.st
  a.warnings.add b.warnings
  if a.stop == stNone and b.stop != stNone:
    a.stop = b.stop
    a.stack = b.stack
