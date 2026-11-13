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
from walkers import FlightMaxima, windowHit, Bounds, cpuWindows, gpuWindows

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
    itfinal*: int32              ## pkTime: number of planes
    zstep*, zfinal*: m           ## pkSpace: distance between planes, last plane
    zplane*: array[1..nPlanes, m]
    izfinal*: int32              ## pkSpace: number of planes

  TimeSums* = object
    ## sums over the electrons crossing each time plane (TPLOUT, TPLOUTG,
    ## TPLOUTH); fgOblique: positions and velocities in the frame of E
    x*, y*, z*: array[1..nPlanes, m]
    xx*, yy*, zz*, yz*, xz*, xy*: array[1..nPlanes, m²]
    e*: array[1..nPlanes, eV]
    t*: array[1..nPlanes, ps]
    vx*, vy*, vz*: array[1..nPlanes, m•ps⁻¹]
    n*: array[1..nPlanes, int32]

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
    born*: array[1..nPlanes + 1, int32]
      ## per ZSTEP slice: electrons created minus electrons attached (NESST)

  Pending* = object
    ## an electron waiting on the stack (or in a wave, `runWaves`)
    pos*: Vector[m]
    st*: ps
    e*: eV
    dir*: Vec3
    plane*: int32   ## pkTime: time planes passed; pkSpace: slice of the creation point
    tag*: int32     ## waves: the chain of primaries it descends from
    primary*: int32 ## waves: 1 for a primary
    count*: int32   ## waves: the chain's collisions since its last 200th
    burn*: int32    ## waves: 1 in a chain's burn-in (not counted)
    seed*: Seed     ## waves: of its random numbers

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
    nmax*: int
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

proc sampleTime[G: static FieldGeometry; P: static PlaneKind; C](c: var C, k: int) =
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

proc sampleSpace[G: static FieldGeometry; P: static PlaneKind; C](c: var C, tau: ps, k: int) =
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

proc reachNextPlane[G: static FieldGeometry; P: static PlaneKind; C](c: var C): bool =
  ## TCALC: the next space plane and the flight time to it (assuming no
  ## collision); false if the electron cannot reach a plane. With isol = 2
  ## the electron first falls back through the plane below (tzstop) and
  ## crosses it again on the way up (tzstop1).
  c.isol = 1
  let a = c.f.efield * c.const2
  let b = sqrt(c.el.e) * c.const3 * unity(m•cm⁻¹) * c.el.dir.z
  let b2 = b * b
  let z = c.pos.z
  template zp: untyped = c.planes.zplane  # (no array copy: GPU code)
  var k = 1
  while k <= nPlanes and not (z < zp[k]): inc k
  c.plane = typeof(c.plane)(k)
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

proc stranded[G: static FieldGeometry; P: static PlaneKind; C](c: var C): bool =
  ## pkSpace: the electron cannot get back to the last plane, or cannot reach
  ## the next one; otherwise prepares the next plane crossing
  if c.pos.z > c.planes.zfinal:
    let epot = c.f.efield * (c.pos.z - c.planes.zfinal) * perElectron  # (E·Δz)·100
    if c.el.e < epot: return true
  not reachNextPlane[G, P, typeof(c)](c)

func slice(planes: Planes, z: m): int =
  ## ZSTEP slice of z (1..9)
  int(clamp(1.0 + trunc(float64(z / planes.zstep)), 1.0, float64(nPlanes + 1)))

# --- electrons ---------------------------------------------------------------

proc finishElectron[C](c: var C) {.inline.} =
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
  c.stack.add Pending(pos: pos, st: st, e: e, dir: dir, plane: int32(plane))
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
        sampleTime[G, P, typeof(c)](c, c.plane)
        if t + c.st >= c.tstop and c.tstop <= c.planes.tfinal: continue
        if t + c.st >= c.planes.tfinal:
          c.finishElectron()
          return
        break
    else:
      while t >= c.tzstop and told < c.tzstop:
        sampleSpace[G, P, typeof(c)](c, c.tzstop, c.plane)
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
      if stranded[G, P, typeof(c)](c):
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
    c.plane = int(p.plane)
    c.zstrt = c.pos.z
    c.tsstrt = c.st
    when P == pkSpace:
      if stranded[G, P, typeof(c)](c):
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

# --- one electron, as GPU code -------------------------------------------------
#
# `follow` with the walkers' flights and angles (the default mode), written
# so that a GPU compiler takes it (no `break`, `continue`, `return`,
# sequences or references): its sums and counts go into a `Chase`, plain
# data; the electrons it creates, the collision counts and the spectrum go
# to a sink (`CascadeSink`: the thread's cascade, on a CPU).

type
  Chase*[G: static FieldGeometry; P: static PlaneKind; T] = object
    ## the electron followed, and the sums of all electrons followed
    f: Field
    const2: m•cm•V⁻¹•ps⁻²
    const3: cm•ps⁻¹•eV⁻¹ᐟ²
    target: T
    estep: eV
    penning: int32
    planes: Planes
    rcs, rsn: float64
    el: Electron
    pos: Vector[m]
    st: ps
    zstrt: m
    tsstrt: ps
    plane: int32
    tstop: ps
    tzstop, tzstop1: ps
    isol: int32
    tag: int32               ## its chain
    i100: int32              ## the chain's collisions since its last 200th
    records: int32           ## 200th collisions of this electron
    record: Electron         ## the electron after the last of them
    full*: int32             ## 1: no room for the electrons it creates
    # (32-bit counts: added up and zeroed after each wave)
    nelec, neion, nclus, iprim: int32
    id*: int32
    nnull: int32
    tmax: ps
    ztot: m
    ttot: ps
    ztots: m
    ttots: ps
    primaryE: eV
    primaryE2: eV²
    cosTheta: float64
    sums: (when P == pkTime: TimeSums else: SpaceSums)

  CascadeSink[G: static FieldGeometry; P: static PlaneKind; T] = object
    c: ptr Cascade[G, P, T]

func initChase*[G: static FieldGeometry; P: static PlaneKind; T](
    c: Cascade[G, P, T]): Chase[G, P, T] =
  ## for the electrons of c's run
  Chase[G, P, T](f: c.f, const2: c.const2, const3: c.const3, target: c.target,
                 estep: c.estep, penning: int32(ord(c.penning)), planes: c.planes,
                 rcs: c.rcs, rsn: c.rsn)

proc push(s: CascadeSink, p: Pending): int {.inline.} =
  ## index of the new electron, -1 if there is no room
  if s.c.stack.len >= stackSize: return -1
  s.c.stack.add p
  s.c.stack.high

func isFull(s: CascadeSink): bool {.inline.} = s.c.stack.len >= stackSize
proc setDir(s: CascadeSink, i: int, dir: Vec3) {.inline.} = s.c.stack[i].dir = dir
proc countFlight(s: CascadeSink, it: int) {.inline.} =
  s.c.times[it] = s.c.times[it] + 1.0
proc countBin(s: CascadeSink, ie: int) {.inline.} =
  s.c.spectrum[ie] = s.c.spectrum[ie] + 1.0
proc countProcess(s: CascadeSink, gas, kind, i: int) {.inline.} = s.c[].count(gas, kind, i)

proc absorb*[G: static FieldGeometry; P: static PlaneKind; T](c: var Cascade[G, P, T],
                                                             ch: Chase[G, P, T]) =
  ## the sums of a chase into c
  c.iprim = c.iprim + int(ch.iprim)
  c.nelec = c.nelec + int(ch.nelec)
  c.neion = c.neion + int(ch.neion)
  c.nclus = c.nclus + int(ch.nclus)
  c.id = c.id + int(ch.id)
  c.nnull = c.nnull + int64(ch.nnull)
  c.tmax = max(c.tmax, ch.tmax)
  c.ztot = c.ztot + ch.ztot
  c.ttot = c.ttot + ch.ttot
  c.ztots = c.ztots + ch.ztots
  c.ttots = c.ttots + ch.ttots
  c.primaryE = c.primaryE + ch.primaryE
  c.primaryE2 = c.primaryE2 + ch.primaryE2
  c.sums.add ch.sums
  c.cosTheta = ch.cosTheta
  c.st = ch.st

proc zeroSums*[G: static FieldGeometry; P: static PlaneKind; T](ch: var Chase[G, P, T]) =
  ## after `absorb`
  ch.iprim = 0
  ch.nelec = 0
  ch.neion = 0
  ch.nclus = 0
  ch.id = 0
  ch.nnull = 0
  ch.tmax = 0.0.ps
  ch.ztot = 0.0.m
  ch.ttot = 0.0.ps
  ch.ztots = 0.0.m
  ch.ttots = 0.0.ps
  ch.primaryE = 0.0.eV
  ch.primaryE2 = 0.0.eV²
  ch.sums = default(typeof(ch.sums))

proc startChase*[G: static FieldGeometry; P: static PlaneKind; T](c: var Chase[G, P, T],
                                                                 p: Pending): bool =
  ## the electron p, to be followed; false if it ends at once (it cannot
  ## reach a space plane)
  result = true
  c.pos = p.pos
  c.st = p.st
  c.el = afterCollision[G](c.f, p.e, p.dir)
  c.plane = p.plane
  c.zstrt = c.pos.z
  c.tsstrt = c.st
  c.tag = p.tag
  c.i100 = p.count
  c.records = 0
  c.full = 0
  if p.primary != 0:
    c.iprim = c.iprim + 1
    c.primaryE2 = c.primaryE2 + c.el.e * c.el.e
    c.primaryE = c.primaryE + c.el.e
    when P == pkSpace:
      c.plane = 0
      c.isol = 1
      c.tzstop = 1000.0.ps
  else:
    when P == pkSpace:
      if stranded[G, P, typeof(c)](c):
        c.nelec = c.nelec + 1
        c.isol = 1
        finishElectron(c)
        result = false
  if result:
    c.nelec = c.nelec + 1
    when P == pkTime:
      c.tstop = c.planes.tstep + float64(c.plane) * c.planes.tstep

proc collide[G: static FieldGeometry; P: static PlaneKind; T, R, S](
    c: var Chase[G, P, T], tabs: MixtureTable | GasTables, h: Hit, rng: var R, sink: S): bool =
  ## the real collision h; false if the electron ends
  mixin push, isFull, setDir, countFlight, countBin, countProcess
  template tab: untyped = tableOf(tabs, h)
  var going = true
  let fl = h.fl
  let t = fl.t
  if t >= c.tmax: c.tmax = t
  discard advance[G](c.f, c.el, fl, c.pos)
  c.st = c.st + t
  let it = int(min(trunc(t.float + 1.0), float64(nTimeBins)))  # 1 ps bins
  sink.countFlight(it)
  sink.countBin(h.ie)
  let r2 = drand48(rng)
  let i = findProcess(tab, h.ie, r2)
  var ei = tab.ein[i]
  if h.e < ei: ei = h.e - 0.0001.eV
  if tab.moderation[i] == 1'i32: ei = h.e - vibrationMark  # MB-58
  var secondary = -1
  var esec = 0.0.eV
  if tab.ipn[i] == -1:
    # attachment: the electron ends
    c.neion = c.neion + 1
    sink.countProcess(h.gas, tab.iarry[i], i)
    c.id = c.id + 1
    sink.countFlight(it)
    finishElectron(c)
    when P == pkSpace:
      c.sums.born[c.planes.slice(c.pos.z)] = c.sums.born[c.planes.slice(c.pos.z)] - 1
    going = false
  else:
    let plane = (when P == pkTime: c.plane else: c.planes.slice(c.pos.z))
    if tab.ipn[i] != 0:
      # ionisation: the secondary takes part of the energy above threshold
      let eistr = ei
      let r9 = drand48(rng)
      esec = tab.wpl[i] * tan(r9 * arctan((h.e - ei) / (2.0 * tab.wpl[i])))
      esec = tab.wpl[i] * pow(esec / tab.wpl[i], 0.9524000287055969)
      ei = esec + ei
      c.nclus = c.nclus + 1
      secondary = sink.push(Pending(pos: c.pos, st: c.st, e: esec, plane: int32(plane), tag: c.tag,
                                    seed: seedOf(rng)))  # direction below
      if secondary < 0:
        c.full = 1
        going = false
      else:
        when P == pkSpace: c.sums.born[plane] = c.sums.born[plane] + 1
        if eistr > 30.0.eV:
          # Auger electrons, isotropic
          let naug = tab.augerCount[i]
          let eavaug = tab.augerEnergy[i] / float64(naug)
          for _ in 1 .. naug:
            if going:
              c.nclus = c.nclus + 1
              let r3 = drand48(rng)
              let f3 = 1.0 - 2.0 * r3
              let theta0 = arccos(f3)
              let f6 = cos(theta0)
              let f5 = sin(theta0)
              let r4 = drand48(rng)
              let phi0 = 2.0 * PI * r4
              let k = sink.push(Pending(pos: c.pos, st: c.st, e: eavaug,
                                        dir: vec3(cos(phi0) * f5, sin(phi0) * f5, f6),
                                        plane: int32(plane), tag: c.tag, seed: seedOf(rng)))
              if k < 0:
                c.full = 1
                going = false
              else:
                when P == pkSpace: c.sums.born[plane] = c.sums.born[plane] + 1
    if going:
      sink.countProcess(h.gas, tab.iarry[i], i)
      c.id = c.id + 1
      # Penning transfer: an electron of 1 eV, displaced in space and time
      if c.penning != 0 and tab.penningFraction[i] != 0.0:
        let ran = drand48(rng)
        if not (ran > tab.penningFraction[i]):
          c.nclus = c.nclus + 1
          if sink.isFull:
            c.full = 1
            going = false
          else:
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
                var looking = true
                kept = false
                for _ in 1 .. c.planes.itfinal:
                  if looking:
                    tstop1 = c.planes.tstep + tstop1
                    if tpen < tstop1:
                      kept = true
                      looking = false
                    else:
                      plane = plane + 1
              else:
                plane = c.planes.slice(p.z)
              if kept:
                discard sink.push(Pending(pos: p, st: tpen, e: 1.0.eV, dir: c.el.dir,
                                          plane: int32(plane), tag: c.tag, seed: seedOf(rng)))
                when P == pkSpace: c.sums.born[plane] = c.sums.born[plane] + 1
            if not kept: c.nclus = c.nclus - 1
    if going:
      # scattering
      ei = lateLossClamp[G, P, T](ei, h.e)
      let s = scatterFast(tab, i, h.ie, h.e, ei, h.din, rng)
      c.cosTheta = s.cosTheta
      if secondary >= 0:
        # the ionisation secondary: transverse momentum balance, opposite azimuth
        var f5s = s.sinTheta * sqrt(float64(s.e / esec))
        if f5s > 1.0: f5s = 1.0
        var f6s = sqrt(1.0 - f5s * f5s)
        if s.cosTheta < 0.0: f6s = -f6s
        sink.setDir(secondary, rotate(h.din, f6s, f5s, -s.cosPhi, -s.sinPhi))
      when T is ColdBound:
        c.el = afterCollision[G](c.f, s.e, s.dir)
      else:
        c.el = leave[G](c.target, c.f, s, h)
      c.i100 = c.i100 + 1
      if c.i100 == 200'i32:
        c.i100 = 0
        c.records = c.records + 1
        c.record = c.el
      when P == pkSpace:
        if stranded[G, P, typeof(c)](c):
          finishElectron(c)
          going = false
  going

proc chase*[G: static FieldGeometry; P: static PlaneKind; T, R, S](
    c: var Chase[G, P, T], tabs: MixtureTable | GasTables, fm: ptr FlightMaxima,
    bounds: static Bounds, rng: var R, sink: S) =
  ## follows the electron (`startChase`) to its end, with the walkers'
  ## flights (window lengths of `bounds`)
  let rstep = 1.0 / c.estep
  var going = true
  while going:
    var nulls = 0'i32
    var h: Hit
    when T is ColdBound:
      h = windowHit[G, bounds](c.f, tabs, fm[], rstep, c.el, rng, nulls)
    else:
      var beyond = 0'i32
      h = windowHit[G, bounds](c.f, tabs, fm[], rstep, c.el, c.target, rng, nulls, beyond)
    c.nnull = c.nnull + nulls
    let t = h.fl.t
    # planes crossed during the flight
    when P == pkTime:
      var crossing = t + c.st >= c.tstop
      while crossing:
        c.plane = c.plane + 1
        c.tstop = c.tstop + c.planes.tstep
        sampleTime[G, P, typeof(c)](c, c.plane)
        if not (t + c.st >= c.tstop and c.tstop <= c.planes.tfinal):
          crossing = false
          if t + c.st >= c.planes.tfinal:
            finishElectron(c)
            going = false
    else:
      var crossing = t >= c.tzstop and 0.0.ps < c.tzstop
      while crossing:
        sampleSpace[G, P, typeof(c)](c, c.tzstop, c.plane)
        if c.plane >= c.planes.izfinal + 1:
          finishElectron(c)
          going = false
          crossing = false
        elif c.isol == 2:
          c.tzstop = c.tzstop1
          c.isol = 1
          crossing = t >= c.tzstop and 0.0.ps < c.tzstop
        else:
          crossing = false
    if going:
      going = collide[G, P, T, R, S](c, tabs, h, rng, sink)

# --- waves --------------------------------------------------------------------
#
# Avalanches electron by electron: all electrons waiting (of one primary or
# of several) form a wave; each is followed to its end (`chase`),
# independently of the others, and the electrons it creates form the next
# wave. The planes' sums and the counts are sums over electrons, so the
# order does not matter: the same estimates as one electron after the
# other. A wave's electrons can be followed side by side, on CPU threads
# (`runWaves`) even within one avalanche, which one primary after the
# other cannot (hi_air: one primary of 6·10⁸ collisions). Each electron has
# random numbers of its own, seeded from its parent's: the numbers do not
# depend on who follows which electron.

func primaryOf(start: Electron, tag, burn: int32, seed: Seed): Pending =
  Pending(e: start.e, dir: start.dir, tag: tag, primary: 1, burn: burn, seed: seed)

type
  Outcome* = object
    ## of an electron followed, for its chain
    tag*: int32
    collisions*, counted*: int32  ## all, and those towards the 200th (no attachments)
    records*: int32
    record*: Electron

proc outcomeOf*[G: static FieldGeometry; P: static PlaneKind; T](ch: Chase[G, P, T],
                                                                 p: Pending, id0: int32): Outcome =
  Outcome(tag: p.tag, collisions: ch.id - id0,
          counted: 200'i32 * ch.records + ch.i100 - p.count,
          records: ch.records, record: ch.record)

type
  WaveJob[G: static FieldGeometry; P: static PlaneKind; T, Tabs, R] = object
    c, b: ptr Cascade[G, P, T]  ## counted; burn-in (its sums are dropped)
    cc, cb: ptr Chase[G, P, T]
    tabs: ptr Tabs
    fm: ptr FlightMaxima
    wave: ptr seq[Pending]
    outcomes: ptr seq[Outcome]
    first, step: int
    created: seq[Pending]

proc waveJob[G: static FieldGeometry; P: static PlaneKind; T, Tabs, R](
    j: ptr WaveJob[G, P, T, Tabs, R]) {.thread.} =
  ## the electrons first, first + step, ... of the wave
  var k = j.first
  while k < j.wave[].len:
    let p = j.wave[][k]
    let c = (if p.burn != 0: j.b else: j.c)
    let ch = (if p.burn != 0: j.cb else: j.cc)
    var rng = initRng(R, p.seed)
    let id0 = ch.id
    if startChase(ch[], p):
      chase[G, P, T, R, CascadeSink[G, P, T]](ch[], j.tabs[], j.fm, cpuWindows, rng,
                                              CascadeSink[G, P, T](c: c))
    j.outcomes[][k] = outcomeOf(ch[], p, id0)
    if ch.full != 0:
      c.stop = stStackFull
      return
    for q in c.stack.mitems: q.burn = p.burn
    j.created.add c.stack
    c.stack.setLen 0
    k = k + j.step

type
  CpuWaves*[G: static FieldGeometry; P: static PlaneKind; T, Tabs, R] = object
    ## waves on CPU threads: one cascade (the sink) and one chase per thread,
    ## and the same for the burn-in (dropped)
    cs, bs: seq[Cascade[G, P, T]]
    ccs, cbs: seq[Chase[G, P, T]]
    tabs: ptr Tabs

proc initCpuWaves*[G: static FieldGeometry; P: static PlaneKind; T, Tabs, R](
    cs: var seq[Cascade[G, P, T]], tabs: Tabs, rng: typedesc[R]): CpuWaves[G, P, T, Tabs, R] =
  ## on `cs.len` threads
  result.cs = move cs
  result.bs = result.cs
  for c in result.cs:
    result.ccs.add initChase(c)
    result.cbs.add initChase(c)
  result.tabs = unsafeAddr tabs

func width*(b: CpuWaves): int = b.cs.len
  ## electrons followed side by side

proc cpuCollisions[G: static FieldGeometry; P: static PlaneKind; T, Tabs, R](
    b: CpuWaves[G, P, T, Tabs, R]): int =
  ## counted real collisions so far
  for c in b.cs: result = result + c.id

proc cpuStopped[G: static FieldGeometry; P: static PlaneKind; T, Tabs, R](
    b: CpuWaves[G, P, T, Tabs, R]): bool =
  for c in b.cs:
    if c.stop != stNone: return true
  for c in b.bs:
    if c.stop != stNone: return true

proc cpuProcess[G: static FieldGeometry; P: static PlaneKind; T, Tabs, R](
    b: var CpuWaves[G, P, T, Tabs, R], wave, next: var seq[Pending], outcomes: var seq[Outcome]) =
  ## follows the wave's electrons: the electrons they create in `next`,
  ## an outcome per electron
  let n = b.cs.len
  outcomes.setLen wave.len
  var jobs = newSeq[WaveJob[G, P, T, Tabs, R]](n)
  for k in 0 ..< n:
    jobs[k] = WaveJob[G, P, T, Tabs, R](c: addr b.cs[k], b: addr b.bs[k], cc: addr b.ccs[k],
                                        cb: addr b.cbs[k], tabs: b.tabs, fm: b.cs[0].fm,
                                        wave: addr wave, outcomes: addr outcomes,
                                        first: k, step: n)
  if n == 1 or wave.len == 1:
    jobs[0].step = 1
    waveJob(addr jobs[0])
  else:
    when compileOption("threads"):
      var threads = newSeq[Thread[ptr WaveJob[G, P, T, Tabs, R]]](n)
      for k in 0 ..< n:
        createThread(threads[k], waveJob[G, P, T, Tabs, R], addr jobs[k])
      joinThreads(threads)
    else:
      jobs[0].step = 1
      waveJob(addr jobs[0])
  next.setLen 0
  for j in jobs.mitems: next.add move(j.created)
  for k in 0 ..< n:
    b.cs[k].absorb(b.ccs[k])
    b.ccs[k].zeroSums()
    b.cbs[k].zeroSums()

proc finish*[G: static FieldGeometry; P: static PlaneKind; T, Tabs, R](
    b: var CpuWaves[G, P, T, Tabs, R]): Cascade[G, P, T] =
  ## the run: the threads' cascades added up
  result = move b.cs[0]
  for k in 1 ..< b.cs.len: result.merge(b.cs[k])

const
  waveTarget* {.intdefine.} = 4096
    ## chains of primaries per electron followed side by side, at most
  minFinished = 16
    ## avalanches ended before new chains (the first to end are the small
    ## ones: with few, their mean says little)
  burnIn = 2000
    ## collisions of a new chain before it counts (ten 200th collisions:
    ## its start, the electron of another chain, is forgotten)

type
  Chain = object
    ## primaries one after the other, as `run` has them: the next starts
    ## when the avalanche of the last has ended, from the electron of the
    ## chain's last 200th collision
    live: int        ## electrons waiting
    collisions: int  ## of the avalanche going
    count: int       ## collisions since the last 200th
    next: Electron
    recorded: bool   ## a 200th collision in the avalanche going: its start
                     ## has had its 200 collisions
    burnLeft: int    ## collisions until it counts
    burning: bool    ## the avalanche going does not count

proc driveWaves*[B](b: var B, nmax: int, start: Electron, seed: var MagRng,
                    process: proc (b: var B, wave, next: var seq[Pending],
                                   outcomes: var seq[Outcome]) {.nimcall.},
                    collisions: proc (b: B): int {.nimcall.},
                    stopped: proc (b: B): bool {.nimcall.}, width: int): (int, Stop) =
  ## Primaries until NMAX real collisions, as `run`, in waves on the
  ## backend b (`process`: CPU threads, `CpuWaves`, or a GPU), random
  ## numbers seeded from `seed`; returns the primaries and how it stopped.
  ##
  ## The primaries come in chains, each as `run`'s: the next primary when
  ## the avalanche of the last has ended, from the electron of the chain's
  ## last 200th collision (an electron used for the primaries of 200
  ## collisions: the primaries' starts are weighted by collisions). Many
  ## chains go side by side. A new chain starts from the first's last start
  ## and runs `burnIn` collisions that do not count (only where avalanches
  ## are small), so that its starts are the gas's; new chains come when
  ## `minFinished` or more finished avalanches leave room in NMAX (those
  ## still going counted at their mean or more), at most as many as have
  ## finished, and NMAX / (10 `burnIn`) in all. A chain goes on while there
  ## is room, or alone as long as there are fewer than NMAX collisions (as
  ## `run`), and in any case until its last start has had its 200
  ## collisions (else long avalanches' starts would weigh more): about NMAX
  ## collisions in all, as `run` has.
  template newSeed(): Seed = seedOf(seed)
  var chains = @[Chain(live: 1, next: start)]
  var wave = @[primaryOf(start, 0, 0, newSeed())]
  var next: seq[Pending]
  var outcomes: seq[Outcome]
  var finished = 0        # counted avalanches that ended
  var done = 0            # their collisions
  var started = 1
  var stop = stNone
  var waves = 0
  while wave.len > 0 and stop == stNone:
    when defined(debugWaves):
      inc waves
      if waves mod 1000 == 1 or wave.len > 100000:
        stderr.writeLine "wave ", waves, " size ", wave.len, " chains ", chains.len,
                         " id ", collisions(b), " nmax ", nmax
    for p in wave.mitems: p.count = int32(chains[p.tag].count)
    process(b, wave, next, outcomes)
    if stopped(b):
      stop = stStackFull
      break
    # the chains' next starts: the electron of the last 200th collision (by
    # the wave's order)
    for k in 0 ..< wave.len:
      let o = outcomes[k]
      if o.records > 0:
        chains[o.tag].recorded = true
        chains[o.tag].next = o.record
      if chains[o.tag].burning:
        chains[o.tag].burnLeft = chains[o.tag].burnLeft - int(o.collisions)
      else:
        chains[o.tag].collisions = chains[o.tag].collisions + int(o.collisions)
      chains[o.tag].count = (chains[o.tag].count + int(o.counted)) mod 200
    for p in wave: dec chains[p.tag].live
    for q in next: inc chains[q.tag].live
    let id = collisions(b)
    var idle: seq[int32]
    var running = 0
    for t, ch in chains.mpairs:
      if ch.live == 0:
        if not ch.burning and (ch.collisions > 0 or t == 0 and finished == 0):
          inc finished
          done = done + ch.collisions
          ch.collisions = 0
        idle.add int32(t)
      elif not ch.burning:
        inc running
    var room = 0
    let mean = max(float64(done) / float64(max(finished, 1)), 1.0)
    if id <= nmax and finished > 0:
      var expected = float64(done)
      for ch in chains:
        if ch.live > 0 and not ch.burning:
          expected = expected + max(float64(ch.collisions), mean)
      room = int((float64(nmax) - expected) / mean)
    template startPrimary(t: int32) =
      if started >= maxPrimaries:
        stop = stTooManyPrimaries
        break
      chains[t].burning = chains[t].burnLeft > 0
      next.add primaryOf(chains[t].next, t, int32(ord(chains[t].burning)), newSeed())
      chains[t].live = 1
      chains[t].recorded = false
      if not chains[t].burning:
        inc started
        inc running
        room = room - 1
    for t in idle:
      if chains[t].burnLeft > 0 and id <= nmax or not chains[t].recorded or
         id <= nmax and (room > 0 or running == 0):
        startPrimary(t)
    if id <= nmax and stop == stNone:
      # new chains (burn-in only for small avalanches: a long one has 200th
      # collisions of its own, its start counts once anyway)
      let most = min(waveTarget * width, nmax div (10 * burnIn))
      var k = min(min(room, finished), most - chains.len)
      if finished < minFinished: k = 0
      while k > 0:
        chains.add Chain(next: chains[0].next, burnLeft: (if mean < burnIn: burnIn else: 0))
        startPrimary(int32(chains.high))
        dec k
    swap(wave, next)
  (started, stop)

proc runWaves*[G: static FieldGeometry; P: static PlaneKind; T, R](
    cs: var seq[Cascade[G, P, T]], tabs: MixtureTable | GasTables, start: Electron,
    seed: var MagRng): Cascade[G, P, T] =
  ## `driveWaves` on `cs.len` CPU threads (one cascade each, all
  ## initialised alike), random numbers R
  let nmax = cs[0].nmax
  var b = initCpuWaves[G, P, T, typeof(tabs), R](cs, tabs, R)
  let (_, stop) = driveWaves(b, nmax, start, seed, cpuProcess[G, P, T, typeof(tabs), R],
                             cpuCollisions[G, P, T, typeof(tabs), R],
                             cpuStopped[G, P, T, typeof(tabs), R], b.width)
  result = b.finish()
  if stop == stTooManyPrimaries: result.stop = stop
  result.nclus = result.nclus + result.iprim - 1  # primaries after the first, as `run`
  result.j1 = 1 + result.id div max(nmax div 10, 1)
