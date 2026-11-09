## Steady-state drift: drift velocity, diffusion tensor, mean energy and
## Townsend coefficients from one long electron trajectory (MONTE, MONTEA,
## MONTEB, MONTEC and their thermal versions MONTET, MONTEAT, MONTEBT,
## MONTECT).
##
## The trajectory runs in 10 segments of NMAX/10 real collisions. Each
## segment gives one estimate of every coefficient; their scatter gives the
## errors. Diffusion is estimated from the spread of the displacement
## between pairs of collisions NCORLN apart (NCORST pairs per collision),
## which needs the positions of the last NCOLM collisions (`Ring`). Only
## B = 0 and B ∥ E estimate the transverse diffusion without pairs when all
## scattering is isotropic.
##
## Formulas and their evaluation order are those of the Fortran, so results
## agree bit for bit. Units as in Magboltz: positions in m, times in ps
## while running; results in cm/s, cm²/s, eV, 1/cm and percent.

import std/math
import units, rm48, tables, vec3, fields, collisions, freeflight

const
  nSegments* = 10      ## ITMAX
  nTimeBins* = 300     ## histogram of flight times, 1 ps bins
  ringSize* = 2000000  ## size of the position ring of Magboltz (NCOLM <= ringSize)

type
  Decorrelation* = object
    ## /DECOR/: diffusion from pairs of collisions
    ncolm*: int    ## collisions kept in the ring
    ncorln*: int   ## distance of the pairs (collisions)
    ncorst*: int   ## pairs per collision

  Mark = object
    pos: Vector[m]
    st: ps

  Ring* = object
    ## positions and times of the last NCOLM collisions (XST, YST, ZST, STO)
    marks: seq[Mark]

  Tensor*[T] = object
    xx*, yy*, zz*, yz*, xy*, xz*: T

  Tally*[N: static int] = object
    ## real collisions of one collision table
    kinds*: array[1..30, int]     ## per process kind IARRY (ICOLL)
    processes*: array[1..N, int]  ## per process (ICOLN)

  GasTallies* = array[1..maxGases, Tally[maxGasProcesses]]

  Sums = object
    ## running sums of the estimators (and their values at a segment start)
    pos: Vector[m]     ## position (OUTPT X, Y, Z)
    st: ps             ## time
    st1, st2: ps       ## time weights of the pair sums (with / without the
                       ## first two segments)
    s: Tensor[m²]      ## pair sums of the squared displacements
    ls, ts: m²         ## pair sums along / across the drift (fgCrossed)
    vx, vy: m²         ## sums for isotropic scattering
    e2: eV•ps          ## ∫ e dt (fgElectric, fgParallel)
    ebar: eV•ps        ## ∫ e dt from the spectrum (other geometries)
    rot: Vector[m]     ## position rotated to the frame of E (fgOblique)
    srot: Tensor[m²]   ## pair sums rotated to the frame of E (fgOblique)

  Drift*[G: static FieldGeometry; T] = object
    ## The state of one trajectory; T is `ColdGas` or `WarmGas`.
    f: Field
    el: Electron
    target*: T
    estep: eV
    gases: int
    penning: bool       ## Penning transfers (IPEN)
    isotropic: bool     ## all scattering isotropic (NISO = 0)
    dec: Decorrelation
    rcs, rsn: float64   ## fgOblique: rotation to the frame of E
    theta*: float64     ## polar angle of the last scattering (THETA)
    nnull*: int64       ## null collisions
    extraIons*: int     ## ionisations by Auger electrons and Penning transfer
    tmax*: ps           ## longest flight
    times*: array[1..nTimeBins, float64]     ## flights per 1 ps bin
    spectrum*: array[1..nEnergyBins, float64]  ## real collisions per energy bin
    when T is ColdGas:
      tally*: Tally[maxMixProcesses]
    else:
      tally*: GasTallies
    ncol, id: int       ## position in the ring, times it was filled
    w: Vector[m•ps⁻¹]   ## pos / st
    cur, prev: Sums
    # estimates (values of the last segment, as Magboltz keeps them)
    vel*: Vector[μm•ns⁻¹]  ## drift velocity
    velRot: Vector[μm•ns⁻¹]  ## fgOblique: drift velocity in the frame of E
    ave*: eV            ## mean energy
    dif*: Tensor[cm²•s⁻¹]  ## diffusion tensor (fgOblique: lab frame)
    difRot: Tensor[cm²•s⁻¹]  ## fgOblique: diffusion tensor in the frame of E
    difln*, diftr*: cm²•s⁻¹
    # per-segment estimates
    wst: array[1..nSegments, Vector[μm•ns⁻¹]]
    avest: array[1..nSegments, eV]
    dst: array[1..nSegments, Tensor[cm²•s⁻¹]]
    dlnst, dtrst: array[1..nSegments, cm²•s⁻¹]

  Transport* = object
    ## results of a trajectory (what MONTE* leave in /VEL/, /DIFLAB/, ...)
    vel*: Vector[cm•s⁻¹]       ## drift velocity
    velErr*: Vec3              ## and its error (%)
    dif*: Tensor[cm²•s⁻¹]      ## diffusion tensor
    difErr*: Tensor[float64]   ## and its errors (%)
    difln*, diftr*: cm²•s⁻¹    ## longitudinal and transverse diffusion
    dfler*, dfter*: float64    ## their errors (%)
    ave*: eV                   ## mean energy
    den*: float64              ## and its error (%)
    alpha*, att*: cm⁻¹         ## ionisation and attachment coefficients
    alper*, atter*: float64    ## their errors (%)

func sq[T](x: T): auto {.inline.} = x * x

func initRing*(): Ring =
  Ring(marks: newSeq[Mark](ringSize))

func isEmpty*(r: Ring): bool = r.marks.len == 0

func initDrift*[G: static FieldGeometry; T](f: Field, start: Electron, target: T,
                                            estep: eV, gases: int, penning,
                                            isotropic: bool, dec: Decorrelation,
                                            btheta: float64): Drift[G, T] =
  doAssert dec.ncolm <= ringSize
  result.f = f
  result.el = start
  result.target = target
  result.estep = estep
  result.gases = gases
  result.penning = penning
  result.isotropic = isotropic
  result.dec = dec
  when G == fgOblique:
    let api = arccos(-1.0)
    result.rcs = cos((btheta - 90.0) * api / 180.0)
    result.rsn = sin((btheta - 90.0) * api / 180.0)

proc record[G: static FieldGeometry; T](d: var Drift[G, T], h: Hit, seg: int, ring: var Ring) =
  ## a real collision at the end of the flight h: positions, histograms,
  ## pair sums
  let t = h.fl.t
  if t >= d.tmax: d.tmax = t
  when G in {fgElectric, fgParallel}:
    d.cur.e2 = d.cur.e2 + energyIntegral[G](d.f, d.el, h.fl)
  inc d.ncol
  let v2 = advance[G](d.f, d.el, h.fl, d.cur.pos)
  d.cur.st = d.cur.st + t
  let it = int(min(trunc(t.float + 1.0), float64(nTimeBins)))  # 1 ps bins
  d.times[it] = d.times[it] + 1.0
  d.spectrum[h.ie] = d.spectrum[h.ie] + 1.0
  d.w.z = d.cur.pos.z / d.cur.st
  when G in {fgCrossed, fgOblique}:
    d.w.y = d.cur.pos.y / d.cur.st
  when G == fgOblique:
    d.w.x = d.cur.pos.x / d.cur.st
  when G != fgOblique:
    d.cur.vx = d.cur.vx + v2.xx
  when G in {fgElectric, fgParallel}:
    d.cur.vy = d.cur.vy + v2.yy
  # pairs with the collisions ncorln, 2 ncorln, ... before
  let pairs = (when G == fgOblique: seg >= 3 else: d.id != 0)
  if pairs:
    let p = d.cur.pos
    var kdum = 0
    for _ in 1 .. d.dec.ncorst:
      var n = d.ncol + kdum
      if n > d.dec.ncolm: n = n - d.dec.ncolm
      let m = ring.marks[n - 1]
      when G == fgOblique:
        d.cur.st1 = d.cur.st1 + t
      else:
        d.cur.st2 = d.cur.st2 + t
      let sdif = d.cur.st - m.st
      when G in {fgElectric, fgParallel}:
        d.cur.s.xx = d.cur.s.xx + sq(p.x - m.pos.x) * t / sdif
        d.cur.s.yy = d.cur.s.yy + sq(p.y - m.pos.y) * t / sdif
        if seg >= 3:
          d.cur.st1 = d.cur.st1 + t
          d.cur.s.zz = d.cur.s.zz + sq(p.z - m.pos.z - d.w.z * sdif) * t / sdif
      elif G == fgCrossed:
        d.cur.s.xx = d.cur.s.xx + sq(p.x - m.pos.x) * t / sdif
        if seg >= 3:
          d.cur.st1 = d.cur.st1 + t
          d.cur.s.zz = d.cur.s.zz + sq(p.z - m.pos.z - d.w.z * sdif) * t / sdif
          d.cur.s.yy = d.cur.s.yy + sq(p.y - m.pos.y - d.w.y * sdif) * t / sdif
          d.cur.s.yz = d.cur.s.yz + (p.z - m.pos.z - d.w.z * sdif) *
                                    (p.y - m.pos.y - d.w.y * sdif) * t / sdif
          # split the y-z spread along and across the drift direction
          let a2 = sq(d.w.z * sdif) + sq(d.w.y * sdif)
          let b2 = sq(p.z - d.w.z * sdif - m.pos.z) + sq(p.y - d.w.y * sdif - m.pos.y)
          let c2 = sq(p.z - m.pos.z) + sq(p.y - m.pos.y)
          let dl2 = sq(a2 + b2 - c2) / (4.0 * a2)
          let dt2 = b2 - dl2
          d.cur.ls = d.cur.ls + dl2 * t / sdif
          d.cur.ts = d.cur.ts + dt2 * t / sdif
      else:
        let dz = p.z - m.pos.z - d.w.z * sdif
        let dy = p.y - m.pos.y - d.w.y * sdif
        let dx = p.x - m.pos.x - d.w.x * sdif
        d.cur.s.zz = d.cur.s.zz + sq(dz) * t / sdif
        d.cur.s.yy = d.cur.s.yy + sq(dy) * t / sdif
        d.cur.s.xx = d.cur.s.xx + sq(dx) * t / sdif
        d.cur.s.yz = d.cur.s.yz + dz * dy * t / sdif
        d.cur.s.xy = d.cur.s.xy + dx * dy * t / sdif
        d.cur.s.xz = d.cur.s.xz + dx * dz * t / sdif
      kdum = kdum + d.dec.ncorln
  ring.marks[d.ncol - 1] = Mark(pos: d.cur.pos, st: d.cur.st)
  if d.ncol >= d.dec.ncolm:
    inc d.id
    d.ncol = 0

proc count[G: static FieldGeometry; T](d: var Drift[G, T], gas, kind, i: int) {.inline.} =
  when T is ColdGas:
    inc d.tally.kinds[kind]
    inc d.tally.processes[i]
  else:
    inc d.tally[gas].kinds[kind]
    inc d.tally[gas].processes[i]

proc runSegment*[G: static FieldGeometry; T](d: var Drift[G, T],
                                             tabs: MixtureTable | GasTables,
                                             seg: int, collisions: int,
                                             ring: var Ring, rng: var MagRng) =
  ## `collisions` real collisions of segment `seg` (1..10). `tabs`:
  ## `MixtureTable` for `ColdGas`, `GasTables` for `WarmGas`.
  for _ in 1 .. collisions:
    let h = nextHit[G](d.f, tabs, d.target, d.estep, d.el, rng, d.nnull)
    d.record(h, seg, ring)
    template tab: untyped = tableOf(tabs, h)
    let r2 = drand48(rng)
    let i = findProcess(tab, h.ie, r2)
    let ei = energyLoss(tab, i, h.e, rng, moderation = true)
    if tab.ipn[i] > 0:
      # Auger electrons or fluorescence: extra ionisations
      d.extraIons = d.extraIons + tab.augerCount[i]
    d.count(h.gas, tab.iarry[i], i)
    # an excitation transfers to ionisation of another gas with probability
    # PENFRA(1, i)
    if d.penning and tab.penningFraction[i] != 0.0:
      let ran = drand48(rng)
      if not (ran > tab.penningFraction[i]): inc d.extraIons
    let s = scatter(tab, i, h.ie, h.e, ei, h.din, rng,
                    r31First = (G == fgElectric and T is ColdGas))
    d.theta = s.theta
    d.el = leave[G](d.target, d.f, s, h)

func meanEnergy(spectrum: array[1..nEnergyBins, float64], energies: array[1..nEnergyBins, eV],
                tcf: array[1..nEnergyBins, ps⁻¹]): eV•ps =
  ## ∫ e dt from the collision spectrum: each collision at energy e stands
  ## for a time 1/tcf
  for ik in 1 .. nEnergyBins:
    result = result + energies[ik] * spectrum[ik] / tcf[ik]

func rotate[T](t: Tensor[T], rcs, rsn: float64): Tensor[T] =
  ## tensor in the frame of E (fgOblique: B along x, E at BTHETA in x-z)
  result.xx = t.xx * rcs * rcs + t.zz * rsn * rsn + 2.0 * rcs * rsn * t.xz
  result.yy = t.yy
  result.zz = t.xx * rsn * rsn + t.zz * rcs * rcs - 2.0 * rcs * rsn * t.xz
  result.xy = rcs * t.xy + rsn * t.yz
  result.yz = rsn * t.xy - rcs * t.yz
  result.xz = (rcs * rcs - rsn * rsn) * t.xz - rsn * rcs * (t.xx - t.zz)

func diffusion(s: m², t: ps): cm²•s⁻¹ {.inline.} =
  ## D = <Δx²> / 2t from a sum of squared displacements weighted by the
  ## times (Magboltz: C · s / t with C = 5e15)
  (s / (2.0 * t)).to(cm²•s⁻¹)

proc endSegment*[G: static FieldGeometry; T](d: var Drift[G, T], seg: int,
                                             energies: array[1..nEnergyBins, eV],
                                             tcf: array[1..nEnergyBins, ps⁻¹]) =
  ## Estimates of segment `seg`. `energies`: centres of the energy bins;
  ## `tcf`: total collision frequency per bin (of all gases).
  let c1 = d.cur
  let p = d.prev
  when G in {fgElectric, fgParallel}:
    d.vel.z = d.w.z.to(μm•ns⁻¹)
    d.ave = c1.e2 / c1.st
    if d.isotropic:
      d.dif.xx = diffusion(c1.vx, c1.st)
      d.dif.yy = diffusion(c1.vy, c1.st)
      d.dst[seg].xx = diffusion((c1.vx - p.vx), (c1.st - p.st))
      d.dst[seg].yy = diffusion((c1.vy - p.vy), (c1.st - p.st))
    else:
      if c1.st2.float != 0.0:
        d.dif.yy = diffusion(c1.s.yy, c1.st2)
        d.dif.xx = diffusion(c1.s.xx, c1.st2)
        d.dst[seg].xx = diffusion((c1.s.xx - p.s.xx), (c1.st2 - p.st2))
        d.dst[seg].yy = diffusion((c1.s.yy - p.s.yy), (c1.st2 - p.st2))
      else:
        d.dst[seg].xx = 0.0.cm²•s⁻¹
        d.dst[seg].yy = 0.0.cm²•s⁻¹
    if c1.st1.float != 0.0:
      d.dif.zz = diffusion(c1.s.zz, c1.st1)
      d.dst[seg].zz = diffusion((c1.s.zz - p.s.zz), (c1.st1 - p.st1))
    else:
      d.dst[seg].zz = 0.0.cm²•s⁻¹
    d.wst[seg].z = ((c1.pos.z - p.pos.z) / (c1.st - p.st)).to(μm•ns⁻¹)
    d.avest[seg] = (c1.e2 - p.e2) / (c1.st - p.st)
  elif G == fgCrossed:
    d.vel.z = d.w.z.to(μm•ns⁻¹)
    d.vel.y = d.w.y.to(μm•ns⁻¹)
    if c1.st2.float != 0.0: d.dif.xx = diffusion(c1.s.xx, c1.st2)
    if c1.st1.float != 0.0:
      d.dif.yy = diffusion(c1.s.yy, c1.st1)
      d.dif.zz = diffusion(c1.s.zz, c1.st1)
      d.dif.yz = -(diffusion(c1.s.yz, c1.st1))
      d.difln = diffusion(c1.ls, c1.st1)
      d.diftr = diffusion(c1.ts, c1.st1)
    if d.isotropic: d.dif.xx = diffusion(c1.vx, c1.st)
    d.cur.ebar = meanEnergy(d.spectrum, energies, tcf)
    d.ave = d.cur.ebar / c1.st
    d.wst[seg].z = ((c1.pos.z - p.pos.z) / (c1.st - p.st)).to(μm•ns⁻¹)
    d.wst[seg].y = ((c1.pos.y - p.pos.y) / (c1.st - p.st)).to(μm•ns⁻¹)
    d.avest[seg] = (d.cur.ebar - p.ebar) / (c1.st - p.st)
    d.dst[seg] = Tensor[cm²•s⁻¹]()
    d.dlnst[seg] = 0.0.cm²•s⁻¹
    d.dtrst[seg] = 0.0.cm²•s⁻¹
    if seg > 2:
      d.dst[seg].zz = diffusion((c1.s.zz - p.s.zz), (c1.st1 - p.st1))
      d.dst[seg].yy = diffusion((c1.s.yy - p.s.yy), (c1.st1 - p.st1))
      d.dst[seg].yz = diffusion((c1.s.yz - p.s.yz), (c1.st1 - p.st1))
      d.dlnst[seg] = diffusion((c1.ls - p.ls), (c1.st1 - p.st1))
      d.dtrst[seg] = diffusion((c1.ts - p.ts), (c1.st1 - p.st1))
    if c1.st2.float != 0.0: d.dst[seg].xx = diffusion((c1.s.xx - p.s.xx), (c1.st2 - p.st2))
    if d.isotropic: d.dst[seg].xx = diffusion((c1.vx - p.vx), (c1.st - p.st))
  else:
    let rcs = d.rcs
    let rsn = d.rsn
    d.vel = d.w.convert(μm•ns⁻¹)
    d.velRot = vec3(d.vel.z * rsn + d.vel.x * rcs, d.vel.y, d.vel.z * rcs - d.vel.x * rsn)
    d.cur.rot = vec3(c1.pos.z * rsn + c1.pos.x * rcs, c1.pos.y,
                     c1.pos.z * rcs - c1.pos.x * rsn)
    d.cur.ebar = meanEnergy(d.spectrum, energies, tcf)
    d.ave = d.cur.ebar / c1.st
    d.wst[seg].z = ((d.cur.rot.z - p.rot.z) / (c1.st - p.st)).to(μm•ns⁻¹)
    d.wst[seg].y = ((d.cur.rot.y - p.rot.y) / (c1.st - p.st)).to(μm•ns⁻¹)
    d.wst[seg].x = ((d.cur.rot.x - p.rot.x) / (c1.st - p.st)).to(μm•ns⁻¹)
    d.avest[seg] = (d.cur.ebar - p.ebar) / (c1.st - p.st)
    if seg >= 3:
      d.dif.xx = diffusion(c1.s.xx, c1.st1)
      d.dif.yy = diffusion(c1.s.yy, c1.st1)
      d.dif.zz = diffusion(c1.s.zz, c1.st1)
      d.dif.yz = diffusion(c1.s.yz, c1.st1)
      d.dif.xz = diffusion(c1.s.xz, c1.st1)
      d.dif.xy = diffusion(c1.s.xy, c1.st1)
      d.difRot = rotate(d.dif, rcs, rsn)
      d.cur.srot = rotate(c1.s, rcs, rsn)
    d.dst[seg] = Tensor[cm²•s⁻¹]()
    if seg > 2:
      let s = d.cur.srot
      let ps = p.srot
      d.dst[seg].zz = diffusion((s.zz - ps.zz), (c1.st1 - p.st1))
      d.dst[seg].yy = diffusion((s.yy - ps.yy), (c1.st1 - p.st1))
      d.dst[seg].xx = diffusion((s.xx - ps.xx), (c1.st1 - p.st1))
      d.dst[seg].xy = diffusion((s.xy - ps.xy), (c1.st1 - p.st1))
      d.dst[seg].yz = diffusion((s.yz - ps.yz), (c1.st1 - p.st1))
      d.dst[seg].xz = diffusion((s.xz - ps.xz), (c1.st1 - p.st1))
  d.prev = d.cur

func segmentVelocity*[G: static FieldGeometry; T](d: Drift[G, T]): Vector[μm•ns⁻¹] =
  ## drift velocity after the last segment, as printed (µm/ns);
  ## fgOblique: in the frame of E
  when G == fgOblique: d.velRot else: d.vel

func segmentDiffusion*[G: static FieldGeometry; T](d: Drift[G, T]): Tensor[cm²•s⁻¹] =
  ## diffusion tensor after the last segment, as printed (cm²/s);
  ## fgOblique: in the frame of E
  when G == fgOblique: d.difRot else: d.dif

func position*[G: static FieldGeometry; T](d: Drift[G, T]): Vector[m] = d.cur.pos
func time*[G: static FieldGeometry; T](d: Drift[G, T]): ps = d.cur.st
func ringFills*[G: static FieldGeometry; T](d: Drift[G, T]): int = d.id
  ## times the position ring was filled (ID)

func relErr[T](xs: array[1..nSegments, T], first: int, v: T): float64 =
  ## error (%) of the mean of the segment estimates xs[first..10] relative
  ## to v (all in the unit T)
  var t, t2 = 0.0
  for k in first .. nSegments:
    t = t + xs[k].float
    t2 = t2 + xs[k].float * xs[k].float
  let n = float64(nSegments - first + 1)
  100.0 * sqrt((t2 - t * t / n) / (n - 1.0)) / v.float / sqrt(n)

template column(a, field: untyped): untyped =
  var r: array[1..nSegments, typeof(a[1].field)]
  for k in 1 .. nSegments: r[k] = a[k].field
  r

func finish*[G: static FieldGeometry; T](d: Drift[G, T]): Transport =
  ## final estimates and errors
  when G in {fgElectric, fgParallel}:
    result.velErr.z = relErr(column(d.wst, z), 1, d.vel.z)
    result.den = relErr(d.avest, 1, d.ave)
    result.difErr.xx = relErr(column(d.dst, xx), 1, d.dif.xx)
    result.difErr.yy = relErr(column(d.dst, yy), 1, d.dif.yy)
    result.difErr.zz = relErr(column(d.dst, zz), 3, d.dif.zz)
    result.dif = d.dif
    result.difln = d.dif.zz
    result.diftr = (d.dif.xx + d.dif.yy) / 2.0
    result.vel.z = d.vel.z.to(cm•s⁻¹)
    result.dfler = result.difErr.zz
    result.dfter = (result.difErr.xx + result.difErr.yy) / 2.0
  elif G == fgCrossed:
    result.velErr.z = relErr(column(d.wst, z), 1, d.vel.z)
    result.velErr.y = relErr(column(d.wst, y), 1, abs(d.vel.y))
    result.den = relErr(d.avest, 1, d.ave)
    result.difErr.xx = relErr(column(d.dst, xx), 1, d.dif.xx)
    result.difErr.yy = relErr(column(d.dst, yy), 3, d.dif.yy)
    result.difErr.zz = relErr(column(d.dst, zz), 3, d.dif.zz)
    result.difErr.yz = relErr(column(d.dst, yz), 3, abs(d.dif.yz))
    result.dfler = relErr(d.dlnst, 3, d.difln)
    result.dfter = relErr(d.dtrst, 3, d.diftr)
    result.dif = d.dif
    result.difln = d.difln
    result.diftr = d.diftr
    result.vel.z = d.vel.z.to(cm•s⁻¹)
    result.vel.y = d.vel.y.to(cm•s⁻¹)
  else:
    let w = d.velRot
    let dr = d.difRot
    result.velErr.z = relErr(column(d.wst, z), 1, w.z)
    result.velErr.y = relErr(column(d.wst, y), 1, abs(w.y))
    result.velErr.x = relErr(column(d.wst, x), 1, abs(w.x))
    result.den = relErr(d.avest, 1, d.ave)
    result.difErr.zz = relErr(column(d.dst, zz), 3, dr.zz)
    result.difErr.yy = relErr(column(d.dst, yy), 3, dr.yy)
    result.difErr.xx = relErr(column(d.dst, xx), 3, dr.xx)
    result.difErr.xy = relErr(column(d.dst, xy), 3, abs(dr.xy))
    result.difErr.yz = relErr(column(d.dst, yz), 3, abs(dr.yz))
    result.difErr.xz = relErr(column(d.dst, xz), 3, abs(dr.xz))
    result.dif = dr
    result.vel = w.convert(cm•s⁻¹)
  result.ave = d.ave
  # Townsend coefficients: attachments and ionisations per cm of drift
  var ancatt, ancion = 0.0
  for g in 1 .. d.gases:
    when T is ColdGas:
      ancatt = ancatt + float64(d.tally.kinds[5 * g - 2])
      ancion = ancion + float64(d.tally.kinds[5 * g - 3])
    else:
      ancatt = ancatt + float64(d.tally[g].kinds[3])
      ancion = ancion + float64(d.tally[g].kinds[2])
  ancion = ancion + float64(d.extraIons)
  if ancatt != 0.0: result.atter = 100.0 * sqrt(ancatt) / ancatt
  result.att = (ancatt / (d.cur.st * result.vel.z)).to(cm⁻¹)
  if ancion != 0.0: result.alper = 100.0 * sqrt(ancion) / ancion
  result.alpha = (ancion / (d.cur.st * result.vel.z)).to(cm⁻¹)
