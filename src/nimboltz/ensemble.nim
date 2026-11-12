## Transport coefficients of an ensemble of walkers (`walkers.nim`).
##
## Magboltz follows one electron and estimates diffusion from displacements
## between pairs of its collisions. Here many independent electrons run
## after forgetting their start, and stop at checkpoints (`Checkpoints`:
## after equal numbers of collisions) where their positions and times are
## kept. A displacement x over a window between two checkpoints, of
## duration t, spreads around the drift as 〈(x − w t)(x − w t)ᵀ〉 = 2 D t + c,
## where the constant c comes from the correlation of the velocities near
## the ends of the window, the same for all windows longer than that
## correlation. Over the windows of one chunk (sums S₁, T₁ of the spreads
## and the durations, n₁ windows) and of two consecutive chunks (S₂, T₂, n₂),
##   w = Σ x_k / Σ t_k,
##   D = (S₂/n₂ − S₁/n₁) / (2 (T₂/n₂ − T₁/n₁)),
## free of c, and the mean energy is Σ ∫e dt / Σ t. With many checkpoints
## each walker gives many windows, so the error of D no longer depends on
## the number of walkers (two checkpoints: one window of each length, a
## mark halfway). The walkers are split into blocks of consecutive walkers;
## the scatter of the block estimates gives the errors, as the scatter of
## the 10 segments does in Magboltz. The results fill a `Transport`, in
## Magboltz's frame and with its conventions per geometry (B ⊥ E: diffusion
## along and across the drift in the y-z plane; oblique B: everything in
## the frame of E).

import std/math
import units, rm48, vec3, fields, drift, walkers

func sq[T](x: T): auto {.inline.} = x * x

type
  Estimate = object
    vel: Vector[m•ps⁻¹]
    dif: Tensor[cm²•s⁻¹]
    difln, diftr: cm²•s⁻¹
    ave: eV

  Checkpoints* = object
    ## positions and times of n walkers at checkpoints 0 .. m (checkpoint j
    ## of walker k at j · n + k)
    n*: int
    pos*: seq[Vector[m]]
    st*: seq[ps]

func initCheckpoints*[T, R](ws: openArray[Walker[T, R]]): Checkpoints =
  ## checkpoint 0: where the walkers are now
  result.n = ws.len
  for w in ws:
    result.pos.add w.pos
    result.st.add w.st

proc add*[T, R](cp: var Checkpoints, ws: openArray[Walker[T, R]]) =
  ## the next checkpoint: where the walkers are now
  doAssert ws.len == cp.n
  for w in ws:
    cp.pos.add w.pos
    cp.st.add w.st

func chunks*(cp: Checkpoints): int = cp.st.len div cp.n - 1

func estimate[G: static FieldGeometry; T, R](ws: openArray[Walker[T, R]], lo, hi: int,
                                          cp: Checkpoints, w: Vector[m•ps⁻¹]): Estimate =
  ## of walkers lo .. hi, spreads around the drift velocity w (of all walkers)
  var x: Vector[m]
  var t = 0.0.ps
  var e2 = 0.0.eV•ps
  type Spread = object
    s: Tensor[m²]
    ls, ts: m²
    t: ps
    n: float64
  var lag: array[1..2, Spread]
  let wyz = sqrt(sq(w.y) + sq(w.z))
  template window(sp: var Spread, a, b: int) =
    let d = (cp.pos[b] - cp.pos[a]) - w * (cp.st[b] - cp.st[a])
    sp.s.xx = sp.s.xx + d.x * d.x
    sp.s.yy = sp.s.yy + d.y * d.y
    sp.s.zz = sp.s.zz + d.z * d.z
    sp.s.xy = sp.s.xy + d.x * d.y
    sp.s.xz = sp.s.xz + d.x * d.z
    sp.s.yz = sp.s.yz + d.y * d.z
    when G == fgCrossed:
      let l = (d.y * w.y + d.z * w.z) / wyz  # along the drift
      sp.ls = sp.ls + l * l
      sp.ts = sp.ts + (d.y * d.y + d.z * d.z - l * l)
    sp.t = sp.t + (cp.st[b] - cp.st[a])
    sp.n = sp.n + 1.0
  let m = cp.chunks
  for k in lo .. hi:
    x = x + ws[k].pos
    t = t + ws[k].st
    e2 = e2 + ws[k].e2
    for j in 1 .. m:
      lag[1].window((j - 1) * cp.n + k, j * cp.n + k)
      if j >= 2: lag[2].window((j - 2) * cp.n + k, j * cp.n + k)
  result.vel = x / t
  result.ave = e2 / t
  let dt = 2.0 * (lag[2].t / lag[2].n - lag[1].t / lag[1].n)
  func growth(a, b: m², na, nb: float64, dt: ps): cm²•s⁻¹ =
    ((b / nb - a / na) / dt).to(cm²•s⁻¹)
  let (a, b) = (lag[1], lag[2])
  result.dif = Tensor[cm²•s⁻¹](xx: growth(a.s.xx, b.s.xx, a.n, b.n, dt),
                               yy: growth(a.s.yy, b.s.yy, a.n, b.n, dt),
                               zz: growth(a.s.zz, b.s.zz, a.n, b.n, dt),
                               xy: growth(a.s.xy, b.s.xy, a.n, b.n, dt),
                               xz: growth(a.s.xz, b.s.xz, a.n, b.n, dt),
                               yz: growth(a.s.yz, b.s.yz, a.n, b.n, dt))
  result.difln = growth(a.ls, b.ls, a.n, b.n, dt)
  result.diftr = growth(a.ts, b.ts, a.n, b.n, dt)

func rotate[T](t: Tensor[T], rcs, rsn: float64): Tensor[T] =
  ## to the frame of E (fgOblique: B along x, E at BTHETA in x-z), as MONTEC
  result.xx = t.xx * rcs * rcs + t.zz * rsn * rsn + 2.0 * rcs * rsn * t.xz
  result.yy = t.yy
  result.zz = t.xx * rsn * rsn + t.zz * rcs * rcs - 2.0 * rcs * rsn * t.xz
  result.xy = rcs * t.xy + rsn * t.yz
  result.yz = rsn * t.xy - rcs * t.yz
  result.xz = (rcs * rcs - rsn * rsn) * t.xz - rsn * rcs * (t.xx - t.zz)

func relErr(xs: openArray[float64], v: float64): float64 =
  ## error (%) of the mean of the block estimates xs relative to v
  var t, t2 = 0.0
  for x in xs:
    t = t + x
    t2 = t2 + x * x
  let n = float64(xs.len)
  if v == 0.0: 0.0
  else: 100.0 * sqrt(max(t2 - t * t / n, 0.0) / (n - 1.0)) / abs(v) / sqrt(n)

type Quantity = enum qVx, qVy, qVz, qAve, qXx, qYy, qZz, qXy, qXz, qYz, qLn, qTr

func values(e: Estimate): array[Quantity, float64] =
  ## in m/ps, eV, cm²/s
  [e.vel.x.float, e.vel.y.float, e.vel.z.float, e.ave.float, e.dif.xx.float,
   e.dif.yy.float, e.dif.zz.float, e.dif.xy.float, e.dif.xz.float, e.dif.yz.float,
   e.difln.float, e.diftr.float]

func transport*[G: static FieldGeometry; T, R](ws: openArray[Walker[T, R]], cp: Checkpoints,
                                            btheta: float64, blocks = 10): Transport =
  ## estimates and errors (blocks of walkers) of the ensemble ws, at least
  ## two chunks between its checkpoints `cp`
  doAssert ws.len >= 2 * blocks and cp.n == ws.len and cp.chunks >= 2
  var x: Vector[m]
  var t = 0.0.ps
  for w in ws:
    x = x + w.pos
    t = t + w.st
  var all = estimate[G, T, R](ws, 0, ws.high, cp, x / t)
  var est = newSeq[Estimate](blocks)
  for b in 0 ..< blocks:
    est[b] = estimate[G, T, R](ws, b * ws.len div blocks, (b + 1) * ws.len div blocks - 1, cp,
                            all.vel)
  when G == fgOblique:
    # to the frame of E, as MONTEC
    let api = arccos(-1.0)
    let rcs = cos((btheta - 90.0) * api / 180.0)
    let rsn = sin((btheta - 90.0) * api / 180.0)
    proc toE(e: var Estimate) =
      e.vel = vec3(e.vel.z * rsn + e.vel.x * rcs, e.vel.y, e.vel.z * rcs - e.vel.x * rsn)
      e.dif = rotate(e.dif, rcs, rsn)
    toE(all)
    for e in est.mitems: toE(e)
  # errors (%) from the scatter of the blocks
  let v = all.values
  var err: array[Quantity, float64]
  for q in Quantity:
    var xs = newSeq[float64](blocks)
    for b in 0 ..< blocks: xs[b] = est[b].values[q]
    err[q] = relErr(xs, v[q])
  result.vel = all.vel.convert(cm•s⁻¹)
  result.velErr = vec3(err[qVx], err[qVy], err[qVz])
  # the drift components Magboltz has, the others zero by symmetry
  when G in {fgElectric, fgParallel, fgCrossed}:
    result.vel.x = 0.0.cm•s⁻¹
    result.velErr.x = 0.0
  when G in {fgElectric, fgParallel}:
    result.vel.y = 0.0.cm•s⁻¹
    result.velErr.y = 0.0
  result.ave = all.ave
  result.den = err[qAve]
  result.dif = all.dif
  result.difErr = Tensor[float64](xx: err[qXx], yy: err[qYy], zz: err[qZz],
                                  xy: err[qXy], xz: err[qXz], yz: err[qYz])
  when G in {fgElectric, fgParallel}:
    result.difln = all.dif.zz
    result.diftr = (all.dif.xx + all.dif.yy) / 2.0
    result.dfler = err[qZz]
    result.dfter = (err[qXx] + err[qYy]) / 2.0
    result.dif.xy = 0.0.cm²•s⁻¹  # by symmetry; Magboltz prints none
    result.dif.xz = 0.0.cm²•s⁻¹
    result.dif.yz = 0.0.cm²•s⁻¹
  elif G == fgCrossed:
    result.difln = all.difln
    result.diftr = all.diftr
    result.dfler = err[qLn]
    result.dfter = err[qTr]
    result.dif.yz = -result.dif.yz  # Magboltz's sign (MONTEB)
    result.dif.xy = 0.0.cm²•s⁻¹
    result.dif.xz = 0.0.cm²•s⁻¹
  # Townsend coefficients: attachments and ionisations per cm of drift
  var ions, atts = 0.0
  for w in ws:
    ions = ions + float64(w.ions)
    atts = atts + float64(w.attachments)
  if atts != 0.0: result.atter = 100.0 * sqrt(atts) / atts
  result.att = (atts / (t * result.vel.z)).to(cm⁻¹)
  if ions != 0.0: result.alper = 100.0 * sqrt(ions) / ions
  result.alpha = (ions / (t * result.vel.z)).to(cm⁻¹)
