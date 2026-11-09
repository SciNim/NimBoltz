## Transport coefficients of an ensemble of walkers (`walkers.nim`).
##
## Magboltz follows one electron and estimates diffusion from displacements
## between pairs of its collisions. Here many independent electrons have
## run for times t_k after forgetting their start; their displacements x_k
## spread as 〈(x − w t)(x − w t)ᵀ〉 = 2 D t + c, where the constant c comes
## from the correlation of the velocities over the first collisions. With
## the positions x1_k at times t1_k of a mark halfway (`walkers.mark`),
##   w = Σ x_k / Σ t_k,
##   D = [Σ (x_k − w t_k)(…)ᵀ − Σ (x1_k − w t1_k)(…)ᵀ] / (2 Σ (t_k − t1_k)),
## free of c (without the mark, D is off by c/2t), and the mean energy is
## Σ ∫e dt / Σ t. The walkers are split into blocks of consecutive walkers;
## the scatter of the block estimates gives the errors, as the scatter of
## the 10 segments does in Magboltz. The results fill a `Transport`, in
## Magboltz's frame and with its conventions per geometry (B ⊥ E: diffusion
## along and across the drift in the y-z plane; oblique B: everything in
## the frame of E).

import std/math
import units, vec3, fields, drift, walkers

func sq[T](x: T): auto {.inline.} = x * x

type
  Estimate = object
    vel: Vector[m•ps⁻¹]
    dif: Tensor[cm²•s⁻¹]
    difln, diftr: cm²•s⁻¹
    ave: eV

func estimate[G: static FieldGeometry; T](ws: openArray[Walker[T]], w: Vector[m•ps⁻¹]): Estimate =
  ## of walkers ws, spreads around the drift velocity w (of all walkers)
  var x: Vector[m]
  var t = 0.0.ps
  var e2 = 0.0.eV•ps
  var s: Tensor[m²]
  var ls, ts = 0.0.m²
  let wyz = sqrt(sq(w.y) + sq(w.z))
  var dt = 0.0.ps
  template spread(p: Vector[m], tk: ps, sign: float64) =
    let d = p - w * tk
    s.xx = s.xx + sign * d.x * d.x
    s.yy = s.yy + sign * d.y * d.y
    s.zz = s.zz + sign * d.z * d.z
    s.xy = s.xy + sign * d.x * d.y
    s.xz = s.xz + sign * d.x * d.z
    s.yz = s.yz + sign * d.y * d.z
    when G == fgCrossed:
      let l = (d.y * w.y + d.z * w.z) / wyz  # along the drift
      ls = ls + sign * l * l
      ts = ts + sign * (d.y * d.y + d.z * d.z - l * l)
  for k in ws.low .. ws.high:
    x = x + ws[k].pos
    t = t + ws[k].st
    dt = dt + (ws[k].st - ws[k].st1)
    e2 = e2 + ws[k].e2
    spread(ws[k].pos, ws[k].st, 1.0)
    spread(ws[k].pos1, ws[k].st1, -1.0)
  result.vel = x / t
  result.ave = e2 / t
  template dif(v: m²): cm²•s⁻¹ = (v / (2.0 * dt)).to(cm²•s⁻¹)
  result.dif = Tensor[cm²•s⁻¹](xx: dif(s.xx), yy: dif(s.yy), zz: dif(s.zz),
                               xy: dif(s.xy), xz: dif(s.xz), yz: dif(s.yz))
  result.difln = dif(ls)
  result.diftr = dif(ts)

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

func transport*[G: static FieldGeometry; T](ws: openArray[Walker[T]], btheta: float64,
                                            blocks = 10): Transport =
  ## estimates and errors (blocks of walkers) of the ensemble ws
  doAssert ws.len >= 2 * blocks
  var x: Vector[m]
  var t = 0.0.ps
  for w in ws:
    x = x + w.pos
    t = t + w.st
  var all = estimate[G, T](ws, x / t)
  var est = newSeq[Estimate](blocks)
  for b in 0 ..< blocks:
    est[b] = estimate[G, T](ws.toOpenArray(b * ws.len div blocks, (b + 1) * ws.len div blocks - 1),
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
