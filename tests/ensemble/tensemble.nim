## Walkers (`walkers.nim`, `ensemble.nim`) on CPU threads against NimBoltz's
## drift run (one trajectory, Magboltz's estimators), for Ar/CO₂ 90/10 at
## 0 K in the four field geometries. Prints the results side by side and
## the z-scores (difference / combined error).
##
##   nim c -r -d:release --path:src tests/ensemble/tensemble.nim [threads]

import std/[math, monotimes, times, strformat, os, strutils, cpuinfo]
import nimboltz/[units, rm48, vec3, deck, conditions, simulation, magboltz, fields,
                 freeflight, drift, walkers, ensemble, api]

const
  total = 40_000_000     # real collisions of either method
  nWalkers = 2000
  warmup = 20_000'i32    # collisions per walker before the sums start

proc setupOf(bfield: float64, angle: float64): Setup =
  Setup(mixture: @[(Gas.argon, 90.0), (Gas.carbonDioxide, 10.0)], efield: 1000.0.V•cm⁻¹,
        bfield: bfield.T, angle: angle, thermal: false, collisions: total)

proc z(a, ea, b, eb: float64): float64 =
  if ea == 0.0 and eb == 0.0: 0.0 else: (a - b) / sqrt(ea * ea + eb * eb)

proc compare(name: string, s: Setup, threads: int) =
  let t0 = getMonoTime()
  let nb = simulate(s)
  let tNb = (getMonoTime() - t0).inMilliseconds
  # the same run, its energy range included, for the walkers' table
  var mixture: ref ColdMixture
  new mixture
  var d = s.deck
  d.efinal = nb.efinal
  var r = initRun(d, mixture, nil)
  r.mix()
  var tr: Transport
  let t1 = getMonoTime()
  withGeometry(geometryOf(d.bmag, d.btheta), G):
    let f = r.field(G)
    let start = initElectron[G](f, r.estart, r.theta, startPhi)
    let target = initColdGas(r.mix.bounds, r.efinal)
    var seeds = initMagRng(nil)
    var ws = newSeq[Walker[ColdGas]](nWalkers)
    for w in ws.mitems: w = initWalker(start, target, int(drand48(seeds) * 9e8))
    walkAll[G, ColdGas](ws, f, r.mix.table, r.estep, warmup, 0, threads)
    for w in ws.mitems: w.reset()
    walkAll[G, ColdGas](ws, f, r.mix.table, r.estep, int32(total div nWalkers), 0, threads)
    tr = transport[G, ColdGas](ws, d.btheta)
  let tW = (getMonoTime() - t1).inMilliseconds
  echo &"== {name}: NimBoltz {tNb} ms, {nWalkers} walkers on {threads} threads {tW} ms"
  proc row(label: string, m: Measurement, v: float64, pct: float64) =
    let a = m.value.float
    let ea = m.error.float
    let eb = abs(v) * pct / 100.0
    let zz = z(a, ea, v, eb)
    echo "  ", alignLeft(label, 8), &" NimBoltz {a:>12.5g} ± {ea:<10.3g} walkers {v:>12.5g} ± {eb:<10.3g} z = {zz:6.2f}"
  let vel = tr.vel.convert(μm•ns⁻¹)
  row("vz", nb.velocity.z, vel.z.float, tr.velErr.z)
  if s.bfield.float != 0.0:
    row("vy", nb.velocity.y, vel.y.float, tr.velErr.y)
    row("vx", nb.velocity.x, vel.x.float, tr.velErr.x)
  row("energy", nb.meanEnergy, tr.ave.float, tr.den)
  row("Dxx", nb.diffusion.xx, tr.dif.xx.float, tr.difErr.xx)
  row("Dyy", nb.diffusion.yy, tr.dif.yy.float, tr.difErr.yy)
  row("Dzz", nb.diffusion.zz, tr.dif.zz.float, tr.difErr.zz)
  if s.bfield.float != 0.0 and s.angle != 0.0:
    row("Dyz", nb.diffusion.yz, tr.dif.yz.float, tr.difErr.yz)
  if s.bfield.float != 0.0 and s.angle != 0.0 and s.angle != 90.0:
    row("Dxy", nb.diffusion.xy, tr.dif.xy.float, tr.difErr.xy)
    row("Dxz", nb.diffusion.xz, tr.dif.xz.float, tr.difErr.xz)
  row("DL", nb.longitudinal, tr.difln.float, tr.dfler)
  row("DT", nb.transverse, tr.diftr.float, tr.dfter)
  row("alpha", nb.ionisation, tr.alpha.float, tr.alper)
  row("eta", nb.attachment, tr.att.float, tr.atter)

let threads = if paramCount() >= 1: parseInt(paramStr(1)) else: countProcessors()
compare("B = 0", setupOf(0.0, 0.0), threads)
compare("B ∥ E, 1 T", setupOf(1.0, 0.0), threads)
compare("B ⊥ E, 1 T", setupOf(1.0, 90.0), threads)
compare("B at 45°, 1 T", setupOf(1.0, 45.0), threads)
