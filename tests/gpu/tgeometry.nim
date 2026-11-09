## GPU walkers in the four field geometries, Ar/CO₂ 90/10, 1000 V/cm,
## B = 1 T, the gas at 0 K and with thermal motion:
## (1) the same walkers on CPU and GPU take the same collisions;
## (2) the drift run on GPU walkers (`simulate` with `wbGpu`) against
##     Magboltz's one trajectory.
##
##   nim c -r -d:release -d:nimboltzGpu --path:src --path:<tattletale> tests/gpu/tgeometry.nim

import std/[math, monotimes, times, strformat, strutils]
import nimboltz/[units, rm48, deck, conditions, simulation, magboltz, fields, freeflight, walkers,
                 gpuwalk, api]

const total = 40_000_000

proc z(a, b: Measurement): float64 =
  let e = sqrt(float(a.error) ^ 2 + float(b.error) ^ 2)
  if e == 0.0: 0.0 else: float(b.value - a.value) / e

proc parity[M](name: string, s: Setup, efinal: eV) =
  ## 256 walkers, 2000 collisions on the CPU and on the GPU (no FMA)
  var mixture: ref M
  new mixture
  var d = s.deck
  d.efinal = efinal
  var r = initRun(d, mixture, nil)
  r.mix()
  withGeometry(geometryOf(d.bmag, d.btheta), G):
    let f = r.field(G)
    let start = initElectron[G](f, r.estart, r.theta, startPhi)
    var seeds = initMagRng(nil)
    when M is WarmMixture:
      type T = WarmGas
      template tabs: untyped = r.mix.tables
      var cpu: seq[Walker[T]]
      for k in 0 ..< 256: cpu.add initWalker(start, r.mix.gas, r.nullBands, r.efinal, int(drand48(seeds) * 9e8))
    else:
      type T = ColdGas
      template tabs: untyped = r.mix.table
      var cpu: seq[Walker[T]]
      let target = initColdGas(r.mix.bounds, r.efinal)
      for k in 0 ..< 256: cpu.add initWalker(start, target, int(drand48(seeds) * 9e8))
    var gpu = cpu
    let fm = flightMaxima(tabs)
    walkAll[G, T](cpu, f, tabs, r.estep, 2000, 0, fm = addr fm[])
    var g = initGpuWalkers[G, T](gpu, f, tabs, fm[], r.estep, exact = true)
    g.walk(2000)
    g.read(gpu)
    var sameNulls = 0
    var worst = 0.0
    for k in 0 ..< cpu.len:
      if cpu[k].nulls == gpu[k].nulls: inc sameNulls
      worst = max(worst, abs(float(gpu[k].pos.z / cpu[k].pos.z) - 1.0))
    echo &"== {name}: {sameNulls}/256 walkers with the same null collisions on CPU and " &
         &"GPU, z differs by ≤ {worst:.1e} (relative)"

proc compare(name: string, s: Setup) =
  let t0 = getMonoTime()
  let one = simulate(s)
  let t1 = getMonoTime()
  var w = s
  w.walkers = WalkerSettings(backend: wbGpu)
  let many = simulate(w)
  let t2 = getMonoTime()
  if s.thermal: parity[WarmMixture](name, s, one.efinal)
  else: parity[ColdMixture](name, s, one.efinal)
  echo &"   one trajectory {(t1 - t0).inMilliseconds} ms, GPU walkers " &
       &"{(t2 - t1).inMilliseconds} ms; {many.warnings}"
  proc row(label: string, a, b: Measurement) =
    echo "   ", alignLeft(label, 8), &" {a.value.float:>12.5g} ± {a.error.float:<10.3g} " &
         &"GPU {b.value.float:>12.5g} ± {b.error.float:<10.3g} z = {z(a, b):6.2f}"
  row("vz", one.velocity.z, many.velocity.z)
  row("vy", one.velocity.y, many.velocity.y)
  row("energy", one.meanEnergy, many.meanEnergy)
  row("DL", one.longitudinal, many.longitudinal)
  row("DT", one.transverse, many.transverse)

proc setupOf(bfield, angle: float64, thermal: bool): Setup =
  Setup(mixture: @[(Gas.argon, 90.0), (Gas.carbonDioxide, 10.0)], efield: 1000.0.V•cm⁻¹,
        bfield: bfield.T, angle: angle, thermal: thermal, collisions: total)

for thermal in [false, true]:
  let t = if thermal: "thermal" else: "0 K"
  compare(&"B = 0, {t}", setupOf(0.0, 0.0, thermal))
  compare(&"B ∥ E, {t}", setupOf(1.0, 0.0, thermal))
  compare(&"B ⊥ E, {t}", setupOf(1.0, 90.0, thermal))
  compare(&"B at 45°, {t}", setupOf(1.0, 45.0, thermal))
