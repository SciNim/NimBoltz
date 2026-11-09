## GPU walkers against CPU walkers and against NimBoltz's drift run, for
## argon at 0 K, 1000 V/cm, EFINAL 12 eV (deck efinal_ar_t0).
##
##   nim c -r -d:release --path:src --path:<tattletale> tests/gpu/twalk.nim

import std/[math, monotimes, times, strformat]
import nimboltz/[units, rm48, deck, conditions, simulation, magboltz, fields,
                 freeflight, walkers, gpuwalk, api]

proc argon(): Deck =
  result.ngas = 1
  result.nmax = 1
  result.efinal = 12.0.eV
  result.gas[1] = 2
  result.percent[1] = 100.0
  result.tempc = 20.0
  result.torr = 760.0.Torr
  result.efield = 1000.0.V•cm⁻¹

proc main() =
  let d = argon()
  var mixture: ref ColdMixture
  new mixture
  var r = initRun(d, mixture, nil)
  r.mix()
  let f = r.field(fgElectric)
  let fm = flightMaxima(r.mix.table)
  let start = initElectron[fgElectric](f, r.estart, r.theta, startPhi)
  let target = initColdGas(r.mix.bounds, r.efinal)
  var seeds = initMagRng(nil)
  proc walkersOf(n: int): seq[Walker[ColdGas]] =
    for k in 0 ..< n:
      result.add initWalker(start, target, int(drand48(seeds) * 9e8))

  # 1. the same walkers on the CPU and the GPU
  block:
    const n = 256
    const collisions = 2000'i32
    var cpu = walkersOf(n)
    var gpu = cpu
    walkAll[fgElectric, ColdGas](cpu, f, r.mix.table, r.estep, collisions, 0, fm = addr fm[])
    var g = initGpuWalkers[fgElectric, ColdGas](gpu, f, r.mix.table, fm[], r.estep, exact = true)
    g.walk(collisions)
    g.read(gpu)
    var same, sameCollisions = 0
    var worst = 0.0
    for k in 0 ..< n:
      # `float`: Unchained's `==` is approximate
      if cpu[k].pos.z.float == gpu[k].pos.z.float and cpu[k].st.float == gpu[k].st.float and
         cpu[k].e2.float == gpu[k].e2.float:
        inc same
      if cpu[k].nulls == gpu[k].nulls: inc sameCollisions
      worst = max(worst, abs(float(gpu[k].pos.z / cpu[k].pos.z) - 1.0))
    echo &"{n} walkers x {collisions} collisions: {same} identical on CPU and GPU, " &
         &"{sameCollisions} with the same null collisions; largest relative " &
         &"difference of z: {worst:.2e}"

  # 2. drift velocity and mean energy, against NimBoltz
  const total = 40_000_000
  let ref0 = getMonoTime()
  let nb = simulate(Setup(mixture: @[(Gas.argon, 100.0)], efield: 1000.0.V•cm⁻¹,
                          efinal: 12.0.eV, thermal: false, collisions: total))
  let tCpu = (getMonoTime() - ref0).inMilliseconds
  echo &"NimBoltz, {total} collisions in one trajectory ({tCpu} ms): " &
       &"v = {nb.velocity.z}, mean energy {nb.meanEnergy}"
  block:
    # In argon an electron loses ~2m/M = 3e-5 of its energy per collision:
    # it takes ~1e5 collisions to forget the start.
    const n = 8192
    const warmup = 300_000
    const perLaunch = 10_000'i32
    const collisions = 20_000
    var ws = walkersOf(n)
    var g = initGpuWalkers[fgElectric, ColdGas](ws, f, r.mix.table, fm[], r.estep)
    let w0 = getMonoTime()
    for _ in 1 .. warmup div perLaunch: g.walk(perLaunch)
    g.read(ws)
    echo &"GPU warm-up, {n} walkers x {warmup} collisions: " &
         &"{(getMonoTime() - w0).inMilliseconds} ms"
    for w in ws.mitems:
      w.pos = default(typeof(w.pos))
      w.st = 0.0.ps
      w.e2 = 0.0.eV•ps
    g.write(ws)
    let t0 = getMonoTime()
    for _ in 1 .. collisions div perLaunch: g.walk(perLaunch)
    g.read(ws)
    let tGpu = (getMonoTime() - t0).inMilliseconds
    # errors from 16 groups of walkers
    const groups = 16
    var vs, es: seq[float]
    for k in 0 ..< groups:
      let part = ws[k * (n div groups) ..< (k + 1) * (n div groups)]
      vs.add part.driftVelocity.float
      es.add part.meanEnergy.float
    func meanErr(xs: seq[float]): (float, float) =
      let m = xs.sum / float(xs.len)
      var s = 0.0
      for x in xs: s += (x - m) ^ 2
      (m, sqrt(s / float(xs.len - 1) / float(xs.len)))
    let (v, ve) = meanErr(vs)
    let (e, ee) = meanErr(es)
    echo &"GPU, {n} walkers x {collisions} collisions ({tGpu} ms): " &
         &"v = {ws.driftVelocity.float:.4f} ± {ve:.4f} μm/ns, " &
         &"mean energy {ws.meanEnergy.float:.4f} ± {ee:.4f} eV"
    let zv = (v - nb.velocity.z.value.float) /
             sqrt(ve ^ 2 + nb.velocity.z.error.float ^ 2)
    let ze = (e - nb.meanEnergy.value.float) /
             sqrt(ee ^ 2 + nb.meanEnergy.error.float ^ 2)
    echo &"difference / error: velocity {zv:.2f}, mean energy {ze:.2f}"

main()
