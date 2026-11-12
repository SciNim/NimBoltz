## WebGPU walkers against CPU walkers, in a browser: 256 walkers, Ar/CO₂
## 90/10 at 1000 V/cm (B = 1 T in the geometries with B), the gas at 0 K and
## with thermal motion. The CPU walkers' random numbers are the GPU's
## (xoshiro128+: 23 bits, exact in a float32), so both walk the same flights, apart from the single
## precision of the GPU: their states after 1, 10 and 100 collisions differ
## by rounding only, until a rounding difference changes a decision.
##
## Built and run by tests/wasm/webgpu.sh (headless Chromium, SwiftShader).

import std/[math, strformat]
import nimboltz/[units, rm48, deck, conditions, simulation, magboltz, fields, freeflight, walkers,
                 webgpuwalk, api]

proc relDiff(a, b: float64): float64 =
  if a == b: 0.0 else: abs(a - b) / max(abs(a), abs(b))

proc parity[M](name: string, s: Setup) =
  var mixture: ref M
  new mixture
  var d = s.deck
  d.efinal = 5.0.eV
  var r = initRun(d, mixture, nil)
  r.mix()
  withGeometry(geometryOf(d.bmag, d.btheta), G):
    let f = r.field(G)
    let start = initElectron[G](f, r.estart, r.theta, startPhi)
    var seeds = initMagRng(nil)
    when M is WarmMixture:
      type T = WarmGas
      template tabs: untyped = r.mix.tables
      var cpu: seq[Walker[T, Xoshiro128Plus]]
      for k in 0 ..< 256: cpu.add initWalker(start, r.mix.gas, r.nullBands, r.efinal, int(drand48(seeds) * 9e8), Xoshiro128Plus)
    else:
      type T = ColdGas
      template tabs: untyped = r.mix.table
      var cpu: seq[Walker[T, Xoshiro128Plus]]
      let target = initColdGas(r.mix.bounds, r.efinal)
      for k in 0 ..< 256: cpu.add initWalker(start, target, int(drand48(seeds) * 9e8), Xoshiro128Plus)
    var gpu = cpu
    let fm = flightMaxima(tabs)
    var g = initGpuWalkers[G, T, Xoshiro128Plus](gpu, f, tabs, fm[], r.estep)
    var done = 0
    for n in [1, 10, 100]:
      for w in cpu.mitems:  # with the GPU's windows, step by step as the GPU
        walkWith[G, T, NoTally, Xoshiro128Plus](w, f, tabs, addr fm[], gpuWindows, r.estep,
                                        int32(n - done), 0, NoTally())
      g.walk(int32(n - done), counting = false)
      done = n
      g.read(gpu)
      var same = 0
      var de, dz, dt, dx = 0.0
      for k in 0 ..< cpu.len:
        if cpu[k].nulls == gpu[k].nulls and cpu[k].collisions == gpu[k].collisions: inc same
        de = max(de, relDiff(cpu[k].el.e.float, gpu[k].el.e.float))
        dz = max(dz, relDiff(cpu[k].pos.z.float, gpu[k].pos.z.float))
        dx = max(dx, relDiff(cpu[k].pos.x.float, gpu[k].pos.x.float))
        dt = max(dt, relDiff(cpu[k].st.float, gpu[k].st.float))
      var se, ge, sz, gz = 0.0
      for k in 0 ..< cpu.len:
        se += cpu[k].el.e.float; ge += gpu[k].el.e.float
        sz += cpu[k].pos.z.float; gz += gpu[k].pos.z.float
      when defined(debugWalkers):
        if n == 1:
          for k in 0 ..< 3:
            when T is WarmGas:
              for (side, w) in [("cpu", cpu[k]), ("gpu", gpu[k])]:
                echo side, " ", k, " rng=", ($w.rng)[^60 .. ^1], " nulls=", w.nulls,
                     " beyond=", w.beyond, " used=", w.target.maxwell, " st=", w.st.float, " e2=", w.e2.float,
                     " pos=", w.pos, " e=", w.el.e.float, " dir=", w.el.dir
      echo &"{name:<18} {n:>4} coll: same nulls {same:>3}/256, max rel diff e {de:.1e} " &
           &"z {dz:.1e} x {dx:.1e} t {dt:.1e}; mean e {se / 256:.5g}/{ge / 256:.5g} " &
           &"z {sz / 256:.5g}/{gz / 256:.5g}"

proc setupOf(bfield, angle: float64, thermal: bool): Setup =
  Setup(mixture: @[(Gas.argon, 90.0), (Gas.carbonDioxide, 10.0)], efield: 1000.0.V•cm⁻¹,
        bfield: bfield.T, angle: angle, thermal: thermal)

proc twebgpuRun() {.exportc, cdecl.} =
  if not webGpuReady():
    echo "no WebGPU device"
    return
  for thermal in [false, true]:
    let t = if thermal: "thermal" else: "0 K"
    for (label, b, a) in (when defined(debugWalkers): [("B = 0", 0.0, 0.0)] else:
                          [("B = 0", 0.0, 0.0), ("B ∥ E", 1.0, 0.0), ("B ⊥ E", 1.0, 90.0),
                           ("B at 45°", 1.0, 45.0)]):
      let s = setupOf(b, a, thermal)
      if thermal: parity[WarmMixture](label & ", " & t, s)
      else: parity[ColdMixture](label & ", " & t, s)
  echo "DONE"
