## Speed of the GPU walkers: collisions per second, after a warm-up, for
## CO₂ and Ar/CO₂ 90/10 at 1000 V/cm, at 0 K and with thermal motion, B = 0.
##
##   nim c -r -d:release -d:nimboltzGpu --path:src --path:<tattletale> tests/gpu/tbench.nim

import std/[monotimes, times, strformat]
import nimboltz/[units, rm48, deck, conditions, simulation, magboltz, fields, freeflight,
                 walkers, gpuwalk, api]

const
  nWalkers = 16384
  warmup = 3000'i32
  collisions = 4000'i32

proc bench[M](name: string, s: Setup) =
  var mixture: ref M
  new mixture
  var q = s
  q.collisions = 1_000_000
  var d = s.deck
  d.efinal = simulate(q).efinal  # Magboltz's energy range
  var r = initRun(d, mixture, nil)
  r.mix()
  const G = fgElectric
  let f = r.field(G)
  let start = initElectron[G](f, r.estart, r.theta, startPhi)
  var seeds = initMagRng(nil)
  when M is WarmMixture:
    type T = WarmGas
    template tabs: untyped = r.mix.tables
    var ws: seq[Walker[T]]
    for k in 0 ..< nWalkers: ws.add initWalker(start, r.mix.gas, r.nullBands, r.efinal, int(drand48(seeds) * 9e8))
  else:
    type T = ColdGas
    template tabs: untyped = r.mix.table
    var ws: seq[Walker[T]]
    let target = initColdGas(r.mix.bounds, r.efinal)
    for k in 0 ..< nWalkers: ws.add initWalker(start, target, int(drand48(seeds) * 9e8))
  let fm = flightMaxima(tabs)
  var g = initGpuWalkers[G, T](ws, f, tabs, fm[], r.estep)
  g.walk(warmup, counting = false)
  g.read(ws)
  var dt = Inf
  for _ in 1 .. 3:  # the best of 3: the GPU's clock varies
    let t0 = getMonoTime()
    g.walk(collisions)
    g.read(ws)
    dt = min(dt, float((getMonoTime() - t0).inMicroseconds))
  var nulls = 0.0
  for w in ws: nulls += float(w.nulls)
  let n = float(nWalkers) * float(warmup + 3 * collisions)
  echo &"{name:<22} {float(nWalkers) * float(collisions) / dt:7.1f} M collisions/s " &
       &"({nulls / n:.2f} null per real, EFINAL {d.efinal.float:.3g} eV)"

proc setupOf(mixture: seq[(Gas, float64)]): Setup =
  Setup(mixture: mixture, efield: 1000.0.V•cm⁻¹, thermal: false)

bench[ColdMixture]("CO2, 0 K", setupOf(@[(Gas.carbonDioxide, 100.0)]))
bench[ColdMixture]("Ar/CO2 90/10, 0 K", setupOf(@[(Gas.argon, 90.0), (Gas.carbonDioxide, 10.0)]))
bench[WarmMixture]("CO2, thermal", setupOf(@[(Gas.carbonDioxide, 100.0)]))
bench[WarmMixture]("Ar/CO2 90/10, thermal", setupOf(@[(Gas.argon, 90.0), (Gas.carbonDioxide, 10.0)]))
