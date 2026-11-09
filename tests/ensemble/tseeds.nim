## One run of Ar/CO₂ 90/10 at 0 K, 1000 V/cm, B = 0, 40 M collisions, for
## seed-to-seed comparisons of walkers and NimBoltz's trajectory:
##   tseeds nimboltz <seed>
##   tseeds walkers <seed> <walkers> <warm-up collisions> [threads]
## prints: vz energy Dxx Dyy Dzz (μm/ns, eV, cm²/s)

import std/[os, strutils, strformat]
import nimboltz/[units, rm48, vec3, deck, conditions, simulation, magboltz, fields,
                 freeflight, drift, walkers, ensemble, api]

const total = 40_000_000
let s = Setup(mixture: @[(Gas.argon, 90.0), (Gas.carbonDioxide, 10.0)],
              efield: 1000.0.V•cm⁻¹, thermal: false, collisions: total, efinal: 4.0.eV,
              seed: parseInt(paramStr(2)))
if paramStr(1) == "nimboltz":
  let r = simulate(s)
  echo &"{r.velocity.z.value.float} {r.meanEnergy.value.float} {r.diffusion.xx.value.float} " &
       &"{r.diffusion.yy.value.float} {r.diffusion.zz.value.float}"
else:
  let n = parseInt(paramStr(3))
  let warmup = int32(parseInt(paramStr(4)))
  let threads = if paramCount() >= 5: parseInt(paramStr(5)) else: 4
  var mixture: ref ColdMixture
  new mixture
  var r = initRun(s.deck, mixture, nil)
  r.mix()
  let f = r.field(fgElectric)
  let start = initElectron[fgElectric](f, r.estart, r.theta, startPhi)
  let target = initColdGas(r.mix.bounds, r.efinal)
  var seeds = initMagRng(nil)
  seeds.rm.rm48in(s.seed, 0, 0)
  var ws = newSeq[Walker[ColdGas]](n)
  for w in ws.mitems: w = initWalker(start, target, int(drand48(seeds) * 9e8))
  walkAll[fgElectric, ColdGas](ws, f, r.mix.table, r.estep, warmup, 0, threads)
  for w in ws.mitems: w.reset()
  let half = int32(total div n div 2)
  walkAll[fgElectric, ColdGas](ws, f, r.mix.table, r.estep, half, 0, threads)
  for w in ws.mitems: w.mark()
  walkAll[fgElectric, ColdGas](ws, f, r.mix.table, r.estep, half, 0, threads)
  let tr = transport[fgElectric, ColdGas](ws, 0.0)
  echo &"{tr.vel.z.to(μm•ns⁻¹).float} {tr.ave.float} {tr.dif.xx.float} {tr.dif.yy.float} " &
       &"{tr.dif.zz.float}"
