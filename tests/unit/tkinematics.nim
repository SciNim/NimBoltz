## Properties of the free-flight and collision kernels (fields.nim,
## collisions.nim, thermal.nim) that must survive any rewrite.

import std/[unittest, math, random]
import nimboltz/[units, vec3, tables, fields, collisions, thermal, rm48, deck, conditions]

proc unitRandom(r: var Rand): Vec3 =
  let ct = r.rand(-1.0 .. 1.0)
  let phi = r.rand(0.0 .. 2 * PI)
  let st = sqrt(1 - ct * ct)
  vec3(st * cos(phi), st * sin(phi), ct)

func dot(a, b: Vec3): float64 = a.x * b.x + a.y * b.y + a.z * b.z

proc field(g: FieldGeometry, efield, bmag, btheta: float64): Field =
  ## (V/cm, kG, degrees) with the conditions of SETUP
  var d = Deck(efield: efield.V•cm⁻¹, bmag: bmag.kGauss, btheta: btheta,
               torr: 760.0.Torr, ngas: 1)
  d.gas[1] = 2
  d.percent[1] = 100.0
  let c = conditions(d)
  initField(g, d.efield, const1, const2, const3, const5, c.wb, c.eovb, btheta)

suite "collisions":
  test "rotate keeps unit length and turns by theta":
    var r = initRand(1)
    for _ in 1 .. 10_000:
      let din = r.unitRandom
      let theta = r.rand(0.0 .. PI)
      let phi = r.rand(0.0 .. 2 * PI)
      let d = rotate(din, cos(theta), sin(theta), cos(phi), sin(phi))
      check abs(norm2(d) - 1) < 1e-12
      check abs(dot(d, din) - cos(theta)) < 1e-12

  test "rotate around the z axis":
    let d = rotate(vec3(0.0, 0.0, 1.0), cos(0.3), sin(0.3), cos(1.1), sin(1.1))
    check d == vec3(cos(1.1) * sin(0.3), sin(1.1) * sin(0.3), cos(0.3))

  test "findProcess is the first process with cf >= r":
    var tab: CollisionTable[64]
    var r = initRand(2)
    for last in 2 .. 64:   # a mixture has at least two processes
      tab.last = last
      tab.isize = 1
      while tab.isize * 2 <= last: tab.isize *= 2   # as MIXER
      for ie in 1 .. 3:
        var c = 0.0
        for i in 1 .. last:
          c += r.rand(1.0)
          tab.cf[ie, i] = c
        for i in 1 .. last: tab.cf[ie, i] = tab.cf[ie, i] / c
        tab.cf[ie, last] = 1.0
        for _ in 1 .. 200:
          let x = if r.rand(1.0) < 0.2: tab.cf[ie, r.rand(1 .. last)] else: r.rand(1.0)
          var naive = 1
          while tab.cf[ie, naive] < x: inc naive
          check findProcess(tab, ie, x) == naive

suite "fields":
  test "E only: unit direction, transverse velocity conserved":
    let f = field(fgElectric, 1000.0, 0, 0)
    let el = initElectron[fgElectric](f, 2.0.eV, 0.7, 0.4)
    let fl = fly[fgElectric](f, el, 3.0.ps)
    let d = directionAt[fgElectric](f, el, fl)
    check abs(norm2(d) - 1) < 1e-12
    # transverse momentum is conserved
    let v0 = el.dir * (f.vfac * sqrt(el.e))
    let v1 = velocityAt[fgElectric](f, el, fl)
    check abs(v1.x - v0.x) < 1e-12 * abs(v0.x)
    check abs(v1.y - v0.y) < 1e-12 * abs(v0.y)

  test "E only: the energy gained is the potential the electron falls through":
    let f = field(fgElectric, 1000.0, 0, 0)
    let el = initElectron[fgElectric](f, 2.0.eV, 0.7, 0.4)
    let fl = fly[fgElectric](f, el, 3.0.ps)
    var pos: Vector[m]
    discard advance[fgElectric](f, el, fl, pos)
    let gain = f.efield * pos.z * perElectron
    check abs(float64((fl.e - el.e - gain) / gain)) < 1e-12

  test "B ⊥ E: E×B motion leaves the energy unchanged after a full turn":
    let f = field(fgCrossed, 1000.0, 10.0, 90)
    var el = initElectron[fgCrossed](f, 1.0.eV, 1.0, 0.5)
    let fl = fly[fgCrossed](f, el, 2 * PI / f.wb)
    check abs(fl.e - el.e) < 1e-9.eV
    let v = velocityAt[fgCrossed](f, el, fl)
    check abs(v.y - el.vel.y) < 1e-15.m•ps⁻¹ and abs(v.z - el.vel.z) < 1e-15.m•ps⁻¹

  test "B ⊥ E: the drift averages to E/B":
    let f = field(fgCrossed, 1000.0, 10.0, 90)
    # 1000 V/cm across 10 kG = 1 T: 1e5 m/s
    check abs(f.eovb.to(m•s⁻¹).float - 1e5) < 1e-9

suite "thermal":
  test "pickGas follows the cumulative shares":
    let g = ThermalGas(tcfmxg: [0.25, 0.75, 1.0, 1.0, 1.0, 1.0])
    check g.pickGas(0.1) == 1
    check g.pickGas(0.25) == 1
    check g.pickGas(0.5) == 2
    check g.pickGas(0.99) == 3

  test "Maxwell velocity components: mean 0, variance 1/2":
    var rng = initMagRng(nil)
    var m = initMaxwellSampler(rng)
    var s, s2: array[3, float64]
    const n = 300_000
    for _ in 1 .. n:
      let v = m.velocity(1.0.m•ps⁻¹, rng)
      for k, c in [v.x.float, v.y.float, v.z.float]:
        s[k] += c
        s2[k] += c * c
    for k in 0 .. 2:
      check abs(s[k] / n) < 0.01
      # GERJAN draws sqrt(-ln u)·cos(2πv): variance 1/2 per component
      check abs(s2[k] / n - 0.5) < 0.01
