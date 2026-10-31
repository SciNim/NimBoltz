## Test of the upper energy limit of the collision tables (ELIMIT, ELIMITB,
## ELIMITC and their thermal versions ELIMITT, ELIMITBT, ELIMITCT).
##
## Magboltz chooses the table range EFINAL itself when the input gives 0:
## it starts low and raises EFINAL by √2 until a short Monte Carlo run
## (NMAX/10 collisions for E only, NMAX/20 with B) never takes an electron
## into the last energy bin.

import std/math
import rm48, tables, vec3, fields, collisions, thermal

func energyBin(e, estep: float64): int {.inline.} =
  clamp(int(trunc(e / estep) + 1.0), 1, nEnergyBins)

func sampleSize(geometry: FieldGeometry, nmax: int): int =
  nmax div (if geometry == fgAlongE: 10 else: 20)

proc energyLimitExceeded*[G: static FieldGeometry](
    f: Field, tab: CollisionTable, bounds: var NullBounds, estep, efinal: float64,
    start: Electron, nmax: int, rng: var MagRng, theta: var float64): bool =
  ## Gas at 0 K. True as soon as an electron, starting from `start`, reaches
  ## the last energy bin. `bounds` are adapted as in the main Monte Carlo;
  ## `theta` receives the polar angle of the last scattering (Magboltz leaves
  ## it in THETA, which MONTE then uses as the start direction).
  const nBands = 8
  let deltae = efinal / float64(nBands)
  var el = start
  var tdash = 0.0
  for _ in 1 .. sampleSize(G, nmax):
    var fl: Flight
    var ie: int
    while true:   # free flight with null collisions
      let r1 = drand48(rng)
      let band = clamp(int(trunc(el.e / deltae) + 1.0), 1, nBands)
      let tlim = bounds[band]
      tdash = -(ln(r1) / tlim) + tdash
      fl = fly[G](f, el, tdash)
      ie = energyBin(fl.e, estep)
      if tab.tcf[ie] > tlim:
        # bound too low: undo this flight, raise the bound
        tdash = tdash + ln(r1) / tlim
        bounds[band] = 1.05 * bounds[band]
        continue
      let r5 = drand48(rng)
      if not (r5 > tab.tcf[ie] / tlim): break
    if ie == nEnergyBins: return true
    tdash = 0.0
    let din = directionAt[G](f, el, fl)
    let r2 = drand48(rng)
    let s = collide(tab, findProcess(tab, ie, r2), ie, fl.e, din, rng)
    theta = s.theta
    el = afterCollision[G](f, s.e, s.dir)

proc energyLimitExceeded*[G: static FieldGeometry](
    f: Field, tabs: array[1..maxGases, CollisionTable], gas: ThermalGas,
    estep: float64, start: Electron, nmax: int, rng: var MagRng,
    theta: var float64): bool =
  ## Same with thermal motion of the gas: collisions are computed in the rest
  ## frame of the molecule that is hit, one table per gas.
  var maxwell = initMaxwellSampler(rng)
  var el = start
  var tdash = 0.0
  for _ in 1 .. sampleSize(G, nmax):
    var fl: Flight
    var ie, kgas: int
    var vgas, vrel: Vec3
    var eok: float64
    while true:   # free flight with null collisions
      let r1 = drand48(rng)
      tdash = -(ln(r1) / gas.tcfmx) + tdash
      fl = fly[G](f, el, tdash)
      let vel = velocityAt[G](f, el, fl)
      kgas = gas.pickGas(drand48(rng))
      vgas = maxwell.velocity(gas.vtmb[kgas], rng)
      vrel = vel - vgas
      eok = norm2(vrel) / f.vfac2   # energy in the rest frame of the molecule
      ie = energyBin(eok, estep)
      let r5 = drand48(rng)
      if not (r5 > tabs[kgas].tcf[ie] / gas.tcfmax[kgas]): break
    if ie == nEnergyBins: return true
    tdash = 0.0
    let din = vrel * (1.0 / (f.vfac * sqrt(eok)))
    let r2 = drand48(rng)
    let s = collide(tabs[kgas], findProcess(tabs[kgas], ie, r2), ie, eok, din, rng)
    theta = s.theta
    # back to the lab frame
    let v = s.dir * (f.vfac * sqrt(s.e)) + vgas
    el.e = norm2(v) / f.vfac2
    el.dir = v * (1.0 / (f.vfac * sqrt(el.e)))
    el.vel = v
