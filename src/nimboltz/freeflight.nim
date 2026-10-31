## Free flights with null collisions: from one real collision to the next.
##
## The total collision frequency depends on the energy, which changes during
## the flight. Magboltz samples flight times with a constant upper bound of
## the frequency and accepts a collision as real with probability
## (real frequency / bound); the rejected ones are "null" collisions.
##
## Two targets, chosen by overloading:
## - `ColdGas`: the gas at 0 K, one collision table for the mixture. The
##   bound depends on the energy at the start of the flight (8 bands) and is
##   raised by 5% whenever it turns out too low (MONTE, ELIMIT).
## - `WarmGas`: thermal motion of the gas, one table per gas. Each collision
##   picks a gas and a molecule velocity; the collision is computed in the
##   rest frame of that molecule (MONTET, ELIMITT).
##
## The random numbers are drawn in the same order as in Magboltz.

import std/math
import rm48, tables, vec3, fields, collisions, thermal

type
  ColdGas* = object
    ## walker state of the null-collision method at 0 K
    bounds*: NullBounds
    bandWidth*: float64   ## width of the 8 bands of starting energy (EFINAL/8)

  WarmGas* = object
    ## thermal motion: per-gas bounds and speeds, and the walker's Maxwell numbers
    gas*: ThermalGas
    maxwell*: MaxwellSampler

  MixtureTable* = CollisionTable[maxMixProcesses]
  GasTables* = array[1..maxGases, CollisionTable[maxGasProcesses]]

  Hit* = object
    ## a real collision
    fl*: Flight
    ie*: int        ## energy bin in the frame of the target
    e*: float64     ## energy in the frame of the target (eV)
    din*: Vec3      ## direction in the frame of the target
    gas*: int       ## gas hit (WarmGas; 1 for ColdGas)
    vgas*: Vec3     ## velocity of the molecule hit (WarmGas)

func energyBin*(e, estep: float64): int {.inline.} =
  clamp(int(trunc(e / estep) + 1.0), 1, nEnergyBins)

func initColdGas*(bounds: NullBounds, efinal: float64): ColdGas =
  ColdGas(bounds: bounds, bandWidth: efinal / float64(bounds.len))

proc initWarmGas*(gas: ThermalGas, rng: var MagRng): WarmGas =
  WarmGas(gas: gas, maxwell: initMaxwellSampler(rng))

template tableOf*(tab: MixtureTable, h: Hit): untyped = tab
template tableOf*(tabs: GasTables, h: Hit): untyped = tabs[h.gas]

proc nextHit*[G: static FieldGeometry](f: Field, tab: MixtureTable, target: var ColdGas,
                                       estep: float64, el: Electron, rng: var MagRng,
                                       nnull: var int): Hit =
  ## Gas at 0 K: next real collision of an electron leaving `el`.
  var tdash = 0.0
  while true:
    let r1 = drand48(rng)
    let band = clamp(int(trunc(el.e / target.bandWidth) + 1.0), 1, target.bounds.len)
    let tlim = target.bounds[band]
    tdash = -(ln(r1) / tlim) + tdash
    result.fl = fly[G](f, el, tdash)
    result.ie = energyBin(result.fl.e, estep)
    if tab.tcf[result.ie] > tlim:
      # bound too low: undo this flight, raise the bound
      tdash = tdash + ln(r1) / tlim
      target.bounds[band] = 1.05 * target.bounds[band]
      continue
    let r5 = drand48(rng)
    if r5 <= tab.tcf[result.ie] / tlim: break
    inc nnull
  finishFlight[G](f, result.fl)
  result.e = result.fl.e
  result.din = directionAt[G](f, el, result.fl)
  result.gas = 1

proc nextHit*[G: static FieldGeometry](f: Field, tabs: GasTables, target: var WarmGas,
                                       estep: float64, el: Electron, rng: var MagRng,
                                       nnull: var int): Hit =
  ## Thermal motion: next real collision; energy and direction in the rest
  ## frame of the molecule hit.
  var tdash = 0.0
  var vrel: Vec3
  while true:
    let r1 = drand48(rng)
    tdash = -(ln(r1) / target.gas.tcfmx) + tdash
    result.fl = fly[G](f, el, tdash)
    finishFlight[G](f, result.fl)
    let vel = velocityAt[G](f, el, result.fl)
    result.gas = target.gas.pickGas(drand48(rng))
    result.vgas = target.maxwell.velocity(target.gas.vtmb[result.gas], rng)
    vrel = vel - result.vgas
    result.e = norm2(vrel) / f.vfac2
    result.ie = energyBin(result.e, estep)
    let r4 = drand48(rng)
    if not (r4 > tabs[result.gas].tcf[result.ie] / target.gas.tcfmax[result.gas]): break
    inc nnull
  result.din = vrel * (1.0 / (f.vfac * sqrt(result.e)))

func leave*[G: static FieldGeometry](target: ColdGas, f: Field, s: Scattered, h: Hit): Electron {.inline.} =
  ## state after the collision, in the lab frame
  afterCollision[G](f, s.e, s.dir)

func leave*[G: static FieldGeometry](target: WarmGas, f: Field, s: Scattered, h: Hit): Electron {.inline.} =
  ## state after the collision, transformed back to the lab frame
  let v = s.dir * (f.vfac * sqrt(s.e)) + h.vgas
  result.e = norm2(v) / f.vfac2
  result.dir = v * (1.0 / (f.vfac * sqrt(result.e)))
  result.vel = v
