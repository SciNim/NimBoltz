## Free flight of an electron between collisions in uniform E and B fields.
##
## Magboltz writes one Monte Carlo routine per field geometry (MONTE,
## MONTEA/B/C, ELIMIT/B/C, ...). Here the geometry is a static parameter,
## so every instantiation is as specialised as the Fortran copies.
##
## Frame and units, as in Magboltz: E points along z; with B ∥ E, B points
## along z too; otherwise B points along x and E lies in the x-z plane.
## Energies in eV, times in ps, velocities in m/ps, positions in m. Formulas
## and their evaluation order are those of the Fortran, so results agree bit
## for bit.

import std/math
import units, vec3

type
  FieldGeometry* = enum
    fgElectric  ## B = 0 (ELIMIT, MONTE)
    fgParallel  ## B parallel to E: the energy changes as for B = 0, the
                ## velocity ⊥ E gyrates (MONTEA)
    fgCrossed   ## B perpendicular to E (ELIMITB, MONTEB)
    fgOblique   ## B at any other angle to E (ELIMITC, MONTEC)

  Field* = object
    ## E and B, with the derived constants of the free-flight formulas
    efield*: V•cm⁻¹
    vfac*: m•ps⁻¹•eV⁻¹ᐟ²    ## speed per sqrt(energy): v = vfac * sqrt(e)  (CONST9)
    vfac2*: m²•ps⁻²•eV⁻¹    ## vfac²  (CONST10)
    # fgElectric, fgParallel: energy gain e(t) - e = (ap + bp t) t
    bp*: eV•ps⁻²
    f2*: eV¹ᐟ²•ps⁻¹         ## ap = cos θ · f2 · √e
    const5*: eV¹ᐟ²•cm•V⁻¹•ps⁻¹  ## change of the direction cosine along E
    # all but fgCrossed: acceleration term of the position along E (fgOblique: along B)
    f1*: m•ps⁻²
    # fgParallel, fgCrossed, fgOblique
    wb*: ps⁻¹              ## cyclotron frequency (rad/ps)
    # fgCrossed, fgOblique
    eovb*: m•ps⁻¹          ## E×B drift speed for the E component ⊥ B
    ez100*, ex100*: eV•m⁻¹  ## energy gain per displacement along z, x

  Electron* = object
    ## state of an electron right after a collision (or at the start)
    e*: eV                 ## kinetic energy
    dir*: Vec3             ## direction cosines
    vel*: Vector[m•ps⁻¹]     ## velocity; kept up to date only when B ≠ 0

  Flight* = object
    ## free flight of duration t from an `Electron` state
    t*: ps                 ## duration
    e*: eV                 ## energy at the end
    cwt*, swt*: float64    ## cos ωt, sin ωt  (B ≠ 0)
    dz*, dx*: m            ## displacement along z, x (fgCrossed: dz; fgOblique: both)

func geometryOf*(bmag: kGauss, btheta: float64): FieldGeometry =
  ## Field geometry; B ∥ E includes the angle -180° (MONTE dispatches -180°
  ## to MONTEC, see simulation.monteGeometry)
  if bmag.float == 0.0: fgElectric
  elif btheta == 0.0 or abs(btheta) == 180.0: fgParallel
  elif btheta == 90.0: fgCrossed
  else: fgOblique

func initField*(geometry: FieldGeometry, efield: V•cm⁻¹, const1: cm²•V⁻¹•ps⁻²,
                const2: m•cm•V⁻¹•ps⁻², const3: cm•ps⁻¹•eV⁻¹ᐟ², const5: eV¹ᐟ²•cm•V⁻¹•ps⁻¹,
                wb: ps⁻¹, eovb: m•ps⁻¹, btheta: float64): Field =
  ## `const1..5` are CNSTS1 of SETUP, `wb`, `eovb` and `btheta` are /BFLD/.
  result.efield = efield
  result.vfac = const3.to(m•ps⁻¹•eV⁻¹ᐟ²)
  result.vfac2 = result.vfac * result.vfac
  result.wb = wb
  case geometry
  of fgElectric, fgParallel:
    result.bp = efield * efield * const1 * perElectron
    result.f1 = efield * const2
    result.f2 = efield * const3 * perElectron
    result.const5 = const5
  of fgCrossed:
    result.eovb = eovb
    result.ez100 = efield.to(V•m⁻¹) * perElectron
  of fgOblique:
    let rtheta = btheta * arccos(-1.0) / 180.0
    result.ez100 = efield.to(V•m⁻¹) * sin(rtheta) * perElectron
    result.ex100 = efield.to(V•m⁻¹) * cos(rtheta) * perElectron
    result.f1 = efield * const2 * cos(rtheta)
    result.eovb = eovb * sin(rtheta)

func initElectron*[G: static FieldGeometry](f: Field, e: eV, theta, phi: float64): Electron =
  ## electron of energy e moving at polar angle theta, azimuth phi (to E)
  result.e = e
  result.dir = vec3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
  when G != fgElectric:
    result.vel = result.dir * (f.vfac * sqrt(e))

func fly*[G: static FieldGeometry](f: Field, el: Electron, t: ps): Flight {.inline.} =
  ## energy after a free flight of duration t; for fgParallel call
  ## `finishFlight` before asking for the velocity or the position
  result.t = t
  when G in {fgElectric, fgParallel}:
    let ap = el.dir.z * f.f2 * sqrt(el.e)
    result.e = el.e + (ap + f.bp * t) * t
  else:
    let wbt = f.wb * t
    result.cwt = cos(wbt)
    result.swt = sin(wbt)
    result.dz = (el.vel.z * result.swt + (f.eovb - el.vel.y) * (1.0 - result.cwt)) / f.wb
    when G == fgCrossed:
      result.e = el.e + result.dz * f.ez100
    else:
      result.dx = el.vel.x * t + f.f1 * t * t
      result.e = el.e + result.dz * f.ez100 + result.dx * f.ex100

func finishFlight*[G: static FieldGeometry](f: Field, fl: var Flight) {.inline.} =
  ## The gyration angle of fgParallel does not change the energy, so `fly`
  ## leaves it out; this adds it once the flight is accepted.
  when G == fgParallel:
    let wbt = f.wb * fl.t
    fl.cwt = cos(wbt)
    fl.swt = sin(wbt)

func velocityAt*[G: static FieldGeometry](f: Field, el: Electron,
                                          fl: Flight): Vector[m•ps⁻¹] {.inline.}
func directionAt*[G: static FieldGeometry](f: Field, el: Electron, fl: Flight): Vec3 {.inline.} =
  ## direction cosines at the end of the flight
  when G == fgElectric:
    result = el.dir * sqrt(el.e / fl.e)
    result.z = result.z + f.efield * fl.t * f.const5 / sqrt(fl.e)
  elif G == fgParallel:
    let vtot = f.vfac * sqrt(fl.e)
    result.x = (el.vel.x * fl.cwt - el.vel.y * fl.swt) / vtot
    result.y = (el.vel.y * fl.cwt + el.vel.x * fl.swt) / vtot
    result.z = el.dir.z * sqrt(el.e / fl.e) + f.efield * fl.t * f.const5 / sqrt(fl.e)
  else:
    let v = velocityAt[G](f, el, fl)
    result = v / sqrt(norm2(v))

func velocityAt*[G: static FieldGeometry](f: Field, el: Electron, fl: Flight): Vector[m•ps⁻¹] =
  ## velocity at the end of the flight
  when G == fgElectric:
    result = directionAt[G](f, el, fl) * f.vfac * sqrt(fl.e)
  elif G == fgParallel:
    result.x = el.vel.x * fl.cwt - el.vel.y * fl.swt
    result.y = el.vel.y * fl.cwt + el.vel.x * fl.swt
    let vtot = f.vfac * sqrt(fl.e)
    result.z = vtot * (el.dir.z * sqrt(el.e / fl.e) + f.efield * fl.t * f.const5 / sqrt(fl.e))
  else:
    when G == fgCrossed:
      result.x = el.vel.x
    else:
      result.x = el.vel.x + 2.0 * f.f1 * fl.t
    result.y = (el.vel.y - f.eovb) * fl.cwt + el.vel.z * fl.swt + f.eovb
    result.z = el.vel.z * fl.cwt - (el.vel.y - f.eovb) * fl.swt

func advance*[G: static FieldGeometry](f: Field, el: Electron, fl: Flight,
                                       pos: var Vector[m]): tuple[xx, yy: m²] {.inline.} =
  ## Moves `pos` over the flight. Returns the flight's contributions to
  ## Magboltz's transverse diffusion estimate for isotropic scattering
  ## (squared x and y displacement; for fgOblique none).
  let t2 = fl.t * fl.t
  when G in {fgElectric, fgParallel}:
    let a = fl.t * (f.vfac * sqrt(el.e))
    when G == fgElectric:
      pos.x = pos.x + el.dir.x * a
      pos.y = pos.y + el.dir.y * a
      let cx = el.dir.x * (f.vfac * sqrt(el.e))
      let cy = el.dir.y * (f.vfac * sqrt(el.e))
      result = (cx * cx * t2, cy * cy * t2)
    else:
      let dx = (el.vel.x * fl.swt - el.vel.y * (1.0 - fl.cwt)) / f.wb
      pos.x = pos.x + dx
      let dy = (el.vel.y * fl.swt + el.vel.x * (1.0 - fl.cwt)) / f.wb
      pos.y = pos.y + dy
      result = (dx * dx, dy * dy)
    pos.z = pos.z + el.dir.z * a + t2 * f.f1
  else:
    when G == fgCrossed:
      pos.x = pos.x + el.vel.x * fl.t
      result = (el.vel.x * el.vel.x * t2, 0.0.m²)
    else:
      pos.x = pos.x + fl.dx
    pos.y = pos.y + f.eovb * fl.t +
            ((el.vel.y - f.eovb) * fl.swt + el.vel.z * (1.0 - fl.cwt)) / f.wb
    pos.z = pos.z + fl.dz

func energyIntegral*[G: static FieldGeometry](f: Field, el: Electron,
                                              fl: Flight): eV•ps {.inline.} =
  ## ∫ e dt over the flight (fgElectric, fgParallel: exact; the energy is
  ## quadratic in t)
  static: doAssert G in {fgElectric, fgParallel}
  let ap = el.dir.z * f.f2 * sqrt(el.e)
  let a = ap * fl.t
  let b = f.bp * (fl.t * fl.t)
  fl.t * (el.e + a / 2.0 + b / 3.0)

func afterCollision*[G: static FieldGeometry](f: Field, e: eV, dir: Vec3): Electron {.inline.} =
  ## lab-frame state leaving a collision with energy e and direction dir
  result.e = e
  result.dir = dir
  when G != fgElectric:
    result.vel = dir * (f.vfac * sqrt(e))

template withGeometry*(geometry: FieldGeometry, G, body: untyped) =
  ## Runs body with `G` bound to `geometry` as a compile-time constant
  ## (one instantiation per geometry).
  case geometry
  of fgElectric:
    const G {.inject.} = fgElectric
    body
  of fgParallel:
    const G {.inject.} = fgParallel
    body
  of fgCrossed:
    const G {.inject.} = fgCrossed
    body
  of fgOblique:
    const G {.inject.} = fgOblique
    body
