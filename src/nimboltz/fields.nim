## Free flight of an electron between collisions in uniform E and B fields.
##
## Magboltz writes one Monte Carlo routine per field geometry (MONTE,
## MONTEA/B/C, ELIMIT/B/C, ...). Here the geometry is a static parameter,
## so every instantiation is as specialised as the Fortran copies.
##
## Frame and units, as in Magboltz: E points along z; with B ⊥ E, B points
## along x. Energies in eV, times in ps, velocities in m/ps. Formulas and
## their evaluation order are those of the Fortran, so results agree bit
## for bit.

import std/math
import vec3

type
  FieldGeometry* = enum
    fgAlongE   ## B = 0, or B parallel to E: only the motion along E changes
               ## the energy (ELIMIT; MONTE and MONTEA)
    fgCrossed  ## B perpendicular to E (ELIMITB, MONTEB)
    fgOblique  ## B at any other angle to E (ELIMITC, MONTEC)

  Field* = object
    ## E and B, with the derived constants of the free-flight formulas
    efield*: float64   ## E (V/cm)
    vfac*: float64     ## speed per sqrt(energy): v = vfac * sqrt(e)  (CONST9)
    vfac2*: float64    ## vfac²  (CONST10)
    # fgAlongE
    bp, f2, const5: float64
    # fgCrossed, fgOblique
    wb*: float64       ## cyclotron frequency (rad/ps)
    eovb*: float64     ## E×B drift speed for the E component ⊥ B (m/ps)
    ez100, ex100, f1: float64

  Electron* = object
    ## state of an electron right after a collision (or at the start)
    e*: float64        ## kinetic energy (eV)
    dir*: Vec3         ## direction cosines
    vel*: Vec3         ## velocity (m/ps); kept up to date only when B ≠ 0

  Flight* = object
    ## free flight of duration t from an `Electron` state
    t*: float64        ## duration (ps)
    e*: float64        ## energy at the end (eV)
    cwt, swt: float64  ## cos ωt, sin ωt  (B ≠ 0)

func geometryOf*(bmag, btheta: float64): FieldGeometry =
  ## Magboltz's choice of routine (main program)
  if bmag == 0.0 or btheta == 0.0 or abs(btheta) == 180.0: fgAlongE
  elif btheta == 90.0: fgCrossed
  else: fgOblique

func initField*(geometry: FieldGeometry, efield, const1, const2, const3,
                const5, wb, eovb, btheta: float64): Field =
  ## `const1..5` are CNSTS1 of SETUP, `wb`, `eovb` and `btheta` are /BFLD/.
  result.efield = efield
  result.vfac = const3 * 0.01
  result.vfac2 = result.vfac * result.vfac
  result.wb = wb
  case geometry
  of fgAlongE:
    result.bp = efield * efield * const1
    result.f2 = efield * const3
    result.const5 = const5
  of fgCrossed:
    result.eovb = eovb
    result.ez100 = efield * 100.0
  of fgOblique:
    let rtheta = btheta * arccos(-1.0) / 180.0
    result.ez100 = efield * 100.0 * sin(rtheta)
    result.ex100 = efield * 100.0 * cos(rtheta)
    result.f1 = efield * const2 * cos(rtheta)
    result.eovb = eovb * sin(rtheta)

func initElectron*[G: static FieldGeometry](f: Field, e, theta, phi: float64): Electron =
  ## electron of energy e moving at polar angle theta, azimuth phi (to E)
  result.e = e
  result.dir = vec3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
  when G != fgAlongE:
    result.vel = result.dir * (f.vfac * sqrt(e))

func fly*[G: static FieldGeometry](f: Field, el: Electron, t: float64): Flight =
  ## energy after a free flight of duration t
  result.t = t
  when G == fgAlongE:
    let ap = el.dir.z * f.f2 * sqrt(el.e)
    result.e = el.e + (ap + f.bp * t) * t
  else:
    let wbt = f.wb * t
    result.cwt = cos(wbt)
    result.swt = sin(wbt)
    let dz = (el.vel.z * result.swt + (f.eovb - el.vel.y) * (1.0 - result.cwt)) / f.wb
    when G == fgCrossed:
      result.e = el.e + dz * f.ez100
    else:
      let dx = el.vel.x * t + f.f1 * t * t
      result.e = el.e + dz * f.ez100 + dx * f.ex100

func velocityAt*[G: static FieldGeometry](f: Field, el: Electron, fl: Flight): Vec3
func directionAt*[G: static FieldGeometry](f: Field, el: Electron, fl: Flight): Vec3 =
  ## direction cosines at the end of the flight
  when G == fgAlongE:
    result = el.dir * sqrt(el.e / fl.e)
    result.z = result.z + f.efield * fl.t * f.const5 / sqrt(fl.e)
  else:
    let v = velocityAt[G](f, el, fl)
    result = v / sqrt(norm2(v))

func velocityAt*[G: static FieldGeometry](f: Field, el: Electron, fl: Flight): Vec3 =
  ## velocity at the end of the flight
  when G == fgAlongE:
    result = directionAt[G](f, el, fl) * f.vfac * sqrt(fl.e)
  else:
    when G == fgCrossed:
      result.x = el.vel.x
    else:
      result.x = el.vel.x + 2.0 * f.f1 * fl.t
    result.y = (el.vel.y - f.eovb) * fl.cwt + el.vel.z * fl.swt + f.eovb
    result.z = el.vel.z * fl.cwt - (el.vel.y - f.eovb) * fl.swt

func afterCollision*[G: static FieldGeometry](f: Field, e: float64, dir: Vec3): Electron =
  ## lab-frame state leaving a collision with energy e and direction dir
  result.e = e
  result.dir = dir
  when G != fgAlongE:
    result.vel = dir * (f.vfac * sqrt(e))

template withGeometry*(geometry: FieldGeometry, G, body: untyped) =
  ## Runs body with `G` bound to `geometry` as a compile-time constant
  ## (one instantiation per geometry).
  case geometry
  of fgAlongE:
    const G {.inject.} = fgAlongE
    body
  of fgCrossed:
    const G {.inject.} = fgCrossed
    body
  of fgOblique:
    const G {.inject.} = fgOblique
    body
