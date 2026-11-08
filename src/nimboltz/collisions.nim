## Real collisions: choice of the process and of the scattering angles.
##
## Shared by all Monte Carlo drivers; works in the rest frame of the target
## (the lab frame for a gas at 0 K). The random numbers are drawn in the
## same order as in Magboltz.

import std/math
import units, rm48, tables, vec3

const small* = 1e-20  ## SMALL of /SETP/: floor of energies (in eV) and of ARG1

type
  Scattered* = object
    e*: eV           ## energy after the collision
    dir*: Vec3       ## direction cosines after the collision
    theta*: float64  ## polar scattering angle (rad, against the incoming direction)
    cosTheta*, sinTheta*: float64  ## of the polar angle, as used (backward: cosTheta < 0)
    phi*: float64    ## azimuth (rad)

func findProcess*[N](tab: CollisionTable[N], ie: int, r: float64): int {.inline.} =
  ## First process i with cf[ie, i] >= r (SORT/SORTT, then a linear scan):
  ## a binary step search narrows it down to 4 processes first.
  ## Needs tab.last >= 2 (as SORT, which reads cf[ie, 0] otherwise).
  ## (Loops without `break`, so the GPU compiler takes it.)
  var istep = int(tab.isize)
  var incr = 0
  var i = 0
  var k = 1
  while k <= 12:
    i = incr
    if istep == 2:
      k = 13
    else:
      i = incr + istep
      if i <= tab.last and tab.cf[ie, i] < r:
        incr = incr + istep
      istep = istep div 2
      inc k
  inc i
  while tab.cf[ie, i] < r:
    inc i
  i

func rotate*(din: Vec3, cosTheta, sinTheta, cosPhi, sinPhi: float64): Vec3 {.inline.} =
  ## Direction at polar angle theta and azimuth phi around the direction din.
  let dz = min(din.z, 1.0)
  let argz = sqrt(din.x * din.x + din.y * din.y)
  if argz == 0.0:
    vec3(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta)
  else:
    vec3(din.x * cosTheta - sinTheta / argz * (din.y * cosPhi + din.x * dz * sinPhi),
         din.y * cosTheta + sinTheta / argz * (din.x * cosPhi - din.y * dz * sinPhi),
         dz * cosTheta + argz * sinTheta * sinPhi)

const vibrationMark* = 0.024999000132083893.eV
  ## EIN of a process whose energy loss is "all but 0.025 eV" (0.024999 as
  ## REAL*4 widened, as Magboltz compares it)

proc energyLoss*[N; R](tab: CollisionTable[N], i: int, e: eV, rng: var R,
                       moderation: static bool = false): eV =
  ## Energy lost in process i at energy e. Ionisation shares the energy above
  ## threshold with the secondary electron (a uniform fraction). With
  ## `moderation` (MONTE*, not ELIMIT*), a process marked by an energy loss
  ## of 0.024999 eV takes all but 0.024999 eV of the energy (MB-38: e is the
  ## energy in the frame of the target, also with thermal motion).
  result = tab.ein[i]
  when moderation:
    if abs(result - vibrationMark) <= 0.0001.eV:
      result = e - vibrationMark
  if tab.ipn[i] > 0:
    let r9 = drand48(rng)
    let extra = r9 * (e - result)
    result = extra + result
  if e < result:
    result = e - 0.0001.eV

proc scatter*[N; R](tab: CollisionTable[N], i, ie: int, e, ei: eV, din: Vec3,
                    rng: var R, r31First: static bool = false): Scattered =
  ## Process i at energy e (bin ie) with energy loss ei, incoming direction
  ## din: scattering angles and the outgoing state. `r31First`: draw the
  ## forward/backward choice of INDEX=1 before the angle (MONTE; every other
  ## Magboltz routine draws the angle first).
  let s1 = tab.rgas[i]
  let s2 = s1 * s1 / (s1 - 1.0)
  var r3, f3: float64
  let index = tab.index[i]
  if index == 1:
    when r31First:
      let r31 = drand48(rng)
      r3 = drand48(rng)
    else:
      r3 = drand48(rng)
      let r31 = drand48(rng)
    f3 = 1.0 - r3 * tab.angct[ie, i]
    if r31 > tab.psct[ie, i]: f3 = -f3
  elif index == 2:
    r3 = drand48(rng)
    let epsi = tab.psct[ie, i]
    f3 = 1.0 - 2.0 * r3 * (1.0 - epsi) / (1.0 + epsi * (1.0 - 2.0 * r3))
  else:
    r3 = drand48(rng)
    f3 = 1.0 - 2.0 * r3
  let theta0 = arccos(f3)
  let r4 = drand48(rng)
  let phi0 = 2.0 * PI * r4
  result.phi = phi0
  let arg1 = max(1.0 - s1 * ei / e, small)
  let d = 1.0 - f3 * sqrt(arg1)
  result.e = max(e * (1.0 - ei / (s1 * e) - 2.0 * d / s2), small.eV)
  let q = min(sqrt(e / result.e * arg1) / s1, 1.0)
  result.theta = arcsin(q * sin(theta0))
  var f6 = cos(result.theta)
  let u = (s1 - 1.0) * (s1 - 1.0) / arg1
  if f3 < 0.0 and f3 * f3 > u: f6 = -f6
  result.cosTheta = f6
  result.sinTheta = sin(result.theta)
  result.dir = rotate(din, f6, result.sinTheta, cos(phi0), sin(phi0))
