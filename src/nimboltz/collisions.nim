## Real collisions: choice of the process and of the scattering angles.
##
## Shared by all Monte Carlo drivers; works in the rest frame of the target
## (the lab frame for a gas at 0 K). The random numbers are drawn in the
## same order as in Magboltz.

import std/math
import rm48, tables, vec3

const small* = 1e-20  ## SMALL of /SETP/: floor of energies and of ARG1

type
  Scattered* = object
    e*: float64      ## energy after the collision (eV)
    dir*: Vec3       ## direction cosines after the collision
    theta*: float64  ## polar scattering angle (rad, against the incoming direction)

func findProcess*[N](tab: CollisionTable[N], ie: int, r: float64): int =
  ## First process i with cf[ie, i] >= r (SORT/SORTT, then a linear scan):
  ## a binary step search narrows it down to 4 processes first.
  ## Needs tab.last >= 2 (as SORT, which reads cf[ie, 0] otherwise).
  var istep = tab.isize
  var incr = 0
  var i = 0
  block search:
    for k in 1 .. 12:
      i = incr
      if istep == 2: break search
      i = incr + istep
      if i <= tab.last and tab.cf[ie, i] < r:
        incr = incr + istep
      istep = istep div 2
  while true:
    inc i
    if not (tab.cf[ie, i] < r): break
  i

func rotate*(din: Vec3, cosTheta, sinTheta, cosPhi, sinPhi: float64): Vec3 =
  ## Direction at polar angle theta and azimuth phi around the direction din.
  let dz = min(din.z, 1.0)
  let argz = sqrt(din.x * din.x + din.y * din.y)
  if argz == 0.0:
    vec3(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta)
  else:
    vec3(din.x * cosTheta - sinTheta / argz * (din.y * cosPhi + din.x * dz * sinPhi),
         din.y * cosTheta + sinTheta / argz * (din.x * cosPhi - din.y * dz * sinPhi),
         dz * cosTheta + argz * sinTheta * sinPhi)

proc collide*[N](tab: CollisionTable[N], i, ie: int, e: float64, din: Vec3,
                 rng: var MagRng): Scattered =
  ## Process i at energy e (bin ie), incoming direction din: energy loss,
  ## scattering angles and the outgoing state.
  let s1 = tab.rgas[i]
  var ei = tab.ein[i]
  if tab.ipn[i] > 0:
    # ionisation: the energy above threshold is shared with the secondary
    let r9 = drand48(rng)
    let extra = r9 * (e - ei)
    ei = extra + ei
  if e < ei:
    ei = e - 0.0001
  let s2 = s1 * s1 / (s1 - 1.0)
  let r3 = drand48(rng)
  var f3: float64
  case tab.index[i]
  of 1:
    let r31 = drand48(rng)
    f3 = 1.0 - r3 * tab.angct[ie, i]
    if r31 > tab.psct[ie, i]: f3 = -f3
  of 2:
    let epsi = tab.psct[ie, i]
    f3 = 1.0 - 2.0 * r3 * (1.0 - epsi) / (1.0 + epsi * (1.0 - 2.0 * r3))
  else:
    f3 = 1.0 - 2.0 * r3
  let theta0 = arccos(f3)
  let r4 = drand48(rng)
  let phi0 = 2.0 * PI * r4
  let arg1 = max(1.0 - s1 * ei / e, small)
  let d = 1.0 - f3 * sqrt(arg1)
  result.e = max(e * (1.0 - ei / (s1 * e) - 2.0 * d / s2), small)
  let q = min(sqrt(e / result.e * arg1) / s1, 1.0)
  result.theta = arcsin(q * sin(theta0))
  var f6 = cos(result.theta)
  let u = (s1 - 1.0) * (s1 - 1.0) / arg1
  if f3 < 0.0 and f3 * f3 > u: f6 = -f6
  result.dir = rotate(din, f6, sin(result.theta), cos(phi0), sin(phi0))
