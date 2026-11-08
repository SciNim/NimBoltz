## Thermal motion of the gas (ITHRM = 1): choice of the target gas and of
## its velocity for each collision.

import std/math
import units, rm48, tables, vec3

type
  ThermalGas* = object
    ## per-gas null-collision bounds and thermal speeds (/THRM/, TCFMAX of /LARGET/)
    count*: int                             ## number of gases
    tcfmx*: ps⁻¹                            ## sum of the per-gas bounds
    tcfmxg*: array[1..maxGases, float64]    ## cumulative share of each gas in tcfmx
    tcfmax*: array[1..maxGases, ps⁻¹]       ## null-collision bound of each gas
    vtmb*: array[1..maxGases, m•ps⁻¹]
      ## sqrt(2kT/M) of each gas; with MaxwellSampler's components
      ## (variance 1/2) each velocity component has the spread sqrt(kT/M)

  MaxwellSampler* = object
    ## Gaussian numbers of variance 1/2, drawn six at a time (GERJAN)
    rn: array[1..6, float64]
    used: int

proc refill(m: var MaxwellSampler, rng: var MagRng) =
  for j in countup(1, 5, 2):
    let ran1 = drand48(rng)
    let ran2 = drand48(rng)
    m.rn[j] = sqrt(-ln(ran1)) * cos(ran2 * (2.0 * PI))
    m.rn[j + 1] = sqrt(-ln(ran1)) * sin(ran2 * (2.0 * PI))
  m.used = 0

proc initMaxwellSampler*(rng: var MagRng): MaxwellSampler =
  result.refill(rng)

proc velocity*(m: var MaxwellSampler, vth: m•ps⁻¹, rng: var MagRng): Vector[m•ps⁻¹] {.inline.} =
  ## velocity of a gas molecule with thermal speed vth
  if m.used == 6: m.refill(rng)
  result = vec3(vth * m.rn[m.used + 1], vth * m.rn[m.used + 2], vth * m.rn[m.used + 3])
  m.used += 3

func pickGas*(g: ThermalGas, r: float64): int {.inline.} =
  ## gas hit by a collision, from a uniform random number r
  result = 1
  while g.tcfmxg[result] < r: inc result
