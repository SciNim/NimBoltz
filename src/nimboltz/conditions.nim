## SETUP and SETUPT: the physical constants and what follows from a deck
## before any cross section is computed (gas densities, decorrelation of
## the diffusion estimate, number of collisions, cyclotron frequency).

import std/math
import deck, tables, drift

const
  nscale* {.intdefine.} = 40_000_000
    ## real collisions per unit of NMAX (the tests use less, see tests/regress)

  # constants of Magboltz (2010 values)
  echarg* = 1.602176565e-19   ## elementary charge (C)
  emass* = 9.10938291e-31     ## electron mass (kg)
  amu* = 1.660538921e-27      ## atomic mass unit (kg)
  boltz* = 0.000086173324     ## Boltzmann constant (eV/K)
  awb = 17588200880.0         ## e/m (C/kg), in units chosen by Magboltz
  alosch* = 2.6867805e+19     ## Loschmidt number (1/cm³)
  abzero* = 273.15            ## 0 °C (K)
  atmos* = 760.0              ## 1 atm (Torr)

  # CNSTS1 of SETUP: factors of the free-flight formulas
  const1* = awb / 2.0 * 1e-19
  const2* = const1 * 0.01
  const3* = sqrt(0.2 * awb) * 1e-9
  const4* = const3 * alosch * 1e-15
  const5* = const3 / 2.0

  naniso* = 2
    ## elastic angular distributions of the gas routines: 2 Okhrimovskyy
    ## (fixed in Magboltz; 1 Capitelli-Longo, 0 isotropic)
  tmax* = 100.0               ## flight time covered by the null-collision bounds (ps)
  startTheta* = 0.785         ## start direction of the first electron (polar angle,
  startPhi* = 0.1             ## azimuth)

type
  SetupError* = enum
    seNone
    seGasInput       ## NGAS not in 1..6, a gas or share missing, shares not 100 %
    seNmaxOverflow   ## NMAX · nscale does not fit an integer

  Conditions* = object
    nmax*: int                         ## real collisions (NMAX · nscale)
    kT*: float64                       ## gas temperature (eV)
    van*: array[1..maxGases, float64]  ## density factor of each gas (VAN1..6)
    wb*: float64                       ## cyclotron frequency (rad/ps)
    eovb*: float64                     ## E/B (m/ps); 0 without B
    dec*: Decorrelation                ## pairs for the diffusion estimate

func corr(d: Deck): float64 =
  abzero * d.torr / (atmos * (abzero + d.tempc) * 100.0)

func setupError*(d: Deck): SetupError =
  if d.ngas == 0 or d.ngas > maxGases: return seGasInput
  var total = 0.0
  for j in 1 .. d.ngas:
    if d.gas[j] == 0 or d.percent[j] == 0.0: return seGasInput
    total = total + d.percent[j]
  if abs(total - 100.0) > 0.000001: return seGasInput
  if d.nmax < 0 or d.nmax > high(int) div nscale: return seNmaxOverflow
  seNone

func decorrelation(d: Deck): Decorrelation =
  ## long decorrelation length, or a short one for mixtures with more than
  ## 3 % molecular gas (not Ar, He-4, He-3, Ne, Kr, Xe) or a large field
  var molecular = 0.0
  for j in 1 .. d.ngas:
    if d.gas[j] notin 2 .. 7: molecular = molecular + d.percent[j]
  if d.efield > 10.0 / corr(d) or molecular > 3.0:
    Decorrelation(ncolm: 400000, ncorln: 50000, ncorst: 4)
  else:
    Decorrelation(ncolm: 2000000, ncorln: 500000, ncorst: 2)

func conditions*(d: Deck): Conditions =
  ## for a deck without `setupError`
  result.nmax = d.nmax * nscale
  result.kT = (abzero + d.tempc) * boltz
  let corr = corr(d)
  for j in 1 .. maxGases:
    result.van[j] = d.percent[j] * corr * const4 * 1000000000000000.0
  result.wb = awb * d.bmag * 1e-12
  if d.bmag != 0.0: result.eovb = d.efield * 1e-9 / d.bmag
  result.dec = decorrelation(d)
