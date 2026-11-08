## SETUP and SETUPT: the physical constants and what follows from a deck
## before any cross section is computed (gas densities, decorrelation of
## the diffusion estimate, number of collisions, cyclotron frequency).

import std/math
import units, deck, tables, drift

const
  nscale* {.intdefine.} = 40_000_000
    ## real collisions per unit of NMAX (the tests use less, see tests/regress)

  # constants of Magboltz (2010 values)
  echarg* = 1.602176565e-19.C   ## elementary charge
  emass* = 9.10938291e-31.kg    ## electron mass
  amu* = 1.660538921e-27.kg     ## atomic mass unit
  boltz* = 0.000086173324.eV•K⁻¹  ## Boltzmann constant
  awb = 17588200880.0.C•kg⁻¹•T•kGauss⁻¹
    ## e/m. Magboltz takes B in kG, so its number is e/m in C/kg times
    ## T/kG = 0.1: the unit says so.
  alosch* = 2.6867805e+19.cm⁻³  ## Loschmidt number
  abzero* = 273.15.K            ## 0 °C
  atmos* = 760.0.Torr           ## 1 atm

  # CNSTS1 of SETUP: factors of the free-flight formulas, from e/m with
  # Magboltz's literal factors as conversions of units (checked below)
  const1* = (awb / 2.0).to(cm²•V⁻¹•ps⁻²)           ## (e/m)/2
  const2* = const1.to(m•cm•V⁻¹•ps⁻²)               ## (e/m)/2
  const3* = (sqrt(2.0 * unity(T•kGauss⁻¹) * awb).to(cm•ps⁻¹•V⁻¹ᐟ²)) / sqrt(perElectron)
    ## sqrt(2 e/m): an electron's speed per √(energy). In Magboltz
    ## sqrt(0.2·AWB): 0.2 = 2 × 0.1 T/kG, which makes AWB e/m.
  const4 = const3 * alosch * 1e-15
    ## CONST4 of SETUP: speed of 1 eV times the Loschmidt number, times
    ## 1e-15 (a scale of Magboltz that `conditions` undoes; kept for the
    ## arithmetic of the densities, VAN)
  const5* = const3 / 2.0 * perElectron             ## sqrt(e/(2m)) · e

static:
  # the conversions above are Magboltz's literals
  doAssert const1.float == 17588200880.0 / 2.0 * 1e-19
  doAssert const2.float == 17588200880.0 / 2.0 * 1e-19 * 0.01
  doAssert const3.float == sqrt(0.2 * 17588200880.0) * 1e-9
  doAssert const5.float == sqrt(0.2 * 17588200880.0) * 1e-9 / 2.0

const
  naniso* = 2
    ## elastic angular distributions of the gas routines: 2 Okhrimovskyy
    ## (fixed in Magboltz; 1 Capitelli-Longo, 0 isotropic)
  tmax* = 100.0.ps            ## flight time covered by the null-collision bounds
  startTheta* = 0.785         ## start direction of the first electron (polar angle,
  startPhi* = 0.1             ## azimuth)

type
  SetupError* = enum
    seNone
    seGasInput       ## NGAS not in 1..6, a gas or share missing, shares not 100 %
    seNmaxOverflow   ## NMAX · nscale does not fit an integer

  Conditions* = object
    nmax*: int                         ## real collisions (NMAX · nscale)
    kT*: eV                            ## gas temperature
    van*: array[1..maxGases, cm⁻²•ps⁻¹•eV⁻¹ᐟ²]
      ## density of each gas times the speed of 1 eV (VAN1..6): times a cross
      ## section and √E, the collision frequency
    wb*: ps⁻¹                          ## cyclotron frequency (rad/ps)
    eovb*: m•ps⁻¹                      ## E/B; 0 without B
    dec*: Decorrelation                ## pairs for the diffusion estimate

func corr(d: Deck): UnitLess =
  ## gas density relative to 0 °C, 1 atm; per cent
  abzero * d.torr / (atmos * (abzero + d.tempc.K) * 100.0)

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
  if d.efield > 10.0.V•cm⁻¹ / corr(d) or molecular > 3.0:
    Decorrelation(ncolm: 400000, ncorln: 50000, ncorst: 4)
  else:
    Decorrelation(ncolm: 2000000, ncorln: 500000, ncorst: 2)

func conditions*(d: Deck): Conditions =
  ## for a deck without `setupError`
  result.nmax = d.nmax * nscale
  result.kT = (abzero + d.tempc.K) * boltz
  let corr = corr(d)
  for j in 1 .. maxGases:
    # (×1e15 undoes the scale of const4)
    result.van[j] = d.percent[j] * corr * const4 * 1000000000000000.0
  result.wb = (awb * d.bmag).to(ps⁻¹)
  if d.bmag != 0.0.kGauss:
    result.eovb = d.efield * unity(m•ps⁻¹•kGauss•cm•V⁻¹) / d.bmag  # E/B: (V/cm)/kG = 1e-9 m/ps
  result.dec = decorrelation(d)
