## The vocabulary of the gas modules (GAS1..GAS80 of Magboltz 11.21).
##
## A gas module turns tabulated measurements and a few analytic formulas
## into the cross sections of one gas on the energy grid (`GasXsec`). The
## Fortran routines repeat the same blocks hundreds of times: interpolate a
## table, fall off as a power of the energy above its end, use the
## Born-Bethe formulas at high energy, and copy the angular distribution of
## the elastic scattering to inelastic processes. Those blocks live here,
## each with the floating point operations in Fortran's order, so that the
## gases still match Magboltz bit for bit.
##
## Precision: Magboltz writes almost all its constants as REAL*4 literals
## (`11.548`, `0.0580`), which are rounded to single precision and then
## widened to double. The tables keep that type (`float32`), and so do the
## parameters below that always receive such literals; a few double
## constants (`1.0D-16`) are `float64`.

import std/math
import ../frt, ../units, ../mixture, ../tables
export math, units, mixture, tables, pow2

const
  # SETUP of Magboltz
  emass* = 9.10938291e-31       ## electron mass (kg)
  amu* = 1.660538921e-27        ## atomic mass unit (kg)
  pir2* = 8.7973554297e-17      ## π a₀² (cm²)
  ary* = 13.605692863464355     ## Rydberg energy (eV), REAL*4

let
  # common to the gas routines (computed at run time: Nim's compile time
  # evaluation does not round single precision values the same way)
  emass2* = float64(1021997.804'f32)   ## 2 m c² (eV)
  a0* = 0.52917720859e-8               ## Bohr radius (cm)
  ry* = float64(13.60569193'f32)       ## Rydberg energy (eV)
  bbconst* = 16.0 * arccos(-1.0) * a0 * a0 * ry * ry / emass2
    ## 16 π a₀² Ry² / (2 m c²), Born-Bethe excitation
  bornConst* = 1.873884e-20            ## Born-Bethe ionisation (cm²)

type
  GasInput* = object
    ## what a gas routine reads besides its own data
    e*: array[1..nEnergyBins, eV]  ## energy grid (bin centres)
    naniso*: int        ## 0 isotropic, 1 and 2 two models of anisotropy
    kT*: eV             ## gas temperature, for the populations of excited
                        ## states
    efinal*: eV         ## upper end of the grid
    estep*: eV          ## bin width; old gas routines step their own
                        ## energies by it
    penning*: bool      ## print the Penning transfers
    output*: File       ## to here (nil: nowhere)
    torr*: Torr         ## gas pressure and temperature (°C), for three-body
    tempc*: float64     ## attachment
    mixture*: array[1..6, tuple[gas: int, percent: float64]]
                        ## the gases of the mixture and their shares, for
                        ## gases that scale with their own density

  Relativistic* = object
    gamma2*, beta2*: float64

func mixtureOf*(gas: array[1..6, int], percent: array[1..6, float64]):
    array[1..6, tuple[gas: int, percent: float64]] =
  for k in 1 .. 6: result[k] = (gas[k], percent[k])

func percent*(inp: GasInput, gas: int): float64 =
  ## share of a gas in the mixture (%; the last slot holding it, as Magboltz
  ## looks it up)
  for m in inp.mixture:
    if m.gas == gas: result = m.percent

iterator bins*(inp: GasInput): (int, eV) =
  for i in 1 .. nEnergyBins: yield (i, inp.e[i])

iterator steps*(inp: GasInput): (int, eV) =
  ## the bins with energies summed up bin by bin, as the older gas routines
  ## do (EN = -ESTEP/2, EN = EN + ESTEP), rather than the grid's
  var en = -inp.estep / 2.0
  for i in 1 .. nEnergyBins:
    en = en + inp.estep
    yield (i, en)

proc relativistic*(en: float64): Relativistic =
  ## γ² and β² of an electron of kinetic energy en (eV)
  let gamma1 = (emass2 + 2.0 * en) / emass2
  result.gamma2 = gamma1 * gamma1
  let beta = sqrt(1.0 - 1.0 / result.gamma2)
  result.beta2 = beta * beta

proc vibrationalDipole*(c: float64, eloss, en: float64): float64 =
  ## Born approximation of a dipole active vibration (excitation or, with
  ## eloss < 0, de-excitation); c in units of 10⁻¹⁶ cm² eV
  let efac = sqrt(1.0 - eloss / en)
  if eloss < 0.0: c * ln((efac + 1.0) / (efac - 1.0)) / en
  else: c * ln((1.0 + efac) / (1.0 - efac)) / en

proc momentumFactor*(en, eloss: float64): float64 =
  ## momentum transfer of Born dipole scattering relative to its total:
  ## 1.5 - forward/(forward + backward)
  let elf = en - eloss
  let fwd = ln((en + elf) / (en + elf - 2.0 * sqrt(en * elf)))
  let bck = ln((en + elf + 2.0 * sqrt(en * elf)) / (en + elf))
  1.5 - fwd / (fwd + bck)

proc power*(x: float64, p: float64): float64 {.inline.} =
  ## x**p as gfortran evaluates it (integer powers by multiplication or
  ## libgcc's __powidf2)
  if p == 1.0: x
  elif p == 2.0: x * x
  elif p == trunc(p): powi(x, int(p))
  else: pow(x, p)

# ------------------------------------------------------------ interpolation

func segment[T](x: openArray[T], en: float64): int {.inline.} =
  ## upper end of the table interval that contains en (the last one above)
  for j in 1 .. x.high:
    if en <= x[j]: return j
  x.high

func lin*[X, Y](x: openArray[X], y: openArray[Y], en: float64): float64 =
  ## linear interpolation, extrapolating the last interval
  let j = segment(x, en)
  let (x0, x1) = (float64(x[j - 1]), float64(x[j]))
  let (y0, y1) = (float64(y[j - 1]), float64(y[j]))
  let a = (y1 - y0) / (x1 - x0)
  let b = (x0 * y1 - x1 * y0) / (x0 - x1)
  a * en + b

proc lin*[X, Y](x: openArray[X], y: openArray[Y], en: float64,
                falloff: float64): float64 =
  ## linear interpolation; above the table y falls off as en^-falloff
  if en > float64(x[^1]): float64(y[^1]) * power(float64(x[^1]) / en, falloff)
  else: lin(x, y, en)

func logLog*[X, Y](x: openArray[X], y: openArray[Y], en: float64): float64 =
  ## interpolation linear in log y over log x
  let j = segment(x, en)
  let yj = ln(float64(y[j]))
  let yj1 = ln(float64(y[j - 1]))
  let xj = ln(float64(x[j]))
  let xj1 = ln(float64(x[j - 1]))
  let a = (yj - yj1) / (xj - xj1)
  let b = (xj1 * yj - xj * yj1) / (xj1 - xj)
  exp(a * ln(en) + b)

func logLin*[X, Y](x: openArray[X], y: openArray[Y], en: float64): float64 =
  ## interpolation linear in log y over x (bremsstrahlung tables)
  let j = segment(x, en)
  let (x0, x1) = (float64(x[j - 1]), float64(x[j]))
  let (y0, y1) = (ln(float64(y[j - 1])), ln(float64(y[j])))
  let a = (y1 - y0) / (x1 - x0)
  let b = (y1 * x0 - y0 * x1) / (x0 - x1)
  exp(a * en + b)

type
  Negative* = enum
    ## what a gas does with Born-Bethe values below zero (near threshold)
    keepNegative, clampNegative, absNegative

func offsets*(inp: GasInput, e: openArray[eV], magnitude = false): seq[int] =
  ## IOFF of each energy loss: number of grid bins at or below it (or its
  ## magnitude). Bin i - IOFF is where an electron of bin i ends up after
  ## the loss.
  result = newSeq[int](e.len)
  for k, ek in e:
    let ek = if magnitude: abs(ek) else: ek
    for i in 1 .. nEnergyBins:
      if inp.e[i] > ek:
        result[k] = i - 1
        break

# ------------------------------------------------------------ gas data

proc elastic*(g: var GasXsec, inp: GasInput, i: int, qela, qmom: cm², pq2: float64) =
  ## elastic cross section of bin i from the total elastic (qela) and the
  ## momentum transfer (qmom) cross sections. Isotropic: the momentum
  ## transfer cross section. Model 1: forward probability
  ## 0.5 + (qela - qmom) / qela; model 2: pq2.
  let pq1 = 0.5 + float64((qela - qmom) / qela)
  case inp.naniso
  of 0: g.peqel[2, i] = 0.5
  of 1: g.peqel[2, i] = pq1
  else: g.peqel[2, i] = pq2
  g.q[2, i] = if inp.naniso == 0: qmom else: qela

proc describe*(g: var GasXsec, inp: GasInput, isotropic, anisotropic: string,
               mass: float32) =
  ## name (isotropic and anisotropic variant) and mass (amu); elastic and
  ## ionisation scattering anisotropic as chosen
  g.name = if inp.naniso == 0: isotropic else: anisotropic
  g.massRatio = 2.0 * emass / (float64(mass) * amu)
  for k in 1 .. 6: g.kel[k] = inp.naniso

type
  Shell* = tuple
    ## an ionisation channel: energy, Opal-Beaty energy, Auger electrons
    ## and their energy, description
    eion, eoby: float32
    nc0: int
    ec0: float32
    text: string

proc shells*(g: var GasXsec, s: openArray[Shell]) =
  g.nion = s.len
  for k, sh in s:
    g.eion[k + 1] = float64(sh.eion).eV
    g.eb[k + 1] = float64(sh.eoby).eV
    g.nc0[k + 1] = sh.nc0
    g.ec0[k + 1] = float64(sh.ec0).eV
    g.scrp[2 + k + 1] = sh.text

proc attachment*(g: var GasXsec, texts: openArray[string]) =
  g.natt = texts.len
  for k, t in texts: g.scrp[2 + g.nion + k + 1] = t

proc levels*(g: var GasXsec, inp: GasInput, l: openArray[(float32, string)]) =
  ## inelastic levels: energy loss (eV, < 0 superelastic), description;
  ## anisotropy as the elastic scattering
  g.nin = l.len
  for j, (e, text) in l:
    g.ei[j + 1] = float64(e).eV
    g.kin[j + 1] = inp.naniso
    g.scrp[4 + g.nion + g.natt + j + 1] = text

proc penning*(g: var GasXsec, inp: GasInput, levels: Slice[int],
              fraction, length, time: float32, fmt: string, all = false) =
  ## Penning transfer from levels: probability, absorption length (µm) and
  ## decay time (ps); printed in Magboltz's format when asked for (only
  ## levels with a transfer, or `all`)
  for j in levels:
    g.penningFraction[j] = fraction
    g.penningRange[j] = float64(length).μm
    g.penningTime[j] = float64(time).ps
  if inp.penning:
    for j in levels:
      if all or g.penningFraction[j] != 0.0:
        inp.output.fwrite(fmt, fstr(g.name, 25), g.ei[j], g.penningFraction[j],
               g.penningRange[j], g.penningTime[j])

proc dropLevelsAbove*(g: var GasXsec, efinal: eV, first = 1) =
  ## levels from the first one at or above the end of the grid are not used
  ## (Magboltz looks from level `first`)
  for k in first .. g.nin:
    if efinal <= g.ei[k]:
      g.nin = k - 1
      return
