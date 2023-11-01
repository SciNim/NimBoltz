## Input cards straight from `magboltz-11.17.f`
## C---------------------------------------------------------------
## C   INPUT CARDS :
## C----------------------------------------------------------
## C  FIRST CARD: 4I10,F10.5  :  NGAS,NMAX,IPEN,ITHRM,EFINAL
## C  NGAS:  NUMBER OF GASES IN MIXTURE
## C  NMAX: NUMBER OF REAL COLLISIONS ( MULTIPLE OF 1*10**7 )
## C  USE NMAX = BETWEEN 2 AND 5 FOR INELASTIC GAS TO OBTAIN 1% ACCURACY
## C      NMAX = ABOVE 10 FOR BETTER THAN 0.5% ACCURACY.
## C      NMAX = AT LEAST 10 FOR PURE ELASTIC GASES LIKE ARGON
## C      HIGHER VALUES THAN NMAX=214 CAN ONLY BE USED ON COMPUTERS SUCH
## C      AS DEC ALPHAS WITH TRUE 64 BIT INTEGERS. PCS ARE LIMITED TO
## C      31 BIT INTEGERS...
## C    IPEN   = 0 PENNING EFFECTS NOT INCLUDED
## C           = 1 PENNING EFFECTS INCLUDED (SEE INSTRUCTIONS ABOVE)
## C    ITHRM  = 0 GAS MOTION ASSUMED TO BE AT O KELVIN (STATIC GAS)
## C    ITHRM  = 1 GAS MOTION TAKEN TO BE AT INPUT TEMPERATURE
## C    EFINAL = UPPER LIMIT OF THE ELECTRON ENERGY IN ELECTRON VOLTS.
## C    EFINAL = 0.0 (PROGRAM AUTOMATICALLY CALCULATES UPPER INTEGRATION
## C                 ENERGY LIMIT)
## C-------------------------------------------------------------
## C  SECOND CARD : 6I5   : NGAS1 , NGAS2, NGAS3 , NGAS4 , NGAS5 , NGAS6
## C       NGAS1,ETC :  GAS NUMBER IDENTIFIERS (BETWEEN 1 AND 80)
## C                   SEE GAS LIST BELOW FOR IDENTIFYING NUMBERS.
## C
## C-------------------------------------------------------------
## C THIRD CARD: 8F10.4  : FRAC1,FRAC2,FRAC3,FRAC4,FRAC5,FRAC6,TEMP,TORR
## C  FRAC1,ETC : PERCENTAGE FRACTION OF GAS1,ETC
## C  TEMP : TEMPERATURE OF GAS IN CENTIGRADE
## C  TORR :  PRESSURE OF GAS IN TORR
## C ------------------------------------------------------------
## C FOURTH CARD : 6F10.3  : EFIELD,BMAG,BTHETA
## C  EFIELD : ELECTRIC FIELD IN VOLTS/ CM.
## C   BMAG  : MAGNITUDE OF THE MAGNETIC FIELD IN KILOGAUSS
## C  BTHETA : ANGLE BETWEEN THE ELECTRIC AND MAGNETIC FIELDS IN DEGREES.
## C-----------------------------------------------------------------------
## C CARD 4*N+1 USES NGAS=0 TO TERMINATE CORRECTLY
## C--------------------------------------------------------------------


## Gas numbers
## Routine Gas     Last update     Notes   Star rating
## GAS1    CF4     2001    anisotropic scattering  5*
## GAS2    Argon   1997            5*
## GAS3    Helium 4        1997            5*
## GAS4    Helium 3        1992            5*
## GAS5    Neon    1992            5*
## GAS6    Krypton 2001            4*
## GAS7    Xenon   2001            4*
## GAS8    methane 1994            5*
## GAS9    ethane  1999            5*
## GAS10   propane 1999            4*
## GAS11   isobutane       1999            3*
## GAS12   CO2     2001            5*
## GAS13   C(CH3)4 neo-pentane     1995            3*
## GAS14   H20     1998            3*
## GAS15   Oxygen  1990    3-body attachment included      4*
## GAS16   Nitrogen        1985    Pitchford and Phelps    4*
## GAS17   nitric oxide    1995    attaching gas   4*
## GAS18   nitrous oxide   1995    attaching gas   4*
## GAS19   C2H4 ethene     1999            4*
## GAS20   C2H2 acetylene  1992            3*
## GAS21   Hydrogen        2001            5*
## GAS22   Deuterium       1998            5*
## GAS23   Carbon monoxide 1998            5*
## GAS24   methylal        1988            2*
## GAS25   DME     1998            4*
## GAS26   Reid step model ?       anisotropic version     -
## GAS27   Maxwell model   ?               -
## GAS28   Reid ramp model ?               -
## GAS29   C2F6    1999    anisotropic     4*
## GAS30   SF6     ?       do not use high percentage      3*
## GAS31   NH3 ammonia     1999            3*
## GAS32   C3H6 propene    1999            4*
## GAS33   C3H6 cyclopropane       1999            4*
## GAS34   CH3OH methanol  1999            2*
## GAS35   C2H5OH ethanol  1999            3*
## GAS36   C3H7OH isopropanol      1999            2*
## GAS37   Cæsium  2001    no dimers       2*
## GAS38   Fluorine        ?       Morgan  2*
## GAS39   CS2     2001    ion drift, dark matter  2*
## GAS40   COS     2001            2*
## GAS41   CD4     2001    TPCs in neutron background environment  3*
## GAS42   BF3 Boron trifluoride   2001    anisotropic     3*
## GAS43   C2HF5 or C2H2F4 ?       estimated no data, anisotropic  1*
## GAS44   Helium 3        2002    anisotropic     5*
## GAS45   Helium 4        2002    anisotropic     5*
## GAS46   Neon    2002    anisotropic     5*
## GAS47   Argon   2002    anisotropic     5*
## GAS48   Krypton 2002    anisotropic     4*
## GAS49   Xenon   2002    anisotropic     4*
## GAS50   methane 2002    anisotropic     5*
## GAS52-80        Dummy routines  ?               -

import std / [math, strutils, strformat, sequtils, parseutils]

import parse_gases

import unchained, measuremancer, npeg

type
  NotImplementedError = object of Defect

defUnit(kV•cm⁻¹)
defUnit(V•cm⁻¹)
type
  Gas* = object
    name*: string
    id*: int
    frac*: float # fraction of total pressure

  Magboltz* = object
    gases*: seq[Gas] ## Each gas
    fractions*: seq[float]
    pressure*: MilliBar ## Pressure of the gas in MilliBar (will be converted to Torr later)
    temp*: Kelvin ## Temperature of the gas (will be converted to centigrade later)
    E*: V•cm⁻¹ ## Field strength
    B*: Tesla  ## Magnetic field in Tesla (will be converted to kilo Gauss later)
    Bθ*: Degree = 0.° ## Angle between electric and magnetic field
    nMax*: int = 10 ## maximum number of collisions to simulate as multiple of 1e7
                   ## MagBoltz documentation:
                   ## NMAX: number of real collisions (multiple of 10**7), use a value between 2 and 5
                   ## for inelastic gas to obtain 1 % accuracy, use a value above 10 for better than
                   ## 0.5 % accuracy and a value of at least 10 for pure elastic gases like Argon;
    eFinal*: float  ## MagBoltz:
                   ## EFINAL: upper limit of the electron energy in electron Volts, if EFINAL = 0.0,
                   ## program automatically calculates upper integration energy limit.
    usePenning*: bool = true ## Whether to include the Penning effect or not
    gasMotionThermal*: bool = true ## Whether gas motion is to be taken at 0 K or input temp
    outpath*: string ## The path in which to store the magboltz input / output files


## Note on diffusion `constant` vs diffusion `coefficient`.
## After a drift time `t` with a diffusion coefficient `D`
## a cloud of electrons will be a distance `σ_T` away from
## the origin
## `σ_T = √(2 N D t)`
## in `N` dimensions.
## In gaseous detector physics `N` is usually taken to be 1
## and each dimension is treated seperately.
##
## It is convenient to introduce a diffusion `constant` D_t such that
## `σ_T = D_t * √x`
## where `x` is the drift length.
## We can therefore relate the diffusion constant to the diffusion
## coefficient by using the drift velocity `v`:
## `D_t = √(2 N D / v)`
##
## Typically the diffusion coefficient `D` is given in `cm²•s⁻¹` and
## the diffusion constant `D_t` in `μm•√cm⁻¹`. The diffusion constant
## therefore says "how much diffusion in μm is expected after a drift
## of `√cm` length?".
##
## In the code below we do not give units to `D_t` (called `Dt_L` and `Dt_T`
## for the 1D transverse and longitudinal components), because `unchained`
## currently does not support square roots of units unfortunately.

defUnit(μm•ns⁻¹)
defUnit(cm²•s⁻¹)
defUnit(cm⁻¹)
type # These come from using `parseMagboltz`
  MagboltzResult* = object
    xDrift*: Measurement[μm•ns⁻¹] ## Drift velocity in `x` direction
    yDrift*: Measurement[μm•ns⁻¹] ## Drift velocity in `y` direction
    zDrift*: Measurement[μm•ns⁻¹] ## Drift velocity in `z` direction. This is the direction of the electric field and therefore drift direction!
    Dt_L*: Measurement[float] ## The longitudinal diffusion *constant* `D_t` in units of `μm•√cm⁻¹`
    Dt_T*: Measurement[float] ## The longitudinal diffusion *constant* `D_t` in units of `μm•√cm⁻¹`
    D_L*: Measurement[cm²•s⁻¹] ## The longitudinal diffusion *coefficient* in units of `cm²•s⁻¹`
    D_T*: Measurement[cm²•s⁻¹] ## The transverse diffusion *coefficient* in units of `cm²•s⁻¹`
    sst*: SteadyStateResult ## May be empty, only for fields strong enough for amplification

  SteadyStateResult* = object
    vDrift*: Measurement[μm•ns⁻¹]
    ws*: Measurement[float] ## XXX: w hat is `WS`?
    D_L*: Measurement[cm²•s⁻¹] ## The longitudinal diffusion *coefficient* in units of `cm²•s⁻¹`
    D_T*: Measurement[cm²•s⁻¹] ## The transverse diffusion *coefficient* in units of `cm²•s⁻¹`
    α*: Measurement[cm⁻¹] ## First Townsend coefficient (Gas gain = `e^{αx}`, `x` drift length). Usually given in `cm⁻¹`.
    att*: Measurement[float] ## The attachment rate I think?

proc initGas*(id: int, fraction: float): Gas =
  let gM = getGas(id)
  result = Gas(name: gM.name, id: gM.id, frac: fraction)

proc initGas*(name: string, fraction: float): Gas =
  let gM = getGas(name)
  result = Gas(name: gM.name, id: gM.id, frac: fraction)

proc initGases*[T: int | string](gases: varargs[(T, float)]): seq[Gas] =
  for arg in gases:
    result.add initGas(arg[0], arg[1])

proc initMagboltz*[F: SomeUnit, P: Pressure](
  gases: seq[Gas], pressure: P, temp: Kelvin, E: F,
  nMax: int,
  usePenning = true, gasMotionThermal = true,
  B: Tesla = 0.T,
  Bθ: Degree = 0.°,
  eFinal = 0.0,
  outpath = "resources"
                                          ): Magboltz =
  var f = gases.mapIt(it.frac)
  let eField = E.to(V•cm⁻¹)
  if gases.len > 6:
    raise newException(ValueError, "Magboltz does not support more than 6 gases in a mixture.")
  f.setLen(6) ## Set to 6, fill with zeros
  result = Magboltz(gases: gases,
                    fractions: f,
                    pressure: pressure.to(MilliBar),
                    temp: temp,
                    E: eField,
                    B: B,
                    Bθ: Bθ,
                    nMax: nMax,
                    eFinal: eFinal,
                    usePenning: usePenning,
                    gasMotionThermal: gasMotionThermal,
                    outpath: outpath)

proc handleInteger[T: SomeNumber | bool](width: int, arg: T): string =
  let argI = when T is SomeFloat: arg.round.int
             elif T is SomeInteger: arg
             else:
               if arg: 1 else: 0
  let strVal = $argI
  if strVal.len > width:
    raise newException(ValueError, "Integer argument: " & strVal & " is too long to fit width " & $width & ".")
  result = alignLeft(strVal, width)

proc handleFloat[T: SomeNumber | bool](widthPrecision: string, arg: T): string =
  ## Input must be of the form `<integer>.<integer>` for `widthPrecision`
  let argF = when T is SomeNumber: arg.float else:
    if arg: 1.0 else: 0.0
  doAssert "." in widthPrecision
  let wp = widthPrecision.split(".")
  doAssert wp.len == 2
  let (width, precision) = (parseInt wp[0], parseInt wp[1])
  let strVal = formatBiggestFloat(arg.float, ffDecimal, precision)
  if strVal.len > width:
    raise newException(ValueError, "Float argument: " & strVal & " is too long to fit width " & $width & ".")
  result = alignLeft(strVal, width)

proc handleExp(): string = doAssert false, "Exp format 'E' not implemented"
proc handleDouble(): string = doAssert false, "Double format 'D' not implemented"

proc formatFortran[T: SomeNumber | bool](f: string, arg: T): string =
  case f[0] # the Fortran format specifier
  of 'I': handleInteger(parseInt f[1 .. ^1], arg)
  of 'F': handleFloat(f[1 .. ^1], arg)
  of 'E': handleExp()
  of 'D': handleDouble()
  else:
    raise newException(NotImplementedError, "The format specifier " & $f[0] & " is not implemented yet.")

iterator iterFormat(fmt: string): string =
  ## Yields all the individual components of the format string
  # 1. split all formats by `,`
  let pieces = fmt.split(",")
  for p in pieces:
    # 2. for each piece parse until the first non digit character
    var numStr = ""
    let idx = parseWhile(p, numStr, {'0' .. '9'}, 0)
    echo "IDX?? ", idx, " for : ", p, " num ", numStr
    if numStr.len > 0: # more than 1 element with same format
      let num = parseInt numStr
      # 3. yield each element
      for i in 0 ..< num:
        yield p[idx .. ^1]
    else: # single element of this format, yield
      yield p

proc formatCard1(): string = "4I10,F10.5"
proc formatCard2(): string = "6I5"
proc formatCard3(): string = "8F10.4"
proc formatCard4(): string = "3F10.3" ## IMPORTANT: The Magboltz documentation here https://magboltz.web.cern.ch/magboltz/usage.html writes
                                      ## "6F10.3" but that seems clearly wrong
                                      ## This is `also` still the case in the current source file of Magboltz!

proc writeCard1(f: var File, nGas: int, nMax: int, eFinal: float, usePenning, gasMotionThermal: bool) =
  ## MagBoltz documentation:
  ## First card
  ## C  NGAS:  NUMBER OF GASES IN MIXTURE
  ## C  NMAX: NUMBER OF REAL COLLISIONS ( MULTIPLE OF 1*10**7 )
  ## C  USE NMAX = BETWEEN 2 AND 5 FOR INELASTIC GAS TO OBTAIN 1% ACCURACY
  ## C      NMAX = ABOVE 10 FOR BETTER THAN 0.5% ACCURACY.
  ## C      NMAX = AT LEAST 10 FOR PURE ELASTIC GASES LIKE ARGON
  ## C      HIGHER VALUES THAN NMAX=214 CAN ONLY BE USED ON COMPUTERS SUCH
  ## C      AS DEC ALPHAS WITH TRUE 64 BIT INTEGERS. PCS ARE LIMITED TO
  ## C      31 BIT INTEGERS...
  ## C    IPEN   = 0 PENNING EFFECTS NOT INCLUDED
  ## C           = 1 PENNING EFFECTS INCLUDED (SEE INSTRUCTIONS ABOVE)
  ## C    ITHRM  = 0 GAS MOTION ASSUMED TO BE AT O KELVIN (STATIC GAS)
  ## C    ITHRM  = 1 GAS MOTION TAKEN TO BE AT INPUT TEMPERATURE
  ## C    EFINAL = UPPER LIMIT OF THE ELECTRON ENERGY IN ELECTRON VOLTS.
  ## C    EFINAL = 0.0 (PROGRAM AUTOMATICALLY CALCULATES UPPER INTEGRATION
  ## C                 ENERGY LIMIT)
  let fmt = formatCard1()
  var idx = 0
  ## Note: we could of course hardcode the format into the code, but this way it is somewhat
  ## generic over the inputs.
  for el in iterFormat(fmt):
    echo "Format: ", el
    case idx
    of 0: f.write(formatFortran(el, nGas))
    of 1: f.write(formatFortran(el, nMax))
    of 2: f.write(formatFortran(el, usePenning))
    of 3: f.write(formatFortran(el, gasMotionThermal))
    of 4: f.write(formatFortran(el, eFinal))
    else: doAssert false, "Unexpected format!"
    inc idx
  f.write("\n")

proc writeCard2(f: var File, gases: seq[Gas]) =
  ## MagBoltz documentation:
  ## Second card
  ##   Format: 6I5, variables: NGAS1,NGAS2,NGAS3,NGAS4,NGAS5,NGAS6
  ##   NGAS1 etc.: gas number identifiers (between 1 and 80) see gas list below for identifying numbers.
  let fmt = formatCard2()
  var idx = 0
  ## Note: we could of course hardcode the format into the code, but this way it is somewhat
  ## generic over the inputs.
  for el in iterFormat(fmt):
    if idx < gases.len:
      f.write(formatFortran(el, gases[idx].id))
    else:
      f.write(formatFortran(el, 80)) ## Weird placeholder
    inc idx
  f.write("\n")

proc writeCard3(f: var File, fractions: seq[float], temp: Kelvin, pressure: MilliBar) =
  ## MagBoltz documentation:
  ## Third card
  ##   Format: 8F10.4, variables: FRAC1,FRAC2,FRAC3,FRAC4,FRAC5,FRAC6,TEMP,TORR
  ##   FRAC1 etc.: percentage fraction of gas1 etc.;
  ##   TEMP: temperature of gas in centigrade;
  ##   TORR: pressure of gas in Torr.
  let fmt = formatCard3()
  var idx = 0
  ## Note: we could of course hardcode the format into the code, but this way it is somewhat
  ## generic over the inputs.
  for el in iterFormat(fmt):
    case idx
    of 0 .. 5:
      doAssert fractions.len == 6
      f.write(formatFortran(el, fractions[idx]))
    of 6:
      f.write(formatFortran(el, temp.float - 273.15)) ## convert to centigrade by hand (unchained does not know °C)
    of 7:
      f.write(formatFortran(el, pressure.to(Torr).float)) ## convert to Torr
    else: doAssert false, "Unexpected format!"
    inc idx
  f.write("\n")

proc writeCard4(f: var File, E: V•cm⁻¹, B: Tesla, Bθ: Degree) =
  ## MagBoltz documentation:
  ## Fourth card
  ##   Format: 6F10.3, variables: EFIELD,BMAG,BTHETA
  ##   EFIELD: electric field in Volt/cm;
  ##   BMAG: magnitude of the magnetic field in kilogauss;
  ##   BTHETA: angle between the electric and magnetic fields in degrees.
  let fmt = formatCard4()
  var idx = 0
  ## Note: we could of course hardcode the format into the code, but this way it is somewhat
  ## generic over the inputs.
  for el in iterFormat(fmt):
    case idx
    of 0: f.write(formatFortran(el, E.to(V•cm⁻¹).float))
    of 1: f.write(formatFortran(el, B.to(KiloGauss).float))
    of 2: f.write(formatFortran(el, Bθ.float))
    else: doAssert false, "Unexpected format!"
    inc idx
  f.write("\n")

proc finishFile(f: var File) =
  ## NOTE: It seems modern Magboltz tries to read the full line to check for
  ## `NGAS=0`, see lines 2740-2743:
  ##
  ##      READ(5,*) NGAS,NMAX,IPEN,ITHRM,EFINAL
  ##C     READ(5,2) NGAS,NMAX,IPEN,ITHRM,EFINAL
  ##C   2 FORMAT(4I10,F10.5)
  ##      IF(NGAS.EQ.0) GO TO 99
  ##
  ## Therefore we emit 5 elements regardless of the fact that they are all 0
  f.write(formatFortran("I10", 0)) # add `NGAS=0`
  f.write(formatFortran("I10", 0))
  f.write(formatFortran("I10", 0))
  f.write(formatFortran("I10", 0))
  f.write(formatFortran("F10.5", 0.0))

proc writeFile(m: Magboltz, fname: string) =
  var f = open(fname, fmWrite)
  f.writeCard1(m.gases.len, m.nMax, m.eFinal, m.usePenning, m.gasMotionThermal)
  f.writeCard2(m.gases)
  f.writeCard3(m.fractions, m.temp, m.pressure)
  f.writeCard4(m.E, m.B, m.Bθ)
  f.finishFile()
  f.close()

proc genSettingsStr(mb: Magboltz): string =
  for g in mb.gases:
    result.add &"{g.name}_{g.frac}"
  result.add &"T_{mb.temp.float}_E_{mb.E.float}_P_{mb.pressure.float}_usePenning_{mb.usePenning}_gasTh_{mb.gasMotionThermal}_nMax_{mb.nMax}"

from std / os import createDir, `/`
proc genInfile(mb: Magboltz): string =
  createDir(mb.outpath)
  var infile = "magboltz_input_" & genSettingsStr(mb) & ".dat"
  result = mb.outpath / infile

proc genOutfile(mb: Magboltz): string =
  createDir(mb.outpath)
  var outfile = "magboltz_output_" & genSettingsStr(mb) & ".dat"
  result = mb.outpath / outfile

import pkg / shell
proc runMagboltz*(mb: Magboltz) =
  let infile = genInfile(mb)
  writeFile(mb, infile)
  let outfile = genOutfile(mb)
  echo "Run:"
  let (res, err) = shellVerbose:
    magboltz "<" ($infile)
  writeFile(outfile, res)

proc runMagboltz*(mbs: openArray[Magboltz]) =
  ## As we don't need to return anything we just create one thread for each task
  proc runCommand(mb: Magboltz) {.thread.} =
    runMagboltz(mb)

  var thrs = newSeq[Thread[Magboltz]](mbs.len)
  for i in 0 ..< mbs.len:
    createThread(thrs[i], runCommand, mbs[i])
  joinThreads(thrs)

proc test() =
  doAssert formatFortran("I6", 123) == &"{123:>6}"
  doAssert formatFortran("I6", 4) == &"{4:>6}"
  try:
    discard formatFortran("I2", 100)
  except ValueError:
    discard
  try:
    discard formatFortran("I2", 100)
  except ValueError:
    discard

  doAssert formatFortran("F6.2", 1.23221) == &"{1.23221:6.2f}"
  try:
    discard formatFortran("F62", 1.23221)
  except AssertionError:
    discard

proc run(usePenning = true, gasMotionThermal = true, drift = false, amp = false) =
  let gases = initGases([("argon", 97.7), ("isobutane", 2.3)])
  if amp:
    let mb = initMagboltz(gases, 1050.mbar, 300.K, 60.kV•cm⁻¹, nMax = 1,
                          usePenning = usePenning, gasMotionThermal = gasMotionThermal)
    # to run a single calculation:
    # runMagboltz(mb)
    # to run multiple in parallel:
    var mb1 = mb
    var mb2 = mb
    var mb3 = mb
    var mb4 = mb
    mb1.temp = 300.K
    mb2.temp = 330.K
    mb3.temp = 360.K
    mb4.temp = 390.K
    let mbs = [mb1, mb2, mb3, mb4]
    runMagboltz(mbs)
  if drift:
    let mb = initMagboltz(gases, 1050.mbar, 300.K, 500.V•cm⁻¹, nMax = 1,
                          usePenning = usePenning, gasMotionThermal = gasMotionThermal)
    # to run a single calculation:
    # runMagboltz(mb)
    # to run multiple in parallel:
    var mb1 = mb
    var mb2 = mb
    var mb3 = mb
    var mb4 = mb
    mb1.temp = 300.K
    mb2.temp = 330.K
    mb3.temp = 360.K
    mb4.temp = 390.K
    let mbs = [mb1, mb2, mb3, mb4]
    runMagboltz(mbs)


## Constants to help with extracting the parts of the Magboltz result files
## that is interesting for us (at the moment)
const sstStart = "SOLUTION FOR STEADY STATE TOWNSEND PARAMETERS"
const sstStop = "SOLUTION FOR PULSED TOWNSEND AND TIME OF FLIGHT PARAMETERS"

const commonStart = "CALCULATED MAX. COLLISION TIME"
const commonStop = "MEAN ELECTRON ENERGY"

## NOTE: Instead of writing a simple parser I ended up using `npeg` instead
## for the parsing...
when false:
  proc parseSteadyState(s: string): SteadyStateResult =
    let idxStart = dat.find(sstStart)
    let idxStop = dat.find(sstStop)
    let lines = dat[idxStart ..< idxStop].splitLines
    for l in lines:
      ## Parse:
      ## ```
      ## SST DRIFT VELOCITIES
      ##
      ##  VD=    257.5 +-   0.29 %   WS=    292.1 +-   1.11 %
      ##
      ##  SST DIFFUSION
      ##
      ##  DL=   2956.7 +-    1.8 %   DT=   3713.0 +-   4.36 %
      ##
      ##  SST TOWNSEND COEFICIENTS
      ##
      ##  ALPHA =   1495.1 +-   0.50 %    ATT=      0.0 +-   0.00 %
      ## ```
      let lS = l.strip
      if lS.startsWith("VD"): discard
        # parse VD and WS
      elif lS.startsWith("DL"): discard
      elif lS.startsWith("ALPHA"): discard

  proc parseDrift(s: string): μm•ns⁻¹ =
    discard

  proc parseCommonFields(s: string): MagboltzResult =
    let idxStart = dat.find(commonStart)
    let idxStop = dat.find(commonStop)
    let lines = dat[idxStart ..< idxStop].splitLines
    for l in lines:
      ## Parsing:
      ## ```
      ##   Z DRIFT VELOCITY = 0.1833E+03 MICRONS/NANOSECOND  +-    0.14%
      ##   Y DRIFT VELOCITY = 0.0000E+00 MICRONS/NANOSECOND  +-    0.00%
      ##   X DRIFT VELOCITY = 0.0000E+00 MICRONS/NANOSECOND  +-    0.00%
      ##
      ##
      ##            DIFFUSION IN CM**2/SEC.
      ##
      ##
      ##   TRANSVERSE DIFFUSION   = 0.2477D+04 +-   10.16%
      ##           =  8.10676 EV. +-  10.163%
      ##           =  164.385 MICRONS/CENTIMETER**0.5  +-    5.08%
      ##
      ##
      ##   LONGITUDINAL DIFFUSION = 0.1761D+04 +-    10.8%
      ##           =   5.7649 EV. +-   10.79%
      ##           =  138.622 MICRONS/CENTIMETER**0.5  +-    5.40%
      ## ```
      let lS = l.strip
      if "DRIFT VELOCITY" in lS:
        if lS.startsWith("Z"):
          result.zDrift = parseDrift(lS)

defUnit(μm²•cm⁻¹)
defUnit(m²•s⁻¹)
defUnit(m•s⁻¹)

proc parseMagboltzNpeg(input: string): MagboltzResult =
  ## Grammar definition.
  ## This is the first time I'm using `npeg` or PEGs in general, so come at me. :P
  ##
  ## Parses:
  ## ```
  ##   Z DRIFT VELOCITY = 0.1833E+03 MICRONS/NANOSECOND  +-    0.14%
  ##   Y DRIFT VELOCITY = 0.0000E+00 MICRONS/NANOSECOND  +-    0.00%
  ##   X DRIFT VELOCITY = 0.0000E+00 MICRONS/NANOSECOND  +-    0.00%
  ##
  ##
  ##            DIFFUSION IN CM**2/SEC.
  ##
  ##
  ##   TRANSVERSE DIFFUSION   = 0.2477D+04 +-   10.16%
  ##           =  8.10676 EV. +-  10.163%
  ##           =  164.385 MICRONS/CENTIMETER**0.5  +-    5.08%
  ##
  ##
  ##   LONGITUDINAL DIFFUSION = 0.1761D+04 +-    10.8%
  ##           =   5.7649 EV. +-   10.79%
  ##           =  138.622 MICRONS/CENTIMETER**0.5  +-    5.40%
  ## ...
  ## SST DRIFT VELOCITIES
  ##
  ##  VD=    257.5 +-   0.29 %   WS=    292.1 +-   1.11 %
  ##
  ##  SST DIFFUSION
  ##
  ##  DL=   2956.7 +-    1.8 %   DT=   3713.0 +-   4.36 %
  ##
  ##  SST TOWNSEND COEFICIENTS
  ##
  ##  ALPHA =   1495.1 +-   0.50 %    ATT=      0.0 +-   0.00 %
  ## ```

  let parser = peg("root", mb: MagboltzResult):
    space <- *{' '}
    eol <- '\n'
    expNumber <- +Digit * ?'.' * *Digit * (i"E" | i"D") * ?{'+', '-'} * +{'0'..'9'}
    number <- +Digit * ?'.' * *Digit
    percentage <- >number * '%'
    junk <- *(' ' | Print) * +eol
    dataLine <- space * >+Alpha * ?space * "=" * space * >number * space * "+-" * space * >number * space * '%' * space * ?eol:
      let val = parseFloat($2)
      let err = parseFloat($3)
      let merr = val * (err / 100.0)
      case $1
      of "VD": mb.sst.vDrift = val.μm•ns⁻¹ ± merr.μm•ns⁻¹
      of "WS": mb.sst.ws = val ± merr
      of "DL": mb.sst.D_L = val.cm²•s⁻¹ ± merr.cm²•s⁻¹
      of "DT": mb.sst.D_T = val.cm²•s⁻¹ ± merr.cm²•s⁻¹
      of "ALPHA": mb.sst.α = val.cm⁻¹ ± merr.cm⁻¹
      of "ATT": mb.sst.att = val ± merr
      else: doAssert false, "Encountered " & ($1)
    velocity <- space * >+Alpha * space * "DRIFT VELOCITY =" * space * >expNumber * space * "MICRONS/NANOSECOND" * space * "+-" * space * >percentage * *eol:
      let val = parseFloat($2)
      let err = parseFloat(($3)[0 .. ^2])
      let mval = val.μm•ns⁻¹ ± (val * err / 100.0).μm•ns⁻¹
      case $1
      of "Z": mb.zDrift = mval
      of "Y": mb.yDrift = mval
      of "X": mb.xDrift = mval
      else: doAssert false, "Encountered " & ($1)
    diffusion <- space * >+Alpha * space * "DIFFUSION" * space * "=" * space * >expNumber * space * "+-" * space * >percentage * eol:
      let val = parseFloat(($2).replace("D", "E"))
      let err = parseFloat(($3)[0 .. ^2])
      let mval = val.cm²•s⁻¹ ± (val * err / 100.0).cm²•s⁻¹
      case $1
      of "TRANSVERSE":   mb.D_L = mval
      of "LONGITUDINAL": mb.D_T = mval
      else: doAssert false, "Encountered " & ($1)
    root <- *(velocity | diffusion | dataLine | junk)

  # And now parse the input file
  doAssert parser.match(input, result).ok

  # now compute the `Dt_L` and `Dt_T` fields using expl. above type defs
  const N = 1 # treated as 1 dimensional!
  proc toUnit[T: SomeUnit; U: SomeUnit](m: Measurement[T], _: typedesc[U]): Measurement[U] =
    # get conversion factor needed
    let factor = m.value.to(U).float / m.value.float
    result = (m * factor).to(U)

  proc toDiffConstant(D: Measurement[cm²•s⁻¹], v: Measurement[μm•ns⁻¹]): Measurement[float] =
    ## NOTE: Because `unchained` does not support sqrt units, we need to manually perform the
    ## unit conversion first before taking the sqrt. Done by calculating the required
    ## conversion factor by hand.
    # 1. First convert each unit to SI units by hand. This makes sure we get the correct
    # conversion factors (something is broken in Measuremancer right now :( )
    let dConstSq = 2 * N * D.toUnit(m²•s⁻¹) / v.toUnit(m•s⁻¹)
    # 2. convert to target unit for conversion factor and then back to float (without factors)
    result = sqrt(dConstSq.toUnit(μm²•cm⁻¹).to(float))

  result.Dt_L = toDiffConstant(result.D_L, result.zDrift)
  result.Dt_T = toDiffConstant(result.D_T, result.zDrift)

proc pretty*(sst: SteadyStateResult, indent: int): string =
  for field, val in fieldPairs(sst):
    result.add repeat(' ', indent) & field & " = " & $val & "\n"
  result = result.strip

proc `$`*(sst: SteadyStateResult): string = pretty(sst, indent = 0)

proc `$`*(mb: MagboltzResult): string =
  for field, val in fieldPairs(mb):
    when typeof(val) is SteadyStateResult:
      result.add field & " = " & pretty(val, 2)
    else:
      result.add field & " = " & $val & "\n"
  result = result.strip

proc cutFile(s: string): string =
  let idxStart = s.find(commonStart)
  let idxStop = if sstStart in s: s.find(sstStop) + sstSTop.len # includes sst
                else: s.find(commonStop) + commonStop.len
  result = s[idxStart ..< idxStop]

proc parseMagboltz(f: string): MagboltzResult =
  let dat = readFile(f)
  ## Magboltz output files depend on the kind of gas & settings we use.
  ## For strong fields with ionization and amplification an additional
  ## steady state solution is written. For pure drift fields this does
  ## not exist.
  ##
  ## First parse the common fields for all Magboltz output results
  ## ⇒ Nope, I used `npeg`.
  ## We still cut the input file from common data to either end of common
  ## or end of SST section.
  let datRead = dat.cutFile()
  result = parseMagboltzNpeg(datRead)

proc read(files: seq[string], alpha = false) =
  ## Note that the files are not very well specified and even differ slightly
  ## between the different Magboltz parameters. Therefore we do a hacky simple
  ## parse so that we might easily adjust it in the future. Not worth writing
  ## some kind of more efficient parser.
  ## The performance of parsing these files is completely irrelevant anyway.
  for f in files:
    let mb = parseMagboltz(f)
    echo "File: ", f
    echo mb
    echo "=================================================="

when isMainModule:
  import cligen
  dispatchMulti([read], [run], [test])
