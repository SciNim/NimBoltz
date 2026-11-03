## ALPCALC, ALPCALCT, ALPCLCA, ALPCLCAT, ALPCLCB, ALPCLCBT, ALPCLCC, ALPCLCCT
## of magboltz-11.21.f and the analyses they call (SST, PT/PTG/PTH,
## TOF/TOFG/TOFH, FRIEDLAND/FRIEDLANDT): bridge to nimboltz/townsend.
##
## Each driver chooses the sampling planes from the Townsend coefficients
## of MONTE, runs the avalanches (space planes for B = 0, then time
## planes), prints the per-plane tables and the results, and replaces the
## Townsend coefficients in /CTOWNS/ and their errors in /CTWNER/ (B = 0:
## from SST; otherwise from TOF). The eight drivers print with slightly
## different formats and differ in a few details; these are kept, marked
## with the routine they belong to.

import ../frt, commons
import colf, colft, output1
import ../[tables, vec3, fields, cascade, townsend, drift]
import montefs

const
  fmt11 = "(/,2X,\'SOLUTION FOR STEADY STATE TOWNSEND PARAMETERS\',/,\'  -------------------------------------------------\')"
  fmt12 = "(1(/),\' SPACE STEP BETWEEN SAMPLING PLANES =\',D12.5,\' MICRONS.\',/)"
  fmt18 = "(/,\' SST DRIFT VELOCITIES\')"
  fmt19 = "(/,\' VD=\',F9.2,\' +- \',F6.2,\' %   WS=\',F9.2,\' +- \',F6.2,\' %\')"
  fmt19t = "(/,\' VD=\',F9.1,\' +- \',F6.2,\' %   WS=\',F9.1,\' +- \',F6.2,\' %\')"
  fmt20 = "(/,\' SST DIFFUSION\')"
  fmt21 = "(/,\' DL=\',F9.2,\' +- \',F6.1,\' %   DT=\',F9.2,\' +- \',F6.2,\' %\')"
  fmt21t = "(/,\' DL=\',F9.1,\' +- \',F6.1,\' %   DT=\',F9.1,\' +- \',F6.2,\' %\')"
  fmt22 = "(/,\' SST TOWNSEND COEFICIENTS\')"
  fmt23 = "(/,\' ALPHA=\',F9.2\' +- \',F6.2,\' %    ATT=\',F9.2,\' +- \',F6.2, \' %\')"
  fmt23t = "(/,\' ALPHA =\',F9.1,\' +- \',F6.2,\' %    ATT=\',F9.1,\' +- \',F6.2,\' %\')"
  fmt25 = "(/,2X,\'SOLUTION FOR PULSED TOWNSEND AND TIME OF FLIGHT PARAMETERS\',/,\'  --------------------------------------------------------------\')"
  fmt26 = "(1(/),\'  TIME STEP BETWEEN SAMPLING PLANES =\',D12.5,\' PICOSECS.\',/)"
  fmt27e = "(/,\' PT IONISATION AND ATTACHMENT RATES *10**12/SEC\',/,\' ALPHA=\',D10.3,\' +- \',F6.2,\' %      ATT=\',D10.3,\' +- \',F6.2,\' %\')"
  fmt27 = "(/,\'PT IONISATION AND ATTACHMENT RATES *10**12/SEC\',/,\'ALPHA=\',D10.3,\' +- \',F6.2,\' %      ATT=\',D10.3,\' +- \',F6.2,\' %\')"
  fmt28 = "(/,\' TOF DIFFUSION\')"
  fmt28b = "(/,\'TOF DIFFUSION\')"
  fmt29 = "(/,\'  DL=\',F8.1,\' +- \',F6.1,\' %       DT=\',F8.1,\' +- \',F6.1,\' %\')"
  fmt29b = "(/,\'DZZ=\',F8.1,\' +- \',F5.1,\' %\',/,\'DXX=\',F8.1,\' +- \',F5.1,\' %\',/,\'DYY=\',F8.1,\' +- \',F5.1,\' %\',/,\'DYZ=\',F8.1,\' +- \',F5.1,\' %\')"
  fmt29c = "(/,\'DZZ=\',F8.1,\' +- \',F5.1,\' %\',/,\'DXX=\',F8.1,\' +- \',F5.1,\' %\',/,\'DYY=\',F8.1,\' +- \',F5.1,\' %\',/,\'DYZ=\',F8.1,\' +- \',F5.1,\' %\',/,\'DXZ=\',F8.1,\' +- \',F5.1,\' %\',/,\'DXY=\',F8.1,\' +- \',F5.1,\' %\')"
  fmt30 = "(/,\' TOF DRIFT VELOCITY\')"
  fmt30b = "(/,\'TOF DRIFT VELOCITY\')"
  fmt31e = "(/,\'    WR=\',F8.2,\' +-\',F6.1,\' %\')"
  fmt31a = "(/,\'    WR=\',F8.2,\' +- \',F6.2,\' %\')"
  fmt31b = "(/,\'WRZ=\',F8.2,\' +-\',F6.1,\' %    WRY=\',F8.2,\' +-\',F6.1,\' %\')"
  fmt31c = "(/,\'WRZ=\',F8.2,\' +-\',F6.1,\' %    WRY=\',F8.2,\' +-\',F6.1,\' %   WRX=\',F8.2,\' +-\',F6.1,\' %\')"
  fmt888 = "(/,\'TOWNSEND COEFICIENT (ALPHA-ATT) CALCULATED FROM TOF RESULTS:\',/,\' ALPHA-ATT /CM.=\',D11.4,\' N.B. APPROXIMATE FORMULA NOT ACCURATE\',/)"
  # SST
  fmt800s = "(2(/),\' STEADY STATE TOWNSEND RESULTS FOR\',I2,\' SEQUENTIAL SPACE PLANES\',2(/),\'PLANE   NEL     VD      WS     DL       DT     EBAR   ALFN    ALFNJ   ALFNE\',/)"
  fmt810s = "(1X,I2,2X,I7,2(1X,F6.1),2F9.1,F6.2,3F8.2)"
  # PT, PTG, PTH
  fmt900p = "(2(/),\' PULSED TOWNSEND RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\' PLANE NO.   (ION-ATT) FREQ.   ENERGY      WV      NO.OF ELECTRONS\',/)"
  fmt910p = "(2X,I2,4X,E15.4,7X,F7.2,4X,F7.2,3X,I8)"
  fmt900pg = "(2(/),\' PULSED TOWNSEND RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\'PLANE   (ION-ATT)FRQ.    ENERGY      WVZ       WVY   NO.OF ELECTRONS\',/)"
  fmt910pg = "(1X,I2,4X,D12.4,4X,F7.2,4X,F6.1,4X,F6.1,4X,I8)"
  fmt900ph = "(2(/),\' PULSED TOWNSEND RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\'PLANE   (ION-ATT)FRQ.    ENERGY      WVZ       WVY       WVX  NO.OF ELECTRNS\',/)"
  fmt910ph = "(1X,I2,4X,D12.4,4X,F7.2,4X,F6.1,4X,F6.1,4X,F6.1,4X,I8)"
  # TOF, TOFG, TOFH
  fmt900t = "(2(/),\' TIME OF FLIGHT RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\' PLANE NO.        DL         DX          DY            WR\',/)"
  fmt910t = "(3X,I3,4X,3F12.1,4X,F8.2)"
  fmt900tg = "(2(/),\' TIME OF FLIGHT RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\'PLANE    DZZ      DXX      DYY      DYZ        WRZ         WRY\',/)"
  fmt910tg = "(1X,I2,2X,4F9.1,4X,F8.2,4X,F8.2)"
  fmt900th = "(2(/),\' TIME OF FLIGHT RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\'PLANE    DZZ     DXX     DYY     DYZ     DXZ     DXY      WRZ     WRY     WRX\',/)"
  fmt910th = "(1X,I2,2X,6F8.1,3F8.2)"
  # FRIEDLAND
  fmt900f = "(2(/),\' ESTIMATE USING FRIEDLAND :\',/,\' AVERAGE ENERGY =\',F8.3,\'EV.\',/,\' AVERAGE IONISATION =\',E11.4,\' *10**12/SEC\',/,\' AVERAGE ATTACHMENT =\',E11.4,\' *10**12/SEC\')"

proc printSst(s: SteadyState) =
  fwrite(fmt800s, s.planes)
  for k in 1 .. s.planes:
    let p = s.plane[k]
    fwrite(fmt810s, k, p.n, p.vd, p.ws, p.dl, p.dt, p.e, p.alfn, p.alfnj, p.alfne)

proc printPt[G: static FieldGeometry](pt: PulsedTownsend, s: TimeSums) =
  when G in {fgElectric, fgParallel}:
    fwrite(fmt900p, pt.planes)
    for k in 1 .. pt.planes:
      let p = pt.plane[k]
      fwrite(fmt910p, k, p.rate, p.e, p.v.z, s.n[k])
  elif G == fgCrossed:
    fwrite(fmt900pg, pt.planes)
    for k in 1 .. pt.planes:
      let p = pt.plane[k]
      fwrite(fmt910pg, k, p.rate, p.e, p.v.z, p.v.y, s.n[k])
  else:
    fwrite(fmt900ph, pt.planes)
    for k in 1 .. pt.planes:
      let p = pt.plane[k]
      fwrite(fmt910ph, k, p.rate, p.e, p.v.z, p.v.y, p.v.x, s.n[k])

proc printTof[G: static FieldGeometry](tof: TimeOfFlight, planes: int) =
  when G in {fgElectric, fgParallel}:
    fwrite(fmt900t, planes)
    for k in 1 .. planes:
      let p = tof.plane[k]
      fwrite(fmt910t, k, p.d.zz, p.d.xx, p.d.yy, p.w.z)
  elif G == fgCrossed:
    fwrite(fmt900tg, planes)
    for k in 1 .. planes:
      let p = tof.plane[k]
      fwrite(fmt910tg, k, p.d.zz, p.d.xx, p.d.yy, p.d.yz, p.w.z, p.w.y)
  else:
    fwrite(fmt900th, planes)
    for k in 1 .. planes:
      let p = tof.plane[k]
      fwrite(fmt910th, k, p.d.zz, p.d.xx, p.d.yy, p.d.yz, p.d.xz, p.d.xy, p.w.z, p.w.y, p.w.x)

proc printFriedland(thermalMotion: bool) =
  ## FRIEDLAND, FRIEDLANDT: from the spectrum of the last run
  var r: Friedland
  if thermalMotion:
    var tcf: array[1..nEnergyBins, float64]
    for i in 1 .. nEnergyBins:
      var tcfsum = 0.0
      for k in 1 .. inpt.ngas: tcfsum = tcfsum + gasTable[k].tcf[i]
      tcf[i] = tcfsum
    r = friedland(outptt.spec, mix2.eg, tcf, fred.fcion, fred.fcatt, ttrm.ttots)
  else:
    r = friedland(outpt.spec, mix2.eg, mixTable.tcf, fred.fcion, fred.fcatt, ttrm.ttots)
  fwrite(fmt900f, r.e, r.ionisation, r.attachment)

proc run[G: static FieldGeometry](thermalMotion: bool) =
  var imax = setp.nmax div 10000000
  if imax < 5: imax = 5
  setp.nmax = imax * 10000000
  let (planes, _, _) = samplingPlanes(ctowns.alpha, ctowns.att, vel.wz, inpt.tempc, inpt.torr)
  when G == fgElectric:
    # steady-state Townsend
    fwrite(fmt11)
    fwrite(fmt12, planes.zstep * 1000000.0)
    let sr = avalanche[G, pkSpace](thermalMotion, planes)
    let sst = steadyState(sr.sums, sr.iprim, planes.zstep, planes.izfinal, sr.ratio)
    printSst(sst)
    ctowns.alpha = sst.alpha
    ctwner.alper = sst.alphaErr
    ctowns.att = sst.att
    ctwner.atter = sst.attErr
    fwrite(fmt18)
    fwrite(if thermalMotion: fmt19t else: fmt19, sst.vd, sst.vdErr, sst.ws, sst.wsErr)
    fwrite(fmt20)
    fwrite(if thermalMotion: fmt21t else: fmt21, sst.dl, sst.dlErr, sst.dt, sst.dtErr)
    fwrite(fmt22)
    fwrite(if thermalMotion: fmt23t else: fmt23, sst.alpha, sst.alphaErr, sst.att, sst.attErr)
    if not thermalMotion: output1()  # ALPCALCT has CALL OUTPUT1T commented out
  # pulsed Townsend and time of flight
  fwrite(fmt25)
  fwrite(fmt26, planes.tstep)
  let tr = avalanche[G, pkTime](thermalMotion, planes)
  printFriedland(thermalMotion)
  let pt = pulsedTownsend[G](tr.sums, tr.iprim, planes.tstep, planes.itfinal)
  printPt[G](pt, tr.sums)
  var ratio = tr.ratio
  when G == fgCrossed:
    # TOFG: the ratio from the collision frequencies of the run
    var dum1, dum2, frion, fratt, dum3: float64
    var idum4: int
    if thrm.ithrm == 1: colft(dum1, dum2, frion, fratt, dum3, idum4)
    if thrm.ithrm == 0: colf(dum1, dum2, frion, fratt, dum3, idum4)
    if frion > 0.0: ratio.ratio = fratt / frion  # MB-40
  let tof = timeOfFlight[G](tr.sums, pt, planes.tstep, ratio)
  printTof[G](tof, pt.planes)
  fwrite(when G == fgElectric: fmt27e else: fmt27, tof.alpha, tof.alphaErr, tof.att, tof.attErr)
  when G in {fgElectric, fgParallel}:
    fwrite(fmt28)
    fwrite(fmt29, tof.dl, tof.dlErr, tof.dt, tof.dtErr)
    fwrite(fmt30)
    fwrite(when G == fgElectric: fmt31e else: fmt31a, tof.w.z, tof.wErr.z)
    let wrzn = tof.w.z * 100000.0
    let fc1 = wrzn / (2.0 * tof.dl)
    let fc2 = (tof.alpha - tof.att) * 1000000000000.0 / tof.dl
    fwrite(fmt888, fc1 - sqrt(fc1*fc1 - fc2))
    when G == fgParallel:
      ctowns.alpha = tof.alpha / tof.w.z * 10000000.0
      ctwner.alper = tof.alphaErr * ctowns.alpha / 100.0
      ctowns.att = tof.att / tof.w.z * 10000000.0
      ctwner.atter = tof.attErr * ctowns.att / 100.0
  else:
    let d = tof.d
    let de = tof.dErr
    fwrite(fmt28b)
    when G == fgCrossed:
      fwrite(fmt29b, d.zz, de.zz, d.xx, de.xx, d.yy, de.yy, d.yz, de.yz)
    else:
      fwrite(fmt29c, d.zz, de.zz, d.xx, de.xx, d.yy, de.yy, d.yz, de.yz, d.xz, de.xz, d.xy, de.xy)
    fwrite(fmt30b)
    when G == fgCrossed:
      fwrite(fmt31b, tof.w.z, tof.wErr.z, tof.w.y, tof.wErr.y)
      let tofwr = sqrt(tof.w.z * tof.w.z + tof.w.y * tof.w.y)
    else:
      fwrite(fmt31c, tof.w.z, tof.wErr.z, tof.w.y, tof.wErr.y, tof.w.x, tof.wErr.x)
      let tofwr = sqrt(tof.w.z * tof.w.z + tof.w.y * tof.w.y + tof.w.x * tof.w.x)
    # as ALPCLCB/C: the drift velocity along E; ALPCLCBT/CT: the full one
    let wrzn = (if thermalMotion: tofwr else: tof.w.z) * 100000.0
    let fc1 = wrzn / (2.0 * d.zz)
    let fc2 = (tof.alpha - tof.att) * 1000000000000.0 / d.zz
    fwrite(fmt888, fc1 - sqrt(fc1*fc1 - fc2))
    ctowns.alpha = tof.alpha / tofwr * 10000000.0
    ctowns.att = tof.att / tofwr * 10000000.0
    when G == fgCrossed:
      if thermalMotion:
        ctwner.alper = tof.alphaErr * ctowns.alpha / 100.0
        ctwner.atter = tof.attErr * ctowns.att * 100.0  # as ALPCLCBT
      # ALPCLCB leaves /CTWNER/ alone
    else:
      ctwner.alper = tof.alphaErr * ctowns.alpha / 100.0
      ctwner.atter = tof.attErr * ctowns.att / 100.0

proc alpcalc*() = run[fgElectric](false)
proc alpcalct*() = run[fgElectric](true)
proc alpclca*() = run[fgParallel](false)
proc alpclcat*() = run[fgParallel](true)
proc alpclcb*() = run[fgCrossed](false)
proc alpclcbt*() = run[fgCrossed](true)
proc alpclcc*() = run[fgOblique](false)
proc alpclcct*() = run[fgOblique](true)
