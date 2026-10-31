## SUBROUTINE ALPCALC of magboltz-11.21.f (line 13997).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import friedland, montefd, monteft, output1, pt, sst, tof

const
  fmt11 = "(/,2X,\'SOLUTION FOR STEADY STATE TOWNSEND PARAMETERS\',/,\'  -------------------------------------------------\')"
  fmt12 = "(1(/),\' SPACE STEP BETWEEN SAMPLING PLANES =\',D12.5,\' MICRONS.\',/)"
  fmt18 = "(/,\' SST DRIFT VELOCITIES\')"
  fmt19 = "(/,\' VD=\',F9.2,\' +- \',F6.2,\' %   WS=\',F9.2,\' +- \',F6.2,\' %\')"
  fmt20 = "(/,\' SST DIFFUSION\')"
  fmt21 = "(/,\' DL=\',F9.2,\' +- \',F6.1,\' %   DT=\',F9.2,\' +- \',F6.2,\' %\')"
  fmt22 = "(/,\' SST TOWNSEND COEFICIENTS\')"
  fmt23 = "(/,\' ALPHA=\',F9.2\' +- \',F6.2,\' %    ATT=\',F9.2,\' +- \',F6.2, \' %\')"
  fmt25 = "(/,2X,\'SOLUTION FOR PULSED TOWNSEND AND TIME OF FLIGHT PARAMETERS\',/,\'  --------------------------------------------------------------\')"
  fmt26 = "(1(/),\'  TIME STEP BETWEEN SAMPLING PLANES =\',D12.5,\' PICOSECS.\',/)"
  fmt27 = "(/,\' PT IONISATION AND ATTACHMENT RATES *10**12/SEC\',/,\' ALPHA=\',D10.3,\' +- \',F6.2,\' %      ATT=\',D10.3,\' +- \',F6.2,\' %\')"
  fmt28 = "(/,\' TOF DIFFUSION\')"
  fmt29 = "(/,\'  DL=\',F8.1,\' +- \',F6.1,\' %       DT=\',F8.1,\' +- \',F6.1,\' %\')"
  fmt30 = "(/,\' TOF DRIFT VELOCITY\')"
  fmt31 = "(/,\'    WR=\',F8.2,\' +-\',F6.1,\' %\')"
  fmt888 = "(/,\'TOWNSEND COEFICIENT (ALPHA-ATT) CALCULATED FROM TOF RESULTS:\',/,\' ALPHA-ATT /CM.=\',D11.4,\' N.B. APPROXIMATE FORMULA NOT ACCURATE\',/)"

proc alpcalc*() =
  var imax: int
  var corr, zstepm, wrzn, fc1, fc2, alptest: float64
  # ESTIMATE TIME STEP AND SPACE STEPS FOR AVALANCHE SIMULATION IN
  # TIME OF FLIGHT AND STEADY STATE TOWNSEND SIMULATIONS.
  # USES ESTIMATED GAIN OF 3.0 BETWEEN PLANES.
  # CALLS SST AND TOF SUBROUTINES AND UPDATES ALPHA AND ATT
  imax = setp.nmax div 10000000
  if imax < 5: imax = 5
  setp.nmax = imax * 10000000
  corr = 760.0 * (inpt.tempc + 273.1499938964844) / (inpt.torr * 293.1499938964844)
  cion.alphast = 0.85 * abs(ctowns.alpha - ctowns.att)
  if ctowns.alpha > 10.0 * cion.alphast or ctowns.att > 10.0 * cion.alphast:
    # LARGE CANCELLATION BETWEEN ALPHA AND NETA REDUCE STEP SIZE BY 30.0
    if ctowns.alpha > 100.0 * corr:
      cion.alphast = cion.alphast * 30.0
    elif ctowns.alpha > 50.0 * corr:
      cion.alphast = cion.alphast * 15.0
    else:
      cion.alphast = cion.alphast * 8.0
  cion.vdst = vel.wz * 0.00001
  cion.tstep = ln(3.0) / (cion.alphast * cion.vdst * 100000.0)
  cion.zstep = ln(3.0) / cion.alphast
  # CONVERT TO METRES AND PICOSECONDS
  cion.tstep = cion.tstep * 1000000000000.0
  cion.zstep = cion.zstep * 0.01
  cion.tfinal = 7.0 * cion.tstep
  cion.itfinal = 7
  cion.zfinal = 8.0 * cion.zstep
  ctcalc.izfinal = 8
  ctcalc.zplane1 = cion.zstep
  ctcalc.zplane2 = 2.0 * cion.zstep
  ctcalc.zplane3 = 3.0 * cion.zstep
  ctcalc.zplane4 = 4.0 * cion.zstep
  ctcalc.zplane5 = 5.0 * cion.zstep
  ctcalc.zplane6 = 6.0 * cion.zstep
  ctcalc.zplane7 = 7.0 * cion.zstep
  ctcalc.zplane8 = 8.0 * cion.zstep
  # CALC SST
  fwrite(fmt11)
  zstepm = cion.zstep * 1000000.0
  fwrite(fmt12, zstepm)
  montefd()
  sst()
  # LOAD NEW ALPHA AND ATTACHMENT INTO COMMON BLOCKS
  ctowns.alpha = sstout.alphsst
  ctwner.alper = sstout.alpherr
  ctowns.att = sstout.attsst
  ctwner.atter = sstout.atterr
  fwrite(fmt18)
  fwrite(fmt19, sstout.vdout, sstout.vderr, sstout.wsout, sstout.wserr)
  fwrite(fmt20)
  fwrite(fmt21, sstout.dlout, sstout.dlerr, sstout.dtout, sstout.dterr)
  fwrite(fmt22)
  fwrite(fmt23, sstout.alphsst, sstout.alpherr, sstout.attsst, sstout.atterr)
  output1()
  # CALC TIME OF FLIGHT AND PT
  fwrite(fmt25)
  fwrite(fmt26, cion.tstep)
  monteft()
  friedland()
  pt()
  tof()
  fwrite(fmt27, tofout.ralpha, tofout.ralper, tofout.rattof, tofout.ratofer)
  fwrite(fmt28)
  fwrite(fmt29, tofout.tofdl, tofout.tofdler, tofout.tofdt, tofout.tofdter)
  fwrite(fmt30)
  fwrite(fmt31, tofout.tofwr, tofout.tofwrer)
  # CALCULATE TOWNSEND SST COEFICIENTS DIRECTLY FROM TOF RESULTS
  wrzn = tofout.tofwr * 100000.0
  fc1 = wrzn / (2.0 * tofout.tofdl)
  fc2 = (tofout.ralpha - tofout.rattof) * 1000000000000.0 / tofout.tofdl
  alptest = fc1 - sqrt(fc1*fc1 - fc2)
  fwrite(fmt888, alptest)
  return
