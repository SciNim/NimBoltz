## SUBROUTINE ALPCLCA of magboltz-11.21.f (line 17210).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import friedland, montefs, pt, tof

const
  fmt25 = "(/,2X,\'SOLUTION FOR PULSED TOWNSEND AND TIME OF FLIGHT PARAMETERS\',/,\'  --------------------------------------------------------------\')"
  fmt26 = "(1(/),\'  TIME STEP BETWEEN SAMPLING PLANES =\',D12.5,\' PICOSECS.\',/)"
  fmt27 = "(/,\'PT IONISATION AND ATTACHMENT RATES *10**12/SEC\',/,\'ALPHA=\',D10.3,\' +- \',F6.2,\' %      ATT=\',D10.3,\' +- \',F6.2,\' %\')"
  fmt28 = "(/,\' TOF DIFFUSION\')"
  fmt29 = "(/,\'  DL=\',F8.1,\' +- \',F6.1,\' %       DT=\',F8.1,\' +- \',F6.1,\' %\')"
  fmt30 = "(/,\' TOF DRIFT VELOCITY\')"
  fmt31 = "(/,\'    WR=\',F8.2,\' +- \',F6.2,\' %\')"
  fmt32 = "(/,\'TOWNSEND COEFICIENT (ALPHA-ATT) CALCULATED FROM TOF RESULTS:\',/,\' ALPHA-ATT /CM.=\',D11.4,\' N.B. APPROXIMATE FORMULA NOT ACCURATE\',/)"

proc alpclca*() =
  var imax: int
  var corr, wrn, fc1, fc2, alpzz: float64
  # ESTIMATE TIME STEP FOR AVALANCHE SIMULATION IN TIME OF FLIGHT SIM.
  # USES ESTIMATED GAIN OF 3.0 BETWEEN PLANES.
  # CALLS TOF AND PT SUBROUTINES AND UPDATES ALPHA AND ATT
  # VERSION WITH BFIELD PARALLEL TO EFIELD
  imax = setp.nmax div 10000000
  if imax < 5: imax = 5
  setp.nmax = imax * 10000000
  corr = 760.0 * (inpt.tempc + 273.1499938964844) / (inpt.torr * 293.1499938964844)
  cion.alphast = 0.85 * abs(ctowns.alpha - ctowns.att)
  if ctowns.alpha > 10.0 * cion.alphast or ctowns.att > 10.0 * cion.alphast:
    # LARGE CANCELLATION BETWEEN ALPHA AND NETA REDUCE STEP BY 30.0
    if ctowns.alpha > 100.0 * corr:
      cion.alphast = cion.alphast * 30.0
    elif ctowns.alpha > 50.0 * corr:
      cion.alphast = cion.alphast * 15.0
    else:
      cion.alphast = cion.alphast * 8.0
  cion.vdst = vel.wz * 0.00001
  cion.tstep = ln(3.0) / (cion.alphast * cion.vdst * 100000.0)
  cion.tstep = cion.tstep * 1000000000000.0
  cion.tfinal = 7.0 * cion.tstep
  cion.itfinal = 7
  # CALC TIME OF FLIGHT AND PT
  fwrite(fmt25)
  fwrite(fmt26, cion.tstep)
  montefta()
  friedland()
  pt()
  tof()
  fwrite(fmt27, tofout.ralpha, tofout.ralper, tofout.rattof, tofout.ratofer)
  fwrite(fmt28)
  fwrite(fmt29, tofout.tofdl, tofout.tofdler, tofout.tofdt, tofout.tofdter)
  fwrite(fmt30)
  fwrite(fmt31, tofout.tofwr, tofout.tofwrer)
  # CALCULATE TOWNSEND SST COEFICIENTS FROM TOF RESULTS
  wrn = tofout.tofwr * 100000.0
  fc1 = wrn / (2.0 * tofout.tofdl)
  fc2 = (tofout.ralpha - tofout.rattof) * 1000000000000.0 / tofout.tofdl
  alpzz = fc1 - sqrt(fc1*fc1 - fc2)
  fwrite(fmt32, alpzz)
  # LOAD NEW ALPHA AND ATTACHMENT INTO COMMON BLOCKS
  ctowns.alpha = tofout.ralpha / tofout.tofwr * 10000000.0
  ctwner.alper = tofout.ralper * ctowns.alpha / 100.0
  ctowns.att = tofout.rattof / tofout.tofwr * 10000000.0
  ctwner.atter = tofout.ratofer * ctowns.att / 100.0
  # WRITE(6,33) ALPHA,ALPER,ATT,ATTER
  # 33  FORMAT(' TOF TOWNSEND IN UNITS OF 1/CM:',/,'   GAIN =',D12.4,' +-'
  # /,D12.4,/,' ATTACH =',D12.4,' +-',D12.4,/)
  return
