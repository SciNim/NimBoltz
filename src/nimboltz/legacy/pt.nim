## SUBROUTINE PT of magboltz-11.21.f (line 6020).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons

const
  fmt900 = "(2(/),\' PULSED TOWNSEND RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\' PLANE NO.   (ION-ATT) FREQ.   ENERGY      WV      NO.OF ELECTRONS\',/)"
  fmt910 = "(2X,I2,4X,E15.4,7X,F7.2,4X,F7.2,3X,I8)"

proc pt*() =
  var antpl: array[1..8, float64]
  var i, ipl: int
  # CALCULATES PULSED TOWNSEND COEFFICIENTS
  antpl[1] = float64(tplout.netpl[1])
  pttof.ri[1] = (ln(antpl[1]) - ln(float64(cion.iprim))) / cion.tstep
  pttof.ept[1] = tplout.etpl[1] / antpl[1]
  pttof.ttest[1] = tplout.ttpl[1] / antpl[1]
  pttof.vzpt[1] = 1000000000.0 * tplout.vztpl[1] / antpl[1]
  block L11:
    fdo i, 2, cion.itfinal:
      if tplout.netpl[i] == 0:
        cion.itfinal = i - 1
        break L11
      antpl[i] = float64(tplout.netpl[i])
      pttof.ri[i] = (ln(antpl[i]) - ln(antpl[i - 1])) / cion.tstep
      pttof.ept[i] = tplout.etpl[i] / antpl[i]
      pttof.ttest[i] = tplout.ttpl[i] / antpl[i]
      pttof.vzpt[i] = 1000000000.0 * tplout.vztpl[i] / antpl[i]
  fwrite(fmt900, cion.itfinal)
  fdo ipl, 1, cion.itfinal:
    fwrite(fmt910, ipl, pttof.ri[ipl], pttof.ept[ipl], pttof.vzpt[ipl], tplout.netpl[ipl])
  return
