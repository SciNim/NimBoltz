## SUBROUTINE FRIEDLAND of magboltz-11.21.f (line 14696).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons

const
  fmt900 = "(2(/),\' ESTIMATE USING FRIEDLAND :\',/,\' AVERAGE ENERGY =\',F8.3,\'EV.\',/,\' AVERAGE IONISATION =\',E11.4,\' *10**12/SEC\',/,\' AVERAGE ATTACHMENT =\',E11.4,\' *10**12/SEC\')"

proc friedland*() =
  var fr: array[1..4000, float64]
  var alfbar, attbar, ebar, fsum: float64
  var i: int
  # CALCULATE DISTRIBUTION FUNCTION USING FRIEDLAND TECHNIQUE
  # CF: J.FRIEDLAND  PHYSICS OF FLUIDS 20(1461)1977
  # USE DITRIBUTION FUNCTION TO CALCULATE AVERAGE ENERGY
  # IONISATION RATE AND ATTACHMENT RATE.
  alfbar = 0.0
  attbar = 0.0
  ebar = 0.0
  fsum = 0.0
  fdo i, 1, 4000:
    # IF(TCF(I).EQ.0.0.OR.TTOTS.EQ.0.0) THEN
    # WRITE(6,888) I,TCF(I),TTOTS
    # 888 FORMAT('IN FRIEDLAND I=',I4,' TCF(I)=',D12.4,' TTOTS=',D12.4)
    # ENDIF
    fr[i] = outpt.spec[i] / mixTable.tcf[i]
    ebar = ebar + mix2.eg[i] * outpt.spec[i] / mixTable.tcf[i]
    alfbar = alfbar + fred.fcion[i] * outpt.spec[i] / mixTable.tcf[i]
    attbar = attbar + fred.fcatt[i] * outpt.spec[i] / mixTable.tcf[i]
    fsum = fsum + fr[i]
  fdo i, 1, 4000:
    fr[i] = fr[i] / fsum
  ebar = ebar / ttrm.ttots
  alfbar = alfbar / ttrm.ttots
  attbar = attbar / ttrm.ttots
  fwrite(fmt900, ebar, alfbar, attbar)
  return
