## SUBROUTINE READIN of magboltz-11.21.f (line 2516).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons

proc readin*(last: var int) =
  # READ IN OUTPUT CONTROL AND INTEGRATION DATA
  # CONVERTED TO FREE FORMAT READIN
  block:
    var rd = listRead()
    rd.get(inpt.ngas)
    rd.get(setp.nmax)
    rd.get(inpt.ipen)
    rd.get(thrm.ithrm)
    rd.get(inpt.efinal)
  # READ(5,2) NGAS,NMAX,IPEN,ITHRM,EFINAL
  # 2 FORMAT(4I10,F10.5)
  if not (inpt.ngas == 0):
    # GAS IDENTIFIERS
    block:
      var rd = listRead()
      rd.get(gasn.ngasn[1])
      rd.get(gasn.ngasn[2])
      rd.get(gasn.ngasn[3])
      rd.get(gasn.ngasn[4])
      rd.get(gasn.ngasn[5])
      rd.get(gasn.ngasn[6])
    # READ(5,3) NGASN(1),NGASN(2),NGASN(3),NGASN(4),NGASN(5),NGASN(6)
    # 3 FORMAT(6I5)
    # GAS PARAMETERS
    block:
      var rd = listRead()
      rd.get(ratio.frac[1])
      rd.get(ratio.frac[2])
      rd.get(ratio.frac[3])
      rd.get(ratio.frac[4])
      rd.get(ratio.frac[5])
      rd.get(ratio.frac[6])
      rd.get(inpt.tempc)
      rd.get(inpt.torr)
    # READ(5,4) FRAC(1),FRAC(2),FRAC(3),FRAC(4),FRAC(5),FRAC(6),TEMPC,
    # /TORR
    # 4 FORMAT(8F10.4)
    # FIELD VALUES
    block:
      var rd = listRead()
      rd.get(setp.efield)
      rd.get(bfld.bmag)
      rd.get(bfld.btheta)
    # READ(5,5) EFIELD,BMAG,BTHETA
    # 5 FORMAT(3F10.3)
    # CHECK INPUT
    # IF(NGAS.EQ.1.AND.IPEN.EQ.1) THEN
    # WRITE(6,992)
    # 992 FORMAT(' PROGRAM STOPPED IT IS NOT POSSIBLE TO HAVE PENNING TRANSF
    # /ERS IN PURE GASES.',/,' ONLY MIXTURES ALLOWED WITH IPEN=1')
    # STOP
    last = 0
    return
  last = 1
  return
