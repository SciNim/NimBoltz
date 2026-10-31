## SUBROUTINE OUTPUT2 of magboltz-11.21.f (line 15482).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import colf

const
  fmt15 = "(\'-----------------------------------------------------------------------------\')"
  fmt110 = "(2(/),\' NUMBER OF COLLISIONS IN FINAL ENERGY BIN =\',F8.1)"
  fmt220 = "(/,6X,\'TOTAL COLL. FREQ. =\',D11.4,\' (*10**12)/SEC.\',/,2X,\'  ELASTIC COLL. FREQ. =\',D11.4,\' (*10**12)/SEC.\',/,2X,\'INELASTIC COLL. FREQ. =\',D11.4,\' (*10**12)/SEC.\',/,1X,\'IONISATION COLL. FREQ. =\',D11.4,\' (*10**12)/SEC.\',/,1X,\'ATTACHMENT COLL. FREQ. =\',D11.4,\' (*10**12/SEC.\',/)"
  fmt301 = "(2(/),10X,\' NORMALISED ENERGY DISTRIBUTION\')"
  fmt302 = "(6X,\'E=\',F11.3,6X,\'SPEC=\',D11.4)"
  fmt1060 = "(/,2X,\'DETAILED COLLISION FREQUENCIES FOR EACH GAS IN UNITS OF 10**12/SEC. :\',2(/))"
  fmt1065 = "(3X,A25,/,\'-----------------------------\',2(/))"
  fmt1070 = "(3X,A50,3X,D11.4,\' +-\',F8.4,\' %\')"

proc output2*() =
  var freqel, freqsp, freine, freatt, freion: array[1..6, float64]
  var specs: array[1..40, float64]
  var freq, freel, freion1, freatt1, frein, frelv, errfre, specn: float64
  var smspec, eplt, ener: float64
  var nreal, i, j, k, j1, j2: int
  fwrite(fmt15)
  fwrite(fmt15)
  fwrite(fmt110, outpt.spec[4000])
  colf(freq, freel, freion1, freatt1, frein, nreal)
  ttrm.ttots = outpt.st
  fwrite(fmt220, freq, freel, frein, freion1, freatt1)
  fwrite(fmt15)
  # ILAST=DINT(TMAX1)+1
  # IF(ILAST.GT.120) ILAST=120
  # WRITE(6,1010) (TIME(I),I=1,ILAST)
  # 1010 FORMAT(/,6X,'DISTRIBUTION OF COLLISION TIMES IN 1 PECOSECOND BINS'
  # /,2(/),20(1X,6(F10.1,2X)/))
  # WRITE(6,15)
  fdo i, 1, inpt.ngas:
    freqel[i] = float64(outpt.icoll[5 * i - 4]) / ttrm.ttots
    freqsp[i] = float64(outpt.icoll[5 * i]) / ttrm.ttots
    freine[i] = float64(outpt.icoll[5 * i - 1]) / ttrm.ttots
    freatt[i] = float64(outpt.icoll[5 * i - 2]) / ttrm.ttots
    freion[i] = float64(outpt.icoll[5 * i - 3]) / ttrm.ttots
  # WRITE(6,1050) (NAMEG(I),FREQEL(I),FREQSP(I),FREINE(I),FREATT(I),
  # /FREION(I),I=1,NGAS)
  # 1050 FORMAT(/,5X,'COLLISION FREQUENCIES SORTED ACCORDING TO GAS AND TYP
  # /E OF COLLISION',/,5X,' IN UNITS OF 10**12/SEC.',2(/),'  GASES USED
  # /     ELASTIC     SUPERELAS   INELASTIC  ATTACHMENT  IONISATION ',2
  # /(/),6(A15,1X,5(D10.3,2X),/))
  # WRITE(6,15)
  fwrite(fmt1060)
  fdo j, 1, inpt.ngas:
    fwrite(fmt1065, names.nameg[j])
    fdo k, 1, mixTable.last:
      if mixTable.iarry[k] <= 5 * j and mixTable.iarry[k] > 5 * (j - 1):
        frelv = freq * float64(outpt.icoln[k]) / float64(nreal)
        if outpt.icoln[k] == 0:
          errfre = 0.0
        else:
          errfre = 100.0 * sqrt(float64(outpt.icoln[k])) / float64(outpt.icoln[k])
        fwrite(fmt1070, scrip.dscrpt[k], frelv, errfre)
  fwrite(fmt15)
  fwrite(fmt301)
  j1 = 0
  j2 = 0
  specn = float64(nreal)
  smspec = 0.0
  fdo k, 1, 4000:
    outpt.spec[k] = outpt.spec[k] / specn
    j1 = j1 + 1
    smspec = smspec + outpt.spec[k]
    if not (j1 < 100):
      j2 = j2 + 1
      specs[j2] = smspec
      smspec = 0.0
      j1 = 0
  eplt = inpt.efinal / 40.0
  fdo i, 1, 40:
    ener = eplt * (float64(i) - 0.5)
    fwrite(fmt302, ener, specs[i])
  return
