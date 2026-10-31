## SUBROUTINE OUTPUT1T of magboltz-11.21.f (line 7219).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import colft

const
  fmt15 = "(\'-----------------------------------------------------------------------------\')"
  fmt110 = "(2(/),\' NUMBER OF COLLISIONS IN FINAL ENERGY BIN =\',F8.1)"
  fmt220 = "(/,6X,\'TOTAL COLL. FREQ. =\',D11.4,\' (*10**12)/SEC.\',/,2X,\'INELASTIC COLL. FREQ. =\',D11.4,\' (*10**12)/SEC.\',/,4X,\'ELASTIC COLL. FREQ. =\',D11.4,\' (*10**12)/SEC.\',/,1X,\'IONISATION COLL. FREQ. =\',D11.4,\' (*10**12)/SEC.\',/,1X,\'ATTACHMENT COLL. FREQ. =\',D11.4,\' (*10**12)/SEC.\',/)"
  fmt1060 = "(/,2X,\'DETAILED COLLISION FREQUENCIES FOR EACH GAS IN UNITS OF 10**12/SEC. :\',2(/))"
  fmt1065 = "(/,3X,A25,/,\'----------------------------\',2(/))"
  fmt1070 = "(3X,A50,3X,D11.4,\' +-\',F8.4,\' %\')"

proc output1t*() =
  var freqel, freqsp, freine, freatt, freion: array[1..6, float64]
  var specs: array[1..40, float64]
  var freq1, freel1, freion1, freatt1, frein1, frelv, errfre: float64
  var nreal, i, j, k: int
  fwrite(fmt15)
  fwrite(fmt15)
  fwrite(fmt110, outptt.spec[4000])
  colft(freq1, freel1, freion1, freatt1, frein1, nreal)
  if ttrm.ttots == 0.0: ttrm.ttots = outptt.st
  fdo i, 1, inpt.ngas:
    freqel[i] = float64(outptt.icoll[i, 1]) / ttrm.ttots
    freqsp[i] = float64(outptt.icoll[i, 5]) / ttrm.ttots
    freine[i] = float64(outptt.icoll[i, 4]) / ttrm.ttots
    freatt[i] = float64(outptt.icoll[i, 3]) / ttrm.ttots
    freion[i] = float64(outptt.icoll[i, 2]) / ttrm.ttots
  fwrite(fmt220, freq1, frein1, freel1, freion1, freatt1)
  fwrite(fmt15)
  fwrite(fmt1060)
  fdo j, 1, inpt.ngas:
    fwrite(fmt1065, names.nameg[j])
    fdo k, 1, larget.last[j]:
      frelv = freq1 * float64(outptt.icoln[j, k]) / float64(nreal)
      if outptt.icoln[j, k] == 0:
        errfre = 0.0
      else:
        errfre = 100.0 * sqrt(float64(outptt.icoln[j, k])) / float64(outptt.icoln[j, k])
      fwrite(fmt1070, script.dscrpt[j, k], frelv, errfre)
  fwrite(fmt15)
  # WRITE(6,301)
  # 301 FORMAT(2(/),10X,' NORMALISED ENERGY DISTRIBUTION')
  # J1=0
  # J2=0
  # SPECN=DFLOAT(NREAL)
  # SMSPEC=0.0D0
  # DO 350 K=1,4000
  # SPEC(K)=SPEC(K)/SPECN
  # J1=J1+1
  # SMSPEC=SMSPEC+SPEC(K)
  # IF(J1.LT.100) GO TO 350
  # J2=J2+1
  # SPECS(J2)=SMSPEC
  # SMSPEC=0.0D0
  # J1=0
  # 350 CONTINUE
  # EPLT=EFINAL/40.0D0
  # DO 420 I=1,40
  # ENER=EPLT*(DFLOAT(I)-0.5D0)
  # WRITE(6,302) ENER,SPECS(I)
  # 302 FORMAT(6X,'E=',F11.3,6X,'SPEC=',D11.4)
  # 420 CONTINUE
  return
