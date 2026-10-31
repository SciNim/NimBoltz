## SUBROUTINE COLF of magboltz-11.21.f (line 15616).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons

proc colf*(freq: var float64, freel: var float64, freion: var float64, freatt: var float64, frein: var float64, ntotal: var int) =
  var ninel, nela, natt, nion, nreal: int
  # CALCULATES REAL COLLISION FREQUENCIES FOR EVENT TYPES
  ninel = outpt.icoll[4] + outpt.icoll[5] + outpt.icoll[9] + outpt.icoll[10] + outpt.icoll[14] + outpt.icoll[15] + outpt.icoll[19] + outpt.icoll[20] + outpt.icoll[24] + outpt.icoll[25] + outpt.icoll[29] + outpt.icoll[30]
  nela = outpt.icoll[1] + outpt.icoll[6] + outpt.icoll[11] + outpt.icoll[16] + outpt.icoll[21] + outpt.icoll[26]
  natt = outpt.icoll[3] + outpt.icoll[8] + outpt.icoll[13] + outpt.icoll[18] + outpt.icoll[23] + outpt.icoll[28]
  nion = outpt.icoll[2] + outpt.icoll[7] + outpt.icoll[12] + outpt.icoll[17] + outpt.icoll[22] + outpt.icoll[27]
  ntotal = nela + ninel + natt + nion
  if ttrm.ttots == 0.0:
    nreal = ntotal
    ttrm.ttots = outpt.st
  else:
    nreal = ntotal
  freq = float64(nreal) / ttrm.ttots
  frein = float64(ninel) / ttrm.ttots
  freel = float64(nela) / ttrm.ttots
  freion = float64(nion) / ttrm.ttots
  freatt = float64(natt) / ttrm.ttots
  return
