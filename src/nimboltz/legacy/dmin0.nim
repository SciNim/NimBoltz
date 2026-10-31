## FUNCTION DMIN0 of magboltz-11.21.f (line 17913).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons

proc dmin0*(ia: int, ib: int): float64 =
  var ione: int
  ione = 1
  if ia > ib:
    result = float64(ib)
  elif ia < ione:
    result = float64(ione)
  else:
    result = float64(ia)
  return
