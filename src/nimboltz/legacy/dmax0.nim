## FUNCTION DMAX0 of magboltz-11.21.f (line 17904).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons

proc dmax0*(ia: int, ib: int): float64 =
  if ia < ib:
    result = float64(ib)
  else:
    result = float64(ia)
  return
