## SUBROUTINE SORT of magboltz-11.21.f (line 11775).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons

proc sort*(i: var int, r2: float64, ie: int) =
  var istep, incr, k: int
  # SELECTS COLLISION TYPE FROM COLLISION ARRAY BY BINARY STEP SAMPLING
  # REDUCES SAMPLING RANGE TO WITHIN 4 POSITIONS IN ARRAY
  # OUTPUT =  I   ( POSITION WITHIN 4 OF CORRECT VALUE)
  istep = large.isize
  incr = 0
  fdo k, 1, 12:
    i = incr
    if istep == 2: return
    i = incr + istep
    if not (i > large.last):
      if large.cf[ie, i] < r2:
        incr = incr + istep
    istep = istep div 2
  return
