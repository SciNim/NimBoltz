## SUBROUTINE GERJAN of magboltz-11.21.f (line 2845).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons

proc gerjan*() =
  var j: int32
  var ran1, rdum, ran2: float64
  # GENERATES 6 RANDOM NUMBERS ACCORDING TO MAXWELLIAN DISTRIBUTION
  # CF GERJAN HAGELLAR
  fdo j, 1, 5, 2:
    ran1 = drand48(rng)
    ran2 = drand48(rng)
    ranm.rnmx[int(j)] = sqrt(-ln(ran1)) * cos(ran2 * ger.twopi)
    ranm.rnmx[int(j + 1)] = sqrt(-ln(ran1)) * sin(ran2 * ger.twopi)
  return
