## SUBROUTINE TPLANEAT of magboltz-11.21.f (line 9913).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons

proc tplaneat*(t: float64, e1: float64, cx1: float64, cy1: float64, dcz1: float64, ap: float64, bp: float64, efld: float64, iplane: int) =
  var timesp, timlft, t2lft, a, b, eplane, wbt, coswt: float64
  var sinwt, const6, dcz2, xplane, yplane, zplane, vzplane: float64
  # STORES POSITION, TIME AND ENERGY AND SUMS REQUIRED
  # TO CALCULATE DEVIATIONS AND MEANS AT PLANE =IPLANE
  timesp = float64(iplane) * cion.tstep
  # CALC TIME LEFT TO ARRIVE AT PLANE
  timlft = timesp - outptt.st
  t2lft = timlft * timlft
  a = ap * timlft
  b = bp * t2lft
  eplane = e1 + a + b
  wbt = bfld.wb * timlft
  coswt = cos(wbt)
  sinwt = sin(wbt)
  const6 = sqrt(e1 / eplane)
  dcz2 = dcz1 * const6 + efld * timlft * cnsts1.const5 / sqrt(eplane)
  xplane = outptt.x + (cx1 * sinwt - cy1 * (1.0 - coswt)) / bfld.wb
  yplane = outptt.y + (cy1 * sinwt + cx1 * (1.0 - coswt)) / bfld.wb
  zplane = outptt.z + dcz1 * timlft * sqrt(e1) * cnsts1.const3 * 0.01 + t2lft * efld * cnsts1.const2
  vzplane = dcz2 * sqrt(eplane) * cnsts1.const3 * 0.01
  tplout.xtpl[iplane] = tplout.xtpl[iplane] + xplane
  tplout.ytpl[iplane] = tplout.ytpl[iplane] + yplane
  tplout.ztpl[iplane] = tplout.ztpl[iplane] + zplane
  tplout.xxtpl[iplane] = tplout.xxtpl[iplane] + xplane * xplane
  tplout.yytpl[iplane] = tplout.yytpl[iplane] + yplane * yplane
  tplout.zztpl[iplane] = tplout.zztpl[iplane] + zplane * zplane
  tplout.etpl[iplane] = tplout.etpl[iplane] + eplane
  tplout.ttpl[iplane] = tplout.ttpl[iplane] + outptt.st + timlft
  tplout.vztpl[iplane] = tplout.vztpl[iplane] + vzplane
  tplout.netpl[iplane] = tplout.netpl[iplane] + 1
  return
