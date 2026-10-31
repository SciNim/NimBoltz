## SUBROUTINE ELIMIT of magboltz-11.21.f (line 13531).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import dmin0, sort

proc elimit*(ielow: var int) =
  var isamp, n4000, intem, j2m, j1, i, ie, ipt: int
  var rdum, e1, tdash, dcz1, dcx1, dcy1, bp, f1: float64
  var f2, f4, deltae, r1, tlim, t, ap, ev: float64
  var r5, const6, dcx2, dcy2, dcz2, r2, s1, ei: float64
  var r9, extra, s2, r3, r31, f3, epsi, theta0: float64
  var r4, phi0, f8, f9, arg1, d, q, f6: float64
  var u, csqd, f5, argz: float64
  # CALCULATES COLLISION EVENTS AND TESTS TO FIND IF THE UPPER ENERGY
  # LIMIT FOR THE ELECTRON ENERGY IS EXCEEDED.
  # IF ENERGY LIMIT IS OK       IELOW = 0
  # IF ENERGY LIMIT IS EXCEEDED IELOW = 1
  # THE TEST IS CARRIED OUT FOR A SAMPLE OF COLLISIONS THAT ARE
  # SMALLER THAN THE FULL SAMPLE BY A FACTOR OF 1/ISAMP
  # USED WITH MAGNETIC FIELD B =0.0   ELECTRIC FIELD IN Z DIRECTION.
  isamp = 10
  setp.small = 1e-20
  rdum = setp.rstart
  e1 = setp.estart
  n4000 = 4000
  intem = 8
  tdash = 0.0
  # INITIAL DIRECTION COSINES
  dcz1 = cos(setp.theta)
  dcx1 = sin(setp.theta) * cos(setp.phi)
  dcy1 = sin(setp.theta) * sin(setp.phi)
  bp = setp.efield * setp.efield * cnsts1.const1
  f1 = setp.efield * cnsts1.const2
  f2 = setp.efield * cnsts1.const3
  f4 = 2.0 * arccos((-1.0))
  deltae = inpt.efinal / float64(intem)
  j2m = setp.nmax div isamp
  # MAIN LOOP
  fdo j1, 1, j2m:
    block R1x:
      while true:
        block R1:
          r1 = drand48(rng)
          i = int(trunc(e1 / deltae) + 1.0)
          i = int(dmin0(i, intem))
          tlim = large.tcfmax[i]
          t = -(ln(r1) / tlim) + tdash
          tdash = t
          ap = dcz1 * f2 * sqrt(e1)
          ev = e1 + (ap + bp * t) * t
          ie = int(trunc(ev / inpt.estep) + 1.0)
          ie = int(dmin0(ie, n4000))
          if mixTable.tcf[ie] > tlim:
            tdash = tdash + ln(r1) / tlim
            large.tcfmax[i] = 1.05 * large.tcfmax[i]
            break R1
          # TEST FOR REAL OR NULL COLLISION
          r5 = drand48(rng)
          tlim = mixTable.tcf[ie] / tlim
          if r5 > tlim: break R1
          break R1x
    # CALCULATE DIRECTION COSINES AT INSTANT BEFORE COLLISION
    if ie == 4000:
      # ELECTRON ENERGY OUT OF RANGE
      ielow = 1
      return
    tdash = 0.0
    const6 = sqrt(e1 / ev)
    dcx2 = dcx1 * const6
    dcy2 = dcy1 * const6
    dcz2 = dcz1 * const6 + setp.efield * t * cnsts1.const5 / sqrt(ev)
    # DETERMINATION OF REAL COLLISION TYPE
    r2 = drand48(rng)
    # FIND LOCATION WITHIN 4 UNITS IN COLLISION ARRAY
    sort(i, r2, ie)
    while true:
      i = i + 1
      if not (mixTable.cf[ie, i] < r2): break
    s1 = mixTable.rgas[i]
    ei = mixTable.ein[i]
    if not (mixTable.ipn[i] <= 0):
      r9 = drand48(rng)
      extra = r9 * (ev - ei)
      ei = extra + ei
    # GENERATE SCATTERING ANGLES AND UPDATE  LABORATORY COSINES AFTER
    # COLLISION ALSO UPDATE ENERGY OF ELECTRON.
    ipt = mixTable.iarry[i]
    if ev < ei:
      ei = ev - 0.0001
    s2 = s1 * s1 / (s1 - 1.0)
    r3 = drand48(rng)
    if mixTable.index[i] == 1:
      r31 = drand48(rng)
      f3 = 1.0 - r3 * mixTable.angct[ie, i]
      if r31 > mixTable.psct[ie, i]: f3 = -f3
    elif mixTable.index[i] == 2:
      epsi = mixTable.psct[ie, i]
      f3 = 1.0 - 2.0 * r3 * (1.0 - epsi) / (1.0 + epsi * (1.0 - 2.0 * r3))
    else:
      f3 = 1.0 - 2.0 * r3
    theta0 = arccos(f3)
    r4 = drand48(rng)
    phi0 = f4 * r4
    f8 = sin(phi0)
    f9 = cos(phi0)
    arg1 = 1.0 - s1 * ei / ev
    arg1 = max(arg1, setp.small)
    d = 1.0 - f3 * sqrt(arg1)
    e1 = ev * (1.0 - ei / (s1 * ev) - 2.0 * d / s2)
    e1 = max(e1, setp.small)
    q = sqrt(ev / e1 * arg1) / s1
    q = min(q, 1.0)
    setp.theta = arcsin(q * sin(theta0))
    f6 = cos(setp.theta)
    u = (s1 - 1.0) * (s1 - 1.0) / arg1
    csqd = f3 * f3
    if f3 < 0.0 and csqd > u: f6 = -(1.0 * f6)
    f5 = sin(setp.theta)
    dcz2 = min(dcz2, 1.0)
    argz = sqrt(dcx2 * dcx2 + dcy2 * dcy2)
    block L210:
      if argz == 0.0:
        dcz1 = f6
        dcx1 = f9 * f5
        dcy1 = f8 * f5
        break L210
      dcz1 = dcz2 * f6 + argz * f5 * f8
      dcy1 = dcy2 * f6 + f5 / argz * (dcx2 * f9 - dcy2 * dcz2 * f8)
      dcx1 = dcx2 * f6 - f5 / argz * (dcy2 * f9 + dcx2 * dcz2 * f8)
    # LOOP
  ielow = 0
  return
