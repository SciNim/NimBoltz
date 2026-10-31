## SUBROUTINE ELIMITT of magboltz-11.21.f (line 4686).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import dmin0, gerjan, sortt

proc elimitt*(ielow: var int) =
  var isamp, n4000, imbpt, j2m, j1, kgas, ie, i: int
  var ipt: int
  var rdum, e1, tdash, const9, const10, dcz1, dcx1, dcy1: float64
  var bp, f1, f2, f4, r1, t, ap, ev: float64
  var const6, dcx2, dcy2, dcz2, r2, vgx, vgy, vgz: float64
  var vex, vey, vez, eok, r5, tlim, const11, dxcom: float64
  var dycom, dzcom, s1, ei, r9, extra, s2, r3: float64
  var r31, f3, epsi, theta0, r4, phi0, f8, f9: float64
  var arg1, d, q, f6, u, csqd, f5, argz: float64
  var const12, vxlab, vylab, vzlab: float64
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
  tdash = 0.0
  const9 = cnsts1.const3 * 0.01
  const10 = const9 * const9
  # GENERATE RANDOM NUMBER FOR MAXWELL BOLTZMAN
  gerjan()
  imbpt = 0
  # INITIAL DIRECTION COSINES
  dcz1 = cos(setp.theta)
  dcx1 = sin(setp.theta) * cos(setp.phi)
  dcy1 = sin(setp.theta) * sin(setp.phi)
  bp = setp.efield * setp.efield * cnsts1.const1
  f1 = setp.efield * cnsts1.const2
  f2 = setp.efield * cnsts1.const3
  f4 = 2.0 * arccos((-1.0))
  j2m = setp.nmax div isamp
  # MAIN LOOP
  fdo j1, 1, j2m:
    while true:
      r1 = drand48(rng)
      t = -(ln(r1) / thrm.tcfmx) + tdash
      tdash = t
      ap = dcz1 * f2 * sqrt(e1)
      ev = e1 + (ap + bp * t) * t
      # CALCULATE DIRECTION COSINES BEFORE COLLISION
      const6 = sqrt(e1 / ev)
      dcx2 = dcx1 * const6
      dcy2 = dcy1 * const6
      dcz2 = dcz1 * const6 + setp.efield * t * cnsts1.const5 / sqrt(ev)
      # FIND IDENTITY OF GAS FOR COLLISION
      kgas = 0
      r2 = drand48(rng)
      while true:
        kgas = kgas + 1
        if not (thrm.tcfmxg[kgas] < r2): break
      # CALCULATE GAS VELOCITY VECTORS VGX,VGY,VGZ
      imbpt = imbpt + 1
      if imbpt > 6:
        gerjan()
        imbpt = 1
      vgx = thrm.vtmb[kgas] * ranm.rnmx[imbpt]
      imbpt = imbpt + 1
      vgy = thrm.vtmb[kgas] * ranm.rnmx[imbpt]
      imbpt = imbpt + 1
      vgz = thrm.vtmb[kgas] * ranm.rnmx[imbpt]
      # CALCULATE ELECTRON VELOCITY VECTORS VEX,VEY,VEZ
      vex = dcx2 * const9 * sqrt(ev)
      vey = dcy2 * const9 * sqrt(ev)
      vez = dcz2 * const9 * sqrt(ev)
      # CALCULATE ENERGY WITH STATIONARY GAS TARGET ,EOK
      eok = (pow2(vex - vgx) + pow2(vey - vgy) + pow2(vez - vgz)) / const10
      ie = int(trunc(eok / inpt.estep) + 1.0)
      ie = int(dmin0(ie, n4000))
      # TEST FOR REAL OR NULL COLLISION
      r5 = drand48(rng)
      tlim = larget.tcf[kgas, ie] / larget.tcfmax[kgas]
      if not (r5 > tlim): break
    if ie == 4000:
      # ELECTRON ENERGY OUT OF RANGE
      ielow = 1
      return
    tdash = 0.0
    # CALCULATE DIRECTION COSINES OF ELECTRON IN 0 KELVIN FRAME
    const11 = 1.0 / (const9 * sqrt(eok))
    dxcom = (vex - vgx) * const11
    dycom = (vey - vgy) * const11
    dzcom = (vez - vgz) * const11
    # DETERMINATION OF REAL COLLISION TYPE
    r2 = drand48(rng)
    # FIND LOCATION WITHIN 4 UNITS IN COLLISION ARRAY
    sortt(kgas, i, r2, ie)
    while true:
      i = i + 1
      if not (larget.cf[kgas, ie, i] < r2): break
    s1 = larget.rgas[kgas, i]
    ei = larget.ein[kgas, i]
    if not (larget.ipn[kgas, i] <= 0):
      r9 = drand48(rng)
      extra = r9 * (eok - ei)
      ei = extra + ei
    # GENERATE SCATTERING ANGLES AND UPDATE  LABORATORY COSINES AFTER
    # COLLISION ALSO UPDATE ENERGY OF ELECTRON.
    ipt = larget.iarry[kgas, i]
    if eok < ei:
      ei = eok - 0.0001
    s2 = s1 * s1 / (s1 - 1.0)
    r3 = drand48(rng)
    if anist.index[kgas, i] == 1:
      r31 = drand48(rng)
      f3 = 1.0 - r3 * anist.angct[kgas, ie, i]
      if r31 > anist.psct[kgas, ie, i]: f3 = -f3
    elif anist.index[kgas, i] == 2:
      epsi = anist.psct[kgas, ie, i]
      f3 = 1.0 - 2.0 * r3 * (1.0 - epsi) / (1.0 + epsi * (1.0 - 2.0 * r3))
    else:
      f3 = 1.0 - 2.0 * r3
    theta0 = arccos(f3)
    r4 = drand48(rng)
    phi0 = f4 * r4
    f8 = sin(phi0)
    f9 = cos(phi0)
    arg1 = 1.0 - s1 * ei / eok
    arg1 = max(arg1, setp.small)
    d = 1.0 - f3 * sqrt(arg1)
    e1 = eok * (1.0 - ei / (s1 * eok) - 2.0 * d / s2)
    e1 = max(e1, setp.small)
    q = sqrt(eok / e1 * arg1) / s1
    q = min(q, 1.0)
    setp.theta = arcsin(q * sin(theta0))
    f6 = cos(setp.theta)
    u = (s1 - 1.0) * (s1 - 1.0) / arg1
    csqd = f3 * f3
    if f3 < 0.0 and csqd > u: f6 = -(1.0 * f6)
    f5 = sin(setp.theta)
    dzcom = min(dzcom, 1.0)
    argz = sqrt(dxcom * dxcom + dycom * dycom)
    block L130:
      if argz == 0.0:
        dcz1 = f6
        dcx1 = f9 * f5
        dcy1 = f8 * f5
        break L130
      dcz1 = dzcom * f6 + argz * f5 * f8
      dcy1 = dycom * f6 + f5 / argz * (dxcom * f9 - dycom * dzcom * f8)
      dcx1 = dxcom * f6 - f5 / argz * (dycom * f9 + dxcom * dzcom * f8)
    # TRANSFORM VELOCITY VECTORS TO LAB FRAME
    const12 = const9 * sqrt(e1)
    vxlab = dcx1 * const12 + vgx
    vylab = dcy1 * const12 + vgy
    vzlab = dcz1 * const12 + vgz
    # CALCULATE ENERGY AND DIRECTION COSINES IN LAB FRAME
    e1 = (vxlab * vxlab + vylab * vylab + vzlab * vzlab) / const10
    const11 = 1.0 / (const9 * sqrt(e1))
    dcx1 = vxlab * const11
    dcy1 = vylab * const11
    dcz1 = vzlab * const11
    # LOOP
  ielow = 0
  return
