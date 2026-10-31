## SUBROUTINE MONTEBT of magboltz-11.21.f (line 3755).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import dmin0, gerjan, sortt

const
  fmt201 = "(/,\'   VELZ    VELY  ENERGY   ID  DIFXX   DIFYY   DIFZZ   DIFYZ  DIFLNG  DIFTRN\',/)"
  fmt202 = "(2(F7.2,1X),F7.4,1X,I4,6(F7.1,1X))"

# NIMBOLTZ FIX MB-37: WAS 1000000; SETUP SETS NCOLM UP TO 2000000 (ALL
# OTHER MONTE* USE 2000000)
var xst: array[1..2000000, float64]
var yst: array[1..2000000, float64]
var zst: array[1..2000000, float64]
var sto: array[1..2000000, float64]
var temp: Arr2[float64, 6, 4000]

proc montebt*() =
  var wzst, wyst, avest, dfzzst, dfyyst, dfxxst, dfyzst, dflnst: array[1..10, float64]
  var dftrst: array[1..10, float64]
  var st1, st2, sumxx, sumyy, sumzz, sumyz, sumls, sumts: float64
  var sumvx, zold, yold, stold, st1old, st2old, szzold, sxxold: float64
  var syyold, syzold, svxold, slnold, strold, ebarold, ef100, rdum: float64
  var e1, tdash, f4, const9, const10, dcz1, dcx1, dcy1: float64
  var vtot, cx1, cy1, cz1, r1, t, wbt, coswt: float64
  var sinwt, dz, ev, cx2, cy2, cz2, r2, vgx: float64
  var vgy, vgz, eok, r4, tlim, const11, dxcom, dycom: float64
  var dzcom, t2, sdif, a2, b2, c2, dl2, dt2: float64
  var s1, ei, r9, extra, ran, s2, r3, r31: float64
  var f3, epsi, theta0, phi0, f8, f9, arg1, d: float64
  var q, f6, u, csqd, f5, argz, ebar, tcfsum: float64
  var twzst, twyst, tave, t2wzst, t2wyst, t2ave, tzzst, tyyst: float64
  var txxst, tyzst, tlnst, ttrst, t2zzst, t2yyst, t2xxst, t2yzst: float64
  var t2lnst, t2trst, ancatt, ancion: float64
  var n4000, n300, itmax, id, ncol, iextra, imbpt, j2m: int
  var j1, j2, kgas, ie, it, kdum, j, ncoldm: int
  var i, ipt, ik, ki, k: int
  # CALCULATES COLLISION EVENTS AND UPDATES DIFFUSION AND VELOCITY.
  # SUBROUTINE HANDLES MAGNETIC FIELD AND ELECTRIC FIELD
  # BFIELD ALONG X-AXIS EFIELD ALONG Z-AXIS (90 DEGREES).
  # BACKGROUND GAS MOTION INCLUDED AT TEMP = TEMPC
  vel.wx = 0.0
  velerr.dwx = 0.0
  diflab.difxz = 0.0
  diflab.difxy = 0.0
  diferb.dxzer = 0.0
  diferb.dxyer = 0.0
  outptt.x = 0.0
  outptt.y = 0.0
  outptt.z = 0.0
  outptt.st = 0.0
  st1 = 0.0
  st2 = 0.0
  sumxx = 0.0
  sumyy = 0.0
  sumzz = 0.0
  sumyz = 0.0
  sumls = 0.0
  sumts = 0.0
  sumvx = 0.0
  zold = 0.0
  yold = 0.0
  stold = 0.0
  st1old = 0.0
  st2old = 0.0
  szzold = 0.0
  sxxold = 0.0
  syyold = 0.0
  syzold = 0.0
  svxold = 0.0
  slnold = 0.0
  strold = 0.0
  ebarold = 0.0
  setp.small = 1e-20
  outptt.tmax1 = 0.0
  ef100 = setp.efield * 100.0
  rdum = setp.rstart
  e1 = setp.estart
  n4000 = 4000
  n300 = 300
  itmax = 10
  id = 0
  ncol = 0
  outptt.nnull = 0
  iextra = 0
  # GENERATE INITIAL RANDOM NUMBER FOR MAXWELL BOLTZMAN
  gerjan()
  # NIMBOLTZ FIX MB-05: WAS IMBTP=0 (TYPO), LEAVING IMBPT UNDEFINED
  imbpt = 0
  tdash = 0.0
  f4 = 2.0 * arccos((-1.0))
  const9 = cnsts1.const3 * 0.01
  const10 = const9 * const9
  # INITIAL DIRECTION COSINES
  dcz1 = cos(setp.theta)
  dcx1 = sin(setp.theta) * cos(setp.phi)
  dcy1 = sin(setp.theta) * sin(setp.phi)
  # INITIAL VELOCITY
  vtot = const9 * sqrt(e1)
  cx1 = dcx1 * vtot
  cy1 = dcy1 * vtot
  cz1 = dcz1 * vtot
  j2m = setp.nmax div itmax
  # MAIN LOOP
  fdo j1, 1, itmax:
    fdo j2, 1, j2m:
      block R1x:
        while true:
          block R1:
            r1 = drand48(rng)
            t = -(ln(r1) / thrm.tcfmx) + tdash
            tdash = t
            wbt = bfld.wb * t
            coswt = cos(wbt)
            sinwt = sin(wbt)
            dz = (cz1 * sinwt + (bfld.eovb - cy1) * (1.0 - coswt)) / bfld.wb
            ev = e1 + dz * ef100
            # CALCULATE ELECTRON VELOCITY IN LAB FRAME
            cx2 = cx1
            cy2 = (cy1 - bfld.eovb) * coswt + cz1 * sinwt + bfld.eovb
            cz2 = cz1 * coswt - (cy1 - bfld.eovb) * sinwt
            # FIND IDENTITY OF GAS FOR COLLISION
            kgas = 0
            block L4:
              if inpt.ngas == 1:
                r2 = drand48(rng)
                kgas = 1
                break L4
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
            # CALCULATE ENERGY WITH STATIONARY GAS TARGET,EOK
            eok = (pow2(cx2 - vgx) + pow2(cy2 - vgy) + pow2(cz2 - vgz)) / const10
            ie = int(trunc(eok / inpt.estep) + 1.0)
            ie = int(dmin0(ie, n4000))
            # TEST FOR REAL OR NULL COLLISION
            r4 = drand48(rng)
            tlim = gasTable[kgas].tcf[ie] / larget.tcfmax[kgas]
            if r4 > tlim:
              outptt.nnull = outptt.nnull + 1
              break R1
            break R1x
      ncol = ncol + 1
      # CALCULATE DIRECTION COSINES OF ELECTRON IN 0 KELVIN FRAME
      const11 = 1.0 / (const9 * sqrt(eok))
      dxcom = (cx2 - vgx) * const11
      dycom = (cy2 - vgy) * const11
      dzcom = (cz2 - vgz) * const11
      # CALCULATE POSITIONS AT INSTANT BEFORE COLLISION
      # ALSO UPDATE DIFFUSION  AND ENERGY CALCULATIONS.
      t2 = t * t
      if t >= outptt.tmax1: outptt.tmax1 = t
      tdash = 0.0
      # CALC NEW POSITION
      outptt.x = outptt.x + cx1 * t
      outptt.y = outptt.y + bfld.eovb * t + ((cy1 - bfld.eovb) * sinwt + cz1 * (1.0 - coswt)) / bfld.wb
      outptt.z = outptt.z + dz
      outptt.st = outptt.st + t
      it = int(trunc(t + 1.0))
      it = int(dmin0(it, n300))
      outptt.time[it] = outptt.time[it] + 1.0
      # ENERGY SPECTRUM FOR O KELVIN FRAME
      outptt.spec[ie] = outptt.spec[ie] + 1.0
      vel.wz = outptt.z / outptt.st
      vel.wy = outptt.y / outptt.st
      sumvx = sumvx + cx1 * cx1 * t2
      if not (id == 0):
        kdum = 0
        fdo j, 1, decor.ncorst:
          st2 = st2 + t
          ncoldm = ncol + kdum
          if ncoldm > decor.ncolm: ncoldm = ncoldm - decor.ncolm
          sdif = outptt.st - sto[ncoldm]
          sumxx = sumxx + pow2(outptt.x - xst[ncoldm]) * t / sdif
          if not (j1 < 3):
            st1 = st1 + t
            sumzz = sumzz + pow2(outptt.z - zst[ncoldm] - vel.wz * sdif) * t / sdif
            sumyy = sumyy + pow2(outptt.y - yst[ncoldm] - vel.wy * sdif) * t / sdif
            sumyz = sumyz + (outptt.z - zst[ncoldm] - vel.wz * sdif) * (outptt.y - yst[ncoldm] - vel.wy * sdif) * t / sdif
            a2 = pow2(vel.wz * sdif) + pow2(vel.wy * sdif)
            b2 = pow2(outptt.z - vel.wz * sdif - zst[ncoldm]) + pow2(outptt.y - vel.wy * sdif - yst[ncoldm])
            c2 = pow2(outptt.z - zst[ncoldm]) + pow2(outptt.y - yst[ncoldm])
            dl2 = pow2(a2 + b2 - c2) / (4.0 * a2)
            dt2 = b2 - dl2
            sumls = sumls + dl2 * t / sdif
            sumts = sumts + dt2 * t / sdif
          kdum = kdum + decor.ncorln
      xst[ncol] = outptt.x
      yst[ncol] = outptt.y
      zst[ncol] = outptt.z
      sto[ncol] = outptt.st
      if ncol >= decor.ncolm:
        id = id + 1
        outptt.xid = float64(id)
        ncol = 0
      # DETERMINATION OF REAL COLLISION TYPE
      r2 = drand48(rng)
      # FIND LOCATION WITHIN 4 UNITS IN COLLISION ARRAY
      sortt(kgas, i, r2, ie)
      while true:
        i = i + 1
        if not (gasTable[kgas].cf[ie, i] < r2): break
      s1 = gasTable[kgas].rgas[i]
      ei = gasTable[kgas].ein[i]
      # NEW CATCH FOR COLLISIONS WITH LARGE VIBRATIONAL ENERGY LOSS:
      # ENERGY LOSS = INCIDENT ENERGY - 0.024999
      if abs(ei - 0.024999000132083893) <= 0.0001:
        ei = ev - 0.024999000132083893
      if not (gasTable[kgas].ipn[i] <= 0):
        # USE FLAT DISTRIBUTION OF  ELECTRON ENERGY BETWEEN E-EION AND 0.0 EV
        # SAME AS IN BOLTZMANN
        r9 = drand48(rng)
        extra = r9 * (eok - ei)
        ei = extra + ei
        # IF FLOUORESCENCE OR AUGER ADD EXTRA ELECTRONS
        iextra = iextra + gasTable[kgas].augerCount[i]
      # GENERATE SCATTERING ANGLES AND UPDATE  LABORATORY COSINES AFTER
      # COLLISION ALSO UPDATE ENERGY OF ELECTRON.
      ipt = gasTable[kgas].iarry[i]
      outptt.icoll[kgas, ipt] = outptt.icoll[kgas, ipt] + 1
      outptt.icoln[kgas, i] = outptt.icoln[kgas, i] + 1
      if eok < ei:
        # FIX ENERGY LOSS SMALLER THAN INCIDENT ENERGY IF ERROR OCCURS
        ei = eok - 0.0001
      block L5:
        # IF EXCITATION THEN ADD PROBABILITY,PENFRAC(1,I), OF TRANSFER TO GIVE
        # IONISATION OF THE OTHER GASES IN THE MIXTURE
        if inpt.ipen == 0: break L5
        if gasTable[kgas].penfra[1, i] != 0.0:
          ran = drand48(rng)
          if ran > gasTable[kgas].penfra[1, i]: break L5
          # ADD EXTRA IONISATION COLLISION
          iextra = iextra + 1
      s2 = s1 * s1 / (s1 - 1.0)
      # ANISOTROPIC SCATTERING
      r3 = drand48(rng)
      if gasTable[kgas].index[i] == 1:
        r31 = drand48(rng)
        f3 = 1.0 - r3 * gasTable[kgas].angct[ie, i]
        if r31 > gasTable[kgas].psct[ie, i]: f3 = -f3
      elif gasTable[kgas].index[i] == 2:
        epsi = gasTable[kgas].psct[ie, i]
        f3 = 1.0 - 2.0 * r3 * (1.0 - epsi) / (1.0 + epsi * (1.0 - 2.0 * r3))
      else:
        # ISOTROPIC SCATTERING
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
      vtot = const9 * sqrt(e1)
      cx1 = dcx1 * vtot + vgx
      cy1 = dcy1 * vtot + vgy
      cz1 = dcz1 * vtot + vgz
      # CALCULATE ENERGY AND DIRECTION COSINES IN LAB FRAME
      e1 = (cx1 * cx1 + cy1 * cy1 + cz1 * cz1) / const10
      const11 = 1.0 / (const9 * sqrt(e1))
      dcx1 = cx1 * const11
      dcy1 = cy1 * const11
      dcz1 = cz1 * const11
    if j1 == 1: fwrite(fmt201)
    vel.wz = vel.wz * 1000000000.0
    vel.wy = vel.wy * 1000000000.0
    if st2 != 0.0: diflab.difxx = 5000000000000000.0 * sumxx / st2
    if st1 != 0.0: diflab.difyy = 5000000000000000.0 * sumyy / st1
    if st1 != 0.0: diflab.difzz = 5000000000000000.0 * sumzz / st1
    if st1 != 0.0: diflab.difyz = -(5000000000000000.0 * sumyz / st1)
    if st1 != 0.0: difvel.difln = 5000000000000000.0 * sumls / st1
    if st1 != 0.0: difvel.diftr = 5000000000000000.0 * sumts / st1
    if anist.niso == 0: diflab.difxx = 5000000000000000.0 * sumvx / outptt.st
    ebar = 0.0
    fdo ik, 1, 4000:
      tcfsum = 0.0
      fdo ki, 1, inpt.ngas:
        tcfsum = tcfsum + gasTable[ki].tcf[ik]
      ebar = ebar + mix2.eg[ik] * outptt.spec[ik] / tcfsum
    outptt.ave = ebar / outptt.st
    wzst[j1] = (outptt.z - zold) / (outptt.st - stold) * 1000000000.0
    wyst[j1] = (outptt.y - yold) / (outptt.st - stold) * 1000000000.0
    avest[j1] = (ebar - ebarold) / (outptt.st - stold)
    ebarold = ebar
    dfzzst[j1] = 0.0
    dfyyst[j1] = 0.0
    dfyzst[j1] = 0.0
    dflnst[j1] = 0.0
    dftrst[j1] = 0.0
    if j1 > 2:
      dfzzst[j1] = 5000000000000000.0 * (sumzz - szzold) / (st1 - st1old)
      dfyyst[j1] = 5000000000000000.0 * (sumyy - syyold) / (st1 - st1old)
      dfyzst[j1] = 5000000000000000.0 * (sumyz - syzold) / (st1 - st1old)
      dflnst[j1] = 5000000000000000.0 * (sumls - slnold) / (st1 - st1old)
      dftrst[j1] = 5000000000000000.0 * (sumts - strold) / (st1 - st1old)
    # NIMBOLTZ FIX MB-23: ST2=ST2OLD=0 IN THE FIRST SEGMENTS GAVE 0/0, I.E.
    # DXXER=NaN; GUARD AS IN MONTE/MONTEA
    dfxxst[j1] = 0.0
    if st2 != 0.0: dfxxst[j1] = 5000000000000000.0 * (sumxx - sxxold) / (st2 - st2old)
    if anist.niso == 0: dfxxst[j1] = 5000000000000000.0 * (sumvx - svxold) / (outptt.st - stold)
    zold = outptt.z
    yold = outptt.y
    stold = outptt.st
    st1old = st1
    st2old = st2
    svxold = sumvx
    szzold = sumzz
    sxxold = sumxx
    syyold = sumyy
    syzold = sumyz
    slnold = sumls
    strold = sumts
    fwrite(fmt202, vel.wz, vel.wy, outptt.ave, id, diflab.difxx, diflab.difyy, diflab.difzz, diflab.difyz, difvel.difln, difvel.diftr)
    # LOOP
  # CALCULATE ERRORS AND CHECK AVERAGES
  twzst = 0.0
  twyst = 0.0
  tave = 0.0
  t2wzst = 0.0
  t2wyst = 0.0
  t2ave = 0.0
  tzzst = 0.0
  tyyst = 0.0
  txxst = 0.0
  tyzst = 0.0
  tlnst = 0.0
  ttrst = 0.0
  t2zzst = 0.0
  t2yyst = 0.0
  t2xxst = 0.0
  t2yzst = 0.0
  t2lnst = 0.0
  t2trst = 0.0
  fdo k, 1, 10:
    twzst = twzst + wzst[k]
    twyst = twyst + wyst[k]
    tave = tave + avest[k]
    t2wzst = t2wzst + wzst[k] * wzst[k]
    t2wyst = t2wyst + wyst[k] * wyst[k]
    t2ave = t2ave + avest[k] * avest[k]
    txxst = txxst + dfxxst[k]
    t2xxst = t2xxst + dfxxst[k] * dfxxst[k]
    if not (k < 3):
      tzzst = tzzst + dfzzst[k]
      tyyst = tyyst + dfyyst[k]
      tyzst = tyzst + dfyzst[k]
      tlnst = tlnst + dflnst[k]
      ttrst = ttrst + dftrst[k]
      t2zzst = t2zzst + dfzzst[k] * dfzzst[k]
      t2yyst = t2yyst + dfyyst[k] * dfyyst[k]
      t2yzst = t2yzst + dfyzst[k] * dfyzst[k]
      t2lnst = t2lnst + dflnst[k] * dflnst[k]
      t2trst = t2trst + dftrst[k] * dftrst[k]
  velerr.dwz = 100.0 * sqrt((t2wzst - twzst * twzst / 10.0) / 9.0) / vel.wz
  velerr.dwy = 100.0 * sqrt((t2wyst - twyst * twyst / 10.0) / 9.0) / abs(vel.wy)
  outptt.den = 100.0 * sqrt((t2ave - tave * tave / 10.0) / 9.0) / outptt.ave
  diferb.dxxer = 100.0 * sqrt((t2xxst - txxst * txxst / 10.0) / 9.0) / diflab.difxx
  diferb.dyyer = 100.0 * sqrt((t2yyst - tyyst * tyyst / 8.0) / 7.0) / diflab.difyy
  diferb.dzzer = 100.0 * sqrt((t2zzst - tzzst * tzzst / 8.0) / 7.0) / diflab.difzz
  diferb.dyzer = 100.0 * sqrt((t2yzst - tyzst * tyzst / 8.0) / 7.0) / abs(diflab.difyz)
  diferl.dfler = 100.0 * sqrt((t2lnst - tlnst * tlnst / 8.0) / 7.0) / difvel.difln
  diferl.dfter = 100.0 * sqrt((t2trst - ttrst * ttrst / 8.0) / 7.0) / difvel.diftr
  velerr.dwz = velerr.dwz / sqrt(10.0)
  velerr.dwy = velerr.dwy / sqrt(10.0)
  outptt.den = outptt.den / sqrt(10.0)
  diferb.dxxer = diferb.dxxer / sqrt(10.0)
  diferb.dyyer = diferb.dyyer / sqrt(8.0)
  diferb.dzzer = diferb.dzzer / sqrt(8.0)
  diferb.dyzer = diferb.dyzer / sqrt(8.0)
  diferl.dfler = diferl.dfler / sqrt(8.0)
  diferl.dfter = diferl.dfter / sqrt(8.0)
  # CONVERT TO CM/SEC
  vel.wz = vel.wz * 100000.0
  vel.wy = vel.wy * 100000.0
  # CALCULATE TOWNSEND COEFICIENTS AND ERRORS
  ancatt = 0.0
  ancion = 0.0
  fdo i, 1, inpt.ngas:
    ancatt = ancatt + float64(outptt.icoll[i, 3])
    ancion = ancion + float64(outptt.icoll[i, 2])
  ancion = ancion + float64(iextra)
  ctwner.atter = 0.0
  if not (ancatt == 0.0):
    ctwner.atter = 100.0 * sqrt(ancatt) / ancatt
  ctowns.att = ancatt / (outptt.st * vel.wz) * 1000000000000.0
  ctwner.alper = 0.0
  if not (ancion == 0.0):
    ctwner.alper = 100.0 * sqrt(ancion) / ancion
  ctowns.alpha = ancion / (outptt.st * vel.wz) * 1000000000000.0
  return
