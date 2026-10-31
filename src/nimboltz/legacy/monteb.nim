## SUBROUTINE MONTEB of magboltz-11.21.f (line 12641).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import dmin0, sort

const
  fmt201 = "(/,\'   VELZ    VELY  ENERGY   ID  DIFXX   DIFYY   DIFZZ   DIFYZ  DIFLNG  DIFTRN\',/)"
  fmt202 = "(2(F7.2,1X),F7.4,1X,I4,6(F7.1,1X))"

var xst: array[1..2000000, float64]
var yst: array[1..2000000, float64]
var zst: array[1..2000000, float64]
var sto: array[1..2000000, float64]

proc monteb*() =
  var wzst, wyst, avest, dfzzst, dfyyst, dfxxst, dfyzst, dflnst: array[1..10, float64]
  var dftrst: array[1..10, float64]
  var temp: array[1..4000, float64]
  var st1, st2, sumxx, sumyy, sumzz, sumyz, sumls, sumts: float64
  var sumvx, zold, yold, stold, st1old, st2old, szzold, sxxold: float64
  var syyold, syzold, svxold, slnold, strold, ebarold, ef100, rdum: float64
  var e1, tdash, const9, dcz1, dcx1, dcy1, vtot, cx1: float64
  var cy1, cz1, f4, deltae, r1, tlim, t, wbt: float64
  var coswt, sinwt, dz, ev, r5, t2, cx2, cy2: float64
  var cz2, dcx2, dcy2, dcz2, sdif, a2, b2, c2: float64
  var dl2, dt2, r2, s1, ei, r9, extra, ran: float64
  var s2, r3, r31, f3, epsi, theta0, r4, phi0: float64
  var f8, f9, arg1, d, q, f6, u, csqd: float64
  var f5, argz, ebar, twzst, twyst, tave, t2wzst, t2wyst: float64
  var t2ave, tzzst, tyyst, txxst, tyzst, tlnst, ttrst, t2zzst: float64
  var t2yyst, t2xxst, t2yzst, t2lnst, t2trst, ancatt, ancion: float64
  var n4000, n300, intem, itmax, id, ncol, iextra, j2m: int
  var j1, j2, i, ie, it, kdum, j, ncoldm: int
  var ipt, ik, k: int
  # CALCULATES COLLISION EVENTS AND UPDATES DIFFUSION AND VELOCITY.
  # SUBROUTINE HANDLES MAGNETIC FIELD AND ELECTRIC FIELD
  # BFIELD ALONG X-AXIS EFIELD ALONG Z-AXIS (90 DEGREES).
  vel.wx = 0.0
  velerr.dwx = 0.0
  diflab.difxz = 0.0
  diflab.difxy = 0.0
  diferb.dxzer = 0.0
  diferb.dxyer = 0.0
  outpt.x = 0.0
  outpt.y = 0.0
  outpt.z = 0.0
  outpt.st = 0.0
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
  outpt.tmax1 = 0.0
  ef100 = setp.efield * 100.0
  rdum = setp.rstart
  e1 = setp.estart
  n4000 = 4000
  n300 = 300
  intem = 8
  itmax = 10
  id = 0
  ncol = 0
  outpt.nnull = 0
  iextra = 0
  tdash = 0.0
  const9 = cnsts1.const3 * 0.01
  # INITIAL DIRECTION COSINES
  dcz1 = cos(setp.theta)
  dcx1 = sin(setp.theta) * cos(setp.phi)
  dcy1 = sin(setp.theta) * sin(setp.phi)
  # INITIAL VELOCITY
  vtot = const9 * sqrt(e1)
  cx1 = dcx1 * vtot
  cy1 = dcy1 * vtot
  cz1 = dcz1 * vtot
  f4 = 2.0 * arccos((-1.0))
  deltae = inpt.efinal / float64(intem)
  j2m = setp.nmax div itmax
  # MAIN LOOP
  fdo j1, 1, itmax:
    fdo j2, 1, j2m:
      block R1x:
        while true:
          block R1:
            r1 = drand48(rng)
            i = int(trunc(e1 / deltae) + 1.0)
            i = int(dmin0(i, intem))
            tlim = large.tcfmax[i]
            t = -(ln(r1) / tlim) + tdash
            tdash = t
            wbt = bfld.wb * t
            coswt = cos(wbt)
            sinwt = sin(wbt)
            dz = (cz1 * sinwt + (bfld.eovb - cy1) * (1.0 - coswt)) / bfld.wb
            ev = e1 + dz * ef100
            # IF(E.LT.0.0) WRITE(6,983) J2,DZ,E1,COSWT,SINWT,WBT,CY1
            # 983  FORMAT(2X,' J2=',I12,' DZ=',D12.3,' E1=',D12.3,' COSWT=',D12.3
            # /,' SINWT=',D12.3,' WBT=',D12.3,' CY1=',D12.3)
            ie = int(trunc(ev / inpt.estep) + 1.0)
            ie = int(dmin0(ie, n4000))
            if mixTable.tcf[ie] > tlim:
              tdash = tdash + ln(r1) / tlim
              large.tcfmax[i] = 1.05 * large.tcfmax[i]
              # WRITE(6,996)
              # 996   FORMAT(/,5X,' WARNING NULL COLLISION TIME INCREASED',/)
              break R1
            # TEST FOR REAL OR NULL COLLISION
            r5 = drand48(rng)
            tlim = mixTable.tcf[ie] / tlim
            if not (r5 <= tlim):
              outpt.nnull = outpt.nnull + 1
              break R1
            break R1x
      # CALCULATE DIRECTION COSINES AND POSITIONS AT INSTANT BEFORE COLLISION
      # ALSO UPDATE DIFFUSION  AND ENERGY CALCULATIONS.
      t2 = t * t
      if t >= outpt.tmax1: outpt.tmax1 = t
      tdash = 0.0
      # CALC VELOCITY
      cx2 = cx1
      cy2 = (cy1 - bfld.eovb) * coswt + cz1 * sinwt + bfld.eovb
      cz2 = cz1 * coswt - (cy1 - bfld.eovb) * sinwt
      # CALC DIRECTION COSINE
      vtot = sqrt(cx2 * cx2 + cy2 * cy2 + cz2 * cz2)
      dcx2 = cx2 / vtot
      dcy2 = cy2 / vtot
      dcz2 = cz2 / vtot
      ncol = ncol + 1
      # CALC NEW POSITION
      outpt.x = outpt.x + cx1 * t
      outpt.y = outpt.y + bfld.eovb * t + ((cy1 - bfld.eovb) * sinwt + cz1 * (1.0 - coswt)) / bfld.wb
      outpt.z = outpt.z + dz
      outpt.st = outpt.st + t
      it = int(trunc(t + 1.0))
      it = int(dmin0(it, n300))
      outpt.time[it] = outpt.time[it] + 1.0
      outpt.spec[ie] = outpt.spec[ie] + 1.0
      vel.wz = outpt.z / outpt.st
      vel.wy = outpt.y / outpt.st
      sumvx = sumvx + cx1 * cx1 * t2
      if not (id == 0):
        kdum = 0
        fdo j, 1, decor.ncorst:
          st2 = st2 + t
          ncoldm = ncol + kdum
          if ncoldm > decor.ncolm: ncoldm = ncoldm - decor.ncolm
          sdif = outpt.st - sto[ncoldm]
          sumxx = sumxx + pow2(outpt.x - xst[ncoldm]) * t / sdif
          if not (j1 < 3):
            st1 = st1 + t
            sumzz = sumzz + pow2(outpt.z - zst[ncoldm] - vel.wz * sdif) * t / sdif
            sumyy = sumyy + pow2(outpt.y - yst[ncoldm] - vel.wy * sdif) * t / sdif
            sumyz = sumyz + (outpt.z - zst[ncoldm] - vel.wz * sdif) * (outpt.y - yst[ncoldm] - vel.wy * sdif) * t / sdif
            a2 = pow2(vel.wz * sdif) + pow2(vel.wy * sdif)
            b2 = pow2(outpt.z - vel.wz * sdif - zst[ncoldm]) + pow2(outpt.y - vel.wy * sdif - yst[ncoldm])
            c2 = pow2(outpt.z - zst[ncoldm]) + pow2(outpt.y - yst[ncoldm])
            dl2 = pow2(a2 + b2 - c2) / (4.0 * a2)
            dt2 = b2 - dl2
            sumls = sumls + dl2 * t / sdif
            sumts = sumts + dt2 * t / sdif
          kdum = kdum + decor.ncorln
      xst[ncol] = outpt.x
      yst[ncol] = outpt.y
      zst[ncol] = outpt.z
      sto[ncol] = outpt.st
      if ncol >= decor.ncolm:
        id = id + 1
        outpt.xid = float64(id)
        ncol = 0
      # DETERMINATION OF REAL COLLISION TYPE
      r2 = drand48(rng)
      # FIND LOCATION WITHIN 4 UNITS IN COLLISION ARRAY
      sort(i, r2, ie)
      while true:
        i = i + 1
        if not (mixTable.cf[ie, i] < r2): break
      s1 = mixTable.rgas[i]
      ei = mixTable.ein[i]
      # NEW CATCH FOR COLLISIONS WITH LARGE VIBRATIONAL ENERGY LOSS:
      # ENERGY LOSS = INCIDENT ENERGY - 0.024999
      if abs(ei - 0.024999000132083893) <= 0.0001:
        ei = ev - 0.024999000132083893
      if not (mixTable.ipn[i] <= 0):
        # USE FLAT DISTRIBUTION OF  ELECTRON ENERGY BETWEEN E-EION AND 0.0 EV
        # SAME AS IN BOLTZMANN
        r9 = drand48(rng)
        extra = r9 * (ev - ei)
        ei = extra + ei
        # IF AUGER OR FLUORESCENCE ADD EXTRA IONISATION ELECTRONS
        iextra = iextra + ionfl1.nc0[i]
      # GENERATE SCATTERING ANGLES AND UPDATE  LABORATORY COSINES AFTER
      # COLLISION ALSO UPDATE ENERGY OF ELECTRON.
      ipt = mixTable.iarry[i]
      outpt.icoll[ipt] = outpt.icoll[ipt] + 1
      outpt.icoln[i] = outpt.icoln[i] + 1
      if ev < ei:
        # WRITE(6,994) E,EI,J2
        # 994  FORMAT(2X,' WARNING ENERGY =',F8.3,' LESS THAN ENERGY LOSS EI=',F8
        # /.3,' AT ITER=',I12,' DUE TO BINNING ERROR')
        # FIX ENERGY LOSS SMALLER THAN INCIDENT ENERGY IF ERROR OCCURS
        ei = ev - 0.0001
      block L5:
        # IF EXCITATION THEN ADD PROBABILITY,PENFRAC(1,I), OF TRANSFER TO GIVE
        # IONISATION OF THE OTHER GASES IN THE MIXTURE
        if inpt.ipen == 0: break L5
        if mixTable.penfra[1, i] != 0.0:
          ran = drand48(rng)
          if ran > mixTable.penfra[1, i]: break L5
          # ADD EXTRA IONISATION COLLISION
          iextra = iextra + 1
      s2 = s1 * s1 / (s1 - 1.0)
      # ANISOTROPIC SCATTERING
      r3 = drand48(rng)
      if mixTable.index[i] == 1:
        r31 = drand48(rng)
        f3 = 1.0 - r3 * mixTable.angct[ie, i]
        if r31 > mixTable.psct[ie, i]: f3 = -f3
      elif mixTable.index[i] == 2:
        epsi = mixTable.psct[ie, i]
        f3 = 1.0 - 2.0 * r3 * (1.0 - epsi) / (1.0 + epsi * (1.0 - 2.0 * r3))
      else:
        # ISOTROPIC SCATTERING
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
      vtot = const9 * sqrt(e1)
      argz = sqrt(dcx2 * dcx2 + dcy2 * dcy2)
      block L130:
        if argz == 0.0:
          # WRITE(6,9232) ITER,ID,E1
          # 9232  FORMAT(3X,'WARNING ARGZ= 0.0 AT ITER =',I10,' ID =',I10,' E1=',E1
          # /2.3)
          dcz1 = f6
          dcx1 = f9 * f5
          dcy1 = f8 * f5
          break L130
        dcz1 = dcz2 * f6 + argz * f5 * f8
        dcy1 = dcy2 * f6 + f5 / argz * (dcx2 * f9 - dcy2 * dcz2 * f8)
        dcx1 = dcx2 * f6 - f5 / argz * (dcy2 * f9 + dcx2 * dcz2 * f8)
      # CALCULATE VELOCITY VECTORS AFTER COLLISION
      cx1 = dcx1 * vtot
      cy1 = dcy1 * vtot
      cz1 = dcz1 * vtot
    if j1 == 1: fwrite(fmt201)
    vel.wz = vel.wz * 1000000000.0
    vel.wy = vel.wy * 1000000000.0
    if st2 != 0.0: diflab.difxx = 5000000000000000.0 * sumxx / st2
    if st1 != 0.0: diflab.difyy = 5000000000000000.0 * sumyy / st1
    if st1 != 0.0: diflab.difzz = 5000000000000000.0 * sumzz / st1
    if st1 != 0.0: diflab.difyz = -(5000000000000000.0 * sumyz / st1)
    if st1 != 0.0: difvel.difln = 5000000000000000.0 * sumls / st1
    if st1 != 0.0: difvel.diftr = 5000000000000000.0 * sumts / st1
    if anis.niso == 0: diflab.difxx = 5000000000000000.0 * sumvx / outpt.st
    ebar = 0.0
    fdo ik, 1, 4000:
      ebar = ebar + mix2.eg[ik] * outpt.spec[ik] / mixTable.tcf[ik]
    outpt.ave = ebar / outpt.st
    wzst[j1] = (outpt.z - zold) / (outpt.st - stold) * 1000000000.0
    wyst[j1] = (outpt.y - yold) / (outpt.st - stold) * 1000000000.0
    avest[j1] = (ebar - ebarold) / (outpt.st - stold)
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
    if anis.niso == 0: dfxxst[j1] = 5000000000000000.0 * (sumvx - svxold) / (outpt.st - stold)
    zold = outpt.z
    yold = outpt.y
    stold = outpt.st
    st1old = st1
    st2old = st2
    svxold = sumvx
    szzold = sumzz
    sxxold = sumxx
    syyold = sumyy
    syzold = sumyz
    slnold = sumls
    strold = sumts
    fwrite(fmt202, vel.wz, vel.wy, outpt.ave, id, diflab.difxx, diflab.difyy, diflab.difzz, diflab.difyz, difvel.difln, difvel.diftr)
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
  outpt.den = 100.0 * sqrt((t2ave - tave * tave / 10.0) / 9.0) / outpt.ave
  diferb.dxxer = 100.0 * sqrt((t2xxst - txxst * txxst / 10.0) / 9.0) / diflab.difxx
  diferb.dyyer = 100.0 * sqrt((t2yyst - tyyst * tyyst / 8.0) / 7.0) / diflab.difyy
  diferb.dzzer = 100.0 * sqrt((t2zzst - tzzst * tzzst / 8.0) / 7.0) / diflab.difzz
  diferb.dyzer = 100.0 * sqrt((t2yzst - tyzst * tyzst / 8.0) / 7.0) / abs(diflab.difyz)
  diferl.dfler = 100.0 * sqrt((t2lnst - tlnst * tlnst / 8.0) / 7.0) / difvel.difln
  diferl.dfter = 100.0 * sqrt((t2trst - ttrst * ttrst / 8.0) / 7.0) / difvel.diftr
  velerr.dwz = velerr.dwz / sqrt(10.0)
  velerr.dwy = velerr.dwy / sqrt(10.0)
  outpt.den = outpt.den / sqrt(10.0)
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
    ancatt = ancatt + float64(outpt.icoll[5 * i - 2])
    ancion = ancion + float64(outpt.icoll[5 * i - 3])
  ancion = ancion + float64(iextra)
  ctwner.atter = 0.0
  if not (ancatt == 0.0):
    ctwner.atter = 100.0 * sqrt(ancatt) / ancatt
  ctowns.att = ancatt / (outpt.st * vel.wz) * 1000000000000.0
  ctwner.alper = 0.0
  if not (ancion == 0.0):
    ctwner.alper = 100.0 * sqrt(ancion) / ancion
  ctowns.alpha = ancion / (outpt.st * vel.wz) * 1000000000000.0
  return
