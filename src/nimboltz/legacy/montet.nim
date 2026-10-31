## SUBROUTINE MONTET of magboltz-11.21.f (line 2859).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import dmin0, gerjan, sortt

const
  fmt50 = "(\' WARNING ENERGY OUT OF RANGE,INCREASE ELECTRON ENERGY INTEGRATION RANGE \')"
  fmt201 = "(/,\'    VEL      POS        TIME      ENERGY   COUNT   DIFXX     DIFYY     DIFZZ\',/)"
  fmt202 = "(1X,F8.2,2(1X,D10.3),F9.4,1X,I5,1X,3(2X,F8.1))"

var xst: array[1..2000000, float64]
var yst: array[1..2000000, float64]
var zst: array[1..2000000, float64]
var sto: array[1..2000000, float64]

proc montet*() =
  var wzst, avest, dfzzst, dfyyst, dfxxst: array[1..10, float64]
  var st1, st2, sume2, sumxx, sumyy, sumzz, sumvx, sumvy: float64
  var zold, stold, st1old, st2old, szzold, sxxold, syyold, svxold: float64
  var svyold, sme2old, rdum, e1, const9, const10, tdash, dcz1: float64
  var dcx1, dcy1, bp, f1, f2, f4, r1, t: float64
  var ap, ev, const6, dcx2, dcy2, dcz2, r2, vgx: float64
  var vgy, vgz, vex, vey, vez, eok, r4, tlim: float64
  var const11, dxcom, dycom, dzcom, t2, a, b, const7: float64
  var cx1, cy1, sdif, r3, s1, ei, r9, extra: float64
  var ran, s2, r31, f3, epsi, theta0, phi0, f8: float64
  var f9, arg1, d, q, f6, u, csqd, f5: float64
  var argz, const12, vxlab, vylab, vzlab, twzst, tave, t2wzst: float64
  var t2ave, tzzst, tyyst, txxst, t2zzst, t2yyst, t2xxst, ancatt: float64
  var ancion: float64
  var n4000, n300, itmax, id, ncol, iextra, imbpt, j2m: int
  var j1, j2, kgas, ie, it, kdum, jdum, ncoldm: int
  var i, ipt, k: int
  # CALCULATES COLLISION EVENTS AND UPDATES DIFFUSION AND VELOCITY.
  # USED WITH MAGNETIC FIELD B =0.0   ELECTRIC FIELD IN Z DIRECTION.
  # BACKGROUND GAS MOTION INCLUDED AT TEMP =  TEMPC
  vel.wx = 0.0
  vel.wy = 0.0
  velerr.dwx = 0.0
  velerr.dwy = 0.0
  diflab.difyz = 0.0
  diflab.difxy = 0.0
  diflab.difxz = 0.0
  diferb.dyzer = 0.0
  diferb.dxyer = 0.0
  diferb.dxzer = 0.0
  outptt.x = 0.0
  outptt.y = 0.0
  outptt.z = 0.0
  outptt.st = 0.0
  st1 = 0.0
  st2 = 0.0
  sume2 = 0.0
  sumxx = 0.0
  sumyy = 0.0
  sumzz = 0.0
  sumvx = 0.0
  sumvy = 0.0
  zold = 0.0
  stold = 0.0
  st1old = 0.0
  st2old = 0.0
  szzold = 0.0
  sxxold = 0.0
  syyold = 0.0
  svxold = 0.0
  svyold = 0.0
  sme2old = 0.0
  setp.small = 1e-20
  outptt.tmax1 = 0.0
  rdum = setp.rstart
  e1 = setp.estart
  const9 = cnsts1.const3 * 0.01
  const10 = const9 * const9
  n4000 = 4000
  n300 = 300
  itmax = 10
  id = 0
  ncol = 0
  outptt.nnull = 0
  iextra = 0
  # GENERATE INITIAL RANDOM MAXWELL BOLTZMAN NUMBERS
  gerjan()
  imbpt = 0
  tdash = 0.0
  # INITIAL DIRECTION COSINES
  dcz1 = cos(setp.theta)
  dcx1 = sin(setp.theta) * cos(setp.phi)
  dcy1 = sin(setp.theta) * sin(setp.phi)
  bp = setp.efield * setp.efield * cnsts1.const1
  f1 = setp.efield * cnsts1.const2
  f2 = setp.efield * cnsts1.const3
  f4 = 2.0 * arccos((-1.0))
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
            ap = dcz1 * f2 * sqrt(e1)
            ev = e1 + (ap + bp * t) * t
            # CALCULATE DIRECTION COSINES BEFORE COLLISION
            const6 = sqrt(e1 / ev)
            dcx2 = dcx1 * const6
            dcy2 = dcy1 * const6
            dcz2 = dcz1 * const6 + setp.efield * t * cnsts1.const5 / sqrt(ev)
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
            # CALCULATE ELECTRON VELOCITY VECTORS VEX,VEY,VEZ
            vex = dcx2 * const9 * sqrt(ev)
            vey = dcy2 * const9 * sqrt(ev)
            vez = dcz2 * const9 * sqrt(ev)
            # CALCULATE ENERGY WITH STATIONARY GAS TARGET, EOK
            eok = (pow2(vex - vgx) + pow2(vey - vgy) + pow2(vez - vgz)) / const10
            ie = int(trunc(eok / inpt.estep) + 1.0)
            ie = int(dmin0(ie, n4000))
            # TEST FOR REAL OR NULL COLLISION
            r4 = drand48(rng)
            tlim = larget.tcf[kgas, ie] / larget.tcfmax[kgas]
            if r4 > tlim:
              outptt.nnull = outptt.nnull + 1
              break R1
            break R1x
      ncol = ncol + 1
      # CALCULATE DIRECTION COSINES OF ELECTRON IN O KELVIN FRAME
      const11 = 1.0 / (const9 * sqrt(eok))
      dxcom = (vex - vgx) * const11
      dycom = (vey - vgy) * const11
      dzcom = (vez - vgz) * const11
      # CALCULATE  POSITIONS AT INSTANT BEFORE COLLISION
      # ALSO UPDATE DIFFUSION  AND ENERGY CALCULATIONS.
      t2 = t * t
      if t >= outptt.tmax1: outptt.tmax1 = t
      tdash = 0.0
      a = ap * t
      b = bp * t2
      sume2 = sume2 + t * (e1 + a / 2.0 + b / 3.0)
      const7 = const9 * sqrt(e1)
      a = t * const7
      cx1 = dcx1 * const7
      cy1 = dcy1 * const7
      outptt.x = outptt.x + dcx1 * a
      outptt.y = outptt.y + dcy1 * a
      outptt.z = outptt.z + dcz1 * a + t2 * f1
      outptt.st = outptt.st + t
      it = int(trunc(t + 1.0))
      it = int(dmin0(it, n300))
      outptt.time[it] = outptt.time[it] + 1.0
      # ENERGY SPECTRUM FOR O KELVIN FRAME
      outptt.spec[ie] = outptt.spec[ie] + 1.0
      vel.wz = outptt.z / outptt.st
      sumvx = sumvx + cx1 * cx1 * t2
      sumvy = sumvy + cy1 * cy1 * t2
      if not (id == 0):
        kdum = 0
        fdo jdum, 1, decor.ncorst:
          st2 = st2 + t
          ncoldm = ncol + kdum
          if ncoldm > decor.ncolm: ncoldm = ncoldm - decor.ncolm
          sdif = outptt.st - sto[ncoldm]
          sumxx = sumxx + pow2(outptt.x - xst[ncoldm]) * t / sdif
          sumyy = sumyy + pow2(outptt.y - yst[ncoldm]) * t / sdif
          if not (j1 < 3):
            st1 = st1 + t
            sumzz = sumzz + pow2(outptt.z - zst[ncoldm] - vel.wz * sdif) * t / sdif
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
      r3 = drand48(rng)
      # FIND LOCATION WITHIN 4 UNITS IN COLLISION ARRAY
      sortt(kgas, i, r3, ie)
      while true:
        i = i + 1
        if not (larget.cf[kgas, ie, i] < r3): break
      s1 = larget.rgas[kgas, i]
      ei = larget.ein[kgas, i]
      # NEW CATCH FOR COLLISIONS WITH LARGE VIBRATIONAL ENERGY LOSS
      # ENERGY LOSS = INCIDENT ENERGY -0.024999EV
      if abs(ei - 0.024999000132083893) <= 0.0001:
        ei = ev - 0.024999000132083893
      if not (larget.ipn[kgas, i] <= 0):
        # USE FLAT DISTRIBUTION OF  ELECTRON ENERGY BETWEEN E-EION AND 0.0 EV
        # SAME AS IN BOLTZMANN
        r9 = drand48(rng)
        extra = r9 * (eok - ei)
        ei = extra + ei
        # IF AUGER OR FLUORESCENCE ADD EXTRA IONISATION COLLISIONS
        iextra = iextra + ionfl.nc0[kgas, i]
      # GENERATE SCATTERING ANGLES AND UPDATE  LABORATORY COSINES AFTER
      # COLLISION ALSO UPDATE ENERGY OF ELECTRON.
      ipt = larget.iarry[kgas, i]
      outptt.icoll[kgas, ipt] = outptt.icoll[kgas, ipt] + 1
      outptt.icoln[kgas, i] = outptt.icoln[kgas, i] + 1
      if eok < ei:
        ei = eok - 0.0001
      block L5:
        # IF EXCITATION THEN ADD PROBABILITY ,PENFRA(1,I), OF TRANSFER TO
        # IONISATION OF THE OTHER GASES IN MIXTURE
        if inpt.ipen == 0: break L5
        if larget.penfra[kgas, 1, i] != 0.0:
          ran = drand48(rng)
          if ran > larget.penfra[kgas, 1, i]: break L5
          # ADD EXTRA IONISATION COLLISION
          iextra = iextra + 1
      s2 = s1 * s1 / (s1 - 1.0)
      # ANISOTROPIC SCATTERING
      r3 = drand48(rng)
      if anist.index[kgas, i] == 1:
        r31 = drand48(rng)
        f3 = 1.0 - r3 * anist.angct[kgas, ie, i]
        if r31 > anist.psct[kgas, ie, i]: f3 = -f3
      elif anist.index[kgas, i] == 2:
        epsi = anist.psct[kgas, ie, i]
        f3 = 1.0 - 2.0 * r3 * (1.0 - epsi) / (1.0 + epsi * (1.0 - 2.0 * r3))
      else:
        # ISOTROPIC  SCATTERING
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
    if j1 == 1: fwrite(fmt201)
    vel.wz = vel.wz * 1000000000.0
    outptt.ave = sume2 / outptt.st
    difvel.difln = 0.0
    if anist.niso == 0:
      diflab.difxx = 5000000000000000.0 * sumvx / outptt.st
      diflab.difyy = 5000000000000000.0 * sumvy / outptt.st
      dfxxst[j1] = 5000000000000000.0 * (sumvx - svxold) / (outptt.st - stold)
      dfyyst[j1] = 5000000000000000.0 * (sumvy - svyold) / (outptt.st - stold)
    else:
      if st2 != 0.0:
        diflab.difyy = 5000000000000000.0 * sumyy / st2
        diflab.difxx = 5000000000000000.0 * sumxx / st2
        dfxxst[j1] = 5000000000000000.0 * (sumxx - sxxold) / (st2 - st2old)
        dfyyst[j1] = 5000000000000000.0 * (sumyy - syyold) / (st2 - st2old)
      else:
        dfxxst[j1] = 0.0
        dfyyst[j1] = 0.0
    if st1 != 0.0:
      diflab.difzz = 5000000000000000.0 * sumzz / st1
      dfzzst[j1] = 5000000000000000.0 * (sumzz - szzold) / (st1 - st1old)
    else:
      dfzzst[j1] = 0.0
    wzst[j1] = (outptt.z - zold) / (outptt.st - stold) * 1000000000.0
    avest[j1] = (sume2 - sme2old) / (outptt.st - stold)
    zold = outptt.z
    stold = outptt.st
    st1old = st1
    st2old = st2
    svxold = sumvx
    svyold = sumvy
    szzold = sumzz
    sxxold = sumxx
    syyold = sumyy
    sme2old = sume2
    fwrite(fmt202, vel.wz, outptt.z, outptt.st, outptt.ave, id, diflab.difxx, diflab.difyy, diflab.difzz)
    # IF LAST ENERGY BIN IS HIGH INCREASE INTEGRATION ENERGY RANGE
    if outptt.spec[4000] > 1000.0 * float64(j1):
      fwrite(fmt50)
      fstop()
    # LOOP
  # CALCULATE ERRORS AND CHECK AVERAGES
  twzst = 0.0
  tave = 0.0
  t2wzst = 0.0
  t2ave = 0.0
  tzzst = 0.0
  tyyst = 0.0
  txxst = 0.0
  t2zzst = 0.0
  t2yyst = 0.0
  t2xxst = 0.0
  fdo k, 1, 10:
    twzst = twzst + wzst[k]
    tave = tave + avest[k]
    t2wzst = t2wzst + wzst[k] * wzst[k]
    t2ave = t2ave + avest[k] * avest[k]
    txxst = txxst + dfxxst[k]
    tyyst = tyyst + dfyyst[k]
    t2yyst = t2yyst + dfyyst[k] * dfyyst[k]
    t2xxst = t2xxst + dfxxst[k] * dfxxst[k]
    if not (k < 3):
      tzzst = tzzst + dfzzst[k]
      t2zzst = t2zzst + dfzzst[k] * dfzzst[k]
  velerr.dwz = 100.0 * sqrt((t2wzst - twzst * twzst / 10.0) / 9.0) / vel.wz
  outptt.den = 100.0 * sqrt((t2ave - tave * tave / 10.0) / 9.0) / outptt.ave
  diferb.dxxer = 100.0 * sqrt((t2xxst - txxst * txxst / 10.0) / 9.0) / diflab.difxx
  diferb.dyyer = 100.0 * sqrt((t2yyst - tyyst * tyyst / 10.0) / 9.0) / diflab.difyy
  diferb.dzzer = 100.0 * sqrt((t2zzst - tzzst * tzzst / 8.0) / 7.0) / diflab.difzz
  velerr.dwz = velerr.dwz / sqrt(10.0)
  outptt.den = outptt.den / sqrt(10.0)
  diferb.dxxer = diferb.dxxer / sqrt(10.0)
  diferb.dyyer = diferb.dyyer / sqrt(10.0)
  diferb.dzzer = diferb.dzzer / sqrt(8.0)
  difvel.difln = diflab.difzz
  difvel.diftr = (diflab.difxx + diflab.difyy) / 2.0
  # CONVERT CM/SEC
  vel.wz = vel.wz * 100000.0
  diferl.dfler = diferb.dzzer
  diferl.dfter = (diferb.dxxer + diferb.dyyer) / 2.0
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
