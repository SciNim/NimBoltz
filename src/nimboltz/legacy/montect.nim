## SUBROUTINE MONTECT of magboltz-11.21.f (line 4194).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import dmin0, gerjan, sortt

const
  fmt201 = "(/,\'   VELZ    VELY    VELX   ENERGY  ID  DIFXX  DIFYY  DIFZZ  DIFYZ  DIFXZ  DIFXY\',/)"
  fmt202 = "(3(F7.2,1X),1X,F7.4,I4,6(F7.1))"

var xst: array[1..2000000, float64]
var yst: array[1..2000000, float64]
var zst: array[1..2000000, float64]
var sto: array[1..2000000, float64]
var temp: Arr2[float64, 6, 4000]

proc montect*() =
  var wzst, wyst, wxst, avest, dfzzst, dfyyst, dfxxst, dfyzst: array[1..10, float64]
  var dfxyst, dfxzst: array[1..10, float64]
  var difxxr, difyyr, difzzr, difyzr, difxzr, difxyr, st1, sumxx: float64
  var sumyy, sumzz, sumyz, sumxy, sumxz, zrold, yrold, xrold: float64
  var szzr, syyr, sxxr, sxyr, syzr, sxzr, stold, st1old: float64
  var st2old, szzold, syyold, sxxold, syzold, sxyold, sxzold, ebarold: float64
  var rcs, rsn, rtheta, efz100, efx100, f1, f4, const9: float64
  var const10, eovbr, rdum, e1, tdash, dcz1, dcx1, dcy1: float64
  var vtot, cx1, cy1, cz1, r1, t, wbt, coswt: float64
  var sinwt, dz, dx, ev, cx2, cy2, cz2, r2: float64
  var vgx, vgy, vgz, eok, r4, tlim, const11, dxcom: float64
  var dycom, dzcom, t2, sdif, s1, ei, r9, extra: float64
  var ran, s2, r3, r31, f3, epsi, theta0, phi0: float64
  var f8, f9, arg1, d, q, f6, u, csqd: float64
  var f5, argz, wzr, wyr, wxr, zr, yr, xr: float64
  var ebar, tcfsum, twzst, twyst, twxst, tave, t2wzst, t2wyst: float64
  var t2wxst, t2ave, tzzst, tyyst, txxst, txyst, tyzst, txzst: float64
  var t2zzst, t2yyst, t2xxst, t2xyst, t2yzst, t2xzst, ancatt, ancion: float64
  var n4000, n300, itmax, id, ncol, iextra, imbpt, j2m: int
  var j1, j2, kgas, ie, it, kdum, j, ncoldm: int
  var i, ipt, ik, ki, k: int
  # CALCULATES COLLISION EVENTS AND UPDATES DIFFUSION AND VELOCITY.
  # SUBROUTINE SOLVES MOTION IN COORDINATE SYSTEM WITH BFIELD
  # ALIGNED ALONG X AXIS AND ELECTRIC FIELD AT AN ANGLE BTHETA IN
  # THE X-Z PLANE.  THE VELOCITY VECTORS AND DIFFUSION ARE THEN
  # ROTATED INTO THE STANDARD COORDINATE FRAME WITH THE ELECTRIC-
  # FIELD ALONG Z-AXIS AND THE  BFIELD AT AN ANGLE BTHETA TO THE
  # ELECTRIC FIELD IN THE X-Z PLANE.
  # BACKGROUND GAS MOTION INCLUDED AT TEMP = TEMPC
  outptt.x = 0.0
  outptt.y = 0.0
  outptt.z = 0.0
  # CALC LONG AND TRANS ONLY FOR 90 DEGREES BETWEEN E AND B
  difvel.difln = 0.0
  difvel.diftr = 0.0
  diferl.dfler = 0.0
  diferl.dfter = 0.0
  difxxr = 0.0
  difyyr = 0.0
  difzzr = 0.0
  difyzr = 0.0
  difxzr = 0.0
  difxyr = 0.0
  outptt.st = 0.0
  st1 = 0.0
  sumxx = 0.0
  sumyy = 0.0
  sumzz = 0.0
  sumyz = 0.0
  sumxy = 0.0
  sumxz = 0.0
  zrold = 0.0
  yrold = 0.0
  xrold = 0.0
  szzr = 0.0
  syyr = 0.0
  sxxr = 0.0
  sxyr = 0.0
  syzr = 0.0
  sxzr = 0.0
  stold = 0.0
  st1old = 0.0
  st2old = 0.0
  szzold = 0.0
  syyold = 0.0
  sxxold = 0.0
  syzold = 0.0
  sxyold = 0.0
  sxzold = 0.0
  ebarold = 0.0
  setp.small = 1e-20
  outptt.tmax1 = 0.0
  setp.api = arccos((-1.0))
  # CALC ROTATION MATRIX ANGLES
  rcs = cos((bfld.btheta - 90.0) * setp.api / 180.0)
  rsn = sin((bfld.btheta - 90.0) * setp.api / 180.0)
  rtheta = bfld.btheta * setp.api / 180.0
  efz100 = setp.efield * 100.0 * sin(rtheta)
  efx100 = setp.efield * 100.0 * cos(rtheta)
  f1 = setp.efield * cnsts1.const2 * cos(rtheta)
  f4 = 2.0 * setp.api
  const9 = cnsts1.const3 * 0.01
  const10 = const9 * const9
  eovbr = bfld.eovb * sin(rtheta)
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
  imbpt = 0
  tdash = 0.0
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
            dz = (cz1 * sinwt + (eovbr - cy1) * (1.0 - coswt)) / bfld.wb
            dx = cx1 * t + f1 * t * t
            ev = e1 + dz * efz100 + dx * efx100
            # CALCULATE ELECTRON VELOCITY IN LAB FRAME
            cx2 = cx1 + 2.0 * f1 * t
            cy2 = (cy1 - eovbr) * coswt + cz1 * sinwt + eovbr
            cz2 = cz1 * coswt - (cy1 - eovbr) * sinwt
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
            tlim = larget.tcf[kgas, ie] / larget.tcfmax[kgas]
            if r4 > tlim:
              outptt.nnull = outptt.nnull + 1
              break R1
            break R1x
      ncol = ncol + 1
      # CALCULATE DIRECTION COSINES OF ELECTRON IN 0 KELVIN FRAME
      const11 = 1.0 / (const9 * sqrt(eok))
      # VTOT=1.0D0/CONST11
      dxcom = (cx2 - vgx) * const11
      dycom = (cy2 - vgy) * const11
      dzcom = (cz2 - vgz) * const11
      # CALCULATE POSITIONS AT INSTANT BEFORE COLLISION
      # ALSO UPDATE DIFFUSION  AND ENERGY CALCULATIONS.
      t2 = t * t
      if t >= outptt.tmax1: outptt.tmax1 = t
      tdash = 0.0
      # CALC NEW POSITION
      outptt.x = outptt.x + dx
      outptt.y = outptt.y + eovbr * t + ((cy1 - eovbr) * sinwt + cz1 * (1.0 - coswt)) / bfld.wb
      outptt.z = outptt.z + dz
      outptt.st = outptt.st + t
      it = int(trunc(t + 1.0))
      it = int(dmin0(it, n300))
      outptt.time[it] = outptt.time[it] + 1.0
      # ENERGY SPECTRUM FOR  0 KELVIN FRAME
      outptt.spec[ie] = outptt.spec[ie] + 1.0
      vel.wz = outptt.z / outptt.st
      vel.wy = outptt.y / outptt.st
      vel.wx = outptt.x / outptt.st
      if not (j1 < 3):
        kdum = 0
        fdo j, 1, decor.ncorst:
          ncoldm = ncol + kdum
          if ncoldm > decor.ncolm: ncoldm = ncoldm - decor.ncolm
          st1 = st1 + t
          sdif = outptt.st - sto[ncoldm]
          sumzz = sumzz + pow2(outptt.z - zst[ncoldm] - vel.wz * sdif) * t / sdif
          sumyy = sumyy + pow2(outptt.y - yst[ncoldm] - vel.wy * sdif) * t / sdif
          sumxx = sumxx + pow2(outptt.x - xst[ncoldm] - vel.wx * sdif) * t / sdif
          sumyz = sumyz + (outptt.z - zst[ncoldm] - vel.wz * sdif) * (outptt.y - yst[ncoldm] - vel.wy * sdif) * t / sdif
          sumxy = sumxy + (outptt.x - xst[ncoldm] - vel.wx * sdif) * (outptt.y - yst[ncoldm] - vel.wy * sdif) * t / sdif
          sumxz = sumxz + (outptt.x - xst[ncoldm] - vel.wx * sdif) * (outptt.z - zst[ncoldm] - vel.wz * sdif) * t / sdif
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
        if not (larget.cf[kgas, ie, i] < r2): break
      s1 = larget.rgas[kgas, i]
      ei = larget.ein[kgas, i]
      # NEW CATCH FOR COLLISIONS WITH LARGE VIBRATIONAL ENERGY LOSS:
      # ENERGY LOSS = INCIDENT ENERGY - 0.024999
      if abs(ei - 0.024999000132083893) <= 0.0001:
        ei = ev - 0.024999000132083893
      if not (larget.ipn[kgas, i] <= 0):
        # USE FLAT DISTRIBUTION OF  ELECTRON ENERGY BETWEEN E-EION AND 0.0 EV
        # SAME AS IN BOLTZMANN
        r9 = drand48(rng)
        extra = r9 * (eok - ei)
        ei = extra + ei
        # IF FLUORESCENCE OR AUGER ADD EXTRA ELECTRONS
        iextra = iextra + ionfl.nc0[kgas, i]
      # GENERATE SCATTERING ANGLES AND UPDATE  LABORATORY COSINES AFTER
      # COLLISION ALSO UPDATE ENERGY OF ELECTRON.
      ipt = larget.iarry[kgas, i]
      outptt.icoll[kgas, ipt] = outptt.icoll[kgas, ipt] + 1
      outptt.icoln[kgas, i] = outptt.icoln[kgas, i] + 1
      if eok < ei:
        ei = eok - 0.0001
      block L5:
        # IF EXCITATION THEN ADD PROBABILITY, PENFRA(1,I), OF TRANSFER TO GIVE
        # IONISATION OF THE OTHER GASES IN THE MIXTURE
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
    vel.wx = vel.wx * 1000000000.0
    # CALCULATE ROTATED VECTORS AND POSITIONS
    wzr = vel.wz * rcs - vel.wx * rsn
    wyr = vel.wy
    wxr = vel.wz * rsn + vel.wx * rcs
    zr = outptt.z * rcs - outptt.x * rsn
    yr = outptt.y
    xr = outptt.z * rsn + outptt.x * rcs
    ebar = 0.0
    fdo ik, 1, 4000:
      tcfsum = 0.0
      fdo ki, 1, inpt.ngas:
        tcfsum = tcfsum + larget.tcf[ki, ik]
      ebar = ebar + mix2.eg[ik] * outptt.spec[ik] / tcfsum
    outptt.ave = ebar / outptt.st
    wzst[j1] = (zr - zrold) / (outptt.st - stold) * 1000000000.0
    wyst[j1] = (yr - yrold) / (outptt.st - stold) * 1000000000.0
    wxst[j1] = (xr - xrold) / (outptt.st - stold) * 1000000000.0
    avest[j1] = (ebar - ebarold) / (outptt.st - stold)
    ebarold = ebar
    if not (j1 < 3):
      diflab.difxx = 5000000000000000.0 * sumxx / st1
      diflab.difyy = 5000000000000000.0 * sumyy / st1
      diflab.difzz = 5000000000000000.0 * sumzz / st1
      diflab.difyz = 5000000000000000.0 * sumyz / st1
      diflab.difxz = 5000000000000000.0 * sumxz / st1
      diflab.difxy = 5000000000000000.0 * sumxy / st1
      # CALCULATE  ROTATED TENSOR .
      difxxr = diflab.difxx * rcs * rcs + diflab.difzz * rsn * rsn + 2.0 * rcs * rsn * diflab.difxz
      difyyr = diflab.difyy
      difzzr = diflab.difxx * rsn * rsn + diflab.difzz * rcs * rcs - 2.0 * rcs * rsn * diflab.difxz
      difxyr = rcs * diflab.difxy + rsn * diflab.difyz
      difyzr = rsn * diflab.difxy - rcs * diflab.difyz
      difxzr = (rcs * rcs - rsn * rsn) * diflab.difxz - rsn * rcs * (diflab.difxx - diflab.difzz)
      sxxr = sumxx * rcs * rcs + sumzz * rsn * rsn + 2.0 * rcs * rsn * sumxz
      syyr = sumyy
      szzr = sumxx * rsn * rsn + sumzz * rcs * rcs - 2.0 * rcs * rsn * sumxz
      sxyr = rcs * sumxy + rsn * sumyz
      syzr = rsn * sumxy - rcs * sumyz
      sxzr = (rcs * rcs - rsn * rsn) * sumxz - rsn * rcs * (sumxx - sumzz)
    dfzzst[j1] = 0.0
    dfyyst[j1] = 0.0
    dfxxst[j1] = 0.0
    dfxyst[j1] = 0.0
    dfyzst[j1] = 0.0
    dfxzst[j1] = 0.0
    if j1 > 2:
      dfzzst[j1] = 5000000000000000.0 * (szzr - szzold) / (st1 - st1old)
      dfyyst[j1] = 5000000000000000.0 * (syyr - syyold) / (st1 - st1old)
      dfxxst[j1] = 5000000000000000.0 * (sxxr - sxxold) / (st1 - st1old)
      dfxyst[j1] = 5000000000000000.0 * (sxyr - sxyold) / (st1 - st1old)
      dfyzst[j1] = 5000000000000000.0 * (syzr - syzold) / (st1 - st1old)
      dfxzst[j1] = 5000000000000000.0 * (sxzr - sxzold) / (st1 - st1old)
    zrold = zr
    yrold = yr
    xrold = xr
    stold = outptt.st
    st1old = st1
    szzold = szzr
    syyold = syyr
    sxxold = sxxr
    sxyold = sxyr
    syzold = syzr
    sxzold = sxzr
    # OUTPUT ROTATED VECTORS AND TENSOR
    fwrite(fmt202, wzr, wyr, wxr, outptt.ave, id, difxxr, difyyr, difzzr, difyzr, difxzr, difxyr)
    # LOOP
  # CALCULATE ERRORS AND CHECK AVERAGES
  twzst = 0.0
  twyst = 0.0
  twxst = 0.0
  tave = 0.0
  t2wzst = 0.0
  t2wyst = 0.0
  t2wxst = 0.0
  t2ave = 0.0
  tzzst = 0.0
  tyyst = 0.0
  txxst = 0.0
  txyst = 0.0
  tyzst = 0.0
  txzst = 0.0
  t2zzst = 0.0
  t2yyst = 0.0
  t2xxst = 0.0
  t2xyst = 0.0
  t2yzst = 0.0
  t2xzst = 0.0
  fdo k, 1, 10:
    twzst = twzst + wzst[k]
    twyst = twyst + wyst[k]
    twxst = twxst + wxst[k]
    tave = tave + avest[k]
    t2wzst = t2wzst + wzst[k] * wzst[k]
    t2wyst = t2wyst + wyst[k] * wyst[k]
    t2wxst = t2wxst + wxst[k] * wxst[k]
    t2ave = t2ave + avest[k] * avest[k]
    if not (k < 3):
      tzzst = tzzst + dfzzst[k]
      tyyst = tyyst + dfyyst[k]
      txxst = txxst + dfxxst[k]
      txyst = txyst + dfxyst[k]
      tyzst = tyzst + dfyzst[k]
      txzst = txzst + dfxzst[k]
      t2zzst = t2zzst + dfzzst[k] * dfzzst[k]
      t2yyst = t2yyst + dfyyst[k] * dfyyst[k]
      t2xxst = t2xxst + dfxxst[k] * dfxxst[k]
      t2xyst = t2xyst + dfxyst[k] * dfxyst[k]
      t2yzst = t2yzst + dfyzst[k] * dfyzst[k]
      t2xzst = t2xzst + dfxzst[k] * dfxzst[k]
  velerr.dwz = 100.0 * sqrt((t2wzst - twzst * twzst / 10.0) / 9.0) / wzr
  velerr.dwy = 100.0 * sqrt((t2wyst - twyst * twyst / 10.0) / 9.0) / abs(wyr)
  velerr.dwx = 100.0 * sqrt((t2wxst - twxst * twxst / 10.0) / 9.0) / abs(wxr)
  outptt.den = 100.0 * sqrt((t2ave - tave * tave / 10.0) / 9.0) / outptt.ave
  diferb.dzzer = 100.0 * sqrt((t2zzst - tzzst * tzzst / 8.0) / 7.0) / difzzr
  diferb.dyyer = 100.0 * sqrt((t2yyst - tyyst * tyyst / 8.0) / 7.0) / difyyr
  diferb.dxxer = 100.0 * sqrt((t2xxst - txxst * txxst / 8.0) / 7.0) / difxxr
  diferb.dxyer = 100.0 * sqrt((t2xyst - txyst * txyst / 8.0) / 7.0) / abs(difxyr)
  diferb.dyzer = 100.0 * sqrt((t2yzst - tyzst * tyzst / 8.0) / 7.0) / abs(difyzr)
  diferb.dxzer = 100.0 * sqrt((t2xzst - txzst * txzst / 8.0) / 7.0) / abs(difxzr)
  velerr.dwz = velerr.dwz / sqrt(10.0)
  velerr.dwy = velerr.dwy / sqrt(10.0)
  velerr.dwx = velerr.dwx / sqrt(10.0)
  outptt.den = outptt.den / sqrt(10.0)
  diferb.dzzer = diferb.dzzer / sqrt(8.0)
  diferb.dyyer = diferb.dyyer / sqrt(8.0)
  diferb.dxxer = diferb.dxxer / sqrt(8.0)
  diferb.dxyer = diferb.dxyer / sqrt(8.0)
  diferb.dyzer = diferb.dyzer / sqrt(8.0)
  diferb.dxzer = diferb.dxzer / sqrt(8.0)
  # LOAD ROTATED VALUES INTO ARRAYS
  vel.wz = wzr
  vel.wx = wxr
  vel.wy = wyr
  diflab.difxx = difxxr
  diflab.difyy = difyyr
  diflab.difzz = difzzr
  diflab.difyz = difyzr
  diflab.difxz = difxzr
  diflab.difxy = difxyr
  # CONVERT TO CM/SEC.
  vel.wz = vel.wz * 100000.0
  vel.wy = vel.wy * 100000.0
  vel.wx = vel.wx * 100000.0
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
