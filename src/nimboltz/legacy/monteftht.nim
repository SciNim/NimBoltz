## SUBROUTINE MONTEFTHT of magboltz-11.21.f (line 8319).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import dmin0, gerjan, sortt, tplaneht

const
  fmt546 = "(2X,\' PROGRAM STOPPED . NPONT=\',I5,\' ITER=\',I10)"
  fmt836 = "(/,2X,\'AVERAGE ENERGY OF PRIMARY ELECTRON =\',F10.3,\' EV.\',/,\'   ENERGY SPREAD OF PRIMARY ELECTRON =\',F10.3,\' EV.\')"
  fmt878 = "(/,\' TOTAL NO OF ELECTRONS=\',I8,/,\' TOTAL NO OF NEG. IONS=\',I8,/,\' TOTAL NO OF PRIMARIES=\',I8)"
  fmt913 = "(3X,\' AFTER STORE ITER=\',I10,\' DZ=\',D12.3,\'E1=\',D12.3,\' COSWT=\',D12.3,\' SINWT=\',D12.3,\' WBT=\',D12.3,\' CY1=\',D12.3)"
  fmt940 = "(2(/),\' PROGRAM STOPPED  (TOO FEW COLLISIONS), DECREASE THE ESTIMATED ALPHA.  NCLUS = \',I7,\'  ITER =\',I9)"
  fmt944 = "(/,2X,\'PROGRAM STOPPED TOO MANY PRIMARIES IPRIM =\',I7)"
  fmt991 = "(2(/),\' PROGRAM STOPPED.  ITER =\',I10,\'    NMAX =\',I10,/,  \' NPONT=\',I4,\' NELEC=\',I8,\' IPRIM=\',I4,\' NMXADD=\',I3)"

var eprm: array[1..10000000, float64]
var temp: Arr2[float64, 6, 4000]

proc monteftht*() =
  var iespecp: array[1..100, int]
  var s, ytot, ytots, xtot, xtots, rdum, e1, zstrt: float64
  var ystrt, xstrt, tsstrt, rtheta, f4, const9, const10, dcz1: float64
  var dcx1, dcy1, vtot, cx1, cy1, cz1, e100, dcz100: float64
  var dcx100, dcy100, tdash, tstop, r1, t, wbt, coswt: float64
  var sinwt, dz, dx, ev, cx2, cy2, cz2, r2: float64
  var vgx, vgy, vgz, eok, r4, ttem, const11, dxcom: float64
  var dycom, dzcom, t2, s1, ei, eistr, r9, esec: float64
  var eavaug, r3, f3, theta0, f6, f5, phi0, f8: float64
  var f9, ran, asign, ran1, tpen, tstop1, s2, r31: float64
  var epsi, arg1, d, q, u, csqd, argz, f5s: float64
  var thsec, f6s, phis, f8s, f9s, wz, wy, wx: float64
  var wzr, wyr, wxr, aneion, anbt, eprmbar, e2prm, ebar: float64
  var eerr: float64
  var i, k, id, n4000, n300, n100, i100, ncol: int
  var nelec, neion, nmxadd, ntmpflg, npont, nclus, j1, imbpt: int
  var jprint, iprint, iter, iplane, idum, kgas, ie, it: int
  var ipt, ncltmp, naug, jfl, iplane1, kdum, jct: int
  # CALCULATES COLLISION EVENTS AND UPDATES DIFFUSION AND VELOCITY.
  # THIS ROUTINE HANDLES TERMINATIONS AT FIXED DRIFT TIMES.
  # SOLVES MOTION IN COORDINATE SYSTEM WITH BFIELD ALIGNED TO X-AXIS
  # ELECTRIC FIELD AT AN ANGLE BTHETA IN THE X-Z PLANE.
  # THE RESULTS FOR THE VELOCITY VECTORS  ARE THEN
  # ROTATED INTO THE STANDARD COORDINATE FRAME WITH THE ELECTRIC FIELD
  # ALONG THE Z-AXIS AND THE BFIELD AT AN ANGLE BTHETA TO THE ELECTRIC
  # FIELD IN THE X-Z PLANE
  s = 0.0
  outptt.st = 0.0
  outptt.x = 0.0
  outptt.y = 0.0
  outptt.z = 0.0
  ttrm.ztot = 0.0
  ttrm.ztots = 0.0
  ytot = 0.0
  ytots = 0.0
  xtot = 0.0
  xtots = 0.0
  ttrm.ttot = 0.0
  ttrm.ttots = 0.0
  setp.small = 1e-20
  outptt.tmax1 = 0.0
  rdum = setp.rstart
  e1 = setp.estart
  fdo i, 1, 300:
    outptt.time[i] = 0.0
  fdo i, 1, 5:
    fdo k, 1, 6:
      outptt.icoll[k, i] = 0
  fdo i, 1, 300:
    fdo k, 1, 6:
      outptt.icoln[k, i] = 0
  fdo i, 1, 4000:
    outptt.spec[i] = 0.0
  fdo i, 1, 100:
    iespecp[i] = 0
  fdo i, 1, 8:
    tplouth.etpl[i] = 0.0
    tplouth.xtpl[i] = 0.0
    tplouth.ytpl[i] = 0.0
    tplouth.ztpl[i] = 0.0
    tplouth.ttpl[i] = 0.0
    tplouth.xxtpl[i] = 0.0
    tplouth.yytpl[i] = 0.0
    tplouth.zztpl[i] = 0.0
    tplouth.yztpl[i] = 0.0
    tplouth.xztpl[i] = 0.0
    tplouth.xytpl[i] = 0.0
    tplouth.vztpl[i] = 0.0
    tplouth.vytpl[i] = 0.0
    tplouth.vxtpl[i] = 0.0
    tplouth.netpl[i] = 0
  id = 0
  n4000 = 4000
  n300 = 300
  n100 = 100
  i100 = 0
  ncol = 0
  outptt.nnull = 0
  nelec = 0
  neion = 0
  nmxadd = 0
  ntmpflg = 0
  npont = 0
  nclus = 0
  j1 = 1
  zstrt = 0.0
  ystrt = 0.0
  xstrt = 0.0
  tsstrt = 0.0
  setp.api = arccos((-1.0))
  # CALC ROTATION MATRIX ANGLES
  rots.rcs = cos((bfld.btheta - 90.0) * setp.api / 180.0)
  rots.rsn = sin((bfld.btheta - 90.0) * setp.api / 180.0)
  rtheta = bfld.btheta * setp.api / 180.0
  rots.efz100 = setp.efield * 100.0 * sin(rtheta)
  rots.efx100 = setp.efield * 100.0 * cos(rtheta)
  rots.f1 = setp.efield * cnsts1.const2 * cos(rtheta)
  f4 = 2.0 * setp.api
  const9 = cnsts1.const3 * 0.01
  const10 = const9 * const9
  rots.eovbr = bfld.eovb * sin(rtheta)
  # GENERATE INITIAL RANDOM NUMBERS FOR MAXWELL BOLTZMAN
  gerjan()
  imbpt = 0
  # INITIAL DIRECTION COSINES
  dcz1 = cos(setp.theta)
  dcx1 = sin(setp.theta) * cos(setp.phi)
  dcy1 = sin(setp.theta) * sin(setp.phi)
  # INITIAL VELOCITY
  vtot = const9 * sqrt(e1)
  cx1 = dcx1 * vtot
  cy1 = dcy1 * vtot
  cz1 = dcz1 * vtot
  e100 = e1
  dcz100 = dcz1
  dcx100 = dcx1
  dcy100 = dcy1
  jprint = setp.nmax div 10
  iprint = 0
  iter = 0
  iplane = 0
  cion.iprim = 0
  block R544x:
    while true:
      block R544:
        # LOOP FOR NEW STARTING ELECTRONS
        cion.iprim = cion.iprim + 1
        block L700:
          if cion.iprim > 1:
            # CHECK IF PROGRAM WILL EXCEED MAXIMUM NUMBER OF ITERATIONS
            # IN THIS CYCLE IF SO OUTPUT CURRENT RESULTS.
            if iter > setp.nmax:
              cion.iprim = cion.iprim - 1
              break L700
            outptt.x = 0.0
            outptt.y = 0.0
            outptt.z = 0.0
            dcz1 = dcz100
            dcx1 = dcx100
            dcy1 = dcy100
            e1 = e100
            vtot = const9 * sqrt(e1)
            cx1 = dcx1 * vtot
            cy1 = dcy1 * vtot
            cz1 = dcz1 * vtot
            nclus = nclus + 1
            outptt.st = 0.0
            tsstrt = 0.0
            zstrt = 0.0
            ystrt = 0.0
            xstrt = 0.0
            iplane = 0
          if cion.iprim > 10000000:
            fwrite(fmt944, cion.iprim)
            break L700
          eprm[cion.iprim] = e1
          idum = int(trunc(e1) + 1.0)
          idum = int(dmin0(idum, n100))
          iespecp[idum] = iespecp[idum] + 1
          block R555x:
            while true:
              block R555:
                # START OF LOOP FOR NEWLY CREATED ELECTRONS
                tdash = 0.0
                nelec = nelec + 1
                tstop = cion.tstep + float64(iplane) * cion.tstep
                block R1x:
                  while true:
                    block R1:
                      # MAIN LOOP
                      # IF(ITER.GT.NMAX)  GO TO 315
                      r1 = drand48(rng)
                      t = -(ln(r1) / thrm.tcfmx) + tdash
                      tdash = t
                      wbt = bfld.wb * t
                      coswt = cos(wbt)
                      sinwt = sin(wbt)
                      block L200:
                        block L20:
                          block R15x:
                            while true:
                              block R15:
                                if t + outptt.st >= tstop:
                                  iplane = iplane + 1
                                  tstop = tstop + cion.tstep
                                  # STORE POSITION AND ENERGY AT TIME PLANE =IPLANE.
                                  tplaneht(t, e1, cx1, cy1, cz1, setp.efield, iplane)
                                  # CHECK IF PASSED THROUGH MORE THAN ONE PLANE IN THIS STEP
                                  if t + outptt.st >= tstop and tstop <= cion.tfinal: break R15
                                  if t + outptt.st >= cion.tfinal:
                                    ttrm.ztot = ttrm.ztot + outptt.z
                                    ytot = ytot + outptt.y
                                    xtot = xtot + outptt.x
                                    ttrm.ttot = ttrm.ttot + outptt.st
                                    ttrm.ztots = ttrm.ztots + outptt.z - zstrt
                                    ytots = ytots + outptt.y - ystrt
                                    xtots = xtots + outptt.x - xstrt
                                    ttrm.ttots = ttrm.ttots + outptt.st - tsstrt
                                    tstop = cion.tstep
                                    # NO MORE ELECTRONS IN CASCADE TRY NEW PRIMARY ELECTRON
                                    if nelec == nclus + 1: break R544
                                    break L20
                                break R15x
                          dz = (cz1 * sinwt + (rots.eovbr - cy1) * (1.0 - coswt)) / bfld.wb
                          dx = cx1 * t + rots.f1 * t * t
                          ev = e1 + dz * rots.efz100 + dx * rots.efx100
                          if ev < 0.0:
                            fwrite(fmt913, iter, dz, ev, coswt, sinwt, wbt, cy1)
                            ev = 0.001
                          # CALCULATE ELECTRON VELOCITY IN LAB FRAME
                          cx2 = cx1 + 2.0 * rots.f1 * t
                          cy2 = (cy1 - rots.eovbr) * coswt + cz1 * sinwt + rots.eovbr
                          cz2 = cz1 * coswt - (cy1 - rots.eovbr) * sinwt
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
                          # CALCULATE ENERGY WITH STATIONARY GAS TARGET, EOK.
                          eok = (pow2(cx2 - vgx) + pow2(cy2 - vgy) + pow2(cz2 - vgz)) / const10
                          ie = int(trunc(eok / inpt.estep) + 1.0)
                          ie = int(dmin0(ie, n4000))
                          # TEST FOR REAL OR NULL COLLISION
                          r4 = drand48(rng)
                          ttem = gasTable[kgas].tcf[ie] / larget.tcfmax[kgas]
                          if not (r4 <= ttem):
                            outptt.nnull = outptt.nnull + 1
                            break R1
                          ncol = ncol + 1
                          # CALCULATE DIRECTION COSINES OF ELECTRON IN 0 KELVIN FRAME
                          const11 = 1.0 / (const9 * sqrt(eok))
                          dxcom = (cx2 - vgx) * const11
                          dycom = (cy2 - vgy) * const11
                          dzcom = (cz2 - vgz) * const11
                          # CALCULATE  POSITIONS AT INSTANT BEFORE COLLISION
                          # ALSO UPDATE DIFFUSION  AND ENERGY CALCULATIONS.
                          t2 = t * t
                          if t >= outptt.tmax1: outptt.tmax1 = t
                          tdash = 0.0
                          # CALC NEW POSITION
                          outptt.x = outptt.x + dx
                          outptt.y = outptt.y + rots.eovbr * t + ((cy1 - rots.eovbr) * sinwt + cz1 * (1.0 - coswt)) / bfld.wb
                          outptt.z = outptt.z + dz
                          outptt.st = outptt.st + t
                          it = int(trunc(t + 1.0))
                          it = int(dmin0(it, n300))
                          outptt.time[it] = outptt.time[it] + 1.0
                          # ENERGY SPECTRUM IN 0 KELVIN FRAME
                          outptt.spec[ie] = outptt.spec[ie] + 1.0
                          # DETERMINATION OF REAL COLLISION TYPE
                          r2 = drand48(rng)
                          # FIND LOCATION WITHIN 4 UNITS IN COLLISION ARRAY
                          sortt(kgas, i, r2, ie)
                          while true:
                            i = i + 1
                            if not (gasTable[kgas].cf[ie, i] < r2): break
                          s1 = gasTable[kgas].rgas[i]
                          ei = gasTable[kgas].ein[i]
                          if eok < ei:
                            # FIX ENERGY LOSS SMALLER THAN INCIDENT ENERGY IF A BINING ERROR OCCURS.
                            ei = eok - 0.0001
                          # CATCH FOR COLLISIONS WITH LARGE VIBRATIONAL ENERGY LOSS :
                          # ENERGY LOSS = INCIDENT ENERGY - 0.024999
                          if abs(ei - 0.024999000132083893) <= 0.0001:
                            ei = ev - 0.024999000132083893
                          if not (gasTable[kgas].ipn[i] == 0):
                            # ATTACHMENT
                            if gasTable[kgas].ipn[i] == (-1):
                              neion = neion + 1
                              ipt = gasTable[kgas].iarry[i]
                              id = id + 1
                              iter = iter + 1
                              iprint = iprint + 1
                              outptt.icoll[kgas, ipt] = outptt.icoll[kgas, ipt] + 1
                              outptt.icoln[kgas, i] = outptt.icoln[kgas, i] + 1
                              it = int(trunc(t + 1.0))
                              it = int(dmin0(it, n300))
                              outptt.time[it] = outptt.time[it] + 1.0
                              ttrm.ztot = ttrm.ztot + outptt.z
                              ytot = ytot + outptt.y
                              xtot = xtot + outptt.x
                              ttrm.ttot = ttrm.ttot + outptt.st
                              ttrm.ztots = ttrm.ztots + outptt.z - zstrt
                              ytots = ytots + outptt.y - ystrt
                              xtots = xtots + outptt.x - xstrt
                              ttrm.ttots = ttrm.ttots + outptt.st - tsstrt
                              # ELECTRON CAPTURED START NEW PRIMARY
                              if nelec == nclus + 1: break R544
                              # ELECTRON CAPTURED TAKE NEXT ELECTRON FROM STORE
                              break L20
                            eistr = ei
                            r9 = drand48(rng)
                            # ESEC=R9*(EOK-EI)
                            # USE OPAL PETERSON AND BEATY SPLITTING FACTOR.
                            esec = gasTable[kgas].wpl[i] * tan(r9 * arctan((eok - ei) / (2.0 * gasTable[kgas].wpl[i])))
                            esec = gasTable[kgas].wpl[i] * pow(esec / gasTable[kgas].wpl[i], 0.9524000287055969)
                            ei = esec + ei
                            # STORE POSITION ,ENERGY, DIRECTION COSINES AND TIME OF GENERATION
                            # OF IONISATION ELECTRON
                            nclus = nclus + 1
                            npont = npont + 1
                            nmxadd = max(npont, nmxadd)
                            if npont > 8000:
                              fwrite(fmt546, npont, iter)
                              fstop()
                            ispt.xs[npont] = outptt.x
                            ispt.ys[npont] = outptt.y
                            ispt.zs[npont] = outptt.z
                            ispt.ts[npont] = outptt.st
                            ispt.es[npont] = esec
                            ispt.ipl[npont] = iplane
                            ncltmp = npont
                            ntmpflg = 1
                            # RANDOMISE SECONDARY ELECTRON DIRECTION
                            # R3=drand48(RDUM)
                            # F3=1.0D0-2.0D0*R3
                            # THETA0=DACOS(F3)
                            # F6=DCOS(THETA0)
                            # F5=DSIN(THETA0)
                            # R4=drand48(rdum)
                            # PHI0=F4*R4
                            # F8=DSIN(PHI0)
                            # F9=DCOS(PHI0)
                            # DCX(NPONT)=F9*F5
                            # DCY(NPONT)=F8*F5
                            # DCZ(NPONT)=F6
                            # STORE POSSIBLE SHELL EMISSIONS AUGER OR FLUORESCENCE
                            if eistr > 30.0:
                              # AUGER EMISSIONS WITHOUT FLUORESCENCE
                              naug = ionfl.nc0[kgas, i]
                              eavaug = ionfl.ec0[kgas, i] / float64(naug)
                              fdo jfl, 1, naug:
                                nclus = nclus + 1
                                npont = npont + 1
                                ispt.xs[npont] = outptt.x
                                ispt.ys[npont] = outptt.y
                                ispt.zs[npont] = outptt.z
                                ispt.ts[npont] = outptt.st
                                ispt.es[npont] = eavaug
                                ispt.ipl[npont] = iplane
                                r3 = drand48(rng)
                                f3 = 1.0 - 2.0 * r3
                                theta0 = arccos(f3)
                                f6 = cos(theta0)
                                f5 = sin(theta0)
                                r4 = drand48(rng)
                                phi0 = f4 * r4
                                f8 = sin(phi0)
                                f9 = cos(phi0)
                                ispt.dcx[npont] = f9 * f5
                                ispt.dcy[npont] = f8 * f5
                                ispt.dcz[npont] = f6
                          # GENERATE SCATTERING ANGLES AND UPDATE  LABORATORY COSINES AFTER
                          # COLLISION ALSO UPDATE ENERGY OF ELECTRON.
                          ipt = gasTable[kgas].iarry[i]
                          id = id + 1
                          iter = iter + 1
                          iprint = iprint + 1
                          outptt.icoll[kgas, ipt] = outptt.icoll[kgas, ipt] + 1
                          outptt.icoln[kgas, i] = outptt.icoln[kgas, i] + 1
                          block L5:
                            # IF EXCITATION THEN ADD PROBABILITY ,PENFRA(KGAS,1,I),OF TRANSFER TO
                            # GIVE IONISATION OF THE OTHER GASES IN THE MIXTURE
                            if inpt.ipen == 0: break L5
                            if gasTable[kgas].penfra[1, i] != 0.0:
                              ran = drand48(rng)
                              if ran > gasTable[kgas].penfra[1, i]: break L5
                              nclus = nclus + 1
                              npont = npont + 1
                              if npont > 8000:
                                fwrite(fmt546, npont, iter)
                                fstop()
                              block L667:
                                # ENTER HERE POSSIBLE DELOCALISATION LENGTH FOR PENNING TRANSFER
                                if gasTable[kgas].penfra[2, i] == 0.0:
                                  ispt.xs[npont] = outptt.x
                                  ispt.ys[npont] = outptt.y
                                  ispt.zs[npont] = outptt.z
                                  break L667
                                asign = 1.0
                                ran = drand48(rng)
                                ran1 = drand48(rng)
                                if ran1 < 0.5: asign = -asign
                                ispt.xs[npont] = outptt.x - ln(ran) * gasTable[kgas].penfra[2, i] * asign
                                ran = drand48(rng)
                                ran1 = drand48(rng)
                                if ran1 < 0.5: asign = -asign
                                ispt.ys[npont] = outptt.y - ln(ran) * gasTable[kgas].penfra[2, i] * asign
                                ran = drand48(rng)
                                ran1 = drand48(rng)
                                if ran1 < 0.5: asign = -asign
                                ispt.zs[npont] = outptt.z - ln(ran) * gasTable[kgas].penfra[2, i] * asign
                              # ENTER POSSIBLE PENNING TRANSFER TIME
                              tpen = outptt.st
                              if not (gasTable[kgas].penfra[3, i] == 0.0):
                                ran = drand48(rng)
                                tpen = outptt.st - ln(ran) * gasTable[kgas].penfra[3, i]
                              ispt.ts[npont] = tpen
                              ispt.es[npont] = 1.0
                              ispt.dcx[npont] = dcx1
                              ispt.dcy[npont] = dcy1
                              ispt.dcz[npont] = dcz1
                              # FIND LAST TIME PLANE BEFORE TPEN
                              tstop1 = 0.0
                              iplane1 = 0
                              block L670:
                                fdo kdum, 1, cion.itfinal:
                                  tstop1 = cion.tstep + tstop1
                                  if tpen < tstop1: break L670
                                  iplane1 = iplane1 + 1
                                # PENNING TRANSFER OCCURS AFTER FINAL TIME PLANE. CLEAR ENTRY
                                npont = npont - 1
                                nclus = nclus - 1
                                break L5
                              ispt.ipl[npont] = iplane1
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
                          # IF(EOK.LT.EI) EI=0.0D0
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
                          block L190:
                            if argz == 0.0:
                              dcz1 = f6
                              dcx1 = f9 * f5
                              dcy1 = f8 * f5
                              if ntmpflg == 1:
                                # USE FREE KINEMATICS FOR IONISATION SECONDARY ANGLE
                                f5s = f5 * sqrt(e1 / ispt.es[ncltmp])
                                if f5s > 1.0: f5s = 1.0
                                thsec = arcsin(f5s)
                                f5s = sin(thsec)
                                f6s = cos(thsec)
                                if f6 < 0.0: f6s = -f6s
                                phis = phi0 + setp.api
                                if phis > f4: phis = phi0 - f4
                                f8s = sin(phis)
                                f9s = cos(phis)
                                ispt.dcz[ncltmp] = f6s
                                ispt.dcy[ncltmp] = f8s * f5s
                                ispt.dcx[ncltmp] = f9s * f5s
                                ntmpflg = 0
                              break L190
                            dcz1 = dzcom * f6 + argz * f5 * f8
                            dcy1 = dycom * f6 + f5 / argz * (dxcom * f9 - dycom * dzcom * f8)
                            dcx1 = dxcom * f6 - f5 / argz * (dycom * f9 + dxcom * dzcom * f8)
                            if ntmpflg == 1:
                              # USE FREE KINEMATICS FOR IONISATION SECONDARY ANGLE
                              f5s = f5 * sqrt(e1 / ispt.es[ncltmp])
                              if f5s > 1.0: f5s = 1.0
                              thsec = arcsin(f5s)
                              f5s = sin(thsec)
                              f6s = cos(thsec)
                              if f6 < 0.0: f6s = -f6s
                              phis = phi0 + setp.api
                              if phis > f4: phis = phi0 - f4
                              f8s = sin(phis)
                              f9s = cos(phis)
                              ispt.dcz[ncltmp] = dzcom * f6s + argz * f5s * f8s
                              ispt.dcy[ncltmp] = dycom * f6s + f5s / argz * (dxcom * f9s - dycom * dzcom * f8s)
                              ispt.dcx[ncltmp] = dxcom * f6s - f5s / argz * (dycom * f9s + dxcom * dzcom * f8s)
                              ntmpflg = 0
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
                          # STORE DIRECTION COSINES AND ENERGY AFTER N COLLISIONS
                          # FOR LATER REUSE IN PRIMARY GENERATION
                          i100 = i100 + 1
                          if i100 == 200:
                            dcz100 = dcz1
                            dcx100 = dcx1
                            dcy100 = dcy1
                            e100 = e1
                            i100 = 0
                          if iprint > jprint: break L200
                          break R1
                        # TAKE ELECTRONS FROM STORE
                        outptt.x = ispt.xs[npont]
                        outptt.y = ispt.ys[npont]
                        outptt.z = ispt.zs[npont]
                        outptt.st = ispt.ts[npont]
                        e1 = ispt.es[npont]
                        dcx1 = ispt.dcx[npont]
                        dcy1 = ispt.dcy[npont]
                        dcz1 = ispt.dcz[npont]
                        vtot = const9 * sqrt(e1)
                        cx1 = dcx1 * vtot
                        cy1 = dcy1 * vtot
                        cz1 = dcz1 * vtot
                        iplane = ispt.ipl[npont]
                        npont = npont - 1
                        zstrt = outptt.z
                        ystrt = outptt.y
                        xstrt = outptt.x
                        tsstrt = outptt.st
                        break R555
                      # INTERMEDIATE PRINTOUT
                      iprint = 0
                      wz = ttrm.ztots / ttrm.ttots
                      wz = wz * 1000000000.0
                      wy = ytots / ttrm.ttots
                      wy = wy * 1000000000.0
                      wx = xtots / ttrm.ttots
                      wx = wx * 1000000000.0
                      jct = id div 100000
                      # IF(J1.EQ.1) WRITE(6,201)
                      # 201 FORMAT(/,7X,'INTERMEDIATE OUTPUT',/,'    VELZ     VELY     VELX
                      # /   TIME       COUNT    ')
                      # ROTATE INTERMEDIATE OUTPUT INTO LAB FRAME
                      wzr = wz * rots.rcs - wx * rots.rsn
                      wyr = wy
                      wxr = wz * rots.rsn + wx * rots.rcs
                      # WRITE(6,202) WZR,WYR,WXR,TTOTS,JCT
                      # 202 FORMAT(3(1X,F8.3),1X,D10.3,4X,I6)
                      j1 = j1 + 1
                      break R1
                      break R1x
                break R555x
        break R544x
  # MAIN LOOP END
  outptt.xid = float64(id)
  if nelec > cion.iprim:
    aneion = float64(neion)
    anbt = float64(nelec - cion.iprim)
    tplouth.attoint = aneion / anbt
    tplouth.attert = sqrt(aneion) / aneion
    tplouth.aioert = sqrt(anbt) / anbt
  else:
    aneion = float64(neion)
    tplouth.attoint = (-1.0)
    tplouth.attert = sqrt(aneion) / aneion
  jct = id div 100000
  if j1 == 1:
    fwrite(fmt940, nclus, iter)
    fstop()
  fwrite(fmt878, nelec, neion, cion.iprim)
  eprmbar = 0.0
  e2prm = 0.0
  if cion.iprim == 1: return
  fdo i, 1, cion.iprim:
    e2prm = e2prm + eprm[i] * eprm[i]
    eprmbar = eprmbar + eprm[i]
  ebar = eprmbar / float64(cion.iprim)
  eerr = sqrt(e2prm / float64(cion.iprim) - ebar*ebar)
  fwrite(fmt836, ebar, eerr)
  # WRITE(6,835) (IESPECP(J),J=1,100)
  # 835  FORMAT(/,2X,'ENERGY SPECTRUM OF PRIMARY ELECTRONS IN 1 EV. BINS',/
  # /,10(2X,10I5,/))
  return
  if iter > setp.nmax:
    fwrite(fmt991, iter, setp.nmax, npont, nelec, cion.iprim, nmxadd)
    fstop()
  return
