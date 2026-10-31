## SUBROUTINE MONTEFTG of magboltz-11.21.f (line 15817).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import dmin0, sort, tplaneg

const
  fmt546 = "(2X,\' PROGRAM STOPPED . NPONT=\',I5,\' ITER=\',I10)"
  fmt836 = "(/,2X,\'AVERAGE ENERGY OF PRIMARY ELECTRON =\',F10.3,\' EV.\',/,\'   ENERGY SPREAD OF PRIMARY ELECTRON =\',F10.3,\' EV.\')"
  fmt878 = "(/,\' TOTAL NO OF ELECTRONS=\',I8,/,\' TOTAL NO OF NEG. IONS=\',I8,/,\' TOTAL NO OF PRIMARIES=\',I8)"
  fmt913 = "(3X,\' AFTER STORE ITER=\',I10,\' DZ=\',D12.3,\'E1=\',D12.3,\' COSWT=\',D12.3,\' SINWT=\',D12.3,\' WBT=\',D12.3,\' CY1=\',D12.3)"
  fmt940 = "(2(/),\' PROGRAM STOPPED  (TOO FEW COLLISIONS), DECREASE THE ESTIMATED ALPHA.  NCLUS = \',I7,\'  ITER =\',I9)"
  fmt944 = "(/,2X,\'PROGRAM STOPPED TOO MANY PRIMARIES IPRIM =\',I7)"
  fmt991 = "(2(/),\' PROGRAM STOPPED.  ITER =\',I10,\'    NMAX =\',I10,/,  \' NPONT=\',I4,\' NELEC=\',I8,\' IPRIM=\',I4,\' NMXADD=\',I3)"

var eprm: array[1..10000000, float64]

proc monteftg*() =
  var iespecp: array[1..100, int]
  var temp: array[1..4000, float64]
  var s, ytot, ytots, rdum, e1, const9, zstrt, ystrt: float64
  var tsstrt, dcz1, dcx1, dcy1, vtot, cx1, cy1, cz1: float64
  var e100, dcz100, dcx100, dcy100, f4, tlim, tdash, tstop: float64
  var r1, t, wbt, coswt, sinwt, dz, ev, r5: float64
  var ttem, t2, cx2, cy2, cz2, dcx2, dcy2, dcz2: float64
  var r2, s1, ei, eistr, r9, esec, eavaug, r3: float64
  var f3, theta0, f6, f5, r4, phi0, f8, f9: float64
  var ran, asign, ran1, tpen, tstop1, s2, r31, epsi: float64
  var arg1, d, q, u, csqd, argz, f5s, thsec: float64
  var f6s, phis, f8s, f9s, w, wy, aneion, anbt: float64
  var eprmbar, e2prm, ebar, eerr: float64
  var i, id, n4000, n300, n100, intem, i100, ncol: int
  var nelec, neion, nmxadd, ntmpflg, npont, nclus, j1, jprint: int
  var iprint, iter, iplane, idum, ie, it, ipt, ncltmp: int
  var naug, jfl, iplane1, kdum, jct: int
  # CALCULATES COLLISION EVENTS AND UPDATES DIFFUSION AND VELOCITY.
  # THIS ROUTINE HANDLES TERMINATIONS AT FIXED DRIFT TIMES.
  # B FIELD AT 90 DEGREES TO EFIELD
  s = 0.0
  outpt.st = 0.0
  outpt.x = 0.0
  outpt.y = 0.0
  outpt.z = 0.0
  ttrm.ztot = 0.0
  ttrm.ztots = 0.0
  ytot = 0.0
  ytots = 0.0
  ttrm.ttot = 0.0
  ttrm.ttots = 0.0
  setp.small = 1e-20
  outpt.tmax1 = 0.0
  rdum = setp.rstart
  e1 = setp.estart
  const9 = cnsts1.const3 * 0.01
  setp.api = arccos((-1.0))
  fdo i, 1, 300:
    outpt.time[i] = 0.0
  fdo i, 1, 30:
    outpt.icoll[i] = 0
  fdo i, 1, 960:
    outpt.icoln[i] = 0
  fdo i, 1, 4000:
    outpt.spec[i] = 0.0
  fdo i, 1, 100:
    iespecp[i] = 0
  fdo i, 1, 8:
    tploutg.etpl[i] = 0.0
    tploutg.xtpl[i] = 0.0
    tploutg.ytpl[i] = 0.0
    tploutg.ztpl[i] = 0.0
    tploutg.ttpl[i] = 0.0
    tploutg.xxtpl[i] = 0.0
    tploutg.yytpl[i] = 0.0
    tploutg.zztpl[i] = 0.0
    tploutg.yztpl[i] = 0.0
    tploutg.vztpl[i] = 0.0
    tploutg.vytpl[i] = 0.0
    tploutg.netpl[i] = 0
  id = 0
  n4000 = 4000
  n300 = 300
  n100 = 100
  intem = 8
  i100 = 0
  ncol = 0
  outpt.nnull = 0
  nelec = 0
  neion = 0
  nmxadd = 0
  ntmpflg = 0
  npont = 0
  nclus = 0
  j1 = 1
  zstrt = 0.0
  ystrt = 0.0
  tsstrt = 0.0
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
  f4 = 2.0 * setp.api
  # SET TO MAXIMUM POSSIBLE COLLISION FREQ.
  tlim = large.tcfmax[1]
  fdo i, 2, intem:
    if tlim < large.tcfmax[i]: tlim = large.tcfmax[i]
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
            outpt.x = 0.0
            outpt.y = 0.0
            outpt.z = 0.0
            dcz1 = dcz100
            dcx1 = dcx100
            dcy1 = dcy100
            e1 = e100
            vtot = const9 * sqrt(e1)
            cx1 = dcx1 * vtot
            cy1 = dcy1 * vtot
            cz1 = dcz1 * vtot
            nclus = nclus + 1
            outpt.st = 0.0
            tsstrt = 0.0
            zstrt = 0.0
            ystrt = 0.0
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
                      t = -(ln(r1) / tlim) + tdash
                      tdash = t
                      wbt = bfld.wb * t
                      coswt = cos(wbt)
                      sinwt = sin(wbt)
                      block L200:
                        block L20:
                          block R15x:
                            while true:
                              block R15:
                                if t + outpt.st >= tstop:
                                  iplane = iplane + 1
                                  tstop = tstop + cion.tstep
                                  # STORE POSITION AND ENERGY AT TIME PLANE =IPLANE.
                                  tplaneg(t, e1, cx1, cy1, cz1, setp.efield, iplane)
                                  # CHECK IF PASSED THROUGH MORE THAN ONE PLANE IN THIS STEP
                                  if t + outpt.st >= tstop and tstop <= cion.tfinal: break R15
                                  if t + outpt.st >= cion.tfinal:
                                    ttrm.ztot = ttrm.ztot + outpt.z
                                    ytot = ytot + outpt.y
                                    ttrm.ttot = ttrm.ttot + outpt.st
                                    ttrm.ztots = ttrm.ztots + outpt.z - zstrt
                                    ytots = ytots + outpt.y - ystrt
                                    ttrm.ttots = ttrm.ttots + outpt.st - tsstrt
                                    tstop = cion.tstep
                                    # NO MORE ELECTRONS IN CASCADE TRY NEW PRIMARY ELECTRON
                                    if nelec == nclus + 1: break R544
                                    break L20
                                break R15x
                          dz = (cz1 * sinwt + (bfld.eovb - cy1) * (1.0 - coswt)) / bfld.wb
                          ev = e1 + dz * setp.efield * 100.0
                          if ev < 0.0:
                            fwrite(fmt913, iter, dz, ev, coswt, sinwt, wbt, cy1)
                            ev = 0.001
                          ie = int(trunc(ev / inpt.estep) + 1.0)
                          ie = int(dmin0(ie, n4000))
                          # TEST FOR REAL OR NULL COLLISION
                          r5 = drand48(rng)
                          ttem = large.tcf[ie] / tlim
                          if not (r5 <= ttem):
                            outpt.nnull = outpt.nnull + 1
                            break R1
                          # CALCULATE DIRECTION COSINES AND POSITIONS AT INSTANT BEFORE COLLISION
                          # ALSO UPDATE DIFFUSION  AND ENERGY CALCULATIONS.
                          t2 = t * t
                          if t >= outpt.tmax1: outpt.tmax1 = t
                          tdash = 0.0
                          # CALC VELOCITY
                          cx2 = cx1
                          cy2 = (cy1 - bfld.eovb) * coswt + cz1 * sinwt + bfld.eovb
                          cz2 = cz1 * coswt - (cy1 - bfld.eovb) * sinwt
                          # CALC DIRECTION COSINES
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
                          # DETERMINATION OF REAL COLLISION TYPE
                          r2 = drand48(rng)
                          # FIND LOCATION WITHIN 4 UNITS IN COLLISION ARRAY
                          sort(i, r2, ie)
                          while true:
                            i = i + 1
                            if not (large.cf[ie, i] < r2): break
                          s1 = large.rgas[i]
                          ei = large.ein[i]
                          if ev < ei:
                            # WRITE(6,994) E,EI,ITER
                            # 994  FORMAT(2X,' WARNING ENERGY =',F8.3,' LESS THAN ENERGY LOSS EI=',F8
                            # /.3,' AT ITER=',I12,' DUE TO BINNING ERROR')
                            # FIX ENERGY LOSS SMALLER THAN INCIDENT ENERGY IF A BINING ERROR OCCURS.
                            ei = ev - 0.0001
                          # NEW CATCH FOR COLLISIONS WITH LARGE VIBRATIONAL ENERGY LOSS:
                          # ENERGY LOSS = INCIDENT ENERGY - 0.024999
                          if abs(ei - 0.024999000132083893) <= 0.0001:
                            ei = ev - 0.024999000132083893
                          if not (large.ipn[i] == 0):
                            # ATTACHMENT
                            if large.ipn[i] == (-1):
                              neion = neion + 1
                              ipt = large.iarry[i]
                              id = id + 1
                              iter = iter + 1
                              iprint = iprint + 1
                              outpt.icoll[ipt] = outpt.icoll[ipt] + 1
                              outpt.icoln[i] = outpt.icoln[i] + 1
                              it = int(trunc(t + 1.0))
                              it = int(dmin0(it, n300))
                              outpt.time[it] = outpt.time[it] + 1.0
                              ttrm.ztot = ttrm.ztot + outpt.z
                              ytot = ytot + outpt.y
                              ttrm.ttot = ttrm.ttot + outpt.st
                              ttrm.ztots = ttrm.ztots + outpt.z - zstrt
                              ytots = ytots + outpt.y - ystrt
                              ttrm.ttots = ttrm.ttots + outpt.st - tsstrt
                              # ELECTRON CAPTURED START NEW PRIMARY
                              if nelec == nclus + 1: break R544
                              # ELECTRON CAPTURED TAKE NEXT ELECTRON FROM STORE
                              break L20
                            eistr = ei
                            r9 = drand48(rng)
                            # ESEC=R9*(E-EI)
                            # USE OPAL PETERSON AND BEATY SPLITTING FACTOR.
                            esec = large.wpl[i] * tan(r9 * arctan((ev - ei) / (2.0 * large.wpl[i])))
                            esec = large.wpl[i] * pow(esec / large.wpl[i], 0.9524000287055969)
                            ei = esec + ei
                            # STORE POSITION ,ENERGY, DIRECTION COSINES AND TIME OF GENERATION
                            # OF IONISATION ELECTRON
                            nclus = nclus + 1
                            npont = npont + 1
                            nmxadd = max(npont, nmxadd)
                            if npont > 8000:
                              fwrite(fmt546, npont, iter)
                              fstop()
                            ispt.xs[npont] = outpt.x
                            ispt.ys[npont] = outpt.y
                            ispt.zs[npont] = outpt.z
                            ispt.ts[npont] = outpt.st
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
                              # AUGER EMISSION WITHOUT FLUORESCENCE
                              naug = ionfl1.nc0[i]
                              eavaug = ionfl1.ec0[i] / float64(naug)
                              fdo jfl, 1, naug:
                                nclus = nclus + 1
                                npont = npont + 1
                                ispt.xs[npont] = outpt.x
                                ispt.ys[npont] = outpt.y
                                ispt.zs[npont] = outpt.z
                                ispt.ts[npont] = outpt.st
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
                          ipt = large.iarry[i]
                          id = id + 1
                          iter = iter + 1
                          iprint = iprint + 1
                          outpt.icoll[ipt] = outpt.icoll[ipt] + 1
                          outpt.icoln[i] = outpt.icoln[i] + 1
                          block L5:
                            # IF EXCITATION THEN ADD PROBABILITY ,PENFRA(1,I),OF TRANSFER TO GIVE
                            # IONISATION OF THE OTHER GASES IN THE MIXTURE
                            if inpt.ipen == 0: break L5
                            # NIMBOLTZ FIX MB-28: WAS .EQ.0.0 (INVERTED; ALL OTHER MONTE ROUTINES
                            # USE .NE.)
                            if large.penfra[1, i] != 0.0:
                              ran = drand48(rng)
                              if ran > large.penfra[1, i]: break L5
                              nclus = nclus + 1
                              npont = npont + 1
                              if npont > 8000:
                                fwrite(fmt546, npont, iter)
                                fstop()
                              block L667:
                                # ENTER HERE POSSIBLE DELOCALISATION LENGTH FOR PENNING TRANSFER
                                if large.penfra[2, i] == 0.0:
                                  ispt.xs[npont] = outpt.x
                                  ispt.ys[npont] = outpt.y
                                  ispt.zs[npont] = outpt.z
                                  break L667
                                asign = 1.0
                                ran = drand48(rng)
                                ran1 = drand48(rng)
                                if ran1 < 0.5: asign = -asign
                                ispt.xs[npont] = outpt.x - ln(ran) * large.penfra[2, i] * asign
                                ran = drand48(rng)
                                ran1 = drand48(rng)
                                if ran1 < 0.5: asign = -asign
                                ispt.ys[npont] = outpt.y - ln(ran) * large.penfra[2, i] * asign
                                ran = drand48(rng)
                                ran1 = drand48(rng)
                                if ran1 < 0.5: asign = -asign
                                ispt.zs[npont] = outpt.z - ln(ran) * large.penfra[2, i] * asign
                              # POSSIBLE PENNING TRANSFER TIME
                              tpen = outpt.st
                              if not (large.penfra[3, i] == 0.0):
                                ran = drand48(rng)
                                tpen = outpt.st - ln(ran) * large.penfra[3, i]
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
                                # PENNING TRANSFER OCCURS AFTER FINAL TIME PLANE CLEAR ENTRY
                                npont = npont - 1
                                nclus = nclus - 1
                                break L5
                              ispt.ipl[npont] = iplane1
                          s2 = s1 * s1 / (s1 - 1.0)
                          # ANISOTROPIC SCATTERING
                          r3 = drand48(rng)
                          if anis.index[i] == 1:
                            r31 = drand48(rng)
                            f3 = 1.0 - r3 * anis.angct[ie, i]
                            if r31 > anis.psct[ie, i]: f3 = -f3
                          elif anis.index[i] == 2:
                            epsi = anis.psct[ie, i]
                            f3 = 1.0 - 2.0 * r3 * (1.0 - epsi) / (1.0 + epsi * (1.0 - 2.0 * r3))
                          else:
                            # ISOTROPIC SCATTERING
                            f3 = 1.0 - 2.0 * r3
                          theta0 = arccos(f3)
                          r4 = drand48(rng)
                          phi0 = f4 * r4
                          f8 = sin(phi0)
                          f9 = cos(phi0)
                          if ev < ei: ei = 0.0
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
                            dcz1 = dcz2 * f6 + argz * f5 * f8
                            dcy1 = dcy2 * f6 + f5 / argz * (dcx2 * f9 - dcy2 * dcz2 * f8)
                            dcx1 = dcx2 * f6 - f5 / argz * (dcy2 * f9 + dcx2 * dcz2 * f8)
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
                              ispt.dcz[ncltmp] = dcz2 * f6s + argz * f5s * f8s
                              # NIMBOLTZ FIX MB-24: WAS DCYS (UNDEFINED SCALAR), AS IN
                              # MONTEFT/MONTEFTA: DCY2
                              ispt.dcy[ncltmp] = dcy2 * f6s + f5s / argz * (dcx2 * f9s - dcy2 * dcz2 * f8s)
                              ispt.dcx[ncltmp] = dcx2 * f6s - f5s / argz * (dcy2 * f9s + dcx2 * dcz2 * f8s)
                              ntmpflg = 0
                          vtot = const9 * sqrt(e1)
                          cx1 = dcx1 * vtot
                          cy1 = dcy1 * vtot
                          cz1 = dcz1 * vtot
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
                        outpt.x = ispt.xs[npont]
                        outpt.y = ispt.ys[npont]
                        outpt.z = ispt.zs[npont]
                        outpt.st = ispt.ts[npont]
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
                        zstrt = outpt.z
                        ystrt = outpt.y
                        tsstrt = outpt.st
                        break R555
                      # INTERMEDIATE PRINTOUT
                      iprint = 0
                      w = ttrm.ztots / ttrm.ttots
                      w = w * 1000000000.0
                      wy = ytots / ttrm.ttots
                      wy = wy * 1000000000.0
                      jct = id div 100000
                      # IF(J1.EQ.1) WRITE(6,201)
                      # 201 FORMAT(/,7X,'INTERMEDIATE OUTPUT',/,'    VELZ     POSZ       TIME
                      # /      VELY        COUNT    ')
                      # WRITE(6,202) W,ZTOTS,TTOTS,WY,JCT
                      # 202 FORMAT(1X,F8.3,3(1X,D10.3),4X,I6)
                      j1 = j1 + 1
                      break R1
                      break R1x
                break R555x
        break R544x
  # MAIN LOOP END
  outpt.xid = float64(id)
  if nelec > cion.iprim:
    aneion = float64(neion)
    anbt = float64(nelec - cion.iprim)
    tploutg.attoint = aneion / anbt
    tploutg.attert = sqrt(aneion) / aneion
    tploutg.aioert = sqrt(anbt) / anbt
  else:
    aneion = float64(neion)
    tploutg.attoint = (-1.0)
    tploutg.attert = sqrt(aneion) / aneion
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
