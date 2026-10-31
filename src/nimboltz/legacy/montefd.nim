## SUBROUTINE MONTEFD of magboltz-11.21.f (line 14740).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import dmin0, sort, splane, tcalc

const
  fmt546 = "(2X,\' PROGRAM STOPPED NPONT=\',I5,\' ITER=\',I10)"
  fmt836 = "(/,2X,\'AVERAGE ENERGY OF PRIMARY ELECTRON =\',F10.3,\' EV.\',/,\'   ENERGY SPREAD OF PRIMARY ELECTRON =\',F10.3,\' EV.\')"
  fmt878 = "(/,\' TOTAL NO OF ELECTRONS=\',I10,/,\' TOTAL NO OF NEG. IONS=\',I10,/,\' TOTAL NO OF PRIMARIES=\',I10)"
  fmt940 = "(2(/),\' PROGRAM STOPPED (TOO FEW COLLISIONS), DECREASE THE  ESTIMATED ALPHA.  NCLUS=\',I7,\'  ITER =\',I9,\' NELEC=\',I9,\' NEION =\',I6)"
  fmt944 = "(2X,\' PROGRAM STOPPED TOO MANY PRIMARIES IPRIM=\',I7)"
  fmt991 = "(2(/),\' PROGRAM STOPPED.  ITER =\',I10,\'    NMAX =\',I10,/,  \' NPONT=\',I4,\' NELEC=\',I8,\' IPRIM=\',I6,\' NMXADD=\',I3)"
  fmt999 = "(2X,\' WARNING ENERGY LT.0. E=\',D12.3,\' E1=\',D12.3,\' AP=\',D12.3,\' BP=\',D12.3,\' T=\',D12.3,/,\' DCZ1=\',D12.3,\' ITER=\',I10)"

var eprm: array[1..10000000, float64]

proc montefd*() =
  var iespecp: array[1..100, int]
  var temp: array[1..4000, float64]
  var s, rdum, e1, const9, zstrt, tsstrt, dcz1, dcx1: float64
  var dcy1, dcx100, dcy100, dcz100, e100, bp, f1, f2: float64
  var f4, tlim, tzstop, tdash, r1, t, told, ap: float64
  var tlft, tzstop1, ev, r5, ttem, t2, const6, dcx2: float64
  var dcy2, dcz2, a, b, const7, cz1, cx1, cy1: float64
  var r2, s1, ei, eistr, r9, esec, eavaug, r3: float64
  var f3, theta0, f6, f5, r4, phi0, f8, f9: float64
  var ran, asign, ran1, tpen, s2, r31, epsi, arg1: float64
  var d, q, u, csqd, argz, f5s, thsec, f6s: float64
  var phis, f8s, f9s, epot, w, aneion, anbt, eprmbar: float64
  var e2prm, ebar, eerr: float64
  var i, id, n4000, n300, n100, intem, i100, ncol: int
  var nelec, neion, nmxadd, ntmpflg, npont, nclus, j1, jprint: int
  var iprint, iter, izplane, idum, isol, ie, it, ipt: int
  var idm1, ncltmp, naug, jfl, jct: int
  # CALCULATES COLLISION EVENTS AND UPDATES DIFFUSION AND VELOCITY.
  # THIS ROUTINE HANDLES TERMINATIONS AT FIXED DRIFT DISTANCES.
  s = 0.0
  outpt.st = 0.0
  outpt.x = 0.0
  outpt.y = 0.0
  outpt.z = 0.0
  ttrm.ztot = 0.0
  ttrm.ztots = 0.0
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
    splout.espl[i] = 0.0
    splout.xspl[i] = 0.0
    splout.yspl[i] = 0.0
    splout.zspl[i] = 0.0
    splout.tspl[i] = 0.0
    splout.xxspl[i] = 0.0
    splout.yyspl[i] = 0.0
    splout.zzspl[i] = 0.0
    splout.vzspl[i] = 0.0
    splout.tssum[i] = 0.0
    splout.tssum2[i] = 0.0
    spl1.tmspl[i] = 0.0
    spl1.ttmspl[i] = 0.0
    spl1.rspl[i] = 0.0
    spl1.rrspl[i] = 0.0
    spl1.rrspm[i] = 0.0
    splout.nesst[i] = 0
  splout.nesst[9] = 0
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
  tsstrt = 0.0
  # INITIAL DIRECTION COSINES
  dcz1 = cos(setp.theta)
  dcx1 = sin(setp.theta) * cos(setp.phi)
  dcy1 = sin(setp.theta) * sin(setp.phi)
  dcx100 = dcx1
  dcy100 = dcy1
  dcz100 = dcz1
  e100 = e1
  bp = setp.efield * setp.efield * cnsts1.const1
  f1 = setp.efield * cnsts1.const2
  f2 = setp.efield * cnsts1.const3
  f4 = 2.0 * setp.api
  # SET TO MAXIMUM POSSIBLE COLLISION FREQ.
  tlim = large.tcfmax[1]
  fdo i, 2, intem:
    if tlim < large.tcfmax[i]: tlim = large.tcfmax[i]
  jprint = setp.nmax div 10
  iprint = 0
  iter = 0
  cion.iprim = 0
  block R544x:
    while true:
      block R544:
        # LOOP FOR NEW STARTING ELECTRONS
        cion.iprim = cion.iprim + 1
        izplane = 0
        tzstop = 1000.0
        block L700:
          if cion.iprim > 1:
            # CHECK IF PROGRAM WILL EXCEED MAXIMUM NUMBER OF ITERATIONS IN THIS
            # CYCLE , IF SO OUTPUT CURRENT RESULTS.
            if iter > setp.nmax:
              cion.iprim = cion.iprim - 1
              break L700
            outpt.x = 0.0
            outpt.y = 0.0
            outpt.z = 0.0
            dcz1 = dcz100
            dcx1 = dcx100
            dcy1 = dcy100
            nclus = nclus + 1
            e1 = e100
            outpt.st = 0.0
            zstrt = 0.0
            tsstrt = 0.0
          if cion.iprim > 10000000:
            fwrite(fmt944, cion.iprim)
            fstop()
          eprm[cion.iprim] = e1
          idum = int(trunc(e1) + 1.0)
          idum = int(dmin0(idum, n100))
          iespecp[idum] = iespecp[idum] + 1
          block R555x:
            while true:
              block R555:
                # START OF LOOP FOR NEWLY CREATED ELECTRONS .
                tdash = 0.0
                nelec = nelec + 1
                block R1x:
                  while true:
                    block R1:
                      # MAIN LOOP
                      # IF(ITER.GT.NMAX)  GO TO 315
                      r1 = drand48(rng)
                      t = -(ln(r1) / tlim) + tdash
                      told = tdash
                      tdash = t
                      ap = dcz1 * f2 * sqrt(e1)
                      block L200:
                        block L20:
                          block R15x:
                            while true:
                              block R15:
                                if t >= tzstop and told < tzstop:
                                  tlft = tzstop
                                  # STORE POSITION AND ENERGY AT Z PLANE = IZPLANE.
                                  # IF(IZPLANE.EQ.0) WRITE(6,8876) IZPLANE
                                  # 8876  FORMAT(' IZPLANE=',I5)
                                  splane(t, e1, dcx1, dcy1, dcz1, ap, bp, setp.efield, tlft, izplane)
                                  # CHANGE  IF STATEMENT FROM (IZFINAL+1) TO (IZFINAL-1)
                                  # FOR ANODE TERMINATION .
                                  if izplane >= ctcalc.izfinal + 1:
                                    ttrm.ztot = ttrm.ztot + outpt.z
                                    ttrm.ttot = ttrm.ttot + outpt.st
                                    ttrm.ztots = ttrm.ztots + outpt.z - zstrt
                                    ttrm.ttots = ttrm.ttots + outpt.st - tsstrt
                                    # NO MORE ELECTRONS IN CASCADE RETURN TO MAIN.
                                    if nelec == nclus + 1: break R544
                                    break L20
                                  # IF TWO SOLUTIONS REPEAT ENTRY FOR SECOND SOLUTION.
                                  if isol == 2:
                                    tzstop = tzstop1
                                    isol = 1
                                    break R15
                                break R15x
                          ev = e1 + (ap + bp * t) * t
                          if ev < 0.0:
                            fwrite(fmt999, ev, e1, ap, bp, t, dcz1, iter)
                            ev = 0.001
                          ie = int(trunc(ev / inpt.estep) + 1.0)
                          ie = int(dmin0(ie, n4000))
                          # TEST FOR REAL OR NULL COLLISION
                          r5 = drand48(rng)
                          ttem = mixTable.tcf[ie] / tlim
                          if not (r5 <= ttem):
                            outpt.nnull = outpt.nnull + 1
                            break R1
                          # CALCULATE DIRECTION COSINES AND POSITIONS AT INSTANT BEFORE COLLISION
                          # ALSO UPDATE DIFFUSION  AND ENERGY CALCULATIONS.
                          t2 = t * t
                          if t >= outpt.tmax1: outpt.tmax1 = t
                          tdash = 0.0
                          const6 = sqrt(e1 / ev)
                          dcx2 = dcx1 * const6
                          dcy2 = dcy1 * const6
                          dcz2 = dcz1 * const6 + setp.efield * t * cnsts1.const5 / sqrt(ev)
                          a = ap * t
                          b = bp * t2
                          const7 = const9 * sqrt(e1)
                          a = t * const7
                          ncol = ncol + 1
                          cz1 = dcz1 * const7
                          outpt.x = outpt.x + dcx1 * a
                          outpt.y = outpt.y + dcy1 * a
                          outpt.z = outpt.z + dcz1 * a + t2 * f1
                          outpt.st = outpt.st + t
                          it = int(trunc(t + 1.0))
                          it = int(dmin0(it, n300))
                          outpt.time[it] = outpt.time[it] + 1.0
                          cx1 = dcx1 * const7
                          cy1 = dcy1 * const7
                          outpt.spec[ie] = outpt.spec[ie] + 1.0
                          # DETERMINATION OF REAL COLLISION TYPE
                          r2 = drand48(rng)
                          # FIND LOCATION WITHIN 4 UNITS IN COLLISION ARRAY
                          sort(i, r2, ie)
                          while true:
                            i = i + 1
                            if not (mixTable.cf[ie, i] < r2): break
                          s1 = mixTable.rgas[i]
                          ei = mixTable.ein[i]
                          if ev < ei:
                            # FIX ENERGY LOSS SMALLER THAN INCIDENT ENERGY IF A BINING ERROR OCCURS.
                            # WRITE(6,994) E,EI,ITER
                            # 994  FORMAT(3X,' WARNING BINNING ERROR ENERGY =',F8.3,' EI=',F8.3,' ITE
                            # /R =',I12)
                            ei = ev - 0.0001
                          # NEW CATCH FOR COLLISIONS WITH LARGE VIBRATIONAL ENERGY LOSS:
                          # ENERGY LOSS = INCIDENT ENERGY - 0.024999
                          if abs(ei - 0.024999000132083893) <= 0.0001:
                            ei = ev - 0.024999000132083893
                          if not (mixTable.ipn[i] == 0):
                            # ATTACHMENT
                            if mixTable.ipn[i] == (-1):
                              neion = neion + 1
                              ipt = mixTable.iarry[i]
                              id = id + 1
                              iter = iter + 1
                              iprint = iprint + 1
                              outpt.icoll[ipt] = outpt.icoll[ipt] + 1
                              outpt.icoln[i] = outpt.icoln[i] + 1
                              it = int(trunc(t + 1.0))
                              it = int(dmin0(it, n300))
                              outpt.time[it] = outpt.time[it] + 1.0
                              ttrm.ztot = ttrm.ztot + outpt.z
                              ttrm.ttot = ttrm.ttot + outpt.st
                              ttrm.ztots = ttrm.ztots + outpt.z - zstrt
                              ttrm.ttots = ttrm.ttots + outpt.st - tsstrt
                              idm1 = int(1.0 + trunc(outpt.z / cion.zstep))
                              if idm1 < 1: idm1 = 1
                              if idm1 > 9: idm1 = 9
                              splout.nesst[idm1] = splout.nesst[idm1] - 1
                              # ELECTRON CAPTURED START NEW PRIMARY
                              if nelec == nclus + 1: break R544
                              # ELECTRON CAPTURED TAKE NEXT ELECTRON FROM STORE
                              break L20
                            eistr = ei
                            r9 = drand48(rng)
                            # ESEC=R9*(E-EI)
                            # USE OPAL PETERSON AND BEATY SPLITTING FACTOR.
                            esec = mixTable.wpl[i] * tan(r9 * arctan((ev - ei) / (2.0 * mixTable.wpl[i])))
                            esec = mixTable.wpl[i] * pow(esec / mixTable.wpl[i], 0.9524000287055969)
                            ei = esec + ei
                            # STORE POSITION ,ENERGY, DIRECTION COSINES AND TIME OF GENERATION
                            # OF IONISATION ELECTRON
                            nclus = nclus + 1
                            npont = npont + 1
                            nmxadd = max(npont, nmxadd)
                            if npont > 8000:
                              fwrite(fmt546, npont, iter)
                              fstop()
                            ips.xss[npont] = outpt.x
                            ips.yss[npont] = outpt.y
                            ips.zss[npont] = outpt.z
                            ips.tss[npont] = outpt.st
                            ips.ess[npont] = esec
                            ntmpflg = 1
                            ncltmp = npont
                            # RANDOMISE SECONDARY ELECTRON DIRECTION
                            # R3=drand48(RDUM)
                            # F3=1.0D0-2.0D0*R3
                            # THETA0=DACOS(F3)
                            # F6=DCOS(THETA0)
                            # F5=DSIN(THETA0)
                            # R4=drand48(RDUM)
                            # PHI0=F4*R4
                            # F8=DSIN(PHI0)
                            # F9=DCOS(PHI0)
                            # DCXS(NPONT)=F9*F5
                            # DCYS(NPONT)=F8*F5
                            # DCZS(NPONT)=F6
                            idm1 = int(1.0 + trunc(outpt.z / cion.zstep))
                            if idm1 < 1: idm1 = 1
                            if idm1 > 9: idm1 = 9
                            ips.ipls[npont] = idm1
                            splout.nesst[ips.ipls[npont]] = splout.nesst[ips.ipls[npont]] + 1
                            # STORE POSSIBLE SHELL EMISSIONS AUGER OR FLUORESCENCE
                            if eistr > 30.0:
                              # AUGER EMISSION WITHOUT FLUORESCENCE
                              naug = ionfl1.nc0[i]
                              eavaug = ionfl1.ec0[i] / float64(naug)
                              fdo jfl, 1, naug:
                                nclus = nclus + 1
                                npont = npont + 1
                                ips.xss[npont] = outpt.x
                                ips.yss[npont] = outpt.y
                                ips.zss[npont] = outpt.z
                                ips.tss[npont] = outpt.st
                                ips.ess[npont] = eavaug
                                r3 = drand48(rng)
                                f3 = 1.0 - 2.0 * r3
                                theta0 = arccos(f3)
                                f6 = cos(theta0)
                                f5 = sin(theta0)
                                r4 = drand48(rng)
                                phi0 = f4 * r4
                                f8 = sin(phi0)
                                f9 = cos(phi0)
                                ips.dcxs[npont] = f9 * f5
                                ips.dcys[npont] = f8 * f5
                                ips.dczs[npont] = f6
                                idm1 = int(1.0 + trunc(outpt.z / cion.zstep))
                                if idm1 < 1: idm1 = 1
                                if idm1 > 9: idm1 = 9
                                ips.ipls[npont] = idm1
                                splout.nesst[ips.ipls[npont]] = splout.nesst[ips.ipls[npont]] + 1
                          # GENERATE SCATTERING ANGLES AND UPDATE  LABORATORY COSINES AFTER
                          # COLLISION ALSO UPDATE ENERGY OF ELECTRON.
                          ipt = mixTable.iarry[i]
                          id = id + 1
                          iter = iter + 1
                          iprint = iprint + 1
                          outpt.icoll[ipt] = outpt.icoll[ipt] + 1
                          outpt.icoln[i] = outpt.icoln[i] + 1
                          block L5:
                            # IF EXCITATION THEN ADD PROBABILITY,PENFRA(1,I), OF TRANSFER TO GIVE
                            # IONISATION OF THE OTHER GASES IN THE MIXTURE.
                            if inpt.ipen == 0: break L5
                            if mixTable.penfra[1, i] != 0.0:
                              ran = drand48(rng)
                              if ran > mixTable.penfra[1, i]: break L5
                              nclus = nclus + 1
                              npont = npont + 1
                              if npont > 8000:
                                fwrite(fmt546, npont, iter)
                                fstop()
                              block L669:
                                block L667:
                                  # ENTER HERE POSSIBLE DELOCALISATION LENGTH FOR PENNING TRANSFER
                                  if mixTable.penfra[2, i] == 0.0:
                                    ips.xss[npont] = outpt.x
                                    ips.yss[npont] = outpt.y
                                    ips.zss[npont] = outpt.z
                                    if ips.zss[npont] > cion.zfinal or ips.zss[npont] < 0.0: break L669
                                    break L667
                                  asign = 1.0
                                  ran = drand48(rng)
                                  ran1 = drand48(rng)
                                  if ran1 < 0.5: asign = -asign
                                  ips.xss[npont] = outpt.x - ln(ran) * mixTable.penfra[2, i] * asign
                                  ran = drand48(rng)
                                  ran1 = drand48(rng)
                                  if ran1 < 0.5: asign = -asign
                                  ips.yss[npont] = outpt.y - ln(ran) * mixTable.penfra[2, i] * asign
                                  ran = drand48(rng)
                                  ran1 = drand48(rng)
                                  if ran1 < 0.5: asign = -asign
                                  ips.zss[npont] = outpt.z - ln(ran) * mixTable.penfra[2, i] * asign
                                  if ips.zss[npont] < 0.0: break L669
                                  if ips.zss[npont] > cion.zfinal or ips.zss[npont] < 0.0: break L669
                                tpen = outpt.st
                                if not (mixTable.penfra[3, i] == 0.0):
                                  ran = drand48(rng)
                                  tpen = outpt.st - ln(ran) * mixTable.penfra[3, i]
                                ips.tss[npont] = tpen
                                ips.ess[npont] = 1.0
                                ips.dcxs[npont] = dcx1
                                ips.dcys[npont] = dcy1
                                ips.dczs[npont] = dcz1
                                # FIND LAST PLANE BEFORE ZSS(NPONT)
                                idm1 = int(1.0 + trunc(ips.zss[npont] / cion.zstep))
                                if idm1 < 1: idm1 = 1
                                if idm1 > 9: idm1 = 9
                                ips.ipls[npont] = idm1
                                splout.nesst[ips.ipls[npont]] = splout.nesst[ips.ipls[npont]] + 1
                                break L5
                              # PENNING TRANSFER OCCURS BEFORE FIRST SPACE PLANE CLEAR ENTRY
                              npont = npont - 1
                              nclus = nclus - 1
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
                                f5s = f5 * sqrt(e1 / ips.ess[ncltmp])
                                if f5s > 1.0: f5s = 1.0
                                thsec = arcsin(f5s)
                                f5s = sin(thsec)
                                f6s = cos(thsec)
                                if f6 < 0.0: f6s = -f6s
                                phis = phi0 + setp.api
                                if phis > f4: phis = phi0 - f4
                                f8s = sin(phis)
                                f9s = cos(phis)
                                ips.dczs[ncltmp] = f6s
                                ips.dcys[ncltmp] = f8s * f5s
                                ips.dcxs[ncltmp] = f9s * f5s
                                ntmpflg = 0
                              break L190
                            dcz1 = dcz2 * f6 + argz * f5 * f8
                            dcy1 = dcy2 * f6 + f5 / argz * (dcx2 * f9 - dcy2 * dcz2 * f8)
                            dcx1 = dcx2 * f6 - f5 / argz * (dcy2 * f9 + dcx2 * dcz2 * f8)
                            if ntmpflg == 1:
                              # USE FREE KINEMATICS FOR IONISATION SECONDARY ANGLE
                              f5s = f5 * sqrt(e1 / ips.ess[ncltmp])
                              if f5s > 1.0: f5s = 1.0
                              thsec = arcsin(f5s)
                              f5s = sin(thsec)
                              f6s = cos(thsec)
                              if f6 < 0.0: f6s = -f6s
                              phis = phi0 + setp.api
                              if phis > f4: phis = phi0 - f4
                              f8s = sin(phis)
                              f9s = cos(phis)
                              ips.dczs[ncltmp] = dcz2 * f6s + argz * f5s * f8s
                              ips.dcys[ncltmp] = dcy2 * f6s + f5s / argz * (dcx2 * f9s - dcy2 * dcz2 * f8s)
                              ips.dcxs[ncltmp] = dcx2 * f6s - f5s / argz * (dcy2 * f9s + dcx2 * dcz2 * f8s)
                              ntmpflg = 0
                          # STORE DIRECTION COSINES AND ENERGY AFTER N COLLISIONS
                          # FOR LATER REUSE IN PRIMARY GENERATION.
                          i100 = i100 + 1
                          if i100 == 200:
                            dcz100 = dcz1
                            dcx100 = dcx1
                            dcy100 = dcy1
                            e100 = e1
                            i100 = 0
                          if outpt.z > cion.zfinal:
                            # CHECK IF ELECTRON HAS ENOUGH ENERGY TO GO BACK TO FINAL PLANE.
                            epot = setp.efield * (outpt.z - cion.zfinal) * 100.0
                            if e1 < epot:
                              ttrm.ztot = ttrm.ztot + outpt.z
                              ttrm.ttot = ttrm.ttot + outpt.st
                              ttrm.ztots = ttrm.ztots + outpt.z - zstrt
                              ttrm.ttots = ttrm.ttots + outpt.st - tsstrt
                              # NO MORE ELECTRONS IN CASCADE RETURN TO MAIN.
                              if nelec == nclus + 1: break R544
                              break L20
                          # CALCULATE TIME TZSTOP TO ARRIVE AT NEXT Z PLANE IZPLANE.
                          tcalc(outpt.z, dcz1, e1, setp.efield, tzstop, tzstop1, isol, izplane)
                          # CATCH RUNAWAY ELECTRONS AT HIGH FIELD
                          if tzstop == (-99.0):
                            ttrm.ztot = ttrm.ztot + outpt.z
                            ttrm.ttot = ttrm.ttot + outpt.st
                            ttrm.ztots = ttrm.ztots + outpt.z - zstrt
                            ttrm.ttots = ttrm.ttots + outpt.st - tsstrt
                            # NO MORE ELECTRONS IN CASCADE RETURN TO MAIN.
                            if nelec == nclus + 1: break R544
                            break L20
                          if iprint > jprint: break L200
                          break R1
                        block R20x:
                          while true:
                            block R20:
                              # TAKE ELECTRONS FROM STORE
                              outpt.x = ips.xss[npont]
                              outpt.y = ips.yss[npont]
                              outpt.z = ips.zss[npont]
                              outpt.st = ips.tss[npont]
                              e1 = ips.ess[npont]
                              dcx1 = ips.dcxs[npont]
                              dcy1 = ips.dcys[npont]
                              dcz1 = ips.dczs[npont]
                              izplane = ips.ipls[npont]
                              npont = npont - 1
                              zstrt = outpt.z
                              tsstrt = outpt.st
                              if outpt.z > cion.zfinal:
                                # CHECK IF ELECTRON HAS ENOUGH ENERGY TO GO BACK TO FINAL PLANE
                                epot = setp.efield * (outpt.z - cion.zfinal) * 100.0
                                if e1 < epot:
                                  nelec = nelec + 1
                                  isol = 1
                                  ttrm.ztot = ttrm.ztot + outpt.z
                                  ttrm.ttot = ttrm.ttot + outpt.st
                                  ttrm.ztots = ttrm.ztots + outpt.z - zstrt
                                  ttrm.ttots = ttrm.ttots + outpt.st - tsstrt
                                  # NO MORE ELECTRONS IN CASCADE RETURN TO MAIN.
                                  if nelec == nclus + 1: break R544
                                  break R20
                              tcalc(outpt.z, dcz1, e1, setp.efield, tzstop, tzstop1, isol, izplane)
                              if tzstop == (-99.0):
                                # CATCH RUNAWAY ELECTRONS AT HIGH FIELD
                                nelec = nelec + 1
                                isol = 1
                                ttrm.ztot = ttrm.ztot + outpt.z
                                ttrm.ttot = ttrm.ttot + outpt.st
                                ttrm.ztots = ttrm.ztots + outpt.z - zstrt
                                ttrm.ttots = ttrm.ttots + outpt.st - tsstrt
                                # NO MORE ELECTRONS IN CASCADE RETURN TO MAIN.
                                if nelec == nclus + 1: break R544
                                break R20
                              break R20x
                        break R555
                      # INTERMEDIATE PRINTOUT
                      iprint = 0
                      w = ttrm.ztots / ttrm.ttots
                      w = w * 1000000000.0
                      outpt.xid = float64(id)
                      jct = id div 100000
                      # IF(J1.EQ.1) WRITE(6,201)
                      # 201 FORMAT(/,7X,'INTERMEDIATE OUTPUT',/,'    VEL       POS       TIME
                      # /     COUNT')
                      # WRITE(6,202) W,ZTOTS,TTOTS,JCT
                      # 202 FORMAT(1X,F8.3,2(1X,D10.3),1X,I6)
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
    splout.attoion = aneion / anbt
    splout.attater = sqrt(aneion) / aneion
    splout.attioer = sqrt(anbt) / anbt
  else:
    splout.attoion = (-1.0)
    aneion = float64(neion)
    splout.attater = sqrt(aneion) / aneion
  jct = id div 10000
  if j1 == 1:
    fwrite(fmt940, nclus, iter, nelec, neion)
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
  # WRITE(6,837) (IESPECP(J),J=1,100)
  # 837 FORMAT(/,2X,'ENERGY SPECTRUM OF PRIMARIES IN 1EV BINS',/,10(2X,10I
  # /5,/))
  return
  if iter > setp.nmax:
    fwrite(fmt991, iter, setp.nmax, npont, nelec, cion.iprim, nmxadd)
    fstop()
  return
