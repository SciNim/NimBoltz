## SUBROUTINE MIXER of magboltz-11.21.f (line 9955).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import angcut, dmax0, dmin0, gasmix

const
  fmt776 = "(\'  WARNING NEGATIVE COLLISION FEQUENCY =\',D12.3,\' IE =\',I6,\' IF =\',I3,\' IARRY=\',I5,\' EIN=\',F7.4)"
  fmt992 = "(/,/,6X,\'WARNING TOO MANY LEVELS IN CALCULATION. CAN INCREASE THE ARRAY SIZES FROM 960 UP TO 1800 MAXIMUM\',/)"

var q1: Arr2[float64, 6, 4000]
var q2: Arr2[float64, 6, 4000]
var q3: Arr2[float64, 6, 4000]
var q4: Arr2[float64, 6, 4000]
var q5: Arr2[float64, 6, 4000]
var q6: Arr2[float64, 6, 4000]
var qatt: Arr2[float64, 6, 4000]
var peqel1: Arr2[float64, 6, 4000]
var peqel2: Arr2[float64, 6, 4000]
var peqel3: Arr2[float64, 6, 4000]
var peqel4: Arr2[float64, 6, 4000]
var peqel5: Arr2[float64, 6, 4000]
var peqel6: Arr2[float64, 6, 4000]
var peqin1: Arr2[float64, 250, 4000]
var peqin2: Arr2[float64, 250, 4000]
var peqin3: Arr2[float64, 250, 4000]
var peqin4: Arr2[float64, 250, 4000]
var peqin5: Arr2[float64, 250, 4000]
var peqin6: Arr2[float64, 250, 4000]
var qion1: Arr2[float64, 30, 4000]
var qion2: Arr2[float64, 30, 4000]
var qion3: Arr2[float64, 30, 4000]
var qion4: Arr2[float64, 30, 4000]
var qion5: Arr2[float64, 30, 4000]
var qion6: Arr2[float64, 30, 4000]
var peqion1: Arr2[float64, 30, 4000]
var peqion2: Arr2[float64, 30, 4000]
var peqion3: Arr2[float64, 30, 4000]
var peqion4: Arr2[float64, 30, 4000]
var peqion5: Arr2[float64, 30, 4000]
var peqion6: Arr2[float64, 30, 4000]
var qatt1: Arr2[float64, 8, 4000]
var qatt2: Arr2[float64, 8, 4000]
var qatt3: Arr2[float64, 8, 4000]
var qatt4: Arr2[float64, 8, 4000]
var qatt5: Arr2[float64, 8, 4000]
var qatt6: Arr2[float64, 8, 4000]
var qnul1: Arr2[float64, 10, 4000]
var qnul2: Arr2[float64, 10, 4000]
var qnul3: Arr2[float64, 10, 4000]
var qnul4: Arr2[float64, 10, 4000]
var qnul5: Arr2[float64, 10, 4000]
var qnul6: Arr2[float64, 10, 4000]

proc mixer*() =
  var name1, name2, name3, name4, name5, name6, dscrptn: string
  var scrp1, scrp2, scrp3, scrp4, scrp5, scrp6: array[1..300, string]
  var scrpn1, scrpn2, scrpn3, scrpn4, scrpn5, scrpn6: array[1..10, string]
  var e1, e2, e3, e4, e5, e6, eion: array[1..6, float64]
  var ei1, ei2, ei3, ei4, ei5, ei6: array[1..250, float64]
  var qqrot, qdrot: array[1..4000, float64]
  var penfra1, penfra2, penfra3, penfra4, penfra5, penfra6: Arr2[float64, 3, 250]
  var kin1, kin2, kin3, kin4, kin5, kin6: array[1..250, int]
  var kel1, kel2, kel3, kel4, kel5, kel6: array[1..6, int]
  var eion1, eion2, eion3, eion4, eion5, eion6, eb1, eb2: array[1..30, float64]
  var eb3, eb4, eb5, eb6, ec01, ec02, ec03, ec04: array[1..30, float64]
  var ec05, ec06, eg11, eg12, eg13, eg14, eg15, eg16: array[1..30, float64]
  var eg21, eg22, eg23, eg24, eg25, eg26, wk1, wk2: array[1..30, float64]
  var wk3, wk4, wk5, wk6, efl1, efl2, efl3, efl4: array[1..30, float64]
  var efl5, efl6: array[1..30, float64]
  var nc01, nc02, nc03, nc04, nc05, nc06, ng11, ng12: array[1..30, int]
  var ng13, ng14, ng15, ng16, ng21, ng22, ng23, ng24: array[1..30, int]
  var ng25, ng26: array[1..30, int]
  var scln1, scln2, scln3, scln4, scln5, scln6: array[1..10, float64]
  var nion1, nion2, nion3, nion4, nion5, nion6, natt1, natt2: int
  var natt3, natt4, natt5, natt6, nul1, nul2, nul3, nul4: int
  var nul5, nul6, n4000, none, j, i, k, ie: int
  var np, l, kion, jj, `if`, kelsum, jlow, jhi: int
  var katt: int
  var ehalf, aj, virial1, virial2, virial3, virial4, virial5, virial6: float64
  var psct1, angc, psct2, rgas1, rgas2, rgas3, rgas4, rgas5: float64
  var rgas6, bp, f2, elow, ehi: float64
  # SUBROUTINE FILLS ARRAYS OF COLLISION FREQUENCY
  # CAN HAVE A MIXTURE OF UP TO 6 GASES
  anis.niso = 0
  mix3.nin1 = 0
  mix3.nin2 = 0
  mix3.nin3 = 0
  mix3.nin4 = 0
  mix3.nin5 = 0
  mix3.nin6 = 0
  nion1 = 0
  nion2 = 0
  nion3 = 0
  nion4 = 0
  nion5 = 0
  nion6 = 0
  natt1 = 0
  natt2 = 0
  natt3 = 0
  natt4 = 0
  natt5 = 0
  natt6 = 0
  nul1 = 0
  nul2 = 0
  nul3 = 0
  nul4 = 0
  nul5 = 0
  nul6 = 0
  n4000 = 4000
  none = 1
  fdo j, 1, 6:
    names.nameg[j] = fstr("-------------------------", 25)
    kel1[j] = 0
    kel2[j] = 0
    kel3[j] = 0
    kel4[j] = 0
    kel5[j] = 0
    kel6[j] = 0
    fdo i, 1, n4000:
      q1[j, i] = 0.0
      q2[j, i] = 0.0
      q3[j, i] = 0.0
      q4[j, i] = 0.0
      q5[j, i] = 0.0
      q6[j, i] = 0.0
    e1[j] = 0.0
    e2[j] = 0.0
    e3[j] = 0.0
    e4[j] = 0.0
    e5[j] = 0.0
    e6[j] = 0.0
  fdo j, 1, 30:
    eion1[j] = 0.0
    eion2[j] = 0.0
    eion3[j] = 0.0
    eion4[j] = 0.0
    eion5[j] = 0.0
    eion6[j] = 0.0
    eb1[j] = 0.0
    eb2[j] = 0.0
    eb3[j] = 0.0
    eb4[j] = 0.0
    eb5[j] = 0.0
    eb6[j] = 0.0
    ec01[j] = 0.0
    ec02[j] = 0.0
    ec03[j] = 0.0
    ec04[j] = 0.0
    ec05[j] = 0.0
    ec06[j] = 0.0
    eg11[j] = 0.0
    eg12[j] = 0.0
    eg13[j] = 0.0
    eg14[j] = 0.0
    eg15[j] = 0.0
    eg16[j] = 0.0
    eg21[j] = 0.0
    eg22[j] = 0.0
    eg23[j] = 0.0
    eg24[j] = 0.0
    eg25[j] = 0.0
    eg26[j] = 0.0
    wk1[j] = 0.0
    wk2[j] = 0.0
    wk3[j] = 0.0
    wk4[j] = 0.0
    wk5[j] = 0.0
    wk6[j] = 0.0
    efl1[j] = 0.0
    efl2[j] = 0.0
    efl3[j] = 0.0
    efl4[j] = 0.0
    efl5[j] = 0.0
    efl6[j] = 0.0
    nc01[j] = 0
    nc02[j] = 0
    nc03[j] = 0
    nc04[j] = 0
    nc05[j] = 0
    nc06[j] = 0
    ng11[j] = 0
    ng12[j] = 0
    ng13[j] = 0
    ng14[j] = 0
    ng15[j] = 0
    ng16[j] = 0
    ng21[j] = 0
    ng22[j] = 0
    ng23[j] = 0
    ng24[j] = 0
    ng25[j] = 0
    ng26[j] = 0
    fdo i, 1, n4000:
      qion1[j, i] = 0.0
      qion2[j, i] = 0.0
      qion3[j, i] = 0.0
      qion4[j, i] = 0.0
      qion5[j, i] = 0.0
      qion6[j, i] = 0.0
  fdo k, 1, 8:
    fdo i, 1, n4000:
      qatt1[k, i] = 0.0
      qatt2[k, i] = 0.0
      qatt3[k, i] = 0.0
      qatt4[k, i] = 0.0
      qatt5[k, i] = 0.0
      qatt6[k, i] = 0.0
  fdo k, 1, 10:
    fdo i, 1, n4000:
      qnul1[k, i] = 0.0
      qnul2[k, i] = 0.0
      qnul3[k, i] = 0.0
      qnul4[k, i] = 0.0
      qnul5[k, i] = 0.0
      qnul6[k, i] = 0.0
  inpt.estep = inpt.efinal / float64(inpt.nstep)
  ehalf = inpt.estep / 2.0
  mix2.eg[1] = ehalf
  fdo i, 2, 4000:
    aj = float64(i - 1)
    mix2.eg[i] = ehalf + inpt.estep * aj
    # LOGARITHMIC SPACING
    # DLOG_MIN=DLOG10(ESTEP)
    # DLOG_MAX=DLOG10(EFINAL)
    # DLOG_STEP=(DLOG_MAX-DLOG_MIN)/DFLOAT(NSTEP-1)
    # DO 3 I=1,NSTEP
    # AJ=DFLOAT(I-1)
    # E(I)=10.0D0**(DLOG_MIN+AJ*DLOG_STEP)
    mix2.eroot[i] = sqrt(mix2.eg[i])
  mix2.eroot[1] = sqrt(ehalf)
  fdo i, 1, 250:
    kin1[i] = 0
    kin2[i] = 0
    kin3[i] = 0
    kin4[i] = 0
    kin5[i] = 0
    kin6[i] = 0
  fdo i, 1, 960:
    anis.index[i] = 0
  fdo i, 1, 8:
    large.tcfmax[i] = 0.0
  # CALL GAS CROSS-SECTIONS
  gasmix(gasn.ngasn[1], q1, mix1.qin1, mix3.nin1, e1, ei1, name1, virial1, eb1, peqel1, peqin1, penfra1, kel1, kin1, qion1, peqion1, eion1, nion1, qatt1, natt1, qnul1, nul1, scln1, nc01, ec01, wk1, efl1, ng11, eg11, ng21, eg21, scrp1, scrpn1)
  block L10:
    if inpt.ngas == 1: break L10
    gasmix(gasn.ngasn[2], q2, mix1.qin2, mix3.nin2, e2, ei2, name2, virial2, eb2, peqel2, peqin2, penfra2, kel2, kin2, qion2, peqion2, eion2, nion2, qatt2, natt2, qnul2, nul2, scln2, nc02, ec02, wk2, efl2, ng12, eg12, ng22, eg22, scrp2, scrpn2)
    if inpt.ngas == 2: break L10
    gasmix(gasn.ngasn[3], q3, mix1.qin3, mix3.nin3, e3, ei3, name3, virial3, eb3, peqel3, peqin3, penfra3, kel3, kin3, qion3, peqion3, eion3, nion3, qatt3, natt3, qnul3, nul3, scln3, nc03, ec03, wk3, efl3, ng13, eg13, ng23, eg23, scrp3, scrpn3)
    if inpt.ngas == 3: break L10
    gasmix(gasn.ngasn[4], q4, mix1.qin4, mix3.nin4, e4, ei4, name4, virial4, eb4, peqel4, peqin4, penfra4, kel4, kin4, qion4, peqion4, eion4, nion4, qatt4, natt4, qnul4, nul4, scln4, nc04, ec04, wk4, efl4, ng14, eg14, ng24, eg24, scrp4, scrpn4)
    if inpt.ngas == 4: break L10
    gasmix(gasn.ngasn[5], q5, mix1.qin5, mix3.nin5, e5, ei5, name5, virial5, eb5, peqel5, peqin5, penfra5, kel5, kin5, qion5, peqion5, eion5, nion5, qatt5, natt5, qnul5, nul5, scln5, nc05, ec05, wk5, efl5, ng15, eg15, ng25, eg25, scrp5, scrpn5)
    if inpt.ngas == 5: break L10
    gasmix(gasn.ngasn[6], q6, mix1.qin6, mix3.nin6, e6, ei6, name6, virial6, eb6, peqel6, peqin6, penfra6, kel6, kin6, qion6, peqion6, eion6, nion6, qatt6, natt6, qnul6, nul6, scln6, nc06, ec06, wk6, efl6, ng16, eg16, ng26, eg26, scrp6, scrpn6)
  # CORRECTION OF NUMBER DENSITY DUE TO VIRIAL COEFFICIENT
  # CAN BE PROGRAMMED HERE NOT YET IMPLEMENTED.
  # CALCULATION OF COLLISION FREQUENCIES FOR AN ARRAY OF
  # ELECTRON ENERGIES IN THE RANGE ZERO TO EFINAL
  # L=5*N-4    ELASTIC NTH GAS
  # L=5*N-3    IONISATION NTH GAS
  # L=5*N-2    ATTACHMENT NTH GAS
  # L=5*N-1    INELASTIC NTH GAS
  # L=5*N      SUPERELASTIC NTH GAS
  fdo ie, 1, 4000:
    fred.fcion[ie] = 0.0
    fred.fcatt[ie] = 0.0
    np = 1
    large.cf[ie, np] = q1[2, ie] * mratio.van1
    anis.psct[ie, np] = 0.5
    anis.angct[ie, np] = 1.0
    anis.index[np] = 0
    # ELASTIC ANG
    if kel1[2] == 1:
      psct1 = peqel1[2, ie]
      angcut(psct1, angc, psct2)
      anis.angct[ie, np] = angc
      anis.psct[ie, np] = psct2
      anis.index[np] = 1
    if kel1[2] == 2:
      anis.psct[ie, np] = peqel1[2, ie]
      anis.index[np] = 2
    if not (ie > 1):
      rgas1 = 1.0 + e1[2] / 2.0
      large.rgas[np] = rgas1
      large.ein[np] = 0.0
      large.ipn[np] = 0
      l = 1
      large.iarry[np] = l
      scrip.dscrpt[np] = fstr(scrp1[2], 50)
      names.nameg[1] = fstr(name1, 25)
      large.penfra[1, np] = 0.0
      large.penfra[2, np] = 0.0
      large.penfra[3, np] = 0.0
    block L30:
      if inpt.efinal < e1[3]: break L30
      if not (nion1 > 1):
        np = np + 1
        large.cf[ie, np] = q1[3, ie] * mratio.van1
        fred.fcion[ie] = fred.fcion[ie] + large.cf[ie, np]
        anis.psct[ie, np] = 0.5
        anis.angct[ie, np] = 1.0
        anis.index[np] = 0
        if kel1[3] == 1:
          psct1 = peqel1[3, ie]
          angcut(psct1, angc, psct2)
          anis.angct[ie, np] = angc
          anis.psct[ie, np] = psct2
          anis.index[np] = 1
        if kel1[3] == 2:
          anis.psct[ie, np] = peqel1[3, ie]
          anis.index[np] = 2
        if ie > 1: break L30
        large.rgas[np] = rgas1
        large.ein[np] = e1[3] / rgas1
        large.wpl[np] = eb1[1]
        ionfl1.nc0[np] = nc01[1]
        ionfl1.ec0[np] = ec01[1]
        ionfl1.ng1[np] = ng11[1]
        ionfl1.eg1[np] = eg11[1]
        ionfl1.ng2[np] = ng21[1]
        ionfl1.eg2[np] = eg21[1]
        ionfl1.wklm[np] = wk1[1]
        ionfl1.efl[np] = efl1[1]
        large.ipn[np] = 1
        l = 2
        large.iarry[np] = l
        scrip.dscrpt[np] = fstr(scrp1[3], 50)
        large.penfra[1, np] = 0.0
        large.penfra[2, np] = 0.0
        large.penfra[3, np] = 0.0
        break L30
      fdo kion, 1, nion1:
        np = np + 1
        large.cf[ie, np] = qion1[kion, ie] * mratio.van1
        fred.fcion[ie] = fred.fcion[ie] + large.cf[ie, np]
        anis.psct[ie, np] = 0.5
        anis.angct[ie, np] = 1.0
        anis.index[np] = 0
        if kel1[3] == 1:
          psct1 = peqion1[kion, ie]
          angcut(psct1, angc, psct2)
          anis.angct[ie, np] = angc
          anis.psct[ie, np] = psct2
          anis.index[np] = 1
        if kel1[3] == 2:
          anis.psct[ie, np] = peqion1[kion, ie]
          anis.index[np] = 2
        if not (ie > 1):
          large.rgas[np] = rgas1
          large.ein[np] = eion1[kion] / rgas1
          large.wpl[np] = eb1[kion]
          ionfl1.nc0[np] = nc01[kion]
          ionfl1.ec0[np] = ec01[kion]
          ionfl1.ng1[np] = ng11[kion]
          ionfl1.eg1[np] = eg11[kion]
          ionfl1.ng2[np] = ng21[kion]
          ionfl1.eg2[np] = eg21[kion]
          ionfl1.wklm[np] = wk1[kion]
          ionfl1.efl[np] = efl1[kion]
          large.ipn[np] = 1
          l = 2
          large.iarry[np] = l
          scrip.dscrpt[np] = fstr(scrp1[2 + kion], 50)
          large.penfra[1, np] = 0.0
          large.penfra[2, np] = 0.0
          large.penfra[3, np] = 0.0
    block L40:
      if inpt.efinal < e1[4]: break L40
      if not (natt1 > 1):
        np = np + 1
        large.cf[ie, np] = q1[4, ie] * mratio.van1
        fred.fcatt[ie] = fred.fcatt[ie] + large.cf[ie, np]
        anis.psct[ie, np] = 0.5
        anis.angct[ie, np] = 1.0
        if ie > 1: break L40
        anis.index[np] = 0
        large.rgas[np] = rgas1
        large.ein[np] = 0.0
        large.ipn[np] = (-1)
        l = 3
        large.iarry[np] = l
        scrip.dscrpt[np] = fstr(scrp1[3 + nion1], 50)
        large.penfra[1, np] = 0.0
        large.penfra[2, np] = 0.0
        large.penfra[3, np] = 0.0
        break L40
      fdo jj, 1, natt1:
        np = np + 1
        large.cf[ie, np] = qatt1[jj, ie] * mratio.van1
        fred.fcatt[ie] = fred.fcatt[ie] + large.cf[ie, np]
        anis.psct[ie, np] = 0.5
        anis.angct[ie, np] = 1.0
        if not (ie > 1):
          scrip.dscrpt[np] = fstr(scrp1[2 + nion1 + jj], 50)
          large.rgas[np] = rgas1
          large.ein[np] = 0.0
          large.ipn[np] = (-1)
          anis.index[np] = 0
          l = 3
          large.iarry[np] = l
          large.penfra[1, np] = 0.0
          large.penfra[2, np] = 0.0
          large.penfra[3, np] = 0.0
    if not (mix3.nin1 == 0):
      fdo j, 1, mix3.nin1:
        np = np + 1
        large.cf[ie, np] = mix1.qin1[j, ie] * mratio.van1
        anis.psct[ie, np] = 0.5
        anis.angct[ie, np] = 1.0
        anis.index[np] = 0
        if kin1[j] == 1:
          psct1 = peqin1[j, ie]
          angcut(psct1, angc, psct2)
          anis.angct[ie, np] = angc
          anis.psct[ie, np] = psct2
          anis.index[np] = 1
        if kin1[j] == 2:
          anis.psct[ie, np] = peqin1[j, ie]
          anis.index[np] = 2
        if not (ie > 1):
          large.rgas[np] = rgas1
          large.ein[np] = ei1[j] / rgas1
          l = 4
          if ei1[j] < 0.0: l = 5
          large.ipn[np] = 0
          large.iarry[np] = l
          scrip.dscrpt[np] = fstr(scrp1[4 + nion1 + natt1 + j], 50)
          large.penfra[1, np] = penfra1[1, j]
          large.penfra[2, np] = penfra1[2, j] * 1e-16 / sqrt(3.0)
          large.penfra[3, np] = penfra1[3, j]
    block L600:
      if inpt.ngas == 1: break L600
      np = np + 1
      large.cf[ie, np] = q2[2, ie] * mratio.van2
      anis.psct[ie, np] = 0.5
      anis.angct[ie, np] = 1.0
      anis.index[np] = 0
      if kel2[2] == 1:
        psct1 = peqel2[2, ie]
        angcut(psct1, angc, psct2)
        anis.angct[ie, np] = angc
        anis.psct[ie, np] = psct2
        anis.index[np] = 1
      if kel2[2] == 2:
        anis.psct[ie, np] = peqel2[2, ie]
        anis.index[np] = 2
      if not (ie > 1):
        rgas2 = 1.0 + e2[2] / 2.0
        large.rgas[np] = rgas2
        large.ein[np] = 0.0
        large.ipn[np] = 0
        l = 6
        large.iarry[np] = l
        scrip.dscrpt[np] = fstr(scrp2[2], 50)
        names.nameg[2] = fstr(name2, 25)
        large.penfra[1, np] = 0.0
        large.penfra[2, np] = 0.0
        large.penfra[3, np] = 0.0
      block L130:
        if inpt.efinal < e2[3]: break L130
        if not (nion2 > 1):
          np = np + 1
          large.cf[ie, np] = q2[3, ie] * mratio.van2
          fred.fcion[ie] = fred.fcion[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kel2[3] == 1:
            psct1 = peqel2[3, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kel2[3] == 2:
            anis.psct[ie, np] = peqel2[3, ie]
            anis.index[np] = 2
          if ie > 1: break L130
          large.rgas[np] = rgas2
          large.ein[np] = e2[3] / rgas2
          large.wpl[np] = eb2[1]
          ionfl1.nc0[np] = nc02[1]
          ionfl1.ec0[np] = ec02[1]
          ionfl1.ng1[np] = ng12[1]
          ionfl1.eg1[np] = eg12[1]
          ionfl1.ng2[np] = ng22[1]
          ionfl1.eg2[np] = eg22[1]
          ionfl1.wklm[np] = wk2[1]
          ionfl1.efl[np] = efl2[1]
          large.ipn[np] = 1
          l = 7
          large.iarry[np] = l
          scrip.dscrpt[np] = fstr(scrp2[3], 50)
          large.penfra[1, np] = 0.0
          large.penfra[2, np] = 0.0
          large.penfra[3, np] = 0.0
          break L130
        fdo kion, 1, nion2:
          np = np + 1
          large.cf[ie, np] = qion2[kion, ie] * mratio.van2
          fred.fcion[ie] = fred.fcion[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kel2[3] == 1:
            psct1 = peqion2[kion, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kel2[3] == 2:
            anis.psct[ie, np] = peqion2[kion, ie]
            anis.index[np] = 2
          if not (ie > 1):
            large.rgas[np] = rgas2
            large.ein[np] = eion2[kion] / rgas2
            large.wpl[np] = eb2[kion]
            ionfl1.nc0[np] = nc02[kion]
            ionfl1.ec0[np] = ec02[kion]
            ionfl1.ng1[np] = ng12[kion]
            ionfl1.eg1[np] = eg12[kion]
            ionfl1.ng2[np] = ng22[kion]
            ionfl1.eg2[np] = eg22[kion]
            ionfl1.wklm[np] = wk2[kion]
            ionfl1.efl[np] = efl1[kion]
            large.ipn[np] = 1
            l = 7
            large.iarry[np] = l
            scrip.dscrpt[np] = fstr(scrp2[2 + kion], 50)
            large.penfra[1, np] = 0.0
            large.penfra[2, np] = 0.0
            large.penfra[3, np] = 0.0
      block L140:
        if inpt.efinal < e2[4]: break L140
        if not (natt2 > 1):
          np = np + 1
          large.cf[ie, np] = q2[4, ie] * mratio.van2
          fred.fcatt[ie] = fred.fcatt[ie] + large.cf[ie, np]
          if ie > 1: break L140
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          large.rgas[np] = rgas2
          large.ein[np] = 0.0
          large.ipn[np] = (-1)
          l = 8
          large.iarry[np] = l
          scrip.dscrpt[np] = fstr(scrp2[3 + nion2], 50)
          large.penfra[1, np] = 0.0
          large.penfra[2, np] = 0.0
          large.penfra[3, np] = 0.0
          break L140
        fdo jj, 1, natt2:
          np = np + 1
          large.cf[ie, np] = qatt2[jj, ie] * mratio.van2
          fred.fcatt[ie] = fred.fcatt[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          if not (ie > 1):
            scrip.dscrpt[np] = fstr(scrp2[2 + nion2 + jj], 50)
            anis.index[np] = 0
            large.rgas[np] = rgas2
            large.ein[np] = 0.0
            large.ipn[np] = (-1)
            l = 8
            large.iarry[np] = l
            large.penfra[1, np] = 0.0
            large.penfra[2, np] = 0.0
            large.penfra[3, np] = 0.0
      if not (mix3.nin2 == 0):
        fdo j, 1, mix3.nin2:
          np = np + 1
          large.cf[ie, np] = mix1.qin2[j, ie] * mratio.van2
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kin2[j] == 1:
            psct1 = peqin2[j, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
            if kin2[j] == 2:
              anis.psct[ie, np] = peqin2[j, ie]
              anis.index[np] = 2
          if not (ie > 1):
            large.rgas[np] = rgas2
            large.ein[np] = ei2[j] / rgas2
            l = 9
            if ei2[j] < 0.0: l = 10
            large.ipn[np] = 0
            large.iarry[np] = l
            scrip.dscrpt[np] = fstr(scrp2[4 + nion2 + natt2 + j], 50)
            large.penfra[1, np] = penfra2[1, j]
            large.penfra[2, np] = penfra2[2, j] * 0.000001 / sqrt(3.0)
            large.penfra[3, np] = penfra2[3, j]
      if inpt.ngas == 2: break L600
      np = np + 1
      large.cf[ie, np] = q3[2, ie] * mratio.van3
      anis.psct[ie, np] = 0.5
      anis.angct[ie, np] = 1.0
      anis.index[np] = 0
      if kel3[2] == 1:
        psct1 = peqel3[2, ie]
        angcut(psct1, angc, psct2)
        anis.angct[ie, np] = angc
        anis.psct[ie, np] = psct2
        anis.index[np] = 1
      if kel3[2] == 2:
        anis.psct[ie, np] = peqel3[2, ie]
        anis.index[np] = 2
      if not (ie > 1):
        rgas3 = 1.0 + e3[2] / 2.0
        large.rgas[np] = rgas3
        large.ein[np] = 0.0
        large.ipn[np] = 0
        l = 11
        large.iarry[np] = l
        scrip.dscrpt[np] = fstr(scrp3[2], 50)
        names.nameg[3] = fstr(name3, 25)
        large.penfra[1, np] = 0.0
        large.penfra[2, np] = 0.0
        large.penfra[3, np] = 0.0
      block L230:
        if inpt.efinal < e3[3]: break L230
        if not (nion3 > 1):
          np = np + 1
          large.cf[ie, np] = q3[3, ie] * mratio.van3
          fred.fcion[ie] = fred.fcion[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kel3[3] == 1:
            psct1 = peqel3[3, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kel3[3] == 2:
            anis.psct[ie, np] = peqel3[3, ie]
            anis.index[np] = 2
          if ie > 1: break L230
          large.rgas[np] = rgas3
          large.ein[np] = e3[3] / rgas3
          large.wpl[np] = eb3[1]
          ionfl1.nc0[np] = nc03[1]
          ionfl1.ec0[np] = ec03[1]
          ionfl1.ng1[np] = ng13[1]
          ionfl1.eg1[np] = eg13[1]
          ionfl1.ng2[np] = ng23[1]
          ionfl1.eg2[np] = eg23[1]
          ionfl1.wklm[np] = wk3[1]
          ionfl1.efl[np] = efl3[1]
          large.ipn[np] = 1
          l = 12
          large.iarry[np] = l
          scrip.dscrpt[np] = fstr(scrp3[3], 50)
          large.penfra[1, np] = 0.0
          large.penfra[2, np] = 0.0
          large.penfra[3, np] = 0.0
          break L230
        fdo kion, 1, nion3:
          np = np + 1
          large.cf[ie, np] = qion3[kion, ie] * mratio.van3
          fred.fcion[ie] = fred.fcion[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kel3[3] == 1:
            psct1 = peqion3[kion, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kel3[3] == 2:
            anis.psct[ie, np] = peqion3[kion, ie]
            anis.index[np] = 2
          if not (ie > 1):
            large.rgas[np] = rgas3
            large.ein[np] = eion3[kion] / rgas3
            large.wpl[np] = eb3[kion]
            ionfl1.nc0[np] = nc03[kion]
            ionfl1.ec0[np] = ec03[kion]
            ionfl1.ng1[np] = ng13[kion]
            ionfl1.eg1[np] = eg13[kion]
            ionfl1.ng2[np] = ng23[kion]
            ionfl1.eg2[np] = eg23[kion]
            ionfl1.wklm[np] = wk3[kion]
            ionfl1.efl[np] = efl3[kion]
            large.ipn[np] = 1
            l = 12
            large.iarry[np] = l
            scrip.dscrpt[np] = fstr(scrp3[2 + kion], 50)
            large.penfra[1, np] = 0.0
            large.penfra[2, np] = 0.0
            large.penfra[3, np] = 0.0
      block L240:
        if inpt.efinal < e3[4]: break L240
        if not (natt3 > 1):
          np = np + 1
          large.cf[ie, np] = q3[4, ie] * mratio.van3
          fred.fcatt[ie] = fred.fcatt[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          if ie > 1: break L240
          anis.index[np] = 0
          large.rgas[np] = rgas3
          large.ein[np] = 0.0
          large.ipn[np] = (-1)
          l = 13
          large.iarry[np] = l
          scrip.dscrpt[np] = fstr(scrp3[3 + nion3], 50)
          large.penfra[1, np] = 0.0
          large.penfra[2, np] = 0.0
          large.penfra[3, np] = 0.0
          break L240
        fdo jj, 1, natt3:
          np = np + 1
          large.cf[ie, np] = qatt3[jj, ie] * mratio.van3
          fred.fcatt[ie] = fred.fcatt[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          if not (ie > 1):
            scrip.dscrpt[np] = fstr(scrp3[2 + nion3 + jj], 50)
            anis.index[np] = 0
            large.rgas[np] = rgas3
            large.ein[np] = 0.0
            large.ipn[np] = (-1)
            l = 13
            large.iarry[np] = l
            large.penfra[1, np] = 0.0
            large.penfra[2, np] = 0.0
            large.penfra[3, np] = 0.0
      if not (mix3.nin3 == 0):
        fdo j, 1, mix3.nin3:
          np = np + 1
          large.cf[ie, np] = mix1.qin3[j, ie] * mratio.van3
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kin3[j] == 1:
            psct1 = peqin3[j, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kin3[j] == 2:
            anis.psct[ie, np] = peqin3[j, ie]
            anis.index[np] = 2
          if not (ie > 1):
            large.rgas[np] = rgas3
            large.ein[np] = ei3[j] / rgas3
            l = 14
            if ei3[j] < 0.0: l = 15
            large.ipn[np] = 0
            large.iarry[np] = l
            scrip.dscrpt[np] = fstr(scrp3[4 + nion3 + natt3 + j], 50)
            large.penfra[1, np] = penfra3[1, j]
            large.penfra[2, np] = penfra3[2, j] * 0.000001 / sqrt(3.0)
            large.penfra[3, np] = penfra3[3, j]
      if inpt.ngas == 3: break L600
      np = np + 1
      large.cf[ie, np] = q4[2, ie] * mratio.van4
      anis.psct[ie, np] = 0.5
      anis.angct[ie, np] = 1.0
      anis.index[np] = 0
      if kel4[2] == 1:
        psct1 = peqel4[2, ie]
        angcut(psct1, angc, psct2)
        anis.angct[ie, np] = angc
        anis.psct[ie, np] = psct2
        anis.index[np] = 1
      if kel4[2] == 2:
        anis.psct[ie, np] = peqel4[2, ie]
        anis.index[np] = 2
      if not (ie > 1):
        rgas4 = 1.0 + e4[2] / 2.0
        large.rgas[np] = rgas4
        large.ein[np] = 0.0
        large.ipn[np] = 0
        l = 16
        large.iarry[np] = l
        scrip.dscrpt[np] = fstr(scrp4[2], 50)
        names.nameg[4] = fstr(name4, 25)
        large.penfra[1, np] = 0.0
        large.penfra[2, np] = 0.0
        large.penfra[3, np] = 0.0
      block L330:
        if inpt.efinal < e4[3]: break L330
        if not (nion4 > 1):
          np = np + 1
          large.cf[ie, np] = q4[3, ie] * mratio.van4
          fred.fcion[ie] = fred.fcion[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kel4[3] == 1:
            psct1 = peqel4[3, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kel4[3] == 2:
            anis.psct[ie, np] = peqel4[3, ie]
            anis.index[np] = 2
          if ie > 1: break L330
          large.rgas[np] = rgas4
          large.ein[np] = e4[3] / rgas4
          large.wpl[np] = eb4[1]
          ionfl1.nc0[np] = nc04[1]
          ionfl1.ec0[np] = ec04[1]
          ionfl1.ng1[np] = ng14[1]
          ionfl1.eg1[np] = eg14[1]
          ionfl1.ng2[np] = ng24[1]
          ionfl1.eg2[np] = eg24[1]
          ionfl1.wklm[np] = wk4[1]
          ionfl1.efl[np] = efl4[1]
          large.ipn[np] = 1
          l = 17
          large.iarry[np] = l
          scrip.dscrpt[np] = fstr(scrp4[3], 50)
          large.penfra[1, np] = 0.0
          large.penfra[2, np] = 0.0
          large.penfra[3, np] = 0.0
          break L330
        fdo kion, 1, nion4:
          np = np + 1
          large.cf[ie, np] = qion4[kion, ie] * mratio.van4
          fred.fcion[ie] = fred.fcion[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kel4[3] == 1:
            psct1 = peqion4[kion, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kel4[3] == 2:
            anis.psct[ie, np] = peqion4[kion, ie]
            anis.index[np] = 2
          if not (ie > 1):
            large.rgas[np] = rgas4
            large.ein[np] = eion4[kion] / rgas4
            large.wpl[np] = eb4[kion]
            ionfl1.nc0[np] = nc04[kion]
            ionfl1.ec0[np] = ec04[kion]
            ionfl1.ng1[np] = ng14[kion]
            ionfl1.eg1[np] = eg14[kion]
            ionfl1.ng2[np] = ng24[kion]
            ionfl1.eg2[np] = eg24[kion]
            ionfl1.wklm[np] = wk4[kion]
            ionfl1.efl[np] = efl4[kion]
            large.ipn[np] = 1
            l = 17
            large.iarry[np] = l
            scrip.dscrpt[np] = fstr(scrp4[2 + kion], 50)
            large.penfra[1, np] = 0.0
            large.penfra[2, np] = 0.0
            large.penfra[3, np] = 0.0
      block L340:
        if inpt.efinal < e4[4]: break L340
        if not (natt4 > 1):
          np = np + 1
          large.cf[ie, np] = q4[4, ie] * mratio.van4
          fred.fcatt[ie] = fred.fcatt[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          if ie > 1: break L340
          anis.index[np] = 0
          large.rgas[np] = rgas4
          large.ein[np] = 0.0
          large.ipn[np] = (-1)
          l = 18
          large.iarry[np] = l
          scrip.dscrpt[np] = fstr(scrp4[3 + nion4], 50)
          large.penfra[1, np] = 0.0
          large.penfra[2, np] = 0.0
          large.penfra[3, np] = 0.0
          break L340
        fdo jj, 1, natt4:
          np = np + 1
          large.cf[ie, np] = qatt4[jj, ie] * mratio.van4
          fred.fcatt[ie] = fred.fcatt[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          if not (ie > 1):
            scrip.dscrpt[np] = fstr(scrp4[2 + nion4 + jj], 50)
            anis.index[np] = 0
            large.rgas[np] = rgas4
            large.ein[np] = 0.0
            large.ipn[np] = (-1)
            l = 18
            large.iarry[np] = l
            large.penfra[1, np] = 0.0
            large.penfra[2, np] = 0.0
            large.penfra[3, np] = 0.0
      if not (mix3.nin4 == 0):
        fdo j, 1, mix3.nin4:
          np = np + 1
          large.cf[ie, np] = mix1.qin4[j, ie] * mratio.van4
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kin4[j] == 1:
            psct1 = peqin4[j, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kin4[j] == 2:
            anis.psct[ie, np] = peqin4[j, ie]
            anis.index[np] = 2
          if not (ie > 1):
            large.rgas[np] = rgas4
            large.ein[np] = ei4[j] / rgas4
            l = 19
            if ei4[j] < 0.0: l = 20
            large.ipn[np] = 0
            large.iarry[np] = l
            scrip.dscrpt[np] = fstr(scrp4[4 + nion4 + natt4 + j], 50)
            large.penfra[1, np] = penfra4[1, j]
            large.penfra[2, np] = penfra4[2, j] * 0.000001 / sqrt(3.0)
            large.penfra[3, np] = penfra4[3, j]
      if inpt.ngas == 4: break L600
      np = np + 1
      large.cf[ie, np] = q5[2, ie] * mratio.van5
      anis.psct[ie, np] = 0.5
      anis.angct[ie, np] = 1.0
      anis.index[np] = 0
      if kel5[2] == 1:
        psct1 = peqel5[2, ie]
        angcut(psct1, angc, psct2)
        anis.angct[ie, np] = angc
        anis.psct[ie, np] = psct2
        anis.index[np] = 1
      if kel5[2] == 2:
        anis.psct[ie, np] = peqel5[2, ie]
        anis.index[np] = 2
      if not (ie > 1):
        rgas5 = 1.0 + e5[2] / 2.0
        large.rgas[np] = rgas5
        large.ein[np] = 0.0
        large.ipn[np] = 0
        l = 21
        large.iarry[np] = l
        scrip.dscrpt[np] = fstr(scrp5[2], 50)
        names.nameg[5] = fstr(name5, 25)
        large.penfra[1, np] = 0.0
        large.penfra[2, np] = 0.0
        large.penfra[3, np] = 0.0
      block L430:
        if inpt.efinal < e5[3]: break L430
        if not (nion5 > 1):
          np = np + 1
          large.cf[ie, np] = q5[3, ie] * mratio.van5
          fred.fcion[ie] = fred.fcion[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kel5[3] == 1:
            psct1 = peqel5[3, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kel5[3] == 2:
            anis.psct[ie, np] = peqel5[3, ie]
            anis.index[np] = 2
          if ie > 1: break L430
          large.rgas[np] = rgas5
          large.ein[np] = e5[3] / rgas5
          large.wpl[np] = eb5[1]
          ionfl1.nc0[np] = nc05[1]
          ionfl1.ec0[np] = ec05[1]
          ionfl1.ng1[np] = ng15[1]
          ionfl1.eg1[np] = eg15[1]
          ionfl1.ng2[np] = ng25[1]
          ionfl1.eg2[np] = eg25[1]
          ionfl1.wklm[np] = wk5[1]
          ionfl1.efl[np] = efl5[1]
          large.ipn[np] = 1
          l = 22
          large.iarry[np] = l
          scrip.dscrpt[np] = fstr(scrp5[3], 50)
          large.penfra[1, np] = 0.0
          large.penfra[2, np] = 0.0
          large.penfra[3, np] = 0.0
          break L430
        fdo kion, 1, nion5:
          np = np + 1
          large.cf[ie, np] = qion5[kion, ie] * mratio.van5
          fred.fcion[ie] = fred.fcion[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kel5[3] == 1:
            psct1 = peqion5[kion, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kel5[3] == 2:
            anis.psct[ie, np] = peqion5[kion, ie]
            anis.index[np] = 2
          if not (ie > 1):
            large.rgas[np] = rgas5
            large.ein[np] = eion5[kion] / rgas5
            large.wpl[np] = eb5[kion]
            ionfl1.nc0[np] = nc05[kion]
            ionfl1.ec0[np] = ec05[kion]
            ionfl1.ng1[np] = ng15[kion]
            ionfl1.eg1[np] = eg15[kion]
            ionfl1.ng2[np] = ng25[kion]
            ionfl1.eg2[np] = eg25[kion]
            ionfl1.wklm[np] = wk5[kion]
            ionfl1.efl[np] = efl5[kion]
            large.ipn[np] = 1
            l = 22
            large.iarry[np] = l
            scrip.dscrpt[np] = fstr(scrp5[2 + kion], 50)
            large.penfra[1, np] = 0.0
            large.penfra[2, np] = 0.0
            large.penfra[3, np] = 0.0
      block L440:
        if inpt.efinal < e5[4]: break L440
        if not (natt5 > 1):
          np = np + 1
          large.cf[ie, np] = q5[4, ie] * mratio.van5
          fred.fcatt[ie] = fred.fcatt[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          if ie > 1: break L440
          anis.index[np] = 0
          large.rgas[np] = rgas5
          large.ein[np] = 0.0
          large.ipn[np] = (-1)
          l = 23
          large.iarry[np] = l
          scrip.dscrpt[np] = fstr(scrp5[3 + nion5], 50)
          large.penfra[1, np] = 0.0
          large.penfra[2, np] = 0.0
          large.penfra[3, np] = 0.0
          break L440
        fdo jj, 1, natt5:
          np = np + 1
          large.cf[ie, np] = qatt5[jj, ie] * mratio.van5
          fred.fcatt[ie] = fred.fcatt[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          if not (ie > 1):
            scrip.dscrpt[np] = fstr(scrp5[2 + nion5 + jj], 50)
            anis.index[np] = 0
            large.rgas[np] = rgas5
            large.ein[np] = 0.0
            large.ipn[np] = (-1)
            l = 23
            large.iarry[np] = l
            large.penfra[1, np] = 0.0
            large.penfra[2, np] = 0.0
            large.penfra[3, np] = 0.0
      if not (mix3.nin5 == 0):
        fdo j, 1, mix3.nin5:
          np = np + 1
          large.cf[ie, np] = mix1.qin5[j, ie] * mratio.van5
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kin5[j] == 1:
            psct1 = peqin5[j, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kin5[j] == 2:
            anis.psct[ie, np] = peqin5[j, ie]
            anis.index[np] = 2
          if not (ie > 1):
            large.rgas[np] = rgas5
            large.ein[np] = ei5[j] / rgas5
            l = 24
            if ei5[j] < 0.0: l = 25
            large.ipn[np] = 0
            large.iarry[np] = l
            scrip.dscrpt[np] = fstr(scrp5[4 + nion5 + natt5 + j], 50)
            large.penfra[1, np] = penfra5[1, j]
            large.penfra[2, np] = penfra5[2, j] * 0.000001 / sqrt(3.0)
            large.penfra[3, np] = penfra5[3, j]
      if inpt.ngas == 5: break L600
      np = np + 1
      large.cf[ie, np] = q6[2, ie] * mratio.van6
      anis.psct[ie, np] = 0.5
      anis.angct[ie, np] = 1.0
      anis.index[np] = 0
      if kel6[2] == 1:
        psct1 = peqel6[2, ie]
        angcut(psct1, angc, psct2)
        anis.angct[ie, np] = angc
        anis.psct[ie, np] = psct2
        anis.index[np] = 1
      if kel6[2] == 2:
        anis.psct[ie, np] = peqel6[2, ie]
        anis.index[np] = 2
      if not (ie > 1):
        rgas6 = 1.0 + e6[2] / 2.0
        large.rgas[np] = rgas6
        large.ein[np] = 0.0
        large.ipn[np] = 0
        l = 26
        large.iarry[np] = l
        scrip.dscrpt[np] = fstr(scrp6[2], 50)
        names.nameg[6] = fstr(name6, 25)
        large.penfra[1, np] = 0.0
        large.penfra[2, np] = 0.0
        large.penfra[3, np] = 0.0
      block L530:
        if inpt.efinal < e6[3]: break L530
        if not (nion6 > 1):
          np = np + 1
          large.cf[ie, np] = q6[3, ie] * mratio.van6
          fred.fcion[ie] = fred.fcion[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kel6[3] == 1:
            psct1 = peqel6[3, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kel6[3] == 2:
            anis.psct[ie, np] = peqel6[3, ie]
            anis.index[np] = 2
          if ie > 1: break L530
          large.rgas[np] = rgas6
          large.ein[np] = e6[3] / rgas6
          large.wpl[np] = eb6[1]
          ionfl1.nc0[np] = nc06[1]
          ionfl1.ec0[np] = ec06[1]
          ionfl1.ng1[np] = ng16[1]
          ionfl1.eg1[np] = eg16[1]
          ionfl1.ng2[np] = ng26[1]
          ionfl1.eg2[np] = eg26[1]
          ionfl1.wklm[np] = wk6[1]
          ionfl1.efl[np] = efl6[1]
          large.ipn[np] = 1
          l = 27
          large.iarry[np] = l
          scrip.dscrpt[np] = fstr(scrp6[3], 50)
          large.penfra[1, np] = 0.0
          large.penfra[2, np] = 0.0
          large.penfra[3, np] = 0.0
          break L530
        fdo kion, 1, nion6:
          np = np + 1
          large.cf[ie, np] = qion6[kion, ie] * mratio.van6
          fred.fcion[ie] = fred.fcion[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kel6[3] == 1:
            psct1 = peqion6[kion, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kel6[3] == 2:
            anis.psct[ie, np] = peqion6[kion, ie]
            anis.index[np] = 2
          if not (ie > 1):
            large.rgas[np] = rgas6
            large.ein[np] = eion6[kion] / rgas6
            large.wpl[np] = eb6[kion]
            ionfl1.nc0[np] = nc06[kion]
            ionfl1.ec0[np] = ec06[kion]
            ionfl1.ng1[np] = ng16[kion]
            ionfl1.eg1[np] = eg16[kion]
            ionfl1.ng2[np] = ng26[kion]
            ionfl1.eg2[np] = eg26[kion]
            ionfl1.wklm[np] = wk6[kion]
            ionfl1.efl[np] = efl6[kion]
            large.ipn[np] = 1
            l = 27
            large.iarry[np] = l
            scrip.dscrpt[np] = fstr(scrp6[2 + kion], 50)
            large.penfra[1, np] = 0.0
            large.penfra[2, np] = 0.0
            large.penfra[3, np] = 0.0
      block L540:
        if inpt.efinal < e6[4]: break L540
        if not (natt6 > 1):
          np = np + 1
          large.cf[ie, np] = q6[4, ie] * mratio.van6
          fred.fcatt[ie] = fred.fcatt[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          if ie > 1: break L540
          anis.index[np] = 0
          large.rgas[np] = rgas6
          large.ein[np] = 0.0
          large.ipn[np] = (-1)
          l = 28
          large.iarry[np] = l
          scrip.dscrpt[np] = fstr(scrp6[3 + nion6], 50)
          large.penfra[1, np] = 0.0
          large.penfra[2, np] = 0.0
          large.penfra[3, np] = 0.0
          break L540
        fdo jj, 1, natt6:
          np = np + 1
          large.cf[ie, np] = qatt6[jj, ie] * mratio.van6
          fred.fcatt[ie] = fred.fcatt[ie] + large.cf[ie, np]
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          if not (ie > 1):
            scrip.dscrpt[np] = fstr(scrp6[2 + nion6 + jj], 50)
            anis.index[np] = 0
            large.rgas[np] = rgas6
            large.ein[np] = 0.0
            large.ipn[np] = (-1)
            l = 28
            large.iarry[np] = l
            large.penfra[1, np] = 0.0
            large.penfra[2, np] = 0.0
            large.penfra[3, np] = 0.0
      if not (mix3.nin6 == 0):
        fdo j, 1, mix3.nin6:
          np = np + 1
          large.cf[ie, np] = mix1.qin6[j, ie] * mratio.van6
          anis.psct[ie, np] = 0.5
          anis.angct[ie, np] = 1.0
          anis.index[np] = 0
          if kin6[j] == 1:
            psct1 = peqin6[j, ie]
            angcut(psct1, angc, psct2)
            anis.angct[ie, np] = angc
            anis.psct[ie, np] = psct2
            anis.index[np] = 1
          if kin6[j] == 2:
            anis.psct[ie, np] = peqin6[j, ie]
            anis.index[np] = 2
          if not (ie > 1):
            large.rgas[np] = rgas6
            large.ein[np] = ei6[j] / rgas6
            l = 29
            if ei6[j] < 0.0: l = 30
            large.ipn[np] = 0
            large.iarry[np] = l
            scrip.dscrpt[np] = fstr(scrp6[4 + nion6 + natt6 + j], 50)
            large.penfra[1, np] = penfra6[1, j]
            large.penfra[2, np] = penfra6[2, j] * 0.000001 / sqrt(3.0)
            large.penfra[3, np] = penfra6[3, j]
    large.last = np
    large.isize = 1
    if large.last >= 2: large.isize = 2
    if large.last >= 4: large.isize = 4
    if large.last >= 8: large.isize = 8
    if large.last >= 16: large.isize = 16
    if large.last >= 32: large.isize = 32
    if large.last >= 64: large.isize = 64
    if large.last >= 128: large.isize = 128
    if large.last >= 256: large.isize = 256
    if large.last >= 512: large.isize = 512
    if large.last >= 1024: large.isize = 1024
    # CAN INCREASE ARRAY SIZE UP TO 1800 IF MORE COMPLEX MIXTURES USED.
    # 1800 = 6 * 300 ( 6 = MAX NO OF GASES. 300 = MAX NO OF LEVELS )
    if large.last > 960: fwrite(fmt992)
    if large.last > 960: fstop()
    # CALCULATION OF TOTAL COLLISION FREQUENCY
    large.tcf[ie] = 0.0
    fdo `if`, 1, large.last:
      large.tcf[ie] = large.tcf[ie] + large.cf[ie, `if`]
      if large.cf[ie, `if`] < 0.0: fwrite(fmt776, large.cf[ie, `if`], ie, `if`, large.iarry[`if`], large.ein[`if`])
    fdo `if`, 1, large.last:
      block L620:
        block L615:
          if large.tcf[ie] == 0.0: break L615
          large.cf[ie, `if`] = large.cf[ie, `if`] / large.tcf[ie]
          break L620
        large.cf[ie, `if`] = 0.0
    fdo `if`, 2, large.last:
      large.cf[ie, `if`] = large.cf[ie, `if`] + large.cf[ie, `if` - 1]
    fred.fcatt[ie] = fred.fcatt[ie] * mix2.eroot[ie]
    fred.fcion[ie] = fred.fcion[ie] * mix2.eroot[ie]
    large.tcf[ie] = large.tcf[ie] * mix2.eroot[ie]
  # WRITE(6,841) (INDEX(J),J, J=1,IPLAST)
  # 841 FORMAT(2X,' INDEX=',I3,' J=',I3)
  # SET ANISOTROPIC FLAG IF ANISOTROPIC SCATTERING DATA IS DETECTED
  kelsum = 0
  fdo j, 1, 6:
    kelsum = kelsum + kel1[j] + kel2[j] + kel3[j] + kel4[j] + kel5[j] + kel6[j]
  fdo j, 1, 250:
    kelsum = kelsum + kin1[j] + kin2[j] + kin3[j] + kin4[j] + kin5[j] + kin6[j]
  if kelsum > 0: anis.niso = 1
  # IF(NISO.EQ.1) WRITE(6,7765) NISO
  # 7765 FORMAT(3X,' ANISOTROPIC SCATTERING DETECTED NISO=',I5)
  # CALCULATE NULL COLLISION FREQUENCY
  bp = setp.efield * setp.efield * cnsts1.const1
  f2 = setp.efield * cnsts1.const3
  elow = setp.tmax * (setp.tmax * bp - f2 * sqrt(0.5 * inpt.efinal)) / inpt.estep - 1.0
  elow = min(elow, setp.small)
  ehi = setp.tmax * (setp.tmax * bp + f2 * sqrt(0.5 * inpt.efinal)) / inpt.estep + 1.0
  if ehi > 10000.0: ehi = 10000.0
  fdo i, 1, 8:
    jlow = int(float64(4000 - 500 * (9 - i) + 1) + trunc(elow))
    jhi = int(float64(4000 - 500 * (8 - i)) + trunc(ehi))
    jlow = int(dmax0(jlow, none))
    jhi = int(dmin0(jhi, n4000))
    fdo j, jlow, jhi:
      if large.tcf[j] >= large.tcfmax[i]: large.tcfmax[i] = large.tcf[j]
  # CROSS SECTION DATA FOR INTEGRALS IN  OUTPUT
  fdo i, 1, inpt.nstep:
    mix2.qt1[i] = ratio.an1 * q1[1, i] + ratio.an2 * q2[1, i] + ratio.an3 * q3[1, i] + ratio.an4 * q4[1, i] + ratio.an5 * q5[1, i] + ratio.an6 * q6[1, i]
    mix2.qt4[i] = ratio.an1 * q1[2, i] + ratio.an2 * q2[2, i] + ratio.an3 * q3[2, i] + ratio.an4 * q4[2, i] + ratio.an5 * q5[2, i] + ratio.an6 * q6[2, i]
    mix1.qion[1, i] = q1[3, i] * ratio.an1
    if nion1 > 1:
      mix1.qion[1, i] = 0.0
      fdo kion, 1, nion1:
        mix1.qion[1, i] = mix1.qion[1, i] + qion1[kion, i] * ratio.an1
    mix1.qion[2, i] = q2[3, i] * ratio.an2
    if nion2 > 1:
      mix1.qion[2, i] = 0.0
      fdo kion, 1, nion2:
        mix1.qion[2, i] = mix1.qion[2, i] + qion2[kion, i] * ratio.an2
    mix1.qion[3, i] = q3[3, i] * ratio.an3
    if nion3 > 1:
      mix1.qion[3, i] = 0.0
      fdo kion, 1, nion3:
        mix1.qion[3, i] = mix1.qion[3, i] + qion3[kion, i] * ratio.an3
    mix1.qion[4, i] = q4[3, i] * ratio.an4
    if nion4 > 1:
      mix1.qion[4, i] = 0.0
      fdo kion, 1, nion4:
        mix1.qion[4, i] = mix1.qion[4, i] + qion4[kion, i] * ratio.an4
    mix1.qion[5, i] = q5[3, i] * ratio.an5
    if nion5 > 1:
      mix1.qion[5, i] = 0.0
      fdo kion, 1, nion5:
        mix1.qion[5, i] = mix1.qion[5, i] + qion5[kion, i] * ratio.an5
    mix1.qion[6, i] = q6[3, i] * ratio.an6
    if nion6 > 1:
      mix1.qion[6, i] = 0.0
      fdo kion, 1, nion6:
        mix1.qion[6, i] = mix1.qion[6, i] + qion6[kion, i] * ratio.an6
    qatt[1, i] = q1[4, i] * ratio.an1
    if natt1 > 1:
      qatt[1, i] = 0.0
      fdo katt, 1, natt1:
        qatt[1, i] = qatt[1, i] + qatt1[katt, i] * ratio.an1
    qatt[2, i] = q2[4, i] * ratio.an2
    if natt2 > 1:
      qatt[2, i] = 0.0
      fdo katt, 1, natt2:
        qatt[2, i] = qatt[2, i] + qatt2[katt, i] * ratio.an2
    qatt[3, i] = q3[4, i] * ratio.an3
    if natt3 > 1:
      qatt[3, i] = 0.0
      fdo katt, 1, natt3:
        qatt[3, i] = qatt[3, i] + qatt3[katt, i] * ratio.an3
    qatt[4, i] = q4[4, i] * ratio.an4
    if natt4 > 1:
      qatt[4, i] = 0.0
      fdo katt, 1, natt4:
        qatt[4, i] = qatt[4, i] + qatt4[katt, i] * ratio.an4
    qatt[5, i] = q5[4, i] * ratio.an5
    if natt5 > 1:
      qatt[5, i] = 0.0
      fdo katt, 1, natt5:
        qatt[5, i] = qatt[5, i] + qatt5[katt, i] * ratio.an5
    qatt[6, i] = q6[4, i] * ratio.an6
    if natt6 > 1:
      qatt[6, i] = 0.0
      fdo katt, 1, natt5:
        qatt[5, i] = qatt[5, i] + qatt5[katt, i] * ratio.an6
    mix2.qt2[i] = 0.0
    mix1.qsatt[i] = 0.0
    mix1.qsum[i] = 0.0
    fdo j, 1, inpt.ngas:
      mix1.qsum[i] = mix1.qsum[i] + mix1.qion[j, i] + qatt[j, i]
      mix1.qsatt[i] = mix1.qsatt[i] + qatt[j, i]
      mix2.qt2[i] = mix2.qt2[i] + mix1.qion[j, i] - qatt[j, i]
    if not (mix3.nin1 == 0):
      fdo j, 1, mix3.nin1:
        mix1.qsum[i] = mix1.qsum[i] + mix1.qin1[j, i] * ratio.an1
    if not (mix3.nin2 == 0):
      fdo j, 1, mix3.nin2:
        mix1.qsum[i] = mix1.qsum[i] + mix1.qin2[j, i] * ratio.an2
    if not (mix3.nin3 == 0):
      fdo j, 1, mix3.nin3:
        mix1.qsum[i] = mix1.qsum[i] + mix1.qin3[j, i] * ratio.an3
    if not (mix3.nin4 == 0):
      fdo j, 1, mix3.nin4:
        mix1.qsum[i] = mix1.qsum[i] + mix1.qin4[j, i] * ratio.an4
    if not (mix3.nin5 == 0):
      fdo j, 1, mix3.nin5:
        mix1.qsum[i] = mix1.qsum[i] + mix1.qin5[j, i] * ratio.an5
    if not (mix3.nin6 == 0):
      fdo j, 1, mix3.nin6:
        mix1.qsum[i] = mix1.qsum[i] + mix1.qin6[j, i] * ratio.an6
  return
