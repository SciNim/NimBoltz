## SUBROUTINE GASMIX of magboltz-11.21.f (line 2135).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import gas28, gas29, gas30, gas31, gas32, gas33, gas34, gas35, gas36, gas37, gas38, gas39, gas40, gas41, gas42, gas43, gas44, gas45, gas46, gas47, gas48, gas49, gas50, gas51, gas52, gas53, gas54, gas55, gas56, gas57, gas58, gas59, gas60, gas61, gas62, gas63, gas64, gas65, gas66, gas67, gas68, gas69, gas70, gas71, gas72, gas73, gas74, gas75, gas76, gas77, gas78, gas79, gas80

proc gasmix*(ngs: int, q: var Arr2[float64, 6, 4000], qin: var Arr2[float64, 250, 4000], nin: var int, ev: var array[1..6, float64], ei: var array[1..250, float64], name: var string, virl: float64, eb: var array[1..30, float64], peqel: var Arr2[float64, 6, 4000], peqin: var Arr2[float64, 250, 4000], penfra: var Arr2[float64, 3, 250], kel: var array[1..6, int], kin: var array[1..250, int], qion: var Arr2[float64, 30, 4000], peqion: var Arr2[float64, 30, 4000], eion: var array[1..30, float64], nion: var int, qatt: var Arr2[float64, 8, 4000], natt: var int, qnull: var Arr2[float64, 10, 4000], nnull: var int, scln: var array[1..10, float64], nc0: var array[1..30, int], ec0: var array[1..30, float64], wk: var array[1..30, float64], efl: var array[1..30, float64], ng1: var array[1..30, int], eg1: var array[1..30, float64], ng2: var array[1..30, int], eg2: var array[1..30, float64], scrpt: var array[1..300, string], scrptn: var array[1..10, string]) =
  var i, j: int
  fdo i, 1, 4000:
    fdo j, 1, 250:
      qin[j, i] = 0.0
      peqin[j, i] = 0.0
    fdo j, 1, 30:
      qion[j, i] = 0.0
      peqion[j, i] = 0.0
    fdo j, 1, 6:
      q[j, i] = 0.0
      peqel[j, i] = 0.0
    fdo j, 1, 8:
      qatt[j, i] = 0.0
    fdo j, 1, 10:
      qnull[j, i] = 0.0
  fdo j, 1, 250:
    ei[j] = 0.0
    kin[j] = 0
    penfra[1, j] = 0.0
    penfra[2, j] = 0.0
    penfra[3, j] = 0.0
  fdo j, 1, 30:
    eion[j] = 0.0
    eb[j] = 0.0
    nc0[j] = 0
    ec0[j] = 0.0
    wk[j] = 0.0
    efl[j] = 0.0
    ng1[j] = 0
    eg1[j] = 0.0
    ng2[j] = 0
    eg2[j] = 0.0
  fdo j, 1, 6:
    ev[j] = 0.0
    kel[j] = 0
  fdo j, 1, 10:
    scln[j] = 0.0
  case ngs
  of 28: gas28(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 29: gas29(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 30: gas30(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 31: gas31(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 32: gas32(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 33: gas33(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 34: gas34(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 35: gas35(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 36: gas36(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 37: gas37(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 38: gas38(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 39: gas39(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 40: gas40(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 41: gas41(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 42: gas42(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 43: gas43(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 44: gas44(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 45: gas45(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 46: gas46(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 47: gas47(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 48: gas48(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 49: gas49(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 50: gas50(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 51: gas51(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 52: gas52(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 53: gas53(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 54: gas54(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 55: gas55(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 56: gas56(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 57: gas57(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 58: gas58(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 59: gas59(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 60: gas60(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 61: gas61(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 62: gas62(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 63: gas63(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 64: gas64(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 65: gas65(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 66: gas66(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 67: gas67(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 68: gas68(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 69: gas69(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 70: gas70(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 71: gas71(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 72: gas72(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 73: gas73(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 74: gas74(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 75: gas75(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 76: gas76(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 77: gas77(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 78: gas78(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 79: gas79(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  of 80: gas80(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  else: raiseAssert "gas " & $ngs & " is not in legacy/"
