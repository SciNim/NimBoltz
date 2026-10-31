## SUBROUTINE GASMIX of magboltz-11.21.f (line 2135).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons
import gas1, gas10, gas11, gas12, gas13, gas14, gas15, gas16, gas17, gas18, gas19, gas2, gas20, gas21, gas22, gas23, gas24, gas25, gas26, gas27, gas28, gas29, gas3, gas30, gas31, gas32, gas33, gas34, gas35, gas36, gas37, gas38, gas39, gas4, gas40, gas41, gas42, gas43, gas44, gas45, gas46, gas47, gas48, gas49, gas5, gas50, gas51, gas52, gas53, gas54, gas55, gas56, gas57, gas58, gas59, gas6, gas60, gas61, gas62, gas63, gas64, gas65, gas66, gas67, gas68, gas69, gas7, gas70, gas71, gas72, gas73, gas74, gas75, gas76, gas77, gas78, gas79, gas8, gas80, gas9

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
  block L80:
    block L79:
      block L78:
        block L77:
          block L76:
            block L75:
              block L74:
                block L73:
                  block L72:
                    block L71:
                      block L70:
                        block L69:
                          block L68:
                            block L67:
                              block L66:
                                block L65:
                                  block L64:
                                    block L63:
                                      block L62:
                                        block L61:
                                          block L60:
                                            block L59:
                                              block L58:
                                                block L57:
                                                  block L56:
                                                    block L55:
                                                      block L54:
                                                        block L53:
                                                          block L52:
                                                            block L51:
                                                              block L50:
                                                                block L49:
                                                                  block L48:
                                                                    block L47:
                                                                      block L46:
                                                                        block L45:
                                                                          block L44:
                                                                            block L43:
                                                                              block L42:
                                                                                block L41:
                                                                                  block L40:
                                                                                    block L39:
                                                                                      block L38:
                                                                                        block L37:
                                                                                          block L36:
                                                                                            block L35:
                                                                                              block L34:
                                                                                                block L33:
                                                                                                  block L32:
                                                                                                    block L31:
                                                                                                      block L30:
                                                                                                        block L29:
                                                                                                          block L28:
                                                                                                            block L27:
                                                                                                              block L26:
                                                                                                                block L25:
                                                                                                                  block L24:
                                                                                                                    block L23:
                                                                                                                      block L22:
                                                                                                                        block L21:
                                                                                                                          block L20:
                                                                                                                            block L19:
                                                                                                                              block L18:
                                                                                                                                block L17:
                                                                                                                                  block L16:
                                                                                                                                    block L15:
                                                                                                                                      block L14:
                                                                                                                                        block L13:
                                                                                                                                          block L12:
                                                                                                                                            block L11:
                                                                                                                                              block L10:
                                                                                                                                                block L9:
                                                                                                                                                  block L8:
                                                                                                                                                    block L7:
                                                                                                                                                      block L6:
                                                                                                                                                        block L5:
                                                                                                                                                          block L4:
                                                                                                                                                            block L3:
                                                                                                                                                              block L2:
                                                                                                                                                                block L1:
                                                                                                                                                                  case ngs
                                                                                                                                                                  of 1: break L1
                                                                                                                                                                  of 2: break L2
                                                                                                                                                                  of 3: break L3
                                                                                                                                                                  of 4: break L4
                                                                                                                                                                  of 5: break L5
                                                                                                                                                                  of 6: break L6
                                                                                                                                                                  of 7: break L7
                                                                                                                                                                  of 8: break L8
                                                                                                                                                                  of 9: break L9
                                                                                                                                                                  of 10: break L10
                                                                                                                                                                  of 11: break L11
                                                                                                                                                                  of 12: break L12
                                                                                                                                                                  of 13: break L13
                                                                                                                                                                  of 14: break L14
                                                                                                                                                                  of 15: break L15
                                                                                                                                                                  of 16: break L16
                                                                                                                                                                  of 17: break L17
                                                                                                                                                                  of 18: break L18
                                                                                                                                                                  of 19: break L19
                                                                                                                                                                  of 20: break L20
                                                                                                                                                                  of 21: break L21
                                                                                                                                                                  of 22: break L22
                                                                                                                                                                  of 23: break L23
                                                                                                                                                                  of 24: break L24
                                                                                                                                                                  of 25: break L25
                                                                                                                                                                  of 26: break L26
                                                                                                                                                                  of 27: break L27
                                                                                                                                                                  of 28: break L28
                                                                                                                                                                  of 29: break L29
                                                                                                                                                                  of 30: break L30
                                                                                                                                                                  of 31: break L31
                                                                                                                                                                  of 32: break L32
                                                                                                                                                                  of 33: break L33
                                                                                                                                                                  of 34: break L34
                                                                                                                                                                  of 35: break L35
                                                                                                                                                                  of 36: break L36
                                                                                                                                                                  of 37: break L37
                                                                                                                                                                  of 38: break L38
                                                                                                                                                                  of 39: break L39
                                                                                                                                                                  of 40: break L40
                                                                                                                                                                  of 41: break L41
                                                                                                                                                                  of 42: break L42
                                                                                                                                                                  of 43: break L43
                                                                                                                                                                  of 44: break L44
                                                                                                                                                                  of 45: break L45
                                                                                                                                                                  of 46: break L46
                                                                                                                                                                  of 47: break L47
                                                                                                                                                                  of 48: break L48
                                                                                                                                                                  of 49: break L49
                                                                                                                                                                  of 50: break L50
                                                                                                                                                                  of 51: break L51
                                                                                                                                                                  of 52: break L52
                                                                                                                                                                  of 53: break L53
                                                                                                                                                                  of 54: break L54
                                                                                                                                                                  of 55: break L55
                                                                                                                                                                  of 56: break L56
                                                                                                                                                                  of 57: break L57
                                                                                                                                                                  of 58: break L58
                                                                                                                                                                  of 59: break L59
                                                                                                                                                                  of 60: break L60
                                                                                                                                                                  of 61: break L61
                                                                                                                                                                  of 62: break L62
                                                                                                                                                                  of 63: break L63
                                                                                                                                                                  of 64: break L64
                                                                                                                                                                  of 65: break L65
                                                                                                                                                                  of 66: break L66
                                                                                                                                                                  of 67: break L67
                                                                                                                                                                  of 68: break L68
                                                                                                                                                                  of 69: break L69
                                                                                                                                                                  of 70: break L70
                                                                                                                                                                  of 71: break L71
                                                                                                                                                                  of 72: break L72
                                                                                                                                                                  of 73: break L73
                                                                                                                                                                  of 74: break L74
                                                                                                                                                                  of 75: break L75
                                                                                                                                                                  of 76: break L76
                                                                                                                                                                  of 77: break L77
                                                                                                                                                                  of 78: break L78
                                                                                                                                                                  of 79: break L79
                                                                                                                                                                  of 80: break L80
                                                                                                                                                                  else: discard
                                                                                                                                                                gas1(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                                                return
                                                                                                                                                              gas2(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                                              return
                                                                                                                                                            gas3(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                                            return
                                                                                                                                                          gas4(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                                          return
                                                                                                                                                        gas5(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                                        return
                                                                                                                                                      gas6(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                                      return
                                                                                                                                                    gas7(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                                    return
                                                                                                                                                  gas8(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                                  return
                                                                                                                                                gas9(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                                return
                                                                                                                                              gas10(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                              return
                                                                                                                                            gas11(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                            return
                                                                                                                                          gas12(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                          return
                                                                                                                                        gas13(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                        return
                                                                                                                                      gas14(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                      return
                                                                                                                                    gas15(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                    return
                                                                                                                                  gas16(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                  return
                                                                                                                                gas17(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                                return
                                                                                                                              gas18(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                              return
                                                                                                                            gas19(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                            return
                                                                                                                          gas20(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                          return
                                                                                                                        gas21(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                        return
                                                                                                                      gas22(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                      return
                                                                                                                    gas23(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                    return
                                                                                                                  gas24(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                  return
                                                                                                                gas25(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                                return
                                                                                                              gas26(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                              return
                                                                                                            gas27(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                            return
                                                                                                          gas28(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                          return
                                                                                                        gas29(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                        return
                                                                                                      gas30(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                      return
                                                                                                    gas31(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                    return
                                                                                                  gas32(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                  return
                                                                                                gas33(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                                return
                                                                                              gas34(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                              return
                                                                                            gas35(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                            return
                                                                                          gas36(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                          return
                                                                                        gas37(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                        return
                                                                                      gas38(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                      return
                                                                                    gas39(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                    return
                                                                                  gas40(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                  return
                                                                                gas41(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                                return
                                                                              gas42(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                              return
                                                                            gas43(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                            return
                                                                          gas44(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                          return
                                                                        gas45(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                        return
                                                                      gas46(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                      return
                                                                    gas47(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                    return
                                                                  gas48(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                  return
                                                                gas49(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                                return
                                                              gas50(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                              return
                                                            gas51(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                            return
                                                          gas52(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                          return
                                                        gas53(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                        return
                                                      gas54(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                      return
                                                    gas55(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                    return
                                                  gas56(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                  return
                                                gas57(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                                return
                                              gas58(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                              return
                                            gas59(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                            return
                                          gas60(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                          return
                                        gas61(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                        return
                                      gas62(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                      return
                                    gas63(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                    return
                                  gas64(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                  return
                                gas65(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                                return
                              gas66(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                              return
                            gas67(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                            return
                          gas68(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                          return
                        gas69(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                        return
                      gas70(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                      return
                    gas71(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                    return
                  gas72(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                  return
                gas73(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
                return
              gas74(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
              return
            gas75(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
            return
          gas76(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
          return
        gas77(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
        return
      gas78(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
      return
    gas79(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
    return
  gas80(q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel, kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0, ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  return
