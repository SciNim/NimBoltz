## SUBROUTINE GAS26 of magboltz-11.21.f (line 64918).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons

proc gas26*(q: var Arr2[float64, 6, 4000], qin: var Arr2[float64, 250, 4000], nin: var int, ev: var array[1..6, float64], ein: var array[1..250, float64], name: var string, virial: float64, eoby: var array[1..30, float64], peqel: var Arr2[float64, 6, 4000], peqin: var Arr2[float64, 250, 4000], penfra: Arr2[float64, 3, 250], kel: var array[1..6, int], kin: var array[1..250, int], qion: Arr2[float64, 30, 4000], peqion: Arr2[float64, 30, 4000], eion: array[1..30, float64], nion: var int, qatt: var Arr2[float64, 8, 4000], natt: var int, qnull: Arr2[float64, 10, 4000], nnull: var int, scln: array[1..10, float64], nc0: array[1..30, int], ec0: array[1..30, float64], wklm: array[1..30, float64], efl: array[1..30, float64], ng1: array[1..30, int], eg1: array[1..30, float64], ng2: array[1..30, int], eg2: array[1..30, float64], scrpt: var array[1..300, string], scrptn: array[1..10, string]) =
  var j, i: int
  var en: float64
  # NEW ANALYSIS UPDATED TO DECEMBER 1994
  name = fstr("REID STEP(ANIS)", 25)
  nion = 1
  natt = 1
  nin = 1
  nnull = 0
  fdo j, 1, 6:
    kel[j] = 0
  fdo j, 1, nin:
    kin[j] = 0
  # ANISOTROPIC SCATTERING FOR LEVEL 1
  kin[1] = 1
  ev[1] = 0.0
  ev[2] = 2.0 * cnsts.emass / (4.0 * cnsts.amu)
  ev[3] = 15.899999618530273
  ev[4] = 0.0
  ev[5] = 0.0
  ev[6] = 0.0
  eoby[1] = 15.899999618530273
  ein[1] = 0.20000000298023224
  scrpt[1] = fstr("                              ", 50)
  scrpt[2] = fstr(" ELASTIC       REID STEP(ANIS)", 50)
  scrpt[3] = fstr(" IONISATION    ELOSS= 15.90   ", 50)
  scrpt[4] = fstr("                              ", 50)
  scrpt[5] = fstr("                              ", 50)
  scrpt[6] = fstr("                              ", 50)
  scrpt[7] = fstr(" EXC (ANIS)    ELOSS=  0.20   ", 50)
  en = -(inpt.estep / 2.0)
  fdo i, 1, inpt.nstep:
    peqel[2, i] = 0.0
    peqin[1, i] = 0.0
    en = en + inpt.estep
    q[2, i] = 1.0000000168623835e-16
    q[3, i] = 0.0
    q[4, i] = 0.0
    qatt[1, i] = q[4, i]
    q[5, i] = 0.0
    q[6, i] = 0.0
    qin[1, i] = 0.0
    if not (en <= ein[1]):
      qin[1, i] = 1.0000000036274937e-15
      peqin[1, i] = 0.5 + (qin[1, i] - 0.699999988079071 * qin[1, i]) / qin[1, i]
    q[1, i] = q[2, i] + q[3, i] + q[4, i] + 0.699999988079071 * qin[1, i]
  # SAVE COMPUTE TIME
  if inpt.efinal <= ein[1]: nin = 0
  return
