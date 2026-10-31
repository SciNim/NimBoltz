## SUBROUTINE GAS28 of magboltz-11.21.f (line 65049).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons

proc gas28*(q: var Arr2[float64, 6, 4000], qin: var Arr2[float64, 250, 4000], nin: var int, ev: var array[1..6, float64], ein: var array[1..250, float64], name: var string, virial: float64, eoby: var array[1..30, float64], peqel: Arr2[float64, 6, 4000], peqin: Arr2[float64, 250, 4000], penfra: Arr2[float64, 3, 250], kel: var array[1..6, int], kin: var array[1..250, int], qion: Arr2[float64, 30, 4000], peqion: Arr2[float64, 30, 4000], eion: array[1..30, float64], nion: var int, qatt: var Arr2[float64, 8, 4000], natt: var int, qnull: Arr2[float64, 10, 4000], nnull: var int, scln: array[1..10, float64], nc0: array[1..30, int], ec0: array[1..30, float64], wklm: array[1..30, float64], efl: array[1..30, float64], ng1: array[1..30, int], eg1: array[1..30, float64], ng2: array[1..30, int], eg2: array[1..30, float64], scrpt: var array[1..300, string], scrptn: array[1..10, string]) =
  var j, i: int
  var en: float64
  # NEW ANALYSIS UPDATED TO DECEMBER 1994
  name = fstr(" REID RAMP S=10", 25)
  nion = 1
  natt = 1
  nin = 1
  nnull = 0
  fdo j, 1, 6:
    kel[j] = 0
  fdo j, 1, nin:
    kin[j] = 0
  ev[1] = 0.0
  ev[2] = 2.0 * cnsts.emass / (4.0 * cnsts.amu)
  ev[3] = 15.899999618530273
  ev[4] = 0.0
  ev[5] = 0.0
  ev[6] = 0.0
  eoby[1] = 15.899999618530273
  ein[1] = 0.20000000298023224
  scrpt[1] = fstr("                              ", 50)
  scrpt[2] = fstr(" ELASTIC       REID RAMP      ", 50)
  scrpt[3] = fstr("                              ", 50)
  scrpt[4] = fstr("                              ", 50)
  scrpt[5] = fstr("                              ", 50)
  scrpt[6] = fstr("                              ", 50)
  scrpt[7] = fstr(" EXC           ELOSS=  0.2    ", 50)
  en = -(inpt.estep / 2.0)
  fdo i, 1, inpt.nstep:
    en = en + inpt.estep
    q[2, i] = 6.000000233523199e-16
    q[3, i] = 0.0
    q[4, i] = 0.0
    qatt[1, i] = q[4, i]
    q[5, i] = 0.0
    q[6, i] = 0.0
    qin[1, i] = 0.0
    if not (en <= ein[1]):
      qin[1, i] = (en - ein[1]) * 1.0000000036274937e-15
    q[1, i] = q[2, i] + q[3, i] + q[4, i] + qin[1, i]
  # SAVE COMPUTE TIME
  if inpt.efinal <= ein[1]: nin = 0
  return
