## SUBROUTINE GAS27 of magboltz-11.21.f (line 64988).
## Translated by f2nim; maintained by hand since (see docs/DESIGN.md).

import ../frt, commons

proc gas27*(q: var Arr2[float64, 6, 4000], qin: Arr2[float64, 250, 4000], nin: var int, ev: var array[1..6, float64], ein: array[1..250, float64], name: var string, virial: float64, eoby: var array[1..30, float64], peqel: Arr2[float64, 6, 4000], peqin: Arr2[float64, 250, 4000], penfra: Arr2[float64, 3, 250], kel: var array[1..6, int], kin: array[1..250, int], qion: Arr2[float64, 30, 4000], peqion: Arr2[float64, 30, 4000], eion: array[1..30, float64], nion: var int, qatt: var Arr2[float64, 8, 4000], natt: var int, qnull: Arr2[float64, 10, 4000], nnull: var int, scln: array[1..10, float64], nc0: array[1..30, int], ec0: array[1..30, float64], wklm: array[1..30, float64], efl: array[1..30, float64], ng1: array[1..30, int], eg1: array[1..30, float64], ng2: array[1..30, int], eg2: array[1..30, float64], scrpt: var array[1..300, string], scrptn: array[1..10, string]) =
  var j, i: int
  var sigc, en: float64
  # MAXWELL MODEL DECEMBER 1994
  name = fstr(" MAXWELL 1994--", 25)
  nion = 1
  natt = 1
  nin = 0
  nnull = 0
  fdo j, 1, 6:
    kel[j] = 0
  # DO 2 J=1,NIN
  # 2 KIN(J)=0
  sigc = 6.000000233523199e-16
  ev[1] = 0.0
  ev[2] = 2.0 * cnsts.emass / (4.0 * cnsts.amu)
  ev[3] = 99.0
  ev[4] = 0.0
  ev[5] = 0.0
  ev[6] = 0.0
  eoby[1] = 99.0
  scrpt[1] = fstr("                              ", 50)
  scrpt[2] = fstr(" ELASTIC       MAXWELL        ", 50)
  scrpt[3] = fstr("                              ", 50)
  scrpt[4] = fstr("                              ", 50)
  scrpt[5] = fstr("                              ", 50)
  scrpt[6] = fstr("                              ", 50)
  en = -(inpt.estep / 2.0)
  fdo i, 1, inpt.nstep:
    en = en + inpt.estep
    block L10:
      if en == 0.0:
        q[2, i] = 9.999999960041972e-12
        break L10
      q[2, i] = sigc / sqrt(en)
    q[3, i] = 0.0
    q[4, i] = 0.0
    qatt[1, i] = q[4, i]
    q[5, i] = 0.0
    q[6, i] = 0.0
    q[1, i] = q[2, i]
  return
