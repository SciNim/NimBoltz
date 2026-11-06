## Model gases of Magboltz 11.21 for tests of the transport: GAS26 ("REID
## STEP(ANIS)", Reid's step model with an anisotropic inelastic level) and
## GAS27 (" MAXWELL 1994--", a constant collision frequency) and GAS28
## (" REID RAMP S=10", Reid's ramp model).

import common

proc reidStep*(inp: GasInput, g: var GasXsec) =
  ## constant elastic 1e-16 cm² (isotropic) and a 0.2 eV level of 1e-15 cm²
  ## of which 70 % counts for momentum transfer
  g.describe(inp, "REID STEP(ANIS)", "REID STEP(ANIS)", mass = 4.0)
  for k in 1 .. 6: g.kel[k] = 0
  g.scrp[2] = " ELASTIC       REID STEP(ANIS)"
  g.e[3] = 15.90'f32
  g.shells [(0.0'f32, 15.90'f32, 0, 0.0'f32, " IONISATION    ELOSS= 15.90")]
  g.attachment [""]
  g.levels(inp, [(0.2'f32, " EXC (ANIS)    ELOSS=  0.20")])
  g.kin[1] = 1
  let share = float64(0.7'f32)
  for i, en in inp.steps:
    g.peqel[2, i] = 0.0
    g.peqin[1, i] = 0.0
    g.q[2, i] = float64(1.0e-16'f32)
    g.q[3, i] = 0.0
    g.q[4, i] = 0.0
    g.qatt[1, i] = 0.0
    g.q[5, i] = 0.0
    g.q[6, i] = 0.0
    g.qin[1, i] = 0.0
    if en > g.ei[1]:
      g.qin[1, i] = float64(10.0e-16'f32)
      g.peqin[1, i] = 0.5 + (g.qin[1, i] - share * g.qin[1, i]) / g.qin[1, i]
    g.q[1, i] = g.q[2, i] + g.q[3, i] + g.q[4, i] + share * g.qin[1, i]
  g.dropLevelsAbove(inp.efinal)

proc maxwell*(inp: GasInput, g: var GasXsec) =
  ## elastic 6e-16 cm² / sqrt(E / eV): constant collision frequency
  g.describe(inp, " MAXWELL 1994--", " MAXWELL 1994--", mass = 4.0)
  for k in 1 .. 6: g.kel[k] = 0
  g.scrp[2] = " ELASTIC       MAXWELL"
  g.e[3] = 99.0'f32
  g.shells [(0.0'f32, 99.0'f32, 0, 0.0'f32, "")]
  g.attachment [""]
  g.nin = 0
  let sigc = float64(6.0e-16'f32)
  for i, en in inp.steps:
    g.q[2, i] = if en == 0.0: float64(100000.0e-16'f32) else: sigc / sqrt(en)
    g.q[3, i] = 0.0
    g.q[4, i] = 0.0
    g.qatt[1, i] = 0.0
    g.q[5, i] = 0.0
    g.q[6, i] = 0.0
    g.q[1, i] = g.q[2, i]

proc reidRamp*(inp: GasInput, g: var GasXsec) =
  ## constant elastic 6e-16 cm² and a 0.2 eV level rising as
  ## (E - 0.2 eV) · 1e-15 cm²/eV
  g.describe(inp, " REID RAMP S=10", " REID RAMP S=10", mass = 4.0)
  for k in 1 .. 6: g.kel[k] = 0
  g.scrp[2] = " ELASTIC       REID RAMP"
  g.e[3] = 15.90'f32
  g.shells [(0.0'f32, 15.9'f32, 0, 0.0'f32, "")]
  g.attachment [""]
  g.levels(inp, [(0.2'f32, " EXC           ELOSS=  0.2")])
  g.kin[1] = 0
  let slope = float64(10.0e-16'f32)
  for i, en in inp.steps:
    g.q[2, i] = float64(6.0e-16'f32)
    g.q[3, i] = 0.0
    g.q[4, i] = 0.0
    g.qatt[1, i] = 0.0
    g.q[5, i] = 0.0
    g.q[6, i] = 0.0
    g.qin[1, i] = 0.0
    if en > g.ei[1]:
      g.qin[1, i] = (en - g.ei[1]) * slope
    g.q[1, i] = g.q[2, i] + g.q[3, i] + g.q[4, i] + g.qin[1, i]
  g.dropLevelsAbove(inp.efinal)
