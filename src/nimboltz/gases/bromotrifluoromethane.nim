## Bromotrifluoromethane: GAS51 of Magboltz 11.21 ("CF3BR -2002--"), an
## older routine: isotropic, ionisation as one channel in Q(3), energies
## stepped by the routine itself.
##
## Elastic (momentum transfer, log-log) and ionisation from tables;
## attachment by a fit below 2 eV; the V1 and V2 vibrations with their
## superelastic partners as Born dipoles (populations at the gas
## temperature), V3 as a dipole, a harmonic and one excitation from tables.

import common

let xel: array[30, float64] = [
  1e-6, 0.001'f32, 0.01'f32, 0.1'f32, 0.4'f32, 0.6'f32, 0.8'f32, 1.0, 1.4'f32,
  2.0, 3.0, 5.0, 7.5, 10.0, 14.0, 20.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0,
  400.0, 700.0, 1000.0, 2000.0, 4000.0, 10000.0, 20000.0, 100000.0]
let yel: array[30, float32] = [
  2800.0, 2800.0, 280.0, 28.0, 7.0, 5.5, 5.0, 5.0, 6.0, 7.0, 9.0, 11.5, 15.5,
  18.0, 18.0, 17.0, 14.0, 10.0, 8.0, 6.0, 4.0, 3.0, 1.45, 0.9, 0.65, 0.3,
  0.15, 0.05, 0.02, 0.004]
let xvibh: array[16, float32] = [
  0.3, 1.0, 3.0, 3.7, 4.7, 5.7, 6.7, 7.7, 9.5, 15.0, 20.0, 50.0, 100.0,
  1000.0, 10000.0, 100000.0]
let yvibh: array[16, float32] = [
  0.0, 0.005, 0.025, 0.1, 0.25, 1.15, 1.75, 1.5, 1.0, 0.08, 0.03, 0.0015,
  0.00005, 5e-6, 5e-7, 5e-8]
let xexc: array[31, float32] = [
  8.7, 9.0, 10.0, 11.0, 12.0, 14.0, 16.0, 21.0, 26.0, 31.0, 36.0, 40.0, 50.0,
  70.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 250.0, 300.0, 400.0, 500.0,
  600.0, 1000.0, 2000.0, 4000.0, 10000.0, 20000.0, 100000.0]
let yexc: array[31, float32] = [
  0.0, 0.07, 0.32, 0.56, 0.85, 1.35, 1.65, 2.15, 2.3, 2.35, 2.4, 2.45, 2.45,
  2.43, 2.41, 2.39, 2.35, 2.3, 2.25, 2.2, 2.05, 1.9, 1.65, 1.4, 1.2, 0.8,
  0.46, 0.27, 0.14, 0.07, 0.02]
let xion: array[40, float32] = [
  11.4, 12.0, 13.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 30.0, 34.0,
  38.0, 42.0, 46.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 175.0,
  200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0,
  2500.0, 3000.0, 5000.0, 10000.0, 20000.0, 100000.0]
let yion: array[40, float32] = [
  0.0, 0.04, 0.13, 0.25, 0.4, 0.66, 1.0, 1.3, 1.7, 2.05, 2.45, 2.85, 3.3, 3.8,
  4.3, 4.6, 5.35, 5.9, 6.4, 6.6, 6.8, 7.0, 7.05, 7.0, 6.9, 6.45, 6.2, 5.5,
  4.75, 4.2, 3.55, 3.05, 2.25, 1.85, 1.6, 1.35, 0.95, 0.6, 0.35, 0.1]

proc bromotrifluoromethane*(inp: GasInput, g: var GasXsec) =
  g.describe(inp, " CF3BR -2002--", " CF3BR -2002--", mass = 148.90991)
  for k in 1 .. 6: g.kel[k] = 0
  g.scrp[2] = " ELASTIC         CF3BR"
  g.e[3] = 11.40'f32
  g.shells [(0.0'f32, 11.40'f32, 0, 0.0'f32, " IONISATION    ELOSS= 11.40")]
  g.attachment [" ATTACHMENT"]
  g.levels(inp, [
    (-0.050'f32, " VIB1          ELOSS= -0.050"),
    (0.050, " VIB1          ELOSS=  0.050"),
    (-0.1345, " VIB2          ELOSS= -0.1345"),
    (0.1345, " VIB2          ELOSS=  0.1345"),
    (0.024999, " VIB3          ELOSS=  0.024999"),
    (0.30, " VIB HAR       ELOSS=  0.30"),
    (8.7, " EXC           ELOSS=  8.7")])
  for j in 1 .. g.nin: g.kin[j] = 0

  # dipole constants of V1, V2, V3
  let amp = [float64(0.100'f32), float64(0.110'f32), float64(0.750'f32)]
  let apop = [exp(g.ei[1] / inp.kT), exp(g.ei[3] / inp.kT)]
  for i, en in inp.steps:
    g.q[2, i] = if en == 0.0: 2800e-16 else: logLog(xel, yel, en) * 1e-16
    g.q[3, i] = if en < g.e[3]: 0.0 else: lin(xion, yion, en) * 1e-16
    g.q[4, i] = 0.0
    if en != 0.0 and en <= 2.0:
      let at1 = float64(0.0353'f32) / sqrt(en)
      let at2 = exp(-en * en / 9e-6)
      let at3 = exp(-en / float64(0.0588'f32))
      g.q[4, i] = (at1 * at2 + at3) * 27.4e-16
    g.qatt[1, i] = g.q[4, i]
    g.q[5, i] = 0.0
    g.q[6, i] = 0.0
    for k in 0 .. 1:
      let (jd, je) = (2 * k + 1, 2 * k + 2)
      g.qin[jd, i] = 0.0
      if en != 0.0:
        g.qin[jd, i] = vibrationalDipole(amp[k], g.ei[jd], en)
        g.qin[jd, i] = g.qin[jd, i] * apop[k] / (1.0 + apop[k]) * 1e-16
      g.qin[je, i] = 0.0
      if en > g.ei[je]:
        g.qin[je, i] = vibrationalDipole(amp[k], g.ei[je], en)
        g.qin[je, i] = g.qin[je, i] * 1.0 / (1.0 + apop[k]) * 1e-16
    g.qin[5, i] = 0.0
    if en > g.ei[5]:
      g.qin[5, i] = vibrationalDipole(amp[2], g.ei[5], en) * 1e-16
    g.qin[6, i] = 0.0
    if en > g.ei[6]: g.qin[6, i] = lin(xvibh, yvibh, en) * 1e-16
    g.qin[7, i] = 0.0
    if en > g.ei[7]: g.qin[7, i] = lin(xexc, yexc, en) * 1e-16
    var q = g.q[2, i] + g.q[3, i] + g.q[4, i]
    for j in 1 .. 7: q = q + g.qin[j, i]
    g.q[1, i] = q
  g.dropLevelsAbove(inp.efinal)
