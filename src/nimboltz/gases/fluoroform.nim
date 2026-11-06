## Fluoroform: GAS50 of Magboltz 11.21 ("CHF3 -2001---"), an older routine:
## isotropic, ionisation as one channel in Q(3), energies stepped by the
## routine itself.
##
## Elastic (momentum transfer, log-log), ionisation and attachment from
## tables; rotation and the V6/V3 vibration with their superelastic partners
## as Born dipoles (populations at the gas temperature), V5/V2/V4 as a
## dipole plus a table, V1 as a dipole, harmonics and one dissociative
## excitation from tables.

import common

let xen: array[32, float64] = [
  1e-6, 0.001'f32, 0.01'f32, 0.1'f32, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 9.0,
  10.0, 15.0, 20.0, 30.0, 35.0, 50.0, 60.0, 75.0, 100.0, 150.0, 200.0, 300.0,
  500.0, 700.0, 1000.0, 2000.0, 4000.0, 10000.0, 20000.0, 100000.0]
let yxsec: array[32, float32] = [
  16300.0, 16300.0, 1630.0, 163.0, 16.3, 11.5, 10.5, 9.5, 9.5, 10.0, 11.0,
  11.0, 10.8, 10.5, 10.0, 9.5, 8.4, 6.5, 5.5, 4.6, 3.54, 2.55, 1.63, 0.96,
  0.54, 0.36, 0.23, 0.1, 0.05, 0.02, 0.01, 0.002]
let xvib3: array[14, float32] = [
  0.1429, 4.0, 4.7, 5.7, 6.7, 7.7, 9.5, 15.0, 20.0, 50.0, 100.0, 1000.0,
  10000.0, 100000.0]
let yvib3: array[14, float32] = [
  0.0, 0.0, 0.109, 0.952, 1.43, 1.22, 0.79, 0.068, 0.023, 0.0001, 0.00001,
  1e-6, 1e-7, 1e-8]
let xvib4: array[16, float32] = [
  0.2858, 1.0, 3.0, 3.7, 4.7, 5.7, 6.7, 7.7, 9.5, 15.0, 20.0, 50.0, 100.0,
  1000.0, 10000.0, 100000.0]
let yvib4: array[16, float32] = [
  0.0, 0.001, 0.005, 0.02, 0.03, 0.23, 0.35, 0.3, 0.2, 0.016, 0.006, 0.0003,
  0.00001, 1e-6, 1e-7, 1e-8]
let xvib6: array[16, float32] = [
  0.4287, 1.0, 3.0, 3.7, 4.7, 5.7, 6.7, 7.7, 9.5, 15.0, 20.0, 50.0, 100.0,
  1000.0, 10000.0, 100000.0]
let yvib6: array[16, float32] = [
  0.0, 0.001, 0.04, 0.08, 0.12, 0.92, 1.4, 1.2, 0.8, 0.064, 0.02, 0.0004,
  0.00004, 4e-6, 4e-7, 4e-8]
let xexc: array[33, float32] = [
  10.7, 11.0, 12.0, 13.0, 14.0, 16.0, 18.0, 23.0, 28.0, 33.0, 38.0, 43.0,
  48.0, 60.0, 70.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 250.0,
  300.0, 400.0, 500.0, 600.0, 1000.0, 2000.0, 4000.0, 10000.0, 20000.0,
  100000.0]
let yexc: array[33, float32] = [
  0.0, 0.05, 0.21, 0.37, 0.57, 0.87, 1.1, 1.42, 1.52, 1.57, 1.6, 1.63, 1.62,
  1.61, 1.6, 1.6, 1.59, 1.57, 1.55, 1.52, 1.49, 1.45, 1.37, 1.27, 1.1, 0.92,
  0.8, 0.53, 0.31, 0.18, 0.09, 0.05, 0.013]
let xion: array[52, float32] = [
  13.86, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0,
  32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0, 46.0, 48.0, 60.0, 70.0, 80.0,
  90.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0,
  500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0,
  2500.0, 3000.0, 5000.0, 10000.0, 20000.0, 40000.0, 100000.0]
let yion: array[52, float32] = [
  0.0, 0.031, 0.072, 0.123, 0.184, 0.266, 0.431, 0.59, 0.843, 1.07, 1.27,
  1.46, 1.65, 1.83, 1.96, 2.14, 2.34, 2.5, 2.68, 2.93, 3.07, 3.57, 3.95, 4.28,
  4.42, 4.61, 4.78, 4.75, 4.78, 4.59, 4.3, 4.13, 3.88, 3.65, 3.45, 3.16, 2.8,
  2.55, 2.35, 2.14, 2.01, 1.7, 1.48, 1.35, 1.21, 1.04, 0.89, 0.63, 0.4, 0.23,
  0.127, 0.065]
let xatt: array[13, float32] = [
  1.0, 2.5, 4.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 100.0, 10000.0]
let yatt: array[13, float32] = [
  0.0, 0.04, 0.17, 0.17, 0.07, 0.07, 0.56, 2.26, 1.3, 0.67, 0.01, 0.0001,
  1e-6]

proc fluoroform*(inp: GasInput, g: var GasXsec) =
  g.describe(inp, " CHF3 -2001---", " CHF3 -2001---", mass = 70.0138)
  for k in 1 .. 6: g.kel[k] = 0
  g.scrp[2] = " ELASTIC         CHF3"
  g.e[3] = 13.86'f32
  g.shells [(0.0'f32, 13.86'f32, 0, 0.0'f32, " IONISATION    ELOSS= 13.86")]
  g.attachment [" ATTACHMENT"]
  g.levels(inp, [
    (-0.035'f32, " ROT           ELOSS= -0.035"),
    (0.035, " ROT           ELOSS=  0.035"),
    (-0.063, " VIB63         ELOSS= -0.063"),
    (0.063, " VIB63         ELOSS=  0.063"),
    (0.1429, " VIB524        ELOSS=  0.1429"),
    (0.2858, " VIB1          ELOSS=  0.2858"),
    (0.3764, " VIB HAR       ELOSS=  0.3764"),
    (0.4287, " VIB HAR       ELOSS=  0.4287"),
    (10.7, " EXC (DISOCTN)  ELOSS= 10.7")])
  for j in 1 .. g.nin: g.kin[j] = 0

  # dipole constants: rotation, V6/V3, V5/V2/V4, V1
  let amp = [float64(0.100'f32), float64(0.152'f32), float64(0.748'f32),
             float64(0.421'f32)]
  let apop = [exp(g.ei[1] / inp.kT), exp(g.ei[3] / inp.kT)]
  let levels = [(6, curve(xvib4, yvib4)), (8, curve(xvib6, yvib6)),
                (9, curve(xexc, yexc))]
  for i, en in inp.steps:
    g.q[2, i] = if en == 0.0: 16300e-16 else: logLog(xen, yxsec, en) * 1e-16
    g.q[3, i] = if en < g.e[3]: 0.0 else: lin(xion, yion, en) * 1e-16
    g.q[4, i] = 0.0
    if en >= float64(xatt[0]) and en <= float64(xatt[^1]):
      g.q[4, i] = lin(xatt, yatt, en) * 1e-19
    g.qatt[1, i] = g.q[4, i]
    g.q[5, i] = 0.0
    g.q[6, i] = 0.0
    # rotation and V6/V3: de-excitation and excitation
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
      g.qin[5, i] = vibrationalDipole(amp[2], g.ei[5], en)
      g.qin[5, i] = (lin(xvib3, yvib3, en) + g.qin[5, i]) * 1e-16
    g.qin[7, i] = 0.0
    if en > g.ei[7]:
      g.qin[7, i] = vibrationalDipole(amp[3], g.ei[7], en)
      g.qin[7, i] = g.qin[7, i] * 1e-16
    for (j, c) in levels:
      g.qin[j, i] = 0.0
      if en > g.ei[j]: g.qin[j, i] = lin(c, en) * 1e-16
    var q = g.q[2, i] + g.q[3, i] + g.q[4, i]
    for j in 1 .. 9: q = q + g.qin[j, i]
    g.q[1, i] = q
  g.dropLevelsAbove(inp.efinal)
