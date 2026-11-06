## Propene and cyclopropane: GAS32 ("PROPENE C3H6 99") and GAS33
## (" CYCLO--C3H6 99") of Magboltz 11.21, two older routines of the same
## shape: isotropic, ionisation as one channel in Q(3) (no QION), energies
## stepped by the routine itself.
##
## Elastic (momentum transfer), ionisation and attachment from linear
## tables; vibrations as pairs of a superelastic and an inelastic level,
## each a Born dipole term amp · ln((1+f)/(1-f)) / E (f = √(1 - EIN/E)) plus
## for most a table, weighted by the population of the excited state;
## harmonics and excitations from tables alone.

import common

let peXen: array[49, float32] = [
  0.0, 0.01, 0.014, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11,
  0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.4, 2.0,
  3.0, 4.0, 5.0, 6.0, 7.5, 8.5, 10.0, 15.0, 20.0, 30.0, 40.0, 70.0, 100.0,
  140.0, 200.0, 250.0, 300.0, 500.0, 1000.0, 1500.0, 10000.0, 100000.0]
let peYxsec: array[49, float32] = [
  54.0, 51.0, 50.0, 49.0, 45.0, 42.0, 39.0, 34.0, 28.5, 22.0, 15.5, 9.4, 6.8,
  4.8, 4.4, 4.8, 6.1, 8.8, 15.5, 19.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.0,
  24.5, 24.0, 25.0, 27.0, 28.0, 30.0, 27.0, 22.0, 15.4, 12.0, 8.31, 6.28,
  3.69, 2.66, 1.57, 0.97, 0.7, 0.57, 0.32, 0.143, 0.092, 0.011, 0.001]
let peXvib1: array[31, float32] = [
  0.114, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 4.0, 5.0,
  6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0, 30.0, 50.0, 100.0, 200.0,
  1000.0, 10000.0, 100000.0]
let peYvib1: array[31, float32] = [
  0.0, 0.001, 0.04, 0.07, 0.14, 0.15, 0.14, 0.1, 0.08, 0.08, 0.08, 0.08, 0.08,
  0.2, 0.28, 0.36, 0.48, 0.64, 0.7, 0.64, 0.59, 0.3, 0.22, 0.17, 0.11, 0.06,
  0.022, 0.008, 0.0008, 0.00003, 3e-6]
let peXvib2: array[31, float32] = [
  0.161, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 4.0, 5.0,
  6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0, 30.0, 50.0, 100.0, 200.0,
  1000.0, 10000.0, 100000.0]
let peYvib2: array[31, float32] = [
  0.0, 0.001, 0.22, 0.75, 3.77, 5.28, 4.65, 3.26, 2.16, 1.21, 0.77, 0.54,
  0.38, 0.42, 0.6, 0.8, 1.11, 1.3, 1.35, 1.2, 1.0, 0.56, 0.44, 0.33, 0.24,
  0.14, 0.06, 0.024, 0.003, 0.0001, 0.00001]
let peXvib3: array[18, float32] = [
  0.322, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 10.0,
  100.0, 1000.0, 10000.0, 100000.0]
let peYvib3: array[18, float32] = [
  0.0, 0.001, 0.15, 0.17, 0.82, 1.5, 1.36, 0.9, 0.52, 0.3, 0.15, 0.08, 0.04,
  0.002, 0.0002, 0.00002, 2e-6, 2e-7]
let peXvib4: array[31, float32] = [
  0.36, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 4.0, 5.0,
  6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0, 30.0, 50.0, 100.0, 200.0,
  1000.0, 10000.0, 100000.0]
let peYvib4: array[31, float32] = [
  0.0, 0.001, 0.052, 0.09, 0.54, 0.86, 0.8, 0.64, 0.46, 0.45, 0.45, 0.45, 0.5,
  0.6, 1.0, 1.4, 1.8, 1.85, 1.7, 1.5, 1.2, 0.65, 0.48, 0.42, 0.28, 0.16, 0.06,
  0.03, 0.004, 0.0001, 0.00001]
let peXvib5: array[21, float32] = [
  0.72, 1.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0,
  30.0, 50.0, 100.0, 200.0, 1000.0, 10000.0, 100000.0]
let peYvib5: array[21, float32] = [
  0.0, 0.0001, 0.008, 0.05, 0.1, 0.16, 0.2, 0.21, 0.18, 0.15, 0.12, 0.06,
  0.05, 0.04, 0.03, 0.015, 0.007, 0.003, 0.0004, 0.00001, 1e-6]
let peXexc1: array[17, float32] = [
  4.18, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0, 30.0,
  100.0, 1000.0, 10000.0, 100000.0]
let peYexc1: array[17, float32] = [
  0.0, 0.11, 0.21, 0.42, 0.84, 0.8, 0.67, 0.61, 0.45, 0.34, 0.27, 0.25, 0.2,
  0.06, 0.006, 0.0006, 0.00006]
let peXexc2: array[23, float32] = [
  7.3, 7.5, 8.0, 8.5, 9.0, 10.0, 11.0, 14.0, 20.0, 25.0, 30.0, 40.0, 60.0,
  80.0, 100.0, 150.0, 200.0, 400.0, 1000.0, 2000.0, 10000.0, 20000.0,
  100000.0]
let peYexc2: array[23, float32] = [
  0.0, 0.026, 0.21, 0.36, 0.65, 1.11, 1.7, 2.38, 2.74, 2.81, 2.86, 2.81, 2.69,
  2.55, 2.38, 2.14, 1.87, 1.46, 0.82, 0.41, 0.09, 0.044, 0.009]
let peXexc3: array[20, float32] = [
  9.0, 10.0, 11.0, 14.0, 16.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 100.0,
  150.0, 200.0, 400.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0]
let peYexc3: array[20, float32] = [
  0.0, 0.43, 1.11, 2.04, 2.3, 2.64, 2.81, 2.86, 2.81, 2.69, 2.55, 2.38, 2.14,
  1.87, 1.46, 0.82, 0.41, 0.09, 0.044, 0.009]
let peXion: array[46, float32] = [
  9.73, 11.0, 12.0, 13.0, 14.0, 16.5, 19.0, 24.0, 30.0, 35.0, 40.0, 45.0,
  50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0,
  300.0, 350.0, 400.0, 450.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0,
  1250.0, 1500.0, 1750.0, 2000.0, 2500.0, 3000.0, 5000.0, 7000.0, 10000.0,
  15000.0, 30000.0, 60000.0, 100000.0]
let peYion: array[46, float32] = [
  0.0, 0.19, 0.43, 0.7, 1.05, 2.12, 3.06, 4.81, 5.97, 6.8, 7.38, 7.88, 8.51,
  9.04, 9.32, 9.42, 9.42, 9.42, 9.14, 8.64, 8.16, 7.71, 7.2, 6.31, 5.77, 5.34,
  4.86, 4.55, 4.0, 3.68, 3.39, 3.02, 2.82, 2.44, 2.1, 1.9, 1.74, 1.5, 1.28,
  0.85, 0.64, 0.47, 0.33, 0.18, 0.097, 0.061]
let peXatt: array[16, float32] = [
  6.85, 7.0, 7.2, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5,
  13.0, 13.2]
let peYatt: array[16, float32] = [
  0.0, 0.67, 1.1, 1.65, 2.8, 4.4, 6.6, 10.3, 14.7, 12.3, 9.7, 6.2, 3.5, 1.3,
  0.5, 0.0]

let cpXen: array[49, float32] = [
  0.0, 0.01, 0.014, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11,
  0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.4, 2.0,
  3.0, 4.0, 5.0, 6.0, 7.5, 8.5, 10.0, 15.0, 20.0, 30.0, 40.0, 70.0, 100.0,
  140.0, 200.0, 250.0, 300.0, 500.0, 1000.0, 1500.0, 10000.0, 100000.0]
let cpYxsec: array[49, float32] = [
  13.0, 11.0, 10.5, 9.8, 7.8, 5.6, 4.2, 2.9, 2.1, 2.0, 2.2, 2.65, 3.25, 3.9,
  5.65, 7.3, 9.15, 10.8, 14.2, 16.8, 20.0, 21.5, 22.0, 22.5, 22.7, 22.8, 22.9,
  23.0, 23.5, 25.5, 27.0, 29.0, 27.0, 22.0, 15.4, 12.0, 8.31, 6.28, 3.69,
  2.66, 1.57, 0.97, 0.7, 0.57, 0.32, 0.143, 0.092, 0.011, 0.001]
let cpXvib1: array[32, float32] = [
  0.107, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0, 4.0, 4.5, 5.0, 5.5,
  6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0, 30.0, 50.0, 100.0,
  200.0, 1000.0, 10000.0, 100000.0]
let cpYvib1: array[32, float32] = [
  0.0, 0.001, 0.022, 0.04, 0.08, 0.08, 0.08, 0.085, 0.085, 0.085, 0.13, 0.22,
  0.7, 1.1, 1.25, 1.15, 0.75, 0.6, 0.71, 0.77, 0.71, 0.64, 0.31, 0.25, 0.18,
  0.12, 0.06, 0.025, 0.01, 0.001, 0.00003, 3e-6]
let cpXvib2: array[31, float32] = [
  0.178, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 4.0, 5.0,
  6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0, 30.0, 50.0, 100.0, 200.0,
  1000.0, 10000.0, 100000.0]
let cpYvib2: array[31, float32] = [
  0.0, 0.001, 0.15, 0.25, 0.4, 0.45, 0.47, 0.5, 0.52, 0.55, 0.57, 0.6, 0.62,
  0.66, 0.74, 0.9, 1.14, 1.33, 1.38, 1.23, 1.01, 0.56, 0.44, 0.34, 0.25, 0.14,
  0.059, 0.025, 0.003, 0.0001, 0.00001]
let cpXvib3: array[15, float32] = [
  0.295, 1.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 10.0, 100.0, 1000.0,
  10000.0, 100000.0]
let cpYvib3: array[15, float32] = [
  0.0, 0.001, 0.01, 0.01, 0.05, 0.1, 0.15, 0.1, 0.05, 0.01, 0.001, 0.0001,
  0.00001, 1e-6, 1e-7]
let cpXvib4: array[28, float32] = [
  0.374, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
  9.0, 10.0, 11.0, 15.0, 20.0, 25.0, 30.0, 50.0, 100.0, 200.0, 1000.0,
  10000.0, 100000.0]
let cpYvib4: array[28, float32] = [
  0.0, 0.001, 0.029, 0.049, 0.3, 0.44, 0.47, 0.5, 0.55, 0.7, 0.75, 1.15, 1.4,
  1.7, 1.8, 1.7, 1.5, 1.4, 0.9, 0.66, 0.57, 0.4, 0.22, 0.092, 0.04, 0.004,
  0.0004, 0.00004]
let cpXvib5: array[21, float32] = [
  0.748, 1.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0,
  30.0, 50.0, 100.0, 200.0, 1000.0, 10000.0, 100000.0]
let cpYvib5: array[21, float32] = [
  0.0, 0.0001, 0.002, 0.03, 0.052, 0.088, 0.11, 0.12, 0.1, 0.084, 0.065,
  0.035, 0.025, 0.02, 0.016, 0.009, 0.004, 0.0014, 0.0002, 5e-6, 5e-7]
let cpXexc1: array[23, float32] = [
  7.3, 7.5, 8.0, 8.5, 9.0, 10.0, 11.0, 14.0, 20.0, 25.0, 30.0, 40.0, 60.0,
  80.0, 100.0, 150.0, 200.0, 400.0, 1000.0, 2000.0, 10000.0, 20000.0,
  100000.0]
let cpYexc1: array[23, float32] = [
  0.0, 0.026, 0.21, 0.36, 0.65, 1.11, 1.7, 2.38, 2.74, 2.81, 2.86, 2.81, 2.69,
  2.55, 2.38, 2.14, 1.87, 1.46, 0.82, 0.41, 0.09, 0.044, 0.009]
let cpXexc2: array[20, float32] = [
  9.0, 10.0, 11.0, 14.0, 16.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 100.0,
  150.0, 200.0, 400.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0]
let cpYexc2: array[20, float32] = [
  0.0, 0.43, 1.11, 2.04, 2.3, 2.64, 2.81, 2.86, 2.81, 2.69, 2.55, 2.38, 2.14,
  1.87, 1.46, 0.82, 0.41, 0.09, 0.044, 0.009]
let cpXion: array[46, float32] = [
  9.86, 11.0, 12.0, 13.0, 14.0, 16.5, 19.0, 24.0, 30.0, 35.0, 40.0, 45.0,
  50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0,
  300.0, 350.0, 400.0, 450.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0,
  1250.0, 1500.0, 1750.0, 2000.0, 2500.0, 3000.0, 5000.0, 7000.0, 10000.0,
  15000.0, 30000.0, 60000.0, 100000.0]
let cpYion: array[46, float32] = [
  0.0, 0.19, 0.43, 0.7, 1.05, 2.12, 3.06, 4.81, 5.97, 6.8, 7.38, 7.88, 8.51,
  9.04, 9.32, 9.42, 9.42, 9.42, 9.14, 8.64, 8.16, 7.71, 7.2, 6.31, 5.77, 5.34,
  4.86, 4.55, 4.0, 3.68, 3.39, 3.02, 2.82, 2.44, 2.1, 1.9, 1.74, 1.5, 1.28,
  0.85, 0.64, 0.47, 0.33, 0.18, 0.097, 0.061]
let cpXatt: array[16, float32] = [
  6.85, 7.0, 7.2, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5,
  13.0, 13.2]
let cpYatt: array[16, float32] = [
  0.0, 0.67, 1.1, 1.65, 2.8, 4.4, 6.6, 10.3, 14.7, 12.3, 9.7, 6.2, 3.5, 1.3,
  0.5, 0.0]

type
  Pair = tuple
    ## superelastic level j, inelastic level j + 1; dipole amplitude;
    ## table (if any) and population of the excited state
    j: int
    amp: float64
    c: Curve
    table: bool
    apop: float64

proc dipoleLog(amp, ei, en: float64; superelastic: bool): float64 {.inline.} =
  let f = sqrt(1.0 - ei / en)
  if superelastic: amp * ln((f + 1.0) / (f - 1.0)) / en
  else: amp * ln((1.0 + f) / (1.0 - f)) / en

proc vibrations(g: var GasXsec, i: int, en: float64, pairs: openArray[Pair]) =
  for p in pairs:
    let (j, k) = (p.j, p.j + 1)
    g.qin[j, i] = 0.0
    if en > 0.0:
      var q = dipoleLog(p.amp, g.ei[j], en, true)
      if p.table:
        let e2 = en + g.ei[k]
        q = q + e2 * lin(p.c, e2) / en
      g.qin[j, i] = q * p.apop / (1.0 + p.apop) * 1e-16
    g.qin[k, i] = 0.0
    if en > g.ei[k]:
      var q = dipoleLog(p.amp, g.ei[k], en, false)
      if p.table: q = q + lin(p.c, en)
      g.qin[k, i] = q / (1.0 + p.apop) * 1e-16

proc common99(g: var GasXsec, i: int, en: float64; xen, yxsec, xion, yion,
              xatt, yatt: openArray[float32]) =
  ## elastic, ionisation, attachment
  g.q[2, i] = lin(xen, yxsec, en) * 1e-16
  g.q[3, i] = 0.0
  if en >= g.e[3]: g.q[3, i] = lin(xion, yion, en) * 1e-16
  g.q[4, i] = 0.0
  if en >= xatt[0] and en <= xatt[^1]: g.q[4, i] = lin(xatt, yatt, en) * 1e-21
  g.qatt[1, i] = g.q[4, i]
  g.q[5, i] = 0.0
  g.q[6, i] = 0.0

proc dipoleTable(g: var GasXsec, i: int, en: float64, j: int, amp: float64,
                 c: Curve) =
  ## an inelastic level with dipole term and table, no population weight
  g.qin[j, i] = 0.0
  if en > g.ei[j]:
    g.qin[j, i] = (dipoleLog(amp, g.ei[j], en, false) + lin(c, en)) * 1e-16

proc table(g: var GasXsec, i: int, en: float64, j: int, c: Curve) =
  g.qin[j, i] = 0.0
  if en > g.ei[j]: g.qin[j, i] = lin(c, en) * 1e-16

proc total(g: var GasXsec, i: int) =
  var q = g.q[2, i] + g.q[3, i] + g.q[4, i]
  for j in 1 .. g.nin: q = q + g.qin[j, i]
  g.q[1, i] = q

proc propene*(inp: GasInput, g: var GasXsec) =
  g.describe(inp, "PROPENE C3H6 99", "PROPENE C3H6 99", mass = 42.08064)
  for k in 1 .. 6: g.kel[k] = 0
  g.scrp[2] = " ELASTIC       PROPENE"
  g.e[3] = 9.73'f32
  g.shells [(0.0'f32, 9.73'f32, 0, 0.0'f32, " IONISATION    ELOSS=  9.73")]
  g.attachment [" ATTACHMENT"]
  g.levels(inp, [
    (-0.0716'f32, " VIB V20       ELOSS= -0.0716"),
    (0.0716, " VIB V20       ELOSS=  0.0716"),
    (-0.114, " VIB V13       ELOSS= -0.114"),
    (0.114, " VIB V13       ELOSS=  0.114"),
    (-0.161, " VIB           ELOSS= -0.161"),
    (0.161, " VIB           ELOSS=  0.161"),
    (0.322, " VIB HAR       ELOSS=  0.322"),
    (0.360, " VIB           ELOSS=  0.360"),
    (0.720, " VIB HAR       ELOSS=  0.720"),
    (4.18, " EXC           ELOSS=  4.18"),
    (7.30, " EXC           ELOSS=  7.30"),
    (9.00, " EXC           ELOSS=  9.00")])
  for j in 1 .. g.nin: g.kin[j] = 0
  let pairs: array[3, Pair] = [
    (1, float64(0.070'f32), curve(peXvib1, peYvib1), false,
     exp(g.ei[1] / inp.kT)),
    (3, float64(0.15'f32), curve(peXvib1, peYvib1), true,
     exp(g.ei[3] / inp.kT)),
    (5, float64(0.15'f32), curve(peXvib2, peYvib2), true,
     exp(g.ei[5] / inp.kT))]
  let amp3 = float64(0.198'f32)
  let v3 = curve(peXvib3, peYvib3)
  let v4 = curve(peXvib4, peYvib4)
  let v5 = curve(peXvib5, peYvib5)
  let x1 = curve(peXexc1, peYexc1)
  let x2 = curve(peXexc2, peYexc2)
  let x3 = curve(peXexc3, peYexc3)
  for i, en in inp.steps:
    g.common99(i, en, peXen, peYxsec, peXion, peYion, peXatt, peYatt)
    g.vibrations(i, en, pairs)
    g.table(i, en, 7, v3)
    g.dipoleTable(i, en, 8, amp3, v4)
    g.table(i, en, 9, v5)
    g.table(i, en, 10, x1)
    g.table(i, en, 11, x2)
    g.table(i, en, 12, x3)
    g.total(i)
  g.dropLevelsAbove(inp.efinal)

proc cyclopropane*(inp: GasInput, g: var GasXsec) =
  g.describe(inp, " CYCLO--C3H6 99", " CYCLO--C3H6 99", mass = 42.08064)
  for k in 1 .. 6: g.kel[k] = 0
  g.scrp[2] = " ELASTIC       CYCLO PROPANE"
  g.e[3] = 9.86'f32
  g.shells [(0.0'f32, 9.86'f32, 0, 0.0'f32, " IONISATION    ELOSS=  9.86")]
  g.attachment [" ATTACHMENT"]
  g.levels(inp, [
    (-0.107'f32, " VIB V11       ELOSS= -0.107"),
    (0.107, " VIB V11       ELOSS=  0.107"),
    (-0.178, " VIB           ELOSS= -0.178"),
    (0.178, " VIB           ELOSS=  0.178"),
    (0.295, " VIB 2V3       ELOSS=  0.295"),
    (0.374, " VIB           ELOSS=  0.374"),
    (0.748, " VIB HAR       ELOSS=  0.748"),
    (7.30, " EXC           ELOSS=  7.30"),
    (9.00, " EXC           ELOSS=  9.00")])
  for j in 1 .. g.nin: g.kin[j] = 0
  let pairs: array[2, Pair] = [
    (1, float64(0.120'f32), curve(cpXvib1, cpYvib1), true,
     exp(g.ei[1] / inp.kT)),
    (3, float64(0.090'f32), curve(cpXvib2, cpYvib2), true,
     exp(g.ei[3] / inp.kT))]
  let amp3 = float64(0.109'f32)
  let v3 = curve(cpXvib3, cpYvib3)
  let v4 = curve(cpXvib4, cpYvib4)
  let v5 = curve(cpXvib5, cpYvib5)
  let x1 = curve(cpXexc1, cpYexc1)
  let x2 = curve(cpXexc2, cpYexc2)
  for i, en in inp.steps:
    g.common99(i, en, cpXen, cpYxsec, cpXion, cpYion, cpXatt, cpYatt)
    g.vibrations(i, en, pairs)
    g.table(i, en, 5, v3)
    g.dipoleTable(i, en, 6, amp3, v4)
    g.table(i, en, 7, v5)
    g.table(i, en, 8, x1)
    g.table(i, en, 9, x2)
    g.total(i)
  g.dropLevelsAbove(inp.efinal)
