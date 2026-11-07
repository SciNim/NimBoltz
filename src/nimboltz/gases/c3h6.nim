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

import engine

const peXen: array[49, float32] = [
  0.0, 0.01, 0.014, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11,
  0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.4, 2.0,
  3.0, 4.0, 5.0, 6.0, 7.5, 8.5, 10.0, 15.0, 20.0, 30.0, 40.0, 70.0, 100.0,
  140.0, 200.0, 250.0, 300.0, 500.0, 1000.0, 1500.0, 10000.0, 100000.0]
const peYxsec: array[49, float32] = [
  54.0, 51.0, 50.0, 49.0, 45.0, 42.0, 39.0, 34.0, 28.5, 22.0, 15.5, 9.4, 6.8,
  4.8, 4.4, 4.8, 6.1, 8.8, 15.5, 19.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.0,
  24.5, 24.0, 25.0, 27.0, 28.0, 30.0, 27.0, 22.0, 15.4, 12.0, 8.31, 6.28,
  3.69, 2.66, 1.57, 0.97, 0.7, 0.57, 0.32, 0.143, 0.092, 0.011, 0.001]
const peXvib1: array[31, float32] = [
  0.114, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 4.0, 5.0,
  6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0, 30.0, 50.0, 100.0, 200.0,
  1000.0, 10000.0, 100000.0]
const peYvib1: array[31, float32] = [
  0.0, 0.001, 0.04, 0.07, 0.14, 0.15, 0.14, 0.1, 0.08, 0.08, 0.08, 0.08, 0.08,
  0.2, 0.28, 0.36, 0.48, 0.64, 0.7, 0.64, 0.59, 0.3, 0.22, 0.17, 0.11, 0.06,
  0.022, 0.008, 0.0008, 0.00003, 3e-6]
const peXvib2: array[31, float32] = [
  0.161, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 4.0, 5.0,
  6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0, 30.0, 50.0, 100.0, 200.0,
  1000.0, 10000.0, 100000.0]
const peYvib2: array[31, float32] = [
  0.0, 0.001, 0.22, 0.75, 3.77, 5.28, 4.65, 3.26, 2.16, 1.21, 0.77, 0.54,
  0.38, 0.42, 0.6, 0.8, 1.11, 1.3, 1.35, 1.2, 1.0, 0.56, 0.44, 0.33, 0.24,
  0.14, 0.06, 0.024, 0.003, 0.0001, 0.00001]
const peXvib3: array[18, float32] = [
  0.322, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 10.0,
  100.0, 1000.0, 10000.0, 100000.0]
const peYvib3: array[18, float32] = [
  0.0, 0.001, 0.15, 0.17, 0.82, 1.5, 1.36, 0.9, 0.52, 0.3, 0.15, 0.08, 0.04,
  0.002, 0.0002, 0.00002, 2e-6, 2e-7]
const peXvib4: array[31, float32] = [
  0.36, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 4.0, 5.0,
  6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0, 30.0, 50.0, 100.0, 200.0,
  1000.0, 10000.0, 100000.0]
const peYvib4: array[31, float32] = [
  0.0, 0.001, 0.052, 0.09, 0.54, 0.86, 0.8, 0.64, 0.46, 0.45, 0.45, 0.45, 0.5,
  0.6, 1.0, 1.4, 1.8, 1.85, 1.7, 1.5, 1.2, 0.65, 0.48, 0.42, 0.28, 0.16, 0.06,
  0.03, 0.004, 0.0001, 0.00001]
const peXvib5: array[21, float32] = [
  0.72, 1.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0,
  30.0, 50.0, 100.0, 200.0, 1000.0, 10000.0, 100000.0]
const peYvib5: array[21, float32] = [
  0.0, 0.0001, 0.008, 0.05, 0.1, 0.16, 0.2, 0.21, 0.18, 0.15, 0.12, 0.06,
  0.05, 0.04, 0.03, 0.015, 0.007, 0.003, 0.0004, 0.00001, 1e-6]
const peXexc1: array[17, float32] = [
  4.18, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0, 30.0,
  100.0, 1000.0, 10000.0, 100000.0]
const peYexc1: array[17, float32] = [
  0.0, 0.11, 0.21, 0.42, 0.84, 0.8, 0.67, 0.61, 0.45, 0.34, 0.27, 0.25, 0.2,
  0.06, 0.006, 0.0006, 0.00006]
const peXexc2: array[23, float32] = [
  7.3, 7.5, 8.0, 8.5, 9.0, 10.0, 11.0, 14.0, 20.0, 25.0, 30.0, 40.0, 60.0,
  80.0, 100.0, 150.0, 200.0, 400.0, 1000.0, 2000.0, 10000.0, 20000.0,
  100000.0]
const peYexc2: array[23, float32] = [
  0.0, 0.026, 0.21, 0.36, 0.65, 1.11, 1.7, 2.38, 2.74, 2.81, 2.86, 2.81, 2.69,
  2.55, 2.38, 2.14, 1.87, 1.46, 0.82, 0.41, 0.09, 0.044, 0.009]
const peXexc3: array[20, float32] = [
  9.0, 10.0, 11.0, 14.0, 16.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 100.0,
  150.0, 200.0, 400.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0]
const peYexc3: array[20, float32] = [
  0.0, 0.43, 1.11, 2.04, 2.3, 2.64, 2.81, 2.86, 2.81, 2.69, 2.55, 2.38, 2.14,
  1.87, 1.46, 0.82, 0.41, 0.09, 0.044, 0.009]
const peXion: array[46, float32] = [
  9.73, 11.0, 12.0, 13.0, 14.0, 16.5, 19.0, 24.0, 30.0, 35.0, 40.0, 45.0,
  50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0,
  300.0, 350.0, 400.0, 450.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0,
  1250.0, 1500.0, 1750.0, 2000.0, 2500.0, 3000.0, 5000.0, 7000.0, 10000.0,
  15000.0, 30000.0, 60000.0, 100000.0]
const peYion: array[46, float32] = [
  0.0, 0.19, 0.43, 0.7, 1.05, 2.12, 3.06, 4.81, 5.97, 6.8, 7.38, 7.88, 8.51,
  9.04, 9.32, 9.42, 9.42, 9.42, 9.14, 8.64, 8.16, 7.71, 7.2, 6.31, 5.77, 5.34,
  4.86, 4.55, 4.0, 3.68, 3.39, 3.02, 2.82, 2.44, 2.1, 1.9, 1.74, 1.5, 1.28,
  0.85, 0.64, 0.47, 0.33, 0.18, 0.097, 0.061]
const peXatt: array[16, float32] = [
  6.85, 7.0, 7.2, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5,
  13.0, 13.2]
const peYatt: array[16, float32] = [
  0.0, 0.67, 1.1, 1.65, 2.8, 4.4, 6.6, 10.3, 14.7, 12.3, 9.7, 6.2, 3.5, 1.3,
  0.5, 0.0]

const cpXen: array[49, float32] = [
  0.0, 0.01, 0.014, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11,
  0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.4, 2.0,
  3.0, 4.0, 5.0, 6.0, 7.5, 8.5, 10.0, 15.0, 20.0, 30.0, 40.0, 70.0, 100.0,
  140.0, 200.0, 250.0, 300.0, 500.0, 1000.0, 1500.0, 10000.0, 100000.0]
const cpYxsec: array[49, float32] = [
  13.0, 11.0, 10.5, 9.8, 7.8, 5.6, 4.2, 2.9, 2.1, 2.0, 2.2, 2.65, 3.25, 3.9,
  5.65, 7.3, 9.15, 10.8, 14.2, 16.8, 20.0, 21.5, 22.0, 22.5, 22.7, 22.8, 22.9,
  23.0, 23.5, 25.5, 27.0, 29.0, 27.0, 22.0, 15.4, 12.0, 8.31, 6.28, 3.69,
  2.66, 1.57, 0.97, 0.7, 0.57, 0.32, 0.143, 0.092, 0.011, 0.001]
const cpXvib1: array[32, float32] = [
  0.107, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0, 4.0, 4.5, 5.0, 5.5,
  6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0, 30.0, 50.0, 100.0,
  200.0, 1000.0, 10000.0, 100000.0]
const cpYvib1: array[32, float32] = [
  0.0, 0.001, 0.022, 0.04, 0.08, 0.08, 0.08, 0.085, 0.085, 0.085, 0.13, 0.22,
  0.7, 1.1, 1.25, 1.15, 0.75, 0.6, 0.71, 0.77, 0.71, 0.64, 0.31, 0.25, 0.18,
  0.12, 0.06, 0.025, 0.01, 0.001, 0.00003, 3e-6]
const cpXvib2: array[31, float32] = [
  0.178, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 4.0, 5.0,
  6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0, 30.0, 50.0, 100.0, 200.0,
  1000.0, 10000.0, 100000.0]
const cpYvib2: array[31, float32] = [
  0.0, 0.001, 0.15, 0.25, 0.4, 0.45, 0.47, 0.5, 0.52, 0.55, 0.57, 0.6, 0.62,
  0.66, 0.74, 0.9, 1.14, 1.33, 1.38, 1.23, 1.01, 0.56, 0.44, 0.34, 0.25, 0.14,
  0.059, 0.025, 0.003, 0.0001, 0.00001]
const cpXvib3: array[15, float32] = [
  0.295, 1.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 10.0, 100.0, 1000.0,
  10000.0, 100000.0]
const cpYvib3: array[15, float32] = [
  0.0, 0.001, 0.01, 0.01, 0.05, 0.1, 0.15, 0.1, 0.05, 0.01, 0.001, 0.0001,
  0.00001, 1e-6, 1e-7]
const cpXvib4: array[28, float32] = [
  0.374, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
  9.0, 10.0, 11.0, 15.0, 20.0, 25.0, 30.0, 50.0, 100.0, 200.0, 1000.0,
  10000.0, 100000.0]
const cpYvib4: array[28, float32] = [
  0.0, 0.001, 0.029, 0.049, 0.3, 0.44, 0.47, 0.5, 0.55, 0.7, 0.75, 1.15, 1.4,
  1.7, 1.8, 1.7, 1.5, 1.4, 0.9, 0.66, 0.57, 0.4, 0.22, 0.092, 0.04, 0.004,
  0.0004, 0.00004]
const cpXvib5: array[21, float32] = [
  0.748, 1.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 15.0, 20.0, 25.0,
  30.0, 50.0, 100.0, 200.0, 1000.0, 10000.0, 100000.0]
const cpYvib5: array[21, float32] = [
  0.0, 0.0001, 0.002, 0.03, 0.052, 0.088, 0.11, 0.12, 0.1, 0.084, 0.065,
  0.035, 0.025, 0.02, 0.016, 0.009, 0.004, 0.0014, 0.0002, 5e-6, 5e-7]
const cpXexc1: array[23, float32] = [
  7.3, 7.5, 8.0, 8.5, 9.0, 10.0, 11.0, 14.0, 20.0, 25.0, 30.0, 40.0, 60.0,
  80.0, 100.0, 150.0, 200.0, 400.0, 1000.0, 2000.0, 10000.0, 20000.0,
  100000.0]
const cpYexc1: array[23, float32] = [
  0.0, 0.026, 0.21, 0.36, 0.65, 1.11, 1.7, 2.38, 2.74, 2.81, 2.86, 2.81, 2.69,
  2.55, 2.38, 2.14, 1.87, 1.46, 0.82, 0.41, 0.09, 0.044, 0.009]
const cpXexc2: array[20, float32] = [
  9.0, 10.0, 11.0, 14.0, 16.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 100.0,
  150.0, 200.0, 400.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0]
const cpYexc2: array[20, float32] = [
  0.0, 0.43, 1.11, 2.04, 2.3, 2.64, 2.81, 2.86, 2.81, 2.69, 2.55, 2.38, 2.14,
  1.87, 1.46, 0.82, 0.41, 0.09, 0.044, 0.009]
const cpXion: array[46, float32] = [
  9.86, 11.0, 12.0, 13.0, 14.0, 16.5, 19.0, 24.0, 30.0, 35.0, 40.0, 45.0,
  50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0,
  300.0, 350.0, 400.0, 450.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0,
  1250.0, 1500.0, 1750.0, 2000.0, 2500.0, 3000.0, 5000.0, 7000.0, 10000.0,
  15000.0, 30000.0, 60000.0, 100000.0]
const cpYion: array[46, float32] = [
  0.0, 0.19, 0.43, 0.7, 1.05, 2.12, 3.06, 4.81, 5.97, 6.8, 7.38, 7.88, 8.51,
  9.04, 9.32, 9.42, 9.42, 9.42, 9.14, 8.64, 8.16, 7.71, 7.2, 6.31, 5.77, 5.34,
  4.86, 4.55, 4.0, 3.68, 3.39, 3.02, 2.82, 2.44, 2.1, 1.9, 1.74, 1.5, 1.28,
  0.85, 0.64, 0.47, 0.33, 0.18, 0.097, 0.061]
const cpXatt: array[16, float32] = [
  6.85, 7.0, 7.2, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5,
  13.0, 13.2]
const cpYatt: array[16, float32] = [
  0.0, 0.67, 1.1, 1.65, 2.8, 4.4, 6.6, 10.3, 14.7, 12.3, 9.7, 6.2, 3.5, 1.3,
  0.5, 0.0]


const propeneModel* = GasModel(
  name: ("PROPENE C3H6 99", "PROPENE C3H6 99"),
  elasticText: (" ELASTIC       PROPENE", " ELASTIC       PROPENE"),
  mass: 42.08064, e3: 9.73,
  elastic: Elastic(interp: linear, xMomentum: widen(peXen), momentum: widen(peYxsec)),
  ionisation: Ionisation(
    scheme: single, counting: table(peXion, peYion), countingFrom: 9.73,
    channels: @[IonChannel(eoby: 9.73, text: " IONISATION    ELOSS=  9.73")]),
  attachment: @[
    Attach(text: " ATTACHMENT", fromStart: true,
           sigma: Sigma(x: widen(peXatt), y: widen(peYatt), unit: 1e-21.cm², beyond: zero))],
  levels: @[
    Level(ein: -0.0716, text: " VIB V20       ELOSS= -0.0716",
          kind: superelastic, tail: noTail, dip: 0.070, population: 1),
    Level(ein: 0.0716, text: " VIB V20       ELOSS=  0.0716",
          dip: 0.070, ground: 1),
    Level(ein: -0.114, text: " VIB V13       ELOSS= -0.114",
          kind: superelastic, tail: noTail, sigma: table(peXvib1, peYvib1), dip: 0.15, population: 2),
    Level(ein: 0.114, text: " VIB V13       ELOSS=  0.114",
          sigma: table(peXvib1, peYvib1), dip: 0.15, ground: 2),
    Level(ein: -0.161, text: " VIB           ELOSS= -0.161",
          kind: superelastic, tail: noTail, sigma: table(peXvib2, peYvib2), dip: 0.15, population: 3),
    Level(ein: 0.161, text: " VIB           ELOSS=  0.161",
          sigma: table(peXvib2, peYvib2), dip: 0.15, ground: 3),
    Level(ein: 0.322, text: " VIB HAR       ELOSS=  0.322",
          sigma: table(peXvib3, peYvib3)),
    Level(ein: 0.360, text: " VIB           ELOSS=  0.360",
          sigma: table(peXvib4, peYvib4), dip: 0.198),
    Level(ein: 0.720, text: " VIB HAR       ELOSS=  0.720",
          sigma: table(peXvib5, peYvib5)),
    Level(ein: 4.18, text: " EXC           ELOSS=  4.18",
          sigma: table(peXexc1, peYexc1)),
    Level(ein: 7.30, text: " EXC           ELOSS=  7.30",
          sigma: table(peXexc2, peYexc2)),
    Level(ein: 9.00, text: " EXC           ELOSS=  9.00",
          sigma: table(peXexc3, peYexc3))],
  isotropic: true, stepped: true,
  populations: @[pop(1, group = 1), pop(3, group = 2), pop(5, group = 3)])

const cyclopropaneModel* = GasModel(
  name: (" CYCLO--C3H6 99", " CYCLO--C3H6 99"),
  elasticText: (" ELASTIC       CYCLO PROPANE", " ELASTIC       CYCLO PROPANE"),
  mass: 42.08064, e3: 9.86,
  elastic: Elastic(interp: linear, xMomentum: widen(cpXen), momentum: widen(cpYxsec)),
  ionisation: Ionisation(
    scheme: single, counting: table(cpXion, cpYion), countingFrom: 9.86,
    channels: @[IonChannel(eoby: 9.86, text: " IONISATION    ELOSS=  9.86")]),
  attachment: @[
    Attach(text: " ATTACHMENT", fromStart: true,
           sigma: Sigma(x: widen(cpXatt), y: widen(cpYatt), unit: 1e-21.cm², beyond: zero))],
  levels: @[
    Level(ein: -0.107, text: " VIB V11       ELOSS= -0.107",
          kind: superelastic, tail: noTail, sigma: table(cpXvib1, cpYvib1), dip: 0.120, population: 1),
    Level(ein: 0.107, text: " VIB V11       ELOSS=  0.107",
          sigma: table(cpXvib1, cpYvib1), dip: 0.120, ground: 1),
    Level(ein: -0.178, text: " VIB           ELOSS= -0.178",
          kind: superelastic, tail: noTail, sigma: table(cpXvib2, cpYvib2), dip: 0.090, population: 2),
    Level(ein: 0.178, text: " VIB           ELOSS=  0.178",
          sigma: table(cpXvib2, cpYvib2), dip: 0.090, ground: 2),
    Level(ein: 0.295, text: " VIB 2V3       ELOSS=  0.295",
          sigma: table(cpXvib3, cpYvib3)),
    Level(ein: 0.374, text: " VIB           ELOSS=  0.374",
          sigma: table(cpXvib4, cpYvib4), dip: 0.109),
    Level(ein: 0.748, text: " VIB HAR       ELOSS=  0.748",
          sigma: table(cpXvib5, cpYvib5)),
    Level(ein: 7.30, text: " EXC           ELOSS=  7.30",
          sigma: table(cpXexc1, cpYexc1)),
    Level(ein: 9.00, text: " EXC           ELOSS=  9.00",
          sigma: table(cpXexc2, cpYexc2))],
  isotropic: true, stepped: true,
  populations: @[pop(1, group = 1), pop(3, group = 2)])

proc propene*(inp: GasInput, g: var GasXsec) = fill(propeneModel, inp, g)

proc cyclopropane*(inp: GasInput, g: var GasXsec) = fill(cyclopropaneModel, inp, g)
