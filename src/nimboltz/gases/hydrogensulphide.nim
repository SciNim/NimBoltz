## Hydrogen sulphide: GAS55 of Magboltz 11.21 ("H2S       2003"), an older
## routine: isotropic, ionisation as one channel in Q(3), energies stepped
## by the routine itself.
##
## Elastic (log-log; the two rotations are taken out of it), ionisation and
## attachment from tables; two rotations with their superelastic partners
## as Born dipoles (populations at the gas temperature), V2 and V13 as
## dipoles plus tables, higher vibrations and one excitation from tables.

import engine

const xel: array[25, float32] = [
  0.0, 0.001, 0.01, 0.1, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.75, 2.0, 2.3, 3.0,
  5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 1000.0, 10000.0, 100000.0]
const yel: array[25, float32] = [
  1900.0, 1900.0, 1600.0, 160.0, 40.0, 23.5, 13.0, 8.1, 7.5, 9.0, 15.0, 18.9,
  21.0, 17.0, 21.8, 21.8, 12.8, 8.3, 5.2, 3.3, 1.75, 0.75, 0.07, 0.007,
  0.0007]
const xvib1: array[16, float32] = [
  0.1466, 0.8, 1.0, 1.5, 2.0, 2.3, 2.7, 3.0, 5.0, 7.0, 10.0, 20.0, 100.0,
  1000.0, 10000.0, 100000.0]
const yvib1: array[16, float32] = [
  0.0, 0.0001, 0.02, 0.22, 0.55, 0.6, 0.6, 0.56, 0.33, 0.2, 0.1, 0.05, 0.001,
  0.0001, 0.00001, 1e-6]
const xvib2: array[16, float32] = [
  0.3242, 0.8, 1.0, 1.5, 2.0, 2.3, 2.7, 3.0, 5.0, 7.0, 10.0, 20.0, 100.0,
  1000.0, 10000.0, 100000.0]
const yvib2: array[16, float32] = [
  0.0, 0.0001, 0.1, 1.2, 2.3, 2.15, 1.5, 1.2, 0.3, 0.1, 0.05, 0.025, 0.001,
  0.0001, 0.00001, 1e-6]
const xvib3: array[15, float32] = [
  0.4708, 0.9, 1.0, 1.5, 2.0, 2.3, 2.7, 3.0, 5.0, 7.0, 10.0, 100.0, 1000.0,
  10000.0, 100000.0]
const yvib3: array[15, float32] = [
  0.0, 0.001, 0.006, 0.24, 0.42, 0.5, 0.42, 0.25, 0.05, 0.01, 0.001, 0.0001,
  0.00001, 1e-6, 1e-7]
const xvib4: array[15, float32] = [
  0.6484, 0.9, 1.0, 1.5, 2.0, 2.3, 2.7, 3.0, 5.0, 7.0, 10.0, 100.0, 1000.0,
  10000.0, 100000.0]
const yvib4: array[15, float32] = [
  0.0, 0.001, 0.01, 1.15, 1.4, 1.45, 1.4, 0.8, 0.16, 0.03, 0.001, 0.0001,
  0.00001, 1e-6, 1e-7]
const xexc: array[18, float32] = [
  7.85, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 100.0,
  150.0, 200.0, 300.0, 1000.0, 10000.0, 100000.0]
const yexc: array[18, float32] = [
  0.0, 0.4, 2.0, 3.0, 3.8, 4.2, 4.6, 4.6, 4.2, 3.8, 2.9, 2.3, 1.7, 1.4, 1.15,
  0.55, 0.055, 0.0055]
const xion: array[66, float32] = [
  10.48, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 16.0, 17.0,
  18.0, 19.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0,
  40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0,
  130.0, 150.0, 170.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0,
  600.0, 650.0, 700.0, 750.0, 800.0, 900.0, 1000.0, 1200.0, 1500.0, 2000.0,
  3000.0, 4000.0, 7000.0, 10000.0, 14000.0, 20000.0, 40000.0, 60000.0,
  100000.0]
const yion: array[66, float32] = [
  0.0, 0.088, 0.175, 0.262, 0.349, 0.435, 0.518, 0.621, 0.724, 0.825, 1.02,
  1.22, 1.42, 1.62, 1.79, 2.11, 2.38, 2.62, 2.82, 3.01, 3.17, 3.3, 3.42, 3.52,
  3.61, 3.77, 3.88, 3.94, 3.98, 3.99, 3.99, 3.95, 3.88, 3.8, 3.71, 3.61, 3.52,
  3.33, 3.16, 2.93, 2.61, 2.36, 2.15, 1.98, 1.83, 1.71, 1.6, 1.51, 1.42, 1.35,
  1.28, 1.22, 1.12, 1.04, 0.903, 0.758, 0.603, 0.433, 0.341, 0.212, 0.155,
  0.116, 0.085, 0.046, 0.032, 0.02]
const xatt: array[26, float32] = [
  1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 4.75, 5.0, 5.25, 5.5, 5.65, 6.0,
  6.25, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5]
const yatt: array[26, float32] = [
  0.0, 0.0023, 0.0088, 0.019, 0.018, 0.012, 0.003, 0.0001, 0.0001, 0.0007,
  0.005, 0.01, 0.01, 0.0048, 0.0012, 0.0001, 0.0001, 0.0018, 0.0024, 0.0018,
  0.0018, 0.003, 0.0033, 0.001, 0.0001, 1e-7]


const hydrogenSulphideModel* = GasModel(
  name: ("H2S       2003", "H2S       2003"),
  elasticText: (" ELASTIC       H2S", " ELASTIC       H2S"),
  mass: 34.08088, e3: 10.48,
  elastic: Elastic(interp: logLog, xMomentum: widen(xel), momentum: widen(yel),
                   low: constant, lowBelow: f64(xel, 1), lowQ: f64(yel, 1) * 1e-16,
                   minusLevels: 1 .. 4, minusEach: true),
  ionisation: Ionisation(
    scheme: single, counting: table(xion, yion), countingFrom: 10.48,
    channels: @[IonChannel(eoby: 10.48, text: " IONISATION    ELOSS= 10.48")]),
  attachment: @[
    Attach(text: " ATTACHMENT", fromStart: true,
           sigma: Sigma(x: widen(xatt), y: widen(yatt), unit: 1e-16.cm², beyond: zero))],
  levels: @[
    Level(ein: -0.025, text: " ROT R1        ELOSS= -0.025",
          kind: superelastic, tail: noTail, dip: 0.1875, population: 1),
    Level(ein: 0.025, text: " ROT R1        ELOSS=  0.025",
          dip: 0.1875, ground: 1),
    Level(ein: -0.075, text: " ROT R2        EL0SS= -0.075",
          kind: superelastic, tail: noTail, dip: 0.1725, population: 2),
    Level(ein: 0.075, text: " ROT R2        ELOSS=  0.075",
          dip: 0.1725, ground: 2),
    Level(ein: 0.1466, text: " VIB V2        ELOSS=  0.1466",
          sigma: table(xvib1, yvib1), dip: 0.075),
    Level(ein: 0.3242, text: " VIB V13       ELOSS=  0.3242",
          sigma: table(xvib2, yvib2), dip: 0.375),
    Level(ein: 0.4708, text: " (V13+V2)+HIGH ELOSS=  0.4708",
          sigma: table(xvib3, yvib3)),
    Level(ein: 0.6484, text: " 2V13+HIGH     ELOSS=  0.6484",
          sigma: table(xvib4, yvib4)),
    Level(ein: 7.85, text: " EXC           ELOSS=  7.85",
          sigma: table(xexc, yexc))],
  isotropic: true, stepped: true,
  populations: @[pop(1, group = 1), pop(3, group = 2)],
  firstCut: 5)

proc hydrogenSulphide*(inp: GasInput, g: var GasXsec) = fill(hydrogenSulphideModel, inp, g)
