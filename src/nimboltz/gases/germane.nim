## Germane: GAS59 of Magboltz 11.21 ("GERMANE 2005"), an older routine:
## ionisation as one channel in Q(3), energies stepped by the routine
## itself.
##
## Elastic momentum transfer (log-log), ionisation and attachment from
## tables; the V2+V4 vibration with its superelastic partner and V1+V3 as
## Born dipoles plus scaled copies of one vibration table (above its start),
## two harmonics as scaled copies of it and one dissociative excitation.

import engine

const xelm: array[69, float32] = [
  0.0, 0.001, 0.01, 0.012, 0.014, 0.017, 0.02, 0.025, 0.03, 0.035, 0.04, 0.05,
  0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.17, 0.2, 0.25, 0.3, 0.35, 0.4,
  0.45, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0,
  6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0,
  100.0, 150.0, 200.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 2000.0,
  4000.0, 6000.0, 8000.0, 10000.0, 20000.0, 100000.0]
const yelm: array[69, float32] = [
  106.0, 106.0, 106.0, 105.0, 104.0, 102.0, 100.0, 95.0, 90.0, 86.0, 83.0,
  76.5, 70.5, 65.5, 60.5, 55.5, 51.0, 43.0, 32.5, 20.5, 13.0, 6.0, 4.0, 3.2,
  3.2, 3.3, 3.6, 4.3, 5.0, 6.5, 9.5, 12.5, 16.5, 21.0, 26.0, 28.5, 30.0, 30.5,
  30.0, 28.5, 26.5, 24.5, 22.5, 20.5, 18.5, 16.0, 13.0, 8.5, 4.7, 3.4, 2.5,
  2.1, 1.55, 1.2, 0.85, 0.65, 0.55, 0.4, 0.35, 0.3, 0.25, 0.2, 0.11, 0.07,
  0.04, 0.03, 0.025, 0.013, 0.002]
const xvibr: array[21, float32] = [
  0.35, 0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,
  6.5, 10.0, 100.0, 1000.0, 10000.0, 100000.0]
const yvibr: array[21, float32] = [
  0.0, 0.16, 0.47, 1.14, 2.25, 3.65, 4.9, 5.4, 4.9, 3.65, 2.25, 1.14, 0.47,
  0.16, 0.04, 0.02, 0.001, 0.0001, 0.00001, 1e-6, 1e-7]
const xion: array[106, float32] = [
  11.33, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5,
  17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5,
  23.0, 23.5, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 45.0,
  50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0,
  110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 160.0, 170.0,
  180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 300.0, 350.0, 400.0,
  450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0,
  1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0,
  1900.0, 2000.0, 2200.0, 2400.0, 2700.0, 3000.0, 4000.0, 5000.0, 6000.0,
  8000.0, 10000.0, 15000.0, 20000.0, 40000.0, 70000.0, 100000.0]
const yion: array[106, float32] = [
  0.0, 0.05, 0.202, 0.359, 0.518, 0.676, 0.832, 0.986, 1.14, 1.28, 1.42, 1.56,
  1.69, 1.82, 1.94, 2.06, 2.17, 2.28, 2.38, 2.49, 2.6, 2.7, 2.8, 2.9, 2.99,
  3.08, 3.16, 3.47, 3.74, 3.96, 4.16, 4.32, 4.46, 4.58, 4.68, 4.88, 5.04,
  5.15, 5.22, 5.26, 5.27, 5.27, 5.26, 5.24, 5.2, 5.16, 5.12, 5.08, 5.03, 4.98,
  4.92, 4.87, 4.82, 4.76, 4.71, 4.66, 4.6, 4.5, 4.4, 4.31, 4.21, 4.12, 4.04,
  3.95, 3.87, 3.8, 3.72, 3.39, 3.12, 2.88, 2.69, 2.52, 2.37, 2.23, 2.12, 2.01,
  1.92, 1.84, 1.76, 1.69, 1.62, 1.56, 1.46, 1.37, 1.29, 1.22, 1.15, 1.1, 1.05,
  1.0, 0.96, 0.922, 0.854, 0.797, 0.724, 0.665, 0.524, 0.435, 0.372, 0.291,
  0.24, 0.169, 0.131, 0.071, 0.043, 0.031]
const xatt: array[10, float32] = [
  7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5]
const yatt: array[10, float32] = [
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
const xdis1: array[28, float32] = [
  7.5, 8.5, 10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 100.0,
  150.0, 200.0, 300.0, 400.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, 3000.0,
  4000.0, 6000.0, 8000.0, 10000.0, 20000.0, 100000.0]
const ydis1: array[28, float32] = [
  0.0, 1.8, 4.5, 5.7, 6.5, 7.0, 7.25, 7.25, 7.0, 6.5, 6.0, 5.5, 4.8, 4.0, 3.0,
  2.5, 2.0, 1.4, 1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.18, 0.16, 0.08, 0.02]


const germaneModel* = GasModel(
  name: ("GERMANE 2005", "GERMANE 2005"),
  elasticText: (" ELASTIC        GERMANE", " ELASTIC        GERMANE"),
  mass: 76.6418, e3: 11.33,
  elastic: Elastic(interp: logLog, xMomentum: widen(xelm), momentum: widen(yelm),
                   low: constant, lowBelow: f64(xelm, 1), lowQ: f64(yelm, 1) * 1e-16),
  ionisation: Ionisation(
    scheme: single, counting: table(xion, yion), countingFrom: 11.33,
    channels: @[IonChannel(eoby: 11.33, text: " IONISATION    ELOSS= 11.33")]),
  attachment: @[
    Attach(text: " ATTACHMENT", fromStart: true,
           sigma: Sigma(x: widen(xatt), y: widen(yatt), unit: 1e-19.cm², beyond: zero))],
  levels: @[
    Level(ein: -0.1016, text: " VIB V2+V4     ELOSS= -0.1016",
          kind: superelastic, tail: noTail, sigma: table(xvibr, yvibr, scale = 1.5, zeroBelow = true), dip: 0.825, population: 1),
    Level(ein: 0.1016, text: " VIB V2+V4     ELOSS=  0.1016",
          sigma: table(xvibr, yvibr, scale = 1.5, zeroBelow = true), dip: 0.825, ground: 1),
    Level(ein: 0.2611, text: " VIB V1+V3     ELOSS=  0.2611",
          sigma: table(xvibr, yvibr, scale = 0.84, zeroBelow = true), dip: 0.530),
    Level(ein: 0.35, text: " VIB  HAR      ELOSS=  0.350",
          sigma: table(xvibr, yvibr, scale = 0.113)),
    Level(ein: 0.50, text: " VIB  HAR      ELOSS=  0.500",
          sigma: table(xvibr, yvibr, scale = 0.074)),
    Level(ein: 7.5, text: " EXC+DISSOC    ELOSS=  7.5",
          sigma: table(xdis1, ydis1))],
  isotropic: true, stepped: true,
  populations: @[pop(1)],
  isotropicPeq: 0.5,
  firstCut: 2)

gas germane, germaneModel
