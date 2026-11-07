## Ozone: GAS53 of Magboltz 11.21 ("OZONE  2002"), an older routine:
## isotropic, ionisation as one channel in Q(3), energies stepped by the
## routine itself.
##
## Elastic (momentum transfer, log-log), ionisation and two dissociative
## attachments (O-, O2-) from tables; rotation and the V2 bend with their
## superelastic partners as Born dipoles (populations at the gas
## temperature), V2 and V3+V1 with tables added, higher vibrations and
## three excitations (Chappuis, Hartley, higher) from tables.

import engine

const xen: array[30, float64] = [
  1e-6, 0.0010000000474974513, 0.009999999776482582, 0.10000000149011612,
  0.20000000298023224, 0.4000000059604645, 0.699999988079071, 1.0, 2.0, 3.0,
  4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0, 300.0,
  600.0, 1000.0, 2000.0, 3000.0, 6000.0, 10000.0, 100000.0]
const yxsec: array[30, float32] = [
  3000.0, 1500.0, 150.0, 19.0, 14.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0,
  12.0, 12.0, 12.0, 12.0, 11.5, 11.0, 10.0, 7.4, 4.7, 1.99, 1.16, 0.41, 0.185,
  0.058, 0.028, 0.0081, 0.0031, 0.0003]
const xion: array[104, float32] = [
  12.75, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0,
  18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0,
  26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 45.0, 50.0, 55.0, 60.0,
  65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0,
  125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0,
  210.0, 220.0, 230.0, 240.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0,
  600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1100.0,
  1250.0, 1500.0, 1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0, 3250.0,
  3500.0, 3750.0, 4000.0, 4500.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0,
  10000.0, 15000.0, 20000.0, 30000.0, 50000.0, 100000.0]
const yion: array[104, float32] = [
  0.0, 0.011, 0.035, 0.059, 0.083, 0.118, 0.162, 0.216, 0.27, 0.324, 0.378,
  0.431, 0.484, 0.535, 0.586, 0.636, 0.685, 0.732, 0.779, 0.824, 0.869, 0.921,
  0.975, 1.03, 1.25, 1.46, 1.66, 1.84, 2.01, 2.17, 2.31, 2.44, 2.72, 2.95,
  3.14, 3.3, 3.42, 3.52, 3.6, 3.66, 3.71, 3.74, 3.77, 3.78, 3.79, 3.8, 3.8,
  3.79, 3.78, 3.77, 3.76, 3.74, 3.72, 3.7, 3.66, 3.61, 3.56, 3.51, 3.46, 3.41,
  3.35, 3.3, 3.25, 3.2, 2.97, 2.76, 2.57, 2.41, 2.27, 2.14, 2.03, 1.93, 1.84,
  1.76, 1.68, 1.62, 1.55, 1.5, 1.44, 1.35, 1.23, 1.07, 0.95, 0.856, 0.78,
  0.717, 0.664, 0.619, 0.58, 0.546, 0.516, 0.489, 0.443, 0.406, 0.348, 0.305,
  0.272, 0.246, 0.225, 0.158, 0.123, 0.086, 0.054, 0.029]
const xvib1: array[18, float32] = [
  0.0869, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0,
  15.0, 100.0, 1000.0, 100000.0]
const yvib1: array[18, float32] = [
  0.0, 0.0005, 0.001, 0.005, 0.01, 0.12, 0.265, 0.24, 0.13, 0.07, 0.03, 0.02,
  0.01, 0.001, 0.0001, 0.00001, 1e-6, 1e-7]
const xvib2: array[18, float32] = [
  0.1292, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0,
  14.0, 100.0, 1000.0, 100000.0]
const yvib2: array[18, float32] = [
  0.0, 0.01, 0.015, 0.02, 0.05, 0.43, 0.83, 0.74, 0.56, 0.44, 0.25, 0.14,
  0.05, 0.01, 0.001, 0.0001, 0.00001, 1e-6]
const xvib3: array[16, float32] = [
  0.2161, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 20.0,
  100.0, 1000.0, 100000.0]
const yvib3: array[16, float32] = [
  0.0, 0.001, 0.01, 0.077, 0.165, 0.14, 0.068, 0.034, 0.022, 0.02, 0.02, 0.01,
  0.005, 0.001, 0.0001, 1e-6]
const xvib4: array[16, float32] = [
  0.266, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 20.0, 100.0,
  1000.0, 100000.0]
const yvib4: array[16, float32] = [
  0.0, 0.003, 0.04, 0.14, 0.278, 0.271, 0.238, 0.21, 0.14, 0.075, 0.029, 0.01,
  0.005, 0.001, 0.0001, 1e-6]
const xvib5: array[16, float32] = [
  0.38, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 20.0, 100.0,
  1000.0, 100000.0]
const yvib5: array[16, float32] = [
  0.0, 0.001, 0.02, 0.07, 0.139, 0.135, 0.119, 0.105, 0.07, 0.038, 0.014,
  0.005, 0.002, 0.001, 0.0001, 1e-6]
const xexc1: array[15, float32] = [
  1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0, 20.0, 40.0, 100.0,
  1000.0, 100000.0]
const yexc1: array[15, float32] = [
  0.0, 0.01, 0.04, 0.12, 0.12, 0.1, 0.08, 0.06, 0.05, 0.04, 0.02, 0.01, 0.004,
  0.0004, 4e-6]
const xexc2: array[14, float32] = [
  4.85, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0, 40.0, 100.0, 1000.0,
  100000.0]
const yexc2: array[14, float32] = [
  0.0, 0.26, 0.63, 0.75, 0.68, 0.65, 0.58, 0.47, 0.37, 0.27, 0.13, 0.05,
  0.005, 5e-6]
const xexc3: array[15, float32] = [
  9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 40.0, 100.0, 150.0, 200.0, 400.0,
  1000.0, 10000.0, 100000.0]
const yexc3: array[15, float32] = [
  0.0, 0.52, 1.5, 1.5, 1.4, 1.4, 1.2, 1.0, 0.8, 0.7, 0.6, 0.3, 0.12, 0.012,
  0.0012]
const xatt: array[25, float32] = [
  0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.6, 3.0, 3.5,
  4.0, 5.0, 6.0, 7.0, 7.5, 8.0, 9.0, 10.0, 20.0, 100000.0]
const yat1: array[25, float32] = [
  0.0, 0.032, 0.07, 0.14, 0.22, 0.29, 0.36, 0.37, 0.36, 0.26, 0.21, 0.12, 0.1,
  0.09, 0.08, 0.06, 0.02, 0.02, 0.05, 0.07, 0.05, 0.01, 0.005, 0.001, 1e-8]
const yat2: array[25, float32] = [
  0.0, 0.0, 0.0, 0.01, 0.08, 0.13, 0.17, 0.15, 0.11, 0.055, 0.025, 0.006,
  0.005, 0.002, 0.002, 0.0015, 0.001, 0.002, 0.003, 0.003, 0.0025, 0.001,
  0.0005, 0.0001, 1e-8]


const ozoneModel* = GasModel(
  name: (" OZONE  2002", " OZONE  2002"),
  elasticText: (" ELASTIC        OZONE", " ELASTIC        OZONE"),
  mass: 47.9982, e3: 12.75,
  elastic: Elastic(interp: logLog, xMomentum: widen(xen), momentum: widen(yxsec)),
  ionisation: Ionisation(
    scheme: single, counting: table(xion, yion), countingFrom: 12.75,
    channels: @[IonChannel(eoby: 12.75, text: " IONISATION    ELOSS= 12.75")]),
  attachment: @[
    Attach(text: " DISOCIATIVE ATTACHMENT  O -", fromStart: true,
           sigma: Sigma(x: widen(xatt), y: widen(yat1), unit: 1e-16.cm², beyond: zero)),
    Attach(text: " DISOCIATIVE ATTACHMENT  O2 -", fromStart: true,
           sigma: Sigma(x: widen(xatt), y: widen(yat2), unit: 1e-16.cm², beyond: zero))],
  levels: @[
    Level(ein: -0.005, text: " ROT           ELOSS= -0.005",
          kind: superelastic, tail: noTail, dip: 0.450, population: 1),
    Level(ein: 0.005, text: " ROT           ELOSS=  0.005",
          dip: 0.450, ground: 1),
    Level(ein: -0.0869, text: " VIB2  BEND    ELOSS= -0.0869",
          kind: superelastic, tail: noTail, dip: 0.0133, population: 2),
    Level(ein: 0.0869, text: " VIB2 BEND     ELOSS=  0.0869",
          sigma: table(xvib1, yvib1), dip: 0.0133, ground: 2),
    Level(ein: 0.1292, text: " VIB3+VIB1     ELOSS=  0.1292",
          sigma: table(xvib2, yvib2), dip: 0.090),
    Level(ein: 0.2161, text: " V12+V23       ELOSS=  0.2161",
          sigma: table(xvib3, yvib3)),
    Level(ein: 0.2660, text: " V13+2V1+2V3   ELOSS=  0.2660",
          sigma: table(xvib4, yvib4)),
    Level(ein: 0.380, text: " SUM HIGH VIB  ELOSS=  0.380",
          sigma: table(xvib5, yvib5)),
    Level(ein: 1.50, text: " EXC CHAPPUIS  ELOSS=  1.50",
          sigma: table(xexc1, yexc1)),
    Level(ein: 4.85, text: " EXC HARTLEY   ELOSS=  4.85",
          sigma: table(xexc2, yexc2)),
    Level(ein: 9.00, text: " EXC           ELOSS=  9.00",
          sigma: table(xexc3, yexc3))],
  isotropic: true, stepped: true,
  populations: @[pop(1, group = 1), pop(3, group = 2)])

proc ozone*(inp: GasInput, g: var GasXsec) = fill(ozoneModel, inp, g)
