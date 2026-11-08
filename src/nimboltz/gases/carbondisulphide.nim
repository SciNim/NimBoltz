## Carbon disulphide: GAS39 of Magboltz 11.21 ("CS2  -2001---"), an older
## routine: isotropic, ionisation as one channel in Q(3), energies stepped by
## the routine itself.
##
## Elastic (momentum transfer) and ionisation from tables; three-body (at
## 40 Torr, scaled by 1.3) and two-body attachment; the V2 and V1
## vibrations with their superelastic partners (detailed balance on the same
## tables), V3 as a Born dipole and one excitation.

import engine

const xen: array[30, float32] = [
  0.0, 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5,
  1.8, 2.2, 3.0, 3.5, 5.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0, 40.0, 100.0,
  1000.0, 10000.0, 100000.0]
const yxsec: array[30, float32] = [
  99.0, 90.0, 80.0, 58.0, 45.0, 36.5, 28.5, 23.0, 16.0, 9.82, 7.62, 8.61,
  11.6, 14.3, 20.0, 23.8, 27.6, 28.1, 26.4, 28.1, 29.2, 29.2, 26.4, 17.1, 9.9,
  6.5, 2.7, 0.27, 0.027, 0.0027]
const xvib1: array[39, float32] = [
  0.049, 0.05, 0.055, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.17, 0.2,
  0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 3.5,
  4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 40.0, 100.0, 1000.0, 10000.0,
  100000.0]
const yvib1: array[39, float32] = [
  0.0, 2.4, 10.0, 18.5, 21.0, 22.0, 21.0, 19.5, 14.0, 10.0, 7.0, 5.0, 3.2,
  2.1, 1.5, 1.2, 0.9, 0.78, 0.6, 0.5, 0.43, 0.39, 0.33, 0.29, 0.25, 0.22,
  0.25, 0.32, 0.4, 0.45, 0.2, 0.1, 0.075, 0.052, 0.032, 0.013, 0.0013,
  0.00013, 0.000013]
const xvib2: array[34, float32] = [
  0.081, 0.09, 0.1, 0.11, 0.12, 0.13, 0.15, 0.17, 0.2, 0.25, 0.3, 0.35, 0.4,
  0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0,
  10.0, 20.0, 100.0, 1000.0, 10000.0, 100000.0]
const yvib2: array[34, float32] = [
  0.0, 2.5, 6.0, 7.5, 8.0, 7.5, 5.8, 4.2, 2.5, 1.5, 1.05, 0.74, 0.58, 0.4,
  0.29, 0.16, 0.12, 0.1, 0.1, 0.125, 0.165, 0.27, 0.43, 0.51, 0.49, 0.2, 0.12,
  0.07, 0.057, 0.033, 0.008, 0.0008, 0.00008, 8e-6]
const xexc: array[18, float32] = [
  6.2, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 17.0, 20.0, 30.0, 40.0, 60.0,
  80.0, 100.0, 1000.0, 10000.0, 100000.0]
const yexc: array[18, float32] = [
  0.0, 0.6, 1.5, 3.3, 5.2, 7.0, 8.0, 8.8, 9.2, 8.9, 8.0, 7.4, 6.3, 5.5, 5.0,
  0.5, 0.05, 0.005]
const xion: array[69, float32] = [
  10.07, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5,
  16.0, 16.5, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 26.0, 28.0,
  30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0,
  75.0, 80.0, 90.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0,
  250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 800.0,
  900.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 5000.0, 10000.0, 20000.0,
  40000.0, 100000.0]
const yion: array[69, float32] = [
  0.0, 0.192, 0.421, 0.652, 0.88, 1.1, 1.32, 1.53, 1.72, 1.92, 2.2, 2.48,
  2.74, 3.03, 3.31, 3.84, 4.34, 4.83, 5.28, 5.69, 6.06, 6.4, 6.99, 7.48, 7.9,
  8.26, 8.58, 8.84, 9.05, 9.23, 9.53, 9.69, 9.75, 9.74, 9.68, 9.59, 9.47,
  9.34, 9.06, 8.76, 8.17, 7.63, 7.14, 6.71, 6.33, 5.99, 5.69, 5.55, 4.96,
  4.48, 4.1, 3.79, 3.52, 3.29, 3.09, 2.92, 2.76, 2.5, 2.29, 2.12, 1.56, 1.28,
  1.1, 0.94, 0.66, 0.42, 0.24, 0.134, 0.069]
const xatt: array[30, float32] = [
  2.41, 2.5, 2.6, 2.7, 2.8, 3.0, 3.2, 3.35, 3.6, 3.7, 3.8, 4.0, 4.2, 4.4, 5.4,
  5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.5, 7.75, 8.0, 8.25, 8.5, 10.0,
  100.0, 100000.0]
const yatt: array[30, float32] = [
  0.0, 0.01, 0.02, 0.04, 0.08, 0.4, 2.0, 3.7, 3.0, 3.1, 2.7, 1.5, 0.4, 0.01,
  0.01, 0.1, 0.5, 1.45, 1.8, 0.9, 0.3, 0.2, 0.3, 0.9, 0.5, 0.1, 0.01, 0.01,
  0.001, 0.0001]
const xat1: array[9, float32] = [
  0.0001, 0.001, 0.01, 0.017, 0.025, 0.03, 0.035, 0.04, 10.0]
const yat1: array[9, float32] = [
  35.0, 35.0, 28.0, 20.0, 10.0, 5.0, 1.5, 0.00001, 1e-7]


const carbonDisulphideModel* = GasModel(
  name: (" CS2  -2001---", " CS2  -2001---"),
  elasticText: (" ELASTIC       CS2", " ELASTIC       CS2"),
  mass: 76.1427, e3: 10.07,
  elastic: Elastic(interp: linear, xMomentum: widen(xen), momentum: widen(yxsec)),
  ionisation: Ionisation(
    scheme: single, counting: table(xion, yion), countingFrom: 10.07,
    channels: @[IonChannel(eoby: 10.07, text: " IONISATION    ELOSS= 10.07")]),
  attachment: @[
    Attach(text: " ATTACHMENT 3-BODY AT 40TORR", fromStart: true,
           sigma: Sigma(x: widen(xat1), y: widen(yat1), unit: 1e-16.cm², beyond: zero, scale: 1.3)),
    Attach(text: " ATTACHMENT 2-BODY", fromStart: true,
           sigma: Sigma(x: widen(xatt), y: widen(yatt), unit: 1e-19.cm², beyond: zero))],
  levels: @[
    Level(ein: -0.049, text: " VIB V2        ELOSS= -0.049",
          kind: superelastic, tail: noTail, sigma: table(xvib1, yvib1), population: 1),
    Level(ein: 0.049, text: " VIB V2        ELOSS=  0.049",
          sigma: table(xvib1, yvib1), ground: 1),
    Level(ein: -0.081, text: " VIB V1        ELOSS= -0.081",
          kind: superelastic, tail: noTail, sigma: table(xvib2, yvib2), population: 2),
    Level(ein: 0.081, text: " VIB V1        ELOSS=  0.081",
          sigma: table(xvib2, yvib2), ground: 2),
    Level(ein: 0.190, text: " VIB V3        ELOSS=  0.190",
          dip: 0.710),
    Level(ein: 6.20, text: " EXC           ELOSS=  6.20",
          sigma: table(xexc, yexc))],
  isotropic: true, stepped: true,
  populations: @[pop(1, group = 1), pop(3, group = 2)],
  firstCut: 5)

gas carbonDisulphide, carbonDisulphideModel
