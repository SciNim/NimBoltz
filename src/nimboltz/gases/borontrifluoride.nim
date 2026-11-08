## Boron trifluoride: GAS42 of Magboltz 11.21 ("BF3  -2001---"), an older
## routine: isotropic elastic, ionisation as one channel in Q(3), energies
## stepped by the routine itself.
##
## Elastic (momentum transfer), ionisation and attachment from tables; the
## V4, V2 and V1 vibrations and their superelastic partners as Born dipoles
## (populations at the gas temperature), V1 and V3 plus tables and
## anisotropic, one harmonic and one excitation from tables.

import engine

const xen: array[62, float32] = [
  0.0, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12,
  0.14, 0.16, 0.18, 0.2, 0.24, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5,
  2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 17.0, 20.0,
  24.0, 28.0, 32.0, 36.0, 40.0, 45.0, 50.0, 60.0, 75.0, 100.0, 150.0, 200.0,
  300.0, 500.0, 700.0, 1000.0, 2000.0, 4000.0, 10000.0, 20000.0, 100000.0]
const yxsec: array[62, float32] = [
  100.0, 80.0, 50.3, 43.0, 39.0, 35.5, 33.0, 31.0, 29.4, 27.8, 26.8, 25.5,
  23.5, 21.0, 18.5, 16.0, 14.0, 10.5, 7.2, 5.65, 4.25, 3.15, 2.7, 2.7, 3.3,
  4.3, 6.2, 7.8, 9.3, 10.4, 11.1, 11.3, 11.3, 10.9, 10.5, 10.0, 9.0, 8.5, 8.0,
  7.5, 7.2, 6.8, 6.5, 6.4, 6.3, 6.2, 6.0, 5.75, 5.05, 4.5, 3.75, 2.7, 1.75,
  1.0, 0.57, 0.38, 0.24, 0.11, 0.05, 0.02, 0.01, 0.002]
const xvibh: array[15, float32] = [
  0.36, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 10.0, 100.0, 1000.0,
  10000.0, 100000.0]
const yvibh: array[15, float32] = [
  0.0, 0.00001, 0.21, 0.54, 0.78, 0.9, 0.78, 0.54, 0.36, 0.18, 0.0024, 0.001,
  0.0001, 0.00001, 1e-6]
const xvib1: array[15, float32] = [
  0.11, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 10.0, 100.0, 1000.0,
  10000.0, 100000.0]
const yvib1: array[15, float32] = [
  0.0, 0.00001, 0.2, 0.48, 0.72, 0.8, 0.72, 0.48, 0.32, 0.12, 0.0016, 0.001,
  0.0001, 0.00001, 1e-6]
const xvib3: array[15, float32] = [
  0.18, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 10.0, 100.0, 1000.0,
  10000.0, 100000.0]
const yvib3: array[15, float32] = [
  0.0, 0.00001, 0.4, 0.96, 1.44, 1.6, 1.44, 0.96, 0.64, 0.24, 0.0032, 0.001,
  0.0001, 0.00001, 1e-6]
const xexc: array[34, float32] = [
  10.0, 10.1, 10.6, 11.1, 11.6, 12.1, 12.6, 13.1, 13.6, 14.1, 14.6, 15.1,
  16.2, 17.2, 18.2, 20.2, 22.2, 24.2, 27.2, 30.3, 40.0, 50.0, 100.0, 200.0,
  300.0, 400.0, 500.0, 600.0, 1000.0, 2000.0, 4000.0, 10000.0, 20000.0,
  100000.0]
const yexc: array[34, float32] = [
  0.0, 0.01, 0.11, 0.21, 0.39, 0.58, 0.65, 0.73, 0.82, 0.89, 0.97, 1.03, 1.15,
  1.24, 1.33, 1.49, 1.61, 1.68, 1.78, 1.82, 1.81, 1.83, 1.88, 1.88, 1.7, 1.4,
  1.1, 0.88, 0.49, 0.22, 0.11, 0.05, 0.03, 0.008]
const xion: array[71, float32] = [
  15.56, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 22.0,
  23.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 45.0, 50.0,
  55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0,
  115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 150.0, 160.0, 170.0, 180.0, 200.0,
  220.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 600.0, 700.0, 800.0, 900.0,
  1000.0, 1200.0, 1400.0, 2000.0, 2500.0, 3000.0, 4000.0, 5000.0, 6000.0,
  8000.0, 10000.0, 20000.0, 100000.0]
const yion: array[71, float32] = [
  0.0, 0.045, 0.064, 0.079, 0.13, 0.183, 0.236, 0.295, 0.356, 0.419, 0.493,
  0.645, 0.8, 0.96, 1.26, 1.54, 1.8, 2.03, 2.25, 2.45, 2.63, 2.79, 3.15, 3.48,
  3.76, 3.99, 4.19, 4.35, 4.48, 4.58, 4.67, 4.74, 4.8, 4.84, 4.88, 4.9, 4.92,
  4.93, 4.93, 4.93, 4.93, 4.92, 4.89, 4.86, 4.81, 4.76, 4.66, 4.54, 4.37,
  4.08, 3.83, 3.59, 3.38, 3.2, 2.88, 2.62, 2.41, 2.23, 2.07, 1.85, 1.66, 1.37,
  1.15, 1.02, 0.82, 0.67, 0.58, 0.45, 0.36, 0.21, 0.06]
const xatt: array[33, float32] = [
  10.0, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4,
  11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6,
  12.7, 12.8, 12.9, 13.0, 20.0, 100.0, 1000.0, 10000.0, 100000.0]
const yatt: array[33, float32] = [
  0.0, 0.0015, 0.0032, 0.0046, 0.0063, 0.0084, 0.01, 0.014, 0.017, 0.02,
  0.022, 0.024, 0.025, 0.025, 0.023, 0.021, 0.018, 0.015, 0.012, 0.0097,
  0.0069, 0.0048, 0.0033, 0.0022, 0.0015, 0.00092, 0.00061, 0.00024, 0.0002,
  0.0001, 0.00001, 1e-6, 1e-7]


const boronTrifluorideModel* = GasModel(
  name: (" BF3  -2001---", " BF3  -2001---"),
  elasticText: (" ELASTIC       BF3", " ELASTIC       BF3"),
  mass: 67.8062, e3: 15.56,
  elastic: Elastic(interp: linear, xMomentum: widen(xen), momentum: widen(yxsec)),
  ionisation: Ionisation(
    scheme: single, counting: table(xion, yion), countingFrom: 15.56,
    channels: @[IonChannel(eoby: 15.56, text: " IONISATION    ELOSS= 15.56")]),
  attachment: @[
    Attach(text: " ATTACHMENT", fromStart: true,
           sigma: Sigma(x: widen(xatt), y: widen(yatt), unit: 1e-16.cm², beyond: zero))],
  levels: @[
    Level(ein: -0.0596, text: " VIB V4        ELOSS= -0.0596",
          kind: superelastic, tail: noTail, dip: 0.018, population: 1),
    Level(ein: -0.086, text: " VIB V2        ELOSS= -0.086",
          kind: superelastic, tail: noTail, dip: 0.045, population: 2),
    Level(ein: -0.110, text: " VIB V1        ELOSS= -0.110",
          kind: superelastic, tail: noTail, dip: 0.37, population: 3),
    Level(ein: 0.0596, text: " VIB V4        ELOSS=  0.0596",
          dip: 0.018, ground: 1),
    Level(ein: 0.086, text: " VIB V2        ELOSS=  0.086",
          dip: 0.045, ground: 2),
    Level(ein: 0.110, text: " VIB V1 (ANIS) ELOSS=  0.110",
          sigma: table(xvib1, yvib1), dip: 0.37, ground: 3, kin: 1, peq: peqBorn),
    Level(ein: 0.180, text: " VIB V3 (ANIS) ELOSS=  0.180",
          sigma: table(xvib3, yvib3), dip: 0.74, kin: 1, peq: peqBorn),
    Level(ein: 0.360, text: " VIB HAR       ELOSS=  0.360",
          sigma: table(xvibh, yvibh)),
    Level(ein: 10.0, text: " EXC           ELOSS= 10.0",
          sigma: table(xexc, yexc))],
  isotropic: true, stepped: true,
  populations: @[pop(1, group = 1), pop(2, group = 2), pop(3, group = 3)],
  momentumRatio: 0.58)

gas boronTrifluoride, boronTrifluorideModel
