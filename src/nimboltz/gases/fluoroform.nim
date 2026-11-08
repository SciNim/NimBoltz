## Fluoroform: GAS50 of Magboltz 11.21 ("CHF3 -2001---"), an older routine:
## isotropic, ionisation as one channel in Q(3), energies stepped by the
## routine itself.
##
## Elastic (momentum transfer, log-log), ionisation and attachment from
## tables; rotation and the V6/V3 vibration with their superelastic partners
## as Born dipoles (populations at the gas temperature), V5/V2/V4 as a
## dipole plus a table, V1 as a dipole, harmonics and one dissociative
## excitation from tables.

import engine

const xen: array[32, float64] = [
  1e-6, 0.0010000000474974513, 0.009999999776482582, 0.10000000149011612, 1.0,
  1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 10.0, 15.0, 20.0, 30.0, 35.0, 50.0, 60.0,
  75.0, 100.0, 150.0, 200.0, 300.0, 500.0, 700.0, 1000.0, 2000.0, 4000.0,
  10000.0, 20000.0, 100000.0]
const yxsec: array[32, float32] = [
  16300.0, 16300.0, 1630.0, 163.0, 16.3, 11.5, 10.5, 9.5, 9.5, 10.0, 11.0,
  11.0, 10.8, 10.5, 10.0, 9.5, 8.4, 6.5, 5.5, 4.6, 3.54, 2.55, 1.63, 0.96,
  0.54, 0.36, 0.23, 0.1, 0.05, 0.02, 0.01, 0.002]
const xvib3: array[14, float32] = [
  0.1429, 4.0, 4.7, 5.7, 6.7, 7.7, 9.5, 15.0, 20.0, 50.0, 100.0, 1000.0,
  10000.0, 100000.0]
const yvib3: array[14, float32] = [
  0.0, 0.0, 0.109, 0.952, 1.43, 1.22, 0.79, 0.068, 0.023, 0.0001, 0.00001,
  1e-6, 1e-7, 1e-8]
const xvib4: array[16, float32] = [
  0.2858, 1.0, 3.0, 3.7, 4.7, 5.7, 6.7, 7.7, 9.5, 15.0, 20.0, 50.0, 100.0,
  1000.0, 10000.0, 100000.0]
const yvib4: array[16, float32] = [
  0.0, 0.001, 0.005, 0.02, 0.03, 0.23, 0.35, 0.3, 0.2, 0.016, 0.006, 0.0003,
  0.00001, 1e-6, 1e-7, 1e-8]
const xvib6: array[16, float32] = [
  0.4287, 1.0, 3.0, 3.7, 4.7, 5.7, 6.7, 7.7, 9.5, 15.0, 20.0, 50.0, 100.0,
  1000.0, 10000.0, 100000.0]
const yvib6: array[16, float32] = [
  0.0, 0.001, 0.04, 0.08, 0.12, 0.92, 1.4, 1.2, 0.8, 0.064, 0.02, 0.0004,
  0.00004, 4e-6, 4e-7, 4e-8]
const xexc: array[33, float32] = [
  10.7, 11.0, 12.0, 13.0, 14.0, 16.0, 18.0, 23.0, 28.0, 33.0, 38.0, 43.0,
  48.0, 60.0, 70.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 250.0,
  300.0, 400.0, 500.0, 600.0, 1000.0, 2000.0, 4000.0, 10000.0, 20000.0,
  100000.0]
const yexc: array[33, float32] = [
  0.0, 0.05, 0.21, 0.37, 0.57, 0.87, 1.1, 1.42, 1.52, 1.57, 1.6, 1.63, 1.62,
  1.61, 1.6, 1.6, 1.59, 1.57, 1.55, 1.52, 1.49, 1.45, 1.37, 1.27, 1.1, 0.92,
  0.8, 0.53, 0.31, 0.18, 0.09, 0.05, 0.013]
const xion: array[52, float32] = [
  13.86, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0,
  32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0, 46.0, 48.0, 60.0, 70.0, 80.0,
  90.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0,
  500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0,
  2500.0, 3000.0, 5000.0, 10000.0, 20000.0, 40000.0, 100000.0]
const yion: array[52, float32] = [
  0.0, 0.031, 0.072, 0.123, 0.184, 0.266, 0.431, 0.59, 0.843, 1.07, 1.27,
  1.46, 1.65, 1.83, 1.96, 2.14, 2.34, 2.5, 2.68, 2.93, 3.07, 3.57, 3.95, 4.28,
  4.42, 4.61, 4.78, 4.75, 4.78, 4.59, 4.3, 4.13, 3.88, 3.65, 3.45, 3.16, 2.8,
  2.55, 2.35, 2.14, 2.01, 1.7, 1.48, 1.35, 1.21, 1.04, 0.89, 0.63, 0.4, 0.23,
  0.127, 0.065]
const xatt: array[13, float32] = [
  1.0, 2.5, 4.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 100.0, 10000.0]
const yatt: array[13, float32] = [
  0.0, 0.04, 0.17, 0.17, 0.07, 0.07, 0.56, 2.26, 1.3, 0.67, 0.01, 0.0001,
  1e-6]


const fluoroformModel* = GasModel(
  name: (" CHF3 -2001---", " CHF3 -2001---"),
  elasticText: (" ELASTIC         CHF3", " ELASTIC         CHF3"),
  mass: 70.0138, e3: 13.86,
  elastic: Elastic(interp: logLog, xMomentum: widen(xen), momentum: widen(yxsec)),
  ionisation: Ionisation(
    scheme: single, counting: table(xion, yion), countingFrom: 13.86,
    channels: @[IonChannel(eoby: 13.86, text: " IONISATION    ELOSS= 13.86")]),
  attachment: @[
    Attach(text: " ATTACHMENT", fromStart: true,
           sigma: table(xatt, yatt, unit = 1e-19.cm², beyond = zero))],
  levels: @[
    Level(ein: -0.035, text: " ROT           ELOSS= -0.035",
          kind: superelastic, tail: noTail, dip: 0.100, population: 1),
    Level(ein: 0.035, text: " ROT           ELOSS=  0.035",
          dip: 0.100, ground: 1),
    Level(ein: -0.063, text: " VIB63         ELOSS= -0.063",
          kind: superelastic, tail: noTail, dip: 0.152, population: 2),
    Level(ein: 0.063, text: " VIB63         ELOSS=  0.063",
          dip: 0.152, ground: 2),
    Level(ein: 0.1429, text: " VIB524        ELOSS=  0.1429",
          sigma: table(xvib3, yvib3), dip: 0.748),
    Level(ein: 0.2858, text: " VIB1          ELOSS=  0.2858",
          sigma: table(xvib4, yvib4)),
    Level(ein: 0.3764, text: " VIB HAR       ELOSS=  0.3764",
          dip: 0.421),
    Level(ein: 0.4287, text: " VIB HAR       ELOSS=  0.4287",
          sigma: table(xvib6, yvib6)),
    Level(ein: 10.7, text: " EXC (DISOCTN)  ELOSS= 10.7",
          sigma: table(xexc, yexc))],
  isotropic: true, stepped: true,
  populations: @[pop(1, group = 1), pop(3, group = 2)])

gas fluoroform, fluoroformModel
