## Bromotrifluoromethane: GAS51 of Magboltz 11.21 ("CF3BR -2002--"), an
## older routine: isotropic, ionisation as one channel in Q(3), energies
## stepped by the routine itself.
##
## Elastic (momentum transfer, log-log) and ionisation from tables;
## attachment by a fit below 2 eV; the V1 and V2 vibrations with their
## superelastic partners as Born dipoles (populations at the gas
## temperature), V3 as a dipole, a harmonic and one excitation from tables.

import engine

const xel: array[30, float64] = [
  1e-6, 0.0010000000474974513, 0.009999999776482582, 0.10000000149011612,
  0.4000000059604645, 0.6000000238418579, 0.800000011920929, 1.0,
  1.399999976158142, 2.0, 3.0, 5.0, 7.5, 10.0, 14.0, 20.0, 30.0, 50.0, 75.0,
  100.0, 150.0, 200.0, 400.0, 700.0, 1000.0, 2000.0, 4000.0, 10000.0, 20000.0,
  100000.0]
const yel: array[30, float32] = [
  2800.0, 2800.0, 280.0, 28.0, 7.0, 5.5, 5.0, 5.0, 6.0, 7.0, 9.0, 11.5, 15.5,
  18.0, 18.0, 17.0, 14.0, 10.0, 8.0, 6.0, 4.0, 3.0, 1.45, 0.9, 0.65, 0.3,
  0.15, 0.05, 0.02, 0.004]
const xvibh: array[16, float32] = [
  0.3, 1.0, 3.0, 3.7, 4.7, 5.7, 6.7, 7.7, 9.5, 15.0, 20.0, 50.0, 100.0,
  1000.0, 10000.0, 100000.0]
const yvibh: array[16, float32] = [
  0.0, 0.005, 0.025, 0.1, 0.25, 1.15, 1.75, 1.5, 1.0, 0.08, 0.03, 0.0015,
  0.00005, 5e-6, 5e-7, 5e-8]
const xexc: array[31, float32] = [
  8.7, 9.0, 10.0, 11.0, 12.0, 14.0, 16.0, 21.0, 26.0, 31.0, 36.0, 40.0, 50.0,
  70.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 250.0, 300.0, 400.0, 500.0,
  600.0, 1000.0, 2000.0, 4000.0, 10000.0, 20000.0, 100000.0]
const yexc: array[31, float32] = [
  0.0, 0.07, 0.32, 0.56, 0.85, 1.35, 1.65, 2.15, 2.3, 2.35, 2.4, 2.45, 2.45,
  2.43, 2.41, 2.39, 2.35, 2.3, 2.25, 2.2, 2.05, 1.9, 1.65, 1.4, 1.2, 0.8,
  0.46, 0.27, 0.14, 0.07, 0.02]
const xion: array[40, float32] = [
  11.4, 12.0, 13.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 30.0, 34.0,
  38.0, 42.0, 46.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 175.0,
  200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0,
  2500.0, 3000.0, 5000.0, 10000.0, 20000.0, 100000.0]
const yion: array[40, float32] = [
  0.0, 0.04, 0.13, 0.25, 0.4, 0.66, 1.0, 1.3, 1.7, 2.05, 2.45, 2.85, 3.3, 3.8,
  4.3, 4.6, 5.35, 5.9, 6.4, 6.6, 6.8, 7.0, 7.05, 7.0, 6.9, 6.45, 6.2, 5.5,
  4.75, 4.2, 3.55, 3.05, 2.25, 1.85, 1.6, 1.35, 0.95, 0.6, 0.35, 0.1]


const bromotrifluoromethaneModel* = GasModel(
  name: (" CF3BR -2002--", " CF3BR -2002--"),
  elasticText: (" ELASTIC         CF3BR", " ELASTIC         CF3BR"),
  mass: 148.90991, e3: 11.40,
  elastic: Elastic(interp: logLog, xMomentum: widen(xel), momentum: widen(yel)),
  ionisation: Ionisation(
    scheme: single, counting: table(xion, yion), countingFrom: 11.40,
    channels: @[IonChannel(eoby: 11.40, text: " IONISATION    ELOSS= 11.40")]),
  attachment: @[
    Attach(text: " ATTACHMENT", formula: cf3brAttachment)],
  levels: @[
    Level(ein: -0.050, text: " VIB1          ELOSS= -0.050",
          kind: superelastic, tail: noTail, dip: 0.100, population: 1),
    Level(ein: 0.050, text: " VIB1          ELOSS=  0.050",
          dip: 0.100, ground: 1),
    Level(ein: -0.1345, text: " VIB2          ELOSS= -0.1345",
          kind: superelastic, tail: noTail, dip: 0.110, population: 2),
    Level(ein: 0.1345, text: " VIB2          ELOSS=  0.1345",
          dip: 0.110, ground: 2),
    Level(ein: 0.024999, text: " VIB3          ELOSS=  0.024999",
          dip: 0.750),
    Level(ein: 0.30, text: " VIB HAR       ELOSS=  0.30",
          sigma: table(xvibh, yvibh)),
    Level(ein: 8.7, text: " EXC           ELOSS=  8.7",
          sigma: table(xexc, yexc))],
  isotropic: true, stepped: true,
  populations: @[pop(1, group = 1), pop(3, group = 2)])

gas bromotrifluoromethane, bromotrifluoromethaneModel
