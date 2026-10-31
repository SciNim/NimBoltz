## Collision tables built from synthetic gases (nimboltz/mixture).

import std/[unittest, math]
import nimboltz/[mixture, tables]

proc fill(result: var GasXsec, elastic, ion, att: float64,
          levels: openArray[float64], eion = 15.0) =
  ## constant cross sections; `levels` are energy losses (< 0 superelastic)
  result.e[2] = 1e-4
  result.e[3] = eion
  result.nion = 1
  result.natt = 1
  result.nin = levels.len
  for ie in 1 .. nEnergyBins:
    result.q[2, ie] = elastic
    result.q[3, ie] = ion
    result.q[4, ie] = att
    for j, l in levels: result.qin[j + 1, ie] = 0.1 * elastic
  for j, l in levels: result.ei[j + 1] = l
  result.kel[2] = 2
  for ie in 1 .. nEnergyBins: result.peqel[2, ie] = 0.3
  result.name = "test gas"

var gases = newSeq[GasXsec](2)   # (too big for the stack)
gases[0].fill(1e-15, 2e-16, 0.0, [11.5, -0.1])
gases[1].fill(5e-16, 0.0, 1e-17, [0.3])
let van = @[2.0, 1.0]
let grid = energyGrid(100.0, nEnergyBins)

suite "mixture":
  test "energy grid: bin centres":
    check grid.estep == 0.025
    check grid.e[1] == 0.0125
    check abs(grid.e[nEnergyBins] - (100.0 - 0.0125)) < 1e-9
    check grid.eroot[400] == sqrt(grid.e[400])

  test "angleCut leaves p <= 1 alone":
    check angleCut(0.7) == (1.0, 0.7)
    let (angc, psct) = angleCut(1.2)
    check angc > 0 and angc < 1
    check psct > 0.5

  test "one table for the mixture":
    var tab = new CollisionTable[maxMixProcesses]
    var descr: array[1..maxMixProcesses, string]
    var fcion, fcatt: array[1..nEnergyBins, float64]
    var negative: seq[NegativeFrequency]
    buildMixture(tab[], descr, gases, van, 100.0, grid, fcion, fcatt, negative)
    # gas 1: elastic, ionisation, attachment, 2 levels; gas 2: 4 processes
    check tab.last == 9
    check tab.isize == 8
    check tab.iarry[1 .. 9] == [1, 2, 3, 4, 5, 6, 7, 8, 9]
    check tab.ipn[1 .. 5] == [0, 1, -1, 0, 0]
    check tab.ein[2] == 15.0 / (1.0 + 0.5e-4)
    check tab.index[1] == 2 and tab.index[2] == 0
    let ie = 1234
    let total = (2.0 * (1e-15 + 2e-16 + 2 * 1e-16) + 1.0 * (5e-16 + 1e-17 + 5e-17)) *
                grid.eroot[ie]
    check abs(tab.tcf[ie] / total - 1) < 1e-14
    check abs(tab.cf[ie, tab.last] - 1) < 1e-14
    for i in 2 .. tab.last: check tab.cf[ie, i] >= tab.cf[ie, i - 1]
    check abs(fcion[ie] / (2.0 * 2e-16 * grid.eroot[ie]) - 1) < 1e-14
    check abs(fcatt[ie] / (1.0 * 1e-17 * grid.eroot[ie]) - 1) < 1e-14
    check negative.len == 0
    check anisotropic(gases)

  test "one table per gas (thermal)":
    var tabs = new array[1..maxGases, CollisionTable[maxGasProcesses]]
    var descr = new array[1..maxGases, array[1..maxGasProcesses, string]]
    var fcion, fcatt: array[1..nEnergyBins, float64]
    var negative: seq[NegativeFrequency]
    buildPerGas(tabs[], descr[], gases, van, 100.0, grid, fcion, fcatt, negative)
    check tabs[1].last == 5 and tabs[2].last == 4
    check tabs[2].iarry[1 .. 4] == [1, 2, 3, 4]
    let ie = 77
    # MB-36: the frequencies are summed over the gases
    check abs(fcion[ie] / (2.0 * 2e-16 * grid.eroot[ie]) - 1) < 1e-14
    check abs(fcatt[ie] / (1.0 * 1e-17 * grid.eroot[ie]) - 1) < 1e-14
    let tg = thermalGas(tabs[], gases, 0.025, 1.602176565e-19, 9.10938291e-31)
    check abs(tg.tcfmxg[2] - 1) < 1e-15
    check tg.tcfmax[1] > tg.tcfmax[2]

  test "Penning delocalisation length: microns to metres per axis (MB-34)":
    var g = newSeq[GasXsec](1)
    g[0].fill(1e-15, 0.0, 0.0, [12.0])
    g[0].penfra[2, 1] = 1.0
    var tab = new CollisionTable[maxMixProcesses]
    var descr: array[1..maxMixProcesses, string]
    var fcion, fcatt: array[1..nEnergyBins, float64]
    var negative: seq[NegativeFrequency]
    buildMixture(tab[], descr, g, [1.0], 100.0, grid, fcion, fcatt, negative)
    check abs(tab.penfra[2, 4] - 1e-6 / sqrt(3.0)) < 1e-20
