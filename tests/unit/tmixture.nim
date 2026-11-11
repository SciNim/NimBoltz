## Collision tables built from synthetic gases (nimboltz/mixture).

import std/[unittest, math]
import nimboltz/[units, mixture, tables, conditions]

proc fill(result: var GasXsec, elastic, ion, att: cm²,
          levels: openArray[eV], eion = 15.0.eV) =
  ## constant cross sections; `levels` are energy losses (< 0 superelastic)
  result.massRatio = 1e-4
  result.ionThreshold = eion
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
gases[0].fill(1e-15.cm², 2e-16.cm², 0.0.cm², [11.5.eV, -0.1.eV])
gases[1].fill(5e-16.cm², 0.0.cm², 1e-17.cm², [0.3.eV])
let van = @[2.0.cm⁻²•ps⁻¹•eV⁻¹ᐟ², 1.0.cm⁻²•ps⁻¹•eV⁻¹ᐟ²]
let grid = energyGrid(100.0.eV, nEnergyBins)

suite "mixture":
  test "energy grid: bin centres":
    check grid.estep.float == 0.025
    check grid.e[1].float == 0.0125
    check abs(grid.e[nEnergyBins] - (100.0.eV - 0.0125.eV)) < 1e-9.eV
    check grid.eroot[400].float == sqrt(grid.e[400]).float

  test "angleCut leaves p <= 1 alone":
    check angleCut(0.7) == (1.0, 0.7)
    let (angc, psct) = angleCut(1.2)
    check angc > 0 and angc < 1
    check psct > 0.5

  test "one table for the mixture":
    var tab = new CollisionTable[maxMixProcesses]
    var descr: array[1..maxMixProcesses, string]
    var fcion, fcatt: array[1..nEnergyBins, ps⁻¹]
    var negative: seq[NegativeFrequency]
    buildMixture(tab[], descr, gases, van, 100.0.eV, grid, fcion, fcatt, negative)
    # gas 1: elastic, ionisation, attachment, 2 levels; gas 2: 4 processes
    check tab.last == 9
    check tab.isize == 8
    check tab.iarry[1 .. 9] == [1'i32, 2, 3, 4, 5, 6, 7, 8, 9]
    check tab.ipn[1 .. 5] == [0'i32, 1, -1, 0, 0]
    check tab.ein[2].float == 15.0 / (1.0 + 0.5e-4)
    check tab.index[1] == 2 and tab.index[2] == 0
    let ie = 1234
    let total = (2.0 * (1e-15 + 2e-16 + 2 * 1e-16) + 1.0 * (5e-16 + 1e-17 + 5e-17)) *
                grid.eroot[ie].float
    check abs(tab.tcf[ie].float / total - 1) < 1e-14
    check abs(tab.cf[ie, tab.last] - 1) < 1e-14
    for i in 2 .. tab.last: check tab.cf[ie, i] >= tab.cf[ie, i - 1]
    check abs(fcion[ie].float / (2.0 * 2e-16 * grid.eroot[ie].float) - 1) < 1e-14
    check abs(fcatt[ie].float / (1.0 * 1e-17 * grid.eroot[ie].float) - 1) < 1e-14
    check negative.len == 0
    check anisotropic(gases)

  test "moderation levels: only the exact mark (MB-58)":
    var one = newSeq[GasXsec](1)
    one[0].fill(1e-15.cm², 0.0.cm², 0.0.cm²,
                [0.025.eV, float64(0.024999'f32).eV, 0.0249991.eV, 0.1.eV])
    var tab = new CollisionTable[maxMixProcesses]
    var descr: array[1..maxMixProcesses, string]
    var fcion, fcatt: array[1..nEnergyBins, ps⁻¹]
    var negative: seq[NegativeFrequency]
    buildMixture(tab[], descr, one, @[1.0.cm⁻²•ps⁻¹•eV⁻¹ᐟ²], 100.0.eV, grid,
                 fcion, fcatt, negative)
    # elastic, ionisation, attachment, then the 4 levels
    check tab.last == 7
    check tab.moderation[1 .. 7] == [0'i32, 0, 0, 0, 1, 0, 0]

  test "one table per gas (thermal)":
    var tabs = new array[1..maxGases, CollisionTable[maxGasProcesses]]
    var descr = new array[1..maxGases, array[1..maxGasProcesses, string]]
    var fcion, fcatt: array[1..nEnergyBins, ps⁻¹]
    var negative: seq[NegativeFrequency]
    buildPerGas(tabs[], descr[], gases, van, 100.0.eV, grid, fcion, fcatt, negative)
    check tabs[1].last == 5 and tabs[2].last == 4
    check tabs[2].iarry[1 .. 4] == [1'i32, 2, 3, 4]
    let ie = 77
    # MB-36: the frequencies are summed over the gases
    check abs(fcion[ie].float / (2.0 * 2e-16 * grid.eroot[ie].float) - 1) < 1e-14
    check abs(fcatt[ie].float / (1.0 * 1e-17 * grid.eroot[ie].float) - 1) < 1e-14
    let tg = thermalGas(tabs[], gases, 0.025.eV, echarg, emass)
    check abs(tg.tcfmxg[2] - 1) < 1e-15
    check tg.tcfmax[1] > tg.tcfmax[2]
    # sqrt(2kT/M) with M = 2m / 1e-4, from Unchained's (2019) constants
    let mgas = 2.0 * 9.1093837015e-31.kg / 1e-4
    let v = sqrt(2.0 * 0.025.eV.to(J) / mgas).to(m•ps⁻¹)
    check abs(float64(tg.vtmb[1] / v) - 1) < 1e-7

  test "Penning delocalisation length: microns to metres per axis (MB-34)":
    var g = newSeq[GasXsec](1)
    g[0].fill(1e-15.cm², 0.0.cm², 0.0.cm², [12.0.eV])
    g[0].penningRange[1] = 1.0.μm
    var tab = new CollisionTable[maxMixProcesses]
    var descr: array[1..maxMixProcesses, string]
    var fcion, fcatt: array[1..nEnergyBins, ps⁻¹]
    var negative: seq[NegativeFrequency]
    buildMixture(tab[], descr, g, [1.0.cm⁻²•ps⁻¹•eV⁻¹ᐟ²], 100.0.eV, grid, fcion, fcatt,
                 negative)
    check abs(tab.penningRange[4].float - 1e-6 / sqrt(3.0)) < 1e-20
