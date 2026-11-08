## The library interface (nimboltz/api): a setup runs as its deck does,
## and the results are the numbers Magboltz prints.
##   nim c -r -d:release --path:src tests/unit/tapi.nim

import std/[unittest, os, strutils]
import nimboltz/api
import nimboltz/[deck, magboltz]

const arco2 = Setup(mixture: @[(argon, 90.0), (carbonDioxide, 10.0)],
                    efield: 1000.0.V•cm⁻¹, collisions: 4_000_000)

proc printout(s: Setup): (string, Results) =
  ## the results of a setup and what it prints
  let path = getTempDir() / "tapi.out"
  var f = open(path, fmWrite)
  var t = s
  t.output = f
  let r = simulate(t)
  f.close
  (readFile(path), r)

proc printed(text, label: string, field = 0): float =
  ## the number after `label` in the printout
  let line = text[text.find(label) + label.len .. ^1].splitLines[0]
  line.splitWhitespace[field].replace('D', 'E').parseFloat

suite "api":
  test "a setup is a deck":
    let d = arco2.deck
    check d.ngas == 2
    check d.gas[1] == 2 and d.gas[2] == 12
    check d.percent[1] == 90.0 and d.percent[2] == 10.0
    check d.tempc == 20.0
    check d.torr == 760.0.Torr
    check d.ithrm == 1 and d.ipen == 0
    check Setup(mixture: arco2.mixture, bfield: 0.5.T).deck.bmag == 5.0.kGauss

  let (text, r) = printout(arco2)

  test "the printout is that of the deck":
    let path = getTempDir() / "tapi-deck.out"
    var f = open(path, fmWrite)
    var s = initSession(f)
    discard s.run(arco2.deck, arco2.collisions)
    f.close
    check readFile(path) == text

  test "the results are the printed numbers":
    check r.gases.len == 2 and r.gases[0].startsWith("ARGON")
    check abs(r.velocity.z.value.float - text.printed("Z DRIFT VELOCITY =")) <
          1e-4 * r.velocity.z.value.float
    check abs(r.meanEnergy.value.float - text.printed("MEAN ELECTRON ENERGY =")) < 1e-4
    check abs(r.longitudinal.value.float - text.printed("LONGITUDINAL DIFFUSION =")) <
          1e-4 * r.longitudinal.value.float
    # the characteristic energy and the spread after 1 cm, transverse: the
    # two lines after TRANSVERSE DIFFUSION (=  value  EV. +- error %)
    let lines = text[text.find("TRANSVERSE DIFFUSION") .. ^1].splitLines
    let energy = lines[1].splitWhitespace
    let spread = lines[2].splitWhitespace
    let ce = r.characteristicEnergy(r.transverse)
    let sp = r.spread(r.transverse)
    check abs(ce.value.float - energy[1].parseFloat) < 1e-4
    check abs(sp.value.float - spread[1].parseFloat) < 1e-3
    check abs(100.0 * float(sp.error / sp.value) - spread[4].strip(chars = {'%'}).parseFloat) < 0.01
    check r.collisions == 4_000_000
    check r.processes.len > 0 and r.spectrum.len == 4000

  test "without output, the same numbers":
    let q = simulate(arco2)
    check q.velocity.z.value == r.velocity.z.value
    check q.transverse.value == r.transverse.value

  test "another seed, other numbers within errors":
    var s = arco2
    s.seed = 1
    let q = simulate(s)
    check q.velocity.z.value != r.velocity.z.value
    let z = (q.velocity.z.value - r.velocity.z.value) / r.velocity.z.error
    check abs(z.float) < 6.0

  test "an input error raises":
    expect MagboltzError:
      discard simulate(Setup(mixture: @[(argon, 90.0)], efield: 1000.0.V•cm⁻¹))
    expect MagboltzError:
      discard simulate(Setup(efield: 1000.0.V•cm⁻¹))

  test "the analyses run where ionisation matters":
    let q = simulate(Setup(mixture: @[(cf4, 100.0)], efield: 20000.0.V•cm⁻¹,
                           thermal: false, collisions: 4_000_000))
    check q.steadyState.isSome and q.timeOfFlight.isSome and q.friedland.isSome
    check q.ionisation.value == q.steadyState.get.alpha
