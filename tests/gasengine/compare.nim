## Compares the generic gas engine (nimboltz/gases) with the per-gas
## modules it replaces (module `refgases`, extracted from git by run.sh;
## those match Magboltz 11.21 of the `fixes` fork bit for bit, see
## tests/gasdump):
##
##   compare [-t tol] gas...
##
## Each gas runs on the 72 cases of tests/gasdump plus EFINAL 1 and 15 (EFINAL, NANISO,
## temperature, alone or 10 % in argon). Every output is compared, floats
## to a relative tolerance (default 1e-12), integers and strings exactly.
## The anisotropy parameter of a process is skipped where its cross section
## is zero in both (the Monte Carlo never picks the process there; the
## routines leave different values there).
## Prints one line per gas and the largest relative difference; exit code 1
## if any gas differs.

import std/[os, strutils, strformat, math]
import nimboltz/[units, deck, conditions, mixture]
import nimboltz/gases as engine
import refgases as reference
import refxsec

var a: refxsec.GasXsec    # the reference (numbers without units)
var b: mixture.GasXsec    # the engine (large: static storage)

func value(x: SomeFloat): float64 = x
func value(x: SomeUnit): float64 = x.float

var tol = 1e-12
var gases: seq[int]
var args = commandLineParams()
var k = 0
while k < args.len:
  if args[k] == "-t":
    tol = parseFloat(args[k + 1])
    k += 2
  else:
    gases.add parseInt(args[k])
    inc k


type Worst = object
  rel: float64
  where: string
  differ, total: int  # values not bit-identical / compared

proc relDiff(x, y: float64): float64 =
  if x == y: 0.0
  elif x.isNaN or y.isNaN: Inf
  else: abs(x - y) / max(abs(x), abs(y))

proc compare(a: refxsec.GasXsec, b: mixture.GasXsec, worst: var Worst,
             fails: var seq[string]) =
  ## the outputs MIXER and MIXERT read (as tests/gasdump/dumpdiff -u)
  template num(name: string, x, y: float64) =
    let r = relDiff(x, y)
    inc worst.total
    if r > 0.0: inc worst.differ
    if r > worst.rel: worst = Worst(rel: r, where: name)
    if r > tol and fails.len < 8:
      fails.add name & ": " & formatFloat(x, ffDefault, 17) & " vs " &
                formatFloat(y, ffDefault, 17)
  template exact(name: string, x, y: untyped) =
    if x != y and fails.len < 8: fails.add name & ": " & $x & " vs " & $y
  template rows(field: untyped, rs: untyped) =
    for r in rs:
      for c in 1 .. nEnergyBins:
        num(astToStr(field) & "(" & $r & "," & $c & ")", a.field[r, c], b.field[r, c].value)
  template angles(field, xs: untyped, rs: untyped, xrow: int = 0) =
    ## anisotropy parameters where the cross section xs is open
    for r in rs:
      let xr = if xrow > 0: xrow else: r
      for c in 1 .. nEnergyBins:
        if a.xs[xr, c] != 0.0 or b.xs[xr, c].value != 0.0:
          num(astToStr(field) & "(" & $r & "," & $c & ")", a.field[r, c], b.field[r, c])
  template elems(field: untyped, rs: untyped) =
    for r in rs:
      when a.field[r] is float64: num(astToStr(field) & "(" & $r & ")", a.field[r], b.field[r].value)
      else: exact(astToStr(field) & "(" & $r & ")", a.field[r], b.field[r])
  let (nin, nion, natt) = (a.nin, a.nion, a.natt)
  exact("nin", a.nin, b.nin)
  exact("nion", a.nion, b.nion)
  exact("natt", a.natt, b.natt)
  exact("name", a.name, b.name)
  var qrows = @[2]
  if nion <= 1: qrows.add 3
  if natt <= 1: qrows.add 4
  rows(q, qrows)
  angles(peqel, q, [2])
  if nion <= 1: angles(peqel, q, [3])
  if nion > 1:
    rows(qion, 1 .. nion)
    angles(peqion, qion, 1 .. nion)
    elems(eion, 1 .. nion)
  if natt > 1: rows(qatt, 1 .. natt)
  rows(qin, 1 .. nin)
  angles(peqin, qin, 1 .. nin)
  for c in 1 .. nin:
    num("penfra(1," & $c & ")", a.penfra[1, c], b.penningFraction[c])
    num("penfra(2," & $c & ")", a.penfra[2, c], b.penningRange[c].float)
    num("penfra(3," & $c & ")", a.penfra[3, c], b.penningTime[c].float)
  num("e(2)", a.e[2], b.massRatio)
  num("e(3)", a.e[3], b.ionThreshold.float)
  num("e(4)", a.e[4], b.attachThreshold.float)
  elems(kel, 1 .. 6)
  elems(eb, 1 .. max(nion, 1)); elems(nc0, 1 .. max(nion, 1)); elems(ec0, 1 .. max(nion, 1))
  elems(ei, 1 .. nin); elems(kin, 1 .. nin)
  var texts = @[2]
  if nion <= 1: texts.add 3
  else:
    for j in 3 .. 2 + nion: texts.add j
  if natt <= 1: texts.add 3 + nion
  else:
    for j in 3 + nion .. 2 + nion + natt: texts.add j
  for j in 5 + nion + natt .. 4 + nion + natt + nin: texts.add j
  elems(scrp, texts)

var failed = 0
for gas in gases:
  var worst: Worst
  var fails: seq[string]
  var n = 0
  for ef in [0.3, 1.0, 4.0, 15.0, 40.0, 400.0, 8000.0, 200000.0]:
    for an in 0 .. 2:
      for t in [20.0, -50.0]:
        for fr in [100.0, 10.0]:
          inc n
          var d = Deck(ngas: 2, nmax: 1, ithrm: 1, efinal: ef.eV, tempc: t,
                       torr: 760.0.Torr, efield: 1000.0.V•cm⁻¹)
          d.gas[1] = gas
          d.percent[1] = fr
          d.gas[2] = 2
          d.percent[2] = 100.0 - fr
          if fr == 100.0:
            d.ngas = 1
            d.gas[2] = 0
            d.percent[2] = 0.0
          let cond = conditions(d)
          let grid = energyGrid(d.efinal, nEnergyBins)
          let mix = engine.mixtureOf(d.gas, d.percent)
          var inpA = reference.GasInput(naniso: an, kT: cond.kT.float, efinal: ef,
            estep: grid.estep.float, torr: d.torr.float, tempc: t, mixture: mix)
          for i in 1 .. nEnergyBins: inpA.e[i] = grid.e[i].float
          let inpB = engine.GasInput(e: grid.e, naniso: an, kT: cond.kT,
            efinal: ef.eV, estep: grid.estep, torr: d.torr, tempc: t, mixture: mix)
          reset a
          reset b
          reference.crossSections(gas, inpA, a)
          engine.crossSections(gas, inpB, b)
          var f: seq[string]
          compare(a, b, worst, f)
          if f.len > 0 and fails.len == 0:
            fails.add &"case EFINAL={ef} NANISO={an} T={t} share={fr}:"
            fails.add f
  if fails.len == 0:
    echo &"ok    gas{gas} ({n} cases, largest relative difference {worst.rel:.2g} in {worst.where}, {worst.differ} of {worst.total} values not bit-identical)"
  else:
    inc failed
    echo &"FAIL  gas{gas}"
    for l in fails: echo "  ", l
quit(if failed > 0: 1 else: 0)
