## Statistical comparison of two Magboltz outputs.
##
##   statcmp [--z=5] [--floor=1] [--errors=TSV --deck=NAME] a.out b.out
##
## Every result printed with an error ("value ... +- p %", "+/- p PERCENT")
## is paired with the result of the same label in the other output, and
##   z = (a - b) / sqrt(σa² + σb²),   σ = |value|·sqrt(p² + floor²)/100.
## The floor (percent) covers errors Magboltz underestimates: attachment
## rates of strongly attaching gases (SF6: ±0.08% printed, ~1% real).
## With --errors, p is the larger of the printed error and the spread
## measured over several random streams (spread.nim output, rows of
## deck NAME), so the floor can be small.
##
## Not compared: the per-process collision frequencies ("DETAILED COLLISION
## FREQUENCIES" up to the next run): their errors are counting errors that
## ignore the correlation along an electron's history and are far too small
## for runs with different random numbers.
##
## Prints a summary line (n, Σz²/n, worst |z|), then every pair with
## |z| > zmax ("z=..."), results found in only one output ("unpaired") and
## runs that Magboltz aborted ("aborted"). Exit status: 0 all pairs within
## zmax, 1 a pair beyond zmax, 2 no pair beyond zmax but an aborted run.

import std/[strutils, math, tables, strformat, parseopt]
import magout

proc main =
  var zmax = 5.0
  var floor = 1.0
  var errorsFile, deck = ""
  var files: seq[string]
  for kind, key, val in getopt():
    case kind
    of cmdLongOption, cmdShortOption:
      if key == "z": zmax = parseFloat(val)
      elif key == "floor": floor = parseFloat(val)
      elif key == "errors": errorsFile = val
      elif key == "deck": deck = val
      else: quit("unknown option " & key)
    of cmdArgument: files.add key
    of cmdEnd: discard
  if files.len != 2: quit("usage: statcmp [--z=5] [--floor=1] a.out b.out")
  var a = parse(files[0])
  var b = parse(files[1])
  if errorsFile != "":
    var spread: Table[string, float]
    for line in lines(errorsFile):
      let c = line.split('\t')
      if c.len >= 7 and c[0] == deck:
        let e = parseFloat(c[4])
        if e.classify notin {fcNan, fcInf}: spread[c[1]] = e
    for p in [addr a, addr b]:
      for r in p.results.mitems:
        r.pct = max(r.pct, spread.getOrDefault(r.key))
  var bByKey: Table[string, Result]
  for r in b.results: bByKey[r.key] = r
  var n = 0
  var chi2 = 0.0
  var worst = 0.0
  var worstLine = ""
  var beyond, unpaired: seq[string]
  var paired: CountTable[string]
  for ra in a.results:
    if ra.key notin bByKey:
      unpaired.add &"unpaired ({files[0]}): {ra.line}"; continue
    paired.inc ra.key
    let rb = bByKey[ra.key]
    let s2 = ra.value ^ 2 * (ra.pct ^ 2 + floor ^ 2) / 1e4 +
             rb.value ^ 2 * (rb.pct ^ 2 + floor ^ 2) / 1e4
    var z: float
    if s2 > 0: z = (ra.value - rb.value) / sqrt(s2)
    elif ra.value == rb.value: z = 0
    else: z = Inf
    inc n
    if z.classify != fcInf: chi2 += z * z
    if abs(z) > worst: (worst = abs(z); worstLine = ra.line)
    if abs(z) > zmax:
      beyond.add &"z={z:+.2f}: {ra.line}  |  {rb.line}"
  for rb in b.results:
    if rb.key notin paired: unpaired.add &"unpaired ({files[1]}): {rb.line}"
  echo &"n={n} chi2/n={chi2 / max(n, 1).float:.2f} max|z|={worst:.2f} " &
       &"unpaired={unpaired.len} [{worstLine}]"
  for l in beyond: echo "  ", l
  for l in unpaired: echo "  ", l
  for (p, f) in [(a, files[0]), (b, files[1])]:
    if p.aborted: echo &"  aborted ({f}): ENERGY OUT OF RANGE"
  quit(if beyond.len > 0: 1 elif a.aborted or b.aborted: 2 else: 0)

main()
