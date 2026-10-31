## fixefinal deck.in reference.out > stat.in
##
## Copies a Magboltz deck, replacing every EFINAL = 0 (automatic energy
## range) by twice the range Magboltz chose in the reference output
## ("INTEGRATION FROM 0.0 TO x EV."). The automatic choice depends on the
## random numbers and can be too low for the main run, which then stops
## with ENERGY OUT OF RANGE; a fixed range keeps statistical comparisons
## from failing on that.

import std/[os, strutils]

proc main =
  if paramCount() != 2: quit("usage: fixefinal deck.in reference.out")
  var ranges: seq[float]
  for line in lines(paramStr(2)):
    let i = line.find("INTEGRATION FROM 0.0 TO")
    if i >= 0:
      ranges.add parseFloat(line[i + 23 .. ^1].strip.splitWhitespace[0])
  let deck = readFile(paramStr(1)).strip(leading = false).splitLines
  var run = 0
  var i = 0
  while i < deck.len:
    var f = deck[i].splitWhitespace
    if f.len >= 5 and f[0] != "0":         # first card of a run
      if parseFloat(f[4]) == 0.0:
        if run >= ranges.len: quit("no energy range for run " & $(run + 1) & " in " & paramStr(2))
        f[4] = formatFloat(2 * ranges[run], ffDecimal, 5)
      echo f.join(" ")
      for k in 1 .. 3:
        if i + k < deck.len: echo deck[i + k]
      inc run
      i += 4
    else:
      echo deck[i]
      inc i

main()
