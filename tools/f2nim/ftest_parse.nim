import std/[os, strutils]
import fast, fparse
let raws = readStatements(paramStr(1))
echo raws.len, " statements"
let units = parseUnits(raws)
echo units.len, " units"
var n = 0
proc count(l: seq[Stmt]) =
  for s in l:
    inc n
    case s.kind
    of sIfBlock:
      for b in s.bodies: count(b)
    of sDo: count(s.body)
    of sIf: inc n
    else: discard
for u in units: count(u.body)
echo n, " executable statements"
