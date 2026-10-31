## Results of a Magboltz output: every value printed with an error
## ("value ... +- p %", "+/- p PERCENT"), keyed by its label (the line
## without value and error) and occurrence, so that the same result can be
## found in another output.
##
## Not included: the per-process collision frequencies ("DETAILED COLLISION
## FREQUENCIES" up to the next run), whose counting errors are meaningless
## (see statcmp.nim).

import std/[strutils, tables, strformat]

type
  Result* = object
    key*: string        # label (line without value and error) and occurrence
    line*: string
    value*, pct*: float

proc isNumber(t: string): bool =
  ## plain number, Fortran exponent allowed: 12, -1.5, .5, 0.1D+04
  var i = 0
  if i < t.len and t[i] in {'+', '-'}: inc i
  var digits = 0
  while i < t.len and t[i].isDigit: (inc i; inc digits)
  if i < t.len and t[i] == '.':
    inc i
    while i < t.len and t[i].isDigit: (inc i; inc digits)
  if digits == 0: return false
  if i < t.len and t[i] in {'D', 'E'}:
    inc i
    if i < t.len and t[i] in {'+', '-'}: inc i
    if i == t.len or not t[i].isDigit: return false
    while i < t.len and t[i].isDigit: inc i
  i == t.len

proc parseNum(s: string): float = parseFloat(s.replace('D', 'E'))

proc findError(line: string, start: int): tuple[first, last: int, pct: string] =
  ## next "+- p %" or "+/- p PERCENT" at or after start; first = -1 if none
  var i = start
  while true:
    var j = line.find("+-", i)
    var k = line.find("+/-", i)
    if j < 0 and k < 0: return (-1, -1, "")
    var p: int
    if j >= 0 and (k < 0 or j < k): (p = j + 2; k = j)
    else: (p = k + 3)
    while p < line.len and line[p] == ' ': inc p
    let q0 = p
    while p < line.len and (line[p].isDigit or line[p] == '.'): inc p
    let pct = line[q0 ..< p]
    while p < line.len and line[p] == ' ': inc p
    if pct.isNumber:
      if line.continuesWith("%", p): return (k, p, pct)
      if line.continuesWith("PERCENT", p): return (k, p + 6, pct)
    i = k + 1

proc lastNumber(s: string): tuple[first, last: int] =
  ## span of the last plain number in s (tokens split at blanks and '=')
  var e = s.len - 1
  while e >= 0:
    while e >= 0 and s[e] in {' ', '='}: dec e
    if e < 0: break
    var b = e
    while b > 0 and s[b - 1] notin {' ', '='}: dec b
    if s[b .. e].isNumber: return (b, e)
    e = b - 1
  (-1, -1)

type Parsed* = object
  results*: seq[Result]
  aborted*: bool

proc parse*(path: string): Parsed =
  var seen: CountTable[string]
  var inFrequencies = false
  for line in lines(path):
    if "PROGRAM MAGBOLTZ" in line: inFrequencies = false
    if "DETAILED COLLISION FREQUENCIES" in line: inFrequencies = true
    if "ENERGY OUT OF RANGE" in line: result.aborted = true
    if inFrequencies: continue
    var pos = 0
    var label = ""
    var found: seq[(float, float)]
    while true:
      let (first, last, pct) = findError(line, pos)
      if first < 0: break
      let seg = line[pos ..< first]
      let (vb, ve) = lastNumber(seg)
      if vb >= 0:
        label.add seg[0 ..< vb]
        label.add seg[ve + 1 .. ^1]
        label.add " ± "
        found.add (parseNum(seg[vb .. ve]), parseNum(pct))
      pos = last + 1
    if found.len == 0: continue
    if pos < line.len: label.add line[pos .. ^1]
    label = label.splitWhitespace.join(" ")
    seen.inc label
    for k, (v, p) in found:
      result.results.add Result(key: &"{label} #{seen[label]} @{k}",
                                line: line.strip, value: v, pct: p)

