## Statistical comparison of two Magboltz outputs.
##
##   statcmp [--z=5] [--floor=1] a.out b.out
##
## Every result printed with an error ("value ... +- p %", "+/- p PERCENT")
## is paired with the result of the same label in the other output, and
##   z = (a - b) / sqrt(σa² + σb²),   σ = |value|·sqrt(p² + floor²)/100.
## The floor (percent) covers errors Magboltz underestimates: attachment
## rates of strongly attaching gases (SF6: ±0.08% printed, ~1% real).
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

type
  Result = object
    key: string        # label (line without value and error) and occurrence
    line: string
    value, pct: float

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

type Parsed = object
  results: seq[Result]
  aborted: bool

proc parse(path: string): Parsed =
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

proc main =
  var zmax = 5.0
  var floor = 1.0
  var files: seq[string]
  for kind, key, val in getopt():
    case kind
    of cmdLongOption, cmdShortOption:
      if key == "z": zmax = parseFloat(val)
      elif key == "floor": floor = parseFloat(val)
      else: quit("unknown option " & key)
    of cmdArgument: files.add key
    of cmdEnd: discard
  if files.len != 2: quit("usage: statcmp [--z=5] [--floor=1] a.out b.out")
  let a = parse(files[0])
  let b = parse(files[1])
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
