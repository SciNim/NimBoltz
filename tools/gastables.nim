## Extracts the DATA tables of a gas routine for src/nimboltz/gases.
##
##   gastables legacy/gas2.nim [names...] > tables.txt
##
## Reads the tables that f2nim wrote (`let x: array[1..n, float64] = [...]`,
## widened from Fortran REAL*4 literals) and prints them as Nim `let` arrays,
## 0-based, ten values a line. A table whose values are all exact single
## precision numbers becomes `float32` (printed in the shortest form that
## reads back to the same float32, i.e. as written in the Fortran DATA
## statement); a table with double precision values (Fortran D exponents)
## becomes `float64`, with its single precision values marked 'f32
## (unless the decimal is exact in both). With names given, only those
## tables are printed, in that order.

import std/[os, strutils, strscans, tables]

type Tab = object
  name: string
  ints: bool
  values: seq[string]

proc parse(path: string): OrderedTable[string, Tab] =
  let src = readFile(path)
  var pos = 0
  while true:
    let start = src.find("\nlet ", pos)
    if start < 0: break
    let eol = src.find('\n', start + 1)
    var name, typ: string
    var n: int
    pos = eol
    if not scanf(src[start + 1 ..< eol], "let $w: array[1..$i, $w] = [", name, n, typ):
      continue
    let close = src.find(']', eol)
    var t = Tab(name: name, ints: typ == "int")
    for v in src[eol + 1 ..< close].split(','):
      let s = v.strip
      if s.len > 0: t.values.add s.strip(chars = {'(', ')'})
    if t.values.len != n:
      quit "table " & t.name & ": " & $t.values.len & " values, expected " & $n
    result[t.name] = t
    pos = close

proc emit(t: Tab): string =
  var single = not t.ints
  var texts: seq[string]
  for s in t.values:
    if t.ints:
      texts.add s
      continue
    let x = parseFloat(s)
    if float64(float32(x)) != x: single = false
  for s in t.values:
    if t.ints: continue
    let x = parseFloat(s)
    texts.add(if single: $float32(x)
              elif parseFloat($float32(x)) == x: $float32(x)
              elif float64(float32(x)) == x: $float32(x) & "'f32"
              else: $x)
  let typ = if t.ints: "int" elif single: "float32" else: "float64"
  if typ == "float64" and texts[0].endsWith("'f32"):
    texts[0] = "float64(" & texts[0] & ")"   # (the first value types the literal)
  result = "let " & t.name & ": array[" & $texts.len & ", " & typ & "] = ["
  var line = "\n  "
  for i, s in texts:
    let item = s & (if i < texts.high: "," else: "]")
    if line.len + item.len + 1 > 80:
      result.add line.strip(leading = false)
      line = "\n  "
    line.add item & " "
  result.add line.strip(leading = false) & "\n"

let tabs = parse(paramStr(1))
if paramCount() > 1:
  for i in 2 .. paramCount():
    let n = paramStr(i)
    if n notin tabs: quit "no table " & n
    stdout.write emit(tabs[n])
else:
  for t in tabs.values: stdout.write emit(t)
