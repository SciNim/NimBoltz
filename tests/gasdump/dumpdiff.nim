## Decodes the differences between two gasdmp dumps (one case each):
##   dumpdiff [-u] ref.bin test.bin [maxShown]
## prints field name, Fortran indices and both values of differing elements.
## -u: only the outputs that Magboltz uses (see `used`); exit code 1 if any
## of them differ.

import std/[os, strutils, strformat, math]

type Field = tuple[name: string, kind: char, dims: seq[int]]  # kind: f, i, c (c: dims[0] = length)

const layout: seq[Field] = @[
  ("Q", 'f', @[6, 4000]), ("QIN", 'f', @[250, 4000]), ("NIN", 'i', @[]),
  ("E", 'f', @[6]), ("EI", 'f', @[250]), ("NAME", 'c', @[25]), ("VIRL", 'f', @[]),
  ("EB", 'f', @[30]), ("PEQEL", 'f', @[6, 4000]), ("PEQIN", 'f', @[250, 4000]),
  ("PENFRA", 'f', @[3, 250]), ("KEL", 'i', @[6]), ("KIN", 'i', @[250]),
  ("QION", 'f', @[30, 4000]), ("PEQION", 'f', @[30, 4000]), ("EION", 'f', @[30]),
  ("NION", 'i', @[]), ("QATT", 'f', @[8, 4000]), ("NATT", 'i', @[]),
  ("QNULL", 'f', @[10, 4000]), ("NNULL", 'i', @[]), ("SCLN", 'f', @[10]),
  ("NC0", 'i', @[30]), ("EC0", 'f', @[30]), ("WK", 'f', @[30]), ("EFL", 'f', @[30]),
  ("NG1", 'i', @[30]), ("EG1", 'f', @[30]), ("NG2", 'i', @[30]), ("EG2", 'f', @[30]),
  ("SCRPT", 'c', @[50, 300]), ("SCRPTN", 'c', @[50, 10]),
  ("EG", 'f', @[4000]), ("EROOT", 'f', @[4000]), ("QT1", 'f', @[4000]),
  ("QT2", 'f', @[4000]), ("QT3", 'f', @[4000]), ("QT4", 'f', @[4000]),
  ("DEN", 'f', @[4000]), ("NANISO", 'i', @[]), ("ESTEP", 'f', @[]), ("AKT", 'f', @[])]

proc idx(dims: seq[int], k: int): string =
  ## column-major linear index -> 1-based Fortran subscripts
  var k = k
  var parts: seq[string]
  for d in dims:
    parts.add $(k mod d + 1)
    k = k div d
  "(" & parts.join(",") & ")"

var args = commandLineParams()
let onlyUsed = args.len > 0 and args[0] == "-u"
if onlyUsed: args.delete(0)
let a = readFile(args[0])
let b = readFile(args[1])
let maxShown = if args.len >= 3: parseInt(args[2]) else: 20

proc intField(d: string, name: string): int =
  var off = 0
  for f in layout:
    let n = (if f.dims.len == 0: 1 else: prod(f.dims))
    if f.name == name:
      var x: int64
      copyMem(addr x, unsafeAddr d[off], 8)
      return int(x)
    off += n * (if f.kind == 'c': 1 else: 8)

let (nion, natt, nin) = (intField(a, "NION"), intField(a, "NATT"), intField(a, "NIN"))

proc used(name: string, k: int): bool =
  ## element k (column-major) of a field is read by MIXER/MIXERT (mixture.nim)
  let r = k mod 6 + 1           # row of the 6 x 4000 arrays
  case name
  of "Q": r == 2 or (r == 3 and nion <= 1) or (r == 4 and natt <= 1)
  of "PEQEL": r == 2 or (r == 3 and nion <= 1)
  of "E": k + 1 in 2 .. 4
  of "NIN", "NION", "NATT", "NAME", "KEL": true
  of "QIN", "PEQIN": k mod 250 + 1 <= nin
  of "EI", "KIN": k + 1 <= nin
  of "PENFRA": k div 3 + 1 <= nin
  of "EB", "NC0", "EC0": k + 1 <= max(nion, 1)
  of "EION": nion > 1 and k + 1 <= nion
  of "QION", "PEQION": nion > 1 and k mod 30 + 1 <= nion
  of "QATT": natt > 1 and k mod 8 + 1 <= natt
  of "SCRPT":
    let j = k + 1
    j == 2 or (if nion <= 1: j == 3 else: j in 3 .. 2 + nion) or
      (if natt <= 1: j == 3 + nion else: j in 3 + nion .. 2 + nion + natt) or
      j in 5 + nion + natt .. 4 + nion + natt + nin
  else: false
if a.len != b.len: echo &"sizes differ: {a.len} vs {b.len}"
var off = 0
var shown = 0
var perField: seq[(string, int)]
for f in layout:
  let n = (if f.dims.len == 0: 1 else: prod(f.dims))
  let size = (if f.kind == 'c': 1 else: 8)
  var cnt = 0
  if f.kind == 'c':
    let w = f.dims[0]
    for k in 0 ..< n div w:
      let o = off + k * w
      if onlyUsed and not used(f.name, k): continue
      if o + w <= min(a.len, b.len) and a[o ..< o+w] != b[o ..< o+w]:
        inc cnt
        if shown < maxShown:
          inc shown
          echo &"{f.name}({k+1}): '{a[o ..< o+w]}' vs '{b[o ..< o+w]}'"
  else:
    for k in 0 ..< n:
      let o = off + k * 8
      if o + 8 > min(a.len, b.len): break
      if onlyUsed and not used(f.name, k): continue
      if a[o ..< o+8] != b[o ..< o+8]:
        inc cnt
        if shown < maxShown:
          inc shown
          if f.kind == 'f':
            var x, y: float64
            copyMem(addr x, unsafeAddr a[o], 8); copyMem(addr y, unsafeAddr b[o], 8)
            echo &"{f.name}{idx(f.dims, k)}: {x:.17g} vs {y:.17g}"
          else:
            var x, y: int64
            copyMem(addr x, unsafeAddr a[o], 8); copyMem(addr y, unsafeAddr b[o], 8)
            echo &"{f.name}{idx(f.dims, k)}: {x} vs {y}"
  if cnt > 0: perField.add (f.name, cnt)
  off += n * size
if off != a.len: echo &"layout covers {off} bytes, dump has {a.len}"
for (n, c) in perField: echo &"  {n}: {c} elements differ"
if onlyUsed and perField.len > 0: quit 1
