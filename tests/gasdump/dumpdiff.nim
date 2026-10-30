## Decodes the differences between two gasdmp dumps (one case each):
##   dumpdiff ref.bin test.bin [maxShown]
## prints field name, Fortran indices and both values of differing elements.

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

let a = readFile(paramStr(1))
let b = readFile(paramStr(2))
let maxShown = if paramCount() >= 3: parseInt(paramStr(3)) else: 20
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
      if o + w <= min(a.len, b.len) and a[o ..< o+w] != b[o ..< o+w]:
        inc cnt
        if shown < maxShown:
          inc shown
          echo &"{f.name}({k+1}): '{a[o ..< o+w]}' vs '{b[o ..< o+w]}'"
  else:
    for k in 0 ..< n:
      let o = off + k * 8
      if o + 8 > min(a.len, b.len): break
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
