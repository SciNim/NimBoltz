## Lists gas-routine outputs that MIXER reads but the gas routine never sets.
##   undef zero.bin poison.bin efinal
## zero.bin / poison.bin: gasdmp dumps of the same case, the driver having
## filled every output with 0 resp. NaN / -1234567 before CALL GASMIX.
## An element MIXER reads (same selection logic as MIXER) that differs between
## the two dumps was left undefined by the gas routine. Angular parameters
## (PEQEL/PEQIN/PEQION) count only where their cross section is nonzero.

import std/[os, strutils, strformat, math, tables]

type Field = tuple[name: string, kind: char, dims: seq[int]]

const layout: seq[Field] = @[
  ("Q", 'f', @[6, 4000]), ("QIN", 'f', @[250, 4000]), ("NIN", 'i', @[]),
  ("E", 'f', @[6]), ("EI", 'f', @[250]), ("NAME", 'c', @[25]), ("VIRL", 'f', @[]),
  ("EB", 'f', @[30]), ("PEQEL", 'f', @[6, 4000]), ("PEQIN", 'f', @[250, 4000]),
  ("PENFRA", 'f', @[3, 250]), ("KEL", 'i', @[6]), ("KIN", 'i', @[250]),
  ("QION", 'f', @[30, 4000]), ("PEQION", 'f', @[30, 4000]), ("EION", 'f', @[30]),
  ("NION", 'i', @[]), ("QATT", 'f', @[8, 4000]), ("NATT", 'i', @[]),
  ("QNULL", 'f', @[10, 4000]), ("NNULL", 'i', @[]), ("SCLN", 'f', @[10]),
  ("NC0", 'i', @[30]), ("EC0", 'f', @[30]), ("WK", 'f', @[30]), ("EFL", 'f', @[30]),
  ("NG1", 'i', @[30]), ("EG1", 'f', @[30]), ("NG2", 'i', @[30]), ("EG2", 'f', @[30])]

let a = readFile(paramStr(1))
let b = readFile(paramStr(2))
let efinal = parseFloat(paramStr(3))
var offs: Table[string, (int, seq[int])]
var off = 0
for f in layout:
  offs[f.name] = (off, f.dims)
  off += (if f.dims.len == 0: 1 else: prod(f.dims)) * (if f.kind == 'c': 1 else: 8)

proc pos(name: string, ix: varargs[int]): int =
  let (o, dims) = offs[name]
  var k = 0
  var m = 1
  for d, i in ix:
    k += (i - 1) * m
    m *= dims[d]
  o + 8 * k

proc geti(s: string, name: string, ix: varargs[int]): int64 =
  copyMem(addr result, unsafeAddr s[pos(name, ix)], 8)
proc getf(s: string, name: string, ix: varargs[int]): float64 =
  copyMem(addr result, unsafeAddr s[pos(name, ix)], 8)

var bad: OrderedTable[string, (int, int, int)]   # key -> (count, first row, last row)
proc chk(name: string, ix: varargs[int]) =
  let p = pos(name, ix)
  if a[p ..< p+8] != b[p ..< p+8]:
    var key = name
    if ix.len == 2: key &= &"({ix[0]},*)"
    elif ix.len == 1: key &= &"({ix[0]})"
    elif ix.len == 3: key &= &"({ix[0]},{ix[1]})"
    let r = (if ix.len == 2: ix[1] else: 0)
    if key in bad:
      let (c, f, _) = bad[key]
      bad[key] = (c + 1, f, r)
    else: bad[key] = (1, r, r)

for n in ["NIN", "NION", "NATT", "NNULL"]:
  if geti(a, n) != geti(b, n): echo &"{n} differs: {geti(a, n)} vs {geti(b, n)}"
let nin = geti(a, "NIN").int
let nion = geti(a, "NION").int
let natt = geti(a, "NATT").int
let nnull = geti(a, "NNULL").int
for k in [2, 3]: chk("KEL", k)
for k in [2, 3, 4]: chk("E", k)
for i in 1..4000:
  chk("Q", 1, i); chk("Q", 2, i)
  if geti(a, "KEL", 2) in [1'i64, 2]: chk("PEQEL", 2, i)
for i in 1..4000:
  if efinal >= getf(a, "E", 3):
    if nion <= 1:
      chk("Q", 3, i)
      if geti(a, "KEL", 3) in [1'i64, 2]: chk("PEQEL", 3, i)
    else:
      for k in 1..nion:
        chk("QION", k, i)
        if geti(a, "KEL", 3) in [1'i64, 2] and getf(a, "QION", k, i) != 0: chk("PEQION", k, i)
  if efinal >= getf(a, "E", 4):
    if natt <= 1: chk("Q", 4, i)
    else:
      for k in 1..natt: chk("QATT", k, i)
  for k in 1..nin:
    chk("QIN", k, i)
    if geti(a, "KIN", k) in [1'i64, 2] and getf(a, "QIN", k, i) != 0: chk("PEQIN", k, i)
  for k in 1..nnull: chk("QNULL", k, i)
if efinal >= getf(a, "E", 3):
  for k in 1..max(nion, 1):
    for n in ["EB", "NC0", "EC0", "NG1", "EG1", "NG2", "EG2", "WK", "EFL"]: chk(n, k)
    if nion > 1: chk("EION", k)
for k in 1..nin:
  chk("KIN", k); chk("EI", k)
  for l in 1..3: chk("PENFRA", l, k)
for k in 1..nnull: chk("SCLN", k)
for key, v in bad:
  if v[1] > 0: echo &"{key} {v[0]} pts, I={v[1]}..{v[2]}"
  else: echo &"{key}"
