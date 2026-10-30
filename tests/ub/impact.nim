## Summarises how two gasdmp dumps (e.g. upstream vs fixed Magboltz) differ:
## one line per array row (cross-section level) with the number of energy
## points that differ, their energy range and the largest relative change.
##   impact a.bin b.bin
import std/[os, strformat, math]
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
  ("NG1", 'i', @[30]), ("EG1", 'f', @[30]), ("NG2", 'i', @[30]), ("EG2", 'f', @[30]),
  ("SCRPT", 'c', @[50, 300]), ("SCRPTN", 'c', @[50, 10]),
  ("EG", 'f', @[4000]), ("EROOT", 'f', @[4000]), ("QT1", 'f', @[4000]),
  ("QT2", 'f', @[4000]), ("QT3", 'f', @[4000]), ("QT4", 'f', @[4000]),
  ("DEN", 'f', @[4000]), ("NANISO", 'i', @[]), ("ESTEP", 'f', @[]), ("AKT", 'f', @[])]
let a = readFile(paramStr(1))
let b = readFile(paramStr(2))
proc f64(s: string, o: int): float64 = copyMem(addr result, unsafeAddr s[o], 8)
var offs: seq[int]
var off = 0
var egOff = -1
for f in layout:
  offs.add off
  if f.name == "EG": egOff = off
  let n = (if f.dims.len == 0: 1 else: prod(f.dims))
  off += n * (if f.kind == 'c': 1 else: 8)
for fi, f in layout:
  if f.kind != 'f' or f.dims.len != 2 or f.dims[1] != 4000: continue
  let rows = f.dims[0]
  for r in 0 ..< rows:
    var cnt = 0
    var elo = Inf
    var ehi = -Inf
    var maxRel = 0.0
    var ex = (0.0, 0.0, 0.0)
    for i in 0 ..< 4000:
      let o = offs[fi] + (i * rows + r) * 8
      let x = f64(a, o)
      let y = f64(b, o)
      if x == y or (x.isNaN and y.isNaN): continue
      inc cnt
      let e = f64(a, egOff + i * 8)
      elo = min(elo, e); ehi = max(ehi, e)
      let rel = (if max(abs(x), abs(y)) > 0: abs(x - y) / max(abs(x), abs(y)) else: 0.0)
      if rel > maxRel or rel.isNaN: maxRel = rel; ex = (e, x, y)
    if cnt > 0:
      echo &"{f.name}({r+1},*): {cnt:4} pts  E {elo:.4g}..{ehi:.4g} eV  max rel {maxRel:.3g} (E={ex[0]:.4g}: {ex[1]:.6g} -> {ex[2]:.6g})"
