## Walkers on a GPU field by field ("structure of arrays"): every number of
## a walker in a slot of its own, slot j of walker k at j · stride + k,
## floats (every unit of Unchained too) in one array, `int32` and `uint32`
## in another. A warp's 32 walkers then read and write each slot in one
## piece; one walker after the other, they would read 32 pieces 200 bytes
## or more apart. A kernel copies its walker into a variable of its own
## first (`slots`), walks, and copies it back.
##
## In single precision (single.nim) the floats are `float32`, as the
## demoted kernel has them: the slots are the same.

import std/macros
import rm48

proc slotCode(path, t: NimNode, fs, us, k, stride: NimNode, store: bool,
              nf, nu: var int, code: NimNode) =
  ## code for the slots of `path` (of type t) from slot nf of the floats
  ## and nu of the words on; `code` nil: counts only
  let impl = getTypeImpl(t)
  template at(j: int): NimNode =
    infix(infix(newLit(int32(j)), "*", stride), "+", k)
  template float() =
    if code != nil:
      if store: code.add newAssignment(nnkBracketExpr.newTree(fs, at(nf)), newCall(bindSym"float64", path))
      else: code.add newAssignment(path, newCall(t, nnkBracketExpr.newTree(fs, at(nf))))
    inc nf
  template word() =
    if code != nil:
      if store: code.add newAssignment(nnkBracketExpr.newTree(us, at(nu)), newCall(bindSym"uint32", path))
      else: code.add newAssignment(path, newCall(t, nnkBracketExpr.newTree(us, at(nu))))
    inc nu
  proc floatBase(t: NimNode): bool =
    ## a float, or a distinct float (Unchained's units, distinct of distinct)
    let i = getTypeImpl(t)
    case i.kind
    of nnkDistinctTy: floatBase(i[0])
    of nnkSym: i.typeKind in {ntyFloat64, ntyFloat, ntyFloat32}
    else: false
  case impl.kind
  of nnkObjectTy:
    for defs in impl[2]:
      expectKind defs, nnkIdentDefs
      for f in defs[0 ..< ^2]:
        let name = if f.kind == nnkPostfix: f[1] else: f
        slotCode(newDotExpr(path, name), defs[^2], fs, us, k, stride, store, nf, nu, code)
  of nnkBracketExpr:
    if not impl[0].eqIdent("array"): error("slots: no slot for " & repr(t), path)
    let r = impl[1]
    var lo, hi: BiggestInt
    if r.kind == nnkInfix:
      lo = r[1].intVal
      hi = r[2].intVal
    else:
      let ri = getTypeImpl(r)
      lo = ri[0].intVal
      hi = ri[1].intVal
    for i in lo .. hi:
      slotCode(nnkBracketExpr.newTree(path, newLit(int(i))), impl[2], fs, us, k, stride,
               store, nf, nu, code)
  elif floatBase(t): float()
  elif impl.kind == nnkSym and impl.typeKind in {ntyInt32, ntyUInt32}: word()
  else:
    error("slots: no slot for " & repr(t) & " (" & $impl.kind & ")", path)

macro load*(x: typed, fs, us, k, stride: typed): untyped =
  ## x from its slots: x's slot j from fs[j · stride + k] (floats) or
  ## us[j · stride + k] (`int32` and `uint32`)
  result = newStmtList()
  var nf, nu = 0
  slotCode(x, x.getTypeInst, fs, us, k, stride, false, nf, nu, result)

macro store*(x: typed, fs, us, k, stride: typed): untyped =
  ## x to its slots (`load` the other way)
  result = newStmtList()
  var nf, nu = 0
  slotCode(x, x.getTypeInst, fs, us, k, stride, true, nf, nu, result)

macro slotCounts*(T: typedesc): untyped =
  ## (floats, words): the slots of a T
  var nf, nu = 0
  let dummy = ident"x"
  slotCode(dummy, T.getTypeInst[1], dummy, dummy, dummy, dummy, false, nf, nu, nil)
  nnkTupleConstr.newTree(newLit(nf), newLit(nu))

proc toSlots*[W](ws: openArray[W], stride: int): (seq[float64], seq[uint32]) =
  ## the walkers' slots, with room for `stride` walkers
  const c = slotCounts(W)
  var fs = newSeq[float64](max(c[0] * stride, 1))
  var us = newSeq[uint32](max(c[1] * stride, 1))
  let pf = cast[ptr UncheckedArray[float64]](addr fs[0])
  let pu = cast[ptr UncheckedArray[uint32]](addr us[0])
  for k in 0 ..< ws.len:
    let kk = int32(k)
    let st = int32(stride)
    store(ws[k], pf, pu, kk, st)
  (fs, us)

proc fromSlots*[W](ws: var openArray[W], fs: openArray[float64], us: openArray[uint32],
                   stride: int) =
  ## the walkers from their slots
  let pf = cast[ptr UncheckedArray[float64]](unsafeAddr fs[0])
  let pu = cast[ptr UncheckedArray[uint32]](unsafeAddr us[0])
  for k in 0 ..< ws.len:
    let kk = int32(k)
    let st = int32(stride)
    load(ws[k], pf, pu, kk, st)

proc slotBytes*[W](ws: openArray[W], stride: int, single: bool): (seq[byte], seq[byte]) =
  ## the walkers' slots, as a GPU has them: floats in single
  ## precision if `single` (and Ranmar's numbers then `toSingle`)
  var (fs, us) = when typeof(ws[0].rng) is Ranmar:
      if single:
        var c = @ws
        for w in c.mitems: w.rng.toSingle()
        toSlots(c, stride)
      else: toSlots(ws, stride)
    else: toSlots(ws, stride)
  result[1] = newSeq[byte](us.len * 4)
  copyMem(addr result[1][0], addr us[0], result[1].len)
  if single:
    var f32 = newSeq[float32](fs.len)
    for i, x in fs: f32[i] = float32(x)
    result[0] = newSeq[byte](f32.len * 4)
    copyMem(addr result[0][0], addr f32[0], result[0].len)
  else:
    result[0] = newSeq[byte](fs.len * 8)
    copyMem(addr result[0][0], addr fs[0], result[0].len)

proc fromSlotBytes*[W](ws: var openArray[W], fb, ub: openArray[byte],
                       stride: int, single: bool) =
  ## the walkers from their slots as the GPU has them (`slotBytes`)
  var us = newSeq[uint32](ub.len div 4)
  copyMem(addr us[0], unsafeAddr ub[0], ub.len)
  var fs: seq[float64]
  if single:
    fs = newSeq[float64](fb.len div 4)
    let p = cast[ptr UncheckedArray[float32]](unsafeAddr fb[0])
    for i in 0 ..< fs.len: fs[i] = float64(p[i])
  else:
    fs = newSeq[float64](fb.len div 8)
    copyMem(addr fs[0], unsafeAddr fb[0], fb.len)
  fromSlots(ws, fs, us, stride)
