## Nim code generation: typing, constant folding, GOTO structuring, emission.
##
## Every construct is translated so that the generated Nim performs exactly
## the floating point operations gfortran (-O0 -ffp-contract=off) performs:
## same operand order and grouping, same kinds (REAL*4 arithmetic stays
## float32), REAL*4 literals emitted as their exact widened double value.

import std/[tables, sets, strutils, strformat, algorithm, math, sequtils]
import fast, fsema

proc c_strtod(s: cstring, endp: pointer): cdouble {.importc: "strtod", header: "<stdlib.h>".}
proc c_strtof(s: cstring, endp: pointer): cfloat {.importc: "strtof", header: "<stdlib.h>".}
proc powidf2(x: cdouble, n: cint): cdouble {.importc: "__powidf2".}
proc powisf2(x: cfloat, n: cint): cfloat {.importc: "__powisf2".}

# ---------------------------------------------------------------- constants

type
  Const = object
    typ: TypeKind
    i: int64
    f: float64      ## for tkReal4 holds the exact float32 value

proc isNum(t: TypeKind): bool = t in {tkInt4, tkInt8, tkReal4, tkReal8}
proc isInt(t: TypeKind): bool = t in {tkInt4, tkInt8}
proc isReal(t: TypeKind): bool = t in {tkReal4, tkReal8}
proc maxNum(a, b: TypeKind): TypeKind =
  if not (a.isNum and b.isNum): serr(&"non-numeric operands {a} {b}")
  if ord(a) >= ord(b): a else: b

proc realLiteral(text: string, dbl: bool): float64 =
  let t = text.replace('D', 'E')
  if dbl: float64(c_strtod(t.cstring, nil))
  else: float64(c_strtof(t.cstring, nil))

proc convConst(c: Const, to: TypeKind): Const =
  ## Fortran conversion of a constant, as gfortran folds it
  result.typ = to
  case to
  of tkInt4, tkInt8:
    if c.typ.isInt: result.i = c.i
    else: result.i = int64(trunc(c.f))
    if to == tkInt4: result.i = int64(int32(result.i))
  of tkReal4:
    result.f = (if c.typ.isInt: float64(float32(c.i)) else: float64(float32(c.f)))
  of tkReal8:
    result.f = (if c.typ.isInt: float64(c.i) else: c.f)
  else: serr("bad constant conversion")

proc arith(op: string, a, b: Const, t: TypeKind): Const =
  let a = convConst(a, t)
  let b = convConst(b, t)
  result.typ = t
  case t
  of tkInt4, tkInt8:
    case op
    of "+": result.i = a.i + b.i
    of "-": result.i = a.i - b.i
    of "*": result.i = a.i * b.i
    of "/": result.i = a.i div b.i
    else: serr("bad int const op " & op)
    if t == tkInt4: result.i = int64(int32(result.i))
  of tkReal4:
    let x = float32(a.f)
    let y = float32(b.f)
    var r: float32
    case op
    of "+": r = x + y
    of "-": r = x - y
    of "*": r = x * y
    of "/": r = x / y
    else: serr("bad real const op " & op)
    result.f = float64(r)
  of tkReal8:
    case op
    of "+": result.f = a.f + b.f
    of "-": result.f = a.f - b.f
    of "*": result.f = a.f * b.f
    of "/": result.f = a.f / b.f
    else: serr("bad real const op " & op)
  else: serr("bad const type")

# ---------------------------------------------------------------- generator state

type
  Gen* = ref object
    prog*: Program
    ui*: UnitInfo
    lines*: seq[string]
    tmpCount*: int
    # goto structuring, per unit
    labelList*: Table[int, int]          ## label -> list id
    labelIdx*: Table[int, int]           ## label -> index in its list
    gotoCount*: Table[int, int]          ## label -> number of jumps to it
    listGotos*: Table[int, seq[tuple[target, anchor: int, st: pointer]]]
    gotoMode*: Table[(pointer, int), string]  ## (jump stmt, label) -> block name to break
    listIds*: Table[pointer, int]
    constFoldNotes*: seq[string]        ## intrinsic calls on constants (verify!)
    warnings*: seq[string]

proc nimType*(t: TypeKind): string =
  case t
  of tkInt4: "int32"
  of tkInt8: "int"
  of tkReal4: "float32"
  of tkReal8: "float64"
  of tkLogical: "bool"
  of tkChar: "string"
  of tkNone: serr("no type")

proc nimDeclType*(s: Sym): string =
  let et = nimType(s.typ)
  case s.dims.len
  of 0: et
  of 1:
    if s.dims[0] < 0: serr(&"assumed-size array {s.name} not supported")
    &"array[1..{s.dims[0]}, {et}]"
  of 2: &"Arr2[{et}, {s.dims[0]}, {s.dims[1]}]"
  of 3: &"Arr3[{et}, {s.dims[0]}, {s.dims[1]}, {s.dims[2]}]"
  else: serr(&"rank {s.dims.len} array {s.name}")

proc fmtF64*(v: float64): string =
  result = $v
  if 'e' notin result and '.' notin result and "inf" notin result and "nan" notin result:
    result.add ".0"

proc literal(c: Const): string =
  case c.typ
  of tkInt4, tkInt8: (if c.i < 0: "(" & $c.i & ")" else: $c.i)
  of tkReal8: (if c.f < 0 or (c.f == 0 and signbit(c.f)): "(" & fmtF64(c.f) & ")" else: fmtF64(c.f))
  of tkReal4:
    let s = fmtF64(c.f) & "'f32"
    if c.f < 0 or (c.f == 0 and signbit(c.f)): "(" & s & ")" else: s
  else: serr("bad literal")

# ---------------------------------------------------------------- intrinsics

type IntrInfo = object
  nim: string          ## Nim function name ("" = special)
  argT: TypeKind       ## required arg type (tkNone = generic)
  resT: TypeKind       ## tkNone = same as (promoted) argument type

const intrinsics = {
  "DSQRT": IntrInfo(nim: "sqrt", argT: tkReal8), "SQRT": IntrInfo(nim: "sqrt"),
  "DEXP": IntrInfo(nim: "exp", argT: tkReal8), "EXP": IntrInfo(nim: "exp"),
  "DLOG": IntrInfo(nim: "ln", argT: tkReal8), "LOG": IntrInfo(nim: "ln"),
  "ALOG": IntrInfo(nim: "ln", argT: tkReal4),
  "DLOG10": IntrInfo(nim: "log10", argT: tkReal8), "LOG10": IntrInfo(nim: "log10"),
  "DSIN": IntrInfo(nim: "sin", argT: tkReal8), "SIN": IntrInfo(nim: "sin"),
  "DCOS": IntrInfo(nim: "cos", argT: tkReal8), "COS": IntrInfo(nim: "cos"),
  "DTAN": IntrInfo(nim: "tan", argT: tkReal8), "TAN": IntrInfo(nim: "tan"),
  "DASIN": IntrInfo(nim: "arcsin", argT: tkReal8), "ASIN": IntrInfo(nim: "arcsin"),
  "DACOS": IntrInfo(nim: "arccos", argT: tkReal8), "ACOS": IntrInfo(nim: "arccos"),
  "DATAN": IntrInfo(nim: "arctan", argT: tkReal8), "ATAN": IntrInfo(nim: "arctan"),
  "DATAN2": IntrInfo(nim: "arctan2", argT: tkReal8), "ATAN2": IntrInfo(nim: "arctan2"),
  "DABS": IntrInfo(nim: "abs", argT: tkReal8), "ABS": IntrInfo(nim: "abs"),
  "IABS": IntrInfo(nim: "abs", argT: tkInt4),
  "DMAX1": IntrInfo(nim: "max", argT: tkReal8), "DMIN1": IntrInfo(nim: "min", argT: tkReal8),
  "AMAX1": IntrInfo(nim: "max", argT: tkReal4), "AMIN1": IntrInfo(nim: "min", argT: tkReal4),
  "MAX0": IntrInfo(nim: "max", argT: tkInt4), "MIN0": IntrInfo(nim: "min", argT: tkInt4),
  "MAX": IntrInfo(nim: "max"), "MIN": IntrInfo(nim: "min"),
  "MOD": IntrInfo(nim: "mod"),
  "DINT": IntrInfo(nim: "trunc", argT: tkReal8), "AINT": IntrInfo(nim: "trunc"),
  "DFLOAT": IntrInfo(nim: "", resT: tkReal8), "DBLE": IntrInfo(nim: "", resT: tkReal8),
  "FLOAT": IntrInfo(nim: "", resT: tkReal4), "REAL": IntrInfo(nim: "", resT: tkReal4),
  "SNGL": IntrInfo(nim: "", argT: tkReal8, resT: tkReal4),
  "INT": IntrInfo(nim: "", resT: tkInt4), "IFIX": IntrInfo(nim: "", argT: tkReal4, resT: tkInt4),
  "IDINT": IntrInfo(nim: "", argT: tkReal8, resT: tkInt4),
  "NINT": IntrInfo(nim: "nint", resT: tkInt4), "IDNINT": IntrInfo(nim: "nint", argT: tkReal8, resT: tkInt4),
  "DSIGN": IntrInfo(nim: "fsign", argT: tkReal8), "SIGN": IntrInfo(nim: "fsign"),
}.toTable

# ---------------------------------------------------------------- typing

proc isUserFunc(g: Gen, name: string): bool =
  name in g.prog.units and g.prog.units[name].u.kind == ukFunction

proc resolve*(g: Gen, e: Expr) =
  let ui = g.ui
  case e.kind
  of eInt: e.typ = tkInt4
  of eReal: e.typ = (if e.dbl: tkReal8 else: tkReal4)
  of eStr: e.typ = tkChar
  of eLog: e.typ = tkLogical
  of eName:
    let s = ui.getOrImplicit(e.name)
    s.used = true
    e.refKind = (if s.dims.len > 0: rkArray else: rkScalar)
    e.typ = s.typ
  of eRef:
    for a in e.args: g.resolve(a)
    let s = ui.lookup(e.name)
    if s != nil and s.dims.len > 0:
      s.used = true
      e.refKind = rkArray
      e.typ = s.typ
      if e.args.len != s.dims.len: serr(&"{ui.u.name}: rank mismatch for {e.name}")
    elif e.name == "DRAND48":
      e.refKind = rkFunc
      e.typ = tkReal8
    elif g.isUserFunc(e.name):
      e.refKind = rkFunc
      e.typ = (if s != nil: s.typ else: implicitType(ui.u, e.name))
      if s != nil: s.storage = stParam   # declared type only; not a variable
      ui.callees.incl e.name
    elif e.name in intrinsics:
      e.refKind = rkIntrinsic
      let info = intrinsics[e.name]
      var t = tkNone
      for a in e.args:
        t = (if t == tkNone: a.typ else: maxNum(t, a.typ))
      if info.argT != tkNone: t = info.argT
      e.typ = (if info.resT != tkNone: info.resT else: t)
    else:
      serr(&"{ui.u.name}: unknown function or undeclared array {e.name}")
  of eUn:
    g.resolve(e.a)
    e.typ = (if e.op == "NOT": tkLogical else: e.a.typ)
  of eBin:
    g.resolve(e.a)
    g.resolve(e.b)
    case e.op
    of "+", "-", "*", "/": e.typ = maxNum(e.a.typ, e.b.typ)
    of "**":
      e.typ = (if e.b.typ.isInt: e.a.typ else: maxNum(e.a.typ, e.b.typ))
    of "EQ", "NE", "LT", "LE", "GT", "GE", "AND", "OR": e.typ = tkLogical
    else: serr("bad op " & e.op)

proc fold(g: Gen, e: Expr): (bool, Const) =
  ## constant value of e (literals, PARAMETERs, + - * / **int)
  case e.kind
  of eInt: (true, Const(typ: tkInt4, i: e.ival))
  of eReal: (true, Const(typ: e.typ, f: realLiteral(e.rtext, e.dbl)))
  of eName:
    let s = g.ui.lookup(e.name)
    if s != nil and s.storage == stParam and s.pval != nil:
      let (ok, c) = g.fold(s.pval)
      if ok: return (true, convConst(c, s.typ))
    (false, Const())
  of eUn:
    if e.op == "NOT": return (false, Const())
    let (ok, c) = g.fold(e.a)
    if not ok: return (false, Const())
    var r = c
    if e.op == "-":
      if r.typ.isInt: r.i = -r.i else: r.f = -r.f
    (true, r)
  of eBin:
    if e.op notin ["+", "-", "*", "/", "**"]: return (false, Const())
    let (oka, a) = g.fold(e.a)
    let (okb, b) = g.fold(e.b)
    if not (oka and okb): return (false, Const())
    if e.op == "**":
      if b.typ.isInt and a.typ.isReal:
        var r = Const(typ: a.typ)
        if b.i == 2: r = arith("*", a, a, a.typ)
        elif a.typ == tkReal8: r.f = powidf2(a.f, cint(b.i))
        else: r.f = float64(powisf2(float32(a.f), cint(b.i)))
        return (true, r)
      g.constFoldNotes.add &"{g.ui.u.name}: constant power {e}"
      return (false, Const())
    (true, arith(e.op, a, b, e.typ))
  else: (false, Const())

# ---------------------------------------------------------------- expression emission

const precAtom = 10
proc opPrec(op: string): int =
  case op
  of "*", "/": 9
  of "+", "-": 8
  of "EQ", "NE", "LT", "LE", "GT", "GE": 5
  of "AND": 4
  of "OR": 3
  else: 10

proc conv(code: string, prec: int, frm, to: TypeKind): string =
  if frm == to: return code
  if to == tkChar or frm == tkChar or to == tkLogical or frm == tkLogical:
    serr(&"cannot convert {frm} to {to}: {code}")
  &"{nimType(to)}({code})"

proc emitE*(g: Gen, e: Expr): tuple[code: string, prec: int]

proc emitAs*(g: Gen, e: Expr, to: TypeKind): string =
  ## e converted to type `to` (Fortran conversion semantics)
  let (ok, c) = g.fold(e)
  if ok and to.isNum: return literal(convConst(c, to))
  let (code, prec) = g.emitE(e)
  if e.typ == tkReal4 and to == tkReal8 and e.kind == eReal:
    return literal(convConst(Const(typ: tkReal4, f: realLiteral(e.rtext, false)), tkReal8))
  conv(code, prec, e.typ, to)

proc paren(x: tuple[code: string, prec: int], minPrec: int): string =
  if x.prec < minPrec: "(" & x.code & ")" else: x.code

proc emitOperand(g: Gen, e: Expr, t: TypeKind, minPrec: int): string =
  let (ok, c) = g.fold(e)
  if ok: return literal(convConst(c, t))
  if e.typ == t:
    return paren(g.emitE(e), minPrec)
  g.emitAs(e, t)    # conversion call -> atom

proc isSimpleAtom(e: Expr): bool =
  e.kind in {eName, eInt, eReal} or (e.kind == eRef and e.refKind == rkArray and
    e.args.allIt(it.kind in {eName, eInt}))

proc arrayRef(g: Gen, e: Expr): string
proc symRef(g: Gen, s: Sym): string

proc emitCallArgs(g: Gen, callee: UnitInfo, args: seq[Expr], pre: var seq[string]): string

proc emitE*(g: Gen, e: Expr): tuple[code: string, prec: int] =
  let (ok, c) = g.fold(e)
  if ok and e.typ.isNum: return (literal(c), precAtom)
  case e.kind
  of eInt, eReal: (literal(c), precAtom)
  of eStr: (escape(e.sval), precAtom)
  of eLog: ((if e.bval: "true" else: "false"), precAtom)
  of eName:
    let s = g.ui.lookup(e.name)
    (g.symRef(s), precAtom)
  of eRef:
    case e.refKind
    of rkArray: (g.arrayRef(e), precAtom)
    of rkFunc:
      if e.name == "DRAND48": return ("drand48(rng)", precAtom)
      let callee = g.prog.units[e.name]
      var pre: seq[string]
      let args = g.emitCallArgs(callee, e.args, pre)
      if pre.len > 0: serr(&"{g.ui.u.name}: function {e.name} needs temporaries")
      (&"{callee.procName}({args})", precAtom)
    of rkIntrinsic:
      let info = intrinsics[e.name]
      var argT = info.argT
      if argT == tkNone:
        for a in e.args:
          argT = (if argT == tkNone: a.typ else: maxNum(argT, a.typ))
      let args = e.args.mapIt(g.emitAs(it, argT))
      if e.args.allIt(g.fold(it)[0]):
        g.constFoldNotes.add &"{g.ui.u.name}: {e}"
      let r = (case e.name
        of "DFLOAT", "DBLE", "FLOAT", "REAL", "SNGL": conv(args[0], precAtom, argT, e.typ)
        of "INT", "IFIX", "IDINT": conv(args[0], precAtom, argT, tkInt4)
        of "MOD": "(" & args[0] & " mod " & args[1] & ")"
        of "MAX", "MIN", "DMAX1", "DMIN1", "AMAX1", "AMIN1", "MAX0", "MIN0":
          (var s = args[0]; (for a in args[1..^1]: s = &"{info.nim}({s}, {a})"); s)
        else: &"{info.nim}({args.join(\", \")})")
      (r, precAtom)
    else: serr("unresolved reference " & e.name)
  of eUn:
    if e.op == "NOT":
      return ("not (" & g.emitE(e.a).code & ")", precAtom)
    let a = g.emitE(e.a)
    if e.op == "+": return a
    if a.prec >= precAtom and not a.code.startsWith("-") and not a.code.startsWith("("):
      ("-" & a.code, precAtom)
    else:
      ("-(" & a.code & ")", precAtom)
  of eBin:
    case e.op
    of "+", "-", "*", "/":
      let p = opPrec(e.op)
      let l = g.emitOperand(e.a, e.typ, p)
      let r = g.emitOperand(e.b, e.typ, p+1)
      let op = (if e.op == "/" and e.typ.isInt: "div" else: e.op)
      (&"{l} {op} {r}", p)
    of "**":
      let (okb, cb) = g.fold(e.b)
      if e.b.typ.isInt:
        if not okb: serr(&"{g.ui.u.name}: non-constant integer exponent {e}")
        if e.a.typ.isReal:
          if cb.i == 2:
            if isSimpleAtom(e.a):
              let x = g.emitOperand(e.a, e.typ, 10)
              return (&"{x}*{x}", 9)
            return (&"pow2({g.emitAs(e.a, e.typ)})", precAtom)
          return (&"powi({g.emitAs(e.a, e.typ)}, {cb.i})", precAtom)
        return (&"ipow({g.emitAs(e.a, e.typ)}, {cb.i})", precAtom)
      (&"pow({g.emitAs(e.a, e.typ)}, {g.emitAs(e.b, e.typ)})", precAtom)
    of "EQ", "NE", "LT", "LE", "GT", "GE":
      let t = (if e.a.typ == tkChar: tkChar else: maxNum(e.a.typ, e.b.typ))
      let p = 5
      let l = (if t == tkChar: g.emitE(e.a).code else: g.emitOperand(e.a, t, p+1))
      let r = (if t == tkChar: g.emitE(e.b).code else: g.emitOperand(e.b, t, p+1))
      let op = {"EQ": "==", "NE": "!=", "LT": "<", "LE": "<=", "GT": ">", "GE": ">="}.toTable[e.op]
      (&"{l} {op} {r}", p)
    of "AND", "OR":
      let p = opPrec(e.op)
      let l = paren(g.emitE(e.a), p)
      let r = paren(g.emitE(e.b), p+1)
      (&"{l} {e.op.toLowerAscii} {r}", p)
    else: serr("bad op " & e.op)

proc commonVar(g: Gen, s: Sym): string =
  let blk = g.prog.commons[s.blk]
  blk.nim & "." & blk.members[s.idx].nim

proc symRef(g: Gen, s: Sym): string =
  case s.storage
  of stCommon: g.commonVar(s)
  of stResult: "result"
  of stParam: s.nim
  else: s.nim

proc arrayRef(g: Gen, e: Expr): string =
  let s = g.ui.lookup(e.name)
  let idx = e.args.mapIt(g.emitAs(it, tkInt8))
  g.symRef(s) & "[" & idx.join(", ") & "]"

# ---------------------------------------------------------------- calls

proc isLvalue(e: Expr): bool =
  e.kind == eName or (e.kind == eRef and e.refKind == rkArray)

proc emitCallArgs(g: Gen, callee: UnitInfo, args: seq[Expr], pre: var seq[string]): string =
  let cu = callee.u
  if args.len != cu.params.len:
    serr(&"{g.ui.u.name}: {cu.name} called with {args.len} args, has {cu.params.len}")
  var res: seq[string]
  for i, a in args:
    let d = callee.syms[cu.params[i]]
    let isOut = callee.outParams[i]
    if d.dims.len > 0:
      # array dummy: actual must be a whole array of the same shape
      if a.kind != eName or a.refKind != rkArray:
        serr(&"{g.ui.u.name}: array element/expression passed to array dummy {cu.name}.{d.name}: {a}")
      let s = g.ui.lookup(a.name)
      if s.dims != d.dims or s.typ != d.typ:
        serr(&"{g.ui.u.name}: shape/type mismatch passing {a.name}{s.dims} to {cu.name}.{d.name}{d.dims}")
      res.add g.symRef(s)
    elif isOut:
      if a.isLvalue and a.typ == d.typ and not (a.kind == eName and a.refKind == rkArray):
        res.add g.emitE(a).code
      else:
        inc g.tmpCount
        let t = &"tmp{g.tmpCount}"
        pre.add &"var {t} = {g.emitAs(a, d.typ)}"
        res.add t
    else:
      if d.typ == tkChar: res.add g.emitE(a).code
      else: res.add g.emitAs(a, d.typ)
  res.join(", ")

# ---------------------------------------------------------------- statements

type
  IvKind = enum ivBlock, ivLoop
  Interval = object
    kind: IvKind
    label: int
    s, e: int
    simple: bool

proc add(g: Gen, ind: int, line: string) =
  g.lines.add repeat("  ", ind) & line

proc emitComments(g: Gen, ind: int, cs: seq[string]) =
  for c in cs:
    let t = c.strip
    if t.len == 0 or t.allCharsInSet({'-', '*', '=', 'C', 'c', '+', '.', '/'}): continue
    g.add ind, "# " & t

proc emitList*(g: Gen, list: seq[Stmt], ind: int)

proc gotoCode(g: Gen, st: Stmt, label: int): string =
  let m = g.gotoMode.getOrDefault((cast[pointer](st), label), "")
  if m == "": serr(&"{g.ui.u.name}: unstructured jump to {label}")
  "break " & m

proc emitAssign(g: Gen, st: Stmt, ind: int) =
  let lhs = st.lhs
  var target: Sym
  if lhs.kind == eName:
    target = g.ui.getOrImplicit(lhs.name)
    if target.dims.len > 0: serr(&"{g.ui.u.name}: whole-array assignment {lhs.name}")
  elif lhs.kind == eRef and lhs.refKind == rkArray:
    target = g.ui.lookup(lhs.name)
  else:
    serr(&"{g.ui.u.name}: statement function or bad assignment at line {st.line}: {lhs}")
  let l = g.emitE(lhs).code
  if target.typ == tkChar:
    g.add ind, &"{l} = fstr({g.emitE(st.rhs).code}, {target.clen})"
  else:
    g.add ind, &"{l} = {g.emitAs(st.rhs, target.typ)}"

proc emitIoItems(g: Gen, items: seq[IoItem], ind: int, sink: string, isRead: bool) =
  for it in items:
    case it.kind
    of ioExpr:
      if isRead: g.add ind, &"rd.get({g.emitE(it.e).code})"
      else: g.add ind, &"{sink}.add toFVal({g.emitE(it.e).code})"
    of ioDo:
      let v = g.ui.getOrImplicit(it.doVar)
      var hdr = &"fdo {g.symRef(v)}, {g.emitAs(it.lo, v.typ)}, {g.emitAs(it.hi, v.typ)}"
      if it.step != nil: hdr.add ", " & g.emitAs(it.step, v.typ)
      g.add ind, hdr & ":"
      g.emitIoItems(it.items, ind+1, sink, isRead)

proc hasDo(items: seq[IoItem]): bool = items.anyIt(it.kind == ioDo)

proc emitStmt(g: Gen, st: Stmt, ind: int) =
  g.emitComments(ind, st.comments)
  case st.kind
  of sAssign: g.emitAssign(st, ind)
  of sIf:
    let c = g.emitE(st.cond).code
    if st.then.kind == sGoto:
      g.add ind, &"if {c}: {g.gotoCode(st, st.then.target)}"
    else:
      let save = g.lines.len
      g.emitStmt(st.then, ind+1)
      let body = g.lines[save..^1]
      g.lines.setLen(save)
      if body.len == 1:
        g.add ind, &"if {c}: " & body[0].strip
      else:
        g.add ind, &"if {c}:"
        g.lines.add body
  of sIfBlock:
    for i, cond in st.conds:
      let kw = (if i == 0: "if" else: "elif")
      g.add ind, &"{kw} {g.emitE(cond).code}:"
      g.emitList(st.bodies[i], ind+1)
    if st.hasElse:
      g.add ind, "else:"
      g.emitList(st.bodies[^1], ind+1)
  of sDo:
    let v = g.ui.getOrImplicit(st.dvar)
    v.written = true
    var hdr = &"fdo {g.symRef(v)}, {g.emitAs(st.dlo, v.typ)}, {g.emitAs(st.dhi, v.typ)}"
    if st.dstep != nil: hdr.add ", " & g.emitAs(st.dstep, v.typ)
    g.add ind, hdr & ":"
    g.emitList(st.body, ind+1)
  of sGoto: g.add ind, g.gotoCode(st, st.target)
  of sCGoto:
    g.add ind, &"case {g.emitAs(st.sel, tkInt8)}"
    for i, t in st.targets:
      g.add ind, &"of {i+1}: {g.gotoCode(st, t)}"
    g.add ind, "else: discard"
  of sCall:
    if st.callee notin g.prog.units: serr(&"{g.ui.u.name}: call of unknown {st.callee}")
    let callee = g.prog.units[st.callee]
    var pre: seq[string]
    let args = g.emitCallArgs(callee, st.cargs, pre)
    if pre.len > 0:
      g.add ind, "block:"
      for p in pre: g.add ind+1, p
      g.add ind+1, &"{callee.procName}({args})"
    else:
      g.add ind, &"{callee.procName}({args})"
  of sReturn: g.add ind, "return"
  of sStop: g.add ind, "fstop()"
  of sContinue: discard
  of sWrite:
    let fmt = (if st.fmtInline.len > 0: escape("(" & st.fmtInline.strip(chars = {'(', ')'}) & ")")
               elif st.fmtLabel == 0: ""
               else: "fmt" & $st.fmtLabel)
    if st.fmtLabel != 0 and st.fmtLabel notin g.ui.formats:
      serr(&"{g.ui.u.name}: missing FORMAT {st.fmtLabel}")
    if hasDo(st.witems):
      g.add ind, "block:"
      g.add ind+1, "var w: seq[FVal]"
      g.emitIoItems(st.witems, ind+1, "w", false)
      if fmt == "": g.add ind+1, "fwriteList(w)"
      else: g.add ind+1, &"fwrite({fmt}, w)"
    else:
      var args: seq[string]
      for it in st.witems: args.add g.emitE(it.e).code
      if fmt == "": g.add ind, &"fwriteList({args.join(\", \")})"
      elif args.len == 0: g.add ind, &"fwrite({fmt})"
      else: g.add ind, &"fwrite({fmt}, {args.join(\", \")})"
  of sRead:
    g.add ind, "block:"
    g.add ind+1, "var rd = listRead()"
    g.emitIoItems(st.ritems, ind+1, "", true)
  of sFormat: discard

# ---------------------------------------------------------------- goto structuring

proc listId(g: Gen, list: seq[Stmt]): int =
  ## a statement list is identified by its first statement (refs are unique)
  if list.len == 0: return -1
  let p = cast[pointer](list[0])
  if p notin g.listIds: g.listIds[p] = g.listIds.len
  g.listIds[p]

proc collectGotos*(g: Gen, list: seq[Stmt], path: var seq[tuple[list, idx: int]]) =
  ## records labels and every jump with its anchor in each enclosing list
  let lid = g.listId(list)
  for i, st in list:
    if st.label != 0:
      if st.label in g.labelList: serr(&"{g.ui.u.name}: duplicate label {st.label}")
      g.labelList[st.label] = lid
      g.labelIdx[st.label] = i
    path.add (lid, i)
    proc jump(g: Gen, target: int, path: seq[tuple[list, idx: int]], st: Stmt) =
      g.gotoCount.mgetOrPut(target, 0).inc
      for (l, k) in path:
        g.listGotos.mgetOrPut(l, @[]).add (target, k, cast[pointer](st))
    case st.kind
    of sGoto: g.jump(st.target, path, st)
    of sCGoto:
      for t in st.targets: g.jump(t, path, st)
    of sIf:
      if st.then.kind == sGoto: g.jump(st.then.target, path, st)
      elif st.then.kind == sCGoto: serr("IF(..) GO TO (..) not supported")
    of sIfBlock:
      for b in st.bodies: g.collectGotos(b, path)
    of sDo: g.collectGotos(st.body, path)
    else: discard
    discard path.pop

proc isIfGoto(st: Stmt, label: int): bool =
  st.kind == sIf and st.then.kind == sGoto and st.then.target == label

proc emitRange(g: Gen, list: seq[Stmt], lo, hi: int, ivs: seq[Interval], ind: int)

proc emitInterval(g: Gen, list: seq[Stmt], iv: Interval, inner: seq[Interval], ind: int) =
  let L = (if iv.kind == ivBlock: "L" else: "R") & $iv.label
  case iv.kind
  of ivBlock:
    if iv.simple:
      # IF(c) GO TO L ; stmts ; L ...   ->   if not c: stmts
      let st = list[iv.s]
      g.emitComments(ind, st.comments)
      if iv.e > iv.s:
        g.add ind, &"if not ({g.emitE(st.cond).code}):"
        g.emitRange(list, iv.s+1, iv.e, inner, ind+1)
    else:
      g.add ind, &"block {L}:"
      g.emitRange(list, iv.s, iv.e, inner, ind+1)
  of ivLoop:
    if iv.simple:
      let last = list[iv.e]
      g.add ind, "while true:"
      g.emitRange(list, iv.s, iv.e-1, inner, ind+1)
      g.emitComments(ind+1, last.comments)
      if last.kind == sIf:
        g.add ind+1, &"if not ({g.emitE(last.cond).code}): break"
    else:
      g.add ind, &"block {L}x:"
      g.add ind+1, "while true:"
      g.add ind+2, &"block {L}:"
      g.emitRange(list, iv.s, iv.e, inner, ind+3)
      g.add ind+3, &"break {L}x"

proc emitRange(g: Gen, list: seq[Stmt], lo, hi: int, ivs: seq[Interval], ind: int) =
  let start = g.lines.len
  var i = lo
  while i <= hi:
    # outermost interval starting here
    var best = -1
    for k, iv in ivs:
      if iv.s == i and (best < 0 or iv.e > ivs[best].e or
          (iv.e == ivs[best].e and iv.kind == ivBlock and ivs[best].kind == ivLoop)):
        best = k
    if best >= 0:
      let iv = ivs[best]
      var inner: seq[Interval]
      for k, x in ivs:
        if k != best and x.s >= iv.s and x.e <= iv.e: inner.add x
      g.emitInterval(list, iv, inner, ind)
      i = iv.e + 1
    else:
      g.emitStmt(list[i], ind)
      inc i
  if g.lines.len == start:
    g.add ind, "discard"

proc emitList*(g: Gen, list: seq[Stmt], ind: int) =
  let lid = g.listId(list)
  var ivs: seq[Interval]
  var fwd, bwd: Table[int, Interval]
  for (target, anchor, stp) in g.listGotos.getOrDefault(lid, @[]):
    if g.labelList.getOrDefault(target, -1) != lid: continue
    let j = g.labelIdx[target]
    g.gotoMode[(stp, target)] = (if anchor < j: "L" else: "R") & $target
    if anchor < j:
      if target in fwd: fwd[target].s = min(fwd[target].s, anchor)
      else: fwd[target] = Interval(kind: ivBlock, label: target, s: anchor, e: j-1)
    else:
      if target in bwd: bwd[target].e = max(bwd[target].e, anchor)
      else: bwd[target] = Interval(kind: ivLoop, label: target, s: j, e: anchor)
  # labels of this list must only be targeted from inside it
  for (l, li) in g.labelList.pairs:
    if li == lid:
      var inside = 0
      for (target, anchor, stp) in g.listGotos.getOrDefault(lid, @[]):
        if target == l: inc inside
      if inside != g.gotoCount.getOrDefault(l, 0):
        serr(&"{g.ui.u.name}: jump into a block/loop to label {l}")
  for t, iv in fwd:
    var iv = iv
    iv.simple = isIfGoto(list[iv.s], t) and g.gotoCount[t] == 1
    ivs.add iv
  for t, iv in bwd:
    var iv = iv
    let last = list[iv.e]
    iv.simple = g.gotoCount[t] == 1 and (isIfGoto(last, t) or
      (last.kind == sGoto and last.target == t))
    ivs.add iv
  # make intervals nest properly; loops first (extending a loop's end is
  # always possible, and moving block starts early can be avoided that way)
  var changed = true
  while changed:
    changed = false
    for a in 0 ..< ivs.len:
      for b in 0 ..< ivs.len:
        let A = ivs[a]
        let B = ivs[b]
        if a != b and A.kind == ivLoop and B.kind == ivLoop and
            A.s < B.s and B.s <= A.e and A.e < B.e:
          ivs[a].e = B.e
          ivs[a].simple = false
          changed = true
  changed = true
  while changed:
    changed = false
    for a in 0 ..< ivs.len:
      for b in 0 ..< ivs.len:
        if a == b: continue
        let A = ivs[a]
        let B = ivs[b]
        if A.s < B.s and B.s <= A.e and A.e < B.e:
          if B.kind == ivBlock:
            ivs[b].s = A.s
            ivs[b].simple = false
          elif A.kind == ivLoop:
            ivs[a].e = B.e
            ivs[a].simple = false
          else:
            var msg = &"{g.ui.u.name}: jump from before into loop range of label {B.label}\n"
            for iv in ivs:
              msg.add &"   {iv.kind} {iv.label}: [{iv.s} (line {list[iv.s].line}), {iv.e} (line {list[iv.e].line})]\n"
            serr(msg)
          changed = true
  # simple forms need exclusive use of their first/last statement
  for k in 0 ..< ivs.len:
    if not ivs[k].simple: continue
    for m in 0 ..< ivs.len:
      if m == k: continue
      if ivs[k].kind == ivBlock and ivs[m].s == ivs[k].s: ivs[k].simple = false
      if ivs[k].kind == ivLoop and ivs[m].s >= ivs[k].s and ivs[m].e >= ivs[k].e and
          ivs[m].s <= ivs[k].e: ivs[k].simple = false
  g.emitRange(list, 0, list.len-1, ivs, ind)

# ---------------------------------------------------------------- dummy argument intents

proc markWritten(g: Gen, ui: UnitInfo, e: Expr) =
  let s = ui.lookup(e.name)
  if s != nil: s.written = true

proc scanWrites*(g: Gen, ui: UnitInfo, list: seq[Stmt]) =
  ## marks symbols written by assignments, DO, READ, and out arguments
  proc items(g: Gen, ui: UnitInfo, its: seq[IoItem], isRead: bool) =
    for it in its:
      if it.kind == ioDo:
        ui.getOrImplicit(it.doVar).written = true
        g.items(ui, it.items, isRead)
      elif isRead: g.markWritten(ui, it.e)
  for st in list:
    case st.kind
    of sAssign: g.markWritten(ui, st.lhs)
    of sDo:
      ui.getOrImplicit(st.dvar).written = true
      g.scanWrites(ui, st.body)
    of sIf: g.scanWrites(ui, @[st.then])
    of sIfBlock:
      for b in st.bodies: g.scanWrites(ui, b)
    of sRead: g.items(ui, st.ritems, true)
    of sWrite: g.items(ui, st.witems, false)
    of sCall:
      if st.callee in g.prog.units:
        let c = g.prog.units[st.callee]
        for i, a in st.cargs:
          if i < c.outParams.len and c.outParams[i] and a.kind in {eName, eRef}:
            g.markWritten(ui, a)
    else: discard
