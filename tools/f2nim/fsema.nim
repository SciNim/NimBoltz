## Symbol tables, typing and whole-program facts for the translator.

import std/[tables, sets, strutils, strformat, algorithm]
import fast

type
  StorageKind* = enum stLocal, stDummy, stCommon, stParam, stResult

  Sym* = ref object
    name*: string            ## Fortran name (uppercase)
    nim*: string             ## emitted identifier
    typ*: TypeKind
    clen*: int
    dims*: seq[int]          ## constant extents (-1 = assumed size)
    storage*: StorageKind
    blk*: string             ## COMMON block name
    idx*: int                ## position in COMMON block / dummy list
    pval*: Expr              ## PARAMETER value
    data*: seq[Expr]         ## DATA values, flattened, in storage order
    saved*: bool
    written*: bool           ## assigned / passed to an out argument
    used*: bool

  UnitInfo* = ref object
    u*: Unit
    syms*: OrderedTable[string, Sym]
    modName*: string         ## Nim module name
    procName*: string
    callees*: OrderedSet[string]
    formats*: Table[int, string]
    outParams*: seq[bool]    ## per dummy: written in this unit or callees

  CommonMember* = object
    name*: string            ## canonical Fortran name
    nim*: string
    typ*: TypeKind
    clen*: int
    dims*: seq[int]

  CommonBlock* = object
    name*: string
    nim*: string             ## global variable name
    typeName*: string
    members*: seq[CommonMember]

  Program* = ref object
    units*: OrderedTable[string, UnitInfo]
    commons*: OrderedTable[string, CommonBlock]
    skip*: HashSet[string]   ## units provided by hand-written Nim code

  SemaError* = object of CatchableError

proc serr*(msg: string) {.noreturn.} =
  raise newException(SemaError, msg)

const nimKeywords* = ["addr", "and", "as", "asm", "bind", "block", "break",
  "case", "cast", "concept", "const", "continue", "converter", "defer",
  "discard", "distinct", "div", "do", "elif", "else", "end", "enum", "except",
  "export", "finally", "for", "from", "func", "if", "import", "in", "include",
  "interface", "is", "isnot", "iterator", "let", "macro", "method", "mixin",
  "mod", "nil", "not", "notin", "object", "of", "or", "out", "proc", "ptr",
  "raise", "ref", "return", "shl", "shr", "static", "template", "try",
  "tuple", "type", "using", "var", "when", "while", "xor", "yield"]

## identifiers the generated code uses unqualified; Fortran names that
## normalise to one of these get renamed
const reservedNames* = ["result", "sqrt", "exp", "ln", "log10", "sin", "cos",
  "tan", "arcsin", "arccos", "arctan", "arctan2", "abs", "max", "min", "trunc",
  "pow", "float64", "float32", "int", "int32", "int64", "fdo", "fstr", "fwrite",
  "fstop", "pow2", "powi", "ipow", "fint", "nint", "rng", "drand48", "echo",
  "len", "high", "low", "ord", "chr", "true", "false", "inc", "dec", "system",
  "math", "frt", "commons", "copySign", "string", "char", "bool", "seq",
  "array", "Arr2", "Arr3", "items", "add", "listRead", "get", "quit", "e",
  "pi", "round", "floor", "ceil", "sgn", "fmod", "isNil", "new", "default"]

proc normId*(s: string): string =
  ## Nim style-insensitive key: first char as is, rest lowercase, no '_'
  if s.len == 0: return s
  result.add s[0]
  for c in s[1..^1]:
    if c != '_': result.add c.toLowerAscii

proc nimName*(fortran: string): string =
  result = fortran.toLowerAscii
  if result in nimKeywords: result = "`" & result & "`"

proc implicitType*(u: Unit, name: string): TypeKind =
  let c = name[0]
  for (lo, hi, t) in u.implicits:
    if c >= lo and c <= hi: return t
  if u.implicitNone: return tkNone
  if c in {'I'..'N'}: tkInt4 else: tkReal4

proc constInt*(e: Expr, params: Table[string, Expr]): int64 =
  case e.kind
  of eInt: e.ival
  of eName:
    if e.name in params: constInt(params[e.name], params)
    else: serr("not a constant: " & e.name)
  of eBin:
    let a = constInt(e.a, params)
    let b = constInt(e.b, params)
    case e.op
    of "+": a+b
    of "-": a-b
    of "*": a*b
    of "/": a div b
    else: serr("bad constant op " & e.op)
  of eUn:
    if e.op == "-": -constInt(e.a, params) else: constInt(e.a, params)
  else: serr("not an integer constant: " & $e)

proc elemSize*(t: TypeKind, clen: int): int =
  case t
  of tkInt4, tkReal4, tkLogical: 4
  of tkInt8, tkReal8: 8
  of tkChar: clen
  of tkNone: 0

proc totalSize*(s: Sym): int =
  result = elemSize(s.typ, s.clen)
  for d in s.dims: result *= max(d, 1)

proc lookup*(ui: UnitInfo, name: string): Sym =
  ui.syms.getOrDefault(name, nil)

proc getOrImplicit*(ui: UnitInfo, name: string): Sym =
  ## symbol for a (scalar) name, creating an implicitly typed local
  result = ui.lookup(name)
  if result == nil:
    let t = implicitType(ui.u, name)
    if t == tkNone: serr(&"{ui.u.name}: undeclared {name} under IMPLICIT NONE")
    result = Sym(name: name, typ: t, storage: stLocal)
    ui.syms[name] = result

proc buildSymbols*(ui: UnitInfo) =
  let u = ui.u
  var params: Table[string, Expr]
  for (n, e) in u.parameters: params[n] = e
  proc declare(name: string): Sym =
    result = ui.lookup(name)
    if result == nil:
      result = Sym(name: name, typ: implicitType(u, name), storage: stLocal)
      ui.syms[name] = result
  proc setDims(s: Sym, dims: seq[Expr]) =
    if dims.len == 0: return
    if s.dims.len > 0: serr(&"{u.name}: dims of {s.name} given twice")
    for d in dims:
      s.dims.add (if d == nil: -1 else: int(constInt(d, params)))
  # dummies first, to keep their order
  for i, p in u.params:
    let s = declare(p)
    s.storage = stDummy
    s.idx = i
  if u.kind == ukFunction:
    let s = declare(u.name)
    s.storage = stResult
    if u.resType != tkNone: s.typ = u.resType
  for (t, clen, ents) in u.typeDecls:
    for e in ents:
      let s = declare(e.name)
      s.typ = t
      s.clen = (if e.clen > 0: e.clen else: clen)
      s.setDims(e.dims)
  for e in u.dims:
    declare(e.name).setDims(e.dims)
  for (blk, ents) in u.commons:
    for i, e in ents:
      let s = declare(e.name)
      if s.storage == stDummy: serr(&"{u.name}: dummy {e.name} in COMMON")
      s.storage = stCommon
      s.blk = blk
      s.setDims(e.dims)
  # COMMON member index: position within all lists of the same block
  var pos: Table[string, int]
  for (blk, ents) in u.commons:
    for e in ents:
      let s = ui.syms[e.name]
      s.idx = pos.getOrDefault(blk, 0)
      pos[blk] = s.idx + 1
  for (n, e) in u.parameters:
    let s = declare(n)
    s.storage = stParam
    s.pval = e
  for n in u.saves:
    declare(n).saved = true
  for d in u.datas:
    var vals: seq[Expr]
    for (rep, v) in d.values:
      for k in 1..rep: vals.add v
    var k = 0
    for n in d.names:
      if n.kind != eName: serr(&"{u.name}: DATA for array elements not supported: {n}")
      let s = declare(n.name)
      var cnt = 1
      for x in s.dims: cnt *= x
      if k + cnt > vals.len: serr(&"{u.name}: DATA too few values for {s.name}")
      s.data = vals[k ..< k+cnt]
      s.saved = true
      k += cnt
    if k != vals.len: serr(&"{u.name}: DATA value count mismatch")
  for s in ui.syms.values:
    if s.typ == tkNone: serr(&"{u.name}: no type for {s.name}")
