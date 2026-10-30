## AST of the Fortran 77 subset used by Magboltz.

type
  TypeKind* = enum
    tkNone, tkInt4, tkInt8, tkReal4, tkReal8, tkLogical, tkChar

  ExprKind* = enum
    eInt, eReal, eStr, eLog, eName, eRef, eUn, eBin

  RefKind* = enum
    rkUnknown, rkScalar, rkArray, rkIntrinsic, rkFunc

  Expr* = ref object
    typ*: TypeKind
    case kind*: ExprKind
    of eInt: ival*: int64
    of eReal:
      rtext*: string        ## literal text as written (uppercase, no blanks)
      dbl*: bool            ## has D exponent -> REAL*8
    of eStr: sval*: string
    of eLog: bval*: bool
    of eName, eRef:
      name*: string
      args*: seq[Expr]
      refKind*: RefKind
    of eUn, eBin:
      op*: string           ## "+", "-", "*", "/", "**", ".EQ.", ..., ".NOT."
      a*, b*: Expr

  IoItemKind* = enum ioExpr, ioDo
  IoItem* = ref object
    case kind*: IoItemKind
    of ioExpr: e*: Expr
    of ioDo:
      doVar*: string
      lo*, hi*, step*: Expr
      items*: seq[IoItem]

  StmtKind* = enum
    sAssign, sIf, sIfBlock, sDo, sGoto, sCGoto, sCall, sReturn, sStop,
    sContinue, sWrite, sRead, sFormat

  Stmt* = ref object
    label*: int
    line*: int
    comments*: seq[string]   ## comment lines preceding the statement
    case kind*: StmtKind
    of sAssign:
      lhs*, rhs*: Expr
    of sIf:
      cond*: Expr
      then*: Stmt
    of sIfBlock:
      conds*: seq[Expr]
      bodies*: seq[seq[Stmt]] ## one per cond, plus else body if hasElse
      hasElse*: bool
    of sDo:
      dvar*: string
      dlo*, dhi*, dstep*: Expr  ## dstep nil = 1
      dlabel*: int
      body*: seq[Stmt]
    of sGoto: target*: int
    of sCGoto:
      targets*: seq[int]
      sel*: Expr
    of sCall:
      callee*: string
      cargs*: seq[Expr]
    of sReturn, sStop, sContinue: discard
    of sWrite:
      fmtLabel*: int          ## 0 = list directed or inline
      fmtInline*: string      ## inline format string, if given
      witems*: seq[IoItem]
    of sRead:
      ritems*: seq[IoItem]
    of sFormat:
      ftext*: string          ## "(...)" as written

  DeclKind* = enum dImplicit, dDimension, dCommon, dType, dParameter, dSave, dData, dEntry

  Entity* = object
    name*: string
    dims*: seq[Expr]        ## empty = scalar; "*" encoded as nil
    clen*: int              ## CHARACTER*n override per entity (0 = none)

  DataSpec* = object
    names*: seq[Expr]       ## variables / array names / elements
    values*: seq[tuple[rep: int, v: Expr]]

  UnitKind* = enum ukProgram, ukSubroutine, ukFunction

  Unit* = ref object
    kind*: UnitKind
    name*: string
    params*: seq[string]
    resType*: TypeKind      ## for functions with explicit type
    line*: int
    implicitNone*: bool
    implicits*: seq[tuple[lo, hi: char, t: TypeKind]]
    dims*: seq[Entity]
    commons*: seq[tuple[blk: string, ents: seq[Entity]]]
    typeDecls*: seq[tuple[t: TypeKind, clen: int, ents: seq[Entity]]]
    parameters*: seq[tuple[name: string, e: Expr]]
    saves*: seq[string]
    datas*: seq[DataSpec]
    entries*: seq[tuple[name: string, params: seq[string]]]
    body*: seq[Stmt]
    headComments*: seq[string]

proc `$`*(e: Expr): string =
  ## Fortran-ish rendering, for messages.
  if e == nil: return "<nil>"
  case e.kind
  of eInt: $e.ival
  of eReal: e.rtext
  of eStr: "'" & e.sval & "'"
  of eLog: (if e.bval: ".TRUE." else: ".FALSE.")
  of eName: e.name
  of eRef:
    var s = e.name & "("
    for i, a in e.args:
      if i > 0: s.add ","
      s.add $a
    s & ")"
  of eUn: e.op & "(" & $e.a & ")"
  of eBin: "(" & $e.a & e.op & $e.b & ")"
