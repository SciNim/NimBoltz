## Fixed-form reader, lexer and parser for the Fortran 77 subset used by
## Magboltz. Produces `Unit`s with nested statement trees (block IF, DO).

import std/[strutils, strformat]
import fast

type
  RawStmt* = object
    label*: int
    text*: string        ## columns 7-72, continuation lines joined
    norm*: string        ## uppercase, blanks removed outside quotes
    line*: int
    comments*: seq[string]

  ParseError* = object of CatchableError

proc perr(line: int, msg: string) {.noreturn.} =
  raise newException(ParseError, &"line {line}: {msg}")

proc normalize*(s: string): string =
  ## Remove blanks and uppercase, except inside quoted strings.
  var inStr = false
  var i = 0
  while i < s.len:
    let c = s[i]
    if inStr:
      result.add c
      if c == '\'':
        if i+1 < s.len and s[i+1] == '\'':
          result.add '\''
          inc i
        else:
          inStr = false
    elif c == '\'':
      inStr = true
      result.add c
    elif c notin {' ', '\t'}:
      result.add c.toUpperAscii
    inc i

proc readStatements*(path: string): seq[RawStmt] =
  var comments: seq[string]
  var lineNo = 0
  for raw in lines(path):
    inc lineNo
    let stripped = raw.strip(leading = false)
    if stripped.len == 0 or stripped[0] in {'C', 'c', '*', '!'}:
      comments.add (if stripped.len > 0: stripped[1..^1] else: "")
      continue
    # gfortran pads fixed-form lines with blanks up to column 72; this matters
    # for character constants continued over several lines
    let l = alignLeft(raw, 72)
    let lab = l[0..4].strip
    let isCont = l[5] notin {' ', '0'}
    let text = l[6 ..< 72]
    if isCont:
      if result.len == 0: perr(lineNo, "continuation without statement")
      if lab.len > 0: perr(lineNo, "label on continuation line")
      result[^1].text.add text
    else:
      var label = 0
      if lab.len > 0:
        try: label = parseInt(lab)
        except ValueError: perr(lineNo, "bad label '" & lab & "'")
      result.add RawStmt(label: label, text: text, line: lineNo, comments: comments)
      comments = @[]
  for s in result.mitems:
    s.norm = normalize(s.text)

# ---------------------------------------------------------------- lexer

type
  TokKind = enum tName, tInt, tReal, tStr, tOp, tDot, tEnd
  Tok = object
    kind: TokKind
    s: string

const dotOps = ["EQ", "NE", "LT", "LE", "GT", "GE", "AND", "OR", "NOT",
                "EQV", "NEQV", "TRUE", "FALSE"]

proc isDotOpAt(s: string, i: int): bool =
  ## s[i] == '.', followed by an operator name and '.'
  if i >= s.len or s[i] != '.': return false
  var j = i+1
  while j < s.len and s[j] in Letters: inc j
  j < s.len and s[j] == '.' and s[i+1 ..< j] in dotOps

proc lex(s: string, line: int): seq[Tok] =
  var i = 0
  while i < s.len:
    let c = s[i]
    if c in Letters:
      var j = i
      while j < s.len and s[j] in IdentChars: inc j
      result.add Tok(kind: tName, s: s[i ..< j])
      i = j
    elif c in Digits or (c == '.' and i+1 < s.len and s[i+1] in Digits):
      var j = i
      while j < s.len and s[j] in Digits: inc j
      var isReal = false
      if j < s.len and s[j] == '.' and not isDotOpAt(s, j):
        isReal = true
        inc j
        while j < s.len and s[j] in Digits: inc j
      if j < s.len and s[j] in {'E', 'D'} and
          (j+1 < s.len and (s[j+1] in Digits or
            (s[j+1] in {'+', '-'} and j+2 < s.len and s[j+2] in Digits))):
        isReal = true
        inc j
        if s[j] in {'+', '-'}: inc j
        while j < s.len and s[j] in Digits: inc j
      result.add Tok(kind: (if isReal: tReal else: tInt), s: s[i ..< j])
      i = j
    elif c == '\'':
      var j = i+1
      var v = ""
      while true:
        if j >= s.len: perr(line, "unterminated string")
        if s[j] == '\'':
          if j+1 < s.len and s[j+1] == '\'':
            v.add '\''; j += 2; continue
          break
        v.add s[j]; inc j
      result.add Tok(kind: tStr, s: v)
      i = j+1
    elif c == '.':
      if not isDotOpAt(s, i): perr(line, "bad '.' in " & s)
      var j = i+1
      while s[j] != '.': inc j
      result.add Tok(kind: tDot, s: s[i+1 ..< j])
      i = j+1
    elif c == '*' and i+1 < s.len and s[i+1] == '*':
      result.add Tok(kind: tOp, s: "**"); i += 2
    elif c == '/' and i+1 < s.len and s[i+1] == '/':
      result.add Tok(kind: tOp, s: "//"); i += 2
    elif c in {'+', '-', '*', '/', '(', ')', ',', '=', ':'}:
      result.add Tok(kind: tOp, s: $c); inc i
    else:
      perr(line, "unexpected char '" & c & "' in " & s)
  result.add Tok(kind: tEnd)

# ---------------------------------------------------------------- expressions

type P = object
  toks: seq[Tok]
  pos: int
  line: int

proc peek(p: P): Tok = p.toks[p.pos]
proc isOp(p: P, s: string): bool = p.peek.kind == tOp and p.peek.s == s
proc isDot(p: P, s: string): bool = p.peek.kind == tDot and p.peek.s == s
proc next(p: var P): Tok =
  result = p.toks[p.pos]
  inc p.pos
proc expect(p: var P, s: string) =
  if not p.isOp(s): perr(p.line, "expected '" & s & "' got '" & p.peek.s & "'")
  inc p.pos

proc parseExpr(p: var P): Expr

proc parsePrimary(p: var P): Expr =
  let t = p.next
  case t.kind
  of tInt: Expr(kind: eInt, ival: parseBiggestInt(t.s))
  of tReal: Expr(kind: eReal, rtext: t.s, dbl: 'D' in t.s)
  of tStr: Expr(kind: eStr, sval: t.s)
  of tDot:
    if t.s == "TRUE": Expr(kind: eLog, bval: true)
    elif t.s == "FALSE": Expr(kind: eLog, bval: false)
    else: perr(p.line, "unexpected ." & t.s & ".")
  of tName:
    if p.isOp("("):
      inc p.pos
      var args: seq[Expr]
      if not p.isOp(")"):
        while true:
          args.add p.parseExpr
          if p.isOp(","): inc p.pos
          else: break
      p.expect ")"
      Expr(kind: eRef, name: t.s, args: args)
    else:
      Expr(kind: eName, name: t.s)
  of tOp:
    if t.s == "(":
      let e = p.parseExpr
      p.expect ")"
      e   # parentheses are kept implicitly by tree structure
    else: perr(p.line, "unexpected '" & t.s & "'")
  of tEnd: perr(p.line, "unexpected end of statement")

proc parsePow(p: var P): Expr

proc parseSignedPow(p: var P): Expr =
  ## extension: sign directly after an operator (A*-B, A**-B)
  if p.isOp("-") or p.isOp("+"):
    let op = p.next.s
    return Expr(kind: eUn, op: op, a: p.parsePow)
  p.parsePow

proc parsePow(p: var P): Expr =
  result = p.parsePrimary
  if p.isOp("**"):
    inc p.pos
    result = Expr(kind: eBin, op: "**", a: result, b: p.parseSignedPow)

proc parseMult(p: var P): Expr =
  result = p.parseSignedPow
  while p.isOp("*") or p.isOp("/"):
    let op = p.next.s
    result = Expr(kind: eBin, op: op, a: result, b: p.parseSignedPow)

proc parseAdd(p: var P): Expr =
  if p.isOp("-") or p.isOp("+"):
    let op = p.next.s
    result = Expr(kind: eUn, op: op, a: p.parseMult)
  else:
    result = p.parseMult
  while p.isOp("+") or p.isOp("-"):
    let op = p.next.s
    result = Expr(kind: eBin, op: op, a: result, b: p.parseMult)

proc parseRel(p: var P): Expr =
  result = p.parseAdd
  if p.peek.kind == tDot and p.peek.s in ["EQ", "NE", "LT", "LE", "GT", "GE"]:
    let op = p.next.s
    result = Expr(kind: eBin, op: op, a: result, b: p.parseAdd)

proc parseNot(p: var P): Expr =
  if p.isDot("NOT"):
    inc p.pos
    return Expr(kind: eUn, op: "NOT", a: p.parseNot)
  p.parseRel

proc parseAnd(p: var P): Expr =
  result = p.parseNot
  while p.isDot("AND"):
    inc p.pos
    result = Expr(kind: eBin, op: "AND", a: result, b: p.parseNot)

proc parseOr(p: var P): Expr =
  result = p.parseAnd
  while p.isDot("OR"):
    inc p.pos
    result = Expr(kind: eBin, op: "OR", a: result, b: p.parseAnd)

proc parseExpr(p: var P): Expr = p.parseOr

proc initP(s: string, line: int): P = P(toks: lex(s, line), line: line)

proc parseExprStr*(s: string, line: int): Expr =
  var p = initP(s, line)
  result = p.parseExpr
  if p.peek.kind != tEnd: perr(line, "trailing tokens in expression: " & s)

# ---------------------------------------------------------------- helpers

proc matchParen(s: string, i: int): int =
  ## index of the ')' matching '(' at s[i]
  var d = 0
  var inStr = false
  for j in i ..< s.len:
    let c = s[j]
    if inStr:
      if c == '\'': inStr = false
    elif c == '\'': inStr = true
    elif c == '(': inc d
    elif c == ')':
      dec d
      if d == 0: return j
  -1

proc splitTop*(s: string, sep = ','): seq[string] =
  ## split at depth-0 separators (outside parens and strings)
  var d = 0
  var inStr = false
  var cur = ""
  for c in s:
    if inStr:
      if c == '\'': inStr = false
    elif c == '\'': inStr = true
    elif c == '(': inc d
    elif c == ')': dec d
    elif c == sep and d == 0:
      result.add cur
      cur = ""
      continue
    cur.add c
  result.add cur

proc topIndex(s: string, c: char): int =
  ## first depth-0 occurrence of c, -1 if none
  var d = 0
  var inStr = false
  for i, ch in s:
    if inStr:
      if ch == '\'': inStr = false
    elif ch == '\'': inStr = true
    elif ch == '(': inc d
    elif ch == ')': dec d
    elif ch == c and d == 0: return i
  -1

proc parseEntity(s: string, line: int): Entity =
  let k = s.find('(')
  if k < 0:
    var name = s
    let st = s.find('*')
    if st > 0:
      name = s[0 ..< st]
      result.clen = parseInt(s[st+1..^1])
    return Entity(name: name, clen: result.clen)
  result.name = s[0 ..< k]
  let e = matchParen(s, k)
  for d in splitTop(s[k+1 ..< e]):
    if d == "*": result.dims.add nil
    elif ':' in d: perr(line, "array bounds with ':' not supported: " & s)
    else: result.dims.add parseExprStr(d, line)

proc parseEntities(s: string, line: int): seq[Entity] =
  for part in splitTop(s):
    if part.len > 0: result.add parseEntity(part, line)

proc parseTypeSpec(s: string): tuple[t: TypeKind, clen: int, rest: string] =
  ## type keyword at start of s
  if s.startsWith("DOUBLEPRECISION"): return (tkReal8, 0, s[15..^1])
  if s.startsWith("REAL*8"): return (tkReal8, 0, s[6..^1])
  if s.startsWith("REAL*4"): return (tkReal4, 0, s[6..^1])
  if s.startsWith("INTEGER*8"): return (tkInt8, 0, s[9..^1])
  if s.startsWith("INTEGER*4"): return (tkInt4, 0, s[9..^1])
  if s.startsWith("LOGICAL"): return (tkLogical, 0, s[7..^1])
  if s.startsWith("CHARACTER*"):
    var j = 10
    while j < s.len and s[j] in Digits: inc j
    return (tkChar, parseInt(s[10 ..< j]), s[j..^1])
  if s.startsWith("INTEGER"): return (tkInt4, 0, s[7..^1])
  if s.startsWith("REAL"): return (tkReal4, 0, s[4..^1])
  (tkNone, 0, s)

proc parseIoItems(s: string, line: int): seq[IoItem] =
  for part in splitTop(s):
    if part.len == 0: continue
    if part[0] == '(' and matchParen(part, 0) == part.len-1 and
        topIndex(part[1..^2], '=') >= 0:
      # implied DO: (items, V=lo,hi[,step])
      let inner = splitTop(part[1..^2])
      var k = -1
      for i, x in inner:
        if topIndex(x, '=') >= 0: k = i; break
      let eq = inner[k].find('=')
      var it = IoItem(kind: ioDo, doVar: inner[k][0 ..< eq],
                      lo: parseExprStr(inner[k][eq+1..^1], line),
                      hi: parseExprStr(inner[k+1], line))
      if k+2 < inner.len: it.step = parseExprStr(inner[k+2], line)
      it.items = parseIoItems(inner[0 ..< k].join(","), line)
      result.add it
    else:
      result.add IoItem(kind: ioExpr, e: parseExprStr(part, line))

# ---------------------------------------------------------------- statements

proc isAssignment(s: string): bool =
  var j = 0
  if j >= s.len or s[j] notin Letters: return false
  while j < s.len and s[j] in IdentChars: inc j
  if j < s.len and s[j] == '(':
    let e = matchParen(s, j)
    if e < 0: return false
    j = e+1
  j < s.len and s[j] == '=' and (j+1 >= s.len or s[j+1] != '=')

proc parseSimple(r: RawStmt, s: string): Stmt =
  ## executable statement other than block-IF/DO structure statements
  let line = r.line
  template mk(k: StmtKind): Stmt = Stmt(kind: k, line: line, label: r.label, comments: r.comments)
  if s.startsWith("IF(") and isAssignment(s) == false:
    let e = matchParen(s, 2)
    let rest = s[e+1..^1]
    result = mk(sIf)
    result.cond = parseExprStr(s[3 ..< e], line)
    var r2 = r
    r2.label = 0
    r2.comments = @[]
    result.then = parseSimple(r2, rest)
    return
  if isAssignment(s):
    let k = topIndex(s, '=')
    result = mk(sAssign)
    result.lhs = parseExprStr(s[0 ..< k], line)
    result.rhs = parseExprStr(s[k+1..^1], line)
    return
  if s.startsWith("GOTO("):
    let e = matchParen(s, 4)
    result = mk(sCGoto)
    for x in splitTop(s[5 ..< e]): result.targets.add parseInt(x)
    var rest = s[e+1..^1]
    if rest.startsWith(","): rest = rest[1..^1]
    result.sel = parseExprStr(rest, line)
    return
  if s.startsWith("GOTO"):
    result = mk(sGoto)
    result.target = parseInt(s[4..^1])
    return
  if s.startsWith("CALL"):
    result = mk(sCall)
    let rest = s[4..^1]
    let k = rest.find('(')
    if k < 0:
      result.callee = rest
    else:
      result.callee = rest[0 ..< k]
      let e = matchParen(rest, k)
      if e > k+1:
        for a in splitTop(rest[k+1 ..< e]): result.cargs.add parseExprStr(a, line)
    return
  if s == "RETURN": return mk(sReturn)
  if s == "STOP" or (s.startsWith("STOP") and s[4..^1].allCharsInSet(Digits)): return mk(sStop)
  if s == "CONTINUE": return mk(sContinue)
  if s.startsWith("WRITE("):
    let e = matchParen(s, 5)
    let ctl = splitTop(s[6 ..< e])
    if ctl.len != 2 or ctl[0] != "6": perr(line, "unsupported WRITE control: " & s)
    result = mk(sWrite)
    if ctl[1] == "*": result.fmtLabel = 0
    elif ctl[1][0] == '\'': result.fmtInline = ctl[1][1..^2].replace("''", "'")
    else: result.fmtLabel = parseInt(ctl[1])
    result.witems = parseIoItems(s[e+1..^1], line)
    return
  if s.startsWith("READ("):
    let e = matchParen(s, 4)
    if s[5 ..< e] != "5,*": perr(line, "unsupported READ control: " & s)
    result = mk(sRead)
    result.ritems = parseIoItems(s[e+1..^1], line)
    return
  if s.startsWith("FORMAT("):
    result = mk(sFormat)
    result.ftext = r.text.strip   # original text keeps blanks in strings
    result.ftext = result.ftext[result.ftext.find('(') .. ^1]
    return
  perr(line, "unrecognised statement: " & r.text)

proc parseDataStmt(s: string, line: int): seq[DataSpec] =
  ## s without the DATA keyword: nlist/clist/[,]nlist/clist/...
  var i = 0
  while i < s.len:
    if s[i] == ',': inc i
    let a = s.find('/', i)
    let b = s.find('/', a+1)
    if a < 0 or b < 0: perr(line, "bad DATA")
    var spec: DataSpec
    for n in splitTop(s[i ..< a]): spec.names.add parseExprStr(n, line)
    for v in splitTop(s[a+1 ..< b]):
      let st = topIndex(v, '*')
      if st > 0 and v[0 ..< st].allCharsInSet(Digits):
        spec.values.add (parseInt(v[0 ..< st]), parseExprStr(v[st+1..^1], line))
      else:
        spec.values.add (1, parseExprStr(v, line))
    result.add spec
    i = b+1

proc parseImplicit(u: Unit, s: string, line: int) =
  if s == "NONE":
    u.implicitNone = true
    return
  let k = s.find('(')
  let (t, _, rest) = parseTypeSpec(s[0 ..< k])
  if t == tkNone or rest.len > 0: perr(line, "bad IMPLICIT " & s)
  for rng in splitTop(s[k+1 ..< matchParen(s, k)]):
    if rng.len == 1: u.implicits.add (rng[0], rng[0], t)
    elif rng.len == 3 and rng[1] == '-': u.implicits.add (rng[0], rng[2], t)
    else: perr(line, "bad IMPLICIT range " & rng)

proc parseHeader(s: string, line: int): Unit =
  var s = s
  var rt = tkNone
  var kind: UnitKind
  if s.startsWith("PROGRAM"):
    return Unit(kind: ukProgram, name: s[7..^1], line: line)
  let (t, _, rest) = parseTypeSpec(s)
  if t != tkNone:
    rt = t
    s = rest
  if s.startsWith("SUBROUTINE"):
    kind = ukSubroutine; s = s[10..^1]
  elif s.startsWith("FUNCTION"):
    kind = ukFunction; s = s[8..^1]
  else:
    return nil
  result = Unit(kind: kind, resType: rt, line: line)
  let k = s.find('(')
  if k < 0: result.name = s
  else:
    result.name = s[0 ..< k]
    let e = matchParen(s, k)
    if e > k+1: result.params = splitTop(s[k+1 ..< e])

type Frame = object
  st: Stmt       ## sIfBlock or sDo; nil for the unit body
  list: ptr seq[Stmt]

proc parseUnits*(raws: seq[RawStmt]): seq[Unit] =
  var i = 0
  var pending: seq[string]
  while i < raws.len:
    let r = raws[i]
    let u = parseHeader(r.norm, r.line)
    if u == nil: perr(r.line, "expected unit header, got: " & r.text)
    u.headComments = r.comments
    inc i
    var stack: seq[Frame] = @[Frame(st: nil, list: addr u.body)]
    var inExec = false
    while true:
      if i >= raws.len: perr(r.line, "missing END for " & u.name)
      let r = raws[i]
      inc i
      let s = r.norm
      if s == "END":
        if stack.len != 1: perr(r.line, "unclosed IF/DO at END of " & u.name)
        break
      # --- declarations
      if not inExec:
        var handled = true
        if s.startsWith("IMPLICIT"): parseImplicit(u, s[8..^1], r.line)
        elif s.startsWith("DIMENSION"): u.dims.add parseEntities(s[9..^1], r.line)
        elif s.startsWith("COMMON/"):
          var j = 6
          while j < s.len:
            let e = s.find('/', j+1)
            let blk = s[j+1 ..< e]
            var nx = s.find('/', e+1)
            if nx < 0: nx = s.len
            var lst = s[e+1 ..< nx]
            if lst.endsWith(","): lst.setLen(lst.len-1)
            u.commons.add (blk, parseEntities(lst, r.line))
            j = nx
        elif s.startsWith("PARAMETER("):
          for a in splitTop(s[10 ..< matchParen(s, 9)]):
            let k = a.find('=')
            u.parameters.add (a[0 ..< k], parseExprStr(a[k+1..^1], r.line))
        elif s.startsWith("SAVE"):
          for n in splitTop(s[4..^1]): u.saves.add n
        elif s.startsWith("DATA") and not isAssignment(s):
          u.datas.add parseDataStmt(s[4..^1], r.line)
        elif s.startsWith("EXTERNAL"): discard
        elif not isAssignment(s) and parseTypeSpec(s).t != tkNone:
          let (t, clen, rest) = parseTypeSpec(s)
          u.typeDecls.add (t, clen, parseEntities(rest, r.line))
        else:
          handled = false
        if handled:
          if r.label != 0: perr(r.line, "label on declaration")
          continue
      if s.startsWith("DATA") and not isAssignment(s):
        # DATA after executable statements (legal F77)
        u.datas.add parseDataStmt(s[4..^1], r.line)
        continue
      if s.startsWith("ENTRY"):
        let k = s.find('(')
        u.entries.add (s[5 ..< k], splitTop(s[k+1 ..< matchParen(s, k)]))
        continue
      inExec = true
      template top: untyped = stack[^1]
      # --- structure statements
      if s.startsWith("IF(") and s.endsWith("THEN") and matchParen(s, 2) == s.len-5:
        let st = Stmt(kind: sIfBlock, label: r.label, line: r.line, comments: r.comments)
        st.conds.add parseExprStr(s[3 ..< s.len-5], r.line)
        st.bodies.add @[]
        top.list[].add st
        stack.add Frame(st: st, list: addr st.bodies[^1])
        continue
      if s.startsWith("ELSEIF(") and s.endsWith("THEN"):
        if top.st == nil or top.st.kind != sIfBlock or top.st.hasElse: perr(r.line, "ELSE IF outside IF")
        if r.label != 0: perr(r.line, "label on ELSE IF")
        let st = top.st
        st.conds.add parseExprStr(s[7 ..< s.len-5], r.line)
        st.bodies.add @[]
        top.list = addr st.bodies[^1]
        continue
      if s == "ELSE":
        if top.st == nil or top.st.kind != sIfBlock or top.st.hasElse: perr(r.line, "ELSE outside IF")
        if r.label != 0: perr(r.line, "label on ELSE")
        let st = top.st
        st.hasElse = true
        st.bodies.add @[]
        top.list = addr st.bodies[^1]
        continue
      if s == "ENDIF":
        if top.st == nil or top.st.kind != sIfBlock: perr(r.line, "ENDIF outside IF")
        discard stack.pop
        if r.label != 0:
          # target of a jump to END IF: becomes a CONTINUE after the IF
          top.list[].add Stmt(kind: sContinue, label: r.label, line: r.line)
        continue
      if s.startsWith("DO") and s.len > 2 and s[2] in Digits:
        var j = 2
        while s[j] in Digits: inc j
        let eq = s.find('=')
        if eq > j and topIndex(s[eq+1..^1], ',') >= 0:
          let st = Stmt(kind: sDo, label: r.label, line: r.line, comments: r.comments)
          st.dlabel = parseInt(s[2 ..< j])
          var v = s[j ..< eq]
          if v.startsWith(","): v = v[1..^1]
          st.dvar = v
          let parts = splitTop(s[eq+1..^1])
          st.dlo = parseExprStr(parts[0], r.line)
          st.dhi = parseExprStr(parts[1], r.line)
          if parts.len > 2: st.dstep = parseExprStr(parts[2], r.line)
          top.list[].add st
          stack.add Frame(st: st, list: addr st.body)
          continue
      # --- simple statement
      let st = parseSimple(r, s)
      top.list[].add st
      # close DO loops terminated by this label
      if r.label != 0:
        while stack.len > 1 and stack[^1].st.kind == sDo and stack[^1].st.dlabel == r.label:
          discard stack.pop
        for f in stack:
          if f.st != nil and f.st.kind == sDo and f.st.dlabel == r.label:
            perr(r.line, "DO terminator " & $r.label & " closes across IF block")
    result.add u
    discard pending
