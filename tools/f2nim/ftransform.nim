## AST transformations applied before structuring.

import std/[tables, strformat, sequtils]
import fast, fsema

proc clone*(st: Stmt): Stmt =
  ## copy of a statement tree (expressions are shared), labels dropped
  result = Stmt(kind: st.kind, line: st.line, comments: st.comments)
  case st.kind
  of sAssign:
    result.lhs = st.lhs; result.rhs = st.rhs
  of sIf:
    result.cond = st.cond; result.then = clone(st.then)
  of sIfBlock:
    result.conds = st.conds; result.hasElse = st.hasElse
    for b in st.bodies:
      var nb: seq[Stmt]
      for s in b: nb.add clone(s)
      result.bodies.add nb
  of sDo:
    result.dvar = st.dvar; result.dlo = st.dlo; result.dhi = st.dhi
    result.dstep = st.dstep; result.dlabel = st.dlabel
    for s in st.body: result.body.add clone(s)
  of sGoto: result.target = st.target
  of sCGoto:
    result.targets = st.targets; result.sel = st.sel
  of sCall:
    result.callee = st.callee; result.cargs = st.cargs
  of sWrite:
    result.fmtLabel = st.fmtLabel; result.fmtInline = st.fmtInline
    result.witems = st.witems
  of sRead: result.ritems = st.ritems
  of sFormat: result.ftext = st.ftext
  of sReturn, sStop, sContinue: discard

proc hasLabels(st: Stmt): bool =
  if st.label != 0: return true
  case st.kind
  of sIf: hasLabels(st.then)
  of sIfBlock:
    for b in st.bodies:
      for s in b:
        if hasLabels(s): return true
    false
  of sDo:
    for s in st.body:
      if hasLabels(s): return true
    false
  else: false

type
  ListRef = ptr seq[Stmt]
  Level = tuple[list: ListRef, idx: int, isDo: bool]   ## list + index of the containing stmt

proc hoistTails*(u: Unit): int =
  ## Removes jumps into a block (illegal F77, accepted by gfortran).
  ##
  ## The jumped-to code TAIL (from the target label to the first unconditional
  ## GO TO / RETURN / STOP in the same list) is never left by falling through,
  ## so it can be moved out of its IF blocks without changing the semantics:
  ##
  ##   S (contains: ... t TAIL ...)   ->   S (contains: ... GO TO t ...)
  ##                                       ...
  ##                                       RETURN
  ##                                     t TAIL
  ##
  ## TAIL goes to the end of the innermost list enclosing both the jump and
  ## the target, so that all jumps to it are forward jumps.
  ## Returns the number of hoisted blocks.
  var nextLabel = 99001
  while true:
    var labels: Table[int, seq[Level]]   # label -> chain, last = its own list
    var bad: seq[tuple[chain: seq[ListRef], target: int]]
    proc scan(list: ListRef, chain: seq[Level]) =
      for i in 0 ..< list[].len:
        let st = list[][i]
        let here = chain & @[(list, i, false)]
        if st.label != 0: labels[st.label] = here
        case st.kind
        of sIfBlock:
          for b in 0 ..< st.bodies.len: scan(addr st.bodies[b], chain & @[(list, i, false)])
        of sDo: scan(addr st.body, chain & @[(list, i, true)])
        else: discard
    proc findBad(list: ListRef, enclosing: seq[ListRef]) =
      let encl = enclosing & list
      for i in 0 ..< list[].len:
        let st = list[][i]
        var ts: seq[int]
        if st.kind == sGoto: ts.add st.target
        elif st.kind == sCGoto: ts.add st.targets
        elif st.kind == sIf and st.then.kind == sGoto: ts.add st.then.target
        for t in ts:
          if t notin labels: serr(&"{u.name}: jump to unknown label {t}")
          if labels[t][^1].list notin encl: bad.add (encl, t)
        case st.kind
        of sIfBlock:
          for b in 0 ..< st.bodies.len: findBad(addr st.bodies[b], encl)
        of sDo: findBad(addr st.body, encl)
        else: discard
    scan(addr u.body, @[])
    findBad(addr u.body, @[])
    if bad.len == 0: return
    if result > 50: serr(&"{u.name}: hoisting does not terminate")
    let (encl, t) = bad[0]
    let chain = labels[t]
    # innermost common list
    var c = -1
    for lv in 0 ..< chain.len:
      if chain[lv].list in encl: c = lv
    if c < 0: serr(&"{u.name}: no common list for jump to {t}")
    for lv in c ..< chain.len-1:
      if chain[lv].isDo: serr(&"{u.name}: jump into DO loop at label {t}")
    let tl = chain[^1].list
    let j = chain[^1].idx
    var k = j
    while true:
      if k >= tl[].len: serr(&"{u.name}: code at label {t} falls off its block")
      if tl[][k].kind in {sGoto, sReturn, sStop}: break
      inc k
    var tail = tl[][j .. k]
    let gotoT = Stmt(kind: sGoto, target: t, line: tail[0].line)
    tl[].delete(j .. k)
    tl[].insert(gotoT, j)
    # TAIL goes where it can only be reached by jumps: after the first
    # unconditional jump in the common list that follows all jumps into it;
    # else at the end (after a RETURN, or behind a skip label)
    let cl = chain[c].list
    var tailLabels: seq[int]
    for x in tail:
      proc labs(st: Stmt) =
        if st.label != 0: tailLabels.add st.label
        case st.kind
        of sIf: labs(st.then)
        of sIfBlock:
          for b in st.bodies:
            for y in b: labs(y)
        of sDo:
          for y in st.body: labs(y)
        else: discard
      labs(x)
    proc jumpsInto(st: Stmt): bool =
      case st.kind
      of sGoto: st.target in tailLabels
      of sCGoto: st.targets.anyIt(it in tailLabels)
      of sIf: jumpsInto(st.then)
      of sIfBlock: st.bodies.anyIt(it.anyIt(jumpsInto(it)))
      of sDo: st.body.anyIt(jumpsInto(it))
      else: false
    var lastAnchor = -1
    for i, x in cl[]:
      if jumpsInto(x): lastAnchor = i
    var p = -1
    for i in max(lastAnchor, 0) ..< cl[].len:
      if cl[][i].kind in {sGoto, sReturn, sStop}:
        p = i
        break
    let last = cl[][^1]
    if p >= 0 and p < cl[].len - 1:
      for m, x in tail: cl[].insert(x, p + 1 + m)
    elif last.kind in {sGoto, sReturn, sStop}:
      cl[].add tail
    elif c == 0:
      cl[].add Stmt(kind: (if u.kind == ukProgram: sStop else: sReturn), line: last.line)
      cl[].add tail
    else:
      let skip = nextLabel
      inc nextLabel
      cl[].add Stmt(kind: sGoto, target: skip, line: tail[0].line)
      cl[].add tail
      cl[].add Stmt(kind: sContinue, label: skip, line: tail[^1].line)
    inc result

proc splitLabel*(u: Unit, x: int): int =
  ## Node splitting for irreducible control flow: every `GO TO x` is replaced
  ## by a copy of the code at label x (up to the next label y, which it falls
  ## into) followed by `GO TO y`. Returns the number of replaced jumps.
  var found: ListRef
  var j = -1
  proc find(list: ListRef) =
    for i in 0 ..< list[].len:
      let st = list[][i]
      if st.label == x: (found = list; j = i)
      case st.kind
      of sIfBlock:
        for b in 0 ..< st.bodies.len: find(addr st.bodies[b])
      of sDo: find(addr st.body)
      else: discard
  find(addr u.body)
  if j < 0: serr(&"{u.name}: split label {x} not found")
  var k = j + 1
  while k < found[].len and found[][k].label == 0: inc k
  if k >= found[].len: serr(&"{u.name}: code at split label {x} is not followed by a label")
  let seg = found[][j ..< k]
  for i, st in seg:
    if (i > 0 and hasLabels(st)) or st.kind in {sGoto, sReturn, sStop}:
      serr(&"{u.name}: code at split label {x} is not straight-line")
  let y = found[][k].label
  proc copyFor(line: int): seq[Stmt] =
    for st in seg: result.add clone(st)
    result.add Stmt(kind: sGoto, target: y, line: line)
  var n = 0
  proc rewrite(list: ListRef) =
    var i = 0
    while i < list[].len:
      let st = list[][i]
      if st.kind == sGoto and st.target == x:
        let c = copyFor(st.line)
        c[0].label = st.label
        list[].delete(i)
        for m, s in c: list[].insert(s, i + m)
        i += c.len
        inc n
        continue
      if st.kind == sIf and st.then.kind == sGoto and st.then.target == x:
        let blk = Stmt(kind: sIfBlock, label: st.label, line: st.line, comments: st.comments)
        blk.conds = @[st.cond]
        blk.bodies = @[copyFor(st.line)]
        list[][i] = blk
        inc n
      case list[][i].kind
      of sIfBlock:
        for b in 0 ..< list[][i].bodies.len: rewrite(addr list[][i].bodies[b])
      of sDo: rewrite(addr list[][i].body)
      else: discard
      inc i
  rewrite(addr u.body)
  n
