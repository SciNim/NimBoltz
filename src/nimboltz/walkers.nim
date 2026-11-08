## Independent electrons drifting in a gas at 0 K with E only (B = 0).
##
## MONTE's collision loop (`drift.runSegment` for fgElectric and `ColdGas`)
## without its bookkeeping: each walker is one electron with its own random
## numbers (`Ranmar`), and keeps only the sums that give the drift velocity
## and the mean energy. Many walkers give the same estimates as one long
## trajectory, not the same numbers; they can run side by side, e.g. one
## per GPU thread. No Penning transfers (IPEN = 0).
##
## Plain Nim the GPU compiler takes too: no `break`, `continue` or `case`,
## no references, `int32` counters.

import units, rm48, tables, vec3, fields, collisions, freeflight

type
  Walker* = object
    el*: Electron        ## state after the last collision
    target*: ColdGas     ## null-collision bounds, raised as needed
    rng*: Ranmar
    pos*: Vector[m]      ## position
    st*: ps              ## time
    e2*: eV•ps           ## ∫ e dt
    collisions*: int32   ## real collisions
    nulls*: int32        ## null collisions

proc initWalker*(start: Electron, target: ColdGas, seed: int): Walker =
  Walker(el: start, target: target, rng: initRanmar(seed))

proc walk*(w: var Walker, f: Field, tab: MixtureTable, estep: eV, collisions: int32) =
  ## `collisions` real collisions, with the random numbers drawn in the
  ## order of MONTE
  for _ in 0'i32 ..< collisions:
    let h = nextHit[fgElectric](f, tab, w.target, estep, w.el, w.rng, w.nulls)
    w.e2 = w.e2 + energyIntegral[fgElectric](f, w.el, h.fl)
    discard advance[fgElectric](f, w.el, h.fl, w.pos)
    w.st = w.st + h.fl.t
    let r2 = drand48(w.rng)
    let i = findProcess(tab, h.ie, r2)
    let ei = energyLoss(tab, i, h.e, w.rng, moderation = true)
    let s = scatter(tab, i, h.ie, h.e, ei, h.din, w.rng, r31First = true)
    w.el = leave[fgElectric](w.target, f, s, h)
  w.collisions = w.collisions + collisions

func driftVelocity*(ws: openArray[Walker]): μm•ns⁻¹ =
  ## Σ z / Σ t over all walkers
  var z = 0.0.m
  var t = 0.0.ps
  for w in ws:
    z = z + w.pos.z
    t = t + w.st
  (z / t).to(μm•ns⁻¹)

func meanEnergy*(ws: openArray[Walker]): eV =
  ## Σ ∫ e dt / Σ t over all walkers
  var e2 = 0.0.eV•ps
  var t = 0.0.ps
  for w in ws:
    e2 = e2 + w.e2
    t = t + w.st
  e2 / t
