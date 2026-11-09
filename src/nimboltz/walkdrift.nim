## The drift run (MONTE) with walkers:
## many independent electrons (`walkers.nim`) instead of one long
## trajectory, on CPU threads (`CpuWalkers`) or a GPU (`gpuwalk.GpuWalkers`),
## for the gas at 0 K or with thermal motion. The same physics, other estimators (`ensemble.nim`):
## the results agree with Magboltz's within errors, not number for number.
##
## A run: the walkers start where Magboltz's electron starts and run until
## they have forgotten it (`warmUp`), then the sums restart; halfway through
## the main run the walkers mark their positions (for the diffusion), and
## every collision of the main run is counted per process and energy bin.
##
## How many walkers: the time to forget the start is a few energy
## relaxation times, and the main run of each walker must last at least as
## long, else the walkers spread apart by their differences in energy
## (ballistically, not by diffusion) and the diffusion comes out too large.
## So noble gases, which lose energy slowly, get longer walkers: on CPUs
## fewer of them, on a GPU as many, which runs them in the same time (and
## then more collisions than asked for).

import std/math
import units, vec3, tables, fields, freeflight, drift, walkers, ensemble

type
  WarmUp* = object
    ## how the walkers forgot their start
    collisions*: int     ## per walker
    converged*: bool     ## false: stopped at the limit
    energy*: eV          ## mean energy of the last two chunks
    previous*: eV

  Chunk = object
    e, ee: float64       ## mean energy and its error (eV)
    v, ve: float64       ## drift velocity along E and its error (m/ps)

func chunk[T](ws: openArray[Walker[T]], blocks = 10): Chunk =
  ## mean energy and drift velocity of the walkers' sums, with errors from
  ## the scatter of blocks
  var es, vs: seq[float64]
  for b in 0 ..< blocks:
    var e2 = 0.0.eV•ps
    var z = 0.0.m
    var t = 0.0.ps
    for k in b * ws.len div blocks ..< (b + 1) * ws.len div blocks:
      e2 = e2 + ws[k].e2
      z = z + ws[k].pos.z
      t = t + ws[k].st
    es.add float64(e2 / t)
    vs.add float64(z / t)
  func meanErr(xs: seq[float64]): (float64, float64) =
    let m = sum(xs) / float64(xs.len)
    var s = 0.0
    for x in xs: s = s + (x - m) * (x - m)
    (m, sqrt(s / float64(xs.len - 1) / float64(xs.len)))
  (result.e, result.ee) = meanErr(es)
  (result.v, result.ve) = meanErr(vs)

const minWalkers* = 100  ## 10 blocks of 10 for the errors

proc drop*[B, T](b: var B, ws: var seq[Walker[T]], n: int) =
  ## keeps the first n walkers
  mixin shrink
  ws.setLen(n)
  b.shrink(n)

proc warmUp*[B, T](b: var B, ws: var seq[Walker[T]], budget: int, first = 1000'i32): WarmUp =
  ## Runs the walkers in chunks of doubling length until the mean energy and
  ## the drift velocity of a chunk agree with those of the chunk before
  ## (within 2 errors); then clears their sums. At most budget / 2 /
  ## `minWalkers` collisions per walker (then not converged). On CPUs, where
  ## time goes with all collisions, a chunk spends at most budget / 2 of
  ## them: half of the walkers are dropped for longer chunks (down to
  ## `minWalkers`); on a GPU, where time goes with the collisions per
  ## walker, all walkers stay.
  mixin walk, read, write, dropsWalkers
  proc run(b: var B, ws: var seq[Walker[T]], c: int32): Chunk =
    b.walk(c, counting = false)
    b.read(ws)
    result = chunk(ws)
    for w in ws.mitems: w.reset()
    b.write(ws)
  let cap = budget div 2 div minWalkers
  template fit(c: int32) =
    while b.dropsWalkers and ws.len div 2 >= minWalkers and ws.len * int(c) > budget div 2:
      b.drop(ws, ws.len div 2)
  var c = first
  fit(c)
  var prev = b.run(ws, c)
  result.collisions = c
  while true:
    c = 2'i32 * c
    fit(c)
    if result.collisions + int(c) > cap: return
    let cur = b.run(ws, c)
    result.collisions = result.collisions + c
    result.energy = cur.e.eV
    result.previous = prev.e.eV
    if abs(cur.e - prev.e) <= 2.0 * sqrt(cur.ee ^ 2 + prev.ee ^ 2) and
       abs(cur.v - prev.v) <= 2.0 * sqrt(cur.ve ^ 2 + prev.ve ^ 2):
      result.converged = true
      return
    prev = cur

type
  WalkerRun* = object
    ## the outcome of the main run
    transport*: Transport
    tallies*: Tallies
    st*: ps        ## time of all walkers
    nulls*: int64
    tmax*: ps

proc driftRun*[G: static FieldGeometry; B, T](b: var B, ws: var seq[Walker[T]],
                                              perWalker: int32, btheta: float64): WalkerRun =
  ## `perWalker` collisions of every walker (after `warmUp`), marked
  ## halfway, counted
  mixin walk, read, write, clearTally, readTally
  b.clearTally()
  let half = perWalker div 2
  b.walk(half)
  b.read(ws)
  for w in ws.mitems: w.mark()
  b.write(ws)
  b.walk(perWalker - half)
  b.read(ws)
  b.readTally(result.tallies)
  result.transport = transport[G, T](ws, btheta)
  for w in ws:
    result.st = result.st + w.st
    result.nulls = result.nulls + int64(w.nulls)
    if w.tmax > result.tmax: result.tmax = w.tmax
