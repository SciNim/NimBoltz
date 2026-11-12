## The drift run (MONTE) with walkers:
## many independent electrons (`walkers.nim`) instead of one long
## trajectory, on CPU threads (`CpuWalkers`) or a GPU (`gpuwalk.GpuWalkers`),
## for the gas at 0 K or with thermal motion. The same physics, other estimators (`ensemble.nim`):
## the results agree with Magboltz's within errors, not number for number.
##
## A run: the walkers start where Magboltz's electron starts and run until
## they have forgotten it (`warmUp`), then the sums restart; the main run
## stops at checkpoints, where the walkers' positions are kept (for the
## diffusion, `ensemble.nim`), and every collision of the main run is
## counted per process and energy bin.
##
## How many walkers: the time to forget the start is a few energy
## relaxation times, and the main run of each walker must last at least as
## long, else the walkers spread apart by their differences in energy
## (ballistically, not by diffusion) and the diffusion comes out too large.
## So noble gases, which lose energy slowly, get longer walkers: on CPUs
## fewer of them, on a GPU as many, which runs them in the same time (and
## then more collisions than asked for).
##
## ELIMIT on walkers (`energyLimitExceeded`): whether, once warm, they get
## a collision in the last energy bin among ELIMIT's number of them.

import std/math
import units, vec3, tables, fields, freeflight, drift, walkers, ensemble

type
  WarmUp* = object
    ## how the walkers forgot their start
    collisions*: int     ## per walker
    converged*: bool     ## false: stopped at the limit
    energy*: eV          ## mean energy of the last two chunks
    previous*: eV
    exceeded*: bool      ## `exceeding`: too many collisions in the last energy bin

  Chunk = object
    e, ee: float64       ## mean energy and its error (eV)
    v, ve: float64       ## drift velocity along E and its error (m/ps)

func chunk[T, R](ws: openArray[Walker[T, R]], blocks = 10): Chunk =
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

proc drop*[B, T, R](b: var B, ws: var seq[Walker[T, R]], n: int) =
  ## keeps the first n walkers
  mixin shrink
  ws.setLen(n)
  b.shrink(n)

proc lastBinShare*[B](b: B): float64 =
  ## share of the collisions counted since `clearTally` in the last energy
  ## bin
  mixin readTally
  var t {.noinit.}: Tallies
  b.readTally(t)
  t.spectrum[nEnergyBins] / max(sum(t.spectrum), 1.0)

proc warmUp*[B, T, R](b: var B, ws: var seq[Walker[T, R]], budget: int, first = 1000'i32,
                   exceeding = 0.0): WarmUp =
  ## Runs the walkers in chunks of doubling length until the mean energy and
  ## the drift velocity of a chunk agree with those of the chunk before
  ## (within 2 errors); then clears their sums. At most budget / 2 /
  ## `minWalkers` collisions per walker (then not converged). On CPUs, where
  ## time goes with all collisions, a chunk spends at most budget / 2 of
  ## them: half of the walkers are dropped for longer chunks (down to
  ## `minWalkers`); on a GPU, where time goes with the collisions per
  ## walker, all walkers stay. `exceeding` > 0: the collisions are counted
  ## (after a `clearTally`), and the warm-up stops (`exceeded`) when more
  ## than this share of them is in the last energy bin.
  mixin walk, read, write, dropsWalkers
  proc run(b: var B, ws: var seq[Walker[T, R]], c: int32): Chunk =
    b.walk(c, counting = exceeding > 0.0)
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
    if exceeding > 0.0 and b.lastBinShare > exceeding:
      result.exceeded = true
      return
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

proc energyLimitExceeded*[B, T, R](b: var B, ws: var seq[Walker[T, R]], budget, sample: int): WarmUp =
  ## ELIMIT on walkers: whether an electron reaches the last energy bin
  ## (`result.exceeded`) in `sample` collisions of all walkers together,
  ## once they have warmed up (`warmUp`). Magboltz's ELIMIT counts its one
  ## electron from the start; the walkers' starts are many, and would
  ## raise EFINAL for nothing. Stops early when the last bin gets more than
  ## 1000 / `sample` of the collisions so far (exceeded for sure). Leaves
  ## the walkers warm, their sums cleared.
  mixin walk, read, write, clearTally
  const probe = 250'i32  # a range far too small shows at once
  let share = 1000.0 / float64(sample)
    # as MONTE's limit: one walker running away gives a burst of
    # collisions in the last bin, which a warm-up of many more collisions
    # than `sample` catches more often than ELIMIT would
  b.clearTally()
  b.walk(probe)
  if b.lastBinShare > share: return WarmUp(collisions: probe, exceeded: true)
  b.clearTally()
  result = b.warmUp(ws, budget, exceeding = share)
  result.collisions = result.collisions + probe
  if result.exceeded: return
  b.clearTally()
  b.walk(int32(max(sample div ws.len, 1)))
  result.exceeded = b.lastBinShare > 0.0
  b.read(ws)
  for w in ws.mitems: w.reset()
  b.write(ws)

type
  WalkerRun* = object
    ## the outcome of the main run
    transport*: Transport
    tallies*: Tallies
    st*: ps        ## time of all walkers
    nulls*: int64
    tmax*: ps
    aborted*: bool ## stopped early: the last energy bin got more than
                   ## `abortAbove` of the collisions (no transport then)

const maxChunks* = 64
  ## checkpoints of the main run (`ensemble.Checkpoints`) at most

func chunksOf*(perWalker, warmUp: int): int =
  ## chunks of the main run between checkpoints: each at least as long as
  ## the warm-up (longer than the correlation of the velocities), at least
  ## two, at most `maxChunks`
  clamp(perWalker div max(warmUp, 1), 2, maxChunks)

proc driftRun*[G: static FieldGeometry; B, T, R](b: var B, ws: var seq[Walker[T, R]],
                                              perWalker: int32, chunks: int,
                                              btheta: float64,
                                              abortAbove = 0.0, anyBefore = 0): WalkerRun =
  ## `perWalker` collisions of every walker (after `warmUp`), in `chunks`
  ## chunks with checkpoints in between, counted. With `abortAbove` > 0,
  ## stops after a chunk with more than this share of its collisions in the
  ## last energy bin, or, while the walkers have made fewer than `anyBefore`
  ## collisions together before it (and in the first chunk), with any
  ## (`aborted`: the energy range is too small, no need to run on). The
  ## second is ELIMIT's test continued: a tail that fills slowly can pass
  ## it in the warm-up.
  mixin walk, read, write, clearTally, readTally
  b.clearTally()
  var cp = initCheckpoints(ws)
  var done = 0'i32
  var last, all = 0.0  # collisions in the last energy bin and in all, so far
  for j in 1 .. chunks:
    let upTo = int32(int64(perWalker) * j div chunks)
    b.walk(upTo - done)
    done = upTo
    if abortAbove > 0.0 and j < chunks:
      var t {.noinit.}: Tallies
      b.readTally(t)
      let (l, a) = (t.spectrum[nEnergyBins], sum(t.spectrum))
      if (l - last) / max(a - all, 1.0) > abortAbove or
         (l > 0.0 and (j == 1 or all < float64(anyBefore))):
        result.aborted = true
        b.read(ws)
        return
      (last, all) = (l, a)
    b.read(ws)
    cp.add ws
  b.readTally(result.tallies)
  result.transport = transport[G, T, R](ws, cp, btheta)
  for w in ws:
    result.st = result.st + w.st
    result.nulls = result.nulls + int64(w.nulls)
    if w.tmax > result.tmax: result.tmax = w.tmax
