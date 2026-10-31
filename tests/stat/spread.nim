## Spread of Magboltz results over independent random streams.
##
##   spread DECK OUT...        one line per result of the deck (TSV):
##                             deck, key, runs, mean, empirical σ (% of
##                             |mean|), printed σ (rms, %), ratio
##   spread --summary TSV...   ratio empirical/printed by kind of result
##
## The outputs are runs of the same deck with different seeds
## (NIMBOLTZ_SEED). Results missing from a run (or an aborted run) are
## skipped; a result needs at least 4 runs.

import std/[strutils, math, tables, strformat, algorithm, os]
import magout

proc perDeck(deck: string, files: seq[string]) =
  var vals: OrderedTable[string, seq[(float, float)]]
  for f in files:
    let p = parse(f)
    if p.aborted: continue
    for r in p.results:
      vals.mgetOrPut(r.key, @[]).add (r.value, r.pct)
  for key, v in vals:
    if v.len < 4: continue
    var mean, pct2 = 0.0
    for (x, p) in v: (mean += x; pct2 += p * p)
    mean /= v.len.float
    var ss = 0.0
    for (x, _) in v: ss += (x - mean) ^ 2
    let sd = sqrt(ss / float(v.len - 1))
    let rel = if mean != 0: 100 * sd / abs(mean) else: (if sd == 0: 0.0 else: Inf)
    let printed = sqrt(pct2 / v.len.float)
    let ratio = if printed > 0: rel / printed else: (if rel == 0: 1.0 else: Inf)
    echo &"{deck}\t{key}\t{v.len}\t{mean:.6g}\t{rel:.4g}\t{printed:.4g}\t{ratio:.4g}"

proc kindOf(key: string): string =
  ## coarse kind of a result, from its label
  const kinds = [
    ("DETAILED", "per-level"), ("ELVL", "per-level"), ("ELOSS", "per-level"),
    ("SST", "SST"), ("TOF", "TOF"), ("PULSED", "PT"), ("VD=", "SST/TOF VD,WS"),
    ("WS=", "SST/TOF VD,WS"), ("DL=", "SST/TOF DL,DT"), ("ALPHA", "alpha/att"),
    ("ATT", "alpha/att"), ("DRIFT VEL", "drift"), ("MEAN ELECTRON", "mean energy"),
    ("LONGITUDINAL", "diffusion L"), ("TRANSVERSE", "diffusion T"),
    ("DIFFUSION", "diffusion"), ("IONISATION RATE", "ion/att rate"),
    ("ATTACHMENT RATE", "ion/att rate"), ("MICRONS/CM", "diffusion"),
    ("EV.", "char. energy"), ("ELASTIC", "per-level"), ("EXC", "per-level"),
    ("ROTATION", "per-level"), ("VIBRATION", "per-level")]
  for (s, k) in kinds:
    if s in key: return k
  "other"

proc summary(files: seq[string]) =
  var byKind: OrderedTable[string, seq[float]]
  for f in files:
    for line in lines(f):
      let c = line.split('\t')
      if c.len < 7: continue
      let r = parseFloat(c[6])
      if r.classify in {fcNan, fcInf}: continue
      if parseFloat(c[5]) == 0: continue
      byKind.mgetOrPut(kindOf(c[1]), @[]).add r
  echo &"""{"kind":<16} {"n":>5} {"median":>7} {"p90":>7} {"max":>7}"""
  for k, v in byKind.mpairs:
    v.sort()
    let med = v[min(v.high, int(0.5 * v.len.float))]
    let p90 = v[min(v.high, int(0.9 * v.len.float))]
    echo &"{k:<16} {v.len:>5} {med:>7.2f} {p90:>7.2f} {v[^1]:>7.2f}"

let args = commandLineParams()
if args.len >= 2 and args[0] == "--summary": summary(args[1 .. ^1])
elif args.len >= 2: perDeck(args[0], args[1 .. ^1])
else: quit("usage: spread DECK OUT... | spread --summary TSV...")
