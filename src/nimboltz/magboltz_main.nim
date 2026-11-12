## Entry point of NimBoltz: runs the input decks on stdin, like Magboltz.
##
## By default the drift run (MONTE, and ELIMIT) follows many independent
## electrons ("walkers", walkdrift.nim) on all cores, and the Townsend
## analyses follow their avalanches with exact bounds of the collision
## frequency: the same results as Magboltz within errors, and faster.
## NIMBOLTZ_FORTRAN=1 runs Magboltz 11.21 bit for bit instead: one electron,
## its null-collision method and formulas, one thread; with the same deck
## the printout is Magboltz's, number for number (for cross-checks).
##
## NIMBOLTZ_SEED=<n> starts the random numbers from seed n instead of
## Magboltz's fixed 54217137 (for statistical comparisons, tests/stat).
## NIMBOLTZ_THREADS=<n>: n threads instead of all cores (walkers and
## avalanches; with NIMBOLTZ_FORTRAN, n > 1 splits the avalanches of the
## analyses over n threads, with random numbers of their own).
## NIMBOLTZ_WALKERS=cpu[:n], gpu[:n] or gpu32[:n]: n walkers (default 100
## per thread), on the CPU or a GPU (a build with -d:nimboltzGpu; gpu32: in
## single precision).
## NIMBOLTZ_RNG=xoshiro (default) or ranmar: the random numbers of the
## walkers and the avalanches (xoshiro256+, or Magboltz's RM48 one number at
## a time; GPUs always use RM48).
## NIMBOLTZ_TIMING=1 prints the times of the steps of each deck to stderr.
import std/[os, strutils]
import rm48, magboltz

var s = initSession(fortran = getEnv("NIMBOLTZ_FORTRAN") notin ["", "0"])
if existsEnv("NIMBOLTZ_SEED"):
  s.rng.rm.rm48in(parseInt(getEnv("NIMBOLTZ_SEED")), 0, 0)
if existsEnv("NIMBOLTZ_THREADS"):
  s.threads = max(parseInt(getEnv("NIMBOLTZ_THREADS")), 1)
if existsEnv("NIMBOLTZ_WALKERS"):
  let w = getEnv("NIMBOLTZ_WALKERS").split(':')
  s.walkers.backend = case w[0]
    of "cpu": wbCpu
    of "gpu", "gpu32": wbGpu
    else: quit("NIMBOLTZ_WALKERS: cpu[:n], gpu[:n] or gpu32[:n], not " &
               getEnv("NIMBOLTZ_WALKERS"))
  s.walkers.single = w[0] == "gpu32"
  if w.len > 1: s.walkers.n = parseInt(w[1])
if existsEnv("NIMBOLTZ_RNG"):
  s.walkers.rng = case getEnv("NIMBOLTZ_RNG")
    of "ranmar": rkRanmar
    of "xoshiro": rkXoshiro
    else: quit("NIMBOLTZ_RNG: ranmar or xoshiro, not " & getEnv("NIMBOLTZ_RNG"))
s.timing = existsEnv("NIMBOLTZ_TIMING")
s.magboltz()
