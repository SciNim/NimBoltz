## Entry point of NimBoltz: runs the input decks on stdin, like Magboltz.
##
## NIMBOLTZ_SEED=<n> starts the random numbers from seed n instead of
## Magboltz's fixed 54217137 (for statistical comparisons, tests/stat).
## NIMBOLTZ_THREADS=<n> runs the avalanches of the Townsend analyses on n
## threads, with random numbers of their own: the same estimates, not the
## same numbers as Magboltz.
## NIMBOLTZ_WALKERS=cpu[:n] or gpu[:n] runs MONTE of the gas at
## 0 K with n independent electrons (walkdrift.nim), on NIMBOLTZ_THREADS
## CPU threads or a GPU (a build with -d:nimboltzGpu); again the same
## estimates within errors.
## NIMBOLTZ_TIMING=1 prints the times of the steps of each deck to stderr.
import std/[os, strutils]
import rm48, magboltz

var s = initSession()
if existsEnv("NIMBOLTZ_SEED"):
  s.rng.rm.rm48in(parseInt(getEnv("NIMBOLTZ_SEED")), 0, 0)
if existsEnv("NIMBOLTZ_THREADS"):
  s.threads = max(parseInt(getEnv("NIMBOLTZ_THREADS")), 1)
if existsEnv("NIMBOLTZ_WALKERS"):
  let w = getEnv("NIMBOLTZ_WALKERS").split(':')
  s.walkers.backend = case w[0]
    of "cpu": wbCpu
    of "gpu": wbGpu
    else: quit("NIMBOLTZ_WALKERS: cpu[:n] or gpu[:n], not " & getEnv("NIMBOLTZ_WALKERS"))
  if w.len > 1: s.walkers.n = parseInt(w[1])
  s.walkers.threads = s.threads
s.timing = existsEnv("NIMBOLTZ_TIMING")
s.magboltz()
