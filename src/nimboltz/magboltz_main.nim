## Entry point of NimBoltz: runs the input decks on stdin, like Magboltz.
##
## NIMBOLTZ_SEED=<n> starts the random numbers from seed n instead of
## Magboltz's fixed 54217137 (for statistical comparisons, tests/stat).
## NIMBOLTZ_THREADS=<n> runs the avalanches of the Townsend analyses on n
## threads, with random numbers of their own: the same estimates, not the
## same numbers as Magboltz.
import std/[os, strutils]
import rm48, magboltz

var s = initSession()
if existsEnv("NIMBOLTZ_SEED"):
  s.rng.rm.rm48in(parseInt(getEnv("NIMBOLTZ_SEED")), 0, 0)
if existsEnv("NIMBOLTZ_THREADS"):
  s.threads = max(parseInt(getEnv("NIMBOLTZ_THREADS")), 1)
s.magboltz()
