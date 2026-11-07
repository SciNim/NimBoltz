## Entry point of NimBoltz: runs the input decks on stdin, like Magboltz.
##
## NIMBOLTZ_SEED=<n> starts the random numbers from seed n instead of
## Magboltz's fixed 54217137 (for statistical comparisons, tests/stat).
import std/[os, strutils]
import rm48, magboltz

var s = initSession()
if existsEnv("NIMBOLTZ_SEED"):
  s.rng.rm.rm48in(parseInt(getEnv("NIMBOLTZ_SEED")), 0, 0)
s.magboltz()
