## Entry point of the translated Magboltz: reads the input deck on stdin.
##
## NIMBOLTZ_SEED=<n> starts the random numbers from seed n instead of
## Magboltz's fixed 54217137 (for statistical comparisons, tests/stat).
import std/[os, strutils]
import frt, legacy/magboltz

if existsEnv("NIMBOLTZ_SEED"):
  rng.rm.rm48in(parseInt(getEnv("NIMBOLTZ_SEED")), 0, 0)
magboltz()
