# NimBoltz design notes

How the Magboltz 11.21 port becomes idiomatic Nim. Agreed 2026-10-31.

## Starting point

`src/nimboltz/legacy/` holds the translation of `magboltz-11.21.f` (with
the fixes MB-01..MB-29 of the `fixes` branch of the Magboltz fork) that
`tools/f2nim` produced. It matches the gfortran build bit for bit. The
translator is retired: a better translator cannot produce idiomatic code.
From now on `legacy/` is edited by hand and shrinks as routine families
are rewritten. f2nim stays for one-off jobs (extracting DATA tables,
diffing a new Magboltz release).

## Order of work

1. ELIMIT family (done first: small, and it settles the patterns below).
2. SETUP / MIXER / OUTPUT: inputs, the tables, the results.
3. The MONTE family (MONTE, MONTEA/B/C, thermal versions, MONTEFT*, ...).
4. Gas cross sections (95 kLOC): data tables as data, shared interpolation
   helpers.

COMMON blocks are replaced family by family. While old and new code
coexist, a small bridge in `legacy/` (e.g. `legacy/elimits.nim`) reads the
COMMON blocks, calls the new code and writes the results back.

## Variants are compile-time parameters

Magboltz copies whole routines for each variant: field geometry (B = 0 or
B ∥ E, B ⊥ E, any other angle), thermal gas motion, time-of-flight
analyses. Each family becomes one generic implementation:

- `FieldGeometry` (`fields.nim`) is a `static` parameter; `when` picks the
  formulas, so each instantiation is as specialised as a Fortran copy.
- `withGeometry(g, G): body` turns the runtime choice into a constant once,
  at the top.
- Where the data differ, overloading does the dispatch: gas at 0 K uses one
  `CollisionTable` for the mixture, thermal motion one table per gas plus
  `ThermalGas` (`elimit.nim` has both overloads; they share the free-flight
  (`fields.nim`) and collision (`collisions.nim`) code).

## State: shared tables, copyable walkers

The Monte Carlo is meant to run on many threads, and later on a GPU.
State therefore falls into four groups:

| group | examples | mutable? | sharing |
|---|---|---|---|
| configuration | `Field`, energy grid, NMAX | no | copied (small) |
| tables | `CollisionTable` (30–170 MB) | no, after MIXER | one copy, shared |
| walker | `Electron`, RNG, `NullBounds`, `MaxwellSampler` | yes | one per thread/lane |
| accumulators | sums of MONTE | yes | per walker, merged at the end |

Rules:

- No globals in new code. Everything is passed in explicitly.
- Plain value objects, no `ref`, no `seq` in hot data, so that a walker is
  copied with `=` and a table can be uploaded as one block.
- Shared state is read-only. Magboltz mutates "tables" during the run
  (TCFMAX, the null-collision bounds, grows by 5% when exceeded); such data
  belongs to the walker (`NullBounds`).
- Random numbers come from a `var MagRng` that belongs to the walker. For
  now that is still the one global stream, so that results stay comparable
  with Fortran.

Known leftovers that matter here: table layout is Magboltz's (`cf[ie, i]`
with the energy bin fastest), which is poor for the per-bin scan over
processes. It will be transposed once MIXER is rewritten.

## Equivalence with Magboltz

Bitwise identity does not scale past this point. Rewriting changes the
evaluation order, which changes the last bits, and the Monte Carlo
amplifies any change into a different random history. The criteria:

- Non-MC code (cross sections, tables, derived constants): equal within a
  relative epsilon, with a dump comparison like `tests/gasdump`.
- MC results: equal within their statistical errors, i.e. Nim and Fortran
  runs with independent random streams agree within a few combined σ
  (`tests/stat`, to be written).
- Bitwise comparison stays as a cheap tripwire while a change keeps the
  evaluation order, as the ELIMIT rewrite does. `tests/regress` (MC decks)
  and `tests/gasdump` (every gas routine) remain the tools for that.

## Layout

```
src/nimboltz/
  legacy/        translated Magboltz, shrinking
  frt.nim        Fortran runtime for legacy (arrays, FORMAT, powi, ...)
  rm48.nim       Magboltz's random numbers, bit-exact
  vec3.nim       small 3-vector
  tables.nim     CollisionTable, NullBounds
  fields.nim     FieldGeometry, Field, Electron, free flight
  collisions.nim process choice, scattering
  thermal.nim    ThermalGas, Maxwell velocities
  elimit.nim     upper energy limit search (ELIMIT*)
  magboltz_main.nim  program: reads the deck on stdin
```
