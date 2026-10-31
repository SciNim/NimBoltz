# NimBoltz design notes

How the Magboltz 11.21 port becomes idiomatic Nim. Agreed 2026-10-31.

## Starting point

`src/nimboltz/legacy/` holds the translation of `magboltz-11.21.f` (with
the fixes MB-01..MB-29 of the `fixes` branch of the Magboltz fork) that
`tools/f2nim` produced. It matches the gfortran build bit for bit. The
translator is retired: a better translator cannot produce idiomatic code.
From now on `legacy/` is edited by hand and shrinks as routine families
are rewritten. f2nim stays for one-off jobs (extracting DATA tables,
diffing a new Magboltz release). Defects found later (MB-30 onwards) are
fixed by hand in both `legacy/` and the fork, so that the bitwise
comparison keeps holding.

## Order of work

1. ELIMIT family (done first: small, and it settles the patterns below).
2. SETUP / MIXER / OUTPUT: inputs, the tables, the results. MIXER and
   MIXERT done (`mixture.nim`, one builder for both table layouts).
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
processes in theory. Tried after the MIXER rewrite (2026-10-31): storing
each bin's processes contiguously made the legacy MONTE 3–8% slower, so
the layout is revisited with the MONTE rewrite.

## Equivalence with Magboltz

Bitwise identity does not scale past this point. Rewriting changes the
evaluation order, which changes the last bits, and the Monte Carlo
amplifies any change into a different random history. The criteria:

- Non-MC code (cross sections, tables, derived constants): equal within a
  relative epsilon, with a dump comparison like `tests/gasdump`.
- MC results: equal within their statistical errors, i.e. Nim and Fortran
  runs with independent random streams agree within a few combined σ
  (`tests/stat`: the Nim build runs with `NIMBOLTZ_SEED`, `statcmp`
  compares every result Magboltz prints with an error, |z| ≤ 5 with a 1%
  error floor added in quadrature).

Calibration of `tests/stat` (2026-10-31, port unchanged, so only the
random numbers differ; fast mode, seed 1): of 96 decks, 90 pass (mean
Σz²/n = 0.67, max |z| ≤ 4.7), 5 are not compared (one of the runs took
more than an hour), 1 fails (penning_neco2_t0). What it showed about
Magboltz's printed errors:

- The per-process collision frequencies carry counting errors (±0.05% for
  4·10⁶ collisions) but differ by several percent between random streams,
  because collisions along one electron's history are correlated. They
  are not compared.
- Attachment rates of strongly attaching gases have the same problem
  (SF₆: ±0.08% printed, 1.3% difference); the 1% floor covers it.
- The SST fits can print errors far below their real scatter: gas15 SST
  drift 15.9 ± 0.13% against 13.9 ± 1.6% from the MC of the same run;
  penning_neco2_t0 SST longitudinal diffusion 4283 ± 0.3% vs 5406 ± 0.8%
  (TOF values agree). These are the expected failures with this criterion.
- In fast mode the MC diffusion estimates have heavy tails (gas04,
  longitudinal: 46 ± 24% with one seed, 440 ± 24% with another).
- EFINAL found by ELIMIT depends on the random numbers; with NMAX/10
  collisions it is often too low for the main run, which then stops with
  "ENERGY OUT OF RANGE". `tests/stat` therefore runs both builds on copies
  of the decks with EFINAL fixed at twice the value Fortran chose
  (`fixefinal.nim`); ELIMIT itself is covered by `tests/regress`.

An empirical error per result, from several seeds of the unchanged port,
would remove the guesswork from the floor; not done yet.
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
  mixture.nim    GasXsec (output of a gas routine), collision tables of a
                 mixture (MIXER, MIXERT)
  magboltz_main.nim  program: reads the deck on stdin
```
