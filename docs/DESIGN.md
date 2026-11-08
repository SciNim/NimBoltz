# NimBoltz design notes

How the Magboltz 11.21 port becomes idiomatic Nim. Agreed 2026-10-31.

## Starting point

`tools/f2nim` translated `magboltz-11.21.f` (with the fixes MB-01..MB-29
of the `fixes` branch of the Magboltz fork) into `src/nimboltz/legacy/`,
bit for bit equal to the gfortran build. The translator is retired: a
better translator cannot produce idiomatic code. `legacy/` was then
rewritten by hand, routine family by routine family, and removed
(2026-11-07). f2nim stays for one-off jobs (extracting DATA tables,
diffing a new Magboltz release). Defects found later (MB-30 onwards) are
fixed by hand in both the port and the fork, so that the bitwise
comparison keeps holding.

## Order of work

1. ELIMIT family (done first: small, and it settles the patterns below).
2. SETUP / MIXER / OUTPUT: inputs, the tables, the results. Done: MIXER
   and MIXERT (`mixture.nim`, one builder for both table layouts), READIN
   (`deck.nim`), SETUP* (`conditions.nim`), PRNTER*, OUTPUT* (`report.nim`)
   and the main program (`magboltz.nim`).
3. The MONTE family. Done: the main Monte Carlo (MONTE, MONTEA/B/C and
   the thermal versions; `drift.nim`) and the avalanche runs behind the
   SST and time-of-flight analyses (MONTEFT, MONTEFTA/G/H, MONTEFD and
   the thermal versions; `cascade.nim`), one generic driver each, and
   their analyses (the eight ALPCALC drivers, SST, PT, TOF, FRIEDLAND;
   `townsend.nim`, driven by `analyses.nim`).
4. Gas cross sections (95 kLOC): first every GASn by hand, bit for bit
   (2026-11-06), then one generic engine (2026-11-08): `gases/model.nim`
   defines the description of a gas, `GasModel`, whose fields are the
   variants of `docs/gas-features.md`; each gas module is its tables plus
   one `const GasModel`, and `engine.fill`, instantiated per model,
   computes the cross sections. All 61 real gases use it; the three model
   gases (26–28) and the empty slots (62–80) stay small hand-written
   routines in `gases/models.nim`. `tests/gasengine` compares the engine
   with the hand-written modules (commit 9e63630) on 96 cases per gas:
   largest relative difference 2.5e-13 (most gases below 1e-15); 15–40% of
   the values differ in the last bits, from a different order of
   operations, so MC decks stay byte-identical only by luck.

COMMON blocks were replaced family by family, with small bridges in
`legacy/` while old and new code coexisted. Now the state of a run is a
`Run[M]` value (`simulation.nim`; M = `ColdMixture` or `WarmMixture`
selects the thermal variants at compile time) and what lives across decks
(random numbers, position rings of MONTE, the table memory) a `Session`
(`magboltz.nim`). No module keeps global state.

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
- Random numbers come from a `var MagRng` that belongs to the walker. By
  default that is the one stream of the `Session`, passed to every run in
  Magboltz's order, so that results stay comparable with Fortran bit for
  bit.
- Threads (2026-11-08, opt-in with `NIMBOLTZ_THREADS=n`): the avalanche
  runs of the Townsend analyses (`cascade.nim`) are split over n walkers,
  each with its share of NMAX and of the limit on primaries, its own
  target state and random numbers seeded from the session's stream
  (`analyses.avalanche`). Primaries are independent apart from the random
  numbers, so their counts, histograms and plane sums just add up
  (`cascade.merge`); the estimates are the same, the numbers are not, and
  the comparison with Fortran is statistical (`tests/stat`). The main
  drift run (MONTE) is still one trajectory on one thread.

Known leftovers that matter here: table layout is Magboltz's (`cf[ie, i]`
with the energy bin fastest), which is poor for the per-bin scan over
processes in theory. Tried after the MIXER rewrite (2026-10-31): storing
each bin's processes contiguously made the legacy MONTE 3–8% slower, so
the layout is revisited with the MONTE rewrite.

`drift.nim` keeps Magboltz's statistics as they are: one trajectory, 10
segments, errors from the scatter of the segment estimates, diffusion from
pairs of collisions NCORLN apart, kept in a ring of NCOLM positions
(64 MB). Each geometry has its own estimators (B ⊥ E splits the y-z spread
along and across the drift; B at other angles rotates the tensor to the
frame of E). The `Session` (`magboltz.nim`) keeps one ring per Magboltz
routine across runs, because the Fortran arrays are static and MONTEC
reads unfilled entries when NCOLM exceeds two segments (reduced NSCALE
only). Merging walkers (threads) needs per-walker rings and a different
error estimate; not done.

`cascade.nim` (MONTEF*) follows avalanches: primaries one after another,
secondaries, Auger and Penning electrons on a stack (8000, as Magboltz),
sampled where they cross planes in time (`pkTime`) or in space along E
(`pkSpace`, B = 0 only), a static parameter like the field geometry. The
ten Fortran routines differ in details beyond their variant (the energy
loss clamp after ionisation, the order of a few products, MONTEFDT's
renormalised direction); these are reproduced under `when` with a comment
naming the routine, so bitwise comparison holds. Thermal runs compute the
ionisation secondary in the rest frame of the molecule and push it
without transforming it back to the lab frame, as Magboltz does (an
error of order v_gas/v_e). Speed equal to the legacy routines (±4%).

Speed (2026-10-31): `drift.nim` is 2–4% slower than the legacy routines.
Nim emits a generic instance into the C file of the module that first
uses it, so calls between modules cannot be inlined; the small hot procs
(`fly`, `velocityAt`, `findProcess`, ...) are therefore `{.inline.}`. With
LTO (`--passC:-flto=auto --passL:-flto=auto`) the gap is 0–2%, at the
cost of much slower builds (95 kLOC of gas routines). About half of the
run time was the random number generator (RM48) and the Maxwell sampler;
RM48 without branches (2026-11-08: the wrap-arounds of each number are
selects, the same numbers) made it 4.5 times faster and whole runs about
25% faster.

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
- Time-of-flight diffusion of strongly attaching gases rests on few
  electrons at the later planes, and its printed error does not reflect
  that: gas45 longitudinal, 483 ± 1.9% (Fortran), 347 ± 0.8%, 404 ± 44.7%,
  411 ± 26.9% (Nim, other random numbers, 1 or 4 threads).
- In fast mode the MC diffusion estimates have heavy tails (gas04,
  longitudinal: 46 ± 24% with one seed, 440 ± 24% with another).
- EFINAL found by ELIMIT depends on the random numbers; with NMAX/10
  collisions it is often too low for the main run, which then stops with
  "ENERGY OUT OF RANGE". `tests/stat` therefore runs both builds on copies
  of the decks with EFINAL fixed at twice the value Fortran chose
  (`fixefinal.nim`); ELIMIT itself is covered by `tests/regress`.

Measured errors (2026-10-31): `tests/stat/spread.sh` runs each deck with
8 more seeds of the unchanged port and records, per printed result, the
spread of the 9 runs (with the Fortran run) in `spread-fast.tsv`; `run.sh`
uses the larger of that and the printed error, with a floor of 0.25%. Ratio
of measured to printed error over 91 decks (5 too slow): drift velocity and mean energy
median 0.93 (worst 1.7); longitudinal diffusion median 1.3, 10% of decks
above 2.4, worst 4.8 (low-field noble gases, fast mode: ±25% printed,
100–125% real); α/attachment median 0.8, worst 16 (CCl₄ attachment ±0.09%
printed, 1.4% real; SF₆ ±0.08% / 0.6%). With these errors seed 1 of the
unchanged port passes all 91 decks (mean Σz²/n = 0.43, max |z| = 3.2),
including the SST cases that failed with the printed errors.
- Bitwise comparison stays as a cheap tripwire while a change keeps the
  evaluation order, as the ELIMIT rewrite does. `tests/regress` (MC decks)
  and `tests/gasdump` (every gas routine) remain the tools for that.

## Units

Every physical quantity has its unit in its type (Unchained, since
2026-11-08): `eV`, `ps`, `m`, `m•ps⁻¹`, `cm²`, `V•cm⁻¹`, `kGauss`, ...
`src/nimboltz/units.nim` defines all compound units once (a compound unit
type defined in two modules would be two different types) and the two
helpers that state Magboltz's conventions:

- `unity(U)`: the number 1 in a dimensionless unit, e.g. `unity(m•cm⁻¹)`
  is 0.01. Magboltz writes its conversions as literals (`*0.01D0`,
  `*1D-9`); the code states them as units, and Unchained computes the
  factor. `conditions.nim` asserts at compile time that each constant of
  SETUP equals Magboltz's literal arithmetic (the 1e-19 in CONST1, for
  instance, is T/kG · m²/cm² · ps²/s² = 10 · 1e4 · 1e-24).
- `perElectron` = 1 eV/V: Magboltz uses a potential in V as an electron's
  energy in eV. Unchained knows the eV as 1.602e-19 J, so the step from
  volts to an electron's energy is written out (multiplying by 1.0 changes
  no bit). kT in J goes through Magboltz's own e (`echarg / perElectron`),
  not Unchained's.

Speeds come from energies as v = vfac·√E, so energies appear with half
powers (`eV¹ᐟ²`, `cm⁻²•ps⁻¹•eV⁻¹ᐟ²` for VAN); this needed rational powers
in Unchained (branch `rational-powers` of ~/src/units, with exact
power-of-ten conversions and `-d:unchainedImplicitScale=warn|error`).

Rules that keep Magboltz's arithmetic:

- The trajectory (fields, flights, tables, null-collision bounds, energy
  loss) keeps every operation of the Fortran: all 96 decks stay
  byte-identical. A conversion at the end of an expression is `.to(U)`
  (Unchained multiplies last, as Magboltz does); Unchained converts on its
  own when a product mixes prefixes of one unit (cm and m), and then
  multiplies right after that product.
- Results that are only printed (diffusion in cm²/s, velocities in cm/s,
  Townsend coefficients) are written in the natural form, D = Σ/2t
  converted at the end, where Magboltz folds its literal into the middle
  (`5D15*SUM/T`). That changes the last bit of a printed 4-5 digit
  number, which Magboltz's output does not show.
- Build with `-d:unchainedImplicitScale=warn` to list every conversion
  Unchained inserts; on 2026-11-08 these are the printed results above
  and Magboltz's own `DZ*EFIELD*100`.
- Unchained's `==` is approximate (`almostEqual`); exact comparisons use
  `.float`. Logarithms of quantities that only enter through their
  differences take `inUnits` (the value in its unit).

## Layout

```
src/nimboltz/
  frt.nim        Fortran FORMAT output, list-directed input, fstop
  farrays.nim    Fortran-style 2D arrays (collision tables, gases)
  rm48.nim       Magboltz's random numbers, bit-exact
  vec3.nim       small 3-vector (Vec3 directions, Vector[T] with units)
  units.nim      the units (Unchained) and Magboltz's conversions
  tables.nim     CollisionTable, NullBounds
  fields.nim     FieldGeometry, Field, Electron, free flight
  collisions.nim process choice, scattering
  thermal.nim    ThermalGas, Maxwell velocities
  freeflight.nim next real collision (null-collision method), gas at 0 K
                 or with thermal motion
  elimit.nim     upper energy limit search (ELIMIT*)
  drift.nim      main Monte Carlo: drift velocity, diffusion, mean energy,
                 Townsend coefficients (MONTE, MONTEA/B/C, MONTE*T)
  cascade.nim    avalanches sampled at time or space planes (MONTEFT*,
                 MONTEFD, MONTEF*T)
  townsend.nim   analyses of the avalanches: sampling planes, steady-state
                 and pulsed Townsend, time of flight, Friedland estimate
  mixture.nim    GasXsec (output of a gas routine), collision tables of a
                 mixture (MIXER, MIXERT)
  gases.nim, gases/  cross sections of the 80 gases (GAS1..GAS80)
  deck.nim       the input deck (READIN)
  conditions.nim constants and what follows from a deck (SETUP, SETUPT)
  simulation.nim state of one run: Run[M], ColdMixture, WarmMixture,
                 collision counts and frequencies (COLF)
  analyses.nim   Townsend analyses after MONTE (ALPCALC*, avalanche runs)
  report.nim     printout (PRNTER*, OUTPUT*, OUTPUT1, OUTPUT2*)
  magboltz.nim   the program: Session, per-deck driver (MIXER, ELIMIT,
                 MONTE calls)
  magboltz_main.nim  entry point: decks on stdin
```
