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

## Library

`nimboltz/api` runs NimBoltz from Nim code (since 2026-11-08):

```nim
import nimboltz/api

let r = simulate(Setup(mixture: @[(argon, 90.0), (carbonDioxide, 10.0)],
                       efield: 500.0.V•cm⁻¹))
echo r.velocity.z, " ", r.spread(r.transverse)
```

- A `Setup` is one deck, typed: gases as the enum `Gas` (Magboltz's gas
  numbers), temperature in K, pressure in Torr, B in T, the angle in
  degrees, and `collisions` (real collisions of the drift run) in place of
  NMAX · nscale. Defaults: 20 °C, 760 Torr, no B, thermal motion, no
  Penning transfers, 2 · nscale collisions, Magboltz's seed, energy range
  found by ELIMIT.
- `simulate` runs it as the program runs a deck (`Session.run`) and
  returns `Results`: values with errors (measuremancer `Measurement`, from
  the errors in per cent of Magboltz) in Magboltz's printed units, the
  analyses' results if they ran (`Option`), the collision frequencies by
  process, the energy spectrum and the warnings. `characteristicEnergy` and
  `spread` give D/µ and √(2D/v) of a diffusion coefficient.
- Nothing is printed unless `output` is set: every `fwrite` goes to the
  run's `output` (nil: nowhere), so with `output = stdout` the printout is
  Magboltz's, bit for bit. Where Magboltz stops the program (input error,
  energy range too small, avalanche limits) the run raises a
  `MagboltzError`; the program catches it and stops. `tests/unit/tapi`
  checks that a setup prints what its deck prints and that the results
  are the printed numbers.
- Runs share no state (each `simulate` has a session of its own), but
  measuremancer numbers its measurements with a global counter, so runs
  on several threads would need that made thread-safe first.

## GPU

Since 2026-11-08 the collision loop of MONTE runs on NVIDIA GPUs:
`walkers.nim` has independent electrons (`Walker[ColdGas]` at 0 K,
`Walker[WarmGas]` with thermal motion of the gas), each with its own
`Ranmar` (RM48 one number at a time, the numbers of `rm48`), and
`gpuwalk.nim` runs one per GPU thread, in all four field geometries. The
kernel is NimBoltz's own code (`nextHit`, `findProcess`, `energyLoss`,
`scatter`, `advance`, ...), compiled to CUDA at Nim compile time by
Crucible, the Nim-to-GPU compiler of tattletale
(github.com/mratsim/tattletale; local fixes on branch `nimboltz` of
~/src/gpu/tattletale), and to PTX by NVRTC at run time. Unchained units
compile to plain `double`s. One kernel per geometry and target (8),
written once as a template.

What this asks of the code the GPU runs:

- no `break`, `continue` or `case` (Crucible's IR has none): loops with a
  flag, `if` chains
- `int32` in what is copied to the GPU (Crucible's `int` is 32 bits), so
  `CollisionTable`, `ThermalGas` and `MaxwellSampler` have `int32` integers
- no references, no `File`, no strings; the random number generator is a
  generic parameter (`MagRng` on the CPU, `Ranmar` on a GPU)
- collision counts per process and energy bin: atomic adds into one tally
  (`GpuTally`)

### The drift run on walkers

`NIMBOLTZ_WALKERS=cpu[:n]`, `gpu[:n]` or `gpu32[:n]` (`Setup.walkers`,
`Session.walkers`) runs MONTE on walkers (`walkdrift.nim`, estimators in
`ensemble.nim`), and ELIMIT too when the deck leaves EFINAL to it; the
Townsend analyses keep Magboltz's one trajectory. The CPU backend runs the
walkers on threads, the GPU one needs a build with `-d:nimboltzGpu` and
`--path:` to tattletale; `gpu32` runs it in single precision (below).
The results fill the same `Transport`, printout and `Results`; one line
more in the printout (walkers, collisions each),
and the collisions counted (`Results.collisions`) may exceed NMAX on a GPU.
`NIMBOLTZ_TIMING=1` (`Session.timing`) prints the time of each step of a
deck (MIXER and ELIMIT, MONTE, analyses) to stderr.

- A walker starts as Magboltz's electron, and runs in chunks of doubling
  length until the mean energy and drift velocity of two chunks agree
  (warm-up), then its sums restart.
- ELIMIT (`magboltz.walkerEnergyRange`, `walkdrift.energyLimitExceeded`):
  EFINAL from 0.5 eV (8 eV above 15 V/cm/Torr) up by √2 until no
  collision reaches the last energy bin. Each try starts the walkers from
  Magboltz's electron (in a range too small their energies run away; the
  walkers keep their random numbers: seeding 16384 generators takes 0.6 s),
  warms them up, then counts ELIMIT's collisions (NMAX/10, NMAX/20 with B)
  among all of them. Not the warm-up's: Magboltz counts its one
  electron's start, many starts would raise EFINAL for nothing. A range
  far too small shows in a probe of 250 collisions each, or when the
  warm-up gets more than 1000 per sample in the last bin (MONTE's limit;
  at 10 per sample, bursts of one walker running away, caught in a
  warm-up of many more collisions than the sample, raised EFINAL in 22 of
  95 decks). The
  walkers of the last try go on to MONTE warm. Where MONTE then has more
  than 1000 collisions per NMAX/10 in the last bin (Magboltz stops there),
  walkers raise EFINAL and run MONTE again. Over the regression decks the
  choice is Magboltz's in 75 of 94, the others one √2 step off, 11 up and
  8 down, as Magboltz's own varies between seeds; a tail that builds up slowly (Ar/CBrF₃ 90/10:
  electrons escaping into argon heat up over 10⁵ collisions) is seen by
  Magboltz's long trajectory, not by short walkers (gas51: 4 eV against 8).
  Full scale (4·10⁷), 16 threads: arco2_nob_t1 7.0 → 3.3 s a deck, CF₄
  (gas01) 23.1 → 5.5 s.
- Drift velocity Σx/Σt, mean energy Σ∫e dt/Σt (B ⊥ E and oblique: Σ e/ν,
  as MONTEB/C). Diffusion from the spread of the walkers around the mean
  drift: with positions marked halfway, D = [Σ(x − wt)² − Σ(x₁ − wt₁)²] /
  2Σ(t − t₁), so the constant from the correlation of the first
  collisions cancels exactly. Errors from the scatter of 10 blocks of
  walkers. Townsend coefficients from the walkers' ionisations and
  attachments.
- How many walkers: the main run of each walker must last at least as long
  as its warm-up (a few energy relaxation times), else walkers of
  different energies spread ballistically and D comes out too large
  (Ar/He 90/10 at 1000 V/cm, 2000 collisions per walker: DL 10× too
  large). On CPUs, where time goes with all collisions, a warm-up chunk
  spends at most half the collisions of the main run, dropping half of the
  walkers when it would spend more (down to 100), and after the warm-up
  walkers are dropped until each runs at least its warm-up length:
  molecular gases keep many short walkers, noble gases get few long ones.
  On a GPU, where time goes with the collisions per walker (below some
  10⁴ walkers the GPU is not full), all walkers stay and each runs
  max(NMAX/n, warm-up) collisions: maybe more than NMAX in all, for the
  same time. The warm-up ends at NMAX/200 collisions per walker.

Equivalence: the same walkers on CPU and GPU take the same collisions
for a while (`tests/gpu/tgeometry.nim`, `exact = true`: no fused
multiply-adds), until a last-bit difference between CUDA's `logf`, `acos`,
`sin` and glibc's sends a walker on another path: of 256 walkers, 115–256
still have the same null collisions after 2000 collisions. Many walkers
estimate what one trajectory does, not the same numbers: Ar/CO₂ 90/10 at
1000 V/cm, B = 1 T, all geometries, 0 K and thermal, 4·10⁷ collisions:
|z| < 2 for drift, energy and diffusion against `simulate`; the
regression decks with walkers on CPU and GPU agree within errors
(`statcmp`). The walkers' diffusion errors are several times smaller
than Magboltz's for the same collisions. Where they disagree, Magboltz was
off: Ne/CO₂ 90/10 at 20 kV/cm (penning_neco2_t1), DL at the fast scale
(4·10⁶) 1421 ± 21 %, at 4·10⁷ 3940–4726 (3 seeds), walkers 4461 and 4644 at
4·10⁷, 4879 at 4·10⁶ — Magboltz's pairs NCORLN apart need long runs.
Ar/He 90/10 at 1000 V/cm, 4·10⁸ collisions: Magboltz DL 413, 458, 485,
walkers 471 ± 19.

Walkers have free flights of their own (`walkers.nextHit`): the null-
collision method of `freeflight.nextHit` with what depends on the start of
a flight computed once, multiplications for divisions, the logarithm of the
flight time and the Maxwell numbers of thermal motion in single precision
(FP64 runs at 1/64 of FP32 on gaming GPUs), and v(t) = v₀ + a t for thermal
motion at B = 0. `freeflight.nextHit` stays Magboltz's, bit for bit.

Null collisions: Magboltz's bounds in 8 bands of starting energy
(`nullBounds`) cover every energy an electron reaches from the band within
TMAX, nearly the largest frequency of all; with thermal motion it has one
bound per gas over all energies. Walkers at B = 0 and B ∥ E use exact
bounds per time window instead (`walkers.windowHit`): the energy over the
next 2 mean free times is known in closed form, the largest collision
frequency over it (a sparse range-maximum table, `FlightMaxima`, two
lookups) bounds it strictly; no trial in the window, the flight goes on
in the next one. Where the bound is well above the frequency at the
window's start (a Ramsauer minimum), the window shrinks to 2 / bound and
gets its own bound. With thermal motion per gas, in the molecule's frame,
with a margin of 6 thermal speeds. Null collisions per real one: Ar/CO₂
90/10 at 1000 V/cm 5.6 → 0.08 (0 K), 9.6 → 0.10 (thermal); Ar/CF₄ 90/10
thermal (gas01) 4.3 → 0.12. (Window lengths of 1, 2, 4 mean free times
cost about the same on CPUs; on the GPU 2 is best.) (Tried first:
each band's own maximum, raised when exceeded: biased, mean energies 2–6 %
low over the regression decks, mean z −0.83 — flights crossing energies
above the bound without a trial there lose collisions.)

Speed (RTX 3090; `tests/gpu/tbench.nim`, 16384 walkers, B = 0, energy range
from ELIMIT), M collisions/s: Magboltz's free flight → the walkers' own →
with the exact window bounds → shorter windows → single precision:

| | 0 K | thermal |
|---|---|---|
| CO₂ | 284 → 330 → 328 → 298 → 506 | 140 → 145 → 158 → 163 → 319 |
| Ar/CO₂ 90/10 | 57 → 118 → 245 → 268 → 487 | 14 → 42 → 116 → 134 → 269 |
| Ar/CF₄ 90/10 | → 142 (4 mean free times) → 200 → 481 | → 63 → 95 → 255 |

against ~3.6 M/s for one trajectory on a CPU core. Tried, no gain: walkers
in local memory (−5 %), the cumulative frequencies transposed so that a
bin's processes are side by side (−8 %: in Magboltz's layout the threads of
a warp at nearby energies read neighbouring words). Whole decks keep MIXER
on the CPU.

### Single precision

`gpu32` runs the same kernel in single precision: Crucible's
`cudaSingle` (and `webgpuSingle`; WGSL has no doubles) demotes every
`float64` of the device code to `float32` (pass `pass_demote_f64`, local
crucible branch). The host hands over walkers, tables and field in single
precision too (`single.nim`: `pack32`/`unpack32` turn an object's floats
into 4-byte ones; with 4-byte fields only an object has no padding, on the
GPU as on the host). The walkers' random numbers become the original
24-bit RANMAR (`rm48.toSingle`: state on a grid of 2⁻²⁴, Marsaglia and
Zaman's constants), whose arithmetic is exact in single precision; the
GPU's tallies are 32-bit counters, added up on the host after each launch
(single-precision counts stop at 2²⁴). 1.7–2.7× faster than double
precision on an RTX 3090 (a memory-bound kernel: the 64× of the
arithmetic does not show).

## WebAssembly

Since 2026-11-09 NimBoltz also builds for WebAssembly with Emscripten
(`web/build.sh`), as the program for Node (`web/build.sh cli DIR`, decks on
stdin) and as a module for web pages (`src/nimboltz/wasm.nim`: a `Setup`
as JSON in, the `Results` as JSON out). `web/` is a small page that runs
it in a Web Worker: a mixture, the conditions, one or several fields, the
results with errors, a table and a plot of all runs, the energy
distribution, the collision frequencies and Magboltz's printout.

- Two modules: `nimboltz` on one thread (`--threads:off`), and
  `nimboltz-mt` with threads (Emscripten pthreads: Web Workers sharing the
  memory, a SharedArrayBuffer, which browsers allow only with the
  server's cross-origin isolation headers; `web/serve.py`). `worker.js`
  picks by `crossOriginIsolated`. The threads start with the module, one
  per core (`PTHREAD_POOL_SIZE`): a run blocks its thread in
  `joinThreads`, so a thread started then could never load; the worker
  caps the setup's `threads` at that. The threads' workers load
  `worker.js` (Emscripten 7 starts them from the script that started the
  module), which then only loads the module. Crucible's CPU thread pool
  needs Linux futexes, so it stays out under Emscripten (`cpuPool`).
  CPU walkers on 32 threads: Ar/CO₂ 8M collisions 1.8 s against 16.3 s.
- `int` has 32 bits on wasm32. Counters that can pass 2³¹ (null
  collisions: billions in long avalanche runs) are `int64`, and floats
  are clamped before they become integers.
- Memory: ~300 MB for a run (tables of 4000 energy bins), growing as
  needed up to 4 GB; repeated runs in one module don't grow it.
- Equivalence: wasm's maths library (musl) differs from glibc's in the
  last bit of some functions, so trajectories diverge after a while. Of
  the 96 regression decks 52 stay byte-identical, the rest agree
  statistically (`tests/wasm/regress.sh`: max |z| 2.6); the cross
  sections of all gases agree within 7e-15 relative.
- Speed: 1.4–1.7 times the native time under Node.

### WebGPU

With `TATTLETALE=<tattletale root> web/build.sh` the page also offers the
drift run on GPU walkers through WebGPU (2026-11-09):
`src/nimboltz/webgpuwalk.nim` has the interface of `gpuwalk.nim`, its
kernels are `walkWith` compiled to WGSL by Crucible (`webgpuSingle`, 8
kernels: 4 geometries, 0 K or thermal), and `web/webgpu.js` (an
Emscripten JS library) holds the device's buffers and pipelines. The worker
opens the device; setup `"walkers": "gpu"`.

- Single precision only (WGSL has no f64): the walkers, tables and field
  packed by `single.nim`, as `gpu32` on CUDA; all fields are 4 bytes, so
  WGSL's struct layout is the packed one.
- WGSL needed more of Crucible (local branch `nimboltz`, commit 5762f9af):
  storage pointers as function parameters (`-d:wgslPointerParams`,
  `requires unrestricted_pointer_parameters, pointer_composite_access`):
  the tables (46 MB) pass by reference instead of being copied into
  function memory; atomics (`array<atomic<u32>>`), struct ternaries,
  aliased pointer arguments, unreachable functions, `_`.
- Precise maths (`-d:wgslPreciseMath`): WGSL's `acos`/`asin` inherit
  `atan2`'s 4096 ULP and `sin`/`cos` are good to 2⁻¹¹ (measured 2e-4 on
  SwiftShader). With them the walkers drifted 0.8 % fast with thermal
  motion (49.84 against 49.42 μm/ns for Ar/CO₂ 90/10 at 1000 V/cm); with
  Crucible's precise ones (to 4e-7) a walker's state after its first
  collision agrees with the CPU's to 1e-6 in all 8 kernels
  (`tests/wasm/webgpu.sh`), and the drift run with Magboltz's (thermal:
  49.384 ± 0.017 against 49.445 ± 0.048; B ⊥ E at 0 K: v_z, v_y and the
  mean energy within 1σ).
- Bindings: 5 storage buffers (walkers, table, range maxima, the scalars
  and field in one struct, the tally: processes then energy bins), as
  WebGPU allows 8 per stage. Reading back (walkers, tally) is
  asynchronous, so the module is built with Asyncify
  (`ASYNCIFY_IGNORE_INDIRECT`: 7.1 MB instead of 6.5, the CPU path 2–7 %
  slower; without it 11.4 MB and 32 % slower).
- One dispatch walks at most 1000 collisions per walker (browsers stop
  long GPU work); errors of the device (`uncapturederror`, a failed
  pipeline) stop the run with an error instead of wrong numbers.
- Tested on SwiftShader only (Chromium's CPU Vulkan, 2–3 M collisions/s);
  speed on a real GPU untested.

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
  gases.nim, gases/  cross sections of the 80 gases (GAS1..GAS80), the
                 enum Gas
  deck.nim       the input deck (READIN)
  conditions.nim constants and what follows from a deck (SETUP, SETUPT)
  simulation.nim state of one run: Run[M], ColdMixture, WarmMixture,
                 collision counts and frequencies (COLF)
  analyses.nim   Townsend analyses after MONTE (ALPCALC*, avalanche runs)
  report.nim     printout (PRNTER*, OUTPUT*, OUTPUT1, OUTPUT2*)
  magboltz.nim   the program: Session, per-deck driver (MIXER, ELIMIT,
                 MONTE calls)
  walkers.nim    independent electrons, any geometry, 0 K or thermal (CPU threads)
  ensemble.nim   transport coefficients of an ensemble of walkers
  walkdrift.nim  MONTE on walkers: warm-up, how many, main run
  gpuwalk.nim    walkers on an NVIDIA GPU (Crucible, NVRTC; -d:nimboltzGpu)
  webgpuwalk.nim walkers on WebGPU in a page (Crucible, WGSL; -d:nimboltzWebGpu)
  magboltz_main.nim  entry point: decks on stdin
  results.nim    the results of a run as values with errors (Results)
  api.nim        the library: Setup, simulate
  wasm.nim       the library for web pages (Emscripten): JSON in and out
web/             build script and page of the WebAssembly build (webgpu.js:
                 the WebGPU side of webgpuwalk.nim)
```
