## NimBoltz as a library: a gas mixture and the fields in, the transport
## coefficients out, as values with errors and units.
##
## .. code-block:: nim
##   import nimboltz/api
##
##   let r = simulate(Setup(mixture: @[(argon, 90.0), (carbonDioxide, 10.0)],
##                          efield: 500.0.V•cm⁻¹, collisions: 4_000_000))
##   echo r.velocity.z               # 41.1 ± 0.122 MicroMeter•NanoSecond⁻¹
##   echo r.spread(r.transverse)     # √(2D/v): 250. ± 33.5 MicroMeter•CentiMeter⁻¹ᐟ²
##
## (By default `collisions` is 2 · 4·10⁷, Magboltz's NMAX = 2, for smaller
## errors.)
##
## A `Setup` is one input deck of Magboltz. `simulate` runs it as Magboltz
## would (energy range, drift run, then the Townsend analyses if ionisation
## or attachment matters) and prints nothing unless `output` is set. By
## default the drift run follows many independent electrons on all cores:
## Magboltz's results within errors, faster. With `fortran = true`,
## `output = stdout` and the same seed, the printout is Magboltz's, bit for
## bit. A run that Magboltz would stop raises a `MagboltzError`.

import units, rm48, deck, conditions, gases, results, magboltz

export units, results, MagboltzError, nscale, WalkerBackend, WalkerSettings, RngKind
export gases.Gas

type
  Setup* = object
    mixture*: seq[(Gas, float64)]  ## up to 6 gases and their shares (%, adding up to 100)
    temperature*: K = 293.15.K
    pressure*: Torr = 760.0.Torr
    efield*: V•cm⁻¹
    bfield*: Tesla = 0.0.T
    angle*: float64 = 0.0         ## between E and B (degrees)
    collisions*: int = 2 * nscale ## real collisions of the drift run (Magboltz:
                                  ## NMAX · 4·10⁷, NMAX = 2..5 for 1 % errors
                                  ## in inelastic gases, ≥ 10 for argon)
    efinal*: eV = 0.0.eV          ## upper end of the energy range; 0: found by
                                  ## raising it until no electron reaches it
    penning*: bool = false        ## Penning transfers
    thermal*: bool = true         ## thermal motion of the gas (else 0 K)
    seed*: int = DefaultSeed      ## of the random numbers (Magboltz's by default)
    fortran*: bool = false        ## Magboltz 11.21 bit for bit: one electron,
                                  ## its methods and numbers, one thread
                                  ## (slower; for cross-checks); `walkers`
                                  ## unused
    threads*: int = 0             ## CPU threads (0: all cores) for the walkers
                                  ## and the avalanches of the analyses; with
                                  ## `fortran`, 0 is one thread, and more split
                                  ## the avalanches, with random numbers of
                                  ## their own
    walkers*: WalkerSettings      ## the drift run (and ELIMIT) with independent
                                  ## electrons on CPU threads (wbCpu) or a GPU
                                  ## (wbGpu; `single`: in single precision):
                                  ## the same results within errors
    output*: File                 ## Magboltz's printout goes here (nil: nowhere)

func deck*(s: Setup): Deck =
  ## the input deck of a setup
  if s.mixture.len notin 1 .. 6:
    raise newException(MagboltzError, "1 to 6 gases, not " & $s.mixture.len)
  result.ngas = s.mixture.len
  result.ipen = ord(s.penning)
  result.ithrm = ord(s.thermal)
  result.efinal = s.efinal
  for k, (gas, percent) in s.mixture:
    result.gas[k + 1] = ord(gas)
    result.percent[k + 1] = percent
  result.tempc = float64(s.temperature - abzero)
  result.torr = s.pressure
  result.efield = s.efield
  result.bmag = s.bfield.to(kGauss)
  result.btheta = s.angle

proc simulate*(s: Setup): Results =
  ## runs a setup as Magboltz would run its deck
  if s.collisions <= 0:
    raise newException(MagboltzError, "no collisions to run")
  var session = initSession(s.output, s.fortran)
  if s.seed != DefaultSeed: session.rng.rm.rm48in(s.seed, 0, 0)
  if s.threads > 0: session.threads = s.threads
  if not s.fortran: session.walkers = s.walkers
  session.run(s.deck, s.collisions)
