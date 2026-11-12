## The state of one run of Magboltz (one input deck), shared by the
## driver (magboltz.nim), the analyses (analyses.nim) and the printout
## (report.nim).
##
## Two kinds of mixture, chosen at compile time: `ColdMixture` for the gas
## at 0 K (MIXER, MONTE, ...), `WarmMixture` for thermal motion of the gas
## (ITHRM = 1: MIXERT, MONTET, ...).

import std/options
import units, xoshiro, deck, conditions, tables, mixture, thermal, freeflight, fields, drift,
       townsend

type
  MagboltzError* = object of CatchableError
    ## a run that cannot go on, where Magboltz stops the program; its
    ## message in Magboltz's words is in the run's output already

  ColdMixture* = object
    ## gas at 0 K: one collision table for the whole mixture
    table*: MixtureTable
    descr*: array[1..maxMixProcesses, string]  ## process descriptions (50 characters)
    bounds*: NullBounds  ## null-collision bounds; ELIMIT and MONTE raise them
    isotropic*: bool     ## no anisotropic process in any gas (NISO = 0)

  WarmMixture* = object
    ## thermal motion of the gas: one collision table per gas
    tables*: GasTables
    descr*: array[1..maxGases, array[1..maxGasProcesses, string]]
    gas*: ThermalGas     ## per-gas null-collision bounds and thermal speeds
    isotropic*: bool

  Mixture* = ColdMixture | WarmMixture

  Collisions*[M] = object
    ## real collisions of the last Monte Carlo run (MONTE or an avalanche
    ## run), for the collision frequencies (COLF, COLFT) and FRIEDLAND
    when M is ColdMixture:
      tally*: Tally[maxMixProcesses]
    else:
      tally*: GasTallies
    st*: ps          ## time followed
    ttots*: ps
      ## avalanche runs: time summed from the start of each electron (TTOTS);
      ## 0 after MONTE, then `st` is used
    spectrum*: array[1..nEnergyBins, float64]  ## real collisions per energy bin

  Frequencies* = object
    ## real collision frequencies of a run, by kind
    total*, elastic*, ionisation*, attachment*, inelastic*: ps⁻¹
    n*: int  ## real collisions

  Run*[M] = object
    deck*: Deck
    cond*: Conditions
    mix*: ref M
    efinal*: eV           ## upper end of the energy range, maybe found by ELIMIT
    grid*: EnergyGrid
    fcion*, fcatt*: array[1..nEnergyBins, ps⁻¹]
      ## total ionisation and attachment frequencies per energy bin
    names*: array[1..maxGases, string]  ## gas names (25 characters)
    theta*: float64
      ## polar angle of the last scattering: each Monte Carlo run starts in
      ## the direction where the previous one ended (THETA of /SETP/)
    collisions*: int      ## MONTE: real collisions counted (NMAX; with
                          ## walkers on a GPU maybe more)
    nnull*: int64         ## MONTE: null collisions
    tmax*: ps             ## MONTE: longest flight
    transport*: Transport
      ## MONTE's results; the analyses replace the Townsend coefficients
    last*: Collisions[M]
    threads*: int
      ## avalanche runs on this many threads, each with random numbers of
      ## its own (1: Magboltz's single stream, bit for bit)
    fortran*: bool
      ## Magboltz's null collisions and formulas in the analyses, for its
      ## numbers bit for bit; else exact bounds per time window
    generator*: RngKind
      ## random numbers of the avalanches with exact bounds
    output*: File          ## Magboltz's printout; nil: none
    warnings*: seq[string]
    # the analyses, if run (`needsAnalyses`)
    steadyState*: Option[SteadyState]       ## B = 0
    pulsedTownsend*: Option[PulsedTownsend]
    timeOfFlight*: Option[TimeOfFlight]
    friedland*: Option[Friedland]

template thermalMotion*(M: typedesc): bool = M is WarmMixture

func initRun*[M: Mixture](d: Deck, mix: ref M, output: File): Run[M] =
  Run[M](deck: d, cond: conditions(d), mix: mix, efinal: d.efinal, theta: startTheta,
         threads: 1, output: output)

proc stop*(r: Run, msg: string) {.noreturn.} =
  ## where Magboltz stops the program
  raise newException(MagboltzError, msg)

func estart*(r: Run): eV =
  ## energy of the first electron
  r.efinal / 50.0

func estep*(r: Run): eV = r.grid.estep

func monteGeometry*(d: Deck): FieldGeometry =
  ## Routine chosen by Magboltz for MONTE and the analyses: unlike ELIMIT
  ## (`geometryOf`), B ∥ E at -180° goes to MONTEC
  if d.bmag.float == 0.0: fgElectric
  elif d.btheta == 0.0 or d.btheta == 180.0: fgParallel
  elif d.btheta == 90.0: fgCrossed
  else: fgOblique

func field*(r: Run, G: FieldGeometry): Field =
  initField(G, r.deck.efield, const1, const2, const3, const5, r.cond.wb,
            r.cond.eovb, r.deck.btheta)

func totalTcf*(m: WarmMixture, ngas: int): array[1..nEnergyBins, ps⁻¹] =
  ## total collision frequency of the mixture per energy bin
  for ik in 1 .. nEnergyBins:
    var tcfsum = 0.0.ps⁻¹
    for ki in 1 .. ngas: tcfsum = tcfsum + m.tables[ki].tcf[ik]
    result[ik] = tcfsum

func count(t: Tally[maxMixProcesses], gas, kind: int): int = t.kinds[5 * (gas - 1) + kind]
func count(t: GasTallies, gas, kind: int): int = t[gas].kinds[kind]

func frequencies*(c: Collisions, ngas: int): Frequencies =
  ## COLF, COLFT. Kinds: 1 elastic, 2 ionisation, 3 attachment, 4 inelastic,
  ## 5 superelastic (counted as inelastic)
  var ninel, nela, natt, nion = 0
  for j in 1 .. ngas:
    ninel = ninel + c.tally.count(j, 4) + c.tally.count(j, 5)
    nela = nela + c.tally.count(j, 1)
    natt = natt + c.tally.count(j, 3)
    nion = nion + c.tally.count(j, 2)
  result.n = nela + ninel + natt + nion
  let t = if c.ttots.float == 0.0: c.st else: c.ttots
  result.total = float64(result.n) / t
  result.inelastic = float64(ninel) / t
  result.elastic = float64(nela) / t
  result.ionisation = float64(nion) / t
  result.attachment = float64(natt) / t

iterator processes*(r: Run[ColdMixture], gas: int): tuple[descr: string, n: int] =
  ## the processes of a gas with their real collisions in the last run
  template tab: untyped = r.mix.table  # (a `let` would copy the table)
  for k in 1 .. tab.last:
    if tab.iarry[k] <= 5 * gas and tab.iarry[k] > 5 * (gas - 1):
      yield (r.mix.descr[k], r.last.tally.processes[k])

iterator processes*(r: Run[WarmMixture], gas: int): tuple[descr: string, n: int] =
  for k in 1 .. r.mix.tables[gas].last:
    yield (r.mix.descr[gas][k], r.last.tally[gas].processes[k])
