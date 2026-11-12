## The program MAGBOLTZ: reads decks from standard input and, for each,
## computes the collision tables of the mixture (MIXER, MIXERT), the energy
## range if the deck leaves it open (ELIMIT), the transport coefficients
## (MONTE) and, when ionisation or attachment matters, the Townsend
## analyses (analyses.nim). Thermal motion of the gas (ITHRM = 1) selects
## the per-gas versions of all steps at compile time (`WarmMixture`).

import std/[math, times, strutils, cpuinfo]
import frt, units, rm48, deck, conditions, tables, mixture, gases, fields, collisions,
       freeflight, elimit, drift, simulation, report, analyses, results, walkers, walkdrift
when defined(nimboltzGpu):
  import gpuwalk
elif defined(nimboltzWebGpu):
  import webgpuwalk

export MagboltzError, Results, RngKind

const
  fmt50 = "(' WARNING ENERGY OUT OF RANGE,INCREASE ELECTRON ENERGY INTEGRATION RANGE ')"
  fmt776 = "('  WARNING NEGATIVE COLLISION FEQUENCY =',D12.3,' IE =',I6,' IF =',I3,' IARRY='," &
           "I5,' EIN=',F7.4)"
  fmt776t = "('  WARNING NEGATIVE COLLISION FEQUENCY =',D12.3,' KGAS =', I2,' IE =',I6," &
            "' IF =',I3,' IARRY=',I5,' EIN=',F7.4)"
  fmt992 = "(/,/,6X," &
           "'WARNING TOO MANY LEVELS IN CALCULATION. CAN INCREASE THE ARRAY SIZES FROM 960 UP TO 1800 MAXIMUM'," &
           "/)"
  # MONTE progress table, one line per segment
  fmt201e = "(/,'    VEL      POS        TIME      ENERGY   COUNT   DIFXX     DIFYY     DIFZZ',/)"
  fmt201a = "(/,'    VEL      POS        TIME      ENERGY   COUNT   DIFXX     DIFYY     DIFZZ  ',/)"
  fmt202a = "(1X,F8.2,2(1X,D10.3),F9.4,1X,I5,1X,3(2X,F8.1))"
  fmt201b = "(/,'   VELZ    VELY  ENERGY   ID  DIFXX   DIFYY   DIFZZ   DIFYZ  DIFLNG  DIFTRN',/)"
  fmt202b = "(2(F7.2,1X),F7.4,1X,I4,6(F7.1,1X))"
  fmt201c = "(/,'   VELZ    VELY    VELX   ENERGY  ID  DIFXX  DIFYY  DIFZZ  DIFYZ  DIFXZ  DIFXY',/)"
  fmt202c = "(3(F7.2,1X),1X,F7.4,I4,6(F7.1))"
  fmtWalkers = "(/,' MONTE WITH',I7,' WALKERS:',I9,' COLLISIONS EACH TO WARM UP,',I9," &
               "' COUNTED',/)"

type
  WalkerBackend* = enum
    wbNone  ## Magboltz's one long trajectory
    wbCpu   ## walkers on CPU threads
    wbGpu   ## walkers on an NVIDIA GPU (built with -d:nimboltzGpu), or with
            ## WebGPU in a web page (-d:nimboltzWebGpu, single precision)

  WalkerSettings* = object
    ## MONTE with many independent electrons (`walkdrift.nim`), and
    ## ELIMIT on them (`walkerEnergyRange`) when EFINAL is not given: the
    ## same estimates within errors, not the same numbers as Magboltz
    backend*: WalkerBackend = wbCpu
    n*: int        ## walkers; 0: `walkersPerThread` per CPU thread, 16384 on the GPU
    threads*: int  ## CPU threads; 0: the session's
    single*: bool  ## GPU: in single precision (see single.nim)
    rng*: RngKind = rkXoshiro
      ## the walkers' random numbers (xoshiro: xoshiro256+ on the CPU,
      ## xoshiro128** on a GPU, xoshiro128+ in single precision), and the
      ## avalanches' with exact bounds

  Session* = object
    ## what lives from one deck to the next
    rng*: MagRng   ## Magboltz's random numbers: one stream for all decks
    fortran*: bool
      ## Magboltz's numbers bit for bit: one electron for the drift run,
      ## Magboltz's null collisions and formulas everywhere, one thread
      ## (`Run.fortran`); else walkers, exact bounds per time window and
      ## threads, the same results within errors and faster
    threads*: int  ## for the walkers and the avalanche runs (`Run.threads`)
    walkers*: WalkerSettings
    timing*: bool  ## times of the steps of each deck to stderr
    output*: File  ## Magboltz's printout; nil: none
    cold: ref ColdMixture
    warm: ref WarmMixture
    rings: array[FieldGeometry, array[bool, Ring]]
      ## The position arrays of the MONTE routines are static, so a run
      ## starts with the positions left by the previous run of the same
      ## routine; MONTEC reads them when NCOLM exceeds the collisions of the
      ## first two segments (reduced nscale only). One ring per routine.

const walkersPerThread* = 100
  ## on CPUs: more walkers are shorter ones, whose warm-up costs more

proc defaultThreads*(): int =
  ## all cores (one without threads)
  when compileOption("threads"): max(countProcessors(), 1) else: 1

proc initSession*(output: File = stdout, fortran = false): Session =
  ## walkers on all cores, or with `fortran` Magboltz bit for bit
  Session(rng: initMagRng(output), fortran: fortran,
          threads: (if fortran: 1 else: defaultThreads()),
          walkers: WalkerSettings(backend: if fortran: wbNone else: wbCpu),
          output: output)

# --- MIXER, MIXERT ------------------------------------------------------------------

proc crossSections(r: var Run): seq[GasXsec] =
  ## energy grid and the cross sections of each gas
  let d = r.deck
  r.grid = energyGrid(r.efinal, nEnergyBins)
  result = newSeq[GasXsec](d.ngas)
  let inp = GasInput(e: r.grid.e, naniso: naniso, kT: r.cond.kT, efinal: r.efinal,
                     estep: r.estep, penning: d.penning, torr: d.torr,
                     output: r.output, tempc: d.tempc,
                     mixture: mixtureOf(d.gas, d.percent))
  for k, g in result.mpairs:
    crossSections(d.gas[k + 1], inp, g)
  for k in 0 ..< d.ngas:
    r.names[k + 1] = fstr(result[k].name, 25)

proc mix*(r: var Run[ColdMixture]) =
  ## MIXER: the collision table of the mixture
  let gases = r.crossSections()
  let m = r.mix
  var descr: array[1..maxMixProcesses, string]
  var negative: seq[NegativeFrequency]
  try:
    buildMixture(m.table, descr, gases, r.cond.van.toOpenArray(1, r.deck.ngas), r.efinal,
                 r.grid, r.fcion, r.fcatt, negative)
  except ValueError:
    r.output.fwrite(fmt992)
    r.stop("too many levels in the mixture")
  for n in negative:
    r.output.fwrite(fmt776, n.cf, n.ie, n.process, m.table.iarry[n.process],
           m.table.ein[n.process])
    r.warnings.add "negative collision frequency, energy bin " & $n.ie
  for i in 1 .. m.table.last:
    m.descr[i] = fstr(descr[i], 50)
  m.isotropic = not anisotropic(gases)
  m.bounds = nullBounds(m.table, r.deck.efield, tmax, small, r.efinal, r.estep,
                        const1, const3)

proc mix*(r: var Run[WarmMixture]) =
  ## MIXERT: one collision table per gas
  let gases = r.crossSections()
  let m = r.mix
  var descr: array[1..maxGases, array[1..maxGasProcesses, string]]
  var negative: seq[NegativeFrequency]
  buildPerGas(m.tables, descr, gases, r.cond.van.toOpenArray(1, r.deck.ngas), r.efinal,
              r.grid, r.fcion, r.fcatt, negative)
  for n in negative:
    r.warnings.add "negative collision frequency, gas " & $n.gas & ", energy bin " & $n.ie
    r.output.fwrite(fmt776t, n.cf, n.gas, n.ie, n.process, m.tables[n.gas].iarry[n.process],
           m.tables[n.gas].ein[n.process])
  for k in 1 .. gases.len:
    for i in 1 .. m.tables[k].last:
      m.descr[k][i] = fstr(descr[k][i], 50)
  m.isotropic = not anisotropic(gases)
  m.gas = thermalGas(m.tables, gases, r.cond.kT, conditions.echarg, conditions.emass)

# --- ELIMIT -------------------------------------------------------------------------

proc energyLimitExceeded[G: static FieldGeometry; M](r: var Run[M], rng: var MagRng): bool =
  let f = r.field(G)
  let start = initElectron[G](f, r.estart, r.theta, startPhi)
  when thermalMotion(M):
    var target = initWarmGas(r.mix.gas, rng)
    result = energyLimitExceeded[G](f, r.mix.tables, target, r.estep, start, r.cond.nmax,
                                    rng, r.theta)
  else:
    var target = initColdGas(r.mix.bounds, r.efinal)
    result = energyLimitExceeded[G](f, r.mix.table, target, r.estep, start, r.cond.nmax,
                                    rng, r.theta)
    r.mix.bounds = target.bounds

func walkerCount(wk: WalkerSettings): int =
  if wk.n > 0: max(wk.n, 20)
  elif wk.backend == wbGpu: 16384
  else: max(walkersPerThread * wk.threads, minWalkers)

func nullBands*(r: Run[WarmMixture]): array[1..maxGases, NullBounds] =
  ## walkers with thermal motion: null-collision bounds of each gas in the
  ## 8 bands of starting energy, as MIXER has them for the mixture at 0 K
  for k in 1 .. r.deck.ngas:
    result[k] = nullBounds(r.mix.tables[k], r.deck.efield, tmax, small, r.efinal, r.estep,
                           const1, const3)

proc initWalkers[G: static FieldGeometry; M; R](r: Run[M], f: Field, n: int,
                                                rng: var MagRng): auto =
  ## walkers starting as MONTE's electron, with random numbers (R) of their own
  let start = initElectron[G](f, r.estart, r.theta, startPhi)
  when thermalMotion(M):
    let bands = r.nullBands
    var ws = newSeq[Walker[WarmGas, R]](n)
    for w in ws.mitems:
      w = initWalker(start, r.mix.gas, bands, r.efinal, int(drand48(rng) * 9e8), R)
  else:
    let target = initColdGas(r.mix.bounds, r.efinal)
    var ws = newSeq[Walker[ColdGas, R]](n)
    for w in ws.mitems: w = initWalker(start, target, int(drand48(rng) * 9e8), R)
  ws

template tablesOf(r: Run): untyped =
  when thermalMotion(typeof(r.mix[])): r.mix.tables else: r.mix.table

template withWalkers(r: Run, wk: WalkerSettings, G: static FieldGeometry, f: Field,
                     ws: seq, b, body: untyped) =
  ## `body` with walker backend `b` of the settings
  type T = typeof(ws[0].target)
  type R = typeof(ws[0].rng)
  let fm = flightMaxima(r.tablesOf)
  case wk.backend
  of wbGpu:
    when defined(nimboltzGpu) and (R is Xoshiro128 or R is Ranmar and defined(gpuRanmar)) or
         defined(nimboltzWebGpu) and (R is Xoshiro128Plus or R is Ranmar and defined(gpuRanmar)):
      var b = initGpuWalkers[G, T, R](ws, f, r.tablesOf, fm[], r.estep, r.deck.penning,
                                      single = wk.single)
      body
    elif defined(nimboltzGpu) or defined(nimboltzWebGpu):
      r.stop("walkers on a GPU: xoshiro's random numbers (Ranmar's with -d:gpuRanmar)")
    else:
      r.stop("walkers on a GPU: NimBoltz was built without -d:nimboltzGpu")
  else:
    var b = initCpuWalkers[G, T, R](ws, f, addr r.tablesOf, r.estep, r.deck.penning,
                                    max(wk.threads, 1), addr fm[])
    body

template withWalkerRng(wk: WalkerSettings, R, body: untyped) =
  ## body with `R` the walkers' generator: on a GPU xoshiro128** (32-bit
  ## arithmetic, as WGSL has it) or, in single precision, xoshiro128+
  when defined(nimboltzGpu) or defined(nimboltzWebGpu):
    if wk.backend == wbGpu and wk.rng == rkXoshiro:
      if wk.single or defined(nimboltzWebGpu):  # (WebGPU: single precision only)
        type R {.inject.} = Xoshiro128Plus
        body
      else:
        when defined(nimboltzGpu):
          type R {.inject.} = Xoshiro128StarStar
          body
    else:
      withRng(wk.rng, R, body)
  else:
    withRng(wk.rng, R, body)

func firstEfinal(d: Deck): eV =
  ## EFINAL not given: ELIMIT starts at 0.5 eV (8 eV if E/p > 15 V/cm/Torr
  ## at 20 °C)
  let eop = d.efield * (d.tempc + 273.15) / (d.torr * 293.15)
  if eop > 15.0.V•cm⁻¹•Torr⁻¹: 8.0.eV else: 0.5.eV

proc findEnergyRange(r: var Run, rng: var MagRng) =
  ## EFINAL not given: raise it by √2 from `firstEfinal` until ELIMIT never
  ## reaches the last energy bin
  let d = r.deck
  r.efinal = firstEfinal(d)
  while true:
    r.mix()
    var exceeded: bool
    withGeometry(geometryOf(d.bmag, d.btheta), G):
      exceeded = energyLimitExceeded[elimitGeometry(G), typeof(r.mix[])](r, rng)
    if not exceeded: break
    r.efinal = r.efinal * sqrt(2.0)

# --- MONTE --------------------------------------------------------------------------

proc progress[G: static FieldGeometry; T](o: File, d: Drift[G, T], seg: int) =
  let w = d.segmentVelocity
  let df = d.segmentDiffusion
  when G in {fgElectric, fgParallel}:
    if seg == 1: o.fwrite(when G == fgElectric: fmt201e else: fmt201a)
    o.fwrite(fmt202a, w.z, d.position.z, d.time, d.ave, d.ringFills,
           df.xx, df.yy, df.zz)
  elif G == fgCrossed:
    if seg == 1: o.fwrite(fmt201b)
    o.fwrite(fmt202b, w.z, w.y, d.ave, d.ringFills, df.xx, df.yy,
           df.zz, df.yz, d.difln, d.diftr)
  else:
    if seg == 1: o.fwrite(fmt201c)
    o.fwrite(fmt202c, w.z, w.y, w.x, d.ave, d.ringFills, df.xx,
           df.yy, df.zz, df.yz, df.xz, df.xy)

proc monte[G: static FieldGeometry; M](r: var Run[M], ring: var Ring, rng: var MagRng) =
  ## MONTE, MONTEA, MONTEB, MONTEC and the thermal versions
  let f = r.field(G)
  let start = initElectron[G](f, r.estart, r.theta, startPhi)
  let dc = r.deck
  let collisions = r.cond.nmax div nSegments
  if ring.isEmpty: ring = initRing()
  when thermalMotion(M):
    let tcf = r.mix[].totalTcf(dc.ngas)
    var d = initDrift[G, WarmGas](f, start, initWarmGas(r.mix.gas, rng), r.estep, dc.ngas,
                                  dc.penning, r.mix.isotropic, r.cond.dec, dc.btheta)
    template tabs: untyped = r.mix.tables
  else:
    template tcf: untyped = r.mix.table.tcf
    var d = initDrift[G, ColdGas](f, start, initColdGas(r.mix.bounds, r.efinal), r.estep,
                                  dc.ngas, dc.penning, r.mix.isotropic, r.cond.dec,
                                  dc.btheta)
    template tabs: untyped = r.mix.table
  for seg in 1 .. nSegments:
    d.runSegment(tabs, seg, collisions, ring, rng)
    d.endSegment(seg, r.grid.e, tcf)
    r.output.progress(d, seg)
    when G == fgElectric:
      # if the last energy bin is reached often, the table range is too small
      if d.spectrum[nEnergyBins] > 1000.0 * float64(seg):
        r.output.fwrite(fmt50)
        r.stop("energy out of range: electrons reach the end of the energy range " &
               "often, increase EFINAL")
  r.transport = d.finish
  when not thermalMotion(M):
    r.mix.bounds = d.target.bounds
  r.theta = d.theta
  r.tmax = d.tmax
  r.nnull = d.nnull
  r.collisions = r.cond.nmax
  r.last = Collisions[M](tally: d.tally, st: d.time, spectrum: d.spectrum)

proc retarget[M](r: Run[M], ws: var seq) =
  ## walkers after MIXER for another EFINAL: its null-collision bounds
  when thermalMotion(M):
    let bands = r.nullBands
    for w in ws.mitems:
      w.target.gas = r.mix.gas
      w.target.bands = bands
      w.target.bandWidth = r.efinal / float64(bands[1].len)
  else:
    let target = initColdGas(r.mix.bounds, r.efinal)
    for w in ws.mitems: w.target = target

proc walkerEnergyRange[G: static FieldGeometry; M; R](r: var Run[M], rng: var MagRng,
                                                       wk: WalkerSettings): auto =
  ## MIXER and ELIMIT on walkers (EFINAL not given): EFINAL raised by √2
  ## from `firstEfinal` until no walker reaches the last energy bin in
  ## ELIMIT's number of collisions (all walkers together, once warm; see
  ## `walkdrift.energyLimitExceeded`). Each try starts from MONTE's
  ## electron (in a range too small the walkers' energies run away).
  ## Returns the warm walkers of the last try, for MONTE, and their warm-up.
  let f = r.field(G)
  r.efinal = firstEfinal(r.deck)
  r.mix()
  let fresh = initWalkers[G, M, R](r, f, walkerCount(wk), rng)
  var ws = fresh
  while true:
    var wu: WarmUp
    withWalkers(r, wk, G, f, ws, b):
      wu = b.energyLimitExceeded(ws, budget = r.cond.nmax,
                                 sample = sampleSize(elimitGeometry(G), r.cond.nmax))
    if not wu.exceeded: return (ws, wu)
    r.efinal = r.efinal * sqrt(2.0)
    r.mix()
    # back to the start, with the random numbers where they are (seeding
    # many generators takes a while); walkers dropped on CPUs come back
    let kept = ws.len
    ws.setLen(fresh.len)
    for k in 0 ..< ws.len:
      let rng = if k < kept: ws[k].rng else: fresh[k].rng
      ws[k] = fresh[k]
      ws[k].rng = rng
    r.retarget(ws)

type
  WalkerDrift = object
    ## MONTE's run on walkers, before its printout
    res: WalkerRun
    perWalker: int32

proc walkerDrift[G: static FieldGeometry; M; T, R](r: var Run[M], wk: WalkerSettings,
                                                   ws: var seq[Walker[T, R]],
                                                   wu: WarmUp, stopEarly = false): WalkerDrift =
  ## the drift run of warm walkers (`wu`); `stopEarly`: as soon as the energy
  ## range shows too small (for a run that raises EFINAL then)
  let f = r.field(G)
  withWalkers(r, wk, G, f, ws, b):
    # the main run at least as long as the warm-up (see walkdrift.nim)
    while b.dropsWalkers and ws.len div 2 >= minWalkers and
          r.cond.nmax div ws.len < wu.collisions:
      b.drop(ws, ws.len div 2)
    result.perWalker = int32(max(max(r.cond.nmax div ws.len, wu.collisions), 2))
    result.res = driftRun[G, typeof(b), T, R](
      b, ws, result.perWalker, chunksOf(result.perWalker, wu.collisions), r.deck.btheta,
      abortAbove = (if stopEarly: 1000.0 / float64(r.cond.nmax div nSegments) else: 0.0),
      anyBefore = sampleSize(elimitGeometry(G), r.cond.nmax))

func lastBinPerSegment(r: Run, d: WalkerDrift, n: int): float64 =
  ## collisions in the last energy bin per NMAX / 10 of the drift run
  ## (MONTE's segments)
  d.res.tallies.spectrum[nEnergyBins] / (float64(n) * float64(d.perWalker)) *
    float64(r.cond.nmax div nSegments)

proc printWalkerDrift[G: static FieldGeometry; M; T, R](r: var Run[M], ws: seq[Walker[T, R]],
                                                     wu: WarmUp, d: WalkerDrift) =
  ## the outcome of `walkerDrift` into the run, as MONTE's
  let res = d.res
  let perWalker = d.perWalker
  r.output.fwrite(fmtWalkers, ws.len, wu.collisions, perWalker)
  if not wu.converged:
    r.warnings.add "the walkers may not have forgotten their start after " &
                   $wu.collisions & " collisions each (mean energy " & $wu.previous.float &
                   " then " & $wu.energy.float & " eV)"
  when G == fgElectric:
    # MONTE's test, for its 10 segments of NMAX/10 collisions
    let counted = float64(ws.len) * float64(perWalker) / float64(r.cond.nmax)
    if res.tallies.spectrum[nEnergyBins] > 1000.0 * float64(nSegments) * counted:
      r.output.fwrite(fmt50)
      r.stop("energy out of range: electrons reach the end of the energy range " &
             "often, increase EFINAL")
  r.transport = res.transport
  r.tmax = res.tmax
  r.nnull = res.nulls
  r.collisions = ws.len * int(perWalker)
  var c = Collisions[M](st: res.st, spectrum: res.tallies.spectrum)
  when thermalMotion(M):
    for g in 1 .. r.deck.ngas:
      for i in 1 .. r.mix.tables[g].last:
        let k = int(res.tallies.processes[(g - 1) * maxGasProcesses + i])
        c.tally[g].processes[i] = k
        c.tally[g].kinds[r.mix.tables[g].iarry[i]] += k
  else:
    # (the analyses keep Magboltz's null-collision bounds, not the walkers')
    for i in 1 .. r.mix.table.last:
      let k = int(res.tallies.processes[i])
      c.tally.processes[i] = k
      c.tally.kinds[r.mix.table.iarry[i]] += k
  r.last = c

proc monte[G: static FieldGeometry; M; R](r: var Run[M], rng: var MagRng, wk: WalkerSettings) =
  ## MONTE, MONTEA, MONTEB, MONTEC and the thermal versions with walkers,
  ## which start as Magboltz's electron and warm up
  let f = r.field(G)
  var ws = initWalkers[G, M, R](r, f, walkerCount(wk), rng)
  var wu: WarmUp
  withWalkers(r, wk, G, f, ws, b):
    wu = b.warmUp(ws, budget = r.cond.nmax)
  type T = typeof(ws[0].target)
  let d = walkerDrift[G, M, T, R](r, wk, ws, wu)
  printWalkerDrift[G, M, T, R](r, ws, wu, d)

# --- one deck -----------------------------------------------------------------------

proc run[M: Mixture](s: var Session, mix: ref M, d: Deck, collisions: int): Results =
  var r = initRun(d, mix, s.output)
  if collisions > 0: r.cond.nmax = collisions
  r.threads = s.threads
  r.fortran = s.fortran
  var wk = s.walkers
  if s.fortran: wk.backend = wbNone
  r.generator = wk.rng
  if wk.threads <= 0: wk.threads = s.threads
  when not compileOption("threads"): wk.threads = 1
  var clock = epochTime()
  template phase(step: string) =
    if s.timing:
      let now = epochTime()
      stderr.writeLine "NimBoltz: ", alignLeft(step, 20), formatFloat(now - clock, ffDecimal, 3), " s"
      clock = now
  if wk.backend != wbNone and d.efinal <= 0.0.eV:
    # ELIMIT on the walkers, which go on to MONTE warm
    withGeometry(monteGeometry(d), G):
      withWalkerRng(wk, R):
        var (ws, wu) = walkerEnergyRange[G, M, R](r, s.rng, wk)
        phase "MIXER and ELIMIT"
        type T = typeof(ws[0].target)
        var dr = walkerDrift[G, M, T, R](r, wk, ws, wu, stopEarly = true)
        # Where MONTE stops (more than 1000 collisions per segment in the
        # last energy bin: a range too small), walkers raise EFINAL and run
        # again. Short walkers can miss a tail that builds up slowly (an
        # electron escaping into argon heats up over 10⁵ collisions), which
        # Magboltz's ELIMIT sees on its one long trajectory. The drift run
        # stops at the first checkpoint with a chunk over the limit, or with
        # any collision in the last bin within ELIMIT's number of
        # collisions (`driftRun`).
        while dr.res.aborted or r.lastBinPerSegment(dr, ws.len) > 1000.0:
          if s.timing: stderr.writeLine "NimBoltz: EFINAL ", r.efinal.float, " eV too small for MONTE"
          r.efinal = r.efinal * sqrt(2.0)
          r.mix()
          r.retarget(ws)
          for w in ws.mitems: w.reset()
          dr = walkerDrift[G, M, T, R](r, wk, ws, wu, stopEarly = true)
        r.printHeader()
        printWalkerDrift[G, M, T, R](r, ws, wu, dr)
  else:
    if d.efinal > 0.0.eV:
      r.mix()
      phase "MIXER"
    else:
      r.findEnergyRange(s.rng)
      phase "MIXER and ELIMIT"
    r.printHeader()
    withGeometry(monteGeometry(d), G):
      if wk.backend == wbNone: monte[G, M](r, s.rings[G][thermalMotion(M)], s.rng)
      else:
        withWalkerRng(wk, R): monte[G, M, R](r, s.rng, wk)
  phase "MONTE"
  r.printTransport()
  if r.needsAnalyses:
    r.analyse(s.rng)
    phase "analyses"
  r.printFrequencies()
  result = results(r)
  phase "printout" 

proc run*(s: var Session, d: Deck, collisions = 0): Results =
  ## one deck; `collisions` > 0 replaces NMAX · nscale. Raises a
  ## `MagboltzError` where Magboltz stops the program (an input error, an
  ## energy range too small, ...).
  let e = setupError(d)
  if e != seNone:
    s.output.printSetupError(d, e)
    raise newException(MagboltzError, case e
      of seGasInput: "error in the gas input: 1 to 6 gases, shares adding up to 100 %"
      else: "NMAX too large")
  if d.thermal:
    if s.warm == nil: new s.warm
    s.run(s.warm, d, collisions)
  else:
    if s.cold == nil: new s.cold
    s.run(s.cold, d, collisions)

proc magboltz*(s: var Session) =
  ## all decks of standard input, until one makes Magboltz stop
  var d: Deck
  try:
    while readDeck(d):
      discard s.run(d)
  except MagboltzError:
    discard
  fstop()
