## The program MAGBOLTZ: reads decks from standard input and, for each,
## computes the collision tables of the mixture (MIXER, MIXERT), the energy
## range if the deck leaves it open (ELIMIT), the transport coefficients
## (MONTE) and, when ionisation or attachment matters, the Townsend
## analyses (analyses.nim). Thermal motion of the gas (ITHRM = 1) selects
## the per-gas versions of all steps at compile time (`WarmMixture`).

import std/math
import frt, units, rm48, deck, conditions, tables, mixture, gases, fields, collisions,
       freeflight, elimit, drift, simulation, report, analyses, results, walkers, walkdrift
when defined(nimboltzGpu):
  import gpuwalk

export MagboltzError, Results

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
    wbGpu   ## walkers on an NVIDIA GPU (built with -d:nimboltzGpu)

  WalkerSettings* = object
    ## MONTE with many independent electrons (`walkdrift.nim`); ELIMIT
    ## keeps Magboltz's one trajectory (cheap, and so the same energy range) (thermal runs keep the one trajectory): the same
    ## estimates within errors, not the same numbers as Magboltz
    backend*: WalkerBackend
    n*: int        ## walkers; 0: 2000 on the CPU, 16384 on the GPU
    threads*: int  ## CPU threads

  Session* = object
    ## what lives from one deck to the next
    rng*: MagRng   ## Magboltz's random numbers: one stream for all decks
    threads*: int  ## for the avalanche runs (`Run.threads`)
    walkers*: WalkerSettings
    output*: File  ## Magboltz's printout; nil: none
    cold: ref ColdMixture
    warm: ref WarmMixture
    rings: array[FieldGeometry, array[bool, Ring]]
      ## The position arrays of the MONTE routines are static, so a run
      ## starts with the positions left by the previous run of the same
      ## routine; MONTEC reads them when NCOLM exceeds the collisions of the
      ## first two segments (reduced nscale only). One ring per routine.

proc initSession*(output: File = stdout): Session =
  Session(rng: initMagRng(output), threads: 1, output: output)

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
  else: 2000

proc initWalkers[G: static FieldGeometry; M](r: Run[M], f: Field, n: int,
                                             rng: var MagRng): auto =
  ## walkers starting as MONTE's electron, with random numbers of their own
  let start = initElectron[G](f, r.estart, r.theta, startPhi)
  when thermalMotion(M):
    var ws = newSeq[Walker[WarmGas]](n)
    for w in ws.mitems: w = initWalker(start, r.mix.gas, int(drand48(rng) * 9e8))
  else:
    let target = initColdGas(r.mix.bounds, r.efinal)
    var ws = newSeq[Walker[ColdGas]](n)
    for w in ws.mitems: w = initWalker(start, target, int(drand48(rng) * 9e8))
  ws

template tablesOf(r: Run): untyped =
  when thermalMotion(typeof(r.mix[])): r.mix.tables else: r.mix.table

template withWalkers(r: Run, wk: WalkerSettings, G: static FieldGeometry, f: Field,
                     ws: seq, b, body: untyped) =
  ## `body` with walker backend `b` of the settings
  type T = typeof(ws[0].target)
  case wk.backend
  of wbGpu:
    when defined(nimboltzGpu):
      var b = initGpuWalkers[G, T](ws, f, r.tablesOf, r.estep, r.deck.penning)
      body
    else:
      r.stop("walkers on a GPU: NimBoltz was built without -d:nimboltzGpu")
  else:
    var b = initCpuWalkers[G, T](ws, f, addr r.tablesOf, r.estep, r.deck.penning,
                                 max(wk.threads, 1))
    body

proc findEnergyRange(r: var Run, rng: var MagRng) =
  ## EFINAL not given: start at 0.5 eV (8 eV if E/p > 15 V/cm/Torr at 20 °C)
  ## and raise it by √2 until ELIMIT never reaches the last energy bin
  let d = r.deck
  let eop = d.efield * (d.tempc + 273.15) / (d.torr * 293.15)
  r.efinal = if eop > 15.0.V•cm⁻¹•Torr⁻¹: 8.0.eV else: 0.5.eV
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
  r.last = Collisions[M](tally: d.tally, st: d.time, spectrum: d.spectrum)

proc monte[G: static FieldGeometry; M](r: var Run[M], rng: var MagRng, wk: WalkerSettings) =
  ## MONTE, MONTEA, MONTEB, MONTEC and the thermal versions with walkers
  let f = r.field(G)
  let n = walkerCount(wk)
  var ws = initWalkers[G, M](r, f, n, rng)
  var wu: WarmUp
  var res: WalkerRun
  var perWalker: int32
  withWalkers(r, wk, G, f, ws, b):
    wu = b.warmUp(ws, budget = r.cond.nmax div 2)
    # the main run at least as long as the warm-up (see walkdrift.nim)
    while ws.len div 2 >= minWalkers and r.cond.nmax div ws.len < wu.collisions:
      b.drop(ws, ws.len div 2)
    perWalker = int32(max(r.cond.nmax div ws.len, 2))
    res = driftRun[G, typeof(b), typeof(ws[0].target)](b, ws, perWalker, r.deck.btheta)
  r.output.fwrite(fmtWalkers, ws.len, wu.collisions, perWalker)
  if not wu.converged:
    r.warnings.add "the walkers may not have forgotten their start after " &
                   $wu.collisions & " collisions each (mean energy " & $wu.previous.float &
                   " then " & $wu.energy.float & " eV)"
  when G == fgElectric:
    # MONTE's test, for its 10 segments of NMAX/10 collisions
    if res.tallies.spectrum[nEnergyBins] > 1000.0 * float64(nSegments):
      r.output.fwrite(fmt50)
      r.stop("energy out of range: electrons reach the end of the energy range " &
             "often, increase EFINAL")
  r.transport = res.transport
  r.tmax = res.tmax
  r.nnull = res.nulls
  var c = Collisions[M](st: res.st, spectrum: res.tallies.spectrum)
  when thermalMotion(M):
    for g in 1 .. r.deck.ngas:
      for i in 1 .. r.mix.tables[g].last:
        let k = int(res.tallies.processes[(g - 1) * maxGasProcesses + i])
        c.tally[g].processes[i] = k
        c.tally[g].kinds[r.mix.tables[g].iarry[i]] += k
  else:
    r.mix.bounds = res.bounds
    for i in 1 .. r.mix.table.last:
      let k = int(res.tallies.processes[i])
      c.tally.processes[i] = k
      c.tally.kinds[r.mix.table.iarry[i]] += k
  r.last = c

# --- one deck -----------------------------------------------------------------------

proc run[M: Mixture](s: var Session, mix: ref M, d: Deck, collisions: int): Results =
  var r = initRun(d, mix, s.output)
  if collisions > 0: r.cond.nmax = collisions
  r.threads = s.threads
  let wk = s.walkers
  if d.efinal > 0.0.eV: r.mix()
  else: r.findEnergyRange(s.rng)  # one trajectory, as Magboltz, with walkers too
  r.printHeader()
  withGeometry(monteGeometry(d), G):
    if wk.backend == wbNone: monte[G, M](r, s.rings[G][thermalMotion(M)], s.rng)
    else: monte[G, M](r, s.rng, wk)
  r.printTransport()
  if r.needsAnalyses: r.analyse(s.rng)
  r.printFrequencies()
  results(r)

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
