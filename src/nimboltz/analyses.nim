## The analyses after MONTE when ionisation or attachment matters (ALPCALC,
## ALPCALCT, ALPCLCA, ALPCLCAT, ALPCLCB, ALPCLCBT, ALPCLCC, ALPCLCCT of
## Magboltz): avalanche runs (MONTEFT*, MONTEFD and their thermal versions,
## nimboltz/cascade) and the analyses of their sampling planes (SST,
## PT/PTG/PTH, TOF/TOFG/TOFH, FRIEDLAND, nimboltz/townsend).
##
## Each driver chooses the sampling planes from MONTE's Townsend
## coefficients, runs the avalanches (space planes for B = 0, then time
## planes), prints the per-plane tables and the results, and replaces the
## Townsend coefficients and their errors in the run's transport results
## (B = 0: from SST; otherwise from TOF). The eight drivers and ten
## avalanche routines print with slightly different formats and differ in
## a few details; these are kept, marked with the routine they belong to.

import std/[math, options]
import frt, units, rm48, xoshiro, tables, vec3, fields, freeflight, cascade, townsend
from walkers import flightMaxima
import conditions, simulation, report

type
  Formats = object
    stackFull: string      ## FORMAT 546
    totals: string         ## FORMAT 878
    negative: string       ## FORMAT 913 or 999
    negativeKind: int      ## 1: (ITER, E, E1, AP, BP, DCZ1); 2: (ITER, DZ, E, COSWT,
                           ## SINWT, WBT, CY1); 3: (E, E1, AP, BP, T, DCZ1, ITER)
    tooFew: string         ## FORMAT 940
    tooMany: string        ## FORMAT 944

const
  fmt836 = "(/,2X,'AVERAGE ENERGY OF PRIMARY ELECTRON =',F10.3,' EV.',/," &
           "'   ENERGY SPREAD OF PRIMARY ELECTRON =',F10.3,' EV.')"
  fmt546t = "(2X,' PROGRAM STOPPED . NPONT=',I5,' ITER=',I10)"
  fmt546a = "(2X,' PROGRAM STOPPED . NPONT=',I4,' ITER=',I10)"
  fmt546s = "(2X,' PROGRAM STOPPED NPONT=',I5,' ITER=',I10)"
  fmt878w10 = "(/,' TOTAL NO OF ELECTRONS=',I10,/,' TOTAL NO OF NEG. IONS=',I10,/," &
              "' TOTAL NO OF PRIMARIES=',I10)"
  fmt878w8 = "(/,' TOTAL NO OF ELECTRONS=',I8,/,' TOTAL NO OF NEG. IONS=',I8,/," &
             "' TOTAL NO OF PRIMARIES=',I8)"
  fmt913 = "(3X,' AFTER STORE ITER=',I10,' E1=',D12.3,' T=',D12.3,' AP=',D12.3,' BP=',D12.3," &
           "' DCZ1=',D12.3)"
  fmt913at = "(3X,' AFTER STORE ITER=',I10,' E=',D12.3,' E1=',D12.3,' AP=',D12.3,' BP='," &
             "D12.3,' DCZ1=',D12.3)"
  fmt913g = "(3X,' AFTER STORE ITER=',I10,' DZ=',D12.3,'E1=',D12.3,' COSWT=',D12.3,' SINWT='," &
            "D12.3,' WBT=',D12.3,' CY1=',D12.3)"
  fmt999 = "(2X,' WARNING ENERGY LT.0. E=',D12.3,' E1=',D12.3,' AP=',D12.3,' BP=',D12.3," &
           "' T=',D12.3,/,' DCZ1=',D12.3,' ITER=',I10)"
  fmt940t = "(2(/)," &
            "' PROGRAM STOPPED  (TOO FEW COLLISIONS), DECREASE THE ESTIMATED ALPHA.  NCLUS = '," &
            "I7,'  ITER =',I9)"
  fmt940s = "(2(/)," &
            "' PROGRAM STOPPED (TOO FEW COLLISIONS), DECREASE THE  ESTIMATED ALPHA.  NCLUS='," &
            "I7,'  ITER =',I9,' NELEC=',I9,' NEION =',I6)"
  fmt944t = "(/,2X,'PROGRAM STOPPED TOO MANY PRIMARIES IPRIM =',I7)"
  fmt944s = "(2X,' PROGRAM STOPPED TOO MANY PRIMARIES IPRIM=',I7)"

  formatsT = Formats(stackFull: fmt546t, totals: fmt878w10, negative: fmt913,
                     negativeKind: 1, tooFew: fmt940t, tooMany: fmt944t)
  formatsA = Formats(stackFull: fmt546a, totals: fmt878w8, negative: fmt913,
                     negativeKind: 1, tooFew: fmt940t, tooMany: fmt944t)
  formatsAT = Formats(stackFull: fmt546t, totals: fmt878w8, negative: fmt913at,
                      negativeKind: 1, tooFew: fmt940t, tooMany: fmt944t)
  formatsG = Formats(stackFull: fmt546t, totals: fmt878w8, negative: fmt913g,
                     negativeKind: 2, tooFew: fmt940t, tooMany: fmt944t)
  formatsD = Formats(stackFull: fmt546s, totals: fmt878w10, negative: fmt999,
                     negativeKind: 3, tooFew: fmt940s, tooMany: fmt944s)
  fmt11 = "(/,2X,'SOLUTION FOR STEADY STATE TOWNSEND PARAMETERS',/," &
          "'  -------------------------------------------------')"
  fmt12 = "(1(/),' SPACE STEP BETWEEN SAMPLING PLANES =',D12.5,' MICRONS.',/)"
  fmt18 = "(/,' SST DRIFT VELOCITIES')"
  fmt19 = "(/,' VD=',F9.2,' +- ',F6.2,' %   WS=',F9.2,' +- ',F6.2,' %')"
  fmt19t = "(/,' VD=',F9.1,' +- ',F6.2,' %   WS=',F9.1,' +- ',F6.2,' %')"
  fmt20 = "(/,' SST DIFFUSION')"
  fmt21 = "(/,' DL=',F9.2,' +- ',F6.1,' %   DT=',F9.2,' +- ',F6.2,' %')"
  fmt21t = "(/,' DL=',F9.1,' +- ',F6.1,' %   DT=',F9.1,' +- ',F6.2,' %')"
  fmt22 = "(/,' SST TOWNSEND COEFICIENTS')"
  fmt23 = "(/,' ALPHA=',F9.2' +- ',F6.2,' %    ATT=',F9.2,' +- ',F6.2, ' %')"
  fmt23t = "(/,' ALPHA =',F9.1,' +- ',F6.2,' %    ATT=',F9.1,' +- ',F6.2,' %')"
  fmt25 = "(/,2X,'SOLUTION FOR PULSED TOWNSEND AND TIME OF FLIGHT PARAMETERS',/," &
          "'  --------------------------------------------------------------')"
  fmt26 = "(1(/),'  TIME STEP BETWEEN SAMPLING PLANES =',D12.5,' PICOSECS.',/)"
  fmt27e = "(/,' PT IONISATION AND ATTACHMENT RATES *10**12/SEC',/,' ALPHA=',D10.3,' +- '," &
           "F6.2,' %      ATT=',D10.3,' +- ',F6.2,' %')"
  fmt27 = "(/,'PT IONISATION AND ATTACHMENT RATES *10**12/SEC',/,'ALPHA=',D10.3,' +- ',F6.2," &
          "' %      ATT=',D10.3,' +- ',F6.2,' %')"
  fmt28 = "(/,' TOF DIFFUSION')"
  fmt28b = "(/,'TOF DIFFUSION')"
  fmt29 = "(/,'  DL=',F8.1,' +- ',F6.1,' %       DT=',F8.1,' +- ',F6.1,' %')"
  fmt29b = "(/,'DZZ=',F8.1,' +- ',F5.1,' %',/,'DXX=',F8.1,' +- ',F5.1,' %',/,'DYY=',F8.1," &
           "' +- ',F5.1,' %',/,'DYZ=',F8.1,' +- ',F5.1,' %')"
  fmt29c = "(/,'DZZ=',F8.1,' +- ',F5.1,' %',/,'DXX=',F8.1,' +- ',F5.1,' %',/,'DYY=',F8.1," &
           "' +- ',F5.1,' %',/,'DYZ=',F8.1,' +- ',F5.1,' %',/,'DXZ=',F8.1,' +- ',F5.1,' %',/," &
           "'DXY=',F8.1,' +- ',F5.1,' %')"
  fmt30 = "(/,' TOF DRIFT VELOCITY')"
  fmt30b = "(/,'TOF DRIFT VELOCITY')"
  fmt31e = "(/,'    WR=',F8.2,' +-',F6.1,' %')"
  fmt31a = "(/,'    WR=',F8.2,' +- ',F6.2,' %')"
  fmt31b = "(/,'WRZ=',F8.2,' +-',F6.1,' %    WRY=',F8.2,' +-',F6.1,' %')"
  fmt31c = "(/,'WRZ=',F8.2,' +-',F6.1,' %    WRY=',F8.2,' +-',F6.1,' %   WRX=',F8.2,' +-'," &
           "F6.1,' %')"
  fmt888 = "(/,'TOWNSEND COEFICIENT (ALPHA-ATT) CALCULATED FROM TOF RESULTS:',/," &
           "' ALPHA-ATT /CM.=',D11.4,' N.B. APPROXIMATE FORMULA NOT ACCURATE',/)"
  # SST
  fmt800s = "(2(/),' STEADY STATE TOWNSEND RESULTS FOR',I2,' SEQUENTIAL SPACE PLANES',2(/)," &
            "'PLANE   NEL     VD      WS     DL       DT     EBAR   ALFN    ALFNJ   ALFNE',/)"
  fmt810s = "(1X,I2,2X,I7,2(1X,F6.1),2F9.1,F6.2,3F8.2)"
  # PT, PTG, PTH
  fmt900p = "(2(/),' PULSED TOWNSEND RESULTS AT',I2,' SEQUENTIAL TIME PLANES',/," &
            "' PLANE NO.   (ION-ATT) FREQ.   ENERGY      WV      NO.OF ELECTRONS',/)"
  fmt910p = "(2X,I2,4X,E15.4,7X,F7.2,4X,F7.2,3X,I8)"
  fmt900pg = "(2(/),' PULSED TOWNSEND RESULTS AT',I2,' SEQUENTIAL TIME PLANES',/," &
             "'PLANE   (ION-ATT)FRQ.    ENERGY      WVZ       WVY   NO.OF ELECTRONS',/)"
  fmt910pg = "(1X,I2,4X,D12.4,4X,F7.2,4X,F6.1,4X,F6.1,4X,I8)"
  fmt900ph = "(2(/),' PULSED TOWNSEND RESULTS AT',I2,' SEQUENTIAL TIME PLANES',/," &
             "'PLANE   (ION-ATT)FRQ.    ENERGY      WVZ       WVY       WVX  NO.OF ELECTRNS'," &
             "/)"
  fmt910ph = "(1X,I2,4X,D12.4,4X,F7.2,4X,F6.1,4X,F6.1,4X,F6.1,4X,I8)"
  # TOF, TOFG, TOFH
  fmt900t = "(2(/),' TIME OF FLIGHT RESULTS AT',I2,' SEQUENTIAL TIME PLANES',/," &
            "' PLANE NO.        DL         DX          DY            WR',/)"
  fmt910t = "(3X,I3,4X,3F12.1,4X,F8.2)"
  fmt900tg = "(2(/),' TIME OF FLIGHT RESULTS AT',I2,' SEQUENTIAL TIME PLANES',/," &
             "'PLANE    DZZ      DXX      DYY      DYZ        WRZ         WRY',/)"
  fmt910tg = "(1X,I2,2X,4F9.1,4X,F8.2,4X,F8.2)"
  fmt900th = "(2(/),' TIME OF FLIGHT RESULTS AT',I2,' SEQUENTIAL TIME PLANES',/," &
             "'PLANE    DZZ     DXX     DYY     DYZ     DXZ     DXY      WRZ     WRY     WRX'," &
             "/)"
  fmt910th = "(1X,I2,2X,6F8.1,3F8.2)"
  # FRIEDLAND
  fmt900f = "(2(/),' ESTIMATE USING FRIEDLAND :',/,' AVERAGE ENERGY =',F8.3,'EV.',/," &
            "' AVERAGE IONISATION =',E11.4,' *10**12/SEC',/,' AVERAGE ATTACHMENT =',E11.4," &
            "' *10**12/SEC')"

proc printNegative(o: File, fm: Formats, f: Field, w: NegativeEnergy) =
  let ap = w.el.dir.z * f.f2 * sqrt(w.el.e)
  case fm.negativeKind
  of 1: o.fwrite(fm.negative, w.iter, w.fl.e, w.el.e, ap, f.bp,
               w.el.dir.z)
  of 2: o.fwrite(fm.negative, w.iter, w.fl.dz, w.fl.e, w.fl.cwt, w.fl.swt,
               float64(f.wb * w.fl.t), w.el.vel.y)
  else: o.fwrite(fm.negative, w.fl.e, w.el.e, ap, f.bp, w.fl.t,
               w.el.dir.z, w.iter)

type
  AvalancheRun[P: static PlaneKind] = object
    ## what the analyses need from a run
    when P == pkTime:
      sums: TimeSums
    else:
      sums: SpaceSums
    iprim: int
    ratio: AttachmentRatio

func formatsOf(G: FieldGeometry, P: PlaneKind, thermal: bool): Formats =
  if P == pkSpace: formatsD
  else:
    case G
    of fgElectric: formatsT
    of fgParallel: (if thermal: formatsAT else: formatsA)
    of fgCrossed, fgOblique: formatsG

type
  Walker[G: static FieldGeometry; P: static PlaneKind; T, Tabs, Rng] = object
    ## an avalanche run on a thread of its own
    c: Cascade[G, P, T]
    rng: Rng
    tabs: ptr Tabs

proc walk[G: static FieldGeometry; P: static PlaneKind; T, Tabs, Rng](
    w: ptr Walker[G, P, T, Tabs, Rng]) {.thread.} =
  w.c.run(w.tabs[], w.rng)

proc avalanche[G: static FieldGeometry; P: static PlaneKind; M](
    r: var Run[M], planes: Planes, nmax: int, rng: var MagRng): AvalancheRun[P] =
  ## one run of avalanches sampled at `planes`; leaves its collisions in
  ## `r.last`. With `r.threads` > 1, that many walkers share NMAX (and the
  ## limit on primaries), each with random numbers of its own seeded from
  ## `rng`, and their counts and sums are added up: the same estimates,
  ## not the same numbers. Unless `r.fortran`, the avalanches take the
  ## exact bounds per window, and random numbers of `r.generator`'s kind
  ## (also on one thread).
  let fm = formatsOf(G, P, thermalMotion(M))
  let f = r.field(G)
  let start = initElectron[G](f, r.estart, r.theta, startPhi)
  let d = r.deck
  when thermalMotion(M):
    type Target = WarmGas
    template tabs: untyped = r.mix.tables
    template cascade(n, limit: int, g: untyped): Cascade[G, P, Target] =
      initCascade[G, P, WarmGas](f, const2, const3, start, initWarmGas(r.mix.gas, g),
                                 r.estep, n, d.penning, planes, d.btheta, limit, windowBounds)
  else:
    type Target = ColdBound
    template tabs: untyped = r.mix.table
    var tlim = r.mix.bounds[1]
    for i in 2 .. 8:
      if tlim < r.mix.bounds[i]: tlim = r.mix.bounds[i]
    template cascade(n, limit: int, g: untyped): Cascade[G, P, Target] =
      initCascade[G, P, ColdBound](f, const2, const3, start, ColdBound(tlim: tlim),
                                   r.estep, n, d.penning, planes, d.btheta, limit, windowBounds)
  # exact bounds per window, unless Magboltz's numbers are wanted
  let maxima = (if r.fortran: nil else: flightMaxima(tabs))
  let windowBounds = (if maxima == nil: nil else: addr maxima[])
  var c: Cascade[G, P, Target]
  template split(Rng: typedesc, n: int, seed: untyped) =
    ## n walkers, each with random numbers `seed`(`rng`) of type Rng
    var walkers = newSeq[Walker[G, P, Target, typeof(tabs), Rng]](n)
    for k, w in walkers.mpairs:
      w.rng = seed
      w.c = cascade(nmax div n + (if k < nmax mod n: 1 else: 0),
                    maxPrimaries div n, w.rng)
      w.tabs = addr tabs
    if n == 1:
      walk(addr walkers[0])
    else:
      when compileOption("threads"):
        var threads = newSeq[Thread[ptr Walker[G, P, Target, typeof(tabs), Rng]]](n)
        for k in 0 ..< n:
          createThread(threads[k], walk[G, P, Target, typeof(tabs), Rng], addr walkers[k])
        joinThreads(threads)
    c = move walkers[0].c
    for k in 1 ..< n: c.merge(walkers[k].c)
  let n = (when compileOption("threads"): max(r.threads, 1) else: 1)
  if r.fortran and n == 1:
    # Magboltz's one stream of random numbers
    c = cascade(nmax, maxPrimaries, rng)
    c.run(tabs, rng)
  elif r.fortran:
    proc seeded(rng: var MagRng): MagRng =
      result = initMagRng(nil)
      result.rm.rm48in(int(drand48(rng) * 900000000.0), 0, 0)
    split(MagRng, n, seeded(rng))
  else:
    withRng(r.generator, R):
      split(R, n, initRng(R, int(drand48(rng) * 900000000.0)))
  r.last = Collisions[M](tally: c.tally, st: c.st, ttots: c.ttots, spectrum: c.spectrum)
  for w in c.warnings:
    r.output.printNegative(fm, f, w)
    r.warnings.add "negative energy after a collision (ITER " & $w.iter & ")"
  case c.stop
  of stStackFull:
    r.output.fwrite(fm.stackFull, c.stacked + 1, c.id)
    r.stop("avalanche: too many electrons stacked (NPONT " & $(c.stacked + 1) & ")")
  of stTooManyPrimaries:
    r.output.fwrite(fm.tooMany, c.iprim)
    const msg = "avalanche: too many primaries (IPRIM "
    when P == pkSpace: r.stop(msg & $c.iprim & ")")
    else: r.warnings.add msg & $c.iprim & ")"
  of stNone: discard
  r.theta = c.lastTheta
  result.sums = c.sums
  result.iprim = c.iprim
  result.ratio = attachmentRatio(c.neion, c.nelec, c.iprim)
  if c.j1 == 1:
    when P == pkTime:
      r.output.fwrite(fm.tooFew, c.nclus, c.id)
      r.stop("avalanche: too few collisions, the estimated α is too large (NCLUS " &
             $c.nclus & ")")
    else:
      r.output.fwrite(fm.tooFew, c.nclus, c.id, c.nelec, c.neion)
      r.warnings.add "avalanche: too few collisions, the estimated α is too large (NCLUS " &
                     $c.nclus & ")"
  r.output.fwrite(fm.totals, c.nelec, c.neion, c.iprim)
  if c.iprim != 1:
    let ebar = c.primaryE / float64(c.iprim)
    let eerr = sqrt(c.primaryE2 / float64(c.iprim) - ebar*ebar)
    r.output.fwrite(fmt836, ebar, eerr)

proc printSst(o: File, s: SteadyState) =
  o.fwrite(fmt800s, s.planes)
  for k in 1 .. s.planes:
    let p = s.plane[k]
    o.fwrite(fmt810s, k, p.n, p.vd, p.ws, p.dl, p.dt, p.e,
           p.alfn, p.alfnj, p.alfne)

proc printPt[G: static FieldGeometry](o: File, pt: PulsedTownsend, s: TimeSums) =
  when G in {fgElectric, fgParallel}:
    o.fwrite(fmt900p, pt.planes)
    for k in 1 .. pt.planes:
      let p = pt.plane[k]
      o.fwrite(fmt910p, k, p.rate, p.e, p.v.z, s.n[k])
  elif G == fgCrossed:
    o.fwrite(fmt900pg, pt.planes)
    for k in 1 .. pt.planes:
      let p = pt.plane[k]
      o.fwrite(fmt910pg, k, p.rate, p.e, p.v.z, p.v.y, s.n[k])
  else:
    o.fwrite(fmt900ph, pt.planes)
    for k in 1 .. pt.planes:
      let p = pt.plane[k]
      o.fwrite(fmt910ph, k, p.rate, p.e, p.v.z, p.v.y, p.v.x,
             s.n[k])

proc printTof[G: static FieldGeometry](o: File, tof: TimeOfFlight, planes: int) =
  when G in {fgElectric, fgParallel}:
    o.fwrite(fmt900t, planes)
    for k in 1 .. planes:
      let p = tof.plane[k]
      o.fwrite(fmt910t, k, p.d.zz, p.d.xx, p.d.yy, p.w.z)
  elif G == fgCrossed:
    o.fwrite(fmt900tg, planes)
    for k in 1 .. planes:
      let p = tof.plane[k]
      o.fwrite(fmt910tg, k, p.d.zz, p.d.xx, p.d.yy, p.d.yz,
             p.w.z, p.w.y)
  else:
    o.fwrite(fmt900th, planes)
    for k in 1 .. planes:
      let p = tof.plane[k]
      o.fwrite(fmt910th, k, p.d.zz, p.d.xx, p.d.yy, p.d.yz,
             p.d.xz, p.d.xy, p.w.z, p.w.y, p.w.x)

proc printFriedland[M](r: var Run[M]) =
  ## FRIEDLAND, FRIEDLANDT: from the spectrum of the last run
  when thermalMotion(M):
    let tcf = r.mix[].totalTcf(r.deck.ngas)
  else:
    template tcf: untyped = r.mix.table.tcf
  let fr = friedland(r.last.spectrum, r.grid.e, tcf, r.fcion, r.fcatt, r.last.ttots)
  r.friedland = some(fr)
  r.output.fwrite(fmt900f, fr.e, fr.ionisation, fr.attachment)

proc analyse[G: static FieldGeometry; M](r: var Run[M], rng: var MagRng) =
  const thermal = thermalMotion(M)
  let nmax = max(r.cond.nmax div 10000000, 5) * 10000000
  let (planes, _, _) = samplingPlanes(r.transport.alpha, r.transport.att,
                                      r.transport.vel.z, r.deck.tempc, r.deck.torr)
  when G == fgElectric:
    # steady-state Townsend
    r.output.fwrite(fmt11)
    r.output.fwrite(fmt12, planes.zstep.to(μm))
    let sr = avalanche[G, pkSpace, M](r, planes, nmax, rng)
    let sst = steadyState(sr.sums, sr.iprim, planes.zstep, planes.izfinal, sr.ratio)
    r.steadyState = some(sst)
    r.output.printSst(sst)
    r.transport.alpha = sst.alpha
    r.transport.alper = sst.alphaErr
    r.transport.att = sst.att
    r.transport.atter = sst.attErr
    r.output.fwrite(fmt18)
    r.output.fwrite(when thermal: fmt19t else: fmt19, sst.vd, sst.vdErr, sst.ws, sst.wsErr)
    r.output.fwrite(fmt20)
    r.output.fwrite(when thermal: fmt21t else: fmt21, sst.dl, sst.dlErr, sst.dt, sst.dtErr)
    r.output.fwrite(fmt22)
    r.output.fwrite(when thermal: fmt23t else: fmt23, sst.alpha, sst.alphaErr, sst.att,
           sst.attErr)
    when not thermal: r.printCollisionCounts()  # ALPCALCT has CALL OUTPUT1T commented out
  # pulsed Townsend and time of flight
  r.output.fwrite(fmt25)
  r.output.fwrite(fmt26, planes.tstep)
  let tr = avalanche[G, pkTime, M](r, planes, nmax, rng)
  r.printFriedland()
  let pt = pulsedTownsend[G](tr.sums, tr.iprim, planes.tstep, planes.itfinal)
  r.pulsedTownsend = some(pt)
  printPt[G](r.output, pt, tr.sums)
  var ratio = tr.ratio
  when G == fgCrossed:
    # TOFG: the ratio from the collision frequencies of the run
    let fr = frequencies(r.last, r.deck.ngas)
    if fr.ionisation > 0.0.ps⁻¹: ratio.ratio = fr.attachment / fr.ionisation  # MB-40
  let tof = timeOfFlight[G](tr.sums, pt, planes.tstep, ratio)
  r.timeOfFlight = some(tof)
  printTof[G](r.output, tof, pt.planes)
  r.output.fwrite(when G == fgElectric: fmt27e else: fmt27, tof.alpha, tof.alphaErr,
         tof.att, tof.attErr)
  when G in {fgElectric, fgParallel}:
    r.output.fwrite(fmt28)
    r.output.fwrite(fmt29, tof.dl, tof.dlErr, tof.dt, tof.dtErr)
    r.output.fwrite(fmt30)
    r.output.fwrite(when G == fgElectric: fmt31e else: fmt31a, tof.w.z, tof.wErr.z)
    # α - η from the drift and growth of a Gaussian swarm
    let wrzn = tof.w.z.to(cm•s⁻¹)
    let fc1 = wrzn / (2.0 * tof.dl)
    let fc2 = (tof.alpha - tof.att).to(s⁻¹) / tof.dl
    r.output.fwrite(fmt888, float64(fc1 - sqrt(fc1*fc1 - fc2)))
    when G == fgParallel:
      r.transport.alpha = (tof.alpha / tof.w.z).to(cm⁻¹)
      r.transport.att = (tof.att / tof.w.z).to(cm⁻¹)
      # NIMBOLTZ FIX MB-42: errors in percent (was errAlpha·α/100)
      r.transport.alper = tof.alphaErr
      r.transport.atter = tof.attErr
  else:
    let d = tof.d
    let de = tof.dErr
    r.output.fwrite(fmt28b)
    when G == fgCrossed:
      r.output.fwrite(fmt29b, d.zz, de.zz, d.xx, de.xx, d.yy, de.yy, d.yz,
             de.yz)
    else:
      r.output.fwrite(fmt29c, d.zz, de.zz, d.xx, de.xx, d.yy, de.yy, d.yz,
             de.yz, d.xz, de.xz, d.xy, de.xy)
    r.output.fwrite(fmt30b)
    when G == fgCrossed:
      r.output.fwrite(fmt31b, tof.w.z, tof.wErr.z, tof.w.y, tof.wErr.y)
      let tofwr = sqrt(tof.w.z * tof.w.z + tof.w.y * tof.w.y)
    else:
      r.output.fwrite(fmt31c, tof.w.z, tof.wErr.z, tof.w.y, tof.wErr.y, tof.w.x,
             tof.wErr.x)
      let tofwr = sqrt(tof.w.z * tof.w.z + tof.w.y * tof.w.y + tof.w.x * tof.w.x)
    # as ALPCLCB/C: the drift velocity along E; ALPCLCBT/CT: the full one
    let wrzn = (when thermal: tofwr else: tof.w.z).to(cm•s⁻¹)
    let fc1 = wrzn / (2.0 * d.zz)
    let fc2 = (tof.alpha - tof.att).to(s⁻¹) / d.zz
    r.output.fwrite(fmt888, float64(fc1 - sqrt(fc1*fc1 - fc2)))
    r.transport.alpha = (tof.alpha / tofwr).to(cm⁻¹)
    r.transport.att = (tof.att / tofwr).to(cm⁻¹)
    # NIMBOLTZ FIX MB-42: errors in percent (was errAlpha·α/100; ALPCLCB
    # stored that in locals, ALPCLCBT multiplied the η error by 100)
    r.transport.alper = tof.alphaErr
    r.transport.atter = tof.attErr

func needsAnalyses*(r: Run): bool =
  ## MAGBOLTZ runs the analyses when the ionisation or attachment
  ## coefficient (normalised to 760 Torr, 20 °C) exceeds 5 /cm and the
  ## two do not cancel to within 5 /cm
  const sstmin = 5.0.cm⁻¹
  let tgas = 273.15 + r.deck.tempc
  let alpp = r.transport.alpha * 760.0.Torr * tgas / (r.deck.torr * 293.15)
  let attp = r.transport.att * 760.0.Torr * tgas / (r.deck.torr * 293.15)
  if abs(alpp - attp) > sstmin: true
  elif alpp > sstmin or attp > sstmin: not (abs(alpp - attp) < sstmin)
  else: false

proc analyse*[M](r: var Run[M], rng: var MagRng) =
  withGeometry(monteGeometry(r.deck), G):
    analyse[G, M](r, rng)
