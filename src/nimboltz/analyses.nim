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

import std/math
import frt, units, rm48, tables, vec3, fields, freeflight, cascade, townsend
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
  fmt836 = "(/,2X,\'AVERAGE ENERGY OF PRIMARY ELECTRON =\',F10.3,\' EV.\',/,\'   ENERGY SPREAD OF PRIMARY ELECTRON =\',F10.3,\' EV.\')"
  fmt546t = "(2X,\' PROGRAM STOPPED . NPONT=\',I5,\' ITER=\',I10)"
  fmt546a = "(2X,\' PROGRAM STOPPED . NPONT=\',I4,\' ITER=\',I10)"
  fmt546s = "(2X,\' PROGRAM STOPPED NPONT=\',I5,\' ITER=\',I10)"
  fmt878w10 = "(/,\' TOTAL NO OF ELECTRONS=\',I10,/,\' TOTAL NO OF NEG. IONS=\',I10,/,\' TOTAL NO OF PRIMARIES=\',I10)"
  fmt878w8 = "(/,\' TOTAL NO OF ELECTRONS=\',I8,/,\' TOTAL NO OF NEG. IONS=\',I8,/,\' TOTAL NO OF PRIMARIES=\',I8)"
  fmt913 = "(3X,\' AFTER STORE ITER=\',I10,\' E1=\',D12.3,\' T=\',D12.3,\' AP=\',D12.3,\' BP=\',D12.3,\' DCZ1=\',D12.3)"
  fmt913at = "(3X,\' AFTER STORE ITER=\',I10,\' E=\',D12.3,\' E1=\',D12.3,\' AP=\',D12.3,\' BP=\',D12.3,\' DCZ1=\',D12.3)"
  fmt913g = "(3X,\' AFTER STORE ITER=\',I10,\' DZ=\',D12.3,\'E1=\',D12.3,\' COSWT=\',D12.3,\' SINWT=\',D12.3,\' WBT=\',D12.3,\' CY1=\',D12.3)"
  fmt999 = "(2X,\' WARNING ENERGY LT.0. E=\',D12.3,\' E1=\',D12.3,\' AP=\',D12.3,\' BP=\',D12.3,\' T=\',D12.3,/,\' DCZ1=\',D12.3,\' ITER=\',I10)"
  fmt940t = "(2(/),\' PROGRAM STOPPED  (TOO FEW COLLISIONS), DECREASE THE ESTIMATED ALPHA.  NCLUS = \',I7,\'  ITER =\',I9)"
  fmt940s = "(2(/),\' PROGRAM STOPPED (TOO FEW COLLISIONS), DECREASE THE  ESTIMATED ALPHA.  NCLUS=\',I7,\'  ITER =\',I9,\' NELEC=\',I9,\' NEION =\',I6)"
  fmt944t = "(/,2X,\'PROGRAM STOPPED TOO MANY PRIMARIES IPRIM =\',I7)"
  fmt944s = "(2X,\' PROGRAM STOPPED TOO MANY PRIMARIES IPRIM=\',I7)"

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
  fmt11 = "(/,2X,\'SOLUTION FOR STEADY STATE TOWNSEND PARAMETERS\',/,\'  -------------------------------------------------\')"
  fmt12 = "(1(/),\' SPACE STEP BETWEEN SAMPLING PLANES =\',D12.5,\' MICRONS.\',/)"
  fmt18 = "(/,\' SST DRIFT VELOCITIES\')"
  fmt19 = "(/,\' VD=\',F9.2,\' +- \',F6.2,\' %   WS=\',F9.2,\' +- \',F6.2,\' %\')"
  fmt19t = "(/,\' VD=\',F9.1,\' +- \',F6.2,\' %   WS=\',F9.1,\' +- \',F6.2,\' %\')"
  fmt20 = "(/,\' SST DIFFUSION\')"
  fmt21 = "(/,\' DL=\',F9.2,\' +- \',F6.1,\' %   DT=\',F9.2,\' +- \',F6.2,\' %\')"
  fmt21t = "(/,\' DL=\',F9.1,\' +- \',F6.1,\' %   DT=\',F9.1,\' +- \',F6.2,\' %\')"
  fmt22 = "(/,\' SST TOWNSEND COEFICIENTS\')"
  fmt23 = "(/,\' ALPHA=\',F9.2\' +- \',F6.2,\' %    ATT=\',F9.2,\' +- \',F6.2, \' %\')"
  fmt23t = "(/,\' ALPHA =\',F9.1,\' +- \',F6.2,\' %    ATT=\',F9.1,\' +- \',F6.2,\' %\')"
  fmt25 = "(/,2X,\'SOLUTION FOR PULSED TOWNSEND AND TIME OF FLIGHT PARAMETERS\',/,\'  --------------------------------------------------------------\')"
  fmt26 = "(1(/),\'  TIME STEP BETWEEN SAMPLING PLANES =\',D12.5,\' PICOSECS.\',/)"
  fmt27e = "(/,\' PT IONISATION AND ATTACHMENT RATES *10**12/SEC\',/,\' ALPHA=\',D10.3,\' +- \',F6.2,\' %      ATT=\',D10.3,\' +- \',F6.2,\' %\')"
  fmt27 = "(/,\'PT IONISATION AND ATTACHMENT RATES *10**12/SEC\',/,\'ALPHA=\',D10.3,\' +- \',F6.2,\' %      ATT=\',D10.3,\' +- \',F6.2,\' %\')"
  fmt28 = "(/,\' TOF DIFFUSION\')"
  fmt28b = "(/,\'TOF DIFFUSION\')"
  fmt29 = "(/,\'  DL=\',F8.1,\' +- \',F6.1,\' %       DT=\',F8.1,\' +- \',F6.1,\' %\')"
  fmt29b = "(/,\'DZZ=\',F8.1,\' +- \',F5.1,\' %\',/,\'DXX=\',F8.1,\' +- \',F5.1,\' %\',/,\'DYY=\',F8.1,\' +- \',F5.1,\' %\',/,\'DYZ=\',F8.1,\' +- \',F5.1,\' %\')"
  fmt29c = "(/,\'DZZ=\',F8.1,\' +- \',F5.1,\' %\',/,\'DXX=\',F8.1,\' +- \',F5.1,\' %\',/,\'DYY=\',F8.1,\' +- \',F5.1,\' %\',/,\'DYZ=\',F8.1,\' +- \',F5.1,\' %\',/,\'DXZ=\',F8.1,\' +- \',F5.1,\' %\',/,\'DXY=\',F8.1,\' +- \',F5.1,\' %\')"
  fmt30 = "(/,\' TOF DRIFT VELOCITY\')"
  fmt30b = "(/,\'TOF DRIFT VELOCITY\')"
  fmt31e = "(/,\'    WR=\',F8.2,\' +-\',F6.1,\' %\')"
  fmt31a = "(/,\'    WR=\',F8.2,\' +- \',F6.2,\' %\')"
  fmt31b = "(/,\'WRZ=\',F8.2,\' +-\',F6.1,\' %    WRY=\',F8.2,\' +-\',F6.1,\' %\')"
  fmt31c = "(/,\'WRZ=\',F8.2,\' +-\',F6.1,\' %    WRY=\',F8.2,\' +-\',F6.1,\' %   WRX=\',F8.2,\' +-\',F6.1,\' %\')"
  fmt888 = "(/,\'TOWNSEND COEFICIENT (ALPHA-ATT) CALCULATED FROM TOF RESULTS:\',/,\' ALPHA-ATT /CM.=\',D11.4,\' N.B. APPROXIMATE FORMULA NOT ACCURATE\',/)"
  # SST
  fmt800s = "(2(/),\' STEADY STATE TOWNSEND RESULTS FOR\',I2,\' SEQUENTIAL SPACE PLANES\',2(/),\'PLANE   NEL     VD      WS     DL       DT     EBAR   ALFN    ALFNJ   ALFNE\',/)"
  fmt810s = "(1X,I2,2X,I7,2(1X,F6.1),2F9.1,F6.2,3F8.2)"
  # PT, PTG, PTH
  fmt900p = "(2(/),\' PULSED TOWNSEND RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\' PLANE NO.   (ION-ATT) FREQ.   ENERGY      WV      NO.OF ELECTRONS\',/)"
  fmt910p = "(2X,I2,4X,E15.4,7X,F7.2,4X,F7.2,3X,I8)"
  fmt900pg = "(2(/),\' PULSED TOWNSEND RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\'PLANE   (ION-ATT)FRQ.    ENERGY      WVZ       WVY   NO.OF ELECTRONS\',/)"
  fmt910pg = "(1X,I2,4X,D12.4,4X,F7.2,4X,F6.1,4X,F6.1,4X,I8)"
  fmt900ph = "(2(/),\' PULSED TOWNSEND RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\'PLANE   (ION-ATT)FRQ.    ENERGY      WVZ       WVY       WVX  NO.OF ELECTRNS\',/)"
  fmt910ph = "(1X,I2,4X,D12.4,4X,F7.2,4X,F6.1,4X,F6.1,4X,F6.1,4X,I8)"
  # TOF, TOFG, TOFH
  fmt900t = "(2(/),\' TIME OF FLIGHT RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\' PLANE NO.        DL         DX          DY            WR\',/)"
  fmt910t = "(3X,I3,4X,3F12.1,4X,F8.2)"
  fmt900tg = "(2(/),\' TIME OF FLIGHT RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\'PLANE    DZZ      DXX      DYY      DYZ        WRZ         WRY\',/)"
  fmt910tg = "(1X,I2,2X,4F9.1,4X,F8.2,4X,F8.2)"
  fmt900th = "(2(/),\' TIME OF FLIGHT RESULTS AT\',I2,\' SEQUENTIAL TIME PLANES\',/,\'PLANE    DZZ     DXX     DYY     DYZ     DXZ     DXY      WRZ     WRY     WRX\',/)"
  fmt910th = "(1X,I2,2X,6F8.1,3F8.2)"
  # FRIEDLAND
  fmt900f = "(2(/),\' ESTIMATE USING FRIEDLAND :\',/,\' AVERAGE ENERGY =\',F8.3,\'EV.\',/,\' AVERAGE IONISATION =\',E11.4,\' *10**12/SEC\',/,\' AVERAGE ATTACHMENT =\',E11.4,\' *10**12/SEC\')"

proc printNegative(fm: Formats, f: Field, w: NegativeEnergy) =
  let ap = w.el.dir.z * f.f2 * sqrt(w.el.e)
  case fm.negativeKind
  of 1: fwrite(fm.negative, w.iter, w.fl.e.float, w.el.e.float, ap.float, f.bp.float,
               w.el.dir.z)
  of 2: fwrite(fm.negative, w.iter, w.fl.dz.float, w.fl.e.float, w.fl.cwt, w.fl.swt,
               float64(f.wb * w.fl.t), w.el.vel.y.float)
  else: fwrite(fm.negative, w.fl.e.float, w.el.e.float, ap.float, f.bp.float, w.fl.t.float,
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
  Walker[G: static FieldGeometry; P: static PlaneKind; T, Tabs] = object
    ## an avalanche run on a thread of its own
    c: Cascade[G, P, T]
    rng: MagRng
    tabs: ptr Tabs

proc walk[G: static FieldGeometry; P: static PlaneKind; T, Tabs](
    w: ptr Walker[G, P, T, Tabs]) {.thread.} =
  w.c.run(w.tabs[], w.rng)

proc avalanche[G: static FieldGeometry; P: static PlaneKind; M](
    r: var Run[M], planes: Planes, nmax: int, rng: var MagRng): AvalancheRun[P] =
  ## one run of avalanches sampled at `planes`; leaves its collisions in
  ## `r.last`. With `r.threads` > 1, that many walkers share NMAX (and the
  ## limit on primaries), each with random numbers of its own seeded from
  ## `rng`, and their counts and sums are added up: the same estimates,
  ## not the same numbers.
  let fm = formatsOf(G, P, thermalMotion(M))
  let f = r.field(G)
  let start = initElectron[G](f, r.estart, r.theta, startPhi)
  let d = r.deck
  when thermalMotion(M):
    type Target = WarmGas
    template tabs: untyped = r.mix.tables
    template cascade(n, limit: int, g: var MagRng): Cascade[G, P, Target] =
      initCascade[G, P, WarmGas](f, const2, const3, start, initWarmGas(r.mix.gas, g),
                                 r.estep, n, d.penning, planes, d.btheta, limit)
  else:
    type Target = ColdBound
    template tabs: untyped = r.mix.table
    var tlim = r.mix.bounds[1]
    for i in 2 .. 8:
      if tlim < r.mix.bounds[i]: tlim = r.mix.bounds[i]
    template cascade(n, limit: int, g: var MagRng): Cascade[G, P, Target] =
      initCascade[G, P, ColdBound](f, const2, const3, start, ColdBound(tlim: tlim),
                                   r.estep, n, d.penning, planes, d.btheta, limit)
  var c: Cascade[G, P, Target]
  if r.threads <= 1:
    c = cascade(nmax, maxPrimaries, rng)
    c.run(tabs, rng)
  else:
    let n = r.threads
    var walkers = newSeq[Walker[G, P, Target, typeof(tabs)]](n)
    for k, w in walkers.mpairs:
      w.rng = initMagRng(nil)
      w.rng.rm.rm48in(int(drand48(rng) * 900000000.0), 0, 0)
      w.c = cascade(nmax div n + (if k < nmax mod n: 1 else: 0),
                    maxPrimaries div n, w.rng)
      w.tabs = addr tabs
    var threads = newSeq[Thread[ptr Walker[G, P, Target, typeof(tabs)]]](n)
    for k in 0 ..< n: createThread(threads[k], walk[G, P, Target, typeof(tabs)], addr walkers[k])
    joinThreads(threads)
    c = move walkers[0].c
    for k in 1 ..< n: c.merge(walkers[k].c)
  r.last = Collisions[M](tally: c.tally, st: c.st, ttots: c.ttots, spectrum: c.spectrum)
  for w in c.warnings: printNegative(fm, f, w)
  case c.stop
  of stStackFull:
    fwrite(fm.stackFull, c.stacked + 1, c.id)
    fstop()
  of stTooManyPrimaries:
    fwrite(fm.tooMany, c.iprim)
    when P == pkSpace: fstop()
  of stNone: discard
  r.theta = c.theta
  result.sums = c.sums
  result.iprim = c.iprim
  result.ratio = attachmentRatio(c.neion, c.nelec, c.iprim)
  if c.j1 == 1:
    when P == pkTime:
      fwrite(fm.tooFew, c.nclus, c.id)
      fstop()
    else:
      fwrite(fm.tooFew, c.nclus, c.id, c.nelec, c.neion)
  fwrite(fm.totals, c.nelec, c.neion, c.iprim)
  if c.iprim != 1:
    let ebar = c.primaryE / float64(c.iprim)
    let eerr = sqrt(c.primaryE2 / float64(c.iprim) - ebar*ebar)
    fwrite(fmt836, ebar.float, eerr.float)

proc printSst(s: SteadyState) =
  fwrite(fmt800s, s.planes)
  for k in 1 .. s.planes:
    let p = s.plane[k]
    fwrite(fmt810s, k, p.n, p.vd.float, p.ws.float, p.dl.float, p.dt.float, p.e.float,
           p.alfn.float, p.alfnj.float, p.alfne.float)

proc printPt[G: static FieldGeometry](pt: PulsedTownsend, s: TimeSums) =
  when G in {fgElectric, fgParallel}:
    fwrite(fmt900p, pt.planes)
    for k in 1 .. pt.planes:
      let p = pt.plane[k]
      fwrite(fmt910p, k, p.rate.float, p.e.float, p.v.z.float, s.n[k])
  elif G == fgCrossed:
    fwrite(fmt900pg, pt.planes)
    for k in 1 .. pt.planes:
      let p = pt.plane[k]
      fwrite(fmt910pg, k, p.rate.float, p.e.float, p.v.z.float, p.v.y.float, s.n[k])
  else:
    fwrite(fmt900ph, pt.planes)
    for k in 1 .. pt.planes:
      let p = pt.plane[k]
      fwrite(fmt910ph, k, p.rate.float, p.e.float, p.v.z.float, p.v.y.float, p.v.x.float,
             s.n[k])

proc printTof[G: static FieldGeometry](tof: TimeOfFlight, planes: int) =
  when G in {fgElectric, fgParallel}:
    fwrite(fmt900t, planes)
    for k in 1 .. planes:
      let p = tof.plane[k]
      fwrite(fmt910t, k, p.d.zz.float, p.d.xx.float, p.d.yy.float, p.w.z.float)
  elif G == fgCrossed:
    fwrite(fmt900tg, planes)
    for k in 1 .. planes:
      let p = tof.plane[k]
      fwrite(fmt910tg, k, p.d.zz.float, p.d.xx.float, p.d.yy.float, p.d.yz.float,
             p.w.z.float, p.w.y.float)
  else:
    fwrite(fmt900th, planes)
    for k in 1 .. planes:
      let p = tof.plane[k]
      fwrite(fmt910th, k, p.d.zz.float, p.d.xx.float, p.d.yy.float, p.d.yz.float,
             p.d.xz.float, p.d.xy.float, p.w.z.float, p.w.y.float, p.w.x.float)

proc printFriedland[M](r: Run[M]) =
  ## FRIEDLAND, FRIEDLANDT: from the spectrum of the last run
  when thermalMotion(M):
    let tcf = r.mix[].totalTcf(r.deck.ngas)
  else:
    template tcf: untyped = r.mix.table.tcf
  let fr = friedland(r.last.spectrum, r.grid.e, tcf, r.fcion, r.fcatt, r.last.ttots)
  fwrite(fmt900f, fr.e.float, fr.ionisation.float, fr.attachment.float)

proc analyse[G: static FieldGeometry; M](r: var Run[M], rng: var MagRng) =
  const thermal = thermalMotion(M)
  let nmax = max(r.cond.nmax div 10000000, 5) * 10000000
  let (planes, _, _) = samplingPlanes(r.transport.alpha, r.transport.att,
                                      r.transport.vel.z, r.deck.tempc, r.deck.torr)
  when G == fgElectric:
    # steady-state Townsend
    fwrite(fmt11)
    fwrite(fmt12, planes.zstep.to(μm).float)
    let sr = avalanche[G, pkSpace, M](r, planes, nmax, rng)
    let sst = steadyState(sr.sums, sr.iprim, planes.zstep, planes.izfinal, sr.ratio)
    printSst(sst)
    r.transport.alpha = sst.alpha
    r.transport.alper = sst.alphaErr
    r.transport.att = sst.att
    r.transport.atter = sst.attErr
    fwrite(fmt18)
    fwrite(when thermal: fmt19t else: fmt19, sst.vd.float, sst.vdErr, sst.ws.float, sst.wsErr)
    fwrite(fmt20)
    fwrite(when thermal: fmt21t else: fmt21, sst.dl.float, sst.dlErr, sst.dt.float, sst.dtErr)
    fwrite(fmt22)
    fwrite(when thermal: fmt23t else: fmt23, sst.alpha.float, sst.alphaErr, sst.att.float,
           sst.attErr)
    when not thermal: r.printCollisionCounts()  # ALPCALCT has CALL OUTPUT1T commented out
  # pulsed Townsend and time of flight
  fwrite(fmt25)
  fwrite(fmt26, planes.tstep.float)
  let tr = avalanche[G, pkTime, M](r, planes, nmax, rng)
  r.printFriedland()
  let pt = pulsedTownsend[G](tr.sums, tr.iprim, planes.tstep, planes.itfinal)
  printPt[G](pt, tr.sums)
  var ratio = tr.ratio
  when G == fgCrossed:
    # TOFG: the ratio from the collision frequencies of the run
    let fr = frequencies(r.last, r.deck.ngas)
    if fr.ionisation > 0.0.ps⁻¹: ratio.ratio = fr.attachment / fr.ionisation  # MB-40
  let tof = timeOfFlight[G](tr.sums, pt, planes.tstep, ratio)
  printTof[G](tof, pt.planes)
  fwrite(when G == fgElectric: fmt27e else: fmt27, tof.alpha.float, tof.alphaErr,
         tof.att.float, tof.attErr)
  when G in {fgElectric, fgParallel}:
    fwrite(fmt28)
    fwrite(fmt29, tof.dl.float, tof.dlErr, tof.dt.float, tof.dtErr)
    fwrite(fmt30)
    fwrite(when G == fgElectric: fmt31e else: fmt31a, tof.w.z.float, tof.wErr.z)
    # α - η from the drift and growth of a Gaussian swarm
    let wrzn = tof.w.z.to(cm•s⁻¹)
    let fc1 = wrzn / (2.0 * tof.dl)
    let fc2 = (tof.alpha - tof.att).to(s⁻¹) / tof.dl
    fwrite(fmt888, float64(fc1 - sqrt(fc1*fc1 - fc2)))
    when G == fgParallel:
      r.transport.alpha = (tof.alpha / tof.w.z).to(cm⁻¹)
      r.transport.att = (tof.att / tof.w.z).to(cm⁻¹)
      # NIMBOLTZ FIX MB-42: errors in percent (was errAlpha·α/100)
      r.transport.alper = tof.alphaErr
      r.transport.atter = tof.attErr
  else:
    let d = tof.d
    let de = tof.dErr
    fwrite(fmt28b)
    when G == fgCrossed:
      fwrite(fmt29b, d.zz.float, de.zz, d.xx.float, de.xx, d.yy.float, de.yy, d.yz.float,
             de.yz)
    else:
      fwrite(fmt29c, d.zz.float, de.zz, d.xx.float, de.xx, d.yy.float, de.yy, d.yz.float,
             de.yz, d.xz.float, de.xz, d.xy.float, de.xy)
    fwrite(fmt30b)
    when G == fgCrossed:
      fwrite(fmt31b, tof.w.z.float, tof.wErr.z, tof.w.y.float, tof.wErr.y)
      let tofwr = sqrt(tof.w.z * tof.w.z + tof.w.y * tof.w.y)
    else:
      fwrite(fmt31c, tof.w.z.float, tof.wErr.z, tof.w.y.float, tof.wErr.y, tof.w.x.float,
             tof.wErr.x)
      let tofwr = sqrt(tof.w.z * tof.w.z + tof.w.y * tof.w.y + tof.w.x * tof.w.x)
    # as ALPCLCB/C: the drift velocity along E; ALPCLCBT/CT: the full one
    let wrzn = (when thermal: tofwr else: tof.w.z).to(cm•s⁻¹)
    let fc1 = wrzn / (2.0 * d.zz)
    let fc2 = (tof.alpha - tof.att).to(s⁻¹) / d.zz
    fwrite(fmt888, float64(fc1 - sqrt(fc1*fc1 - fc2)))
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
