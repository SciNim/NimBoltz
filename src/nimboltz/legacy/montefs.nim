## MONTEFT, MONTEFTA, MONTEFTG, MONTEFTH, MONTEFD and the thermal versions
## MONTEFTT, MONTEFTAT, MONTEFTGT, MONTEFTHT, MONTEFDT of magboltz-11.21.f:
## bridge from the COMMON blocks to nimboltz/cascade. A run returns the
## plane sums for the analyses (legacy/alpcalcs.nim) and leaves the
## collision counts in /OUTPT/ or /OUTPTT/ and /TTRM/ for OUTPUT1 and COLF.
## The routines print the same messages with slightly different formats.

import ../frt, commons
import ../[tables, fields, collisions, freeflight, thermal, cascade, townsend]

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

proc printNegative(fm: Formats, f: Field, w: NegativeEnergy) =
  let ap = w.el.dir.z * f.f2 * sqrt(w.el.e)
  case fm.negativeKind
  of 1: fwrite(fm.negative, w.iter, w.fl.e, w.el.e, ap, f.bp, w.el.dir.z)
  of 2: fwrite(fm.negative, w.iter, w.fl.dz, w.fl.e, w.fl.cwt, w.fl.swt, f.wb * w.fl.t,
               w.el.vel.y)
  else: fwrite(fm.negative, w.fl.e, w.el.e, ap, f.bp, w.fl.t, w.el.dir.z, w.iter)

proc store(o: var (CommonOutpt | CommonOutptt), c: Cascade) =
  o.st = c.st
  o.x = c.pos.x
  o.y = c.pos.y
  o.z = c.pos.z
  o.tmax1 = c.tmax
  o.time = c.times
  o.spec = c.spectrum
  o.nnull = c.nnull
  o.xid = float64(c.id)

type
  AvalancheRun*[P: static PlaneKind] = object
    ## what the analyses need from a run
    when P == pkTime:
      sums*: TimeSums
    else:
      sums*: SpaceSums
    iprim*: int
    ratio*: AttachmentRatio

func formatsOf(G: FieldGeometry, P: PlaneKind, thermalMotion: bool): Formats =
  if P == pkSpace: formatsD
  else:
    case G
    of fgElectric: formatsT
    of fgParallel: (if thermalMotion: formatsAT else: formatsA)
    of fgCrossed, fgOblique: formatsG

proc avalanche*[G: static FieldGeometry; P: static PlaneKind](thermalMotion: bool,
                                                               planes: Planes): AvalancheRun[P] =
  ## one run of avalanches sampled at `planes`
  let fm = formatsOf(G, P, thermalMotion)
  let f = initField(G, setp.efield, cnsts1.const1, cnsts1.const2, cnsts1.const3,
                    cnsts1.const5, bfld.wb, bfld.eovb, bfld.btheta)
  let start = initElectron[G](f, setp.estart, setp.theta, setp.phi)
  setp.small = small
  setp.api = arccos(-1.0)

  template finish(c: untyped) =
    for w in c.warnings: printNegative(fm, f, w)
    case c.stop
    of stStackFull:
      fwrite(fm.stackFull, c.stacked + 1, c.id)
      fstop()
    of stTooManyPrimaries:
      fwrite(fm.tooMany, c.iprim)
      when P == pkSpace: fstop()
    of stNone: discard
    setp.theta = c.theta
    ttrm.ztot = c.ztot
    ttrm.ttot = c.ttot
    ttrm.ztots = c.ztots
    ttrm.ttots = c.ttots
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
      fwrite(fmt836, ebar, eerr)

  if thermalMotion:
    let gas = ThermalGas(tcfmx: thrm.tcfmx, tcfmxg: thrm.tcfmxg,
                         tcfmax: larget.tcfmax, vtmb: thrm.vtmb, count: inpt.ngas)
    var c = initCascade[G, P, WarmGas](f, cnsts1.const2, cnsts1.const3, start,
                                       initWarmGas(gas, rng), inpt.estep, setp.nmax,
                                       inpt.ipen != 0, planes, bfld.btheta)
    c.run(gasTable, rng)
    outptt.store(c)
    for g in 1 .. maxGases:
      for k in 1 .. 5: outptt.icoll[g, k] = c.tally[g].kinds[k]
      for i in 1 .. maxGasProcesses: outptt.icoln[g, i] = c.tally[g].processes[i]
    finish(c)
  else:
    var tlim = large.tcfmax[1]
    for i in 2 .. 8:
      if tlim < large.tcfmax[i]: tlim = large.tcfmax[i]
    var c = initCascade[G, P, ColdBound](f, cnsts1.const2, cnsts1.const3, start,
                                         ColdBound(tlim: tlim), inpt.estep, setp.nmax,
                                         inpt.ipen != 0, planes, bfld.btheta)
    c.run(mixTable, rng)
    outpt.store(c)
    outpt.icoll = c.tally.kinds
    outpt.icoln = c.tally.processes
    finish(c)
