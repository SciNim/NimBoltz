## MONTEFT, MONTEFTA, MONTEFTG, MONTEFTH, MONTEFD and the thermal versions
## MONTEFTT, MONTEFTAT, MONTEFTGT, MONTEFTHT, MONTEFDT of magboltz-11.21.f:
## bridge from the COMMON blocks to nimboltz/cascade. The routines print
## the same messages with slightly different formats; each entry point
## passes its own.

import ../frt, commons
import ../[tables, fields, collisions, freeflight, thermal, cascade]

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

template storeTime(o: untyped, s: TimeSums) =
  o.etpl = s.e
  o.xtpl = s.x
  o.ytpl = s.y
  o.ztpl = s.z
  o.ttpl = s.t
  o.xxtpl = s.xx
  o.yytpl = s.yy
  o.zztpl = s.zz
  o.vztpl = s.vz
  o.netpl = s.n

proc storePlanes[G: static FieldGeometry; T](c: Cascade[G, pkTime, T]) =
  let aneion = float64(c.neion)
  template ratios(o: untyped) =
    if c.nelec > c.iprim:
      let anbt = float64(c.nelec - c.iprim)
      o.attoint = aneion / anbt
      o.attert = sqrt(aneion) / aneion
      o.aioert = sqrt(anbt) / anbt
    else:
      o.attoint = -1.0
      o.attert = sqrt(aneion) / aneion
  when G in {fgElectric, fgParallel}:
    tplout.storeTime(c.sums)
    ratios(tplout)
  elif G == fgCrossed:
    tploutg.storeTime(c.sums)
    tploutg.yztpl = c.sums.yz
    tploutg.vytpl = c.sums.vy
    ratios(tploutg)
  else:
    tplouth.storeTime(c.sums)
    tplouth.yztpl = c.sums.yz
    tplouth.xztpl = c.sums.xz
    tplouth.xytpl = c.sums.xy
    tplouth.vytpl = c.sums.vy
    tplouth.vxtpl = c.sums.vx
    ratios(tplouth)

proc storePlanes[G: static FieldGeometry; T](c: Cascade[G, pkSpace, T]) =
  let s = c.sums
  splout.espl = s.e
  splout.xspl = s.x
  splout.yspl = s.y
  splout.zspl = s.z
  splout.tspl = s.rate
  splout.xxspl = s.xx
  splout.yyspl = s.yy
  splout.zzspl = s.zz
  splout.vzspl = s.vz
  splout.tssum = s.w
  splout.tssum2 = s.w2
  splout.nesst = s.born
  spl1.tmspl = s.t
  spl1.ttmspl = s.tt
  spl1.rspl = s.r
  spl1.rrspl = default(array[1..8, float64])
  spl1.rrspm = s.rr
  let aneion = float64(c.neion)
  if c.nelec > c.iprim:
    let anbt = float64(c.nelec - c.iprim)
    splout.attoion = aneion / anbt
    splout.attater = sqrt(aneion) / aneion
    splout.attioer = sqrt(anbt) / anbt
  else:
    splout.attoion = -1.0
    splout.attater = sqrt(aneion) / aneion

proc run[G: static FieldGeometry; P: static PlaneKind](thermalMotion: bool, fm: Formats) =
  let f = initField(G, setp.efield, cnsts1.const1, cnsts1.const2, cnsts1.const3,
                    cnsts1.const5, bfld.wb, bfld.eovb, bfld.btheta)
  let start = initElectron[G](f, setp.estart, setp.theta, setp.phi)
  setp.small = small
  setp.api = arccos(-1.0)
  let planes = Planes(tstep: cion.tstep, tfinal: cion.tfinal, itfinal: cion.itfinal,
                      zstep: cion.zstep, zfinal: cion.zfinal, izfinal: ctcalc.izfinal,
                      zplane: [ctcalc.zplane1, ctcalc.zplane2, ctcalc.zplane3, ctcalc.zplane4,
                               ctcalc.zplane5, ctcalc.zplane6, ctcalc.zplane7, ctcalc.zplane8])

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
    cion.iprim = c.iprim
    ttrm.ztot = c.ztot
    ttrm.ttot = c.ttot
    ttrm.ztots = c.ztots
    ttrm.ttots = c.ttots
    storePlanes(c)
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

proc monteft*() = run[fgElectric, pkTime](false, formatsT)
proc montefta*() = run[fgParallel, pkTime](false, formatsA)
proc monteftg*() = run[fgCrossed, pkTime](false, formatsG)
proc montefth*() = run[fgOblique, pkTime](false, formatsG)
proc montefd*() = run[fgElectric, pkSpace](false, formatsD)
proc monteftt*() = run[fgElectric, pkTime](true, formatsT)
proc monteftat*() = run[fgParallel, pkTime](true, formatsAT)
proc monteftgt*() = run[fgCrossed, pkTime](true, formatsG)
proc monteftht*() = run[fgOblique, pkTime](true, formatsG)
proc montefdt*() = run[fgElectric, pkSpace](true, formatsD)
