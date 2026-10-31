## MONTE, MONTEA, MONTEB, MONTEC and the thermal versions MONTET, MONTEAT,
## MONTEBT, MONTECT of magboltz-11.21.f: bridge from the COMMON blocks to
## nimboltz/drift.

import ../frt, commons
import ../[tables, fields, collisions, freeflight, thermal, drift]

const
  fmt50 = "(\' WARNING ENERGY OUT OF RANGE,INCREASE ELECTRON ENERGY INTEGRATION RANGE \')"
  # progress table, one line per segment
  fmt201e = "(/,\'    VEL      POS        TIME      ENERGY   COUNT   DIFXX     DIFYY     DIFZZ\',/)"
  fmt201a = "(/,\'    VEL      POS        TIME      ENERGY   COUNT   DIFXX     DIFYY     DIFZZ  \',/)"
  fmt202a = "(1X,F8.2,2(1X,D10.3),F9.4,1X,I5,1X,3(2X,F8.1))"
  fmt201b = "(/,\'   VELZ    VELY  ENERGY   ID  DIFXX   DIFYY   DIFZZ   DIFYZ  DIFLNG  DIFTRN\',/)"
  fmt202b = "(2(F7.2,1X),F7.4,1X,I4,6(F7.1,1X))"
  fmt201c = "(/,\'   VELZ    VELY    VELX   ENERGY  ID  DIFXX  DIFYY  DIFZZ  DIFYZ  DIFXZ  DIFXY\',/)"
  fmt202c = "(3(F7.2,1X),1X,F7.4,I4,6(F7.1))"

var rings: array[FieldGeometry, array[bool, Ring]]
  ## The position arrays of the Fortran routines are static, so a run
  ## starts with the positions left by the previous run of the same routine;
  ## MONTEC reads them when NCOLM exceeds the collisions of the first two
  ## segments (reduced NSCALE only). One ring per routine keeps that.

func monteGeometry(bmag, btheta: float64): FieldGeometry =
  ## Routine chosen by Magboltz: unlike ELIMIT, B ∥ E at -180° goes to MONTEC
  if bmag == 0.0: fgElectric
  elif btheta == 0.0 or btheta == 180.0: fgParallel
  elif btheta == 90.0: fgCrossed
  else: fgOblique

proc progress[G: static FieldGeometry; T](d: Drift[G, T], seg: int) =
  let w = d.segmentVelocity
  let df = d.segmentDiffusion
  when G in {fgElectric, fgParallel}:
    if seg == 1: fwrite(when G == fgElectric: fmt201e else: fmt201a)
    fwrite(fmt202a, w.z, d.position.z, d.time, d.ave, d.ringFills, df.xx, df.yy, df.zz)
  elif G == fgCrossed:
    if seg == 1: fwrite(fmt201b)
    fwrite(fmt202b, w.z, w.y, d.ave, d.ringFills, df.xx, df.yy, df.zz, df.yz,
           d.difln, d.diftr)
  else:
    if seg == 1: fwrite(fmt201c)
    fwrite(fmt202c, w.z, w.y, w.x, d.ave, d.ringFills, df.xx, df.yy, df.zz, df.yz,
           df.xz, df.xy)

proc store(o: var (CommonOutpt | CommonOutptt), d: Drift, r: Transport) =
  o.time = d.times
  o.spec = d.spectrum
  o.tmax1 = d.tmax
  o.nnull = d.nnull
  o.ave = r.ave
  o.den = r.den
  if d.ringFills > 0: o.xid = float64(d.ringFills)
  let p = d.position
  o.x = p.x
  o.y = p.y
  o.z = p.z
  o.st = d.time

proc run[G: static FieldGeometry](thermalMotion: bool) =
  let f = initField(G, setp.efield, cnsts1.const1, cnsts1.const2, cnsts1.const3,
                    cnsts1.const5, bfld.wb, bfld.eovb, bfld.btheta)
  let start = initElectron[G](f, setp.estart, setp.theta, setp.phi)
  let dec = Decorrelation(ncolm: decor.ncolm, ncorln: decor.ncorln, ncorst: decor.ncorst)
  let collisions = setp.nmax div nSegments
  setp.small = small
  template ring: untyped = rings[G][thermalMotion]
  if ring.isEmpty: ring = initRing()

  template simulate(d, tabs, tcf: untyped) =
    for seg in 1 .. nSegments:
      d.runSegment(tabs, seg, collisions, ring, rng)
      d.endSegment(seg, mix2.eg, tcf)
      progress(d, seg)
      when G == fgElectric:
        # if the last energy bin is reached often, the table range is too small
        if d.spectrum[nEnergyBins] > 1000.0 * float64(seg):
          fwrite(fmt50)
          fstop()

  var r: Transport
  if thermalMotion:
    let gas = ThermalGas(tcfmx: thrm.tcfmx, tcfmxg: thrm.tcfmxg,
                         tcfmax: larget.tcfmax, vtmb: thrm.vtmb, count: inpt.ngas)
    var tcf: array[1..nEnergyBins, float64]
    for ik in 1 .. nEnergyBins:
      var tcfsum = 0.0
      for ki in 1 .. inpt.ngas: tcfsum = tcfsum + gasTable[ki].tcf[ik]
      tcf[ik] = tcfsum
    var d = initDrift[G, WarmGas](f, start, initWarmGas(gas, rng), inpt.estep, inpt.ngas,
                                  inpt.ipen != 0, anist.niso == 0, dec, bfld.btheta)
    simulate(d, gasTable, tcf)
    r = d.finish
    outptt.store(d, r)
    for g in 1 .. maxGases:
      for k in 1 .. 5: outptt.icoll[g, k] = d.tally[g].kinds[k]
      for i in 1 .. maxGasProcesses: outptt.icoln[g, i] = d.tally[g].processes[i]
    setp.theta = d.theta
  else:
    var d = initDrift[G, ColdGas](f, start, initColdGas(large.tcfmax, inpt.efinal),
                                  inpt.estep, inpt.ngas, inpt.ipen != 0, anis.niso == 0,
                                  dec, bfld.btheta)
    simulate(d, mixTable, mixTable.tcf)
    r = d.finish
    outpt.store(d, r)
    outpt.icoll = d.tally.kinds
    outpt.icoln = d.tally.processes
    large.tcfmax = d.target.bounds
    setp.theta = d.theta
  vel.wx = r.vel.x
  vel.wy = r.vel.y
  vel.wz = r.vel.z
  velerr.dwx = r.velErr.x
  velerr.dwy = r.velErr.y
  velerr.dwz = r.velErr.z
  diflab.difxx = r.dif.xx
  diflab.difyy = r.dif.yy
  diflab.difzz = r.dif.zz
  diflab.difyz = r.dif.yz
  diflab.difxy = r.dif.xy
  diflab.difxz = r.dif.xz
  diferb.dxxer = r.difErr.xx
  diferb.dyyer = r.difErr.yy
  diferb.dzzer = r.difErr.zz
  diferb.dyzer = r.difErr.yz
  diferb.dxyer = r.difErr.xy
  diferb.dxzer = r.difErr.xz
  difvel.difln = r.difln
  difvel.diftr = r.diftr
  diferl.dfler = r.dfler
  diferl.dfter = r.dfter
  ctowns.alpha = r.alpha
  ctowns.att = r.att
  ctwner.alper = r.alper
  ctwner.atter = r.atter

proc monte*(thermalMotion: bool) =
  ## the main Monte Carlo: transport coefficients of the mixture
  withGeometry(monteGeometry(bfld.bmag, bfld.btheta), G):
    run[G](thermalMotion)
