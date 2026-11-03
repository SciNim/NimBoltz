## Swarm parameters from the avalanche runs of `cascade.nim` (ALPCALC and
## its variants, SST, PT/PTG/PTH, TOF/TOFG/TOFH, FRIEDLAND):
##
## - steady-state Townsend (SST): drift velocities, diffusion and Townsend
##   coefficients from the space planes (B = 0 only),
## - pulsed Townsend (PT): rate of change of the number of electrons, mean
##   energy and velocity at each time plane,
## - time of flight (TOF): drift and diffusion from the change of the mean
##   and spread of the positions between time planes,
## - Friedland's estimate: mean energy and ionisation and attachment
##   frequencies weighted by the energy spectrum of the avalanche.
##
## Each estimate also comes per plane; the final value is taken where the
## avalanche is best developed: planes 2 and 3 when it shrinks, the last
## planes when it grows, with half the difference between the two as the
## error. Formulas and their order of evaluation are those of the Fortran.
## Units as printed by Magboltz: µm/ns, cm²/s, 1/cm, eV, 10¹²/s.

import std/math
import vec3, fields, tables, cascade, drift

type
  AttachmentRatio* = object
    ## attachments per ionisation in an avalanche run, with the relative
    ## errors of both counts (ATTOINT, ATTERT, AIOERT); ratio = -1 marks a
    ## run without secondary electrons
    ratio*, attErr*, ionErr*: float64

  SstPlane* = object
    n*: int                       ## electrons created up to the plane (NESST)
    vd*, ws*: float64             ## mean velocity at the plane, from arrival times
    dl*, dt*: float64             ## longitudinal, transverse diffusion
    e*: float64                   ## mean energy
    alfn*, alfnj*, alfne*: float64  ## (α−η) from the flux, the current, the count

  SteadyState* = object
    ## SST (SSTOUT, and the table printed by SST)
    planes*: int                  ## planes with electrons (≤ 8)
    plane*: array[1..nPlanes, SstPlane]
    vd*, vdErr*, ws*, wsErr*, dl*, dlErr*, dt*, dtErr*: float64
    alpha*, alphaErr*, att*, attErr*: float64

  PtPlane* = object
    rate*: float64        ## (ionisation − attachment) frequency (10¹² /s)
    e*: float64           ## mean energy
    t*: float64           ## mean time
    v*: Vec3              ## mean velocity (fgOblique: in the frame of E)

  PulsedTownsend* = object
    ## PT (PTTOF, PTTOFG, PTTOFH)
    planes*: int          ## planes with electrons (≤ ITFINAL)
    plane*: array[1..nPlanes, PtPlane]

  TofPlane* = object
    w*: Vec3              ## drift velocity between the previous and this plane
    d*: Tensor            ## diffusion

  TimeOfFlight* = object
    ## TOF (TOFOUT, TOFGOUT, TOFHOUT); for B = 0 and B ∥ E the transverse
    ## diffusion is (DXX + DYY)/2 in `dt`, the drift velocity WR in `w.z`
    plane*: array[1..nPlanes, TofPlane]
    e*, eErr*: float64                  ## mean energy
    wv*, wvErr*: Vec3               ## mean velocity
    d*, dErr*: Tensor               ## diffusion
    dl*, dlErr*, dt*, dtErr*: float64   ## B = 0, B ∥ E
    w*, wErr*: Vec3                 ## drift velocity
    alpha*, alphaErr*: float64          ## ionisation rate (10¹² /s)
    att*, attErr*: float64              ## attachment rate

  Friedland* = object
    e*, ionisation*, attachment*: float64  ## eV, 10¹² /s

func cube(x: float64): float64 {.inline.} =
  ## x**3 as gfortran evaluates it (libgcc __powidf2)
  x * (x * x)

func sq(x: float64): float64 {.inline.} = x * x

func attachmentRatio*(neion, nelec, iprim: int): AttachmentRatio =
  ## from the counts of a run (MONTEFT*, MONTEFD)
  let aneion = float64(neion)
  if nelec > iprim:
    let anbt = float64(nelec - iprim)
    result.ratio = aneion / anbt
    result.attErr = sqrt(aneion) / aneion
    result.ionErr = sqrt(anbt) / anbt
  else:
    result.ratio = -1.0
    result.attErr = sqrt(aneion) / aneion

func samplingPlanes*(alpha, att, wz, tempc, torr: float64):
    tuple[planes: Planes, alphast, vdst: float64] =
  ## Planes for an expected gain of 3 between neighbours, from the Townsend
  ## coefficients of MONTE (1/cm) and the drift velocity wz (cm/s). A large
  ## cancellation of α and η shortens the spacing.
  var alphast = 0.85 * abs(alpha - att)
  # NIMBOLTZ FIX MB-41: ALPCALCT (thermal, B = 0) lacked the shorter spacing
  let corr = 760.0 * (tempc + 273.1499938964844) / (torr * 293.1499938964844)
  if alpha > 10.0 * alphast or att > 10.0 * alphast:
    if alpha > 100.0 * corr:
      alphast = alphast * 30.0
    elif alpha > 50.0 * corr:
      alphast = alphast * 15.0
    else:
      alphast = alphast * 8.0
  let vdst = wz * 0.00001
  var tstep = ln(3.0) / (alphast * vdst * 100000.0)
  var zstep = ln(3.0) / alphast
  tstep = tstep * 1000000000000.0
  zstep = zstep * 0.01
  var p = Planes(tstep: tstep, tfinal: 7.0 * tstep, itfinal: 7,
                 zstep: zstep, zfinal: 8.0 * zstep, izfinal: 8)
  p.zplane[1] = zstep
  for k in 2 .. nPlanes: p.zplane[k] = float64(k) * zstep
  (p, alphast, vdst)

# --- steady-state Townsend --------------------------------------------------------

func steadyState*(s: SpaceSums, iprim: int, zstep: float64, izfinal: int,
                  ar: AttachmentRatio): SteadyState =
  ## SST from the sums of an `pkSpace` run
  var nesst: array[1..nPlanes + 1, int] = s.born
  var esst, vdsst, wtemp, wssst, dxsst, dysst, dlsst, alfn, alfnj, alfne:
    array[1..nPlanes, float64]
  var jprint = izfinal
  # electrons per plane: primaries plus those created minus attached before it
  nesst[1] = iprim + s.born[1]
  for k in 2 .. jprint: nesst[k] = nesst[k - 1] + s.born[k]
  for i in 1 .. jprint:
    if nesst[i] == 0:
      jprint = i - 1
      break
  let ts = s.w
  esst[1] = s.e[1] / ts[1]
  vdsst[1] = s.vz[1] / ts[1]
  wtemp[1] = zstep * ts[1] / s.t[1]
  wssst[1] = wtemp[1]
  dxsst[1] = (s.xx[1] / ts[1] - sq(s.x[1] / ts[1])) * wssst[1] / (2.0 * zstep)
  dysst[1] = (s.yy[1] / ts[1] - sq(s.y[1] / ts[1])) * wssst[1] / (2.0 * zstep)
  dlsst[1] = (s.tt[1] / ts[1] - sq(s.t[1] / ts[1])) * cube(wssst[1]) / (2.0 * zstep)
  if nesst[1] != 0:
    alfne[1] = (ln(float64(nesst[1])) - ln(float64(iprim))) / zstep
  for i in 2 .. jprint:
    esst[i] = s.e[i] / ts[i]
    vdsst[i] = s.vz[i] / ts[i]
    wtemp[i] = zstep * float64(i) * ts[i] / s.t[i]
    wssst[i] = wtemp[i] * wtemp[i - 1] / (float64(i) * wtemp[i - 1] - float64(i - 1) * wtemp[i])
    dxsst[i] = (s.xx[i] / ts[i] - sq(s.x[i] / ts[i]) - s.xx[i - 1] / ts[i - 1] +
                sq(s.x[i - 1] / ts[i - 1])) * wssst[i] / (2.0 * zstep)
    dysst[i] = (s.yy[i] / ts[i] - sq(s.y[i] / ts[i]) - s.yy[i - 1] / ts[i - 1] +
                sq(s.y[i - 1] / ts[i - 1])) * wssst[i] / (2.0 * zstep)
    dlsst[i] = (s.tt[i] / ts[i] - sq(s.t[i] / ts[i]) - s.tt[i - 1] / ts[i - 1] +
                sq(s.t[i - 1] / ts[i - 1])) * cube(wssst[i]) / (2.0 * zstep)
    alfn[i] = (ln(ts[i]) - ln(ts[i - 1])) / zstep
    alfnj[i] = (ln(ts[i] * vdsst[i]) - ln(ts[i - 1] * vdsst[i - 1])) / zstep
    alfne[i] = 0.0
    if not (nesst[i] == 0 or nesst[i - 1] == 0):
      alfne[i] = (ln(float64(nesst[i])) - ln(float64(nesst[i - 1]))) / zstep
  # over all planes
  let j = jprint
  var dxfin = (s.xx[j] / ts[j] - sq(s.x[j] / ts[j])) * wssst[j] / (float64(j) * 2.0 * zstep)
  dxfin = dxfin * 10000000000000000.0
  var dyfin = (s.yy[j] / ts[j] - sq(s.y[j] / ts[j])) * wssst[j] / (float64(j) * 2.0 * zstep)
  dyfin = dyfin * 10000000000000000.0
  var dlfin = (s.tt[j] / ts[j] - sq(s.t[j] / ts[j])) * cube(wssst[j]) / (float64(j) * 2.0 * zstep)
  dlfin = dlfin * 10000000000000000.0
  for k in 1 .. jprint:
    vdsst[k] = vdsst[k] * 1000000000.0
    wssst[k] = wssst[k] * 1000000000.0
    dxsst[k] = dxsst[k] * 10000000000000000.0
    dysst[k] = dysst[k] * 10000000000000000.0
    dlsst[k] = dlsst[k] * 10000000000000000.0
    alfn[k] = alfn[k] * 0.01
    alfnj[k] = alfnj[k] * 0.01
    alfne[k] = alfne[k] * 0.01
  var drsst: array[1..nPlanes, float64]
  result.planes = jprint
  for k in 1 .. jprint:
    drsst[k] = (dxsst[k] + dysst[k]) / 2.0
    result.plane[k] = SstPlane(n: nesst[k], vd: vdsst[k], ws: wssst[k], dl: dlsst[k],
                               dt: drsst[k], e: esst[k], alfn: alfn[k], alfnj: alfnj[k],
                               alfne: alfne[k])
  func relErr(a, b: float64): float64 = 100.0 * abs((a - b) / (2.0 * a))
  # error of the growth rate from the counting statistics of planes 2 and 3
  func countErr(n2, n3: int): float64 =
    let anst2 = float64(n2)
    let anst3 = float64(n3)
    let anst4 = anst3 - sqrt(anst3)
    let anst5 = ln(anst2 / anst3)
    let anst6 = ln(anst2 / anst4)
    let anst7 = anst6 / anst5
    anst7 - 1.0
  if nesst[1] > nesst[5]:
    # shrinking avalanche: planes 2 and 3
    result.vd = vdsst[2]
    result.vdErr = relErr(vdsst[2], vdsst[3])
    result.ws = wssst[2]
    result.wsErr = relErr(wssst[2], wssst[3])
    result.dl = dlsst[2]
    result.dlErr = relErr(dlsst[2], dlsst[3])
    result.dt = drsst[2]
    result.dtErr = relErr(drsst[2], drsst[3])
    let anst8 = countErr(nesst[2], nesst[3])
    if ar.ratio == -1.0:
      result.attErr = 100.0 * sqrt(anst8*anst8 + ar.attErr*ar.attErr)
      result.att = -((alfn[2] + alfnj[2] + alfne[2]) / 3.0)
    else:
      let atmp = (alfn[2] + alfnj[2] + alfne[2]) / 3.0
      result.alpha = atmp / (1.0 - ar.ratio)
      result.alphaErr = 100.0 * sqrt(anst8*anst8 + ar.ionErr*ar.ionErr)
      result.att = ar.ratio * atmp / (1.0 - ar.ratio)
      result.attErr = 100.0 * sqrt(anst8*anst8 + ar.attErr*ar.attErr)
  else:
    # growing avalanche: planes 8 and 7, and all planes
    result.vd = vdsst[8]
    result.vdErr = 100.0 * abs((vdsst[8] - vdsst[7]) / (2.0 * vdsst[8]))
    result.ws = wssst[8]
    result.wsErr = 100.0 * abs((wssst[8] - wssst[7]) / (2.0 * wssst[8]))
    result.dl = dlfin
    result.dlErr = 100.0 * abs((dlfin - dlsst[8]) / (2.0 * dlfin))
    result.dt = (dxfin + dyfin) / 2.0
    result.dtErr = 100.0 * abs((result.dt - drsst[8]) / (2.0 * result.dt))
    let atmp = (alfn[8] + alfnj[8] + alfne[8]) / 3.0
    let atmp2 = (alfn[7] + alfnj[7] + alfne[7]) / 3.0
    let ater = abs((atmp - atmp2) / (2.0 * atmp))
    result.alpha = atmp / (1.0 - ar.ratio)
    result.alphaErr = 100.0 * sqrt(ater*ater + ar.ionErr*ar.ionErr)
    result.att = ar.ratio * atmp / (1.0 - ar.ratio)
    if ar.ratio != 0.0:
      result.attErr = 100.0 * sqrt(ater*ater + ar.attErr*ar.attErr)

# --- pulsed Townsend and time of flight -------------------------------------------

func pulsedTownsend*[G: static FieldGeometry](s: TimeSums, iprim: int, tstep: float64,
                                              itfinal: int): PulsedTownsend =
  ## PT from the sums of a `pkTime` run; stops at the first empty plane
  result.planes = itfinal
  var antpl: array[1..nPlanes, float64]
  template plane(i: int, previous: float64) =
    antpl[i] = float64(s.n[i])
    result.plane[i].rate = (ln(antpl[i]) - ln(previous)) / tstep
    result.plane[i].e = s.e[i] / antpl[i]
    result.plane[i].t = s.t[i] / antpl[i]
    result.plane[i].v.z = 1000000000.0 * s.vz[i] / antpl[i]
    when G in {fgCrossed, fgOblique}:
      result.plane[i].v.y = 1000000000.0 * s.vy[i] / antpl[i]
    when G == fgOblique:
      result.plane[i].v.x = 1000000000.0 * s.vx[i] / antpl[i]
  plane(1, float64(iprim))
  for i in 2 .. itfinal:
    if s.n[i] == 0:
      result.planes = i - 1
      break
    plane(i, antpl[i - 1])

func timeOfFlight*[G: static FieldGeometry](s: TimeSums, pt: PulsedTownsend, tstep: float64,
                                            ar: AttachmentRatio): TimeOfFlight =
  ## TOF from the sums of a `pkTime` run and its PT; `ar` as TOFG uses it
  ## (with FRATT/FRION, MB-40)
  let itfinal = pt.planes
  var antpl: array[1..nPlanes, float64]
  var w: array[1..nPlanes, Vec3]
  var d: array[1..nPlanes, Tensor]
  antpl[1] = float64(s.n[1])
  let n1 = antpl[1]
  w[1].z = s.z[1] / (n1 * tstep)
  when G in {fgCrossed, fgOblique}:
    w[1].y = s.y[1] / (n1 * tstep)
  when G == fgOblique:
    w[1].x = s.x[1] / (n1 * tstep)
  d[1].zz = (s.zz[1] / n1 - sq(s.z[1] / n1)) / (2.0 * tstep)
  d[1].xx = (s.xx[1] / n1 - sq(s.x[1] / n1)) / (2.0 * tstep)
  d[1].yy = (s.yy[1] / n1 - sq(s.y[1] / n1)) / (2.0 * tstep)
  when G in {fgCrossed, fgOblique}:
    d[1].yz = (s.yz[1] / n1 - s.y[1] * s.z[1] / (n1*n1)) / (2.0 * tstep)
  when G == fgOblique:
    d[1].xz = (s.xz[1] / n1 - s.x[1] * s.z[1] / (n1*n1)) / (2.0 * tstep)
    d[1].xy = (s.xy[1] / n1 - s.x[1] * s.y[1] / (n1*n1)) / (2.0 * tstep)
  for i in 2 .. itfinal:
    antpl[i] = float64(s.n[i])
    let a = antpl[i]
    let b = antpl[i - 1]
    w[i].z = (s.z[i] / a - s.z[i - 1] / b) / tstep
    when G in {fgCrossed, fgOblique}:
      w[i].y = (s.y[i] / a - s.y[i - 1] / b) / tstep
    when G == fgOblique:
      w[i].x = (s.x[i] / a - s.x[i - 1] / b) / tstep
    d[i].zz = (s.zz[i] / a - sq(s.z[i] / a) - s.zz[i - 1] / b + sq(s.z[i - 1] / b)) / (2.0 * tstep)
    d[i].xx = (s.xx[i] / a - sq(s.x[i] / a) - s.xx[i - 1] / b + sq(s.x[i - 1] / b)) / (2.0 * tstep)
    d[i].yy = (s.yy[i] / a - sq(s.y[i] / a) - s.yy[i - 1] / b + sq(s.y[i - 1] / b)) / (2.0 * tstep)
    when G == fgCrossed:
      d[i].yz = (s.yz[i] / a - s.y[i] * s.z[i] / (a*a) - s.yz[i - 1] / b +
                 s.y[i - 1] * s.z[i - 1] / sq(b)) / (2.0 * tstep)
    elif G == fgOblique:
      d[i].yz = (s.yz[i] / a - s.y[i] * s.z[i] / (a*a) - s.yz[i - 1] / b +
                 s.y[i - 1] * s.z[i - 1] / (b*b)) / (2.0 * tstep)
      d[i].xz = (s.xz[i] / a - s.x[i] * s.z[i] / (a*a) - s.xz[i - 1] / b +
                 s.x[i - 1] * s.z[i - 1] / (b*b)) / (2.0 * tstep)
      d[i].xy = (s.xy[i] / a - s.x[i] * s.y[i] / (a*a) - s.xy[i - 1] / b +
                 s.x[i - 1] * s.y[i - 1] / (b*b)) / (2.0 * tstep)
  for i in 1 .. itfinal:
    w[i].z = w[i].z * 1000000000.0
    when G in {fgCrossed, fgOblique}:
      w[i].y = w[i].y * 1000000000.0
    when G == fgOblique:
      w[i].x = w[i].x * 1000000000.0
    d[i].zz = d[i].zz * 10000000000000000.0
    d[i].xx = d[i].xx * 10000000000000000.0
    d[i].yy = d[i].yy * 10000000000000000.0
    when G in {fgCrossed, fgOblique}:
      d[i].yz = d[i].yz * 10000000000000000.0
    when G == fgOblique:
      d[i].xz = d[i].xz * 10000000000000000.0
      d[i].xy = d[i].xy * 10000000000000000.0
  for i in 1 .. itfinal:
    result.plane[i] = TofPlane(w: w[i], d: d[i])
  # final values: planes 2 and 3 when the avalanche shrinks, the last two
  # when it grows; error = half the difference (TOFG growing: the full one)
  let shrinking = s.n[1] > s.n[itfinal]
  let (i1, i2) = (if shrinking: (2, 3) else: (itfinal, itfinal - 1))
  let half = (when G == fgCrossed: (if shrinking: 2.0 else: 1.0) else: 2.0)
  template err(a, b: float64): float64 =
    (when G == fgCrossed: 100.0 * abs((a - b) / (half * a))
     else: 100.0 * abs((a - b) / (2.0 * a)))
  let p1 = pt.plane[i1]
  let p2 = pt.plane[i2]
  result.e = p1.e
  result.eErr = err(p1.e, p2.e)
  result.wv.z = p1.v.z
  result.wvErr.z = err(p1.v.z, p2.v.z)
  when G in {fgCrossed, fgOblique}:
    result.wv.y = p1.v.y
    result.wvErr.y = err(p1.v.y, p2.v.y)
  when G == fgOblique:
    result.wv.x = p1.v.x
    result.wvErr.x = err(p1.v.x, p2.v.x)
  when G in {fgElectric, fgParallel}:
    result.dl = d[i1].zz
    result.dlErr = err(d[i1].zz, d[i2].zz)
    let tdt1 = (d[i1].xx + d[i1].yy) / 2.0
    let tdt2 = (d[i2].xx + d[i2].yy) / 2.0
    result.dt = tdt1
    result.dtErr = err(tdt1, tdt2)
  else:
    result.d.zz = d[i1].zz
    result.dErr.zz = err(d[i1].zz, d[i2].zz)
    result.d.xx = d[i1].xx
    result.dErr.xx = err(d[i1].xx, d[i2].xx)
    result.d.yy = d[i1].yy
    result.dErr.yy = err(d[i1].yy, d[i2].yy)
    result.d.yz = d[i1].yz
    result.dErr.yz = err(d[i1].yz, d[i2].yz)
    when G == fgOblique:
      result.d.xz = d[i1].xz
      result.dErr.xz = err(d[i1].xz, d[i2].xz)
      result.d.xy = d[i1].xy
      result.dErr.xy = err(d[i1].xy, d[i2].xy)
  result.w.z = w[i1].z
  when G in {fgCrossed, fgOblique}:
    result.w.y = w[i1].y
  when G == fgOblique:
    result.w.x = w[i1].x
  result.wErr.z = err(w[i1].z, w[i2].z)
  when G in {fgCrossed, fgOblique}:
    result.wErr.y = err(w[i1].y, w[i2].y)
  when G == fgOblique:
    result.wErr.x = err(w[i1].x, w[i2].x)
  # ionisation and attachment rates
  if shrinking:
    let anst2 = float64(s.n[2])
    let anst3 = float64(s.n[3])
    let anst4 = anst3 - sqrt(anst3)
    let anst5 = ln(anst2 / anst3)
    let anst6 = ln(anst2 / anst4)
    let anst7 = anst6 / anst5
    let anst8 = anst7 - 1.0
    if ar.ratio == -1.0:
      result.att = -pt.plane[2].rate
      result.attErr = 100.0 * sqrt(anst8*anst8 + ar.attErr*ar.attErr)
    else:
      result.alpha = pt.plane[2].rate / (1.0 - ar.ratio)
      result.alphaErr = 100.0 * sqrt(anst8*anst8 + ar.ionErr*ar.ionErr)
      result.att = ar.ratio * pt.plane[2].rate / (1.0 - ar.ratio)
      result.attErr = 100.0 * sqrt(anst8*anst8 + ar.attErr*ar.attErr)
  else:
    let ater = (when G == fgCrossed: abs((pt.plane[i1].rate - pt.plane[i2].rate) / pt.plane[i1].rate)
                else: abs((pt.plane[i1].rate - pt.plane[i2].rate) / (2.0 * pt.plane[i1].rate)))
    result.alpha = pt.plane[i1].rate / (1.0 - ar.ratio)
    result.alphaErr = 100.0 * sqrt(ater*ater + ar.ionErr*ar.ionErr)
    result.att = ar.ratio * pt.plane[i1].rate / (1.0 - ar.ratio)
    if ar.ratio != 0.0:
      result.attErr = 100.0 * sqrt(ater*ater + ar.attErr*ar.attErr)

# --- Friedland ---------------------------------------------------------------------

func friedland*(spectrum, energies, tcf, fcion, fcatt: array[1..nEnergyBins, float64],
                ttots: float64): Friedland =
  ## mean energy and ionisation and attachment frequencies over the time
  ## the avalanche spent in each energy bin (spectrum / tcf)
  var ebar, alfbar, attbar = 0.0
  for i in 1 .. nEnergyBins:
    ebar = ebar + energies[i] * spectrum[i] / tcf[i]
    alfbar = alfbar + fcion[i] * spectrum[i] / tcf[i]
    attbar = attbar + fcatt[i] * spectrum[i] / tcf[i]
  Friedland(e: ebar / ttots, ionisation: alfbar / ttots, attachment: attbar / ttots)
