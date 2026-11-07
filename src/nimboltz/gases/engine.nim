## One routine for the cross sections of every gas: `fill` reads a
## `GasModel` (model.nim) and fills a `GasXsec` on the energy grid, as the
## gas routines GASn of Magboltz do.
##
## `fill` takes the model as a static parameter: each gas gets its own
## instance, and the variants it does not use (`when`) are not compiled
## in. The physics is the same as Magboltz's; the floating point operations
## are not always in the same order, so the results agree with the Fortran
## routines to a few units of the last digit (tests/gasengine).

import model
export model

const
  bornIonisation = 1.873884e-20.cm²   ## Born-Bethe ionisation constant

type
  Context = object
    ## what the cross sections of a bin depend on besides their tables
    r: Relativistic
    am2, c: float64      ## Bethe parameters of the total ionisation
    eion: eV             ## ionisation energy, E(3)
    negative: Negative
    qcount: cm²          ## counting ionisation (`countingShare`)

proc betheIonisation(c: Context, fraction = 1.0'f32): cm² =
  ## Bethe's ionisation cross section times a share of it
  let x2 = 1.0 / c.r.beta2
  let x1 = x2 * ln(c.r.beta2 / (1.0 - c.r.beta2)) - 1.0
  bornIonisation * (c.am2 * x1 + c.c * x2) * fraction

proc bornExcitation(c: Context, f: float32, ein, en: eV): cm² =
  ## Born-Bethe cross section of an allowed transition with oscillator
  ## strength f and energy loss ein
  let bb = bbconst.cm²•eV
  let r = c.r
  let d = float64(f) / (ein * r.beta2) *
    (ln(r.beta2 * r.gamma2 * emass2 / (4.0 * ein.float)) - r.beta2) * bb *
    (en / (en + ein + c.eion))
  case c.negative
  of keepNegative: d
  of clampNegative: max(d, 0.0.cm²)
  of absNegative: abs(d)

func interpolate(s: Sigma, en: eV): float64 {.inline.} =
  case s.interp
  of linear: lin(s.x, s.y, en.float)
  of logLog: logLog(s.x, s.y, en.float)
  of linLog: raiseAssert "linLog only for elastic tables"

proc at(s: Sigma, en: eV, ein: eV, c: Context): cm² =
  ## the cross section at energy en (ein: energy loss for `born`)
  let x = en.float
  if s.x.len > 0 and (s.beyond in {extrapolate, power} or x <= s.x[^1]):
    var v = if s.beyond == power and x > s.x[^1]:
              s.y[^1] * power(s.x[^1] / x, s.n)
            else: s.interpolate(en)
    if s.clamp and v < 0.0: v = 0.0
    return v * s.unit * float64(s.scale)
  case s.beyond
  of born: c.bornExcitation(s.f, ein, en) * float64(s.bornScale)
  of bethe: c.betheIonisation(s.fraction)
  of countingShare: c.qcount * float64(s.fraction)
  of extrapolate, power, zero: 0.0.cm²

# ------------------------------------------------------------ elastic

proc mert(p: Mert, en: eV): tuple[qela, qmom: cm²] =
  ## elastic and momentum transfer cross sections from the phase shifts of
  ## the modified effective range theory (partial waves up to l = 99; the
  ## higher ones in single precision, as Magboltz)
  const lmax = 100
  let (apol, a, d, f) = (float64(p.apol), float64(p.a), float64(p.d), float64(p.f))
  let (a1, b1, a2) = (float64(p.a1), float64(p.b1), float64(p.a2))
  let api = arccos(-1.0)
  let ak = sqrt(en / ary.eV)
  let ak2 = ak * ak
  let ak3 = ak2 * ak
  let ak4 = ak3 * ak
  let ak5 = ak4 * ak
  template phase(x: float64): float64 =
    (if p.tangent: arctan(x) else: x)
  let an0 = phase(-(a * ak * (1.0 + 4.0 * apol / 3.0 * ak2 * ln(ak))) -
                  api * apol / 3.0 * ak2 + d * ak3 + f * ak4)
  let an1 = phase((api / 15.0 * apol * ak2 - a1 * ak3) / (1.0 + b1 * ak2))
  let an2 = phase(api * apol * ak2 / 105.0 - a2 * ak5)
  var anhigh = an2
  var sum = pow2(sin(an0 - an1)) + 2.0 * pow2(sin(an1 - an2))
  var sigel = pow2(sin(an0)) + 3.0 * pow2(sin(an1))
  for j in 2 .. lmax - 1:
    let l = float32(j)
    let anlow = anhigh
    let dl = (2.0'f32 * l + 5.0'f32) * (2.0'f32 * l + 3.0'f32) * (2.0'f32 * l + 1.0'f32)
    let sumi = float64(6.0'f32 / (dl * (2.0'f32 * l - 1.0'f32)))
    sum = sum + float64(l + 1.0'f32) * pow2(sin(phase(api * apol * ak2 * sumi)))
    anhigh = phase(api * apol * ak2 / float64(dl))
    sigel = sigel + float64(2.0'f32 * l + 1.0'f32) * pow2(sin(anlow))
  let q0 = 4.0 * pir2.cm² / ak2
  (sigel * q0, sum * q0)

proc elasticAt(m: static Elastic, en: eV): tuple[qela, qmom: cm², pq2: float64] =
  ## total and momentum transfer cross sections, anisotropy parameter of
  ## model 2
  let x = en.float
  template interp(kind: Interp, xs, ys: seq[float64]): float64 =
    when kind == logLog: logLog(xs, ys, x)
    elif kind == linLog:
      if x <= m.switchAt: lin(xs, ys, x) else: logLog(xs, ys, x)
    else: lin(xs, ys, x)
  template eps: float64 = 1.0 - interp(m.epsInterp, m.xEps, m.eps)
  when m.low == phaseShifts:
    if x <= m.lowBelow:
      let (qela, qmom) = mert(m.mert, en)
      return (qela, qmom, eps)
  elif m.low == constant:
    if x <= m.lowBelow:
      return (m.lowQ.cm², m.lowQ.cm², (when m.lowIsotropic: 0.0 else: eps))
  template q(xs, ys: seq[float64]): cm² = interp(m.interp, xs, ys) * 1e-16.cm²
  (q(m.xTotal, m.total), q(m.xMomentum, m.momentum), eps)

# ------------------------------------------------------------ the gas

proc describe(m: static GasModel, inp: GasInput, g: var GasXsec) =
  ## names, energies and descriptions of the processes
  let aniso = inp.naniso != 0
  g.name = if aniso: m.name.anisotropic else: m.name.isotropic
  g.scrp[2] = if aniso: m.elasticText.anisotropic else: m.elasticText.isotropic
  g.e[2] = 2.0 * emass / (float64(m.mass) * amu)
  g.e[3] = m.e3
  g.e[4] = m.e4
  for k in 1 .. 6: g.kel[k] = inp.naniso
  const chans = m.ionisation.channels
  g.nion = chans.len
  for k, ch in chans:
    g.eion[k + 1] = ch.eion
    if ch.eobyScale.len > 0:
      var eb = float64(ch.eion)
      for f in ch.eobyScale: eb = eb * float64(f)
      g.eb[k + 1] = eb
    else: g.eb[k + 1] = ch.eoby
    g.nc0[k + 1] = ch.nc0
    g.ec0[k + 1] = ch.ec0
    g.scrp[2 + k + 1] = ch.text
  g.natt = m.attachment.len
  for k, a in m.attachment: g.scrp[2 + g.nion + k + 1] = a.text
  g.nin = m.levels.len
  for j, lv in m.levels:
    g.ei[j + 1] = lv.ein
    g.kin[j + 1] = if lv.kin < 0: inp.naniso else: lv.kin
    g.scrp[4 + g.nion + g.natt + j + 1] = lv.text
  let first = if m.penning.first > 0: m.penning.first else: 1
  let last = if m.penning.last > 0: m.penning.last else: g.nin
  g.penning(inp, first .. last, m.penning.fraction, m.penning.length,
            m.penning.time, m.penning.format, m.penning.all)

proc ionisation(m: static Ionisation, g: var GasXsec, i: int, en: eV,
                c: Context, iso: float64, ioffion: seq[int]) =
  when m.scheme == grossCounting:
    var gross, counting = 0.0.cm²
    g.peqion[1, i] = iso
    g.peqion[2, i] = iso
    if en > g.eion[1].eV:
      if en.float > m.counting.x[^1]:
        counting = c.betheIonisation()
        gross = counting / float64(m.grossBethe)
      else:
        gross = m.gross.interpolate(en) * m.gross.unit
        counting = m.counting.interpolate(en) * m.counting.unit
      for k in 1 .. 2:
        if en > m.channels[k - 1].gate * g.eion[k].eV:
          g.peqion[k, i] = g.peqel[2, i - ioffion[k - 1]]
    g.qion[1, i] = float(2.0 * counting - gross)
    g.qion[2, i] = float(max(gross - counting, 0.0.cm²))
  else:
    for k, ch in m.channels:
      let k1 = k + 1
      g.qion[k1, i] = 0.0
      g.peqion[k1, i] = iso
      let eion = ch.eion.float64.eV
      if (if ch.openInclusive: en < eion else: en <= eion): continue
      var q = 0.0.cm²
      for t in ch.terms:
        if t.series.len > 0:
          var s = 0.0.cm²
          for n, x in t.series:
            if en.float <= float64(t.seriesOnset[n]): break
            if en.float > x.x[^1]:
              s = c.betheIonisation(t.sigma.fraction)
              break
            s = s + x.interpolate(en) * x.unit * t.seriesFactor[n]
          q = q + s
        elif en.float > float64(t.onset):
          q = q + t.sigma.at(en, eion, c)
      g.qion[k1, i] = q.float
      if (when m.gateInclusive: en >= ch.gate * eion else: en > ch.gate * eion):
        g.peqion[k1, i] = g.peqel[2, i - ioffion[k]]
    when m.renorm == subtractNext:
      for k in 1 ..< m.channels.len:
        g.qion[k, i] = g.qion[k, i] - g.qion[k + 1, i]
    elif m.renorm != noRenorm:
      var outer, shells = 0.0
      for k, ch in m.channels:
        if ch.shell: shells = shells + g.qion[k + 1, i]
        else: outer = outer + g.qion[k + 1, i]
      when m.renorm == counting:
        let qcount =
          if en.float <= float64(m.countingFrom): 0.0
          else: float(c.qcount)
        let fac = if qcount == 0.0: 1.0 else: (qcount - shells) / qcount
      else:
        let fac = if outer == 0.0: 1.0 else: (outer - shells) / outer
      for k, ch in m.channels:
        if not ch.shell: g.qion[k + 1, i] = g.qion[k + 1, i] * fac

type
  Bin = object
    ## what the levels of a bin depend on
    en: eV
    c: Context
    pq1: float64              ## elastic anisotropy parameter of model 1
    iso: float64
    naniso: int

proc vibrationalDipole(lv: Level, en, ein: eV): cm² =
  ## Born cross section of a dipole active vibration
  var c = float64(lv.dip)
  if lv.dipScale == cf4Dipole and en > 0.4.eV:
    c = float64((14.4'f32 - min(en.float, 5.0)) / 14.0'f32) * c
  vibrationalDipole(c, ein.float, en.float) * 1e-16.cm²

type
  Populations = object
    excited: seq[float64]    ## of each population (its ladder summed)
    ground: seq[float64]     ## of the ground state of its group

proc populations(m: static GasModel, g: GasXsec, kT: float64): Populations =
  ## Boltzmann populations of the excited vibrational states, normalised
  ## per group with the ground state
  const ps = m.populations
  for p in ps:
    var a = 0.0
    for n in 1 .. p.ladder:
      a = a + p.deg * exp(float64(float32(n)) * g.ei[p.level] / kT)
    result.excited.add a
  var sums: seq[float64]
  for k, p in ps:
    while sums.len <= p.group: sums.add 1.0
    sums[p.group] = sums[p.group] + result.excited[k]
  for k, p in ps:
    result.excited[k] = result.excited[k] / sums[p.group]
    result.ground.add 1.0 / sums[p.group]

proc level(m: static GasModel, j: int, b: Bin, ein: eV, pops: Populations,
           g: var GasXsec, i: int, ioffn: seq[int]) =
  ## cross section and anisotropy of level j (1-based) in bin i
  let lv = m.levels[j - 1]
  let en = b.en
  g.peqin[j, i] = case lv.peq
                  of peqElastic, peqScreened:
                    (if lv.peqBelow == belowZero: 0.0 else: b.iso)
                  else: 0.5
  if en <= ein or (lv.sigma.beyond == zero and en.float > lv.sigma.x[^1]):
    g.qin[j, i] = 0.0    # (closed)
    return
  var q, dip = 0.0.cm²
  case lv.kind
  of tabulated:
    q = lv.sigma.at(en, ein, b.c)
    if lv.dip != 0.0: dip = lv.vibrationalDipole(en, ein)
    if lv.ground > 0:
      let w = pops.ground[lv.ground - 1]
      q = q * w
      dip = dip * w
  of superelastic:
    let s = lv.sigma
    let e2 = en - ein
    if s.x.len > 0:
      q = if e2.float > s.x[^1]:
            case lv.tail
            of squareTail: s.y[^1] * pow2(s.x[^1] / e2.float) * s.unit
            of inverseTail: s.y[^1] * (s.x[^1] / en.float) * s.unit
            of cf4Tail: s.y[^1] * (s.x[^1] / pow2(float(en * e2))) * s.unit
          else: s.interpolate(e2) * s.unit * float(e2 / en)
    if lv.dip != 0.0: dip = lv.vibrationalDipole(en, ein)
    let p = lv.population - 1
    let w = pops.excited[p] / m.populations[p].deg
    q = q * w
    dip = dip * w
  of scaled:
    let s = m.levels[lv.level - 1].sigma
    let enp = float(en / (ein / g.ei[lv.level].eV))
    q = if enp > s.x[^1]:
          if lv.sigma.f == 0.0:
            float64(lv.factor) * s.y[^1] * power(s.x[^1] / enp, s.n) * s.unit
          else: b.c.bornExcitation(lv.sigma.f, ein, en)
        else:
          # (the ratio of the oscillator strengths in single precision)
          let factor = if lv.sigma.f == 0.0: float64(lv.factor)
                       else: float64(lv.sigma.f / s.f)
          factor * lin(s.x, s.y, enp) * s.unit
  let qin = q + dip
  g.qin[j, i] = qin.float
  case lv.peq
  of peqElastic:
    if (if lv.gateEnergy > 0.0: en.float > lv.gateEnergy else: en > lv.gate * ein):
      g.peqin[j, i] = g.peqel[2, i - ioffn[j - 1]]
  of peqScreened:
    if b.naniso == 2 and en.float >= lv.gate * abs(ein.float):
      g.peqin[j, i] = g.peqel[2, i - ioffn[j - 1]]
  of peqHalf: discard
  of peqBorn:
    let xmt = momentumFactor(en.float, ein.float) * dip +
              float64(m.momentumRatio) * q
    g.peqin[j, i] = 0.5 + float((qin - xmt) / qin)
  of peqRatio:
    g.peqin[j, i] = 0.5 + (1.0 - float64(m.momentumRatio))
  if lv.peq in {peqHalf, peqBorn, peqRatio} and en.float > lv.pq1Above:
    g.peqin[j, i] = b.pq1

proc attachment(m: static GasModel, g: var GasXsec, i: int, en: eV) =
  for k, a in m.attachment:
    let s = a.sigma
    g.qatt[k + 1, i] =
      if s.x.len > 0 and en.float > s.x[0] and en.float <= s.x[^1]:
        float(s.interpolate(en) * s.unit * float64(s.scale))
      else: 0.0
  when m.attachment.len == 1:
    g.q[4, i] = g.qatt[1, i]

proc fill*(m: static GasModel, inp: GasInput, g: var GasXsec) =
  ## the cross sections of the gas described by m (g zeroed before)
  describe(m, inp, g)
  let ioffion = inp.offsets(g.eion[1 .. g.nion])
  let ioffn = inp.offsets(g.ei[1 .. g.nin], magnitude = m.magnitudeOffsets)
  let pops = populations(m, g, inp.kT)
  var b = Bin(iso: if inp.naniso == 2: 0.0 else: 0.5, naniso: inp.naniso)
  b.c = Context(am2: float64(m.ionisation.am2), c: float64(m.ionisation.c),
                eion: float64(m.e3).eV, negative: m.negative)
  for i, e in inp.bins:
    b.en = e.eV
    b.c.r = relativistic(e)
    let (qela, qmom, pq2) = elasticAt(m.elastic, b.en)
    g.elastic(inp, i, qela.float, qmom.float, pq2)
    b.pq1 = 0.5 + float((qela - qmom) / qela)
    when m.ionisation.counting.x.len > 0:
      b.c.qcount = m.ionisation.counting.at(b.en, 0.0.eV, b.c)
    ionisation(m.ionisation, g, i, b.en, b.c, b.iso, ioffion)
    attachment(m, g, i, b.en)
    for j in 1 .. m.levels.len:
      level(m, j, b, g.ei[j].eV, pops, g, i, ioffn)
  g.dropLevelsAbove(inp.efinal, m.firstCut)
