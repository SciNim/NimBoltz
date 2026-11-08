## One routine for the cross sections of every gas: `fill` reads a
## `GasModel` (model.nim) and fills a `GasXsec` on the energy grid, as the
## gas routines GASn of Magboltz do.
##
## A gas module declares its routine with `gas name, description`; `fill`
## is instantiated per gas, and the variants it does not use (`when`) are
## not compiled in. The physics is the same as Magboltz's; the floating point operations
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
    countingBethe: bool  ## qcount is Born-Bethe (above its table)
    density: float64     ## gas density relative to 760 Torr and 0 °C
    bb: cm²•eV           ## Born-Bethe constant
    torr, tempc: float64

proc betheIonisation(c: Context, fraction = 1.0'f32): cm² =
  ## Bethe's ionisation cross section times a share of it
  let x2 = 1.0 / c.r.beta2
  let x1 = x2 * ln(c.r.beta2 / (1.0 - c.r.beta2)) - 1.0
  bornIonisation * (c.am2 * x1 + c.c * x2) * fraction

proc bornExcitation(c: Context, f: float64, ein, en: eV,
                    form = bornStandard): cm² =
  ## Born-Bethe cross section of an allowed transition with oscillator
  ## strength f and energy loss ein
  let bb = c.bb
  let r = c.r
  let b = f / (ein * r.beta2) *
    (ln(r.beta2 * r.gamma2 * emass2 / (4.0 * ein.float)) - r.beta2) * bb
  let d = case form
          of bornStandard: b * (en / (en + ein + c.eion))
          of bornShifted: b * (en + 2.0 * ein) / en
  case c.negative
  of keepNegative: d
  of clampNegative: max(d, 0.0.cm²)
  of absNegative: abs(d)

func interpolate(s: Sigma, en: eV): float64 {.inline.} =
  case s.interp
  of linear: lin(s.x, s.y, en.float)
  of logLog: logLog(s.x, s.y, en.float)
  of logLin: logLin(s.x, s.y, en.float)
  of linLog:
    if en.float <= s.switchAt: lin(s.x, s.y, en.float) else: logLog(s.x, s.y, en.float)

proc at(s: Sigma, en: eV, ein: eV, c: Context): cm² =
  ## the cross section at energy en (ein: energy loss for `born`)
  let x = en.float
  if s.zeroBelow and x <= s.x[0]: return 0.0.cm²
  if s.zeroUnder and x < s.x[0]: return 0.0.cm²
  if s.zeroAbove and x > s.x[^1]: return 0.0.cm²
  if s.logZeros:
    var j = s.x.high
    for k in 1 .. s.x.high:
      if x <= s.x[k]:
        j = k
        break
    if s.y[j] == 0.0 or s.y[j - 1] == 0.0: return 0.0.cm²
  let tableEnd = if s.upTo > 0.0: s.upTo elif s.x.len > 0: s.x[^1] else: 0.0
  if s.x.len > 0 and (s.beyond in {extrapolate, power} or x <= tableEnd):
    var v = if s.beyond == power and x > s.x[^1]:
              s.y[^1] * power(s.x[^1] / x, s.n)
            else: s.interpolate(en)
    if s.clamp and v < 0.0: v = 0.0
    return v * s.unit * (if s.scaleD != 0.0: s.scaleD else: float64(s.scale))
  case s.beyond
  of born:
    c.bornExcitation(if s.fd != 0.0: s.fd else: float64(s.f), ein, en, s.bornForm) *
      float64(s.tailScale) * float64(s.tailScale2)
  of bethe: c.betheIonisation(s.fraction) * float64(s.tailScale)
  of countingShare:
    let f = if c.countingBethe and s.betheFraction != 0.0: s.betheFraction else: s.fraction
    c.qcount * float64(f) * float64(s.tailScale)
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

proc elasticAt[M](en: eV, bend: float64):
    tuple[qela, qmom: cm², pq2: float64] =
  ## total and momentum transfer cross sections, anisotropy parameter of
  ## model 2 (bend: population of the bending states)
  mixin descriptionOf
  const m = descriptionOf(M).elastic
  let x = en.float
  template interp(kind: Interp, xs, ys: seq[float64]): float64 =
    when kind == logLog: logLog(xs, ys, x)
    elif kind == linLog:
      if x <= m.switchAt: lin(xs, ys, x) else: logLog(xs, ys, x)
    else: lin(xs, ys, x)
  template eps: float64 =
    when m.xEps.len > 0:
      when m.epsConstantBelow:
        1.0 - (if x <= m.xEps[0]: m.eps[0] else: interp(m.epsInterp, m.xEps, m.eps))
      else: 1.0 - interp(m.epsInterp, m.xEps, m.eps)
    else: 0.0
  template q(xs, ys: seq[float64]): cm² =
    when m.scale != 1.0: interp(m.interp, xs, ys) * 1e-16.cm² * m.scale
    else: interp(m.interp, xs, ys) * 1e-16.cm²
  when m.bend.populations.len > 0:
    # ground state and bending states
    template mix(q, qb: cm²): cm² =
      (1.0 - bend) * (q - m.bend.share * qb) / (1.0 - m.bend.share) + bend * qb
    if x <= m.lowBelow:
      let qela = m.lowQ.cm²
      let qmom = (if m.lowQMomentum != 0.0: m.lowQMomentum else: m.lowQ).cm²
      return (mix(qela, qela), mix(qmom, qmom), 0.0)
    return (mix(q(m.xTotal, m.total), q(m.xTotal, m.bend.total)),
            mix(q(m.xMomentum, m.momentum), q(m.xTotal, m.bend.momentum)), eps)
  else:
    when m.low == phaseShifts:
      if x <= m.lowBelow:
        let (qela, qmom) = mert(m.mert, en)
        return (qela, qmom, eps)
    elif m.low == constant:
      if x <= m.lowBelow:
        let qmom = if m.lowQMomentum != 0.0: m.lowQMomentum else: m.lowQ
        return (m.lowQ.cm² * m.scale, qmom.cm² * m.scale,
                (when m.lowIsotropic: 0.0 else: eps))
    when m.xTotal.len > 0 and m.xMomentum.len > 0:
      (q(m.xTotal, m.total), q(m.xMomentum, m.momentum), eps)
    elif m.xTotal.len > 0:
      let qela = q(m.xTotal, m.total)
      (qela, qela, eps)
    else:
      let qmom = q(m.xMomentum, m.momentum)
      (qmom, qmom, eps)

# ------------------------------------------------------------ the gas

proc describe[M](inp: GasInput, g: var GasXsec) =
  ## names, energies and descriptions of the processes
  mixin descriptionOf
  template m: untyped = descriptionOf(M)
  let aniso = inp.naniso != 0
  g.name = if aniso: m.name.anisotropic else: m.name.isotropic
  when m.unnamedModel1:
    if inp.naniso == 1: g.name = ""
  g.scrp[2] = if aniso: m.elasticText.anisotropic else: m.elasticText.isotropic
  g.massRatio = 2.0 * emass / (float64(m.mass) * amu)
  g.ionThreshold = float64(m.e3).eV
  g.attachThreshold = float64(m.e4).eV
  for k in 1 .. 6:
    g.kel[k] = when m.forceModel2: 2
               elif m.fixedModel1: (if k in 2 .. 3: 1 else: 0)
               elif m.isotropic: 0
               elif m.elastic.alwaysAnisotropic: 2
               else: inp.naniso
  const chans = m.ionisation.channels
  g.nion = chans.len
  for k, ch in chans:
    g.eion[k + 1] = float64(ch.eion).eV
    if ch.eobyScale.len > 0:
      var eb = float64(if ch.eobyBase > 0: chans[ch.eobyBase - 1].eion else: ch.eion)
      for f in ch.eobyScale: eb = eb * float64(f)
      g.eb[k + 1] = eb.eV
    else: g.eb[k + 1] = float64(ch.eoby).eV
    g.nc0[k + 1] = ch.nc0
    g.ec0[k + 1] = float64(ch.ec0).eV
    g.scrp[2 + k + 1] = ch.text
  g.natt = m.attachment.len
  for k, a in m.attachment: g.scrp[2 + g.nion + k + 1] = a.text
  g.nin = m.levels.len
  for j, lv in m.levels:
    let b0 = m.rotor.b0
    g.ei[j + 1] = eV(
      case lv.kind
      of rotation:              # E(J ± 2) - E(J), E(J) = B0 J (J + 1)
        if lv.up: b0 * float64(4 * lv.j + 6) else: -(b0 * float64(4 * lv.j - 2))
      of tabulatedRotation:     # between J and J + 2
        if lv.up: b0 * float64(4 * lv.j + 6) else: -(b0 * float64(4 * lv.j + 6))
      of dipoleRotation:
        if lv.state > 0: lv.eloss              # E(J ± 1) - E(J)
        elif lv.up: b0 * 2.0 * float64(lv.j + 1)
        else: -(b0 * 2.0 * float64(lv.j))
      elif lv.eloss != 0.0: lv.eloss
      else: lv.ein)
    g.kin[j + 1] = when m.forceModel2: 2
                   elif m.isotropic: max(lv.kin, 0)
                   else: (if lv.kin < 0: inp.naniso else: lv.kin)
    g.scrp[4 + g.nion + g.natt + j + 1] = lv.text
  let first = if m.penning.first > 0: m.penning.first else: 1
  let last = if m.penning.last > 0: m.penning.last else: g.nin
  g.penning(inp, first .. last, m.penning.fraction, m.penning.length,
            m.penning.time, m.penning.format, m.penning.all)

proc ionisation[M](g: var GasXsec, i: int, en: eV,
                c: Context, iso: float64, ioffion: seq[int]) =
  mixin descriptionOf
  const m = descriptionOf(M).ionisation
  when m.scheme == single:
    g.q[3, i] = if (when m.countingExclusive: en.float <= float64(m.countingFrom)
                    else: en.float < float64(m.countingFrom)): 0.0.cm²
                else: c.qcount
  elif m.scheme == grossCounting:
    var gross, counting = 0.0.cm²
    g.peqion[1, i] = iso
    g.peqion[2, i] = iso
    if en > g.eion[1]:
      if en.float > m.counting.x[^1]:
        counting = c.betheIonisation()
        gross = counting / float64(m.grossBethe)
      else:
        gross = m.gross.interpolate(en) * m.gross.unit
        counting = m.counting.interpolate(en) * m.counting.unit
      for k in 1 .. 2:
        if en > m.channels[k - 1].gate * g.eion[k]:
          g.peqion[k, i] = g.peqel[2, i - ioffion[k - 1]]
    g.qion[1, i] = 2.0 * counting - gross
    g.qion[2, i] = max(gross - counting, 0.0.cm²)
  else:
    var firstTerm: array[m.channels.len, cm²]  # (`Subtraction.firstTerm`)
    for k in 0 ..< m.channels.len:
      let ch {.cursor.} = m.channels[k]
      let k1 = k + 1
      g.qion[k1, i] = 0.0.cm²
      g.peqion[k1, i] = when m.peqZero: 0.0 else: iso
      let eion = ch.eion.float64.eV
      when m.fromTables:
        if en.float <= ch.terms[0].sigma.x[0]:
          when m.stopAtClosed: break
          else: continue
      else:
        if (if ch.openInclusive: en < eion else: en <= eion):
          when m.stopAtClosed: break
          else: continue
      var q = 0.0.cm²
      for n in 0 ..< ch.terms.len:
        let t {.cursor.} = ch.terms[n]
        if n == 1: firstTerm[k] = q
        if t.series.len > 0:
          var s = 0.0.cm²
          for n in 0 ..< t.series.len:
            let x {.cursor.} = t.series[n]
            if en.float <= float64(t.seriesOnset[n]): break
            if en.float > x.x[^1]:
              s = c.betheIonisation(t.sigma.fraction)
              break
            s = s + x.interpolate(en) * x.unit * t.seriesFactor[n]
          q = q + s
        elif t.shiftFrom > 0:
          let e0 = en - (eion - m.channels[t.shiftFrom - 1].eion.float64.eV)
          q = q + (if en.float > t.sigma.x[^1]: t.sigma.at(en, eion, c)
                   else: t.sigma.interpolate(e0) * t.sigma.unit * float64(t.sigma.scale))
        elif t.rampTo > 0.0:
          q = q + (if en.float < t.rampTo:
                     (float((en - eion) / (t.rampTo.eV - eion)) * float64(t.rampQ) * 1e-16).cm²
                   else: g.qion[t.shareOf, i] * float64(t.share))
        elif en.float > float64(t.onset):
          q = q + t.sigma.at(en, eion, c)
      if ch.terms.len == 1: firstTerm[k] = q
      g.qion[k1, i] = q
      if (if m.gateInclusive or ch.gateInclusive: en >= ch.gate * eion
          else: en > ch.gate * eion):
        g.peqion[k1, i] = g.peqel[2, i - ioffion[k]]
    when m.split.fractions.len > 0:
      let q = if en > g.eion[1]: m.split.sigma.at(en, g.eion[1], c)
              else: 0.0.cm²
      var open = 1
      for k in 2 .. m.split.fractions.len + 1:
        if en > g.eion[k]: open = k
      g.qion[1, i] = q
      if open > 1:
        for k, f in m.split.fractions[open - 2]: g.qion[k + 1, i] = q * float64(f)
    when m.renorm == subtractNext:
      for k in 1 ..< m.channels.len:
        g.qion[k, i] = g.qion[k, i] - g.qion[k + 1, i]

    elif m.renorm != noRenorm:
      var outer, shells = 0.0.cm²
      for k in 0 ..< m.channels.len:
        if m.channels[k].shell: shells = shells + g.qion[k + 1, i]
        else: outer = outer + g.qion[k + 1, i]
      when m.renorm == counting:
        let qcount =
          if en.float <= float64(m.countingFrom): 0.0.cm²
          else: c.qcount
        let fac = if qcount.float == 0.0: 1.0 else: float64((qcount - shells) / qcount)
      else:
        let fac = if outer.float == 0.0: 1.0 else: float64((outer - shells) / outer)
      for k in 0 ..< m.channels.len:
        if not m.channels[k].shell: g.qion[k + 1, i] = g.qion[k + 1, i] * fac
    for d in m.subtract:
      if en.float <= d.below:
        let v = if d.firstTerm: firstTerm[d.channel - 1] else: g.qion[d.channel, i]
        g.qion[d.target, i] = g.qion[d.target, i] - d.factor * v

type
  Bin = object
    ## what the levels of a bin depend on
    en: eV
    c: Context
    pq1: float64              ## elastic anisotropy parameter of model 1
    dimer: float64            ## density of `GasModel.dimer`
    qmom: float64             ## momentum transfer cross section (cm²)
    iso: float64
    naniso: int

func inside(w: Sigma, x: float64, open, lowOpen: bool): bool {.inline.} =
  ## x inside the table's range (ends excluded if open, the lower one if
  ## lowOpen)
  if open: x > w.x[0] and x < w.x[^1]
  elif lowOpen: x > w.x[0] and x <= w.x[^1]
  else: x >= w.x[0] and x <= w.x[^1]

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
    rotation: seq[float64]   ## of the rotational states J = 0, 2, ..
    rotationEnergy: seq[float64]  ## (`Rotor.ladder`)
    bend: float64            ## of the bending states (Elastic.bend)

proc populations[M](g: GasXsec, kT: eV): Populations =
  ## Boltzmann populations of the excited vibrational states, normalised
  ## per group with the ground state
  mixin descriptionOf
  template m: untyped = descriptionOf(M)
  const ps = m.populations
  for p in ps:
    var a = 0.0
    for n in 1 .. p.ladder:
      a = a + p.deg * exp(float64(float32(n)) * g.ei[p.level] / kT)
    result.excited.add a
  var sums: seq[float64]
  for k, p in ps:
    while sums.len <= p.group: sums.add 1.0
    if not p.unnormalised: sums[p.group] = sums[p.group] + result.excited[k]
  for k, p in ps:
    if not p.unnormalised: result.excited[k] = result.excited[k] / sums[p.group]
    result.ground.add 1.0 / sums[p.group]
  for k in m.elastic.bend.populations: result.bend = result.bend + result.excited[k - 1]
  when m.rotor.ladder.len > 0:
    const r = m.rotor
    var e = newSeq[float64](r.jMax + 1)
    e[1] = r.oddBase
    for j in 2 .. r.jMax: e[j] = e[j - 2] + g.ei[r.ladder[j - 2]].float
    var p = newSeq[float64](r.jMax + 1)
    for j in 0 .. r.jMax: p[j] = r.weights[j mod 2] * float64(2 * j + 1) * exp(-e[j] / kT.float)
    var sum = 0.0
    if r.fixBelow > 0.0 and kT.float < r.fixBelow:
      var odd, even = 0.0
      for j in countup(1, r.jMax, 2): odd = odd + p[j]
      for j in countup(0, r.jMax, 2): even = even + p[j]
      let fix = odd * r.fixRatio.mul / r.fixRatio.den
      let orpafix = fix / even
      for j in countup(2, r.jMax, 2): p[j] = p[j] * orpafix
      sum = odd + fix
      p[0] = p[0] * orpafix
    else:
      for v in p: sum = sum + v
    for v in p: result.rotation.add v / sum
    result.rotationEnergy = e
  elif m.rotor.jMax > 0 or m.rotor.states.len > 0:
    var sum = 0.0
    when m.rotor.states.len > 0:
      for (e, j, w) in m.rotor.states:
        result.rotation.add w * (2.0 * j + 1.0) * exp(-(e / kT.float))
        sum = sum + result.rotation[^1]
    else:
      for j in countup(m.rotor.jStart, m.rotor.jMax, m.rotor.jStep):
        result.rotation.add m.rotor.weights[j mod 2] * float64(2 * j + 1) *
          exp(-(float64(j * (j + 1)) * m.rotor.b0 / kT.float))
        sum = sum + result.rotation[^1]
    for p in result.rotation.mitems: p = p / sum

proc level[M](j: int, b: Bin, ein: eV, pops: Populations,
           g: var GasXsec, i: int, ioffn: seq[int]) =
  ## cross section and anisotropy of level j (1-based) in bin i
  mixin descriptionOf
  template m: untyped = descriptionOf(M)
  let lv {.cursor.} = m.levels[j - 1]
  let en = b.en
  when m.isotropic:
    if lv.kin > 0: g.peqin[j, i] = 0.5
  else:
    g.peqin[j, i] = case lv.peq
                    of peqElastic, peqScreened, peqDipole, peqMomentumRatio, peqSameAs:
                      (if lv.peqBelow == belowZero: 0.0 else: b.iso)
                    else: 0.5
  if (if lv.openAbove > 0.0: en.float <= lv.openAbove else: en <= ein) or
     (lv.sigma.beyond == zero and en.float > lv.sigma.x[^1]) or
     (lv.kind in {dipoleRotation, symmetricTop} and en.float > m.rotor.dipoleCutoff) or
     en.float > lv.closedAbove or
     (lv.kind == tabulatedRotation and lv.up and en.float < lv.sigma.x[0]):
    g.qin[j, i] = 0.0.cm²    # (closed)
    return
  var q, dip = 0.0.cm²
  case lv.kind
  of tabulated:
    let s {.cursor.} = lv.sigma
    when m.stepped:
      # (the old routines add table and dipole term before the unit)
      const fold = true
    else:
      const fold = false
    if lv.windows.len > 0:
      # (dipole term and windows summed before the unit)
      var t = vibrationalDipole(float64(lv.dip), ein.float, en.float) * lv.dipFactor
      for w in lv.windows:
        if inside(w, en.float, lv.windowOpen, lv.windowLowOpen):
          t = t + (if w.scaleD != 0.0: w.scaleD else: float64(w.scale)) * w.interpolate(en)
      q = s.at(en, ein, b.c) + t * lv.windows[0].unit
    elif fold and lv.dip != 0.0 and lv.peq != peqBorn and s.x.len > 0 and
       s.beyond == extrapolate:
      let t = if s.zeroBelow and en.float <= s.x[0]: 0.0
              else: float64(s.scale) * s.interpolate(en)
      q = (vibrationalDipole(float64(lv.dip), ein.float, en.float) + t) * s.unit
    else:
      q = s.at(en, ein, b.c)
      if lv.dip != 0.0: dip = lv.vibrationalDipole(en, ein)
    if lv.extra.x.len > 0: q = q + lv.extra.at(en, ein, b.c)
    if lv.ground > 0:
      let w = pops.ground[lv.ground - 1]
      q = q * w
      dip = dip * w
  of superelastic:
    let s {.cursor.} = lv.sigma
    let e2 = en - ein
    if s.x.len > 0 and not (s.zeroBelow and e2.float <= s.x[0]):
      q = if e2.float > s.x[^1]:
            case lv.tail
            of squareTail: s.y[^1] * pow2(s.x[^1] / e2.float) * s.unit
            of inverseTail: s.y[^1] * (s.x[^1] / en.float) * s.unit
            of crossTail: s.y[^1] * s.x[^1] * float(e2 / (en * en).float) * s.unit
            of mixedTail: s.y[^1] * (s.x[^1] / e2.float) / en.float * s.unit
            of e2Tail: s.y[^1] * s.x[^1] / e2.float * s.unit
            of powerTail: s.y[^1] * s.unit * pow(s.x[^1] / e2.float, s.n)
            of zeroTail: 0.0.cm²
            of cf4Tail: s.y[^1] * (s.x[^1] / pow2(float(en * e2))) * s.unit
            of noTail: float(e2) * s.interpolate(e2) / en.float * s.unit
          else: float(e2) * s.interpolate(e2) / en.float * s.unit
      q = q * float64(s.scale)
    if lv.dip != 0.0: dip = lv.vibrationalDipole(en, ein)
    if lv.windows.len > 0:
      var t = vibrationalDipole(float64(lv.dip), ein.float, en.float) * lv.dipFactor
      for w in lv.windows:
        if inside(w, e2.float, lv.windowOpen, lv.windowLowOpen):
          t = t + (if w.scaleD != 0.0: w.scaleD else: float64(w.scale)) * float(e2) *
                  w.interpolate(e2) / en.float
      q = q + t * lv.windows[0].unit
      dip = 0.0.cm²
    when m.stepped:
      # (the old routines add the dipole term before the unit)
      if dip != 0.0.cm² and lv.tail == noTail:
        var t = 0.0
        if s.x.len > 0 and not (s.zeroBelow and e2.float <= s.x[0]):
          t = float64(s.scale) * (float(e2) * s.interpolate(e2) / en.float)
        q = (vibrationalDipole(float64(lv.dip), ein.float, en.float) + t) * 1e-16.cm²
        dip = 0.0.cm²
    if lv.population > 0:
      let p = lv.population - 1
      let w = pops.excited[p] / m.populations[p].deg
      q = q * w
      dip = dip * w
  of rotation:
    let r {.cursor.} = m.rotor
    let aj = float64(lv.j)
    let p = pops.rotation[(lv.j - r.jStart) div r.jStep]
    let born = p * r.quadrupole.cm² * sqrt(1.0 - ein / en)
    q = if lv.up: born * (aj + 2.0) * (aj + 1.0) / ((2.0 * aj + 3.0) * (2.0 * aj + 1.0))
        else: born * aj * (aj - 1.0) / ((2.0 * aj + 1.0) * (2.0 * aj - 1.0))
    if r.scaleAbove > 0.0 and en.float >= r.scaleAbove:
      q = q * (r.scaleAbove / en.float)
    when m.rotor.enhance.x.len > 0:
      let y = if lv.up: r.enhance.interpolate(en)
              else: float(en - ein) * r.enhance.interpolate(en - ein) / en.float
      q = q * (1.0 + y * float64(r.enhanceNum) / float64(r.enhanceDen))
    when m.rotor.momentumAbove > 0.0:
      if en.float >= r.momentumAbove: q = q * (b.qmom / r.momentumRef)
  of dipoleRotation:
    let r {.cursor.} = m.rotor
    let (p, aj) =
      if lv.state > 0: (pops.rotation[lv.state - 1], lv.aj)
      else: (pops.rotation[(lv.j - r.jStart) div r.jStep], float64(lv.j))
    let strength = if lv.strength != 0.0: lv.strength
                   elif lv.up: aj + 1.0
                   else: aj
    let enrt = sqrt(en.float)
    let s = sqrt(float(en - ein))
    let lg = if lv.up: ln((enrt + s) / (enrt - s)) else: ln((enrt + s) / (s - enrt))
    # operation order as in Magboltz: a tabulated line strength (H2O) goes
    # with the constant, the J factor of a linear rotor at the end
    let num = if lv.strength != 0.0: r.dipole.cm²•eV * strength * p * lg
              else: p * r.dipole.cm²•eV * lg * strength
    q = num / ((2.0 * aj + 1.0) * en) * float64(r.dipoleScale)
  of symmetricTop:
    let r {.cursor.} = m.rotor
    let enrt = sqrt(en.float)
    let s = sqrt(float(en - ein))
    let lg = if ein.float < 0.0: ln((enrt + s) / (s - enrt)) else: ln((enrt + s) / (enrt - s))
    let ajg2 = lv.ajg * lv.ajg
    let den = lv.ajg * (2.0 * lv.aj + 1.0) * en.float
    q = if lv.bracket:
          var w = 0.0
          for (p, k2) in lv.top: w = w + pops.rotation[p - 1] * (ajg2 - k2)
          (r.dipole * lg * w / den).cm²
        else: (r.dipole * lg * pops.rotation[lv.top[0].state - 1] * ajg2 / den).cm²
  of tabulatedRotation:
    let r {.cursor.} = m.rotor
    let s {.cursor.} = lv.sigma
    let x = if lv.up: en else: en - ein    # energy of the excitation
    let jl = float64(lv.j)
    let (p, f) =
      if lv.up: (pops.rotation[(lv.j - r.jStart) div r.jStep], 1.0)
      else: (pops.rotation[(lv.j + 2 - r.jStart) div r.jStep],
             (2.0 * jl + 1.0) / (2.0 * jl + 5.0))
    q = if x.float > s.x[^1]: p * f * s.y[^1] * s.x[^1] / x.float * s.unit
        elif lv.up: p * s.interpolate(x) * s.unit
        else: p * f * float(x / en) * s.interpolate(x) * s.unit
    let rs {.cursor.} = r.resonance
    if x.float <= rs.x[^1] and x.float > rs.x[0]:
      q = q + (if lv.up: rs.interpolate(x) * p
               else: float(x / en) * rs.interpolate(x) * p) * rs.unit
    if en.float >= r.scaleAbove: q = q * (r.scaleAbove / en.float)
    q = q * r.tabulatedScale
  of scaled:
    let s = m.levels[lv.level - 1].sigma
    let enp = float(en / (ein / g.ei[lv.level]))
    q = if enp > s.x[^1]:
          if lv.sigma.f == 0.0:
            float64(lv.factor) * s.y[^1] * power(s.x[^1] / enp, s.n) * s.unit
          else: b.c.bornExcitation(float64(lv.sigma.f), ein, en)
        else:
          # (the ratio of the oscillator strengths in single precision)
          let factor = if lv.sigma.f == 0.0: float64(lv.factor)
                       else: float64(lv.sigma.f / s.f)
          factor * lin(s.x, s.y, enp) * s.unit
  if lv.clampTotal and q < 0.0.cm²: q = 0.0.cm²
  if lv.dampAbove > 0.0 and en.float > lv.dampAbove:
    q = q * sqrt(lv.dampAbove / en.float)
    if q < 0.0.cm²: q = 0.0.cm²
  if lv.rotState >= 0:
    let w = pops.rotation[lv.rotState] * lv.weight
    q = q * w
    dip = dip * w
  if lv.dimer:
    q = q * b.dimer
    dip = dip * b.dimer
  let qin = q + dip
  g.qin[j, i] = qin
  when m.isotropic:
    if lv.kin <= 0: return
  case lv.peq
  of peqElastic:
    let e0 = if lv.gateAbs: abs(ein.float) else: ein.float
    if (if lv.gateEnergy > 0.0:
          (if lv.gateInclusive: en.float >= lv.gateEnergy else: en.float > lv.gateEnergy)
        elif lv.gateInclusive: en.float >= lv.gate * e0
        else: en.float > lv.gate * e0):
      g.peqin[j, i] = g.peqel[2, i - ioffn[j - 1]]
  of peqScreened:
    if b.naniso == 2 and en.float >= lv.gate * abs(ein.float):
      g.peqin[j, i] = g.peqel[2, i - ioffn[j - 1]]
  of peqDipole:
    if lv.elasticAbove > 0.0:
      let e0 = abs(ein.float)
      if en.float > lv.elasticAbove * e0:
        g.peqin[j, i] = g.peqel[2, i - ioffn[j - 1]]
      elif en.float < lv.dipoleUpTo * e0:
        g.peqin[j, i] = 1.0 - lin(m.dipoleAngles.xUp, m.dipoleAngles.eps, en.float / e0)
    elif lv.gateEnergy <= 0.0 or en.float >= lv.gateEnergy:
      let ratio = en.float / abs(ein.float)
      g.peqin[j, i] =
        if ratio > lv.dipoleUpTo: g.peqel[2, i - ioffn[j - 1]]
        elif ein.float < 0.0: 1.0 - lin(m.dipoleAngles.xDown, m.dipoleAngles.eps, ratio)
        else: 1.0 - lin(m.dipoleAngles.xUp, m.dipoleAngles.eps, ratio)
  of peqMomentumRatio:
    if en.float > lv.gate * ein.float:
      if b.naniso == 1:
        let x {.cursor.} = lv.sigma.x
        let ratio = if en.float > x[^1]: lv.ratio[^1] * (x[^1] / en.float)
                    else: lin(x, lv.ratio, en.float)
        g.peqin[j, i] = 1.5 - ratio
      elif b.naniso == 2: g.peqin[j, i] = g.peqel[2, i - ioffn[j - 1]]
  of peqSameAs: g.peqin[j, i] = g.peqin[lv.level, i]
  of peqFixed: g.peqin[j, i] = lv.peqValue
  of peqHalf: discard
  of peqBorn:
    let xmt = momentumFactor(en.float, ein.float) * dip +
              float64(m.momentumRatio) * q
    g.peqin[j, i] = 0.5 + float((qin - xmt) / qin)
  of peqRatio:
    g.peqin[j, i] = 0.5 + (1.0 - float64(m.momentumRatio))
  if lv.peq in {peqHalf, peqBorn, peqRatio} and en.float > lv.pq1Above:
    g.peqin[j, i] = b.pq1
  if en.float > lv.peqUpTo: g.peqin[j, i] = 0.0

proc h2Attachment(a: Attach, en: eV, c: Context, pops: Populations): float64 =
  let h = a.h2
  let s {.cursor.} = a.sigma
  let x = en.float
  let e = pops.rotationEnergy
  let n = h.strengths.len
  block resonances:
    for k in countdown(n, 1):
      if x < h.threshold - e[k]: break resonances
      let t = h.amp * h.strengths[k - 1] * pops.rotation[k] *
              exp(-((x - h.threshold + e[k]) / h.width))
      result = if k == n: t else: result + t
    if x < h.threshold: break resonances
    result = result + h.amp * pops.rotation[0] * exp(-((x - h.threshold) / h.width))
    if x < s.x[0]: break resonances
    if x <= s.x[^1]: result = result + float(s.at(en, 0.0.eV, c))
    if x <= h.highAt: break resonances
    result = result + h.highAmp * exp(-((x - h.highAt) / h.highWidth))

proc attachment[M](g: var GasXsec, i: int, en: eV, c: Context,
                pops: Populations, iso: float64) =
  mixin descriptionOf
  template m: untyped = descriptionOf(M)
  for k in 0 ..< m.attachment.len:
    let a {.cursor.} = m.attachment[k]
    let s {.cursor.} = a.sigma
    let x = en.float
    if a.isoPeq: g.peqel[4, i] = iso
    g.qatt[k + 1, i] = cm²(
      case a.formula
      of h2Attachment: h2Attachment(a, en, c, pops)
      of cf3brAttachment:
        if x > 2.0: 0.0
        else:
          let at1 = float64(0.0353'f32) / sqrt(x)
          let at2 = exp(-x * x / 9e-6)
          let at3 = exp(-x / float64(0.0588'f32))
          (at1 * at2 + at3) * 27.4e-16
      of attachTable:
        if s.x.len > 0 and
           (a.openBelow or (if a.fromStart: x >= s.x[0] else: x > s.x[0])) and
           (s.beyond != zero or (if a.beforeEnd: x < s.x[^1] else: x <= s.x[^1])):
          var v = float(s.at(en, 0.0.eV, c)) *
            (if not a.threeBody: 1.0
             elif a.densityRef > 0.0:
               a.densityRef * c.torr / ((c.tempc + 273.1499938964844) * 760.0)
             else: c.density)
          if a.saturation > 0.0: v = v * (1.0 / (a.saturation / c.torr + 1.0))
          for f in a.factors: v = v * float64(f)
          v
        else: 0.0)
  when m.attachment.len == 1:
    g.q[4, i] = g.qatt[1, i]

proc fill*[M](inp: GasInput, g: var GasXsec) =
  ## the cross sections of the gas described by m (g zeroed before)
  mixin descriptionOf
  template m: untyped = descriptionOf(M)
  var inp = inp
  when m.alwaysModel2: inp.naniso = 2
  when m.forceModel2:
    if inp.naniso != 0: inp.naniso = 2
  describe[M](inp, g)
  when m.roundedOffsets:
    # (int(float32(0.5 + e / ESTEP)), as these routines)
    var ioffion = @[int(float32(0.5 + float64(g.ionThreshold / inp.estep)))]
    var ioffn: seq[int]
    for k in 1 .. g.nin: ioffn.add int(float32(0.5 + float64(g.ei[k] / inp.estep)))
  else:
    let ioffion = inp.offsets(g.eion[1 .. g.nion])
    let ioffn = inp.offsets(g.ei[1 .. g.nin], magnitude = m.magnitudeOffsets)
  let pops = populations[M](g, inp.kT)
  var b = Bin(iso: if inp.naniso == 2: 0.0 else: 0.5, naniso: inp.naniso)
  b.c = Context(am2: float64(m.ionisation.am2), c: float64(m.ionisation.c),
                eion: float64(m.e3).eV, negative: m.negative,
                bb: bbconst.cm²•eV)
  when m.bbBohr > 0.0:
    const (a0, ry) = (m.bbBohr, m.bbRydberg)
    b.c.bb = (16.0 * arccos(-1.0) * a0 * a0 * ry * ry / emass2).cm²•eV
  when m.dimer.gas > 0:
    b.dimer = inp.percent(m.dimer.gas) / 100.0 *
      (float64(m.dimer.temperature) / (273.1499938964844 + inp.tempc)) * inp.torr / 760.0.Torr
    b.dimer = b.dimer * float64(m.dimer.fraction)
  b.c.torr = inp.torr.float
  b.c.tempc = inp.tempc
  # (273.15 single precision in Magboltz)
  b.c.density = 273.1499938964844 * inp.torr / ((inp.tempc + 273.1499938964844) * 760.0.Torr)
  iterator energies(inp: GasInput): (int, eV) =
    when m.stepped:
      for i, e in inp.steps: yield (i, e)
    else:
      for i, e in inp.bins: yield (i, e)
  for i, e in energies(inp):
    b.en = e
    b.c.r = relativistic(e.float)
    let (qela, qmom, pq2) = elasticAt[M](b.en, pops.bend)
    when m.forceModel2:
      if inp.naniso == 0:
        g.q[2, i] = qmom
        g.peqel[2, i] = 0.0
      else:
        g.q[2, i] = qela
        g.peqel[2, i] = pq2
    elif m.fixedModel1:
      g.q[2, i] = qela
      g.peqel[2, i] = 0.5 + float((qela - qmom) / qela)
    elif m.isotropic:
      g.q[2, i] = qmom
      when m.isotropicPeq != 0.0:
        g.peqel[2, i] = m.isotropicPeq
        g.peqel[3, i] = m.isotropicPeq
    elif m.elastic.alwaysAnisotropic:
      g.q[2, i] = qela
      g.peqel[2, i] = pq2
    else: g.elastic(inp, i, qela, qmom, pq2)
    when m.elastic.isotropicBelow > 0.0:
      if inp.naniso > 0 and e.float < m.elastic.isotropicBelow:
        g.q[2, i] = qmom
        g.peqel[2, i] = b.iso
    when not m.isotropic: b.pq1 = 0.5 + float((qela - qmom) / qela)
    b.qmom = qmom.float
    when m.ionisation.counting.x.len > 0:
      b.c.qcount = m.ionisation.counting.at(b.en, 0.0.eV, b.c)
      b.c.countingBethe = e.float > m.ionisation.counting.x[^1]
    ionisation[M](g, i, b.en, b.c, b.iso, ioffion)
    when m.fixedModel1:
      g.peqel[3, i] = 0.5
      if e > 2.0 * g.ionThreshold: g.peqel[3, i] = g.peqel[2, i - ioffion[0]]
    attachment[M](g, i, b.en, b.c, pops, b.iso)
    for j in 1 .. m.levels.len:
      level[M](j, b, g.ei[j], pops, g, i, ioffn)
    when m.elastic.minusLevels.b > 0:
      when m.elastic.minusEach:
        for j in m.elastic.minusLevels: g.q[2, i] = g.q[2, i] - g.qin[j, i]
      else:
        var sum = 0.0.cm²
        for j in m.elastic.minusLevels: sum = sum + g.qin[j, i]
        g.q[2, i] = g.q[2, i] - sum
      when m.elastic.minusFloor > 0.0:
        if g.q[2, i] <= 0.0.cm²: g.q[2, i] = m.elastic.minusFloor.cm²
  when m.appended.n > 0:
    g.nin = m.levels.len - m.appended.n
  when m.ninBy.len > 0:
    for (e, n) in m.ninBy:
      if inp.efinal.float > e:
        g.nin = n
        return
  when m.keepLastAt:
    if inp.efinal.float == g.ei[g.nin].float: return
  when m.cutStrict:
    for k in m.firstCut .. g.nin:
      if inp.efinal < g.ei[k]:
        g.nin = k - 1
        break
  else: g.dropLevelsAbove(inp.efinal, m.firstCut)
  when m.minLevels > 0:
    if g.nin < m.minLevels: g.nin = m.minLevels
  when m.appended.n > 0:
    if inp.efinal.float > m.appended.above: g.nin = g.nin + m.appended.n

template gas*(name: untyped, description: GasModel) =
  ## Declares the gas routine `name` for the constant `description`. The
  ## engine reaches the description through a type of its own (`Tag`), not
  ## as a static parameter: a static value is copied into every use and
  ## every nested instantiation, which made the compiler about 20 times
  ## slower.
  type Tag = object
  template descriptionOf(_: type Tag): GasModel = description
  proc name*(inp: GasInput, g: var GasXsec) = fill[Tag](inp, g)
