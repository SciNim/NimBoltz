## The results of a run (one deck) as values with errors (measuremancer),
## in the units Magboltz prints them in. The frame is Magboltz's: E along
## z; B along z when parallel to E, otherwise along x with E in the x-z
## plane, so that E × B points along y for B ⊥ E.

import std/[math, options, strutils]
import measuremancer
import units, vec3, tables, drift, townsend, simulation

export options, measuremancer

type
  Process* = object
    ## one collision process of a gas, in the last Monte Carlo run
    gas*: int              ## the gas: its index in the mixture (1..6)
    description*: string   ## Magboltz's
    frequency*: Measurement[ps⁻¹]

  Results* = object
    gases*: seq[string]               ## Magboltz's names of the gases
    efield*: V•cm⁻¹
    efinal*: eV                       ## upper end of the energy range, given or found
    collisions*: int                  ## real collisions of the drift run (MONTE)
    nullCollisions*: int64
    velocity*: Vector[Measurement[μm•ns⁻¹]]   ## drift velocity
    diffusion*: Tensor[Measurement[cm²•s⁻¹]]  ## diffusion tensor
    longitudinal*, transverse*: Measurement[cm²•s⁻¹]
      ## diffusion along and across E (B = 0, B ∥ E, B ⊥ E: across both E
      ## and B)
    meanEnergy*: Measurement[eV]
    ionisation*, attachment*: Measurement[cm⁻¹]
      ## Townsend coefficients α and η: from the analyses when they run
      ## (SST for B = 0, TOF otherwise), else from the drift run
    steadyState*: Option[SteadyState]       ## the analyses: SST (B = 0 only),
    pulsedTownsend*: Option[PulsedTownsend] ## PT, TOF and Friedland's estimate
    timeOfFlight*: Option[TimeOfFlight]
    friedland*: Option[Friedland]
    frequencies*: Frequencies         ## real collisions of the last run, by kind
    processes*: seq[Process]          ## and by process
    spectrum*: seq[tuple[energy: eV, fraction: float64]]
      ## energy distribution of the real collisions of the last run, per
      ## energy bin
    warnings*: seq[string]            ## what Magboltz warns about, once each

proc withError[T](x: T, percent: float64): Measurement[T] =
  ## x with an error in per cent of it
  x ± abs(x) * (percent / 100.0)

proc withErrors[T](v: Vector[T], percent: Vec3): Vector[Measurement[T]] =
  vec3(v.x.withError(percent.x), v.y.withError(percent.y), v.z.withError(percent.z))

proc withErrors[T](d: Tensor[T], percent: Tensor[float64]): Tensor[Measurement[T]] =
  Tensor[Measurement[T]](xx: d.xx.withError(percent.xx), yy: d.yy.withError(percent.yy),
                         zz: d.zz.withError(percent.zz), yz: d.yz.withError(percent.yz),
                         xy: d.xy.withError(percent.xy), xz: d.xz.withError(percent.xz))

func percentOf[T](m: Measurement[T]): float64 = 100.0 * float64(m.error / m.value)

proc characteristicEnergy*(r: Results, d: Measurement[cm²•s⁻¹]): Measurement[eV] =
  ## D/µ = D·E/v of a diffusion coefficient d (`r.longitudinal`,
  ## `r.transverse`), as Magboltz prints it for B = 0 and B ∥ E
  let v = r.velocity.z
  let e = (d.value * r.efield / v.value.to(cm•s⁻¹) * perElectron).to(eV)
  e.withError(sqrt(d.percentOf^2 + v.percentOf^2))

proc spread*(r: Results, d: Measurement[cm²•s⁻¹]): Measurement[μm•cm⁻¹ᐟ²] =
  ## √(2D/v): the spread of a diffusion coefficient d after a drift of
  ## 1 cm, as Magboltz prints it for B = 0 and B ∥ E
  let v = r.velocity.z
  let s = sqrt(2.0 * d.value / v.value.to(cm•s⁻¹)).to(μm•cm⁻¹ᐟ²)
  s.withError(sqrt(d.percentOf^2 + v.percentOf^2) / 2.0)

proc results*[M](r: Run[M]): Results =
  let t = r.transport
  for k in 1 .. r.deck.ngas: result.gases.add r.names[k].strip(leading = false)
  result.efield = r.deck.efield
  result.efinal = r.efinal
  result.collisions = r.cond.nmax
  result.nullCollisions = r.nnull
  result.velocity = t.vel.convert(μm•ns⁻¹).withErrors(t.velErr)
  result.diffusion = t.dif.withErrors(t.difErr)
  result.longitudinal = t.difln.withError(t.dfler)
  result.transverse = t.diftr.withError(t.dfter)
  result.meanEnergy = t.ave.withError(t.den)
  result.ionisation = t.alpha.withError(t.alper)
  result.attachment = t.att.withError(t.atter)
  result.steadyState = r.steadyState
  result.pulsedTownsend = r.pulsedTownsend
  result.timeOfFlight = r.timeOfFlight
  result.friedland = r.friedland
  let f = frequencies(r.last, r.deck.ngas)
  result.frequencies = f
  for j in 1 .. r.deck.ngas:
    for (descr, n) in r.processes(j):
      let rate = f.total * float64(n) / float64(f.n)
      result.processes.add Process(gas: j, description: descr.strip(leading = false),
                                   frequency: rate ± rate / sqrt(max(float64(n), 1.0)))
  for k in 1 .. nEnergyBins:
    result.spectrum.add (r.grid.e[k], r.last.spectrum[k] / float64(f.n))
  for w in r.warnings:
    if w notin result.warnings: result.warnings.add w
