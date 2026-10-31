## Test of the upper energy limit of the collision tables (ELIMIT, ELIMITB,
## ELIMITC and their thermal versions ELIMITT, ELIMITBT, ELIMITCT).
##
## Magboltz chooses the table range EFINAL itself when the input gives 0:
## it starts low and raises EFINAL by √2 until a short Monte Carlo run
## (NMAX/10 collisions for E only, NMAX/20 with B) never takes an electron
## into the last energy bin.

import rm48, tables, fields, collisions, freeflight

func elimitGeometry*(g: FieldGeometry): FieldGeometry =
  ## Magboltz runs ELIMIT (B = 0) for B ∥ E, which does not change the energy
  if g == fgParallel: fgElectric else: g

func sampleSize(geometry: FieldGeometry, nmax: int): int =
  nmax div (if geometry == fgElectric: 10 else: 20)

proc energyLimitExceeded*[G: static FieldGeometry](
    f: Field, tabs: MixtureTable | GasTables, target: var (ColdGas | WarmGas), estep: float64,
    start: Electron, nmax: int, rng: var MagRng, theta: var float64): bool =
  ## True as soon as an electron, starting from `start`, reaches the last
  ## energy bin. `tabs`: `MixtureTable` for `ColdGas`, `GasTables` for
  ## `WarmGas`. `theta` receives the polar angle of the last scattering
  ## (Magboltz leaves it in THETA, which MONTE then uses as the start
  ## direction).
  static: doAssert G != fgParallel, "use elimitGeometry"
  var el = start
  var nnull = 0
  for _ in 1 .. sampleSize(G, nmax):
    let h = nextHit[G](f, tabs, target, estep, el, rng, nnull)
    if h.ie == nEnergyBins: return true
    let r2 = drand48(rng)
    let i = findProcess(tableOf(tabs, h), h.ie, r2)
    let ei = energyLoss(tableOf(tabs, h), i, h.e, rng)
    let s = scatter(tableOf(tabs, h), i, h.ie, h.e, ei, h.din, rng)
    theta = s.theta
    el = leave[G](target, f, s, h)
