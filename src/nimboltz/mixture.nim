## Collision tables of a gas mixture (MIXER and MIXERT of Magboltz).
##
## The gas routines return the cross sections of one gas on the energy grid
## (`GasXsec`). From those, one `CollisionTable` is built for the mixture
## when the gas is at 0 K (MIXER), or one table per gas when thermal motion
## is simulated (MIXERT). Both cases list the same processes in the same
## order, gas by gas:
##
##   elastic, ionisation (one, or one per shell), attachment (one or
##   several), inelastic and superelastic levels.
##
## The process kind L (`iarry`) is 1..5 (elastic, ionisation, attachment,
## inelastic, superelastic); in a mixture table it is offset by 5 per gas,
## so that the Monte Carlo can count collisions per gas and kind.

import std/math
import farrays, tables, thermal

const
  maxLevels* = 250    ## inelastic levels of one gas
  maxShells* = 30     ## ionisation shells of one gas
  maxAttach* = 8      ## attachment processes of one gas

type
  GasXsec* = object
    ## Cross sections of one gas on the energy grid, as the gas routines
    ## (GASn, called through GASMIX) return them. Cross sections in cm²;
    ## an anisotropy kind (`kel`, `kin`) is 0 isotropic, 1 or 2 with the
    ## parameter in `peqel`, `peqion`, `peqin` (see `CollisionTable.index`).
    q*: Arr2[float64, 6, nEnergyBins]
      ## 1 momentum transfer, 2 elastic, 3 ionisation, 4 attachment
    e*: array[1..6, float64]
      ## 2: 2m/M, 3: ionisation energy, 4: attachment threshold (eV)
    peqel*: Arr2[float64, 6, nEnergyBins]  ## anisotropy of q[2], q[3]
    kel*: array[1..6, int]                 ## anisotropy kind of q[2], q[3]
    nion*: int                             ## > 1: ionisation by shells
    qion*, peqion*: Arr2[float64, maxShells, nEnergyBins]
    eion*: array[1..maxShells, float64]    ## ionisation energy of a shell
    eb*: array[1..maxShells, float64]      ## binding energy of a shell
    nc0*: array[1..maxShells, int]         ## Auger electrons after shell ionisation
    ec0*: array[1..maxShells, float64]     ## their total energy (eV)
    natt*: int                             ## > 1: several attachment processes
    qatt*: Arr2[float64, maxAttach, nEnergyBins]
    nin*: int                              ## inelastic and superelastic levels
    qin*, peqin*: Arr2[float64, maxLevels, nEnergyBins]
    ei*: array[1..maxLevels, float64]      ## energy loss, < 0 for superelastic
    kin*: array[1..maxLevels, int]
    penfra*: Arr2[float64, 3, maxLevels]
      ## Penning transfer: fraction, delocalisation length (µm), time (ps)
    name*: string
    scrp*: array[1..300, string]           ## process descriptions
    # returned by the gas routines but not used by the tables
    qnul*: Arr2[float64, 10, nEnergyBins]
    nul*: int
    scln*: array[1..10, float64]
    wk*, efl*, eg1*, eg2*: array[1..maxShells, float64]
    ng1*, ng2*: array[1..maxShells, int]
    scrpn*: array[1..10, string]

  EnergyGrid* = object
    estep*: float64                        ## bin width (eV)
    e*, eroot*: array[1..nEnergyBins, float64]  ## bin centres and their roots

  NegativeFrequency* = object
    ## a process with a negative collision frequency (Magboltz warns)
    gas*: int            ## gas of a per-gas table, 0 for a mixture table
    ie*, process*: int
    cf*: float64

  Builder[N: static int] = object
    ## state of one table while its processes are appended
    np: int                                ## processes so far
    ie: int                                ## energy bin
    l0: int                                ## offset of the process kind
    gas: int                               ## per-gas table: its gas, else 0

func energyGrid*(efinal: float64, nstep: int): EnergyGrid =
  result.estep = efinal / float64(nstep)
  let ehalf = result.estep / 2.0
  result.e[1] = ehalf
  for i in 2 .. nEnergyBins:
    result.e[i] = ehalf + result.estep * float64(i - 1)
    result.eroot[i] = sqrt(result.e[i])
  result.eroot[1] = sqrt(ehalf)

func angleCut*(p: float64): tuple[angc, psct: float64] =
  ## ANGCUT: for p > 1, the cut angle (in units of pi/2) of the angular
  ## distribution and the renormalised forward scattering probability
  result = (1.0, p)
  if p <= 1.0: return
  let rads = 2.0 / arccos(-1.0)
  let cns = p - 0.5
  let thetac = arcsin(2.0 * sqrt(cns - cns * cns))
  let fac = (1.0 - cos(thetac)) / (sin(thetac) * sin(thetac))
  result = (thetac * rads, cns * fac + 0.5)

func isotropic(g: GasXsec): bool =
  for k in g.kel:
    if k != 0: return false
  for j in 1 .. g.nin:
    if g.kin[j] != 0: return false
  true

func anisotropic*(gases: openArray[GasXsec]): bool =
  ## NISO: any anisotropic scattering data in the mixture
  for k in 0 ..< gases.len:   # (by index: a GasXsec is too big to copy)
    if not gases[k].isotropic: return true

proc add[N: static int](b: var Builder[N], tab: var CollisionTable[N],
                        cf: float64, kind: int, peq: float64) =
  ## next process in this energy bin, with its scattering parameters
  inc b.np
  if b.np > N:
    raise newException(ValueError, "too many collision processes (more than " &
                       $N & ")")
  tab.cf[b.ie, b.np] = cf
  tab.psct[b.ie, b.np] = 0.5
  tab.angct[b.ie, b.np] = 1.0
  tab.index[b.np] = 0
  if kind == 1:
    let (angc, psct) = angleCut(peq)
    tab.angct[b.ie, b.np] = angc
    tab.psct[b.ie, b.np] = psct
    tab.index[b.np] = 1
  if kind == 2:
    tab.psct[b.ie, b.np] = peq
    tab.index[b.np] = 2

proc describe[N: static int](b: Builder[N], tab: var CollisionTable[N],
                             descr: var array[1..N, string], kind: int,
                             rgas, ein: float64, ipn: int, text: string) =
  ## energy-independent data of the process just added
  let i = b.np
  tab.iarry[i] = b.l0 + kind
  tab.rgas[i] = rgas
  tab.ein[i] = ein
  tab.ipn[i] = ipn
  tab.penfra[1, i] = 0.0
  tab.penfra[2, i] = 0.0
  tab.penfra[3, i] = 0.0
  descr[i] = text

proc addGas[N: static int](b: var Builder[N], tab: var CollisionTable[N],
                           descr: var array[1..N, string], g: GasXsec,
                           van, efinal: float64; fcion, fcatt: var float64) =
  ## appends the processes of gas g in energy bin b.ie; van converts cross
  ## sections into collision frequencies (density · speed of 1 eV)
  let ie = b.ie
  let first = ie == 1
  let rgas = 1.0 + g.e[2] / 2.0
  b.add(tab, g.q[2, ie] * van, g.kel[2], g.peqel[2, ie])
  if first: b.describe(tab, descr, 1, rgas, 0.0, 0, g.scrp[2])
  if not (efinal < g.e[3]):
    if g.nion <= 1:
      b.add(tab, g.q[3, ie] * van, g.kel[3], g.peqel[3, ie])
      fcion = fcion + tab.cf[ie, b.np]
      if first:
        b.describe(tab, descr, 2, rgas, g.e[3] / rgas, 1, g.scrp[3])
        tab.wpl[b.np] = g.eb[1]
        tab.augerCount[b.np] = g.nc0[1]
        tab.augerEnergy[b.np] = g.ec0[1]
    else:
      for k in 1 .. g.nion:
        b.add(tab, g.qion[k, ie] * van, g.kel[3], g.peqion[k, ie])
        fcion = fcion + tab.cf[ie, b.np]
        if first:
          b.describe(tab, descr, 2, rgas, g.eion[k] / rgas, 1, g.scrp[2 + k])
          tab.wpl[b.np] = g.eb[k]
          tab.augerCount[b.np] = g.nc0[k]
          tab.augerEnergy[b.np] = g.ec0[k]
  if not (efinal < g.e[4]):
    if g.natt <= 1:
      b.add(tab, g.q[4, ie] * van, 0, 0.0)
      fcatt = fcatt + tab.cf[ie, b.np]
      if first: b.describe(tab, descr, 3, rgas, 0.0, -1, g.scrp[3 + g.nion])
    else:
      for k in 1 .. g.natt:
        b.add(tab, g.qatt[k, ie] * van, 0, 0.0)
        fcatt = fcatt + tab.cf[ie, b.np]
        if first:
          b.describe(tab, descr, 3, rgas, 0.0, -1, g.scrp[2 + g.nion + k])
  for j in 1 .. g.nin:
    b.add(tab, g.qin[j, ie] * van, g.kin[j], g.peqin[j, ie])
    if first:
      let kind = if g.ei[j] < 0.0: 5 else: 4
      b.describe(tab, descr, kind, rgas, g.ei[j] / rgas, 0,
                 g.scrp[4 + g.nion + g.natt + j])
      tab.penfra[1, b.np] = g.penfra[1, j]
      tab.penfra[2, b.np] = g.penfra[2, j] * 0.000001 / sqrt(3.0)  # µm → m, per axis
      tab.penfra[3, b.np] = g.penfra[3, j]

proc finish[N: static int](b: Builder[N], tab: var CollisionTable[N],
                           eroot: float64, negative: var seq[NegativeFrequency]) =
  ## total collision frequency of the bin, cumulative fractions of the processes
  let ie = b.ie
  tab.last = b.np
  tab.isize = 1
  while 2 * tab.isize <= tab.last: tab.isize *= 2
  tab.tcf[ie] = 0.0
  for i in 1 .. tab.last:
    tab.tcf[ie] = tab.tcf[ie] + tab.cf[ie, i]
    if tab.cf[ie, i] < 0.0:
      negative.add NegativeFrequency(gas: b.gas, ie: ie, process: i, cf: tab.cf[ie, i])
  for i in 1 .. tab.last:
    if tab.tcf[ie] == 0.0: tab.cf[ie, i] = 0.0
    else: tab.cf[ie, i] = tab.cf[ie, i] / tab.tcf[ie]
  for i in 2 .. tab.last:
    tab.cf[ie, i] = tab.cf[ie, i] + tab.cf[ie, i - 1]
  tab.tcf[ie] = tab.tcf[ie] * eroot

proc clear[N: static int](tab: var CollisionTable[N]) =
  zeroMem(addr tab, sizeof(tab))

proc buildMixture*(tab: var CollisionTable[maxMixProcesses],
                   descr: var array[1..maxMixProcesses, string],
                   gases: openArray[GasXsec], van: openArray[float64],
                   efinal: float64, grid: EnergyGrid,
                   fcion, fcatt: var array[1..nEnergyBins, float64],
                   negative: var seq[NegativeFrequency]) =
  ## MIXER: one table for the whole mixture at 0 K. fcion, fcatt: total
  ## ionisation and attachment frequencies of each bin.
  tab.clear()
  for ie in 1 .. nEnergyBins:
    fcion[ie] = 0.0
    fcatt[ie] = 0.0
    var b = Builder[maxMixProcesses](ie: ie)
    for k in 0 ..< gases.len:
      b.l0 = 5 * k
      b.addGas(tab, descr, gases[k], van[k], efinal, fcion[ie], fcatt[ie])
    b.finish(tab, grid.eroot[ie], negative)
    fcatt[ie] = fcatt[ie] * grid.eroot[ie]
    fcion[ie] = fcion[ie] * grid.eroot[ie]

proc buildPerGas*(tabs: var array[1..maxGases, CollisionTable[maxGasProcesses]],
                  descr: var array[1..maxGases, array[1..maxGasProcesses, string]],
                  gases: openArray[GasXsec], van: openArray[float64],
                  efinal: float64, grid: EnergyGrid,
                  fcion, fcatt: var array[1..nEnergyBins, float64],
                  negative: var seq[NegativeFrequency]) =
  ## MIXERT: one table per gas, for thermal motion of the gas
  for k in 1 .. gases.len: tabs[k].clear()
  for ie in 1 .. nEnergyBins:
    fcion[ie] = 0.0
    fcatt[ie] = 0.0
    for k in 0 ..< gases.len:
      var b = Builder[maxGasProcesses](ie: ie, gas: k + 1)
      b.addGas(tabs[k + 1], descr[k + 1], gases[k], van[k], efinal, fcion[ie], fcatt[ie])
      b.finish(tabs[k + 1], grid.eroot[ie], negative)
    fcatt[ie] = fcatt[ie] * grid.eroot[ie]
    fcion[ie] = fcion[ie] * grid.eroot[ie]

func nullBounds*(tab: CollisionTable, efield, tmax, small, efinal, estep,
                 const1, const3: float64): NullBounds =
  ## TCFMAX of MIXER: the largest total collision frequency an electron can
  ## reach within tmax in each of the 8 bands of starting energy
  let bp = efield * efield * const1
  let f2 = efield * const3
  let elow = min(tmax * (tmax * bp - f2 * sqrt(0.5 * efinal)) / estep - 1.0, small)
  var ehi = tmax * (tmax * bp + f2 * sqrt(0.5 * efinal)) / estep + 1.0
  if ehi > 10000.0: ehi = 10000.0
  for i in 1 .. 8:
    let jlow = max(int(float64(4000 - 500 * (9 - i) + 1) + trunc(elow)), 1)
    let jhi = min(int(float64(4000 - 500 * (8 - i)) + trunc(ehi)), nEnergyBins)
    for j in jlow .. jhi:
      if tab.tcf[j] >= result[i]: result[i] = tab.tcf[j]

func thermalGas*(tabs: openArray[CollisionTable[maxGasProcesses]],
                 gases: openArray[GasXsec], kT, echarg, emass: float64): ThermalGas =
  ## null-collision bounds and thermal speeds of each gas (end of MIXERT);
  ## kT in eV
  result.count = gases.len
  for k in 1 .. gases.len:
    for ie in 1 .. nEnergyBins:
      if tabs[k - 1].tcf[ie] >= result.tcfmax[k]: result.tcfmax[k] = tabs[k - 1].tcf[ie]
  for k in 1 .. gases.len:
    result.tcfmx = result.tcfmx + result.tcfmax[k]
  for k in 1 .. gases.len:
    result.tcfmxg[k] = result.tcfmax[k] / result.tcfmx
  for k in 2 .. gases.len:
    result.tcfmxg[k] = result.tcfmxg[k] + result.tcfmxg[k - 1]
  for k in 1 .. gases.len:
    let amgas = 2.0 * emass / gases[k - 1].e[2]
    result.vtmb[k] = sqrt(2.0 * echarg * kT / amgas) * 1e-12
