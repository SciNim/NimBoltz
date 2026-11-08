## The cross-section record of the reference gas modules (tests/gasengine):
## GasXsec as it was before units (commit f336489), so that the frozen
## per-gas modules compile unchanged.

import nimboltz/[farrays, tables]
export farrays, tables

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

