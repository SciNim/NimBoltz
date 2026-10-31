## Collision tables: what MIXER/MIXERT compute from the cross sections and
## the Monte Carlo reads.
##
## One `CollisionTable` describes one collision target, binned in energy:
## the whole gas mixture when the gas is at 0 K (COMMON /LARGE/ + /ANIS/ in
## Magboltz), or a single gas of the mixture when thermal motion is
## simulated (/LARGET/ + /ANIST/, one table per gas).
##
## Tables are filled once and then only read, so all walkers of a
## simulation can share one copy (threads, or a GPU upload). They are plain
## value objects without references. The adaptive null-collision bounds
## (TCFMAX) change during the Monte Carlo and are not part of the table.

import frt

const
  nEnergyBins* = 4000      ## energy grid of all tables (NSTEP in MIXER)
  maxMixProcesses* = 960   ## collision processes of a mixture at 0 K
  maxGasProcesses* = 300   ## collision processes of one gas (thermal)
  maxGases* = 6

type
  CollisionTable*[N: static int] = object
    ## Up to N collision processes per energy bin. Index order follows
    ## Magboltz (`cf[ie, i]`: energy bin ie, process i).
    tcf*: array[1..nEnergyBins, float64]
      ## total collision frequency in each energy bin (real collisions)
    cf*: Arr2[float64, nEnergyBins, N]
      ## cumulative collision frequency of processes 1..i, divided by tcf
    psct*, angct*: Arr2[float64, nEnergyBins, N]
      ## anisotropic scattering parameters (meaning depends on `index`)
    ein*: array[1..N, float64]    ## energy loss of the process (eV)
    rgas*: array[1..N, float64]   ## mass ratio factor of the target
    wpl*: array[1..N, float64]    ## binding energy of the ionised shell (eV)
    iarry*: array[1..N, int]      ## process kind L (1 elastic, 2 ionisation, ...)
    ipn*: array[1..N, int]        ## > 0: ionisation, energy shared with the secondary
    index*: array[1..N, int]      ## 0 isotropic, 1 forward/backward cut, 2 screened
    penfra*: Arr2[float64, 3, N]  ## Penning transfer fraction, distance, time
    last*: int                    ## number of processes in use
    isize*: int                   ## start step of the binary search (power of 2)

  NullBounds* = array[1..8, float64]
    ## Upper bounds of the total collision frequency for the null-collision
    ## method at 0 K, in 8 bands of the starting energy (TCFMAX of /LARGE/).
    ## The Monte Carlo raises a bound by 5% whenever it is exceeded, so this
    ## is mutable per-walker state, not part of the shared table.
