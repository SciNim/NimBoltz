## ELIMIT, ELIMITB, ELIMITC, ELIMITT, ELIMITBT, ELIMITCT of magboltz-11.21.f:
## bridge from the COMMON blocks to nimboltz/elimit.

import ../frt, commons
import ../[fields, collisions, thermal, freeflight, elimit]

proc run[G: static FieldGeometry](thermalMotion: bool): bool =
  let f = initField(G, setp.efield, cnsts1.const1, cnsts1.const2, cnsts1.const3,
                    cnsts1.const5, bfld.wb, bfld.eovb, bfld.btheta)
  let start = initElectron[G](f, setp.estart, setp.theta, setp.phi)
  setp.small = small
  if thermalMotion:
    let gas = ThermalGas(tcfmx: thrm.tcfmx, tcfmxg: thrm.tcfmxg,
                         tcfmax: larget.tcfmax, vtmb: thrm.vtmb, count: inpt.ngas)
    var target = initWarmGas(gas, rng)
    result = energyLimitExceeded[G](f, gasTable, target, inpt.estep, start,
                                    setp.nmax, rng, setp.theta)
  else:
    var target = initColdGas(large.tcfmax, inpt.efinal)
    result = energyLimitExceeded[G](f, mixTable, target, inpt.estep, start,
                                    setp.nmax, rng, setp.theta)
    large.tcfmax = target.bounds

proc elimit*(thermalMotion: bool, ielow: var int) =
  ## IELOW = 1 if the energy limit EFINAL was exceeded, else 0
  withGeometry(geometryOf(bfld.bmag, bfld.btheta), G):
    ielow = ord(run[elimitGeometry(G)](thermalMotion))
