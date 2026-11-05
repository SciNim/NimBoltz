## MIXER and MIXERT of magboltz-11.21.f: bridge from the COMMON blocks to
## nimboltz/mixture.

import ../frt, commons, gasmix
import ../mixture, ../gases

const
  fmt776 = "(\'  WARNING NEGATIVE COLLISION FEQUENCY =\',D12.3,\' IE =\',I6,\' IF =\',I3,\' IARRY=\',I5,\' EIN=\',F7.4)"
  fmt776t = "(\'  WARNING NEGATIVE COLLISION FEQUENCY =\',D12.3,\' KGAS =\', I2,\' IE =\',I6,\' IF =\',I3,\' IARRY=\',I5,\' EIN=\',F7.4)"
  fmt992 = "(/,/,6X,\'WARNING TOO MANY LEVELS IN CALCULATION. CAN INCREASE THE ARRAY SIZES FROM 960 UP TO 1800 MAXIMUM\',/)"

proc crossSections(): (EnergyGrid, seq[GasXsec], seq[float64]) =
  ## energy grid (/MIX2/), cross sections of each gas and their VAN factors
  let grid = energyGrid(inpt.efinal, inpt.nstep)
  inpt.estep = grid.estep
  mix2.eg = grid.e
  mix2.eroot = grid.eroot
  var gases = newSeq[GasXsec](inpt.ngas)
  let inp = GasInput(e: grid.e, naniso: inpt.naniso, kT: inpt.akt,
                     efinal: inpt.efinal, penning: inpt.ipen != 0)
  for k, g in gases.mpairs:
    if not crossSections(gasn.ngasn[k + 1], inp, g):
      gasmix(gasn.ngasn[k + 1], g.q, g.qin, g.nin, g.e, g.ei, g.name, 0.0, g.eb,
             g.peqel, g.peqin, g.penfra, g.kel, g.kin, g.qion, g.peqion, g.eion,
             g.nion, g.qatt, g.natt, g.qnul, g.nul, g.scln, g.nc0, g.ec0, g.wk,
             g.efl, g.ng1, g.eg1, g.ng2, g.eg2, g.scrp, g.scrpn)
  let van = [mratio.van1, mratio.van2, mratio.van3, mratio.van4, mratio.van5,
             mratio.van6]
  (grid, gases, van[0 ..< inpt.ngas])

proc mixer*() =
  ## collision table of the mixture at 0 K
  let (grid, gases, van) = crossSections()
  var descr: array[1..maxMixProcesses, string]
  var negative: seq[NegativeFrequency]
  try:
    buildMixture(mixTable, descr, gases, van, inpt.efinal, grid, fred.fcion,
                 fred.fcatt, negative)
  except ValueError:
    fwrite(fmt992)
    fstop()
  for n in negative:
    fwrite(fmt776, n.cf, n.ie, n.process, mixTable.iarry[n.process],
           mixTable.ein[n.process])
  for j in 1 .. maxGases:
    names.nameg[j] = fstr("-------------------------", 25)
  for k in 0 ..< gases.len:
    names.nameg[k + 1] = fstr(gases[k].name, 25)
  for i in 1 .. mixTable.last:
    scrip.dscrpt[i] = fstr(descr[i], 50)
  anis.niso = ord(anisotropic(gases))
  large.tcfmax = nullBounds(mixTable, setp.efield, setp.tmax, setp.small,
                            inpt.efinal, inpt.estep, cnsts1.const1, cnsts1.const3)

proc mixert*() =
  ## one collision table per gas, for thermal motion of the gas
  let (grid, gases, van) = crossSections()
  var descr: array[1..maxGases, array[1..maxGasProcesses, string]]
  var negative: seq[NegativeFrequency]
  buildPerGas(gasTable, descr, gases, van, inpt.efinal, grid, fred.fcion,
              fred.fcatt, negative)
  for n in negative:
    fwrite(fmt776t, n.cf, n.gas, n.ie, n.process, gasTable[n.gas].iarry[n.process],
           gasTable[n.gas].ein[n.process])
  for j in 1 .. maxGases:
    names.nameg[j] = fstr("---------------", 25)
  for k in 0 ..< gases.len:
    names.nameg[k + 1] = fstr(gases[k].name, 25)
    for i in 1 .. gasTable[k + 1].last:
      script.dscrpt[k + 1, i] = fstr(descr[k + 1][i], 50)
  anist.niso = ord(anisotropic(gases))
  let tg = thermalGas(gasTable, gases, inpt.akt, cnsts.echarg, cnsts.emass)
  larget.tcfmax = tg.tcfmax
  thrm.tcfmx = tg.tcfmx
  thrm.tcfmxg = tg.tcfmxg
  thrm.vtmb = tg.vtmb
