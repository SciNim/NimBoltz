## Cross sections of the gases Magboltz knows, by its gas number.

import gases/common
import gases/[argon, cf4, helium, neon, krypton, xenon, methane, ethane, propane, isobutane, carbondioxide, neopentane, water, oxygen, nitrogen, nitricoxide, nitrousoxide, ethylene, acetylene, hydrogen, deuterium, c3h6, hexafluoroethane, sulphurhexafluoride, methanol, ethanol, methylal, models, carbonmonoxide, dimethylether]
export common

proc crossSections*(gas: int, inp: GasInput, g: var GasXsec): bool =
  ## fills g (zeroed) with the cross sections of gas number `gas`; false
  ## for a gas that is still in legacy/
  case gas
  of 1: cf4(inp, g)
  of 2: argon(inp, g)
  of 3: helium4(inp, g)
  of 4: helium3(inp, g)
  of 5: neon(inp, g)
  of 6: krypton(inp, g)
  of 7: xenon(inp, g)
  of 8: methane(inp, g)
  of 9: ethane(inp, g)
  of 10: propane(inp, g)
  of 11: isobutane(inp, g)
  of 12: carbonDioxide(inp, g)
  of 13: neopentane(inp, g)
  of 14: water(inp, g)
  of 15: oxygen(inp, g)
  of 16: nitrogen(inp, g)
  of 17: nitricOxide(inp, g)
  of 18: nitrousOxide(inp, g)
  of 19: ethylene(inp, g)
  of 20: acetylene(inp, g)
  of 21: hydrogen(inp, g)
  of 22: deuterium(inp, g)
  of 23: carbonMonoxide(inp, g)
  of 24: methylal(inp, g)
  of 25: dimethylEther(inp, g)
  of 26: reidStep(inp, g)
  of 27: maxwell(inp, g)
  of 28: reidRamp(inp, g)
  of 29: hexafluoroethane(inp, g)
  of 30: sulphurHexafluoride(inp, g)
  of 32: propene(inp, g)
  of 33: cyclopropane(inp, g)
  of 34: methanol(inp, g)
  of 35: ethanol(inp, g)
  else: return false
  true
