## Cross sections of the gases Magboltz knows, by its gas number.

import gases/common
import gases/[argon, cf4, helium, neon, krypton, xenon, methane, ethane, propane, isobutane]
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
  else: return false
  true
