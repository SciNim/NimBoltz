## The table here was simply converted to an Org table using emacs
## `org-table-create-or-convert-from-region`

const tab = """
| Routine  | Gas                   | Last update | Notes                                  | Star rating |
| GAS1     | CF4                   |        2001 | anisotropic scattering                 | 5*          |
| GAS2     | Argon                 |        1997 |                                        | 5*          |
| GAS3     | Helium 4              |        1997 |                                        | 5*          |
| GAS4     | Helium 3              |        1992 |                                        | 5*          |
| GAS5     | Neon                  |        1992 |                                        | 5*          |
| GAS6     | Krypton               |        2001 |                                        | 4*          |
| GAS7     | Xenon                 |        2001 |                                        | 4*          |
| GAS8     | methane               |        1994 |                                        | 5*          |
| GAS9     | ethane                |        1999 |                                        | 5*          |
| GAS10    | propane               |        1999 |                                        | 4*          |
| GAS11    | isobutane             |        1999 |                                        | 3*          |
| GAS12    | CO2                   |        2001 |                                        | 5*          |
| GAS13    | C(CH3)4 neo-pentane   |        1995 |                                        | 3*          |
| GAS14    | H20                   |        1998 |                                        | 3*          |
| GAS15    | Oxygen                |        1990 | 3-body attachment included             | 4*          |
| GAS16    | Nitrogen              |        1985 | Pitchford and Phelps                   | 4*          |
| GAS17    | nitric oxide          |        1995 | attaching gas                          | 4*          |
| GAS18    | nitrous oxide         |        1995 | attaching gas                          | 4*          |
| GAS19    | C2H4 ethene           |        1999 |                                        | 4*          |
| GAS20    | C2H2 acetylene        |        1992 |                                        | 3*          |
| GAS21    | Hydrogen              |        2001 |                                        | 5*          |
| GAS22    | Deuterium             |        1998 |                                        | 5*          |
| GAS23    | Carbon monoxide       |        1998 |                                        | 5*          |
| GAS24    | methylal              |        1988 |                                        | 2*          |
| GAS25    | DME                   |        1998 |                                        | 4*          |
| GAS26    | Reid step model       |           ? | anisotropic version                    | -           |
| GAS27    | Maxwell model         |           ? |                                        | -           |
| GAS28    | Reid ramp model       |           ? |                                        | -           |
| GAS29    | C2F6                  |        1999 | anisotropic                            | 4*          |
| GAS30    | SF6                   |           ? | do not use high percentage             | 3*          |
| GAS31    | NH3 ammonia           |        1999 |                                        | 3*          |
| GAS32    | C3H6 propene          |        1999 |                                        | 4*          |
| GAS33    | C3H6 cyclopropane     |        1999 |                                        | 4*          |
| GAS34    | CH3OH methanol        |        1999 |                                        | 2*          |
| GAS35    | C2H5OH ethanol        |        1999 |                                        | 3*          |
| GAS36    | C3H7OH isopropanol    |        1999 |                                        | 2*          |
| GAS37    | Cæsium                |        2001 | no dimers                              | 2*          |
| GAS38    | Fluorine              |           ? | Morgan                                 | 2*          |
| GAS39    | CS2                   |        2001 | ion drift, dark matter                 | 2*          |
| GAS40    | COS                   |        2001 |                                        | 2*          |
| GAS41    | CD4                   |        2001 | TPCs in neutron background environment | 3*          |
| GAS42    | BF3 Boron trifluoride |        2001 | anisotropic                            | 3*          |
| GAS43    | C2HF5 or C2H2F4       |           ? | estimated no data, anisotropic         | 1*          |
| GAS44    | Helium 3              |        2002 | anisotropic                            | 5*          |
| GAS45    | Helium 4              |        2002 | anisotropic                            | 5*          |
| GAS46    | Neon                  |        2002 | anisotropic                            | 5*          |
| GAS47    | Argon                 |        2002 | anisotropic                            | 5*          |
| GAS48    | Krypton               |        2002 | anisotropic                            | 4*          |
| GAS49    | Xenon                 |        2002 | anisotropic                            | 4*          |
| GAS50    | methane               |        2002 | anisotropic                            | 5*          |
| GAS52-80 | Dummy routines        |           ? |                                        | -           |
"""
import std / [strutils, tables, sequtils, sugar]
type
  GasMagboltz* = object
    id*: int
    name*: string
    lastUpdate*: string
    note*: string
    rating*: string

proc initMagboltzGases*(): seq[GasMagboltz] =
  result = newSeq[GasMagBoltz]()
  let lines = tab.strip.splitLines
  for i, line in lines:
    if i == 0: continue # skip header
    let data = line.strip(chars = {'|'}).split("|").mapIt(it.strip)
    let (idStr, gas, upd, note, rating) = (data[0], data[1], data[2], data[3], data[4])
    let idNoP = idStr.dup(removePrefix("GAS"))
    if "-" in idStr:
      let (strt, stop) = (parseInt idNoP.split("-")[0], parseInt idNoP.split("-")[1])
      for j in strt .. stop:
        result.add GasMagboltz(id: j, name: gas, lastUpdate: upd, note: note, rating: rating)
    else:
      result.add GasMagboltz(id: parseInt idNoP, name: gas, lastUpdate: upd, note: note, rating: rating)

var MagboltzGases = initMagboltzGases()

proc getGas*(id: int): GasMagboltz =
  for g in MagboltzGases:
    if g.id == id:
      return g
  raise newException(KeyError, "The gas with id = " & $id & " does not exist in Magboltz!")

proc getGas*(name: string): GasMagboltz =
  for g in MagboltzGases:
    if g.name.normalize == name.normalize:
      return g
  raise newException(KeyError, "The gas with name = " & $name & " does not exist in Magboltz!")
