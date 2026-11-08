## Units of NimBoltz (Unchained). Every quantity is kept in the unit
## Magboltz computes it in, so that the arithmetic stays that of the Fortran:
##
## - energies in eV, times in ps, positions in m, speeds in m/ps while an
##   electron is followed;
## - fields in V/cm and kG, cross sections in cm², collision frequencies in
##   1/ps (Magboltz's "10¹²/s");
## - results in cm/s, cm²/s, 1/cm (and µm/ns where Magboltz prints those).
##
## Magboltz writes its conversions as literal factors (`* 1D-9`, `* 0.01`).
## Here each is a `unity` of the unit it converts between, i.e. the number 1
## written in that unit, so the code says what it converts and Unchained
## checks it; `conversions` below asserts that each equals the literal
## Magboltz uses. Speeds come from energies as v = vfac·√E, hence the
## half-integer powers of eV.

import unchained
export unchained

# energies and speeds
defUnit(eV¹ᐟ², toExport = true)            ## √E, as in v = vfac·√E
defUnit(eV⁻¹ᐟ², toExport = true)
defUnit(eV⁻¹, toExport = true)
defUnit(eV•V⁻¹, toExport = true)           ## `perElectron`
defUnit(eV•ps, toExport = true)            ## ∫ e dt
defUnit(eV•m⁻¹, toExport = true)           ## energy gain per displacement
defUnit(cm•ps⁻¹•eV⁻¹ᐟ², toExport = true)    ## CONST3
defUnit(m•ps⁻¹, toExport = true)           ## speed while following an electron
defUnit(m²•ps⁻², toExport = true)
defUnit(m•ps⁻¹•eV⁻¹ᐟ², toExport = true)     ## vfac
defUnit(m²•ps⁻²•eV⁻¹, toExport = true)      ## vfac²
defUnit(cm•ps⁻¹•V⁻¹ᐟ², toExport = true)     ## CONST3: speed per √(energy per charge)
defUnit(μm•ns⁻¹, toExport = true)          ## drift velocities as printed
defUnit(cm•s⁻¹, toExport = true)           ## drift velocities of the results
defUnit(m•s⁻¹, toExport = true)
# fields
defUnit(V•cm⁻¹, toExport = true)           ## E as given
defUnit(V•m⁻¹, toExport = true)
defUnit(V•cm⁻¹•Torr⁻¹, toExport = true)    ## E/p
defUnit(ps⁻¹, toExport = true)             ## cyclotron and collision frequencies
defUnit(m•ps⁻², toExport = true)           ## acceleration term of the position
defUnit(eV•ps⁻², toExport = true)          ## acceleration term of the energy
defUnit(eV¹ᐟ²•ps⁻¹, toExport = true)
defUnit(cm²•V⁻¹•ps⁻², toExport = true)      ## CONST1: (e/m)/2
defUnit(m•cm•V⁻¹•ps⁻², toExport = true)     ## CONST2
defUnit(eV¹ᐟ²•cm•V⁻¹•ps⁻¹, toExport = true) ## CONST5
# the gas
defUnit(cm⁻³, toExport = true)             ## number density
defUnit(cm², toExport = true)              ## cross sections
defUnit(cm²•eV, toExport = true)           ## Born-Bethe constant
defUnit(cm²•eV¹ᐟ², toExport = true)        ## cross section ∝ 1/√E (GAS27)
defUnit(cm²•eV⁻¹, toExport = true)         ## cross section ∝ E (GAS28)
defUnit(eV•K⁻¹, toExport = true)           ## Boltzmann's constant
defUnit(cm⁻²•ps⁻¹•eV⁻¹ᐟ², toExport = true)  ## gas density times the speed of 1 eV (VAN): the
                                            ## collision frequency per cm² of cross section
                                            ## and √eV of energy
defUnit(ps⁻¹•eV⁻¹ᐟ², toExport = true)       ## collision frequency / √E (MIXER)
# results
defUnit(cm²•s⁻¹, toExport = true)          ## diffusion
defUnit(μm•cm⁻¹ᐟ², toExport = true)        ## spread per √(drift distance)
defUnit(m², toExport = true)
defUnit(m²•ps⁻¹, toExport = true)
defUnit(cm⁻¹, toExport = true)             ## Townsend coefficients
defUnit(m⁻¹, toExport = true)

# sums over the electrons crossing a space plane, weighted by 1/|vz|
defUnit(ps•m⁻¹, toExport = true)
defUnit(ps²•m⁻¹, toExport = true)
defUnit(ps³•m⁻¹, toExport = true)
defUnit(ps²•m⁻², toExport = true)
defUnit(m•ps, toExport = true)
defUnit(eV•ps•m⁻¹, toExport = true)
defUnit(eV², toExport = true)

# dimensionless units of the conversions
defUnit(m•cm⁻¹, toExport = true)
defUnit(cm•m⁻¹, toExport = true)
defUnit(T•kGauss⁻¹, toExport = true)
defUnit(m•ps⁻¹•kGauss•cm•V⁻¹, toExport = true)
defUnit(C•kg⁻¹•T•kGauss⁻¹, toExport = true)
defUnit(cm⁻²•ps⁻¹•V⁻¹ᐟ², toExport = true)

template unity*(U: untyped): untyped =
  ## The number 1 in the dimensionless unit U, e.g. `unity(m•cm⁻¹)` is
  ## 0.01 m/cm: multiplying by it converts cm to m.
  1.0.UnitLess.to(U)

const perElectron* = 1.0.eV•V⁻¹
  ## The (magnitude of the) electron's charge in eV per volt: the energy an
  ## electron gains through a potential difference. Multiplying by it
  ## changes no bit.

template inUnits*(x: SomeUnit): float64 =
  ## the value of x in its unit (for printing, and for the logarithms of
  ## Magboltz that only enter through their differences)
  x.float
