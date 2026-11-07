## The description of a gas: what `engine.fill` needs to compute its cross
## sections on the energy grid.
##
## A gas is data: the tables of its measured cross sections, the energy
## losses of its levels and the variant of each process model it uses
## (docs/gas-features.md lists the variants by gas). Each gas module defines
## one `const GasModel`; `engine.fill` is instantiated per model, so the
## variants a gas does not use are not compiled into its routine.
##
## Units: energies in eV, cross sections in cm². Tables keep the values of
## Magboltz's DATA statements (single precision where Magboltz has them),
## widened to double; a table's `unit` (cm²) is the power of ten it is
## given in.
##
## Constants that Magboltz writes as single precision literals stay
## `float32` fields and are only widened or combined at run time: Nim's
## compile time evaluation of arithmetic on them does not round as the
## run time does (a widened `0.9'f32` times `1e-18` comes out as 0.9e-18).

import unchained
import common
export unchained, common

defUnit(cm², toExport = true)
defUnit(cm²•eV, toExport = true)

type
  Interp* = enum
    linear       ## linear in y over x
    logLog       ## linear in log y over log x
    linLog       ## linear up to `Elastic.switchAt`, log-log above

  Beyond* = enum
    ## what a cross section does above the end of its table
    extrapolate  ## the last table interval continued
    power        ## falls off as (x_end / E)^n from the last value
    born         ## Born-Bethe cross section of an allowed transition (`f`)
    bethe        ## Bethe ionisation cross section times `fraction`
    countingShare  ## `fraction` of the counting ionisation
                   ## (`Ionisation.counting`)
    zero         ## nothing

  Sigma* = object
    ## a cross section: a table, interpolated, and its continuation. Without
    ## table (`x` empty) only the continuation (`born`, `bethe`).
    x*, y*: seq[float64]   ## eV, and y in units of `unit`
    unit*: cm²
    scale*: float32 = 1.0  ## factor of the table values
    interp*: Interp
    beyond*: Beyond
    n*: float64            ## `power`: the exponent
    f*: float32            ## `born`: oscillator strength
    bornScale*: float32 = 1.0  ## `born`: factor of the Born-Bethe value
    fraction*: float32     ## `bethe`: share of the total
    clamp*: bool           ## table values below 0 are taken as 0

  LowEnergy* = enum
    lowNone
    phaseShifts    ## modified effective range theory below `lowBelow`
    constant       ## `lowQ` (total = momentum transfer) up to `lowBelow`

  Mert* = object
    ## parameters of the phase shifts of the modified effective range
    ## theory (polarisability in a₀³; scattering length, effective range
    ## terms). `tangent`: the phase shifts are arctan of the expressions.
    apol*, a*, d*, f*, a1*, b1*, a2*: float32
    tangent*: bool

  Elastic* = object
    interp*: Interp
    xTotal*, total*: seq[float64]        ## total elastic (1e-16 cm²)
    xMomentum*, momentum*: seq[float64]  ## momentum transfer (1e-16 cm²)
    xEps*, eps*: seq[float64]            ## 1 - forward/backward parameter
    epsInterp*: Interp
    switchAt*: float64   ## `linLog`: eV
    low*: LowEnergy
    lowBelow*: float64   ## eV
    lowQ*: float64       ## cm²
    lowIsotropic*: bool  ## `constant`: model 2 parameter 0 below lowBelow
    mert*: Mert

  IonTerm* = object
    ## part of an ionisation channel: a cross section times a factor
    ## (charge or number of equivalent atoms), from an energy on
    sigma*: Sigma
    onset*: float32 ## 0: from the channel's ionisation energy (eV)
    series*: seq[Sigma]
      ## charge states counted together (xenon 4, 5, 6): each table from its
      ## onset; above any table's end the share `sigma` (Bethe) replaces the
      ## whole sum
    seriesOnset*: seq[float32]
    seriesFactor*: seq[float64]     ## (charge / 3, double in Magboltz)

  IonChannel* = object
    eion*, eoby*: float32   ## ionisation energy, Opal-Beaty energy (eV)
    eobyScale*: seq[float32]  ## if given, Opal-Beaty energy = eion times
                              ## these factors, in order
    nc0*: int               ## Auger electrons
    ec0*: float32           ## their energy (eV)
    text*: string
    terms*: seq[IonTerm]
    shell*: bool            ## inner shell (not renormalised, subtracted)
    openInclusive*: bool    ## open from eion on (not above)
    gate*: float64 = 2.0    ## anisotropic above gate × eion

  IonScheme* = enum
    channels      ## each channel from its terms, then renormalised
    grossCounting ## helium: charge 1 = 2 counting - gross, charge 2 =
                  ## gross - counting

  Renormalisation* = enum
    noRenorm
    counting      ## outer channels × (counting - shells) / counting
    fragments     ## outer channels × (Σ outer - shells) / Σ outer
    subtractNext  ## each channel less the next one (isobutane)

  Ionisation* = object
    scheme*: IonScheme
    am2*, c*: float32        ## Bethe parameters of the total
    channels*: seq[IonChannel]
    renorm*: Renormalisation
    gateInclusive*: bool     ## anisotropic from gate × eion on (not above)
    counting*: Sigma         ## counting (or, for helium, gross) cross section
    countingFrom*: float32   ## zero up to this energy (eV)
    gross*: Sigma            ## helium: gross ionisation
    grossBethe*: float32     ## helium: Bethe of the gross = Bethe / this

  LevelKind* = enum
    tabulated     ## `sigma`, plus a vibrational dipole term (`dip`)
    superelastic  ## de-excitation of a vibration: the excitation's table
                  ## (`sigma`) at the energy after the gain, by detailed
                  ## balance, plus the dipole term; times the population of
                  ## the excited state over its degeneracy
    scaled        ## another level's table at energies scaled by the ratio
                  ## of the energy losses, times `factor` (helium n > 4)

  SuperelasticTail* = enum
    ## a superelastic level whose energy after the gain lies above the table
    squareTail    ## y_end (x_end / E')²
    inverseTail   ## y_end x_end / E
    cf4Tail       ## y_end x_end / (E E')² (CF4, see bugs/QUESTIONS.md Q-10)

  DipoleScale* = enum
    constantDipole
    cf4Dipole     ## times (14.4 - min(E, 5)) / 14 above 0.4 eV (CF4 V3)

  PeqKind* = enum
    ## the anisotropy parameter of a level
    peqElastic    ## isotropic (`peqBelow`) up to gate × ein, above it the
                  ## elastic scattering's at the energy after the loss
    peqHalf       ## 0.5; above `pq1Above` the elastic model 1 parameter
    peqBorn       ## from the momentum transfer of the dipole term and
                  ## `GasModel.momentumRatio` of the table part (if open);
                  ## above `pq1Above` the elastic model 1 parameter
    peqRatio      ## 0.5 + (1 - momentumRatio) (if open); above `pq1Above`
                  ## the elastic model 1 parameter
    peqScreened   ## as `peqElastic`, from gate × |ein| on, in anisotropy
                  ## model 2 only

  PeqBelow* = enum
    belowIso      ## 0.5, or 0 in anisotropy model 2
    belowZero

  Level* = object
    ein*: float32           ## energy loss (eV), < 0 superelastic
    text*: string
    kin*: int = -1          ## anisotropy kind; -1: NANISO
    kind*: LevelKind
    sigma*: Sigma
    dip*: float32           ## vibrational dipole constant (1e-16 cm² eV)
    dipScale*: DipoleScale
    tail*: SuperelasticTail ## `superelastic`
    population*: int        ## `superelastic`: index in GasModel.populations
    level*: int             ## `scaled`: the level whose table is used
    factor*: float32        ## `scaled`: factor (Born-Bethe with sigma.f if 0)
    peq*: PeqKind
    peqBelow*: PeqBelow
    gate*: float64 = 2.0    ## `peqElastic`: anisotropic above gate × ein
    gateEnergy*: float64    ## `peqElastic`: > 0: above this energy instead
    ground*: int            ## > 0: weighted by the ground state population
                            ## of this mode (`popPairs`)
    pq1Above*: float64 = Inf  ## eV

  Population* = object
    ## an excited vibrational state, by its superelastic level: degeneracy;
    ## the states of a group are normalised together with their ground
    ## state, deg exp(-E/kT) / (1 + Σ group); a ladder of n states (E, 2E,
    ## ... nE) counts as one, the sum of their populations
    level*: int
    deg*: float64 = 1.0
    group*: int
    ladder*: int = 1

  Attach* = object
    text*: string
    sigma*: Sigma           ## zero outside the table

  Penning* = object
    fraction*, length*, time*: float32
    format*: string
    all*: bool              ## print the levels without transfer as well
    first*, last*: int      ## levels with transfer (0: all)

  Names* = tuple[isotropic, anisotropic: string]

  GasModel* = object
    name*, elasticText*: Names
    mass*: float32           ## amu
    e3*: float32             ## ionisation energy (eV)
    e4*: float64             ## E(4) (helium 4)
    elastic*: Elastic
    ionisation*: Ionisation
    attachment*: seq[Attach]
    levels*: seq[Level]
    populations*: seq[Population]
      ## Boltzmann populations of excited vibrational states
    momentumRatio*: float32  ## momentum transfer over total of vibrational
                             ## tables (`peqBorn`, `peqRatio`)
    magnitudeOffsets*: bool  ## bins after a loss from |ein| (superelastic)
    negative*: Negative      ## Born-Bethe values below 0
    penning*: Penning
    firstCut*: int = 1       ## EFINAL cuts the levels from this one on

# ------------------------------------------------------------ builders

func widen*[T](a: openArray[T]): seq[float64] =
  for v in a: result.add float64(v)

func table*[X, Y](x: openArray[X], y: openArray[Y], unit = 1e-16.cm²,
                  interp = linear, clamp = false, scale = 1.0'f32): Sigma =
  ## a table, the last interval continued above its end
  Sigma(x: widen(x), y: widen(y), unit: unit, scale: scale, interp: interp,
        beyond: extrapolate, clamp: clamp)

func falling*[X, Y](x: openArray[X], y: openArray[Y], n: float64,
                    unit = 1e-18.cm², scale = 1.0'f32): Sigma =
  ## a table, falling off as E^-n above its end
  Sigma(x: widen(x), y: widen(y), unit: unit, scale: scale, beyond: power,
        n: n)

func withBorn*[X, Y](x: openArray[X], y: openArray[Y], f: float32,
                     unit = 1e-18.cm², scale = 1.0'f32): Sigma =
  ## a table, Born-Bethe with oscillator strength f above its end
  Sigma(x: widen(x), y: widen(y), unit: unit, scale: scale, beyond: born,
        f: f, bornScale: 1.0)

func bornOnly*(f: float32): Sigma =
  ## Born-Bethe alone
  Sigma(beyond: born, f: f, bornScale: 1.0)

func withBethe*[X, Y](x: openArray[X], y: openArray[Y], fraction: float32,
                      unit = 1e-16.cm², clamp = false): Sigma =
  ## a table, Bethe ionisation times fraction above its end
  Sigma(x: widen(x), y: widen(y), unit: unit, scale: 1.0, beyond: bethe,
        fraction: fraction, clamp: clamp)

func term*(s: Sigma, onset = 0.0'f32): IonTerm = IonTerm(sigma: s, onset: onset)

func channel*(eion, eoby: float32, nc0: int, ec0: float32, text: string,
              terms: varargs[IonTerm]): IonChannel =
  ## an outer channel, anisotropic above twice its energy
  IonChannel(eion: eion, eoby: eoby, nc0: nc0, ec0: ec0, text: text,
             terms: @terms, gate: 2.0)

func level*(ein: float32, text: string, s: Sigma, gate = 2.0): Level =
  Level(ein: ein, text: text, kind: tabulated, sigma: s, gate: gate)

func attach*(text: string): Attach = Attach(text: text, sigma: Sigma(beyond: zero))

func attach*[X, Y](text: string, x: openArray[X], y: openArray[Y],
                   unit = 1e-16.cm², scale = 1.0'f32): Attach =
  Attach(text: text, sigma: Sigma(x: widen(x), y: widen(y), unit: unit,
                                  scale: scale, beyond: zero))

func shell*(eion, eoby: float32, nc0: int, ec0: float32, text: string,
            s: Sigma, gate = 1.0): IonChannel =
  ## an inner shell: its table, not renormalised; anisotropic above gate ×
  ## eion (as soon as open)
  IonChannel(eion: eion, eoby: eoby, nc0: nc0, ec0: ec0, text: text,
             terms: @[term(s)], shell: true, gate: gate)

func pop*(level: int, deg = 1.0, group = 0, ladder = 1): Population =
  Population(level: level, deg: deg, group: group, ladder: ladder)
