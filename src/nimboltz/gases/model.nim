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
    linLog       ## linear up to `switchAt` (Elastic, Sigma), log-log above

  Beyond* = enum
    ## what a cross section does above the end of its table
    extrapolate  ## the last table interval continued
    power        ## falls off as (x_end / E)^n from the last value
    born         ## Born-Bethe cross section of an allowed transition (`f`)
    bethe        ## Bethe ionisation cross section times `fraction`
    countingShare  ## `fraction` (times `tailScale`) of the counting ionisation
                   ## (`Ionisation.counting`)
    zero         ## nothing

  BornForm* = enum
    bornStandard  ## times E / (E + ein + E(3))
    bornShifted   ## at E + 2 ein, over E (neutral dissociation, N2)

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
    fd*: float64           ## `born`: if not 0, the strength in double precision
    tailScale*: float32 = 1.0  ## `born`, `bethe`: factor of the continuation
    fraction*: float32     ## `bethe`: share of the total
    betheFraction*: float32
      ## `countingShare`: if not 0, the share while the counting ionisation
      ## is Born-Bethe (above its table)
    clamp*: bool           ## table values below 0 are taken as 0
    upTo*: float64         ## > 0: the table only up to this energy (eV)
    bornForm*: BornForm    ## `born`: its energy factor
    switchAt*: float64     ## `linLog`: eV
    zeroBelow*: bool       ## nothing at or below the table's first energy
    zeroUnder*: bool       ## nothing below the table's first energy
    zeroAbove*: bool       ## nothing above the table (the level stays open)
    logZeros*: bool        ## `logLog`: nothing in intervals with a zero end

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
    lowQMomentum*: float64  ## cm² (0: lowQ)
    lowIsotropic*: bool  ## `constant`: model 2 parameter 0 below lowBelow
    alwaysAnisotropic*: bool
      ## model 2 (1 - eps, total cross section, KEL = 2) whatever NANISO
    epsConstantBelow*: bool  ## eps = its first value below the table
    scale*: float64 = 1.0    ## factor of the cross sections (D2 0.988)
    mert*: Mert
    bend*: Bend
    isotropicBelow*: float64
      ## eV; below it momentum transfer and isotropic (if anisotropic)
    minusLevels*: Slice[int]
      ## levels (rotations) whose sum is taken out of the elastic cross
      ## section
    minusFloor*: float64
      ## cm²; > 0: the elastic cross section if not above 0 after that

  Bend* = object
    ## elastic scattering of the molecules in excited bending states
    ## (CO2): their tables, the share of them in the tables of the gas, and
    ## the populations whose sum gives their share in the gas
    total*, momentum*: seq[float64]   ## on Elastic.xTotal (1e-16 cm²)
    share*: float64
    populations*: seq[int]

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
    shiftFrom*: int
      ## > 0: the table at the energy less the difference of this channel's
      ## ionisation energy and that of channel `shiftFrom` (O2+ A, B from
      ## the X table)
    rampTo*: float64
      ## > 0 (eV): linear from 0 at the ionisation energy to `rampQ` (1e-16
      ## cm²) here, from here on `share` of channel `shareOf` (N2)
    rampQ*, share*: float32
    shareOf*: int

  IonChannel* = object
    eion*, eoby*: float32   ## ionisation energy, Opal-Beaty energy (eV)
    eobyScale*: seq[float32]  ## if given, Opal-Beaty energy = eion times
                              ## these factors, in order
    eobyBase*: int            ## > 0: the eion of this channel instead
    nc0*: int               ## Auger electrons
    ec0*: float32           ## their energy (eV)
    text*: string
    terms*: seq[IonTerm]
    shell*: bool            ## inner shell (not renormalised, subtracted)
    openInclusive*: bool    ## open from eion on (not above)
    gateInclusive*: bool    ## anisotropic from gate × eion on (not above)
    gate*: float64 = 2.0    ## anisotropic above gate × eion

  IonScheme* = enum
    channels      ## each channel from its terms, then renormalised
    grossCounting ## helium: charge 1 = 2 counting - gross, charge 2 =
                  ## gross - counting
    single        ## one cross section (`counting`) in Q(3), from
                  ## `countingFrom` on

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
    split*: Split            ## channels 1 .. sharing one cross section
    subtract*: seq[Subtraction]  ## after the renormalisation, in order
    fromTables*: bool        ## channels open above their table's first
                             ## energy (not their ionisation energy)
    stopAtClosed*: bool      ## `fromTables`: and the first one closed
                             ## closes the following ones
    peqZero*: bool           ## anisotropy parameter 0 (not isotropic)
                             ## below the gates
    counting*: Sigma         ## counting (or, for helium, gross) cross section
    countingFrom*: float32   ## zero up to this energy (eV)
    countingExclusive*: bool ## and at it
    gross*: Sigma            ## helium: gross ionisation
    grossBethe*: float32     ## helium: Bethe of the gross = Bethe / this

  Split* = object
    ## one cross section (from channel 1's energy on) shared by the open
    ## ones of channels 1 .. fractions.len + 1: with channels 1 .. n open,
    ## `fractions[n - 2]` (channel 1 alone: all of it) (N2+ X, A, B, C)
    sigma*: Sigma
    fractions*: seq[seq[float32]]

  Subtraction* = object
    ## a channel taken out of another (excited ion states, Auger electrons
    ## of inner shells): target -= factor × channel, up to `below` (eV)
    target*, channel*: int
    factor*: float64 = 1.0
    below*: float64 = Inf
    firstTerm*: bool       ## only the first term of the channel

  LevelKind* = enum
    tabulated     ## `sigma`, plus a vibrational dipole term (`dip`)
    superelastic  ## de-excitation of a vibration: the excitation's table
                  ## (`sigma`) at the energy after the gain, by detailed
                  ## balance, plus the dipole term; times the population of
                  ## the excited state over its degeneracy
    scaled        ## another level's table at energies scaled by the ratio
                  ## of the energy losses, times `factor` (helium n > 4)
    rotation      ## J -> J ± 2 of a linear rotor, Born quadrupole
                  ## (`GasModel.rotor`)
    dipoleRotation  ## a rotational transition, Born dipole: from state
                    ## `state` (or J = `j`) with line strength `strength`
                    ## (or J + 1, J); `GasModel.rotor`
    symmetricTop
      ## J -> J ± 1 of a symmetric top, Born dipole summed over K: Σ pop
      ## (J'² - K²) over the states `top` (or pop J'² of the one state if not
      ## `bracket`), over J' (2 aj + 1) E (`ajg` = J'); `Rotor.states`,
      ## `Rotor.dipole` (NH3)
    tabulatedRotation
      ## J <-> J + 2 (`j` the lower J) from a table of the excitation
      ## (`sigma`), plus `Rotor.resonance`; the de-excitation by detailed
      ## balance with the ratio of the degeneracies; populations of
      ## `GasModel.rotor`; above `Rotor.scaleAbove` times scaleAbove / E, all
      ## times `Rotor.tabulatedScale` (O2)

  SuperelasticTail* = enum
    ## a superelastic level whose energy after the gain lies above the table
    squareTail    ## y_end (x_end / E')²
    inverseTail   ## y_end x_end / E
    crossTail     ## y_end x_end E' / E²
    mixedTail     ## y_end x_end / (E' E)
    e2Tail        ## y_end x_end / E'
    powerTail     ## y_end (x_end / E')^n (`Sigma.n`)
    zeroTail      ## nothing
    noTail        ## the table's last interval continued
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
    peqMomentumRatio
                  ## isotropic up to gate × ein; above, model 1: 1.5 - the
                  ## ratio of momentum transfer to total (`Level.ratio` on
                  ## the table's energies, falling as 1/E above), model 2:
                  ## the elastic scattering's at the energy after the loss
    peqSameAs     ## that of level `Level.level` (if open)
    peqFixed      ## `Level.peqValue`
    peqDipole     ## Born dipole scattering: 1 - eps from
                  ## `GasModel.dipoleAngles` at E / |ein| (up to
                  ## `dipoleUpTo`; above, the elastic scattering's at the
                  ## energy after the loss); from `gateEnergy` on

  PeqBelow* = enum
    belowIso      ## 0.5, or 0 in anisotropy model 2
    belowZero

  Level* = object
    ein*: float32           ## energy loss (eV), < 0 superelastic
    text*: string
    kin*: int = -1          ## anisotropy kind; -1: NANISO (0 in an
                            ## isotropic routine)
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
    gateAbs*: bool          ## gate × |ein|
    gateInclusive*: bool    ## from gate × ein on (not above)
    gateEnergy*: float64    ## `peqElastic`: > 0: above this energy instead
    dampAbove*: float64     ## > 0: times sqrt(dampAbove / E) above it
    extra*: Sigma           ## `tabulated`: added (e.g. a triplet table to
                            ## a Born-Bethe level); `beyond` set if used
    clampTotal*: bool       ## the sum taken as 0 if below
    openAbove*: float64     ## > 0: open above this energy (eV) instead of ein
    dimer*: bool            ## times the density of `GasModel.dimer`
    windows*: seq[Sigma]
      ## tables added to the dipole term inside their range (both ends
      ## included), each times its scale, at the energy after the gain for
      ## a superelastic level (E' y(E') / E): resonances (N2O, C2H4)
    windowOpen*: bool       ## the ends of the windows excluded
    rotState*: int = -1     ## >= 0: times the population of rotational
                            ## state J = rotState (`GasModel.rotor`) and
                            ## `weight`
    weight*: float64 = 1.0
    closedAbove*: float64 = Inf  ## eV
    ratio*: seq[float64]    ## `peqMomentumRatio`
    peqValue*: float64      ## `peqFixed`
    dipoleUpTo*: float64 = Inf  ## `peqDipole`
    eloss*: float64         ## != 0: energy loss in double precision
    j*: int                 ## `rotation`, `dipoleRotation`: initial J
    up*: bool               ## `rotation`, `dipoleRotation`: J up (else down)
    state*: int             ## `dipoleRotation`: initial state (1-based) of
                            ## `Rotor.states` (0: by J)
    strength*: float64      ## `dipoleRotation`: line strength (0: J + 1, J)
    top*: seq[tuple[state: int, k2: float64]]  ## `symmetricTop` (1-based)
    ajg*: float64           ## `symmetricTop`
    bracket*: bool          ## `symmetricTop`
    aj*: float64            ## `dipoleRotation` from a state: the J of the
                            ## 2J + 1 of the cross section (Magboltz's AJIN;
                            ## not always the state's)
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

  AttachFormula* = enum
    attachTable
    h2Attachment
      ## H2: resonances from the rotational states, then the table and a
      ## high resonance (`Attach.h2`)
    cf3brAttachment
      ## (0.0353 / √E exp(-E²/9e-6) + exp(-E / 0.0588)) 27.4e-16 cm² up to
      ## 2 eV (CF3Br)

  H2Attachment* = object
    ## Σ_J amp strength_J pop_J exp(-(E - threshold + E_J) / width) from
    ## J = n down to 1, open from threshold - E_J on; amp pop_0 exp(-(E -
    ## threshold) / width) above threshold; the table; highAmp exp(-(E -
    ## highAt) / highWidth) above highAt (single precision constants of
    ## Magboltz as their doubles)
    threshold*, width*, amp*: float64
    strengths*: seq[float64]
    highAt*, highWidth*, highAmp*: float64

  Dimer* = object
    ## a level of the dimers of a gas (Hg2 in mercury): tables for pure
    ## gas `gas` at `temperature` (K) and 760 Torr, times the dimer
    ## fraction `fraction`
    gas*: int
    temperature*, fraction*: float32

  Attach* = object
    text*: string
    sigma*: Sigma           ## zero below the table (and above, if `zero`)
    fromStart*: bool        ## open at the table's first energy already
    beforeEnd*: bool        ## `zero`: closed at the table's last energy
    threeBody*: bool        ## table for 760 Torr, 20 °C: times the density
    openBelow*: bool        ## not closed below the table (extrapolated)
    formula*: AttachFormula
    saturation*: float64    ## > 0: times 1 / (saturation / TORR + 1) (three-
                            ## body attachment saturating with the pressure)
    factors*: seq[float32]  ## further factors, in order (single precision)
    isoPeq*: bool           ## sets the attachment's (unused) anisotropy
                            ## parameter to the isotropic value
    h2*: H2Attachment

  Rotor* = object
    ## rotational states: of a linear molecule, energies B0 J (J + 1) for
    ## J = jStart, jStart + jStep .. jMax; or tabulated states (energy,
    ## J, statistical weight). Boltzmann populations w (2J + 1) exp(-E/kT),
    ## normalised. Quadrupole (J ± 2) and dipole (J ± 1 or tabulated)
    ## transitions in the Born approximation.
    b0*: float64            ## eV
    jStart*, jMax*: int
    jStep*: int = 1         ## 2: even (CO2) or odd (O2) J only
    states*: seq[RotState]
    resonance*: Sigma       ## `tabulatedRotation`: inside its table only
    tabulatedScale*: float64 = 1.0
    quadrupole*: float64    ## cm²
    scaleAbove*: float64    ## quadrupole: above it times scaleAbove / E
    dipole*: float64        ## cm² eV
    dipoleScale*: float32 = 1.0   ## factor of the dipole transitions
    dipoleCutoff*: float64 = Inf  ## dipole transitions up to this energy
    ladder*: seq[int]
      ## states J = 0 .. jMax with energies built from the losses of the
      ## levels J -> J + 2 (ladder[J]): E(J + 2) = E(J) + ein; E(1) = oddBase
      ## (H2, D2)
    oddBase*: float64
    fixBelow*: float64
      ## kT (eV); below it the ortho/para ratio stays at its room temperature
      ## value: the even states' total = odd total × fixRatio.mul /
      ## fixRatio.den
    fixRatio*: tuple[mul, den: float64]
    weights*: array[2, float64] = [1.0, 1.0]
      ## statistical weights of even and odd J (nuclear spin; N2 6 and 3)
    enhance*: Sigma
      ## quadrupole rotations times 1 + y enhanceNum / enhanceDen, y this
      ## table at E (excitation) or E' y(E') / E (de-excitation): the
      ## resonance of the vibrations (N2)
    enhanceNum*, enhanceDen*: float32
    momentumAbove*: float64
      ## > 0 (eV): from here on the rotations times the momentum transfer
      ## cross section over `momentumRef` (cm²)
    momentumRef*: float64

  RotState* = tuple
    e: float64              ## eV
    j: float64
    weight: float64

  DipoleAngles* = object
    ## angular distribution of Born dipole scattering, 1 - eps over E / |ein|
    ## (excitation and de-excitation)
    xUp*, xDown*, eps*: seq[float64]

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
    rotor*: Rotor
    dipoleAngles*: DipoleAngles
    momentumRatio*: float32  ## momentum transfer over total of vibrational
                             ## tables (`peqBorn`, `peqRatio`)
    magnitudeOffsets*: bool  ## bins after a loss from |ein| (superelastic)
    negative*: Negative      ## Born-Bethe values below 0
    bbBohr*, bbRydberg*: float64
      ## > 0: the routine's own Born-Bethe constant 16 π a₀² Ry² / (2 m c²)
      ## from these a₀ (cm) and Ry (eV)
    penning*: Penning
    firstCut*: int = 1       ## EFINAL cuts the levels from this one on
    minLevels*: int          ## but keeps at least these
    isotropic*: bool
      ## an old routine: isotropic whatever NANISO (anisotropy kinds 0 and
      ## parameters left 0), the momentum transfer cross section as elastic
    stepped*: bool
      ## an old routine: energies summed up bin by bin (EN = EN + ESTEP)
      ## instead of the grid's
    alwaysModel2*: bool
      ## NANISO taken as 2 whatever the input (C2H4)
    forceModel2*: bool
      ## anisotropic means model 2 (NANISO 1 taken as 2), all KEL and KIN
      ## 2; isotropic: the momentum transfer cross section with parameter 0
    fixedModel1*: bool
      ## elastic and ionisation anisotropic, model 1, whatever NANISO (KEL
      ## 1 for Q(2) and Q(3)), the rest isotropic; the ionisation's angles
      ## those of the elastic scattering above 2 E(3)
    roundedOffsets*: bool
      ## bins after a loss: the loss / ESTEP rounded (single precision), as
      ## some old routines
    isotropicPeq*: float64   ## the (unused) anisotropy parameter an
                             ## isotropic routine stores for Q(2), Q(3)
    cutStrict*: bool         ## EFINAL cuts the levels strictly above it
    keepLastAt*: bool        ## no cut if EFINAL is the last level's energy
    dimer*: Dimer
    ninBy*: seq[tuple[efinal: float64, nin: int]]
      ## instead of the cut: the number of levels used, the first one whose
      ## EFINAL is below the input's

# ------------------------------------------------------------ builders

func widen*[T](a: openArray[T]): seq[float64] =
  for v in a: result.add float64(v)

func f64*[T](a: openArray[T], k: int): float64 =
  ## element k of a table, widened. (Not `float64(a[k])`: for a `const`
  ## float32 table Nim folds that from the decimal literal, unrounded, even
  ## in code run at run time.)
  for i, v in a:
    if i == k: return float64(v)

func table*[X, Y](x: openArray[X], y: openArray[Y], unit = 1e-16.cm²,
                  interp = linear, clamp = false, scale = 1.0'f32,
                  zeroBelow = false): Sigma =
  ## a table, the last interval continued above its end
  Sigma(x: widen(x), y: widen(y), unit: unit, scale: scale, interp: interp,
        beyond: extrapolate, clamp: clamp, zeroBelow: zeroBelow)

func falling*[X, Y](x: openArray[X], y: openArray[Y], n: float64,
                    unit = 1e-18.cm², scale = 1.0'f32): Sigma =
  ## a table, falling off as E^-n above its end
  Sigma(x: widen(x), y: widen(y), unit: unit, scale: scale, beyond: power,
        n: n)

func withBorn*[X, Y](x: openArray[X], y: openArray[Y], f: float32,
                     unit = 1e-18.cm², scale = 1.0'f32): Sigma =
  ## a table, Born-Bethe with oscillator strength f above its end
  Sigma(x: widen(x), y: widen(y), unit: unit, scale: scale, beyond: born,
        f: f, tailScale: 1.0)

func bornOnly*(f: float32): Sigma =
  ## Born-Bethe alone
  Sigma(beyond: born, f: f, tailScale: 1.0)

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
