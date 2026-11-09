## The printout of a run: PRNTER (conditions), OUTPUT (transport
## coefficients), OUTPUT1 and OUTPUT2 (collision frequencies, energy
## spectrum) and their thermal versions PRNTERT, OUTPUTT, OUTPUT2T. The
## thermal versions differ in a few formats.

import std/math
import frt, units, deck, conditions, tables, simulation

const
  fmtRule = "('-----------------------------------------------------------------------------')"
  # SETUP, SETUPT
  fmt87 = "(3(/),4X,' ERROR IN GAS INPUT : NGAS=',I5,6(/,2X,' N=',I3,' NGAS=',I5,' FRAC=',F8.3))"
  fmt91 = "(2(/)," &
          "' ERROR IN INPUT : NMAX TOO LARGE OVERFLOWED THE ALLOWED INTEGER RANGE OF THE COMPUTER OR COMPILER')"
  # PRNTER, PRNTERT
  fmt1 = "(2(/),10X,'PROGRAM MAGBOLTZ 2 VERSION 11.20',/)"
  fmt10 = "(10X,'MONTE CARLO SOLUTION FOR MIXTURE OF ',I2,' GASES.',/, 5X," &
          "'------------------------------------------------------')"
  fmt30 = "(/,5X,'  GASES  USED ',15X,' PERCENTAGE USED ',2(/),6(6X,A25,5X,F9.4,/))"
  fmt30t = "(/,5X,'  GASES  USED ',16X,'PERCENTAGE USED ',2(/),6(6X,A25,5X,F9.4,/))"
  fmt50 = "(/,2X,'GAS TEMPERATURE =',F6.1,' DEGREES CENTIGRADE.',/,2X,'GAS PRESSURE = ',F9.4," &
          "' TORR.')"
  fmt51 = "(1(/),2X,'INTEGRATION FROM 0.0 TO ',F8.2,' EV.  IN ',I4,' STEPS. ')"
  fmt52 = "(/,' PENNING EFFECTS NOT INCLUDED')"
  fmt53 = "(/,' PENNING EFFECTS INCLUDED')"
  fmt61 = "(/,' ISOTROPIC SCATTERING X-SECTIONS USED')"
  fmt62 = "(/,' ANISOTROPIC SCATTERING TYPE 1 (CAPITELLI/LONGO) USED IF AVAILABLE')"
  fmt63 = "(/,' ANISOTROPIC SCATTERING TYPE 2 (OKHRIMOVSKYY) USED IF AVAILABLE')"
  fmt65 = "(2(/),' THERMAL MOTION OF GAS NOT INCLUDED')"
  fmt66 = "(2(/),' THERMAL MOTION OF GAS INCLUDED')"
  fmt90 = "(1(/),'  ELECTRIC FIELD =',F12.4,' VOLTS/CM.',/'  MAGNETIC FIELD =',F11.4," &
          "' KILOGAUSS.',/,'  ANGLE BETWEEN ELECTRIC AND MAGNETIC FIELD =',F10.3,' DEGREES.'," &
          "/,'  CYCLOTRON FREQ. =',D12.3,' RADIANS/PICOSECOND')"
  fmt95 = "(1(/),'  INITIAL ELECTRON ENERGY =',F8.3,' EV.')"
  fmt100 = "(1(/),'  TOTAL NUMBER OF REAL COLLISIONS =',I11)"
  fmt110 = "(1(/)," &
           "'  NULL COLLISION FREQUENCY AT  8 EQUALLY SPACED ENERGY INTERVALS (*10**12/SEC)'," &
           "/,2(4(5X,D10.3)/))"
  fmt110t = "(1(/)," &
            "'  NULL COLLISION FREQUENCY FOR EACH GAS COMPONENT IN UNITS OF  (*10**12/SEC)'," &
            "/,2(3(8X,D10.3)/))"
  fmt111 = "('  REAL COLLISION FREQUENCY AT  8 EQUALLY SPACED ENERGY INTERVALS (*10**12/SEC)'," &
           "/,2(4(5X,D10.3)/))"
  fmt166 = "(/,' LONG DECORRELATION LENGTH =',I7,' COLLISIONS.')"
  fmt167 = "(/,' SHORT DECORRELATION LENGTH =',I7,' COLLISIONS.')"
  # OUTPUT, OUTPUTT
  fmt109 = "(/,2X,'CALCULATED MAX. COLLISION TIME =',F7.2,' PICOSECONDS.',2(/),2X," &
           "'NUMBER OF NULL COLLISIONS =',I11,/,2X,'NUMBER OF REAL COLLISIONS =',I11)"
  fmt333 = "(2(/),'  IONISATION RATE /CM.=',D11.4,' +/-',F6.2,' PERCENT.',/," &
           "'  ATTACHMENT RATE /CM.=',D11.4,' +/-',F6.2,' PERCENT.',2(/))"
  fmt940 = "(/,2X,'Z DRIFT VELOCITY =',E11.4,' MICRONS/NANOSECOND  +-',F8.2,'% ',/,2X," &
           "'Y DRIFT VELOCITY =',E11.4,' MICRONS/NANOSECOND  +-',F8.2,'%',/,2X," &
           "'X DRIFT VELOCITY =',E11.4,' MICRONS/NANOSECOND  +-',F8.2,'%',/)"
  fmt950 = "(/,2X,'TRANSVERSE DIFFUSION   =',D11.4,' +-',F8.2,'%',/,10X,'=',F9.4,' EV. +-'," &
           "F8.2,'%',/,10X,'=',F9.3,' MICRONS/CENTIMETER**0.5  +-',F8.2,'%',/)"
  fmt950t = "(/,2X,'TRANSVERSE DIFFUSION   =',D11.4,' +-',F8.2,'%',/,10X,'=',F9.5,' EV. +-'," &
            "F8.3,'%',/,10X,'=',F9.3,' MICRONS/CENTIMETER**0.5  +-',F8.2,'%',/)"
  fmt954 = "(/,10X,' DIFFUSION IN CM**2/SEC.',/)"
  fmt955 = "(/,2X,'DIFFUSION TENSOR :',/,6X,' DIFXX =',D11.4,' DIFYY =',D11.4,' DIFZZ ='," &
           "D11.4,/,6X,' DIFYZ =',D11.4,' DIFXY =',D11.4,' DIFXZ =',D11.4,/)"
  fmt956 = "(/,2X,'ERROR ON DIFFUSION TENSOR :',/,6X,' DIFXX =',F8.2,'%  DIFYY =',F8.2," &
           "'%  DIFZZ =',F8.2,'%',/,6X,' DIFYZ =',F8.2,'%  DIFXY =',F8.2,'%  DIFXZ =',F8.2," &
           "'%',/)"
  fmt957 = "(/,8X,' LONGITUDINAL DIFFUSION =',D11.4,' +-',F8.2,'%',/,10X," &
           "' TRANSVERSE DIFFUSION =',D11.4,' +-',F8.2,'%',/,2X," &
           "'TRANSVERSE DIFFUSION (PARALLEL TO B-FIELD) DIFXX=',D11.4,' +-',F8.2,'%',/)"
  fmt960 = "(/,2X,'MEAN ELECTRON ENERGY =',F9.4,' EV. ERROR =  +-',F8.2,'%',/)"
  fmt992 = "(/,2X,'LONGITUDINAL DIFFUSION =',D11.4,' +-',F8.1,'%',/,10X,'=',F9.4,' EV. +-'," &
           "F8.2,'%',/,10X,'=',F9.3,' MICRONS/CENTIMETER**0.5  +-',F8.2,'%',/)"
  # OUTPUT1, OUTPUT2, OUTPUT2T
  fmtLastBin = "(2(/),' NUMBER OF COLLISIONS IN FINAL ENERGY BIN =',F8.1)"
  fmt220o1 = "(/,6X,'TOTAL COLL. FREQ. =',D11.4,' (*10**12)/SEC.',/,2X," &
             "'INELASTIC COLL. FREQ. =',D11.4,' (*10**12)/SEC.',/,4X,'ELASTIC COLL. FREQ. ='," &
             "D11.4,' (*10**12)/SEC.',/,1X,'IONISATION COLL. FREQ. =',D11.4," &
             "' (*10**12)/SEC.',/,1X,'ATTACHMENT COLL. FREQ. =',D11.4,' (*10**12)/SEC.',/)"
  fmt220 = "(/,6X,'TOTAL COLL. FREQ. =',D11.4,' (*10**12)/SEC.',/,2X," &
           "'  ELASTIC COLL. FREQ. =',D11.4,' (*10**12)/SEC.',/,2X,'INELASTIC COLL. FREQ. ='," &
           "D11.4,' (*10**12)/SEC.',/,1X,'IONISATION COLL. FREQ. =',D11.4,' (*10**12)/SEC.'," &
           "/,1X,'ATTACHMENT COLL. FREQ. =',D11.4,' (*10**12/SEC.',/)"
  fmt220t = "(/,6X,'TOTAL COLL. FREQ. =',D11.4,' (*10**12)/SEC.',/,2X," &
            "'  ELASTIC COLL. FREQ. =',D11.4,' (*10**12)/SEC.',/,2X," &
            "'INELASTIC COLL. FREQ. =',D11.4,' (*10**12)/SEC.',/,1X," &
            "'IONISATION COLL. FREQ. =',D11.4,' (*10**12)/SEC.',/,1X," &
            "'ATTACHMENT COLL. FREQ. =',D11.4,' (*10**12)/SEC.',/)"
  fmt301 = "(2(/),10X,' NORMALISED ENERGY DISTRIBUTION')"
  fmt302 = "(6X,'E=',F11.3,6X,'SPEC=',D11.4)"
  fmt1060 = "(/,2X,'DETAILED COLLISION FREQUENCIES FOR EACH GAS IN UNITS OF 10**12/SEC. :',2(/))"
  fmt1065 = "(3X,A25,/,'-----------------------------',2(/))"
  fmt1065t = "(/,3X,A25,/,'----------------------------',2(/))"
  fmt1070 = "(3X,A50,3X,D11.4,' +-',F8.4,' %')"

proc printSetupError*(o: File, d: Deck, e: SetupError) =
  case e
  of seGasInput:
    var w: seq[FVal]
    w.add toFVal(d.ngas)
    for j in 1 .. maxGases:
      w.add toFVal(j)
      w.add toFVal(d.gas[j])
      w.add toFVal(d.percent[j])
    o.fwrite(fmt87, w)
  of seNmaxOverflow: o.fwrite(fmt91)
  of seNone: discard

proc printHeader*[M](r: Run[M]) =
  ## PRNTER, PRNTERT
  const thermal = thermalMotion(M)
  let d = r.deck
  r.output.fwrite(fmt1)
  r.output.fwrite(fmt10, d.ngas)
  var w: seq[FVal]
  for j in 1 .. d.ngas:
    w.add toFVal(r.names[j])
    w.add toFVal(d.percent[j])
  r.output.fwrite(when thermal: fmt30t else: fmt30, w)
  r.output.fwrite(fmt50, d.tempc, d.torr)
  r.output.fwrite(fmt51, r.efinal, nEnergyBins)
  r.output.fwrite(if d.penning: fmt53 else: fmt52)
  r.output.fwrite(case naniso
         of 0: fmt61
         of 1: fmt62
         else: fmt63)
  if r.cond.dec.ncolm > 1000000: r.output.fwrite(fmt166, r.cond.dec.ncolm)
  else: r.output.fwrite(fmt167, r.cond.dec.ncolm)
  r.output.fwrite(when thermal: fmt66 else: fmt65)
  r.output.fwrite(fmt90, d.efield, d.bmag, d.btheta, r.cond.wb)
  r.output.fwrite(fmt95, r.estart)
  r.output.fwrite(fmt100, r.cond.nmax)
  w.setLen 0
  when thermal:
    for l in 1 .. maxGases: w.add toFVal(r.mix.gas.tcfmax[l])
    r.output.fwrite(fmt110t, w)
  else:
    for l in 1 .. 8: w.add toFVal(r.mix.bounds[l])
    r.output.fwrite(fmt110, w)
    w.setLen 0
    for l in countup(250, 3750, 500): w.add toFVal(r.mix.table.tcf[l])
    r.output.fwrite(fmt111, w)

proc printTransport*[M](r: Run[M]) =
  ## OUTPUT, OUTPUTT: results of MONTE
  let d = r.deck
  let t = r.transport
  r.output.fwrite(fmtRule)
  r.output.fwrite(fmtRule)
  r.output.fwrite(fmt109, r.tmax, r.nnull, r.collisions)
  r.output.fwrite(fmt940, t.vel.z.to(μm•ns⁻¹), t.velErr.z, t.vel.y.to(μm•ns⁻¹),
         t.velErr.y, t.vel.x.to(μm•ns⁻¹), t.velErr.x)
  if d.bmag.float > 0.0 and (d.btheta > 0.0 and d.btheta < 180.0):
    r.output.fwrite(fmt954)
    r.output.fwrite(fmt955, t.dif.xx, t.dif.yy, t.dif.zz, t.dif.yz,
           t.dif.xy, t.dif.xz)
    r.output.fwrite(fmt956, t.difErr.xx, t.difErr.yy, t.difErr.zz, t.difErr.yz, t.difErr.xy,
           t.difErr.xz)
    if d.btheta == 90.0:
      r.output.fwrite(fmt957, t.difln, t.dfler, t.diftr, t.dfter, t.dif.xx,
             t.difErr.xx)
  else:
    # characteristic energies D/µ, and the spread per √(drift distance)
    let dtovmb = t.diftr * d.efield / t.vel.z * perElectron
    let dtmn = sqrt(2.0 * t.diftr / t.vel.z).to(μm•cm⁻¹ᐟ²)
    let dfter1 = sqrt(t.dfter * t.dfter + t.velErr.z * t.velErr.z)
    r.output.fwrite(fmt954)
    r.output.fwrite(when thermalMotion(M): fmt950t else: fmt950, t.diftr, t.dfter,
           dtovmb, dfter1, dtmn, dfter1 / 2.0)
    let dlovmb = t.difln * d.efield / t.vel.z * perElectron
    let dlmn = sqrt(2.0 * t.difln / t.vel.z).to(μm•cm⁻¹ᐟ²)
    let dfler1 = sqrt(t.dfler * t.dfler + t.velErr.z * t.velErr.z)
    r.output.fwrite(fmt992, t.difln, t.dfler, dlovmb, dfler1, dlmn, dfler1 / 2.0)
  r.output.fwrite(fmt333, t.alpha, t.alper, t.att, t.atter)
  r.output.fwrite(fmt960, t.ave, t.den)

proc printProcesses[M](r: Run[M], f: Frequencies) =
  r.output.fwrite(fmt1060)
  for j in 1 .. r.deck.ngas:
    r.output.fwrite(when thermalMotion(M): fmt1065t else: fmt1065, r.names[j])
    for (descr, n) in r.processes(j):
      let frelv = f.total * float64(n) / float64(f.n)
      let errfre = if n == 0: 0.0 else: 100.0 * sqrt(float64(n)) / float64(n)
      r.output.fwrite(fmt1070, descr, frelv, errfre)
  r.output.fwrite(fmtRule)

proc printCollisionCounts*(r: Run[ColdMixture]) =
  ## OUTPUT1: collision frequencies of the steady-state Townsend run
  r.output.fwrite(fmtRule)
  r.output.fwrite(fmtRule)
  r.output.fwrite(fmtLastBin, r.last.spectrum[nEnergyBins])
  let f = frequencies(r.last, r.deck.ngas)
  r.output.fwrite(fmt220o1, f.total, f.inelastic, f.elastic, f.ionisation,
         f.attachment)
  r.output.fwrite(fmtRule)
  r.printProcesses(f)

proc printFrequencies*[M](r: Run[M]) =
  ## OUTPUT2, OUTPUT2T: collision frequencies and energy spectrum of the
  ## last run
  r.output.fwrite(fmtRule)
  r.output.fwrite(fmtRule)
  r.output.fwrite(fmtLastBin, r.last.spectrum[nEnergyBins])
  let f = frequencies(r.last, r.deck.ngas)
  r.output.fwrite(when thermalMotion(M): fmt220t else: fmt220, f.total, f.elastic,
         f.inelastic, f.ionisation, f.attachment)
  r.output.fwrite(fmtRule)
  r.printProcesses(f)
  # spectrum normalised to the collisions, in 40 bins of 100
  r.output.fwrite(fmt301)
  var specs: array[1..40, float64]
  let specn = float64(f.n)
  var smspec = 0.0
  for k in 1 .. nEnergyBins:
    smspec = smspec + r.last.spectrum[k] / specn
    if k mod 100 == 0:
      specs[k div 100] = smspec
      smspec = 0.0
  let eplt = r.efinal / 40.0
  for i in 1 .. 40:
    r.output.fwrite(fmt302, (eplt * (float64(i) - 0.5)), specs[i])
