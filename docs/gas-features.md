# Gas cross sections: feature matrix

Every gas routine (GAS1..GAS80 of Magboltz 11.21) computes the same things on the
energy grid: elastic (total, momentum transfer, angular parameter), ionisation
channels, attachment, inelastic levels, null collisions. The physics is shared
(`gases/common.nim`); the gases differ in which of the variants below they use.
The per-gas modules in `gases/` reproduce Fortran bit for bit. Once all gases
are ported, these columns become the fields of one compile-time configuration
and a single generic routine replaces the modules.

Legend: `x` used, `-` not used, numbers are counts or factors.

## Elastic

| gas | interpolation | low-energy model | anisotropy (NANISO=1 / 2) |
|---|---|---|---|
| 1 CF4 | linear | - | PQ1 from QELA-QMOM / 1-EPS table |
| 2 Ar | linear | phase shifts (LMAX 100, single precision) < 1 eV | PQ1 / 1-EPS |
| 3,4 He | linear | - | PQ1 / 1-EPS |
| 5 Ne | linear | phase shifts < 1 eV | PQ1 / 1-EPS |
| 6 Kr | log-log | constant at EN=0 | PQ1 / 1-EPS |
| 7 Xe | log-log | constant at EN=0 | PQ1 / 1-EPS |
| 8 CH4 | log-log | - | PQ1 / 1-EPS |
| 9 C2H6, 10 C3H8 | linear < 10 eV, log-log above | - | PQ1 / 1-EPS |
| 11 iC4H10 | linear | - | PQ1 / 1-EPS |
| 12 CO2 | log-log | ground/bend-excited mix (APBEND); minus rotations | PQ1 / 1-EPS, isotropic < 10 eV |
| 13 neo-C5H12 | log-log (momentum transfer only) | - | isotropic only |
| 14 H2O | log-log (total, without rotation) | - | always anisotropic: 1-EPS (KEL=KIN=2 for any NANISO) |
| 15 O2 | linear (total, momentum transfer, eps) | - | PQ1 / 1-EPS |
| 16 N2 | linear (total, momentum transfer, eps), rotations subtracted | - | PQ1 / 1-EPS |
| 17 NO | linear (momentum transfer) | - | isotropic only |
| 18 N2O | log-log (total, momentum transfer), eps linear | - | 1-EPS only (NANISO=1 becomes 2); KEL=KIN=2 for any NANISO |
| 19 C2H4 | log-log | - | always 1-EPS (sets NANISO=2) |
| 20 C2H2 | linear (momentum transfer) | - | isotropic only |
| 21 H2 | log-log (total, momentum transfer), eps linear | constant below XELM(2) | PQ1 / 1-EPS |
| 22 D2 | as H2, scaled by AMPEL 0.988 | constant below XELM(2) | PQ1 / 1-EPS |
| 23 CO | log-log, eps linear | constant below XEL(1) | always 1-EPS (forces NANISO=2) |
| 24 methylal | linear (momentum transfer) | - | isotropic only |
| 25 DME | log-log, eps linear | constant below XEL(1) | always 1-EPS (forces NANISO=2) |
| 26 Reid step | constant 1e-16 | - | isotropic only |
| 27 Maxwell | 6e-16/sqrt(E) | - | isotropic only |

## Ionisation

| gas | channels | high-energy tail | inner shells | renormalisation |
|---|---|---|---|---|
| 1 CF4 | 12 | Born-Bethe x fraction | C-K, F-K (x4) | CRFAC=(tot-shells)/tot on 1..10 |
| 2 Ar | 3 charge + 4 shells | Born-Bethe x .9466/.04448/.00987 | L3 L2 L1 K | QCORR |
| 3,4 He | 2 | Born-Bethe | - | split from total (QTEMP) |
| 5 Ne | 3 charge + K | Born-Bethe x .9594/.0388/.00215 | K | QCORR |
| 6 Kr | 4 charge + 7 shells | Born-Bethe x .9009/.0613/.0291/.0082 | M3..K | QCORR=(count-shells)/count |
| 7 Xe | 3 charge + Q456 + shells | Born-Bethe, nested Q4..6 | M5..K | QCORR |
| 8 CH4 | 8 fragments + K | Born-Bethe x fraction | C-K | (sum-K)/sum |
| 9 C2H6 | 15 fragments + K (x2) | Born-Bethe x fraction | C-K | (sum-K)/sum |
| 10 C3H8 | 23 fragments + K (x3) | Born-Bethe x fraction | C-K | (sum-K)/sum |
| 11 iC4H10 | 3 | Born-Bethe; ch2 = dipole/(EN+EION+E3) | C-K (x4) | ch2, ch3 subtracted from ch1 |
| 12 CO2 | 9 + 2 K | Born-Bethe x fraction | C-K, O-K (x2) | ch2, ch3 subtracted from ch1 |
| 13 neo-C5H12 | 1 in Q(3), table only | - | - | - |
| 14 H2O | 8 + O-K | counting ionisation (table, Born-Bethe > 5 keV) x fraction above each table | O-K | - (a channel below its table end skips the rest) |
| 15 O2 | 3 O2+ (one table, shifted) + 4 dissociative + K | Born-Bethe x .6475 x (.558/.308/.136), .2993/.0446/.0061/.0025 | O-K (x2) | O2+ X minus A, B and 2 K |
| 16 N2 | 7 N2+ states + 4 dissociative + K | N2+ table, fractions by channels open; Born-Bethe x .7973/.197/.0338/.0057 | N-K (x2) | N2+ X minus 1.99 K; N+ minus excited dissociative |
| 17 NO | 1 in Q(3), table only | - | - | - |
| 18 N2O | 9 fragments + N-K (x2) + O-K | Born-Bethe > 1 keV x .4692/.2088/.0623/.0135/.0498/.0201/.0399/.0932/.0432 | N-K, O-K | - (stops at the first channel below its table) |
| 19 C2H4 | 12 fragments + C-K (x2) | counting ionisation (log-log table, Born-Bethe > 3 keV) x fraction above each table | C-K | (sum-K)/sum |
| 20 C2H2 | 1 in Q(3), table only | - | - | - |
| 21 H2 | autoionisation + 18 H2+ vibrational states + dissociative | tables to 3 keV; above, shares AMP of a total table (to 16 keV), then Born-Bethe | - | - (stops at the first closed channel) |
| 22 D2 | autoionisation + 26 D2+ vibrational states + dissociative | as H2 (28 channels) | - | - (stops at the first closed channel) |
| 23 CO | CO+ X/A/B, C+, O+, CO++, C++, O++ (log-log tables) | shares of counting ionisation (x .78654 for the CO+ states), Born-Bethe > 5 keV | C, O K (linear tables, from 285/532 eV) | - (Q(3) = counting, left unset outside 14 eV..5 keV) |
| 24 methylal | 1 in Q(3), table only | - | - | - |
| 25 DME | 25 fragments (log-log tables to 1 keV) | shares of counting ionisation, Born-Bethe > 1 keV | C K (x 2e-16), O K | - |
| 26 Reid step | none (E(3)=15.9 unused) | - | - | - |
| 27 Maxwell | none | - | - | - |

## Attachment, inelastic, other

| gas | att | superelastic vib (detailed balance) | dipole vib (Born) | excitation forms | PEQIN gate | Penning |
|---|---|---|---|---|---|---|
| 1 CF4 | 1 | V1-V4 (APOP/DEG) | x (ADIP, RAT) | 3 tables falloff 1/E^2 + 33 Born | 3 EIN (dipole <= 100 eV) | - |
| 2 Ar | - | - | - | tables 1/E^n + Born | 2 EIN | all levels (.2, 1, 1) |
| 3,4 He | - | - | - | 18 tables (1/E,1/E^3,1/E^4, Born) + 31 scaled from n=4 | 2 EIN | all levels (1, 1, 1) |
| 5 Ne | - | - | - | tables (X/E)^1.5 or ^2 + Born | 2 EIN | all levels (.5, 1, 1) |
| 6 Kr | - | - | - | tables 1/E,1/E^3 + Born, group factors AN* | 2 EIN | all levels (0, 1, 1) |
| 7 Xe | - | - | - | tables 1/E^1,1.5,2,3 + 29 Born | 2 EIN | all levels (0, 1, 1) |
| 8 CH4 | 2 (H- x .02) | V4, V2 | V4, V3 (.076) | tables + Born | 3 EIN | - |
| 9 C2H6 | 2 | torsion + 6 modes | torsion (.0045) | 3 triplet tables + 45 Born | 3 EIN (<=34), 2 EIN | - |
| 10 C3H8 | 2 | 8 modes | j1/2 (.00536) | 4 triplet tables + Born | 3 EIN (<=32, 62), 2 EIN | - |
| 11 iC4H10 | 1 | ladder of 10 torsion states + pairs | torsion (.009) | 2 triplets 1/E^3 + 13 Born x ASING | 3 EIN / 2 EIN | 10..24 |
| 12 CO2 | 1 (1/E^3) | rotations J->J+-2 (Born quadrupole, 6/E > 6 eV) + (010) (020) (100) (030) (001) | (010), (001) | 9 polyad tables 1/E + 2 triplets 1/E^2 + 63 Born | 4 EIN rot, 3 EIN vib, 2 EIN diss | 82..NIN (0, 1, 1) |
| 13 neo-C5H12 | 1 (table) | 3 pairs | - | 2 vib + 3 exc tables, no tail | - | - |
| 14 H2O | 3 (log-log, 1e-21) | 105 rotation pairs (Born dipole, para/ortho populations, 0 > 2 keV) + (010) bend | - | vib tables 1/E, 1/E^1.5; 4 triplets log-log 1/E^1.5; 32 Born dipoles | rotations, vib, dipoles: dipole angular table (E/EIN), PEQEL > 500 EIN; triplets 2 EIN | 215..NIN (0, 1, 1), printed although 0 |
| 15 O2 | 2-body lin, 1/E^3 above; 3-body lin x 273.15 P/(760 (T+273.15)) | 24+24 rotation pairs (log-log per pair, Phelps resonance, x3/E > 3 eV, x0.75) + V1 | V1 (APOP2) | 20 V1 harmonics, a1/b1 (1/E^2), Herzberg, 4 triplet tables (1/E, 1/E^1.5), 69 Born | rotations, vib 3/6 eV; 59 3 EIN; Herzberg, dissoc 2 EIN | 76..NIN (0, 1, 1), printed although 0 |
| 16 N2 | - (zero) | 38+38 rotations (Born quadrupole, para/ortho 6/3, resonance, x QMOM/8.9e-16 > 5 eV; subtracted from Q(2)) + V1 | V1 (APOPV1) | 14 V1 harmonics; 14 triplet + 15 singlet tables (1/E^2, 1/E^1.5 or Born above) with momentum ratio; 4 Born + 2 Born x (E+2EIN)/E above E(3) | 1..92: 3 EIN (NANISO>0); states: 1.5-ratio (1) / PEQEL (2) > 3 EIN; 122-127 copy 121 | 107..NIN (0, 1, 1), printed A15 |
| 17 NO | 2-body lin (1e-18) + 3-body lin x 273.15 P/(760 (T+273.15)), in QATT only (Q(4)=0) | - | - | rot, 2 vib, 1 exc tables, no tail | - | - |
| 18 N2O | 1 (lin x 0.67 x 0.1) | 60+60 rotations (Born dipole DBK, J populations, 0 > 500 eV, angular tables) + (010) (020) (100) | (010), (100), (001) + 2 resonances | 10 vib tables 1/E; 4 triplets log-log 1/E^1.5; 36 Born dipoles (EN+E3+EIN) | vib >= 6 eV; triplets 2 EIN; dipoles angular table, PEQEL > 500 EIN | 129..NIN (0, 1, 1), all printed |
| 19 C2H4 | 1 (lin 1e-16) | V7, V3, V2 (APOP/(1+APOP)) | V7 V3 V2 V1 + 1.85 eV resonance (x amplitude per mode) | 4 resonance-only vib (log-log), 2V1 table; 2 triplets log-log + 2 lin (1/E^2 x amplitude); 36 Born dipoles | vib > 20 eV; triplets > 50/100 eV; dipoles 3 EIN, angular table 1..20 EIN | - |
| 20 C2H2 | 1 (table, 1e-21) | V5 (half the detailed balance) | - | 2V5, V3, 3 exc tables, no tail | - | - |
| 21 H2 | dissociative: exp. resonances from J = 0..7 (FROT-weighted) + table + 14 eV resonance | rotations J->J+2 (8 tables, ortho/para, ratio frozen at 3 below 273 K), de-excitation x (2J+1)/(2J+5) in single precision, 0 > 100 eV | - | 4 vib tables 1/E; B3Σ tables (3 cut at table end, 1 1/E^3); c3Π a3Σ (shares) e3Σ 1/E^3; EF table x .4/.6 then Born; 84 Born dipoles (x1.08, x1.20 groups) | rot/vib > 8 eV; then 2 EIN for 19..114 | none printed (PENFRA 0) |
| 22 D2 | 1 (table, 1e-21, linear extrapolation) | rotations as H2 (para J odd weight 3, ortho 6; ratio frozen at 2 below 273 K) | - | as H2: 4 vib, B3Σ, c3Π a3Σ e3Σ, EF, 111 Born dipoles (+ level 105 special) | rot/vib > 8 eV; then 2 EIN for 19..142 | none printed (PENFRA 0) |
| 23 CO | 1 (table 1e-16, inside the table) | rot J->J±1 Born dipole (tabulated ang. dist.) and J->J±2 quadrupole (< 2 eV, isotropic), all < 100 eV; V1 superelastic | - | 10 vib log-log; a3Π 6 tables (< 1 keV) and 6 more triplet tables 1/E^1.5 (1e-18); 70 Born dipoles (own A0 in BBCONST) | vib >= 4 eV; triplets always; dipoles 2 EIN | 61..NIN all printed (0, 1, 1) |
| 24 methylal | none | - | - | 3 vib, 2 exc tables, linear extrapolation | - | - |
| 25 DME | 1 (table 1e-20, open interval) | rot J->J+1 Born dipole (tabulated ang. dist.) x AMPROT .105 < 1 keV; 3 vib superelastic | 4 Born vib + resonance table | 6 triplets log-log 1/E^2; 44 Born dipoles (own A0) with ang. dist. to 500 EIN | vib >= 2 eV; triplets 2 EIN; dipoles always | 61..NIN all printed (0, 1, 1) |
| 26 Reid step | - | - | - | one level 0.2 eV, 1e-15, PEQIN .5+.3, 70 % in Q(1) | - | - |
| 27 Maxwell | - | - | - | none | - | - |

All gases: EFINAL cuts NIN at the first level above EFINAL (searched from a
gas-specific first level). Bremsstrahlung levels (C2H6, CO2, ...) are computed
by Magboltz but lie beyond NIN, so they are never used; the port skips them.
Older routines (GAS13, ...) step their own energies (EN = EN + ESTEP) instead
of using the grid, put the ionisation in Q(3) and ignore NANISO.
Some gases keep a minimum NIN after the EFINAL cut (O2 52, N2 77, C2H4 12,
N2O 120, CO 113); H2 and D2 set NIN to 19, 23 or 114/142 by EFINAL (8 and 11 eV).: the
vibrational ground levels stay. Three-body attachment (O2) needs pressure and
temperature, so GasInput carries TORR and TEMPC (also NO).
Some routines overwrite NANISO in /INPT/ (C2H4, CO, DME always 2, N2O 2 unless 0, and
more in the legacy gases); SETUP already sets NANISO=2 for every run, so in
Magboltz this only matters to direct calls. The modules override their own copy.
CO and DME use A0=0.5291772083e-8 (others 0.52917720859e-8) before computing
BBCONST, so their Born dipoles pass their own constant (`dipole(..., bb)`).

## Upstream quirks kept bit for bit (to report as MB-44+)

- GAS1 CF4 superelastic V4/V1/V3 tail above the table: `Y*(X/(EN*(EN-EIN))**2)`,
  i.e. Y*X/(E*E')^2, drops ~1e-5 at the table end (table branch is detailed balance).
- GAS1 CF4 level 46 uses F=0.00198, description says 0.00189.
- GAS4 He3 DATA Y41F(53)=0.0007827; He4 has 0.0000727 (neighbours .0000786, .0000678): typo.
- GAS7 Xe `3710 QIN(2,I)=Y4DSUM(N4DSUM)*...` should be QIN(23,I): above the 4D-sum
  table, level 23 is 0 and level 2 (1S4 resonance) is overwritten by the 4D tail.
- GAS7 Xe levels 25-27 (2P3/2P2/2P1) never copy PEQEL: isotropic, unlike the other 47.
- GAS8 CH4 level 33 EIN=20.5, description 20.25.
- GAS9 C2H6 level 20 gated on EIN(21) instead of EIN(20).
- GAS11 iC4H10 `IF(NANISO.EQ.2)PEQIN(2,I)=0.0` should be PEQION(2,I).
- GAS14 H2O SUMTRP sums QIN(230) (a dipole) instead of QIN(241) (higher triplets); the sums are not used, no effect.
- GAS15 O2 level 59 (b1Σg+) sets PEQIN(61,I) instead of PEQIN(59,I): level 59 stays
  isotropic, level 61 gets level 59's distribution between 3 EIN(59) and 6 eV.
- GAS15 O2 level 72 (V1 sum 21-30) takes its tail above the table from YVIB20, not YVIB21.
- GAS16 N2 level 119 takes the ratio's energy above its table from XSNG12(NSNG13), not XSNG13
  (both end at 100 eV, no effect).
- GAS16 N2 QION(9), QION(10): at exactly EN = EION(10) both stay 0, and QION(10) at exactly
  120 eV (`.LT.120` / `.GT.120`).
- GAS18 N2O level 156 (high triplets) gates its angular distribution on 2 EIN(147), not 2 EIN(156).
- GAS18 N2O level 134 (300) tests the table end XV101(NV300) instead of XV300(NV300)
  (both 4.25 eV, no effect); `QELM=YELMT(1)` for QELMT below XEN(1) (1e-9 eV, never reached).
- GAS19 C2H4 PEQEL(2)=0.0 for NANISO=0 (0.5 elsewhere); unreachable, the routine sets NANISO=2.
  Description of the H+ channel says ELOSS=18.46, EION is 18.66.
- GAS20 C2H2 cuts the last level with `EFINAL.LT.EIN(7)`, the others with `.LE.`.
- GAS21 H2 level 98 (D band) starts at EIN(91), not EIN(98) (`727 IF(EN.LE.EIN(91))`).
- GAS22 D2 levels 115 and 132 start at EIN(114) and EIN(131) (`744 IF(EN.LE.EIN(114))`, `761`).
- GAS22 D2 level 105 uses EIN(104) in the logarithm and in the denominator (f=.003391).
- GAS23 CO vibration 115 (3V1) has the angular-distribution test inverted:
  `IF(EN.GT.4.0) GO TO 425`, anisotropic only below 4 eV (the other vibrations `.LT.`).
- GAS23 CO level 199 starts at EIN(194) (`IF(EN.LE.EIN(194)) GO TO 825`).
- GAS23 CO level 128 (a3Π v=5-10) never copies PEQEL (isotropic), unlike 123-127.
- GAS23 CO Q(3) (counting ionisation) is not set below XION(1) or above 5 keV: stale values.
- GAS25 DME resonance of level 57 uses `AMPV3R*(EN*EIN(58))` where `(EN+EIN(58))` is meant
  (cf. level 55).
- GAS25 DME rotational populations weight (2J-1), CO (2J+1).
- GAS25 DME total Q(1) leaves out the excitations 61..110 (SUMINEL has rot, vib, att, ion only).
- GAS25 DME level 58 gate `EN.LT.EIN(58)`, the others `.LE.`.
