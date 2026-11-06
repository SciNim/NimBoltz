# Gas cross sections: feature matrix

Every gas routine (GAS1..GAS80 of Magboltz 11.21) computes the same things on the
energy grid: elastic (total, momentum transfer, angular parameter), ionisation
channels, attachment, inelastic levels, null collisions. The physics is shared
(`gases/common.nim`); the gases differ in which of the variants below they use.
The per-gas modules in `gases/` (all 80 routines; 62-80 are one `dummy`)
reproduce Fortran bit for bit. Next, these columns become the fields of one
compile-time configuration and a single generic routine replaces the modules.

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
| 28 Reid ramp | constant 6e-16 | - | isotropic only |
| 29 C2F6 | log-log, eps linear | constant below XEN(1) | always 1-EPS (forces NANISO=2) |
| 30 SF6 | log-log (total, MT, 1-eps) | 1434e-16, eps 0 below XEN(2) | .5+(QEL-QMT)/QEL / 1-EPS |
| 31 NH3 | log-log, eps linear | constant below XEL(1) | always 1-EPS (forces NANISO=2) |
| 32 propene | linear (momentum transfer) | - | isotropic only |
| 33 cyclopropane | linear (momentum transfer) | - | isotropic only |
| 34 methanol | log-log, eps linear | constant below XEL(1) | always 1-EPS (forces NANISO=2) |
| 35 ethanol | as methanol, x ANORM (1.0) | constant below XEL(1) | always 1-EPS (forces NANISO=2) |
| 36 iso-propanol | as ethanol | constant below XEL(1) | always 1-EPS (forces NANISO=2) |
| 37 Cs | linear | - | isotropic (KEL 0) |
| 38 F2 | linear | - | isotropic |
| 39 CS2 | linear | - | isotropic |
| 40 COS | log-log | 19e-16 at 0 | isotropic |
| 41 CD4 | log-log | constant below XEN(2) | PQ1 |
| 42 BF3 | linear | - | isotropic |
| 43 C2H2F4 | log-log, eps | - | always 1-EPS (forces NANISO=2) |
| 44 TMA | log-log, eps linear | - | PQ1 / 1-EPS |
| 45 C3F4H2 | log-log, eps linear | - | 1-EPS (forces NANISO=2) |
| 46 n-propanol | as iso-propanol | constant below XEL(1) | always 1-EPS |
| 47 para-H2 | as H2 | | |
| 48 ortho-D2 | as D2 | | |
| 49 C4F7N | as C3F4H2 | - | 1-EPS (forces NANISO=2) |
| 50 CHF3 | log-log | 16300e-16 at 0 | isotropic |
| 51 CF3Br | log-log | 2800e-16 at 0 | isotropic |
| 52 C3F8 | log-log, eps | - | 1-EPS (forces NANISO=2) |
| 53 O3 | log-log | 3000e-16 at 0 | isotropic |
| 54 Hg2 | log-log (MT; total unused) | - | isotropic |
| 55 H2S | log-log from XEL(2) | YEL(2) below | isotropic; rotations R1/R2 subtracted from Q(2) |
| 56 n-butane | linear | - | isotropic |
| 57 n-pentane | log-log | constant below XEN(2) | isotropic |
| 58 N2 (Phelps) | linear | - | PQ1 |
| 59 GeH4 | log-log | - | PEQEL 0.5 |
| 60 SiH4 | log-log (total, MT, eps) | - | always 1-EPS (forces NANISO=2) |
| 61 CCl4 | log-log (total, MT, eps) | constant below XEN(1), eps 0 | always 1-EPS (forces NANISO=2) |
| 62-80 dummy | zero | - | - |

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
| 28 Reid ramp | none | - | - | - |
| 29 C2F6 | 10 fragments (log-log tables to 900 eV) | shares of counting ionisation to 3 keV, Born-Bethe above | C K (x 2), F K (x 6), linear | channels 1-10 x (tot-K)/tot |
| 30 SF6 | 7 outer (linear tables to XION(N)) | shares of Born-Bethe | S L3, L2, L1, K, F K (x 6) | outer x (tot-inner)/tot |
| 31 NH3 | 9 (linear tables) | shares of counting ionisation (linear to 1 keV, log-log to 40 keV), Born-Bethe above | N K | - |
| 32 propene | 1 in Q(3), table only | - | - | - |
| 33 cyclopropane | 1 in Q(3), table only | - | - | - |
| 34 methanol | 15 fragments (log-log tables) | shares of counting ionisation, Born-Bethe > 1 keV | C K, O K | - |
| 35 ethanol | 25 fragments (log-log tables) | shares of counting ionisation, Born-Bethe > 1 keV | C K (x 2e-16), O K | - |
| 36 iso-propanol | 27 fragments (log-log tables) | shares of counting ionisation, Born-Bethe > 1 keV | C K, O K | - |
| 37 Cs | 1 linear | - | - | - |
| 38 F2 | 1 linear | - | - | - |
| 39 CS2 | 1 linear | - | - | - |
| 40 COS | 1 linear | - | - | - |
| 41 CD4 | 1 linear | - | - | - |
| 42 BF3 | 1 linear | - | - | - |
| 43 C2H2F4 | 1 log-log (counting x AMPION) | Born-Bethe x 1.0518 > 3 keV | C K, F K | CRFAC on 1 |
| 44 TMA | 1 log-log | Born-Bethe | C K, N K | AUGK subtraction (inside the 2 EION branch) |
| 45 C3F4H2 | 2 log-log x 1.2 | Born-Bethe split .9484/.0516 > 3 keV | C K x3, F K x4 (linear) | CRFAC on 1, 2 |
| 46 n-propanol | 28 fragments | as iso-propanol | C K, O K | - |
| 47 para-H2 | as H2 | | | |
| 48 ortho-D2 | as D2 | | | |
| 49 C4F7N | 2 log-log | Born-Bethe x 1.0486 split .9484/.0516 | C K x3, F K x4, N K (unused) | CRFAC on 1, 2 |
| 50 CHF3 | 1 linear | - | - | - |
| 51 CF3Br | 1 | - | - | - |
| 52 C3F8 | 7 log-log (shares of QIONC) | Born-Bethe x 1.052 | C K, F K (linear, not subtracted) | - |
| 53 O3 | 1 linear | - | - | - |
| 54 Hg2 | 1 linear | - | - | - |
| 55 H2S | 1 | - | - | - |
| 56 n-butane | 1 linear | - | - | - |
| 57 n-pentane | 1 linear | - | - | - |
| 58 N2 (Phelps) | 1 linear | - | - | - |
| 59 GeH4 | 1 | - | - | - |
| 60 SiH4 | 6 linear tables to 100 eV | shares of counting ionisation (log-log), Born-Bethe > 1 keV | Si K | - |
| 61 CCl4 | 9 log-log tables to 1 keV | shares of Born-Bethe > 1 keV | C K, Cl K (x4) | CRFAC on 1..9 |
| 62-80 dummy | 1, zero | - | - | - |

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
| 28 Reid ramp | - | - | - | one level 0.2 eV, (E-0.2) x 1e-15 | - | - |
| 29 C2F6 | F-, CF3- (one log-log table x .574, split .77/.23) | V6, V5, V7 Born, populations of 5 quanta | V6, V5, V7 Born + resonance table (log-log) | 2V7, VHIGH resonance only; 4 triplets linear 1/E^2; 6 Born dipoles | vib Okhrimovskyy (super ENRTS, inel ENROT), resonances 1-EPSRES; exc 3 EIN | none |
| 30 SF6 | 7 ions (SF6- log-log, rest linear), Q(4)=0 | V4, V1, V3 (deg 3, 1, 3) | V4 Born; V3 Born + table; PEQIN .5+(Q-XMT)/Q | V1 harmonics tables 1/E^2; 3 triplets; 22 Born dipoles x 1.013-1.025 | harmonics always; exc 3 EIN | none |
| 31 NH3 | 2 linear (H-, NH2-) | V2 (APOP/(1+APOP)) Born + resonance | V2, V4, V1+V3 Born + resonance; harmonics resonance only | 28 rotations J->J±1 of the symmetric top (Born dipole, K sums with Boltzmann populations, ang. dist., <= 1 keV); 6 triplets log-log 1/E^2; 64 Born dipoles, ang. dist. to 500 EIN | 2 EIN | 38..NIN all printed (0, 1, 1) |
| 32 propene | 1 linear 1e-21 | V20 Born, V13 and V Born + table | same | harmonics, 3 exc tables | - | - |
| 33 cyclopropane | 1 linear 1e-21 | V11 and V Born + table | same | 2V3, harmonics, 2 exc tables | - | - |
| 34 methanol | 3 (H-, O-/OH-, CH3O-) linear 1e-20 | rot J->J+1 Born dipole (tabulated ang. dist.) x AMPROT .1 < 1 keV; V1 | V1, V2, V3 Born + resonance table | 6 triplets log-log 1/E^2; 22 Born dipoles (own A0) with ang. dist. to 500 EIN | vib always; triplets, ions 1.5 EIN | 57..NIN all printed (0, 1, 1) |
| 35 ethanol | 3 (H-, O-/OH-, C2H5O-) x AMPATT 1.0 | as methanol | as methanol | 6 triplets; 31 Born dipoles (own A0) | vib >= 4 eV; triplets, ions 2 EIN | 57..NIN all printed (0, 1, 1) |
| 36 iso-propanol | 3 | as ethanol | as ethanol | 6 triplets; Born dipoles; bremsstrahlung (EFINAL gate) | 2 EIN | as ethanol |
| 37 Cs | 1 | - | - | 5 linear tables | - | - |
| 38 F2 | 1 | rotation/vib APOP | vib tables | 4 linear tables | - | - |
| 39 CS2 | 3-body (x1.3, fixed) + 2-body | V2, V1 (APOP) | - | tables | - | - |
| 40 COS | 1 x 1e-16 | V2, V1, 2V2 (APOP) | dipole V3 | 1 exc | - | - |
| 41 CD4 | - | - | V2+V4, V1+V3 dipole + table | 4 dissociation tables | 2 E(3) | - |
| 42 BF3 | 1 | - | - | 9 levels linear | - | - |
| 43 C2H2F4 | 3-body (density) + 2-body log-log | 4 (APOP, dipole + resonance) | x | 4 triplets 1/E^2, 6 Born dipoles; 3 bremsstrahlung levels > 1 keV | 2 EIN | - |
| 44 TMA | 1 | torsion, V1 (APOP) | x | 3 triplets, 25 Born dipoles (per-level factors) | 3/4/5 EIN | - |
| 45 C3F4H2 | 5 (thermal dimers, 3 x 2-body, 3-body density) | 6 (APOP) Born + L/H resonances | x | 4 triplets, 10 Born dipoles, 2 B3PI; 3 bremsstrahlung > 1 keV; debug print 9011 | 2/3 EIN | - |
| 46 n-propanol | 3 | as iso-propanol | | | | |
| 47 para-H2 | as H2, odd-J population 0 | | | | | |
| 48 ortho-D2 | as D2, odd-J population 0 | | | | | |
| 49 C4F7N | 2 (resonance, 2-body) | as C3F4H2 | x | 5 triplets, 8 Born dipoles; bremsstrahlung N/C/F > 1 keV | 2/3 EIN, 50 eV | - |
| 50 CHF3 | 1 x 1e-19 | rotation/V6, V3/V1 (APOP) | x | vib tables, exc/dissociation | - | - |
| 51 CF3Br | 1 | APOP | x | tables | - | - |
| 52 C3F8 | resonance + 3-body (1/(31/TORR+1)) | 3 (APOP) Born + resonances | x | 4 triplets linear 1/E^2, Born dissociation | 2 EIN | - |
| 53 O3 | 2-body | rotation/V (APOP) dipoles | x | 3 exc linear | - | - |
| 54 Hg2 | 2 x density (GFRAC, T, TORR) | - | - | vib .04 + excitations linear | - | - |
| 55 H2S | 1 | rotations R1, R2 | V2, V13 dipole | tables | - | - |
| 56 n-butane | 1 x 1e-21 | 2 pairs (APOP) | - | 3 vib, 3 exc linear | - | - |
| 57 n-pentane | 1 x 1e-21 | 3 pairs (APOP) | - | 2 vib, 3 exc | - | - |
| 58 N2 (Phelps) | - | rotations J±2 (ortho/para Boltzmann) | - | 8 vib, 7 exc tables | 2 E(3) | - |
| 59 GeH4 | 1 x 1e-19 | V2+V4, V1+V3 (APOP) dipole + log-log resonance | x | 1 dissociation | - | - |
| 60 SiH4 | 4 linear 1e-18 | V4, V2 (APOP/DEG) Born + resonance | V13, harmonics | 3 triplets, 38 Born dipoles (ang. dist. 1-20 EIN); 2 bremsstrahlung beyond NIN | 3 EIN | - |
| 61 CCl4 | 3 log-log | V2, V4, V1, V3 (APOP/DEG) tables, V4/V3 + Born dipole | x | harmonics, 3 triplets linear 1/E^2 x 6.12, 54 Born dipoles; 2 bremsstrahlung levels (NIN 68 if EFINAL > 1 keV) | 2 EIN | - |
| 62-80 dummy | 1, zero | - | - | - | - | - |

All gases: EFINAL cuts NIN at the first level above EFINAL (searched from a
gas-specific first level). Bremsstrahlung levels (C2H6, CO2, ...) are computed
by Magboltz but lie beyond NIN, so they are never used; the port skips them.
Older routines (GAS13, ...) step their own energies (EN = EN + ESTEP) instead
of using the grid, put the ionisation in Q(3) and ignore NANISO.
Some gases keep a minimum NIN after the EFINAL cut (O2 52, N2 77, C2H4 12,
N2O 120, CO 113); H2 and D2 set NIN to 19, 23 or 114/142 by EFINAL (8 and 11 eV).: the
vibrational ground levels stay. Three-body attachment (O2) needs pressure and
temperature, so GasInput carries TORR and TEMPC (also NO).
Some routines overwrite NANISO in /INPT/ (C2H4, CO, DME, C2F6, methanol, ethanol always 2, N2O 2 unless 0, and
more in the legacy gases); SETUP already sets NANISO=2 for every run, so in
Magboltz this only matters to direct calls. The modules override their own copy.
CO, DME, methanol and ethanol use A0=0.5291772083e-8 (others 0.52917720859e-8) before computing
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
- GAS29 C2F6 carbon K-shell: `IF(EN.LE.XION11(J)) GO TO 151` jumps to the loop's own
  CONTINUE (152 meant), so QION(11) always extrapolates the table's last interval.
- GAS29 C2F6 QION(5) (C2F4+) share of QIONC is .0045 to 3 keV, .00045 in the Born-Bethe branch.
- GAS29 C2F6 `DO 395 L=4,6` sets PEQIN(L) where PEQIN(2L-1) is meant: the inelastic V6/V5/V7
  (7, 9, 11) stay isotropic, levels 4, 5, 6 get their values; the gate tests EIN(2L).
- GAS29 C2F6 level 8 opens with EIN(7) (same energy); PEQIN = PEQEL above 20 eV only inside the
  VHIGH branch, which ends at XVIBR(NVIBR) = 12.5 eV: dead code. Q(1) is never set.
- GAS30 SF6 level 23 text says F=0.1211, the code DIP(0.01211, ...).
- GAS30 SF6 Q(4)=0 after the attachment (QATT only); PEQION opens at 2 EION inclusive for
  channels 6, 7 (`.LT.`), above it for the others; F K-shell gated `.LE.EION`, others `.LT.`.
- GAS30 SF6 APOPGS is normalised and then reset to 1.0 (excited states normalised, ground not).
- GAS35 ethanol: the sums are commented out, Q(1) is not set.
- GAS31 NH3 (MB-44 candidates, rotational populations): level 2 (J 0->1) weights PJ(2) (population of J=1) where the ground state PJ(1) is meant; level 3 (2->1) has PJ(5)*(AJG2+1.0) for K=1 (should be AJG2-1.0, cf. all other levels); level 4 (1->2) keeps AJ=2.0 from level 3 (AJG set twice, AJ not reset -> (2AJ+1)=5 not 3); PEQIN(32,I)=PEQIN(32,I-IOFFN(32)) copies PEQIN (always 0), not PEQEL; harmonics level 33 not gated by EIN(33) (resonance 4.5-10.5 eV only); rotations and their angular distributions skipped above 1 keV; Q(1) not set.
- GAS38 F2: SCRPT(2) ' ELASTIC       FLOURINE' (spelling, cosmetic).
- GAS39 CS2: QATT(1) is a three-body attachment "AT 40TORR" scaled by a fixed 1.3, not by the actual pressure (O2/NO scale with TORR/TEMPC).
- GAS41 CD4: EOBY(1)=7.3 while E(3)=EION=12.99 (Opal-Beaty split 0.56·EION; maybe intended).
- GAS43 C2H2F4: EFINAL>1000 adds 3 bremsstrahlung levels (NIN+3) whose EIN stay 0 and KIN unset (loop ran to the old NIN); Q(1) never set; K-shell C uses GO TO 151 quirk (last interval extrapolated) as GAS29; penning block commented out; FAC=(273.15+20.0)·TORR/((TEMPC+273.15)·760) (3-body).
- GAS44 TMA: AUGKC/AUGKN K-shell subtraction from QION(1) sits inside the `EN > 2*EION(3)` angular branch, so between 285 eV (C-K) / 401.6 eV (N-K) and 803 eV the K-shells are double counted; NAME not set for NANISO=1; Bethe uses DEN(1) not DEN(I) (DEN=0 anyway); bremsstrahlung QIN(38..40) computed but NIN never reaches them.
- GAS45 C3F4H2: active debug WRITE(6,9011) inside the energy loop prints 4000 CSV lines (EN,QELA,QMMT,QSUP,...) every call — still in fork44 (candidate MB-44: comment out). Penning loop DO 61 K=61,NIN never runs (NIN=29). QIONC/QIONG computed, unused; Q(1) never set; brems levels 30-32 EIN=0 when EFINAL>1000; QIN(1..12) zeroed above 1 keV.
- GAS46 n-propanol: SCRPT(111) 'ELOSS= 13.25' but EIN(75)=13.35 (and comment 13.35); XTRP4(4) 18.8 vs isopropanol 17.8.
- GAS49 C4F7N: NVIBL=NVIBH=19 while XVIBL/YVIBL have 21 and XVIBH/YVIBH 20 points (tails unused); mass E(2) uses 114.0416 (C3F4H2, copied; C4F7N is 195.04); N K-shell QION(5)/EION(5) beyond NION=4 (never used, not subtracted by CRFAC, only set above the F K-edge 685.4 eV, otherwise stale); SCRPT(4) says ELOSS=35.0 but EION(2)=38.0; debug WRITE 9011 per bin + two WRITE(6,*) headers AFTER the table (candidate MB-44 with GAS45).
- GAS52 C3F8: K-shells (C, F) not subtracted from the other channels (no CRFAC, unlike C2F6/C2H2F4) and C not ×3 / F not ×8; same DO 395 PEQIN(L) for PEQIN(2L-1) quirk as GAS29; SCRPT(2) says CARBON HEXAFLUORIDE (C3F8 is octafluoropropane).
- GAS53 OZONE: Q(4,I)=0.0 after the attachment block, so QATT(1..2) (O-, O2-) are neither in Q(4) nor in the total Q(1). MIXER reads QATT, so transport is fine; only Q(1)/printout of the total is low.
- GAS54 Hg dimer: NIN cut uses EFINAL.LT.EIN (others .LE.); QELA (XEL/YEL table) computed every bin and never used; GFRAC from the LAST slot with NGASN=54 (two slots of Hg = only one counted).
- GAS55 H2S: SCRPT(9) 'ROT R2 EL0SS' (zero for O); below XEL(2) elastic is YEL(2), not YEL(1); rotations R1/R2 (QIN1-4) subtracted from Q(2) (elastic was total incl. rotations — intended?).
- GAS60 SiH4: QIN(25) gated on EN>EIN(26) and QIN(27) on EIN(28) (next level's threshold); PEQIN(24)=PEQEL(2,I-IOFFN(25)) uses level 25's offset; brems QIN(50/51) beyond NIN.
- GAS61 CCl4 (MB-44 candidate, real effect): QIN(7,I)=QIN(7,1)*APOPV3/DEGV3*1.D-16 (index 1 for I) -> V3 superelastic = first-bin value *1e-16 again (~1e-33 cm2), detailed balance with V3 lost; V3 super tail uses (EN-EIN(8)) not (EN+EIN(8)); QIN(19) gated EN>EIN(18); active debug WRITE(6,392) ATRP every call; SCRPT(3) CCL3+ says ELOSS=13.35 but EION(1)=11.35; ion gates mix .LE./.LT.; Q(3) zero; brems levels 67/68 KIN unset.
