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

## Upstream defects and open questions

The irregularities found while writing these modules are reported in the
Magboltz fork (`~/src/Magboltz`, branch `fixes`): `bugs/README.md` MB-44 …
MB-56 are fixed in the fork and here in the same way, so the bitwise
comparison holds against `fixes`. The cases whose intent is unclear are
kept as in Magboltz 11.21 and listed in `bugs/QUESTIONS.md` (Q-01 … Q-13:
descriptions that disagree with the code, CF4 superelastic tails, C2F6 QION(5)
share, C3F8 K-shells, CS2 three-body attachment, CD4 EOBY, C4F7N E(3) and
nitrogen K-shell, HFO-1234ze Penning loop, Hg2 GFRAC, the totals Q(1)).
Harmless inconsistencies (`.LT.`/`.LE.` thresholds, unread quantities) are
in the README's "Inconsistencies left unchanged".
