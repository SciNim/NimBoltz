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

All gases: EFINAL cuts NIN at the first level above EFINAL (searched from a
gas-specific first level). Bremsstrahlung levels (C2H6, CO2, ...) are computed
by Magboltz but lie beyond NIN, so they are never used; the port skips them.
Older routines (GAS13, ...) step their own energies (EN = EN + ESTEP) instead
of using the grid, put the ionisation in Q(3) and ignore NANISO.
Some gases keep a minimum NIN after the EFINAL cut (O2 52, N2 77): the
vibrational ground levels stay. Three-body attachment (O2) needs pressure and
temperature, so GasInput carries TORR and TEMPC.

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
