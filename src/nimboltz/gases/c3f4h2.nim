## HFO-1234ze (C3F4H2): GAS45 of Magboltz 11.21 ("C3F4H2 2024").
##
## Always anisotropic (the routine sets NANISO=2): elastic log-log with
## 1 - ε; single and double ionisation (log-log tables, Born-Bethe shares
## above 3 keV, all scaled by 1.2) less the carbon and fluorine K-shells;
## thermal attachment to dimers, three dissociative two-body channels and a
## density scaled three-body one. Six vibrations as Born dipoles plus low
## and high energy resonance shapes, superelastic and inelastic (cut above
## 1 keV), a moderation level, two B3Π triplet parts (log-log, falloff
## 1/E²), Born dipole singlets and dissociation; bremsstrahlung on H, C and
## F (levels 30..32, above 1 keV). No total cross section (Q(1) is not
## set).

import engine

const xel: array[131, float32] = [
  0.0001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0,
  5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0,
  80.0, 90.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 400.0, 500.0,
  600.0, 800.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 2500.0, 3000.0,
  3500.0, 4000.0, 4500.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0,
  12500.0, 15000.0, 17500.0, 20000.0, 25000.0, 30000.0, 35000.0, 40000.0,
  45000.0, 50000.0, 60000.0, 70000.0, 80000.0, 90000.0, 100000.0, 125000.0,
  150000.0, 175000.0, 200000.0, 250000.0, 300000.0, 350000.0, 400000.0,
  450000.0, 500000.0, 600000.0, 700000.0, 800000.0, 900000.0, 1000000.0,
  1250000.0, 1500000.0, 1750000.0, 2000000.0, 2500000.0, 3000000.0, 3500000.0,
  4000000.0, 4500000.0, 5000000.0, 6000000.0, 7000000.0, 8000000.0, 9000000.0,
  10000000.0, 12500000.0, 15000000.0, 17500000.0, 20000000.0, 25000000.0,
  30000000.0, 35000000.0, 40000000.0, 45000000.0, 50000000.0, 60000000.0,
  70000000.0, 80000000.0, 90000000.0, 100000000.0, 125000000.0, 150000000.0,
  175000000.0, 200000000.0, 250000000.0, 300000000.0, 350000000.0,
  400000000.0, 450000000.0, 500000000.0, 600000000.0, 700000000.0,
  800000000.0, 900000000.0, 1e+9]
const ymt: array[131, float64] = [
  69000.0, 6900.0, 690.0, 69.0, 34.5, 23.0, 17.399999618530273, 15.5, 15.0,
  14.699999809265137, 14.800000190734863, 14.899999618530273, 15.0,
  15.199999809265137, 15.300000190734863, 15.399999618530273, 15.5,
  15.600000381469727, 15.699999809265137, 15.699999809265137,
  14.800000190734863, 12.800000190734863, 10.199999809265137,
  8.949999809265137, 7.900000095367432, 6.949999809265137, 6.199999809265137,
  5.599999904632568, 5.150000095367432, 4.699999809265137, 3.759999990463257,
  3.319999933242798, 2.9800000190734863, 2.5999999046325684,
  2.0399999618530273, 1.7699999809265137, 1.2699999809265137,
  0.9950000047683716, 0.828000009059906, 0.5640000104904175,
  0.4090000092983246, 0.28299999237060547, 0.20900000631809235,
  0.16099999845027924, 0.12800000607967377, 0.08749999850988388,
  0.06379999965429306, 0.048700001090765, 0.03849999979138374,
  0.031300000846385956, 0.02590000070631504, 0.018699999898672104,
  0.0142000000923872, 0.01119999960064888, 0.009069999679923058,
  0.007499999832361937, 0.0050200000405311584, 0.0036100000143051147,
  0.002739999908953905, 0.00215000007301569, 0.00144, 0.00104, 0.000785,
  0.000618, 0.000501, 0.000415, 0.0003, 0.000229, 0.000181, 0.000147, 0.000123,
  0.0000836, 0.0000614, 0.0000474, 0.000038, 0.0000263, 0.0000196, 0.0000154,
  0.0000125, 0.0000104, 8.82e-6, 6.66e-6, 5.26e-6, 4.3e-6, 3.59e-6, 3.06e-6,
  2.19e-6, 1.65e-6, 1.3e-6, 1.05e-6, 7.37e-7, 5.49e-7, 4.26e-7, 3.41e-7,
  2.8e-7, 2.34e-7, 1.71e-7, 1.31e-7, 1.04e-7, 8.45e-8, 7.01e-8, 4.71e-8,
  3.39e-8, 2.56e-8, 2.01e-8, 1.33e-8, 9.51e-9, 7.14e-9, 5.56e-9, 4.45e-9,
  3.65e-9, 2.58e-9, 1.92e-9, 1.48e-9, 1.18e-9, 9.59e-10, 6.19e-10, 4.31e-10,
  3.18e-10, 2.44e-10, 1.56e-10, 1.09e-10, 7.98e-11, 6.11e-11, 4.83e-11,
  3.91e-11, 2.72e-11, 2e-11, 1.53e-11, 1.21e-11, 9.79e-12]
const yel: array[131, float32] = [
  69000.0, 6900.0, 690.0, 69.0, 34.5, 23.0, 17.4, 15.5, 15.0, 14.7, 14.8,
  14.9, 15.0, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.7, 14.8, 12.8, 10.2,
  8.95, 7.9, 6.95, 6.2, 5.6, 5.15, 8.41, 7.54, 7.4, 7.31, 6.93, 6.38, 6.32,
  5.64, 5.27, 5.08, 4.37, 3.81, 3.18, 2.73, 2.39, 2.13, 1.75, 1.49, 1.29,
  1.14, 1.03, 0.93, 0.785, 0.679, 0.599, 0.536, 0.485, 0.393, 0.331, 0.287,
  0.253, 0.206, 0.174, 0.152, 0.135, 0.122, 0.111, 0.0949, 0.0836, 0.076,
  0.0684, 0.063, 0.0535, 0.0471, 0.0426, 0.0392, 0.0345, 0.0315, 0.0293,
  0.0277, 0.0265, 0.0255, 0.0241, 0.0231, 0.0224, 0.0218, 0.0214, 0.0207,
  0.0203, 0.02, 0.0198, 0.0195, 0.0194, 0.0193, 0.0192, 0.0192, 0.0191,
  0.0191, 0.0191, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019,
  0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019,
  0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019,
  0.019, 0.019, 0.019]
const yeps: array[131, float64] = [
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
  0.39250001311302185, 0.3257000148296356, 0.27459999918937683,
  0.23549999296665192, 0.20479999482631683, 0.1606999933719635,
  0.13099999725818634, 0.09421999752521515, 0.07281000167131424,
  0.058959998190402985, 0.04238000139594078, 0.03285999968647957,
  0.025550000369548798, 0.02085000090301037, 0.017580000683665276,
  0.015180000104010105, 0.011909999884665012, 0.009782000444829464,
  0.00829199980944395, 0.007193000055849552, 0.006351000163704157,
  0.005681999959051609, 0.00469199987128377, 0.003994000144302845,
  0.003473999910056591, 0.0030739998910576105, 0.0027550000231713057,
  0.0021860001143068075, 0.0018100000452250242, 0.0015429999912157655,
  0.0013429999817162752, 0.0010659999679774046, 0.0008812, 0.0007499,
  0.0006516, 0.0005756, 0.0005148, 0.0004238, 0.0003589, 0.0003057, 0.0002728,
  0.0002428, 0.0001892, 0.0001537, 0.0001285, 0.0001098, 0.00008401,
  0.00006713, 0.00005528, 0.00004658, 0.00003994, 0.00003473, 0.00002712,
  0.00002188, 0.00001809, 0.00001524, 0.00001305, 9.368e-6, 7.043e-6, 5.5e-6,
  4.422e-6, 3.043e-6, 2.226e-6, 1.701e-6, 1.342e-6, 1.087e-6, 8.981e-7,
  6.433e-7, 4.836e-7, 3.768e-7, 3.02e-7, 2.473e-7, 1.617e-7, 1.139e-7,
  8.447e-8, 6.514e-8, 4.211e-8, 2.942e-8, 2.169e-8, 1.664e-8, 1.315e-8,
  1.065e-8, 7.374e-9, 5.394e-9, 4.107e-9, 3.225e-9, 2.596e-9, 1.636e-9,
  1.119e-9, 8.11e-10, 6.14e-10, 3.85e-10, 2.63e-10, 1.9e-10, 1.44e-10,
  1.12e-10, 9e-11, 6.2e-11, 4.5e-11, 3.4e-11, 2.6e-11, 2.1e-11]
const xion: array[47, float32] = [
  11.0, 13.7, 14.7, 15.7, 16.7, 17.7, 19.7, 21.7, 23.7, 25.7, 27.7, 29.7,
  31.7, 33.7, 35.7, 37.7, 39.7, 41.7, 43.7, 45.7, 47.7, 57.7, 67.7, 77.7,
  87.7, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0,
  500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0,
  2500.0, 3000.0]
const yionc: array[47, float32] = [
  0.00001, 0.006, 0.05, 0.11, 0.19, 0.28, 0.69, 1.16, 1.51, 1.95, 2.25, 2.65,
  3.0, 3.3, 3.65, 4.05, 4.34, 4.6, 4.84, 5.12, 5.41, 5.98, 6.36, 6.55, 6.76,
  6.88, 7.02, 7.07, 6.96, 6.83, 6.45, 6.2, 5.91, 5.59, 5.35, 5.0, 4.59, 4.12,
  3.78, 3.43, 3.23, 2.73, 2.39, 2.13, 1.95, 1.64, 1.44]
const yiong: array[47, float32] = [
  0.00001, 0.006, 0.05, 0.11, 0.19, 0.28, 0.69, 1.16, 1.51, 1.95, 2.25, 2.65,
  3.0, 3.3, 3.65, 4.05, 4.35, 4.61, 4.87, 5.15, 5.46, 6.08, 6.51, 6.76, 7.07,
  7.25, 7.44, 7.5, 7.38, 7.25, 6.82, 6.57, 6.26, 5.89, 5.64, 5.27, 4.84, 4.33,
  3.97, 3.61, 3.4, 2.88, 2.51, 2.24, 2.05, 1.72, 1.51]
const xion1: array[47, float32] = [
  11.0, 13.7, 14.7, 15.7, 16.7, 17.7, 19.7, 21.7, 23.7, 25.7, 27.7, 29.7,
  31.7, 33.7, 35.7, 37.7, 39.7, 41.7, 43.7, 45.7, 47.7, 57.7, 67.7, 77.7,
  87.7, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0,
  500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0,
  2500.0, 3000.0]
const yion1: array[47, float32] = [
  0.00001, 0.006, 0.05, 0.11, 0.19, 0.28, 0.69, 1.2, 1.6, 2.02, 2.35, 2.75,
  3.05, 3.35, 3.65, 4.05, 4.34, 4.59, 4.82, 5.08, 5.36, 5.89, 6.24, 6.34,
  6.39, 6.51, 6.61, 6.67, 6.56, 6.46, 6.05, 5.81, 5.56, 5.28, 5.04, 4.74,
  4.36, 3.9, 3.57, 3.25, 3.07, 2.59, 2.26, 2.02, 1.84, 1.55, 1.36]
const xion2: array[32, float32] = [
  38.0, 40.0, 42.0, 44.0, 46.0, 48.0, 58.0, 68.0, 78.0, 88.0, 100.0, 125.0,
  150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 600.0, 700.0,
  800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 2500.0, 3000.0]
const yion2: array[32, float32] = [
  0.0001, 0.006, 0.012, 0.025, 0.037, 0.05, 0.093, 0.136, 0.186, 0.341, 0.372,
  0.415, 0.415, 0.409, 0.397, 0.384, 0.378, 0.35, 0.304, 0.298, 0.267, 0.242,
  0.217, 0.198, 0.18, 0.167, 0.143, 0.124, 0.112, 0.105, 0.087, 0.074]
const xion3: array[81, float32] = [
  285.0, 298.0, 307.0, 316.0, 325.0, 335.0, 345.0, 365.0, 398.0, 422.0, 447.0,
  473.0, 501.0, 531.0, 613.0, 668.0, 708.0, 750.0, 817.0, 917.0, 1000.0,
  1122.0, 1296.0, 1496.0, 1679.0, 1884.0, 2054.0, 2238.0, 2512.0, 2985.0,
  3981.0, 5012.0, 7079.0, 10000.0, 15000.0, 20500.0, 25100.0, 30700.0,
  41000.0, 50100.0, 61300.0, 70800.0, 81800.0, 100000.0, 154000.0, 205000.0,
  299000.0, 410000.0, 501000.0, 613000.0, 708000.0, 818000.0, 1000000.0,
  1250000.0, 1500000.0, 2050000.0, 3070000.0, 4100000.0, 5010000.0, 6130000.0,
  7080000.0, 8180000.0, 10000000.0, 15000000.0, 20500000.0, 30700000.0,
  41000000.0, 50100000.0, 61300000.0, 70800000.0, 81800000.0, 100000000.0,
  150000000.0, 205000000.0, 307000000.0, 410000000.0, 501000000.0,
  613000000.0, 708000000.0, 818000000.0, 1e+9]
const yion3: array[81, float64] = [
  1e-6, 0.000166, 0.000348, 0.000525, 0.000696, 0.000863, 0.00102, 0.00133,
  0.00175, 0.00201, 0.00224, 0.00246, 0.00266, 0.00284, 0.00321, 0.00338,
  0.00347, 0.00355, 0.00365, 0.00372, 0.00375, 0.00374, 0.00368, 0.00357,
  0.00345, 0.00331, 0.00319, 0.00307, 0.00291, 0.00266, 0.00225, 0.00195,
  0.00155, 0.00121, 0.000897, 0.000707, 0.000607, 0.000521, 0.000421,
  0.000363, 0.000314, 0.000284, 0.000257, 0.000225, 0.000174, 0.00015,
  0.000128, 0.000115, 0.000109, 0.000105, 0.000103, 0.000102, 0.000101,
  0.0001005, 0.000101, 0.000103, 0.000107, 0.000111, 0.000114, 0.000117,
  0.00012, 0.000122, 0.000125, 0.000132, 0.000138, 0.000145, 0.00015,
  0.000154, 0.000158, 0.00016, 0.000163, 0.000167, 0.000174, 0.00018,
  0.000187, 0.000192, 0.000196, 0.0002, 0.000202, 0.000205, 0.000209]
const xion4: array[79, float32] = [
  685.4, 705.0, 726.0, 747.0, 770.0, 792.0, 816.0, 840.0, 865.0, 890.0, 916.0,
  944.0, 1000.0, 1090.0, 1188.0, 1296.0, 1496.0, 1679.0, 1884.0, 2054.0,
  2238.0, 2512.0, 2985.0, 3758.0, 4467.0, 5158.0, 5957.0, 7079.0, 10000.0,
  12600.0, 15000.0, 20500.0, 25100.0, 30700.0, 41000.0, 50100.0, 61300.0,
  70800.0, 81800.0, 100000.0, 150000.0, 205000.0, 299000.0, 410000.0,
  501000.0, 613000.0, 708000.0, 818000.0, 1000000.0, 1260000.0, 1500000.0,
  2050000.0, 3070000.0, 4100000.0, 5010000.0, 6130000.0, 7080000.0, 8180000.0,
  10000000.0, 12600000.0, 15000000.0, 20500000.0, 30700000.0, 41000000.0,
  50100000.0, 61300000.0, 70800000.0, 81800000.0, 100000000.0, 126000000.0,
  150000000.0, 205000000.0, 307000000.0, 410000000.0, 501000000.0,
  613000000.0, 708000000.0, 818000000.0, 1e+9]
const yion4: array[79, float64] = [
  1e-6, 0.0000339, 0.0000677, 0.0001, 0.000132, 0.000163, 0.000192, 0.000221,
  0.000248, 0.000275, 0.0003, 0.000325, 0.000371, 0.000433, 0.000487,
  0.000534, 0.000596, 0.000632, 0.000657, 0.000669, 0.000677, 0.000679,
  0.000668, 0.000633, 0.000597, 0.000562, 0.000525, 0.00048, 0.000393,
  0.000341, 0.000304, 0.000245, 0.000213, 0.000185, 0.000151, 0.000131,
  0.000114, 0.000104, 0.0000946, 0.0000832, 0.0000658, 0.000056, 0.000048,
  0.0000435, 0.0000415, 0.00004, 0.0000393, 0.0000389, 0.0000385, 0.0000386,
  0.0000389, 0.0000398, 0.0000417, 0.0000433, 0.0000445, 0.0000458, 0.0000468,
  0.0000478, 0.0000492, 0.0000509, 0.0000521, 0.0000545, 0.0000575, 0.0000596,
  0.0000612, 0.0000627, 0.0000638, 0.0000649, 0.0000664, 0.0000682, 0.0000695,
  0.0000719, 0.000075, 0.0000772, 0.0000788, 0.0000803, 0.0000814, 0.0000825,
  0.000084]
const xatt1: array[9, float32] = [
  0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25]
const yatt1: array[9, float32] = [
  5e-6, 0.000025, 0.00005, 0.00006, 0.000075, 0.00006, 0.00005, 0.000025,
  5e-6]
const xatt2: array[29, float32] = [
  0.25, 0.5, 0.75, 0.9, 1.0, 1.25, 1.5, 1.6, 1.7, 1.8, 1.9, 2.1, 2.4, 2.6,
  2.8, 3.2, 4.0, 8.0, 8.5, 9.0, 9.3, 9.5, 9.8, 10.0, 10.2, 10.6, 10.8, 11.0,
  12.0]
const yatt2: array[29, float32] = [
  2.66e-6, 0.0000113, 0.0000766, 0.000201, 0.000345, 0.000911, 0.00141,
  0.00146, 0.0014, 0.00124, 0.00103, 0.000601, 0.000187, 0.0000697, 0.0000223,
  2.94e-6, 2.41e-6, 0.0000108, 8.24e-6, 0.0000204, 0.0000307, 0.000034,
  0.0000302, 0.0000233, 0.0000157, 5.55e-6, 3.42e-6, 2.48e-6, 2e-6]
const xatt3: array[29, float32] = [
  0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.7, 1.8, 1.9, 2.1, 2.4, 2.55, 2.8, 4.0,
  5.0, 6.0, 6.6, 7.0, 7.4, 7.5, 7.6, 7.8, 8.3, 8.5, 9.0, 9.5, 10.2, 10.8,
  12.0]
const yatt3: array[29, float32] = [
  3.67e-6, 4.88e-6, 0.000013, 0.0000547, 0.000186, 0.000414, 0.000563,
  0.000587, 0.000568, 0.000427, 0.00017, 0.0000938, 0.0000433, 0.0000369,
  0.0000284, 0.0000791, 0.00008, 0.000115, 0.000194, 0.0002, 0.000196,
  0.000161, 0.0000607, 0.0000536, 0.0000762, 0.000087, 0.0000516, 0.0000183,
  3.38e-6]
const xatt4: array[26, float32] = [
  0.25, 2.0, 3.0, 4.0, 6.0, 6.7, 7.0, 7.2, 7.3, 7.5, 7.7, 7.9, 8.0, 8.2, 8.4,
  8.6, 9.0, 9.3, 9.4, 9.64, 9.9, 10.1, 10.5, 10.7, 11.0, 12.0]
const yatt4: array[26, float32] = [
  1.36e-6, 7.01e-6, 0.0000152, 0.0000374, 0.000147, 0.000204, 0.000317,
  0.000438, 0.000495, 0.000556, 0.000512, 0.000393, 0.00033, 0.000236,
  0.000205, 0.000239, 0.000485, 0.000742, 0.000808, 0.000871, 0.000762,
  0.000598, 0.000296, 0.000214, 0.000163, 0.000163]
const xatt5: array[14, float32] = [
  0.18, 0.23, 0.28, 0.34, 0.5, 0.85, 1.0, 1.25, 1.8, 2.4, 3.0, 4.5, 6.5, 10.0]
const yatt5: array[14, float32] = [
  0.00005, 0.0005, 0.01, 0.028, 0.036, 0.034, 0.018, 0.007, 0.007, 0.008,
  0.005, 0.0035, 0.0015, 0.0007]
const xvibl: array[19, float32] = [
  0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.4,
  2.6, 2.8, 3.0, 3.2]
const yvibl: array[19, float32] = [
  0.0001, 0.0044, 0.0175, 0.054, 0.1295, 0.242, 0.3011, 0.3521, 0.3867,
  0.3989, 0.3867, 0.3521, 0.3011, 0.242, 0.1295, 0.054, 0.0175, 0.0044,
  0.0001]
const xvibh: array[19, float32] = [
  0.5, 1.0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5,
  13.5, 14.5, 15.5, 16.0, 16.5]
const yvibh: array[19, float32] = [
  0.029, 0.04, 0.054, 0.092, 0.144, 0.208, 0.277, 0.339, 0.383, 0.399, 0.383,
  0.339, 0.277, 0.208, 0.144, 0.092, 0.054, 0.04, 0.029]
const xvibm: array[18, float32] = [
  0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
  1.9, 2.0, 2.1]
const yvibm: array[18, float32] = [
  0.00003, 0.0663, 0.133, 0.199, 0.265, 0.293, 0.3, 0.293, 0.255, 0.182,
  0.122, 0.0732, 0.0453, 0.0315, 0.0174, 0.0087, 0.0036, 0.00003]
const x3bp1: array[16, float32] = [
  3.8, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 15.0, 20.0, 25.0, 30.0,
  40.0, 50.0]
const y3bp1: array[16, float32] = [
  0.002, 0.1, 0.152, 0.152, 0.136, 0.126, 0.112, 0.094, 0.08, 0.066, 0.048,
  0.0266, 0.0174, 0.0122, 0.0074, 0.0048]
const x3bp2: array[14, float32] = [
  5.2, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.0,
  50.0]
const y3bp2: array[14, float32] = [
  0.0004, 0.028, 0.104, 0.138, 0.12, 0.102, 0.088, 0.074, 0.052, 0.0294,
  0.0194, 0.0134, 0.0082, 0.0052]
const xtrp1: array[12, float32] = [
  8.8, 9.8, 11.3, 12.3, 13.3, 14.3, 16.3, 18.3, 21.3, 24.3, 28.3, 31.3]
const ytrp1: array[12, float32] = [
  0.0, 0.0068, 0.0224, 0.0344, 0.0384, 0.0408, 0.0464, 0.0464, 0.0408, 0.0344,
  0.0252, 0.02]
const xtrp2: array[11, float32] = [
  10.0, 11.5, 12.5, 13.5, 14.5, 16.5, 18.5, 21.5, 24.5, 28.5, 31.5]
const ytrp2: array[11, float32] = [
  0.0, 0.04, 0.074, 0.0844, 0.0896, 0.0992, 0.1016, 0.0896, 0.0752, 0.0556,
  0.0436]
const xtrp3: array[11, float32] = [
  11.4, 12.4, 13.4, 14.4, 16.4, 18.4, 21.4, 24.4, 28.4, 31.4, 36.4]
const ytrp3: array[11, float32] = [
  0.0, 0.0552, 0.114, 0.134, 0.15, 0.158, 0.158, 0.15, 0.124, 0.0916, 0.0678]
const xtrp4: array[11, float32] = [
  14.4, 15.4, 16.4, 17.4, 19.4, 21.4, 24.4, 27.4, 31.4, 34.4, 39.4]
const ytrp4: array[11, float32] = [
  0.0, 0.23, 0.475, 0.556, 0.627, 0.66, 0.66, 0.627, 0.494, 0.383, 0.283]
const z1t: array[25, float32] = [
  11.3, 6.18, 2.8, 1.54, 0.858, 0.407, 0.251, 0.176, 0.145, 0.15, 0.167,
  0.178, 0.187, 0.193, 0.198, 0.205, 0.21, 0.218, 0.222, 0.228, 0.231, 0.233,
  0.234, 0.235, 0.235]
const z6t: array[25, float32] = [
  298.0, 178.0, 85.2, 47.5, 26.3, 12.2, 7.06, 4.45, 3.06, 2.82, 2.89, 2.99,
  3.08, 3.13, 3.18, 3.25, 3.31, 3.39, 3.44, 3.49, 3.52, 3.54, 3.55, 3.57,
  3.57]
const z9t: array[25, float32] = [
  573.0, 358.0, 179.0, 101.6, 57.3, 26.5, 15.4, 9.63, 6.52, 5.92, 6.01, 6.18,
  6.35, 6.43, 6.52, 6.65, 6.75, 6.87, 6.95, 7.02, 7.07, 7.1, 7.12, 7.13, 7.14]
const ebrm: array[25, float32] = [
  1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 50000.0, 100000.0, 200000.0,
  500000.0, 1000000.0, 2000000.0, 3000000.0, 4000000.0, 5000000.0, 6000000.0,
  8000000.0, 10000000.0, 15000000.0, 20000000.0, 30000000.0, 40000000.0,
  50000000.0, 60000000.0, 80000000.0, 100000000.0]


const c3f4h2Model* = GasModel(
  name: (" C3F4H2 2024    ISOTROPIC", " C3F4H2 2024  ANISOTROPIC"),
  elasticText: (" ELASTIC   ANISOTROPIC HFO 1234ZE", " ELASTIC   ANISOTROPIC HFO 1234ZE"),
  mass: 114.0416, e3: 11.0,
  elastic: Elastic(interp: logLog, xTotal: widen(xel), total: widen(yel),
                   xMomentum: widen(xel), momentum: widen(ymt),
                   xEps: widen(xel), eps: widen(yeps), epsConstantBelow: true,
                   low: constant, lowBelow: f64(xel, 0), lowQ: f64(yel, 0) * 1e-16,
                   lowQMomentum: f64(ymt, 0) * 1e-16),
  # single and double ionisation (shares of Born-Bethe above 3 keV), the K
  # shells (carbon three times, fluorine four times) taken out of them
  ionisation: Ionisation(
    scheme: channels, am2: 12.755, c: 146.2, renorm: fragments,
    channels: @[
      IonChannel(eion: 11.0, eoby: 11.50, nc0: 0, ec0: 0.0, openInclusive: true,
                 text: " IONISATION SUM OF SINGLE CHARGE   ELOSS=  11.0",
                 terms: @[term(Sigma(x: widen(xion1), y: widen(yion1), unit: 1e-16.cm²,
                   interp: logLog, scale: 1.2, beyond: bethe, fraction: 0.9484, tailScale: 1.2))]),
      IonChannel(eion: 38.0, eoby: 11.50, nc0: 1, ec0: 2.0,
                 text: " IONISATION DOUBLE CHARGE STATES   ELOSS=  38.0",
                 terms: @[term(Sigma(x: widen(xion2), y: widen(yion2), unit: 1e-16.cm²,
                   interp: logLog, scale: 1.2, beyond: bethe, fraction: 0.0516, tailScale: 1.2))]),
      IonChannel(eion: 285.0, eoby: 11.50, nc0: 2, ec0: 253.0, shell: true,
                 text: " IONISATION    CARBON   K-SHELL    ELOSS= 285.0", terms: @[term(table(xion3, yion3, scale = 3.0))]),
      IonChannel(eion: 685.4, eoby: 11.50, nc0: 3, ec0: 625.2, shell: true,
                 text: " IONISATION    FLUORINE K-SHELL    ELOSS= 685.4", terms: @[term(table(xion4, yion4, scale = 4.0))])]),
  # thermal attachment to dimers, three dissociative two-body and the
  # three-body attachment (table for 760 Torr and 20 °C)
  attachment: @[
    Attach(text: " ATTACHMENT THERMAL TO DIMERS", fromStart: true, factors: @[1.0'f32],
           sigma: Sigma(x: widen(xatt1), y: widen(yatt1), unit: 1e-16.cm², interp: logLog,
                        beyond: zero)),
    Attach(text: " ATTACHMENT  2 BODY DISSOCIATIVE F2H -", fromStart: true, factors: @[1.0'f32],
           sigma: Sigma(x: widen(xatt2), y: widen(yatt2), unit: 1e-16.cm², interp: logLog,
                        beyond: zero)),
    Attach(text: " ATTACHMENT  2 BODY DISSOCIATIVE C3F3H -", fromStart: true, factors: @[1.0'f32],
           sigma: Sigma(x: widen(xatt3), y: widen(yatt3), unit: 1e-16.cm², interp: logLog,
                        beyond: zero)),
    Attach(text: " ATTACHMENT  2 BODY DISSOCIATIVE F -", fromStart: true, factors: @[1.0'f32],
           sigma: Sigma(x: widen(xatt4), y: widen(yatt4), unit: 1e-16.cm², interp: logLog,
                        beyond: zero)),
    Attach(text: " ATTACHMENT 3 BODY", fromStart: true, threeBody: true, densityRef: 293.1499938964844,
           factors: @[1.0'f32, 1.0],
           sigma: Sigma(x: widen(xatt5), y: widen(yatt5), unit: 1e-16.cm², interp: logLog,
                        beyond: zero))],
  # six vibrations: Born dipole plus low and high energy resonance shapes,
  # the superelastic ones weighted by the populations, none above 1 keV;
  # moderation; B3Π, triplets, Born dipoles; bremsstrahlung on H, C, F
  # above 1 keV
  levels: @[
    Level(ein: -0.050, text: " VIB    C-F BEND                   ELOSS= -0.050",
          kind: superelastic, dip: 1.0, dipFactor: f64([0.130'f32], 0), windows: @[table(xvibl, yvibl, interp = logLog, scale = 0.4),
                    table(xvibh, yvibh, interp = logLog, scale = 2.0)],
          windowLowOpen: true,
          population: 1, closedAbove: 1000.0, peq: peqFixed, peqValue: 0.0),
    Level(ein: -0.1095, text: " VIB    C-C STRETCH                ELOSS= -0.1095",
          kind: superelastic, dip: 1.0, dipFactor: f64([0.546'f32], 0), windows: @[table(xvibl, yvibl, interp = logLog, scale = 2.0),
                    table(xvibh, yvibh, interp = logLog, scale = 2.4)],
          windowLowOpen: true,
          population: 2, closedAbove: 1000.0, peq: peqFixed, peqValue: 0.0),
    Level(ein: -0.150, text: " VIB    C-F STRETCH                ELOSS= -0.150",
          kind: superelastic, dip: 1.0, dipFactor: f64([1.209'f32], 0), windows: @[table(xvibl, yvibl, interp = logLog, scale = 1.2),
                    table(xvibh, yvibh, interp = logLog, scale = 3.0)],
          windowLowOpen: true,
          population: 3, closedAbove: 1000.0, peq: peqFixed, peqValue: 0.0),
    Level(ein: -0.2105, text: " VIB    C=C STRETCH                ELOSS= -0.2105",
          kind: superelastic, dip: 1.0, dipFactor: f64([0.585'f32], 0), windows: @[table(xvibl, yvibl, interp = logLog, scale = 1.6),
                    table(xvibh, yvibh, interp = logLog, scale = 3.2)],
          windowLowOpen: true,
          population: 4, closedAbove: 1000.0, peq: peqFixed, peqValue: 0.0),
    Level(ein: -0.310, text: " VIB    2 * C-F STRETCH            ELOSS= -0.310",
          kind: superelastic, dip: 1.0, dipFactor: f64([0.429'f32], 0), windows: @[table(xvibl, yvibl, interp = logLog, scale = 0.8),
                    table(xvibh, yvibh, interp = logLog, scale = 1.4)],
          windowLowOpen: true,
          population: 5, closedAbove: 1000.0, peq: peqFixed, peqValue: 0.0),
    Level(ein: -0.3856, text: " VIB    C-H STRETCH                ELOSS= -0.3856",
          kind: superelastic, dip: 1.0, dipFactor: f64([1.625'f32], 0), windows: @[table(xvibl, yvibl, interp = logLog, scale = 0.48),
                    table(xvibh, yvibh, interp = logLog, scale = 2.4)],
          windowLowOpen: true,
          population: 6, closedAbove: 1000.0, peq: peqFixed, peqValue: 0.0),
    Level(ein: 0.050, text: " VIB    C-F BEND                   ELOSS=  0.050",
          dip: 1.0, dipFactor: f64([0.130'f32], 0), windows: @[table(xvibl, yvibl, interp = logLog, scale = 0.4),
                    table(xvibh, yvibh, interp = logLog, scale = 2.0)],
          windowLowOpen: true,
          closedAbove: 1000.0, peq: peqElastic, gate: 0.0, peqBelow: belowZero, peqUpTo: 100.0),
    Level(ein: 0.1095, text: " VIB    C-C STRETCH                ELOSS=  0.1095",
          dip: 1.0, dipFactor: f64([0.546'f32], 0), windows: @[table(xvibl, yvibl, interp = logLog, scale = 2.0),
                    table(xvibh, yvibh, interp = logLog, scale = 2.4)],
          windowLowOpen: true,
          closedAbove: 1000.0, peq: peqElastic, gate: 0.0, peqBelow: belowZero, peqUpTo: 100.0),
    Level(ein: 0.150, text: " VIB    C-F STRETCH                ELOSS=  0.150",
          dip: 1.0, dipFactor: f64([1.209'f32], 0), windows: @[table(xvibl, yvibl, interp = logLog, scale = 1.2),
                    table(xvibh, yvibh, interp = logLog, scale = 3.0)],
          windowLowOpen: true,
          closedAbove: 1000.0, peq: peqElastic, gate: 0.0, peqBelow: belowZero, peqUpTo: 100.0),
    Level(ein: 0.2105, text: " VIB    C=C STRETCH                ELOSS=  0.2105",
          dip: 1.0, dipFactor: f64([0.585'f32], 0), windows: @[table(xvibl, yvibl, interp = logLog, scale = 1.6),
                    table(xvibh, yvibh, interp = logLog, scale = 3.2)],
          windowLowOpen: true,
          closedAbove: 1000.0, peq: peqElastic, gate: 0.0, peqBelow: belowZero, peqUpTo: 100.0),
    Level(ein: 0.310, text: " VIB    2 * C-F STRETCH            ELOSS=  0.310",
          dip: 1.0, dipFactor: f64([0.429'f32], 0), windows: @[table(xvibl, yvibl, interp = logLog, scale = 0.8),
                    table(xvibh, yvibh, interp = logLog, scale = 1.4)],
          windowLowOpen: true,
          closedAbove: 1000.0, peq: peqElastic, gate: 0.0, peqBelow: belowZero, peqUpTo: 100.0),
    Level(ein: 0.3856, text: " VIB    C-H STRETCH                ELOSS=  0.3856",
          dip: 1.0, dipFactor: f64([1.625'f32], 0), windows: @[table(xvibl, yvibl, interp = logLog, scale = 0.48),
                    table(xvibh, yvibh, interp = logLog, scale = 2.4)],
          windowLowOpen: true,
          closedAbove: 1000.0, peq: peqElastic, gate: 0.0, peqBelow: belowZero, peqUpTo: 100.0),
    Level(ein: 0.024999, text: " VIB   MODERATION        RESIDUAL ENERGY=  0.024999",
          sigma: Sigma(x: widen(xvibm), y: widen(yvibm), unit: 1e-16.cm², interp: logLog, scale: 10.0,
                       zeroBelow: true, beyond: zero), peq: peqFixed, peqValue: 0.0),
    Level(ein: 3.80, text: " 3B PI TRIPLET PART 1              ELOSS=  3.80",
          sigma: Sigma(x: widen(x3bp1), y: widen(y3bp1), unit: 1e-16.cm², interp: logLog, scale: 1.6,
                       beyond: power, n: 2.0), peq: peqElastic, gateEnergy: 50.0, peqBelow: belowZero),
    Level(ein: 5.20, text: " 3B PI TRIPLET PART 2              ELOSS=  5.20",
          sigma: Sigma(x: widen(x3bp2), y: widen(y3bp2), unit: 1e-16.cm², interp: logLog, scale: 1.6,
                       beyond: power, n: 2.0), peq: peqElastic, gateEnergy: 50.0, peqBelow: belowZero),
    Level(ein: 6.90, text: " 1B PI SINGLET PART 1  F=0.023     ELOSS=  6.90",
          sigma: Sigma(beyond: born, f: 0.0230000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 7.30, text: " 1B PI SINGLET PART 2  F=0.103     ELOSS=  7.30",
          sigma: Sigma(beyond: born, f: 0.1030000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 7.88, text: " 1B PI SINGLET PART 3  F=0.178     ELOSS=  7.88",
          sigma: Sigma(beyond: born, f: 0.1780000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 8.62, text: " 1B PI SINGLET PART 4  F=0.047     ELOSS=  8.62",
          sigma: Sigma(beyond: born, f: 0.0470000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 8.80, text: " NEUTRAL DISS   NON-DIPOLE         ELOSS=  8.8",
          sigma: Sigma(x: widen(xtrp1), y: widen(ytrp1), unit: 1e-16.cm², scale: 1.2,
                       beyond: power, n: 2.0), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 8.90, text: " NEUTRAL DISS   DIPOLE   F=0.06    ELOSS=  8.9",
          sigma: Sigma(beyond: born, f: 0.0600000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 10.0, text: " NEUTRAL DISS   NON-DIPOLE         ELOSS= 10.0",
          sigma: Sigma(x: widen(xtrp2), y: widen(ytrp2), unit: 1e-16.cm², scale: 1.2,
                       beyond: power, n: 2.0), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 10.1, text: " NEUTRAL DISS   DIPOLE   F=0.070   ELOSS= 10.1",
          sigma: Sigma(beyond: born, f: 0.0700000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 11.4, text: " NEUTRAL DISS   NON-DIPOLE         ELOSS= 11.4",
          sigma: Sigma(x: widen(xtrp3), y: widen(ytrp3), unit: 1e-16.cm², scale: 1.2,
                       beyond: power, n: 2.0), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 12.0, text: " NEUTRAL DISS   DIPOLE    F=0.13   ELOSS= 12.0",
          sigma: Sigma(beyond: born, f: 0.1300000, tailScale: 1.2), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 13.5, text: " NEUTRAL DISS   DIPOLE    F=0.48   ELOSS= 13.5",
          sigma: Sigma(beyond: born, f: 0.4800000, tailScale: 1.2), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 14.4, text: " NEUTRAL DISS   NON-DIPOLE         ELOSS= 14.4",
          sigma: Sigma(x: widen(xtrp4), y: widen(ytrp4), unit: 1e-16.cm², scale: 1.2,
                       beyond: power, n: 2.0), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 15.5, text: " NEUTRAL DISS   DIPOLE    F=5.27   ELOSS= 15.5",
          sigma: Sigma(beyond: born, f: 5.270000, tailScale: 1.2), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 20.5, text: " NEUTRAL DISS   DIPOLE    F=2.66   ELOSS= 20.5",
          sigma: Sigma(beyond: born, f: 2.660000, tailScale: 1.2), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(text: " BREMSSTRAHLUNG FROM HYDROGEN", kin: 0, openAbove: 1000.0,
          sigma: Sigma(x: widen(ebrm), y: widen(z1t), unit: 2e-24.cm², interp: logLin),
          peq: peqFixed, peqValue: 0.0),
    Level(text: " BREMSSTRAHLUNG FROM CARBON", kin: 0, openAbove: 1000.0,
          sigma: Sigma(x: widen(ebrm), y: widen(z6t), unit: 3e-24.cm², interp: logLin),
          peq: peqFixed, peqValue: 0.0),
    Level(text: " BREMSSTRAHLUNG FROM FLUORINE", kin: 0, openAbove: 1000.0,
          sigma: Sigma(x: widen(ebrm), y: widen(z9t), unit: 4e-24.cm², interp: logLin),
          peq: peqFixed, peqValue: 0.0)],
  populations: @[pop(1, group = 1), pop(2, group = 2), pop(3, group = 3), pop(4, group = 4), pop(5, group = 5), pop(6, group = 6)],
  bbBohr: 0.5291772083e-8, bbRydberg: 13.605691909790039,
  magnitudeOffsets: true,
  negative: clampNegative,
  alwaysModel2: true,
  minLevels: 6,
  appended: (3, 1000.0))

proc c3f4h2*(inp: GasInput, g: var GasXsec) = fill(c3f4h2Model, inp, g)
