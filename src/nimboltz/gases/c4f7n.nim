## C4F7N (3M Novec 4710): GAS49 of Magboltz 11.21 ("C4F7N 2025").
##
## Built on the HFO-1234ze routine (GAS45): always anisotropic, elastic
## log-log with 1 - ε; single and double ionisation (log-log tables,
## Born-Bethe shares above 3 keV) less the carbon and fluorine K-shells (the
## nitrogen K-shell is computed into channel 5, beyond NION, and not
## subtracted); attachment from Rankovic 2020 below 2 eV and an F- resonance
## (the three-body channel has no data, MB-15). The same six vibrations with
## their resonances, the moderation level, five triplets (log-log, falloff
## 1/E²) and Born dipole dissociation; bremsstrahlung on N, C and F (levels
## 27..29, above 1 keV; MB-21). No total cross section (Q(1) is not set).
## Mass and resonance table lengths per MB-47 (upstream has those of GAS45).

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
  104000.0, 10400.0, 1040.0, 104.0, 51.79999923706055, 34.5,
  23.100000381469727, 20.200000762939453, 18.0, 16.0, 16.200000762939453,
  17.299999237060547, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.0, 24.0, 27.0,
  22.0, 15.5, 12.800000190734863, 11.0, 10.0, 9.0, 8.0, 7.5, 7.0, 6.25, 5.12,
  4.21, 3.74, 2.73, 2.0, 1.6, 1.3, 1.0, 0.702, 0.511, 0.4, 0.3, 0.261, 0.224,
  0.153, 0.111, 0.085, 0.0672, 0.0546, 0.0453, 0.0327, 0.0248, 0.0196, 0.0158,
  0.0131, 0.00876, 0.00631, 0.00478, 0.00391, 0.00265, 0.00193, 0.00148,
  0.00117, 0.000956, 0.000796, 0.00058, 0.000443, 0.000352, 0.000286, 0.000239,
  0.000162, 0.000118, 0.0000902, 0.0000715, 0.0000485, 0.0000354, 0.000027,
  0.0000214, 0.0000175, 0.0000146, 0.0000106, 8.11e-6, 6.43e-6, 5.24e-6,
  4.36e-6, 2.96e-6, 2.16e-6, 1.65e-6, 1.31e-6, 8.88e-7, 6.47e-7, 4.95e-7,
  3.92e-7, 3.2e-7, 2.66e-7, 1.94e-7, 1.48e-7, 1.18e-7, 9.58e-8, 7.98e-8,
  5.41e-8, 3.94e-8, 3.02e-8, 2.39e-8, 1.62e-8, 1.18e-8, 9.04e-9, 7.17e-9,
  5.84e-9, 4.87e-9, 3.54e-9, 2.71e-9, 2.15e-9, 1.75e-9, 1.46e-9, 9.9e-10,
  7.21e-10, 5.52e-10, 4.37e-10, 2.97e-10, 2.16e-10, 1.65e-10, 1.31e-10,
  1.07e-10, 8.9e-11, 6.48e-11, 4.96e-11, 3.93e-11, 3.2e-11, 2.67e-11]
const yel: array[131, float64] = [
  104000.0, 10400.0, 1040.0, 104.0, 51.79999923706055, 34.5,
  23.100000381469727, 20.200000762939453, 18.0, 16.0, 16.200000762939453,
  17.299999237060547, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.0, 24.0, 27.0,
  22.0, 15.5, 12.800000190734863, 11.0, 10.0, 9.0, 8.0, 7.5, 17.9, 16.4, 15.3,
  14.4, 13.7, 12.6, 11.7, 10.5, 9.63, 8.98, 8.04, 7.0, 5.929999828338623,
  4.960000038146973, 4.269999980926514, 3.75, 3.009999990463257,
  2.5199999809265137, 2.180000066757202, 1.9299999475479126,
  1.7300000190734863, 1.5700000524520874, 1.3200000524520874,
  1.149999976158142, 1.0099999904632568, 0.9089999794960022,
  0.8240000009536743, 0.6690000295639038, 0.5640000104904175,
  0.48899999260902405, 0.4320000112056732, 0.35600000619888306,
  0.30399999022483826, 0.26600000262260437, 0.2370000034570694,
  0.21400000154972076, 0.19499999284744263, 0.16699999570846558,
  0.1459999978542328, 0.12999999523162842, 0.11699999868869781,
  0.10700000077486038, 0.09309999644756317, 0.08290000259876251,
  0.07519999891519547, 0.06909999996423721, 0.05999999865889549,
  0.05350000038743019, 0.04919999837875366, 0.04580000042915344,
  0.0430000014603138, 0.0406000018119812, 0.03689999878406525,
  0.033900000154972076, 0.031599998474121094, 0.02969999983906746,
  0.02800000086426735, 0.026399999856948853, 0.025200000032782555,
  0.02419999986886978, 0.02329999953508377, 0.02199999988079071,
  0.020899999886751175, 0.020800000056624413, 0.02070000022649765,
  0.020600000396370888, 0.020500000566244125, 0.0203000009059906,
  0.02019999921321869, 0.019999999552965164, 0.019899999722838402,
  0.01979999989271164, 0.01979999989271164, 0.01979999989271164,
  0.01979999989271164, 0.01979999989271164, 0.01979999989271164,
  0.01979999989271164, 0.01979999989271164, 0.01979999989271164,
  0.01979999989271164, 0.01979999989271164, 0.01979999989271164,
  0.01979999989271164, 0.01979999989271164, 0.01979999989271164,
  0.01979999989271164, 0.01979999989271164, 0.01979999989271164,
  0.01979999989271164, 0.01979999989271164, 0.01979999989271164,
  0.01979999989271164, 0.01979999989271164, 0.01979999989271164,
  0.01979999989271164, 0.01979999989271164, 0.01979999989271164,
  0.01979999989271164, 0.01979999989271164, 0.01979999989271164,
  0.01979999989271164]
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
const xion: array[44, float32] = [
  12.7, 15.0, 16.0, 17.0, 18.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0, 31.0,
  33.0, 35.0, 37.0, 39.0, 44.0, 47.0, 49.0, 54.0, 59.0, 64.0, 69.0, 79.0,
  89.0, 99.0, 119.0, 149.0, 199.0, 249.0, 299.0, 399.0, 499.0, 599.0, 699.0,
  799.0, 899.0, 999.0, 1249.0, 1499.0, 1749.0, 1999.0, 2499.0, 2999.0]
const yionc: array[44, float64] = [
  0.00015, 0.0165, 0.132, 0.27, 0.465, 0.69, 1.3, 2.04, 2.67, 3.36, 3.93, 4.5,
  5.0, 5.47, 5.88, 6.3, 7.57, 8.28, 8.77, 9.55, 10.2, 10.6, 11.0, 11.6, 11.8,
  11.9, 11.9, 11.7, 11.7, 11.0, 10.4, 9.2, 8.16, 7.26, 6.6, 5.99, 5.64, 5.38,
  4.56, 3.99, 3.48, 3.25, 2.69, 2.39]
const yiong: array[44, float64] = [
  0.00015, 0.0165, 0.132, 0.27, 0.465, 0.69, 1.3, 2.04, 2.67, 3.36, 3.93, 4.5,
  5.0, 5.47, 5.91, 6.36, 7.68, 8.42, 8.96, 9.79, 10.5, 11.0, 11.4, 12.1, 12.5,
  12.6, 12.5, 12.3, 12.3, 11.7, 10.9, 9.69, 8.58, 7.63, 6.95, 6.3, 5.93, 5.66,
  4.79, 4.19, 3.66, 3.42, 2.82, 2.5]
const xion1: array[44, float32] = [
  12.7, 15.0, 16.0, 17.0, 18.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0, 31.0,
  33.0, 35.0, 37.0, 39.0, 44.0, 47.0, 49.0, 54.0, 59.0, 64.0, 69.0, 79.0,
  89.0, 99.0, 119.0, 149.0, 199.0, 249.0, 299.0, 399.0, 499.0, 599.0, 699.0,
  799.0, 899.0, 999.0, 1249.0, 1499.0, 1749.0, 1999.0, 2499.0, 2999.0]
const yion1: array[44, float64] = [
  0.00015, 0.0165, 0.132, 0.27, 0.465, 0.69, 1.3, 2.04, 2.67, 3.36, 3.93, 4.5,
  5.0, 5.47, 5.85, 6.24, 7.46, 8.14, 8.59, 9.31, 9.87, 10.2, 10.6, 11.0, 11.2,
  11.3, 11.3, 11.1, 11.1, 10.4, 9.85, 8.71, 7.73, 6.87, 6.25, 5.67, 5.35,
  5.11, 4.31, 3.76, 3.27, 3.06, 2.53, 2.27]
const xion2: array[26, float32] = [
  35.0, 37.0, 39.0, 44.0, 47.0, 49.0, 54.0, 59.0, 64.0, 69.0, 79.0, 89.0,
  99.0, 119.0, 149.0, 199.0, 249.0, 299.0, 399.0, 499.0, 599.0, 699.0, 799.0,
  899.0, 999.0, 3000.0]
const yion2: array[26, float64] = [
  0.00015, 0.03, 0.06, 0.112, 0.135, 0.18, 0.24, 0.3, 0.39, 0.42, 0.54, 0.615,
  0.615, 0.615, 0.615, 0.6, 0.57, 0.54, 0.48, 0.435, 0.39, 0.345, 0.315,
  0.291, 0.27, 0.12]
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
const xion5: array[89, float32] = [
  401.6, 407.0, 419.0, 431.0, 444.0, 457.0, 471.0, 485.0, 499.0, 515.0, 530.0,
  546.0, 562.0, 596.0, 631.0, 668.0, 708.0, 750.0, 794.0, 841.0, 891.0, 944.0,
  1000.0, 1090.0, 1188.0, 1296.0, 1496.0, 1679.0, 1884.0, 2054.0, 2238.0,
  2512.0, 2985.0, 3758.0, 4467.0, 5158.0, 5957.0, 7079.0, 10000.0, 12600.0,
  15000.0, 20500.0, 25100.0, 30700.0, 41000.0, 50100.0, 61300.0, 70800.0,
  81800.0, 100000.0, 154000.0, 205000.0, 299000.0, 410000.0, 501000.0,
  613000.0, 708000.0, 818000.0, 1000000.0, 1250000.0, 1500000.0, 2050000.0,
  3070000.0, 4100000.0, 5010000.0, 6130000.0, 7080000.0, 8180000.0,
  10000000.0, 12200000.0, 15000000.0, 18800000.0, 22400000.0, 28200000.0,
  37600000.0, 48700000.0, 63100000.0, 86600000.0, 100000000.0, 122000000.0,
  150000000.0, 188000000.0, 224000000.0, 282000000.0, 376000000.0,
  487000000.0, 631000000.0, 866000000.0, 1e+9]
const yion5: array[89, float64] = [
  0.0, 0.0000411, 0.000138, 0.000232, 0.000323, 0.000412, 0.000498, 0.000581,
  0.000661, 0.000738, 0.000813, 0.000885, 0.000954, 0.00108, 0.0012, 0.00131,
  0.00141, 0.0015, 0.00158, 0.00165, 0.00172, 0.00177, 0.00182, 0.00187,
  0.00191, 0.00194, 0.00195, 0.00193, 0.00189, 0.00186, 0.00181, 0.00174,
  0.00163, 0.00145, 0.00132, 0.00122, 0.00112, 0.001, 0.000795, 0.000677,
  0.000598, 0.000474, 0.000409, 0.000352, 0.000285, 0.000247, 0.000214,
  0.000194, 0.000176, 0.000154, 0.000119, 0.000103, 0.0000877, 0.0000792,
  0.0000754, 0.0000727, 0.0000714, 0.0000704, 0.0000697, 0.0000696, 0.00007,
  0.0000716, 0.0000747, 0.0000774, 0.0000795, 0.0000817, 0.0000834, 0.0000851,
  0.0000875, 0.00009, 0.0000925, 0.0000954, 0.0000976, 0.000101, 0.000104,
  0.000108, 0.000111, 0.000115, 0.000117, 0.00012, 0.000122, 0.000125,
  0.000127, 0.000131, 0.000134, 0.000138, 0.000141, 0.000145, 0.000147]
const xatt1: array[9, float64] = [
  5.0, 5.58, 6.22, 6.94, 7.75, 8.64, 9.64, 10.8, 12.0]
const yatt1: array[9, float64] = [
  0.00381, 0.00456, 0.00865, 0.0143, 0.0161, 0.00908, 0.00455, 0.00271,
  0.00287]
const xatt2: array[29, float32] = [
  0.001, 0.0125, 0.0156, 0.0195, 0.0243, 0.0303, 0.0379, 0.0473, 0.059,
  0.0737, 0.092, 0.115, 0.143, 0.179, 0.224, 0.279, 0.349, 0.435, 0.543,
  0.678, 0.847, 1.06, 1.32, 1.65, 2.06, 2.57, 3.21, 4.0, 5.0]
const yatt2: array[29, float32] = [
  25.5, 24.0, 21.5, 14.1, 11.2, 10.3, 9.15, 7.74, 5.94, 3.96, 3.17, 2.17,
  1.55, 0.999, 0.649, 0.431, 0.231, 0.115, 0.0735, 0.105, 0.105, 0.217, 0.381,
  0.438, 0.253, 0.0884, 0.0884, 0.0224, 0.0224]
const xvibl: array[21, float32] = [
  0.025, 0.175, 0.325, 0.475, 0.625, 0.775, 0.925, 1.075, 1.225, 1.375, 1.525,
  1.675, 1.825, 1.975, 2.125, 2.275, 2.425, 2.575, 2.725, 2.875, 3.025]
const yvibl: array[21, float32] = [
  0.0601, 0.1032, 0.1619, 0.232, 0.3039, 0.3639, 0.3981, 0.3981, 0.3639,
  0.3039, 0.232, 0.1619, 0.1032, 0.0601, 0.032, 0.0156, 0.0069, 0.0028, 0.001,
  0.0004, 0.0001]
const xvibh: array[20, float32] = [
  2.5, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
  10.5, 11.0, 11.5, 12.0, 12.7]
const yvibh: array[20, float32] = [
  0.001, 0.0273, 0.0531, 0.0948, 0.1553, 0.2331, 0.3211, 0.4056, 0.4698,
  0.4991, 0.4863, 0.4346, 0.3562, 0.2677, 0.1846, 0.1167, 0.0677, 0.036,
  0.0176, 0.001]
const xvibm: array[18, float32] = [
  0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
  1.9, 2.0, 2.1]
const yvibm: array[18, float32] = [
  0.00003, 0.0663, 0.133, 0.199, 0.265, 0.293, 0.3, 0.293, 0.255, 0.182,
  0.122, 0.0732, 0.0453, 0.0315, 0.0174, 0.0087, 0.0036, 0.00003]
const xtrp1: array[12, float32] = [
  6.0, 6.68, 7.71, 8.39, 9.07, 9.75, 11.1, 12.5, 14.5, 16.5, 19.3, 21.3]
const ytrp1: array[12, float32] = [
  0.001, 0.0068, 0.0224, 0.0344, 0.0384, 0.0408, 0.0464, 0.0464, 0.0408,
  0.0344, 0.0252, 0.02]
const xtrp2: array[12, float32] = [
  7.5, 8.35, 9.64, 10.5, 11.3, 12.2, 13.9, 15.6, 18.1, 20.6, 24.1, 26.6]
const ytrp2: array[12, float32] = [
  0.0012, 0.00816, 0.0269, 0.0413, 0.0461, 0.049, 0.0557, 0.0557, 0.049,
  0.0413, 0.0302, 0.024]
const xtrp3: array[12, float32] = [
  8.5, 9.46, 10.9, 11.9, 12.8, 13.8, 15.7, 17.7, 20.5, 23.4, 27.3, 30.2]
const ytrp3: array[12, float32] = [
  0.001, 0.0068, 0.0224, 0.0344, 0.0384, 0.0408, 0.0464, 0.0464, 0.0408,
  0.0344, 0.0252, 0.02]
const xtrp4: array[12, float32] = [
  10.2, 11.4, 13.1, 14.3, 15.4, 16.6, 18.9, 21.2, 24.6, 28.1, 32.8, 36.2]
const ytrp4: array[12, float32] = [
  0.001, 0.0068, 0.0224, 0.0344, 0.0384, 0.0408, 0.0464, 0.0464, 0.0408,
  0.0344, 0.0252, 0.02]
const xtrp5: array[12, float32] = [
  12.5, 13.9, 16.1, 17.5, 18.9, 20.3, 23.1, 26.0, 30.2, 34.4, 40.2, 44.4]
const ytrp5: array[12, float32] = [
  0.001, 0.0068, 0.0224, 0.0344, 0.0384, 0.0408, 0.0464, 0.0464, 0.0408,
  0.0344, 0.0252, 0.02]
const z6t: array[25, float32] = [
  298.0, 178.0, 85.2, 47.5, 26.3, 12.2, 7.06, 4.45, 3.06, 2.82, 2.89, 2.99,
  3.08, 3.13, 3.18, 3.25, 3.31, 3.39, 3.44, 3.49, 3.52, 3.54, 3.55, 3.57,
  3.57]
const z7t: array[25, float32] = [
  385.0, 234.0, 113.0, 63.5, 35.5, 16.4, 9.48, 5.96, 4.07, 3.73, 3.81, 3.93,
  4.04, 4.11, 4.18, 4.26, 4.33, 4.42, 4.48, 4.55, 4.59, 4.61, 4.63, 4.64,
  4.65]
const z9t: array[25, float32] = [
  573.0, 358.0, 179.0, 101.6, 57.3, 26.5, 15.4, 9.63, 6.52, 5.92, 6.01, 6.18,
  6.35, 6.43, 6.52, 6.65, 6.75, 6.87, 6.95, 7.02, 7.07, 7.1, 7.12, 7.13, 7.14]
const ebrm: array[25, float32] = [
  1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 50000.0, 100000.0, 200000.0,
  500000.0, 1000000.0, 2000000.0, 3000000.0, 4000000.0, 5000000.0, 6000000.0,
  8000000.0, 10000000.0, 15000000.0, 20000000.0, 30000000.0, 40000000.0,
  50000000.0, 60000000.0, 80000000.0, 100000000.0]


const c4f7nModel* = GasModel(
  name: (" C4F7N 2025    ISOTROPIC", " C4F7N 2025  ANISOTROPIC"),
  elasticText: (" ELASTIC   ANISOTROPIC C4F7N", " ELASTIC   ANISOTROPIC C4F7N"),
  mass: 195.0383, e3: 11.0,
  elastic: Elastic(interp: logLog, xTotal: widen(xel), total: widen(yel),
                   xMomentum: widen(xel), momentum: widen(ymt),
                   xEps: widen(xel), eps: widen(yeps), epsConstantBelow: true,
                   low: constant, lowBelow: f64(xel, 0), lowQ: f64(yel, 0) * 1e-16,
                   lowQMomentum: f64(ymt, 0) * 1e-16),
  # single and double ionisation (shares of Born-Bethe above 3 keV), the K
  # shells (carbon three times, fluorine four times) taken out of them
  ionisation: Ionisation(
    scheme: channels, am2: 19.0, c: 232.0, renorm: fragments,
    channels: @[
      IonChannel(eion: 12.7, eoby: 11.50, nc0: 0, ec0: 0.0, openInclusive: true,
                 text: " IONISATION SUM OF SINGLE CHARGE   ELOSS=  12.7",
                 terms: @[term(Sigma(x: widen(xion1), y: widen(yion1), unit: 1e-16.cm²,
                   interp: logLog, beyond: bethe, fraction: 0.9484))]),
      IonChannel(eion: 38.0, eoby: 11.50, nc0: 1, ec0: 2.0,
                 text: " IONISATION DOUBLE CHARGE STATES   ELOSS=  35.0",
                 terms: @[term(Sigma(x: widen(xion2), y: widen(yion2), unit: 1e-16.cm²,
                   interp: logLog, beyond: bethe, fraction: 0.0516))]),
      IonChannel(eion: 285.0, eoby: 11.50, nc0: 2, ec0: 253.0, shell: true,
                 text: " IONISATION    CARBON   K-SHELL    ELOSS= 285.0", terms: @[term(table(xion3, yion3, scale = 3.0))]),
      IonChannel(eion: 685.4, eoby: 11.50, nc0: 3, ec0: 625.2, shell: true,
                 text: " IONISATION    FLUORINE K-SHELL    ELOSS= 685.4", terms: @[term(table(xion4, yion4, scale = 4.0))])]),
  # total attachment up to 2 eV (Rankovic 2020) and the F- resonance at
  # 10 eV; the other three channels are empty (the first text is that of
  # the nitrogen K shell, which lies beyond NION in Magboltz)
  attachment: @[
    Attach(text: " IONISATION    NITROGEN K-SHELL    ELOSS= 401.6", fromStart: true, factors: @[39.0'f32],
           sigma: Sigma(x: widen(xatt1), y: widen(yatt1), unit: 1e-16.cm², interp: logLog,
                        beyond: zero)),
    Attach(text: " ATTACHMENT TOTAL FROM RANKOVIC 2020 UP TO 2 EV", fromStart: true, factors: @[5.3'f32],
           sigma: Sigma(x: widen(xatt2), y: widen(yatt2), unit: 1e-16.cm², interp: logLog,
                        beyond: zero)),
    attach(" ATTACHMENT F- RESONANCE AT 10 EV"),
    attach(""),
    attach("")],
  # six vibrations: Born dipole plus low and high energy resonance shapes,
  # the superelastic ones weighted by the populations, none above 1 keV;
  # moderation; triplets, Born dipoles; bremsstrahlung on N, C, F
  # above 1 keV
  levels: @[
    Level(ein: -0.065, text: " VIB    CCN BEND                   ELOSS= -0.065",
          kind: superelastic, dip: 1.0, dipFactor: f64([0.2'f32], 0) * 2.0, windows: @[Sigma(x: widen(xvibl), y: widen(yvibl), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([0.4'f32], 0) * f64([2.0'f32], 0)),
                    Sigma(x: widen(xvibh), y: widen(yvibh), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([2.0'f32], 0) * f64([1.1'f32], 0))],
          windowLowOpen: true,
          population: 1, closedAbove: 1000.0, peq: peqFixed, peqValue: 0.0),
    Level(ein: -0.088, text: " VIB    CF3 ASYM UMBRELLA          ELOSS= -0.088",
          kind: superelastic, dip: 1.0, dipFactor: f64([0.2'f32], 0) * 2.0, windows: @[Sigma(x: widen(xvibl), y: widen(yvibl), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([2.0'f32], 0) * f64([2.0'f32], 0)),
                    Sigma(x: widen(xvibh), y: widen(yvibh), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([2.4'f32], 0) * f64([1.1'f32], 0))],
          windowLowOpen: true,
          population: 2, closedAbove: 1000.0, peq: peqFixed, peqValue: 0.0),
    Level(ein: -0.122, text: " VIB    CF STRETCH                 ELOSS= -0.122",
          kind: superelastic, dip: 1.0, dipFactor: f64([0.3'f32], 0) * 2.0, windows: @[Sigma(x: widen(xvibl), y: widen(yvibl), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([1.2'f32], 0) * f64([2.0'f32], 0)),
                    Sigma(x: widen(xvibh), y: widen(yvibh), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([3.0'f32], 0) * f64([1.1'f32], 0))],
          windowLowOpen: true,
          population: 3, closedAbove: 1000.0, peq: peqFixed, peqValue: 0.0),
    Level(ein: -0.155, text: " VIB    CC STRETCH                 ELOSS= -0.155",
          kind: superelastic, dip: 1.0, dipFactor: f64([1.0'f32], 0) * 2.0, windows: @[Sigma(x: widen(xvibl), y: widen(yvibl), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([1.6'f32], 0) * f64([2.0'f32], 0)),
                    Sigma(x: widen(xvibh), y: widen(yvibh), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([3.2'f32], 0) * f64([1.1'f32], 0))],
          windowLowOpen: true,
          population: 4, closedAbove: 1000.0, peq: peqFixed, peqValue: 0.0),
    Level(ein: -0.279, text: " VIB    CN  and 2 CC STRETCH       ELOSS= -0.279",
          kind: superelastic, dip: 1.0, dipFactor: f64([0.12'f32], 0) * 2.0, windows: @[Sigma(x: widen(xvibl), y: widen(yvibl), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([0.8'f32], 0) * f64([2.0'f32], 0)),
                    Sigma(x: widen(xvibh), y: widen(yvibh), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([1.4'f32], 0) * f64([1.1'f32], 0))],
          windowLowOpen: true,
          population: 5, closedAbove: 1000.0, peq: peqFixed, peqValue: 0.0),
    Level(ein: -0.443, text: " VIB    CN + CC AND 3CC STRETCH    ELOSS= -0.443",
          kind: superelastic, dip: 1.0, dipFactor: f64([0.1'f32], 0) * 2.0, windows: @[Sigma(x: widen(xvibl), y: widen(yvibl), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([0.48'f32], 0) * f64([2.0'f32], 0)),
                    Sigma(x: widen(xvibh), y: widen(yvibh), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([2.4'f32], 0) * f64([1.1'f32], 0))],
          windowLowOpen: true,
          population: 6, closedAbove: 1000.0, peq: peqFixed, peqValue: 0.0),
    Level(ein: 0.065, text: " VIB    CCN BEND                   ELOSS=  0.065",
          dip: 1.0, dipFactor: f64([0.2'f32], 0) * 2.0, windows: @[Sigma(x: widen(xvibl), y: widen(yvibl), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([0.4'f32], 0) * f64([2.0'f32], 0)),
                    Sigma(x: widen(xvibh), y: widen(yvibh), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([2.0'f32], 0) * f64([1.1'f32], 0))],
          windowLowOpen: true,
          closedAbove: 1000.0, peq: peqElastic, gate: 0.0, peqBelow: belowZero, peqUpTo: 100.0),
    Level(ein: 0.088, text: " VIB    CF3 ASYM UMBRELLA          ELOSS=  0.088",
          dip: 1.0, dipFactor: f64([0.2'f32], 0) * 2.0, windows: @[Sigma(x: widen(xvibl), y: widen(yvibl), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([2.0'f32], 0) * f64([2.0'f32], 0)),
                    Sigma(x: widen(xvibh), y: widen(yvibh), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([2.4'f32], 0) * f64([1.1'f32], 0))],
          windowLowOpen: true,
          closedAbove: 1000.0, peq: peqElastic, gate: 0.0, peqBelow: belowZero, peqUpTo: 100.0),
    Level(ein: 0.122, text: " VIB    CF STRETCH                 ELOSS=  0.122",
          dip: 1.0, dipFactor: f64([0.3'f32], 0) * 2.0, windows: @[Sigma(x: widen(xvibl), y: widen(yvibl), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([1.2'f32], 0) * f64([2.0'f32], 0)),
                    Sigma(x: widen(xvibh), y: widen(yvibh), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([3.0'f32], 0) * f64([1.1'f32], 0))],
          windowLowOpen: true,
          closedAbove: 1000.0, peq: peqElastic, gate: 0.0, peqBelow: belowZero, peqUpTo: 100.0),
    Level(ein: 0.155, text: " VIB    CC STRETCH                 ELOSS=  0.155",
          dip: 1.0, dipFactor: f64([1.0'f32], 0) * 2.0, windows: @[Sigma(x: widen(xvibl), y: widen(yvibl), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([1.6'f32], 0) * f64([2.0'f32], 0)),
                    Sigma(x: widen(xvibh), y: widen(yvibh), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([3.2'f32], 0) * f64([1.1'f32], 0))],
          windowLowOpen: true,
          closedAbove: 1000.0, peq: peqElastic, gate: 0.0, peqBelow: belowZero, peqUpTo: 100.0),
    Level(ein: 0.279, text: " VIB    CN and 2 CC STRETCH        ELOSS=  0.279",
          dip: 1.0, dipFactor: f64([0.12'f32], 0) * 2.0, windows: @[Sigma(x: widen(xvibl), y: widen(yvibl), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([0.8'f32], 0) * f64([2.0'f32], 0)),
                    Sigma(x: widen(xvibh), y: widen(yvibh), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([1.4'f32], 0) * f64([1.1'f32], 0))],
          windowLowOpen: true,
          closedAbove: 1000.0, peq: peqElastic, gate: 0.0, peqBelow: belowZero, peqUpTo: 100.0),
    Level(ein: 0.443, text: " VIB    CN + CC and 3 CC STRETCH   ELOSS=  0.443",
          dip: 1.0, dipFactor: f64([0.1'f32], 0) * 2.0, windows: @[Sigma(x: widen(xvibl), y: widen(yvibl), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([0.48'f32], 0) * f64([2.0'f32], 0)),
                    Sigma(x: widen(xvibh), y: widen(yvibh), unit: 1e-16.cm², interp: logLog,
                          scaleD: f64([2.4'f32], 0) * f64([1.1'f32], 0))],
          windowLowOpen: true,
          closedAbove: 1000.0, peq: peqElastic, gate: 0.0, peqBelow: belowZero, peqUpTo: 100.0),
    Level(ein: 0.024999, text: "",
          sigma: Sigma(x: widen(xvibm), y: widen(yvibm), unit: 1e-16.cm², interp: logLog, scale: 10.0,
                       zeroBelow: true, beyond: zero), peq: peqFixed, peqValue: 0.0),
    Level(ein: 6.00, text: " NEUTRAL DISS   NON-DIPOLE         ELOSS=  6.0",
          sigma: Sigma(x: widen(xtrp1), y: widen(ytrp1), unit: 1e-16.cm², interp: logLog,
                       scale: 5.2, beyond: power, n: 2.0), peq: peqElastic, gateEnergy: 50.0, peqBelow: belowZero),
    Level(ein: 7.50, text: " NEUTRAL DISS   NON-DIPOLE         ELOSS=  7.5",
          sigma: Sigma(x: widen(xtrp2), y: widen(ytrp2), unit: 1e-16.cm², interp: logLog,
                       scale: 5.2, beyond: power, n: 2.0), peq: peqElastic, gateEnergy: 50.0, peqBelow: belowZero),
    Level(ein: 8.50, text: " NEUTRAL DISS   NON-DIPOLE         ELOSS=  8.5",
          sigma: Sigma(x: widen(xtrp3), y: widen(ytrp3), unit: 1e-16.cm², interp: logLog,
                       scale: 5.2, beyond: power, n: 2.0), peq: peqElastic, gateEnergy: 50.0, peqBelow: belowZero),
    Level(ein: 8.55, text: " NEUTRAL DISS   DIPOLE   F=0.02    ELOSS=  8.55",
          sigma: Sigma(beyond: born, f: 0.0200000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 9.10, text: " NEUTRAL DISS   DIPOLE   F=0.02    ELOSS=  9.1",
          sigma: Sigma(beyond: born, f: 0.0200000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 10.2, text: " NEUTRAL DISS   NON-DIPOLE         ELOSS= 10.2",
          sigma: Sigma(x: widen(xtrp4), y: widen(ytrp4), unit: 1e-16.cm², interp: logLog,
                       scale: 5.2, beyond: power, n: 2.0), peq: peqElastic, gateEnergy: 50.0, peqBelow: belowZero),
    Level(ein: 10.5, text: " NEUTRAL DISS   DIPOLE   F=0.48    ELOSS= 10.5",
          sigma: Sigma(beyond: born, f: 0.4800000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 12.0, text: " NEUTRAL DISS   DIPOLE   F=0.2     ELOSS= 12.0",
          sigma: Sigma(beyond: born, f: 0.2000000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 12.5, text: " NEUTRAL DISS   NON-DIPOLE         ELOSS= 12.5",
          sigma: Sigma(x: widen(xtrp5), y: widen(ytrp5), unit: 1e-16.cm², interp: logLog,
                       scale: 5.2, beyond: power, n: 2.0), peq: peqElastic, gateEnergy: 50.0, peqBelow: belowZero),
    Level(ein: 13.5, text: " NEUTRAL DISS   DIPOLE   F=0.3     ELOSS= 13.5",
          sigma: Sigma(beyond: born, f: 0.3000000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 14.5, text: " NEUTRAL DISS   DIPOLE   F=0.2     ELOSS= 14.5",
          sigma: Sigma(beyond: born, f: 0.2000000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 16.5, text: " NEUTRAL DISS   DIPOLE   F=2.0     ELOSS= 16.5",
          sigma: Sigma(beyond: born, f: 7.0000000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(ein: 20.5, text: " NEUTRAL DISS   DIPOLE   F=10.0    ELOSS= 20.5",
          sigma: Sigma(beyond: born, f: 14.000000), peq: peqElastic, gate: 3.0, peqBelow: belowZero),
    Level(text: " BREMSSTRAHLUNG FROM NITROGEN", kin: 0, openAbove: 1000.0,
          sigma: Sigma(x: widen(ebrm), y: widen(z7t), unit: 1e-24.cm², interp: logLin),
          peq: peqFixed, peqValue: 0.0),
    Level(text: " BREMSSTRAHLUNG FROM CARBON", kin: 0, openAbove: 1000.0,
          sigma: Sigma(x: widen(ebrm), y: widen(z6t), unit: 4e-24.cm², interp: logLin),
          peq: peqFixed, peqValue: 0.0),
    Level(text: " BREMSSTRAHLUNG FROM FLUORINE", kin: 0, openAbove: 1000.0,
          sigma: Sigma(x: widen(ebrm), y: widen(z9t), unit: 7e-24.cm², interp: logLin),
          peq: peqFixed, peqValue: 0.0)],
  populations: @[pop(1, group = 1), pop(2, group = 2), pop(3, group = 3), pop(4, group = 4), pop(5, group = 5), pop(6, group = 6)],
  bbBohr: 0.5291772083e-8, bbRydberg: 13.605691909790039,
  magnitudeOffsets: true,
  negative: clampNegative,
  alwaysModel2: true,
  minLevels: 6,
  appended: (3, 1000.0))

proc c4f7n*(inp: GasInput, g: var GasXsec) = fill(c4f7nModel, inp, g)
