## Generates the input decks of the bitwise regression suite (`decks/*.in`).
## Run from this directory: `nim r gendecks.nim`.
##
## Deck format (free format, see READIN in magboltz-11.21.f):
##   NGAS NMAX IPEN ITHRM EFINAL / NGASN(1..6) / FRAC(1..6) TEMPC TORR /
##   EFIELD BMAG BTHETA, repeated; a card with NGAS=0 terminates.

import std/[strformat, strutils, os]

type Deck = object
  gases: seq[(int, float)]   # (gas number, percentage)
  nmax = 1
  ipen = 0
  ithrm = 0
  efinal = 0.0
  tempc = 20.0
  torr = 760.0
  efield = 1000.0
  bmag = 0.0
  btheta = 0.0

proc card(d: Deck): string =
  var ids, fracs: seq[string]
  for i in 0 ..< 6:
    if i < d.gases.len:
      ids.add $d.gases[i][0]
      fracs.add &"{d.gases[i][1]:.4f}"
    else:
      ids.add "0"
      fracs.add "0.0"
  result = &"{d.gases.len} {d.nmax} {d.ipen} {d.ithrm} {d.efinal:.5f}\n"
  result.add ids.join(" ") & "\n"
  result.add fracs.join(" ") & &" {d.tempc:.4f} {d.torr:.4f}\n"
  result.add &"{d.efield:.3f} {d.bmag:.3f} {d.btheta:.3f}\n"

proc write(name: string, decks: varargs[Deck]) =
  var s = ""
  for d in decks: s.add d.card
  s.add "0 0 0 0 0.0\n"
  writeFile("decks" / name & ".in", s)

proc bname(ang: float): string =
  ## deck name part for the angle between E and B
  if ang == 0.0: "par" elif ang == 90.0: "perp" else: $ang.int

const nGases = 61   # GAS62..80 are dummy routines

createDir "decks"
# 1. Every gas once, as 10% admixture to argon (pure argon for GAS2);
#    alternating static / thermal gas so both MIXER and MIXERT see every gas.
for g in 1 .. nGases:
  let mix = if g == 2: @[(2, 100.0)] else: @[(2, 90.0), (g, 10.0)]
  write(&"gas{g:02}", Deck(gases: mix, ithrm: g mod 2))

# 2. Field configurations: B = 0, B parallel E, B perpendicular E, B at 45°,
#    each for static and thermal gas (MONTE, MONTEA, MONTEB, MONTEC + T versions).
let arco2 = @[(2, 90.0), (12, 10.0)]
for th in 0 .. 1:
  write(&"arco2_nob_t{th}", Deck(gases: arco2, ithrm: th))
  for ang in [0.0, 90.0, 45.0]:
    write(&"arco2_b{bname(ang)}_t{th}",
          Deck(gases: arco2, ithrm: th, bmag: 10.0, btheta: ang))

# 3. High field: Townsend coefficients above the SST threshold -> ALPCALC/ALPCLCx
#    and the MONTEFT/MONTEFD time-of-flight / steady-state routines.
for th in 0 .. 1:
  write(&"hi_arco2_nob_t{th}", Deck(gases: arco2, ithrm: th, efield: 30000.0))
  for ang in [0.0, 90.0, 45.0]:
    write(&"hi_arco2_b{bname(ang)}_t{th}",
          Deck(gases: arco2, ithrm: th, efield: 30000.0, bmag: 50.0, btheta: ang))
  # attachment dominated
  write(&"hi_sf6_t{th}", Deck(gases: @[(30, 100.0)], ithrm: th, efield: 50000.0))
  write(&"hi_air_t{th}", Deck(gases: @[(16, 79.0), (15, 21.0)], ithrm: th,
                              efield: 30000.0))

# 4. Penning transfer, explicit EFINAL, six-gas mixture, pure molecular gases.
for th in 0 .. 1:
  write(&"penning_aribut_t{th}",
        Deck(gases: @[(2, 97.7), (11, 2.3)], ipen: 1, ithrm: th, efield: 500.0))
  write(&"penning_neco2_t{th}",
        Deck(gases: @[(5, 90.0), (12, 10.0)], ipen: 1, ithrm: th, efield: 20000.0))
  write(&"efinal_ar_t{th}", Deck(gases: @[(2, 100.0)], ithrm: th, efinal: 12.0))
  write(&"mix6_t{th}", Deck(gases: @[(2, 50.0), (12, 20.0), (8, 10.0), (1, 10.0),
                                     (5, 5.0), (16, 5.0)], ithrm: th, efield: 2000.0))
  write(&"pure_ch4_t{th}", Deck(gases: @[(8, 100.0)], ithrm: th, efield: 300.0,
                                torr: 1000.0, tempc: -10.0))
  write(&"pure_co2_t{th}", Deck(gases: @[(12, 100.0)], ithrm: th, efield: 3000.0))
  write(&"pure_h2o_t{th}", Deck(gases: @[(14, 100.0)], ithrm: th, efield: 1000.0))

# 5. Several decks in one run: state carried over between decks (RNG stream,
#    COMMON blocks) must match too.
write("chain", Deck(gases: @[(2, 100.0)]),
      Deck(gases: @[(2, 80.0), (12, 20.0)], ithrm: 1, bmag: 5.0, btheta: 30.0),
      Deck(gases: @[(1, 100.0)], efield: 20000.0))
