## The input of Magboltz (READIN): a sequence of decks of four cards, read
## list-directed from standard input. A deck with no gases ends the input.
##
##   card 1: NGAS NMAX IPEN ITHRM EFINAL
##   card 2: the six gas numbers (1..80)
##   card 3: their six percentages, TEMPC, TORR
##   card 4: EFIELD BMAG BTHETA

import frt, units

type
  Deck* = object
    ## One run of Magboltz as read. Values that a card leaves out (`/`) keep
    ## their values from the previous deck, as in Magboltz.
    ngas*: int                      ## gases in the mixture; 0 ends the input
    nmax*: int                      ## real collisions, in units of `nscale`
    ipen*: int                      ## 0: no Penning transfers; 1: Penning transfers
    ithrm*: int                     ## 0: gas at 0 K; 1: thermal motion of the gas
    efinal*: eV                     ## upper end of the energy range; <= 0:
                                    ## found by ELIMIT
    gas*: array[1..6, int]          ## gas numbers (Magboltz's list in gases.nim)
    percent*: array[1..6, float64]  ## their shares (%)
    tempc*: float64                 ## gas temperature (°C)
    torr*: Torr                     ## gas pressure
    efield*: V•cm⁻¹                 ## electric field
    bmag*: kGauss                   ## magnetic field
    btheta*: float64                ## angle between E and B (degrees)

func penning*(d: Deck): bool = d.ipen != 0
func thermal*(d: Deck): bool = d.ithrm != 0

proc get(r: var ListReader, x: var SomeUnit) =
  ## a number in the unit of x (a value left out keeps x)
  var v = x.float
  r.get(v)
  x = typeof(x)(v)

proc readDeck*(d: var Deck): bool =
  ## READIN: reads the next deck into d; false at the end of the input
  ## (NGAS = 0, the other cards are not read then)
  var rd = listRead()
  rd.get(d.ngas)
  rd.get(d.nmax)
  rd.get(d.ipen)
  rd.get(d.ithrm)
  rd.get(d.efinal)
  if d.ngas == 0: return false
  rd = listRead()
  for k in 1 .. 6: rd.get(d.gas[k])
  rd = listRead()
  for k in 1 .. 6: rd.get(d.percent[k])
  rd.get(d.tempc)
  rd.get(d.torr)
  rd = listRead()
  rd.get(d.efield)
  rd.get(d.bmag)
  rd.get(d.btheta)
  true
