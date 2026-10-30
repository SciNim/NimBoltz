## Nim counterpart of gasdmp.f: same input, same raw output layout.
## Build: nim c -d:release --path:../../port gasdmp.nim

import std/[os, strutils]
import frt, gen/[commons, setupt, gasmix]

var
  q, peqel: Arr2[float64, 6, 4000]
  qin, peqin: Arr2[float64, 250, 4000]
  qion, peqion: Arr2[float64, 30, 4000]
  qatt: Arr2[float64, 8, 4000]
  qnull: Arr2[float64, 10, 4000]
  penfra: Arr2[float64, 3, 250]
  ev: array[1..6, float64]
  ei: array[1..250, float64]
  kin: array[1..250, int]
  kel: array[1..6, int]
  eb, eion, ec0, wk, efl, eg1, eg2: array[1..30, float64]
  nc0, ng1, ng2: array[1..30, int]
  scln: array[1..10, float64]
  scrpt: array[1..300, string]
  scrptn: array[1..10, string]
  name: string
  nin, nion, natt, nnull, last: int
  virl: float64

proc put[T](f: File, x: T) =
  when T is string:
    raiseAssert "use putStr"
  else:
    var y = x
    if f.writeBuffer(addr y, sizeof(y)) != sizeof(y): quit "write failed"

proc putStr(f: File, s: string, n: int) =
  if s.len != n: quit "string of length " & $s.len & ", expected " & $n & ": " & s
  f.write s

let f = open(paramStr(1), fmWrite)
for line in stdin.lines:
  let w = line.splitWhitespace
  if w.len == 0: continue
  let ngs = parseInt(w[0])
  let fracg = parseFloat(w[5])
  inpt.ngas = 2
  gasn.ngasn[1] = ngs
  ratio.frac[1] = fracg
  gasn.ngasn[2] = 2
  ratio.frac[2] = 100.0 - fracg
  if fracg == 100.0: inpt.ngas = 1
  for j in inpt.ngas + 1 .. 6:
    gasn.ngasn[j] = 0
    ratio.frac[j] = 0.0
  setp.nmax = 1
  inpt.ipen = 0
  thrm.ithrm = 1
  inpt.efinal = parseFloat(w[2])
  inpt.tempc = parseFloat(w[3])
  inpt.torr = parseFloat(w[4])
  setp.efield = 1000.0
  bfld.bmag = 0.0
  bfld.btheta = 0.0
  setupt(last)
  inpt.naniso = parseInt(w[1])
  # grid exactly as in MIXERT
  inpt.estep = inpt.efinal / float64(inpt.nstep)
  let ehalf = inpt.estep / 2.0
  mix2.eg[1] = ehalf
  for i in 2 .. 4000:
    let aj = float64(i - 1)
    mix2.eg[i] = ehalf + inpt.estep * aj
    mix2.eroot[i] = sqrt(mix2.eg[i])
  mix2.eroot[1] = sqrt(ehalf)
  # all outputs start out zero
  reset q; reset peqel; reset qin; reset peqin; reset qion; reset peqion
  reset qatt; reset qnull; reset penfra; reset ev; reset ei; reset kin
  reset kel; reset eb; reset eion; reset ec0; reset wk; reset efl; reset eg1
  reset eg2; reset nc0; reset ng1; reset ng2; reset scln
  for s in scrpt.mitems: s = repeat(' ', 50)
  for s in scrptn.mitems: s = repeat(' ', 50)
  name = repeat(' ', 25)
  nin = 0; nion = 0; natt = 0; nnull = 0; virl = 0.0
  gasmix(ngs, q, qin, nin, ev, ei, name, virl, eb, peqel, peqin, penfra, kel,
         kin, qion, peqion, eion, nion, qatt, natt, qnull, nnull, scln, nc0,
         ec0, wk, efl, ng1, eg1, ng2, eg2, scrpt, scrptn)
  f.put q; f.put qin; f.put nin; f.put ev; f.put ei; f.putStr(name, 25)
  f.put virl; f.put eb; f.put peqel; f.put peqin; f.put penfra; f.put kel
  f.put kin; f.put qion; f.put peqion; f.put eion; f.put nion; f.put qatt
  f.put natt; f.put qnull; f.put nnull; f.put scln; f.put nc0; f.put ec0
  f.put wk; f.put efl; f.put ng1; f.put eg1; f.put ng2; f.put eg2
  for s in scrpt: f.putStr(s, 50)
  for s in scrptn: f.putStr(s, 50)
  f.put mix2.eg; f.put mix2.eroot; f.put mix2.qt1; f.put mix2.qt2
  f.put mix2.qt3; f.put mix2.qt4; f.put dens.den
  f.put inpt.naniso; f.put inpt.estep; f.put inpt.akt
f.close()
