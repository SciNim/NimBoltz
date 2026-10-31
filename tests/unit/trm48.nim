## Compares the Nim RNG bitwise against the Fortran one (harness/rng/rngdump
## in the Magboltz reference repo, built from magboltz-11.21.f).
import std/[os, osproc, strutils]
import nimboltz/rm48

const N = 10_000_000

let mbDir = getEnv("MAGBOLTZ_DIR", getHomeDir() / "src/Magboltz")
let dump = getTempDir() / "rng_fortran.bin"
let (msg, code) = execCmdEx(quoteShellCommand([mbDir / "build/rngdump", $N, dump]))
doAssert code == 0, msg
doAssert msg == " RM48 INITIALIZED:  54217137           0         0\n", msg

var f = open(dump)
var buf = newSeq[float64](N)
doAssert f.readBuffer(buf[0].addr, N*sizeof(float64)) == N*sizeof(float64)
f.close()

var g = initMagRng(nil)
for i in 0 ..< N:
  let x = g.drand48()
  if cast[uint64](x) != cast[uint64](buf[i]):
    quit "mismatch at draw " & $(i+1) & ": nim " & cast[uint64](x).toHex & " fortran " &
      cast[uint64](buf[i]).toHex
echo "first ", N, " drand48 values match bitwise"

# init message goes to the given File
let tmp = getTempDir() / "rm48_msg.txt"
var o = open(tmp, fmWrite)
var g2 = initMagRng(o)
discard g2.drand48()
o.close()
doAssert readFile(tmp) == msg
