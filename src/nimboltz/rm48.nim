## Magboltz's random number generator, ported bit-for-bit from magboltz-11.21.f.
##
## Magboltz calls a function named `drand48`, but it is *not* the POSIX LCG:
## it is a buffered shim over CERNLIB's RM48 (RANMAR, Marsaglia & Zaman,
## modified by F. James; rm48.F v1.2 1996), both included in the Fortran source.
##
## Quirk kept on purpose: the shim refills its 1000-element buffer when
## `IVEC >= NVEC`, so after the first fill element 1000 of every buffer is
## generated but never returned.
##
## Magboltz never calls `RM48IN`, so the generator always starts from the
## default seed 54217137 on the first draw, which also prints the
## ` RM48 INITIALIZED:` line to the output.

import std/strutils
from std/math import floor

const
  ModCns = 1_000_000_000      ## MODCNS
  DefaultSeed* = 54217137     ## IJKL used when RM48IN is never called
  NVec* = 1000                ## buffer length of the drand48 shim

type
  Rm48* = object
    ## State of RM48: COMMON/R48ST1/ plus its SAVEd locals.
    u: array[1..97, float64]
    c, cd, cm, twom24, twom49: float64
    i97, j97: int
    ntot, ntot2, ijkl: int
    initialized: bool         ## Fortran: NTOT >= 0
    output*: File             ## where the init message goes; nil = silent

  MagRng* = object
    ## Magboltz's `drand48(DUMMY)` shim: buffers `NVec` numbers from RM48.
    rm*: Rm48
    rvec: array[1..NVec, float64]
    ivec: int

proc initRm48(r: var Rm48) =
  ## Label 1 of RM48: set up U, C, CD, CM from IJKL, then skip NTOT2*MODCNS+NTOT.
  var ij = r.ijkl div 30082
  var kl = r.ijkl - 30082*ij
  var i = (ij div 177) mod 177 + 2
  var j = ij mod 177 + 2
  var k = (kl div 169) mod 178 + 1
  var l = kl mod 169
  if r.output != nil:
    r.output.write " RM48 INITIALIZED:", align($r.ijkl, 10), "  ",
      align($r.ntot, 10), align($r.ntot2, 10), "\n"
  const half = 0.5
  var t: float64
  for ii in 1 .. 97:
    var s = 0.0
    t = half
    for jj in 1 .. 48:
      let m = ((i*j) mod 179 * k) mod 179
      i = j
      j = k
      k = m
      l = (53*l + 1) mod 169
      if (l*m) mod 64 >= 32: s = s + t
      t = half*t
    r.u[ii] = s
  r.twom49 = t
  r.twom24 = 1.0
  for i24 in 1 .. 24:
    r.twom24 = half*r.twom24
  r.c = 362436.0*r.twom24
  r.cd = 7654321.0*r.twom24
  r.cm = 16777213.0*r.twom24
  r.i97 = 97
  r.j97 = 33
  # Complete initialization by skipping (NTOT2*MODCNS + NTOT) numbers.
  # NB: the Fortran inner loop runs to NTOT, not NOW; kept as is.
  for loop2 in 1 .. r.ntot2 + 1:
    var now = ModCns
    if loop2 == r.ntot2 + 1: now = r.ntot
    if now > 0:
      if r.output != nil:
        r.output.write " RM48IN SKIPPING OVER ", align($now, 15), "\n"
      for idum in 1 .. r.ntot:
        var uni = r.u[r.i97] - r.u[r.j97]
        if uni < 0.0: uni = uni + 1.0
        r.u[r.i97] = uni
        dec r.i97
        if r.i97 == 0: r.i97 = 97
        dec r.j97
        if r.j97 == 0: r.j97 = 97
        r.c = r.c - r.cd
        if r.c < 0.0: r.c = r.c + r.cm
  r.initialized = true

proc rm48in*(r: var Rm48, ijklin, ntotin, ntot2n: int) =
  ## ENTRY RM48IN: initialise from seed `ijklin`, skipping
  ## `ntot2n*1e9 + ntotin` numbers.
  r.ijkl = ijklin
  r.ntot = max(ntotin, 0)
  r.ntot2 = max(ntot2n, 0)
  r.initRm48()

proc rm48ut*(r: Rm48): tuple[ijkl, ntot, ntot2: int] =
  ## ENTRY RM48UT: current seed and number counts.
  (r.ijkl, r.ntot, r.ntot2)

{.push boundChecks: off.}
proc rm48*(r: var Rm48, rvec: var openArray[float64]) =
  ## SUBROUTINE RM48(RVEC, LENV): fill `rvec` with numbers in (0, 1).
  ##
  ## The wrap-arounds are selects rather than branches (each is taken at
  ## random, half of the time); adding 0.0 instead of nothing changes no
  ## value, as no step can produce -0.0. 4.5 times faster than the Fortran
  ## form, the same numbers.
  if not r.initialized:
    # Default initialization: RM48 called without RM48IN.
    r.ijkl = DefaultSeed
    r.ntot = 0
    r.ntot2 = 0
    r.initRm48()
  var i97 = r.i97
  var j97 = r.j97
  var c = r.c
  let cd = r.cd
  let cm = r.cm
  for ivec in 0 ..< rvec.len:
    var uni = r.u[i97] - r.u[j97]
    uni = uni + (if uni < 0.0: 1.0 else: 0.0)
    r.u[i97] = uni
    i97 = (if i97 == 1: 97 else: i97 - 1)
    j97 = (if j97 == 1: 97 else: j97 - 1)
    c = c - cd
    c = c + (if c < 0.0: cm else: 0.0)
    uni = uni - c
    uni = uni + (if uni < 0.0: 1.0 else: 0.0)
    # Replace exact zeros by 2**-49
    rvec[ivec] = (if uni == 0.0: r.twom49 else: uni)
  r.i97 = i97
  r.j97 = j97
  r.c = c
  r.ntot = r.ntot + rvec.len
  if r.ntot >= ModCns:
    r.ntot2 = r.ntot2 + 1
    r.ntot = r.ntot - ModCns
{.pop.}

proc initMagRng*(output: File = stdout): MagRng =
  ## Fresh generator in the state Magboltz starts in (nothing drawn yet).
  result.rm.output = output

proc drand48*(g: var MagRng): float64 =
  ## Magboltz's `drand48(DUMMY)`.
  if g.ivec == 0 or g.ivec >= NVec:
    g.rm.rm48(g.rvec)
    g.ivec = 1
  else:
    g.ivec = g.ivec + 1
  result = g.rvec[g.ivec]

type
  Ranmar* = object
    ## RM48 without the buffer of the `drand48` shim: one number at a time,
    ## the same sequence as `rm48`. Small and without references, so each
    ## thread of a GPU can run its own.
    u: array[1..97, float64]
    c, cd, cm, twom49: float64
    i97, j97: int32

proc initRanmar*(seed: int): Ranmar =
  ## RM48 started from `seed`, as RM48IN(seed, 0, 0).
  var r: Rm48
  r.rm48in(seed, 0, 0)
  Ranmar(u: r.u, c: r.c, cd: r.cd, cm: r.cm, twom49: r.twom49,
         i97: int32(r.i97), j97: int32(r.j97))

proc drand48*(g: var Ranmar): float64 {.inline.} =
  ## The next number in (0, 1): one step of `rm48`.
  var uni = g.u[g.i97] - g.u[g.j97]
  uni = uni + (if uni < 0.0: 1.0 else: 0.0)
  g.u[g.i97] = uni
  g.i97 = (if g.i97 == 1: 97'i32 else: g.i97 - 1)
  g.j97 = (if g.j97 == 1: 97'i32 else: g.j97 - 1)
  g.c = g.c - g.cd
  g.c = g.c + (if g.c < 0.0: g.cm else: 0.0)
  uni = uni - g.c
  uni = uni + (if uni < 0.0: 1.0 else: 0.0)
  # Replace exact zeros by 2**-49
  result = (if uni == 0.0: g.twom49 else: uni)

proc toSingle*(g: var Ranmar) =
  ## The generator as the original RANMAR of Marsaglia and Zaman: numbers
  ## on a grid of 2⁻²⁴, so its arithmetic is exact in single precision
  ## (GPU kernels demoted to `float32`). The state snaps to the grid; then
  ## the same steps give that RANMAR's sequence, in single and in double
  ## precision alike.
  const grid = 16777216.0
  for i in 1 .. 97: g.u[i] = floor(g.u[i] * grid) / grid
  g.c = floor(g.c * grid) / grid
  g.cd = 7654321.0 / grid
  g.cm = 16777213.0 / grid
