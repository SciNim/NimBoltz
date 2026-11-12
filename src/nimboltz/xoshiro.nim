## xoshiro256+ (Blackman and Vigna, 2018): a small, fast generator for the
## walkers and the avalanches when Magboltz's numbers are not wanted
## (`WalkerRng`). 32 bytes of state against Ranmar's 800, about 1 ns per
## number. Its lowest bits are weak; doubles take the top 53, which are
## not. Seeded through splitmix64, as its authors recommend.

import rm48

type
  Xoshiro256Plus* = object
    s: array[4, uint64]

  RngKind* = enum
    rkXoshiro  ## xoshiro256+
    rkRanmar   ## Magboltz's RM48 (RANMAR), one number at a time

  WalkerRng* = Ranmar | Xoshiro256Plus
    ## the generators of walkers and of the avalanches with exact bounds

func rotl(x: uint64, k: int): uint64 {.inline.} = (x shl k) or (x shr (64 - k))

proc splitmix64(x: var uint64): uint64 =
  x = x + 0x9e3779b97f4a7c15'u64
  var z = x
  z = (z xor (z shr 30)) * 0xbf58476d1ce4e5b9'u64
  z = (z xor (z shr 27)) * 0x94d049bb133111eb'u64
  z xor (z shr 31)

proc initXoshiro*(seed: int): Xoshiro256Plus =
  var x = uint64(seed)
  for i in 0 .. 3: result.s[i] = splitmix64(x)

proc drand48*(g: var Xoshiro256Plus): float64 {.inline.} =
  ## the next number in (0, 1): the top 53 bits, centred in their cell of
  ## 2⁻⁵³ (never 0, as −ln r needs)
  let r = g.s[0] + g.s[3]
  let t = g.s[1] shl 17
  g.s[2] = g.s[2] xor g.s[0]
  g.s[3] = g.s[3] xor g.s[1]
  g.s[1] = g.s[1] xor g.s[2]
  g.s[0] = g.s[0] xor g.s[3]
  g.s[2] = g.s[2] xor t
  g.s[3] = rotl(g.s[3], 45)
  (float64(r shr 11) + 0.5) * (1.0 / 9007199254740992.0)

proc initRng*(R: typedesc[Ranmar], seed: int): Ranmar = initRanmar(seed)
proc initRng*(R: typedesc[Xoshiro256Plus], seed: int): Xoshiro256Plus = initXoshiro(seed)

template withRng*(kind: RngKind, R, body: untyped) =
  ## body with `R` the generator type of `kind`
  case kind
  of rkXoshiro:
    type R {.inject.} = Xoshiro256Plus
    body
  of rkRanmar:
    type R {.inject.} = Ranmar
    body
