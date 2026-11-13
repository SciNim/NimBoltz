## xoshiro256+ (Blackman and Vigna, 2018): a small, fast generator for the
## walkers and the avalanches when Magboltz's numbers are not wanted
## (`WalkerRng`). 32 bytes of state against Ranmar's 800, about 1 ns per
## number. Its lowest bits are weak; doubles take the top 53, which are
## not. Seeded through splitmix64, as its authors recommend. GPUs take
## xoshiro128** or xoshiro128+ (16 bytes, 32-bit arithmetic).

import rm48

type
  Xoshiro256Plus* = object
    s: array[4, uint64]

  Xoshiro128StarStar* = object
    ## xoshiro128** for GPU walkers in double precision: 32-bit arithmetic
    ## only (WGSL has no 64-bit integers), 16 bytes of state; numbers with
    ## 32 random bits
    s0, s1, s2, s3: uint32

  Xoshiro128Plus* = object
    ## xoshiro128+ for GPU walkers in single precision: numbers with 23
    ## random bits, exact in a float32, so a GPU demoted to single precision
    ## draws the CPU's numbers (as RANMAR's 24 bits, `rm48.toSingle`)
    s0, s1, s2, s3: uint32

  RngKind* = enum
    rkXoshiro  ## xoshiro256+ (GPUs: xoshiro128** or, in single precision, xoshiro128+)
    rkRanmar   ## Magboltz's RM48 (RANMAR), one number at a time

  Xoshiro128* = Xoshiro128StarStar | Xoshiro128Plus

  WalkerRng* = Ranmar | Xoshiro256Plus | Xoshiro128
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

# (GPU code: shifts by constants only, no arrays, so the state stays in registers)

proc step(g: var Xoshiro128) {.inline.} =
  let t = g.s1 shl 9
  g.s2 = g.s2 xor g.s0
  g.s3 = g.s3 xor g.s1
  g.s1 = g.s1 xor g.s2
  g.s0 = g.s0 xor g.s3
  g.s2 = g.s2 xor t
  g.s3 = (g.s3 shl 11) or (g.s3 shr 21)

proc drand48*(g: var Xoshiro128StarStar): float64 {.inline.} =
  ## the next number in (0, 1): 32 bits, centred in their cell of 2⁻³²
  let x = g.s1 * 5'u32
  let r = ((x shl 7) or (x shr 25)) * 9'u32
  g.step()
  (float64(r) + 0.5) * (1.0 / 4294967296.0)

proc drand48*(g: var Xoshiro128Plus): float64 {.inline.} =
  ## the next number in (0, 1): the top 23 bits, centred in their cell of
  ## 2⁻²³ (k + 1/2 with k < 2²³ is exact in a float32)
  let r = (g.s0 + g.s3) shr 9
  g.step()
  (float64(r) + 0.5) * (1.0 / 8388608.0)

proc initXoshiro128[X: Xoshiro128](seed: int): X =
  var x = uint64(seed)
  let a = splitmix64(x)
  let b = splitmix64(x)
  X(s0: uint32(a and 0xffffffff'u64), s1: uint32(a shr 32),
    s2: uint32(b and 0xffffffff'u64), s3: uint32(b shr 32))

proc seedOf*(g: var (Ranmar | Xoshiro256Plus)): int {.inline.} =
  ## a seed for another generator, from g's numbers
  int(drand48(g) * 9007199254740992.0)

proc seedOf*(g: var Xoshiro128): int {.inline.} =
  ## from two numbers (one has 32 or 23 bits only)
  let a = int(drand48(g) * 4294967296.0)
  let b = int(drand48(g) * 4294967296.0)
  (a shl 21) xor b

proc initRng*(R: typedesc[Ranmar], seed: int): Ranmar = initRanmar(seed)
proc initRng*(R: typedesc[Xoshiro256Plus], seed: int): Xoshiro256Plus = initXoshiro(seed)
proc initRng*(R: typedesc[Xoshiro128StarStar], seed: int): Xoshiro128StarStar =
  initXoshiro128[Xoshiro128StarStar](seed)
proc initRng*(R: typedesc[Xoshiro128Plus], seed: int): Xoshiro128Plus =
  initXoshiro128[Xoshiro128Plus](seed)

template withRng*(kind: RngKind, R, body: untyped) =
  ## body with `R` the generator type of `kind`
  case kind
  of rkXoshiro:
    type R {.inject.} = Xoshiro256Plus
    body
  of rkRanmar:
    type R {.inject.} = Ranmar
    body
