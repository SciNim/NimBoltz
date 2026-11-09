## Data in single precision, for GPU kernels demoted to `float32`
## (Crucible's `cudaSingle`, `webgpuSingle`): the same objects with every
## `float64` (and every unit of Unchained) as a `float32`. With 4-byte
## fields only, an object has no padding on either side, so a packed
## object is the device's.

import std/typetraits
import rm48

proc pack32*[T](x: T, buf: var seq[byte])
proc unpack32*[T](x: var T, buf: openArray[byte], pos: var int)

proc pack32*(x: Ranmar, buf: var seq[byte]) =
  var g = x
  g.toSingle()  # exact in single precision
  for _, f in g.fieldPairs: pack32(f, buf)

proc add4(buf: var seq[byte], p: pointer) =
  let n = buf.len
  buf.setLen(n + 4)
  copyMem(addr buf[n], p, 4)

proc pack32*[T](x: T, buf: var seq[byte]) =
  ## appends x in single precision
  when T is Ranmar:
    pack32(Ranmar(x), buf)
  elif distinctBase(T) is float64:
    var s = float32(distinctBase(T)(x))
    buf.add4(addr s)
  elif T is int32 | uint32 | float32:
    var s = x
    buf.add4(addr s)
  elif T is array:
    for e in x: pack32(e, buf)
  elif T is object:
    for _, f in x.fieldPairs: pack32(f, buf)
  else:
    {.error: "pack32: no single-precision form of " & $T.}

proc unpack32*[T](x: var T, buf: openArray[byte], pos: var int) =
  ## x from single precision at buf[pos], then pos after it
  when distinctBase(T) is float64:
    var s: float32
    copyMem(addr s, unsafeAddr buf[pos], 4)
    x = T(float64(s))
    pos = pos + 4
  elif T is int32 | uint32 | float32:
    copyMem(addr x, unsafeAddr buf[pos], 4)
    pos = pos + 4
  elif T is array:
    for e in x.mitems: unpack32(e, buf, pos)
  elif T is object:
    for _, f in x.fieldPairs: unpack32(f, buf, pos)
  else:
    {.error: "unpack32: no single-precision form of " & $T.}

proc single*[T](x: T): seq[byte] =
  ## x in single precision
  pack32(x, result)

proc single*[T](xs: openArray[T]): seq[byte] =
  ## xs in single precision, one after the other
  for x in xs: pack32(x, result)

proc singleSize*(T: typedesc): int =
  ## bytes of a T in single precision
  var x: T
  single(x).len
