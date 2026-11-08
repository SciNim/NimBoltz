## Minimal 3-vector for electron and gas-molecule kinematics: directions
## (`Vec3`, plain numbers) and vectors with a unit (`Vector[m•ps⁻¹]`, ...).
##
## Every operation is component-wise and evaluated in the order written, so
## results match the scalar Fortran expressions they replace bit for bit
## (e.g. `v / s` divides each component; it does not multiply by 1/s).

type
  Vector*[T] = object
    x*, y*, z*: T
  Vec3* = Vector[float64]

func vec3*[T](x, y, z: T): Vector[T] {.inline.} = Vector[T](x: x, y: y, z: z)

func `+`*[T](a, b: Vector[T]): Vector[T] {.inline.} =
  Vector[T](x: a.x + b.x, y: a.y + b.y, z: a.z + b.z)
func `-`*[T](a, b: Vector[T]): Vector[T] {.inline.} =
  Vector[T](x: a.x - b.x, y: a.y - b.y, z: a.z - b.z)
func plain[U](x: U): auto {.inline.} =
  ## a ratio of equal units (Unchained's `UnitLess`) as a plain number
  when $U == "UnitLess": float64(x) else: x

func `*`*[T, S](a: Vector[T], s: S): auto {.inline.} =
  vec3(plain(a.x * s), plain(a.y * s), plain(a.z * s))
func `/`*[T, S](a: Vector[T], s: S): auto {.inline.} =
  vec3(plain(a.x / s), plain(a.y / s), plain(a.z / s))

func norm2*[T](a: Vector[T]): auto {.inline.} =
  ## squared length, summed x, y, z
  a.x * a.x + a.y * a.y + a.z * a.z

template convert*[T](v: Vector[T], U: untyped): untyped =
  ## each component converted to the unit U
  let w = v
  vec3(w.x.to(U), w.y.to(U), w.z.to(U))
