## Minimal 3-vector for electron and gas-molecule kinematics.
##
## Every operation is component-wise and evaluated in the order written, so
## results match the scalar Fortran expressions they replace bit for bit
## (e.g. `v / s` divides each component; it does not multiply by 1/s).

type
  Vec3* = object
    x*, y*, z*: float64

func vec3*(x, y, z: float64): Vec3 {.inline.} = Vec3(x: x, y: y, z: z)

func `+`*(a, b: Vec3): Vec3 {.inline.} = Vec3(x: a.x + b.x, y: a.y + b.y, z: a.z + b.z)
func `-`*(a, b: Vec3): Vec3 {.inline.} = Vec3(x: a.x - b.x, y: a.y - b.y, z: a.z - b.z)
func `*`*(a: Vec3, s: float64): Vec3 {.inline.} = Vec3(x: a.x * s, y: a.y * s, z: a.z * s)
func `/`*(a: Vec3, s: float64): Vec3 {.inline.} = Vec3(x: a.x / s, y: a.y / s, z: a.z / s)

func norm2*(a: Vec3): float64 {.inline.} =
  ## squared length, summed x, y, z
  a.x * a.x + a.y * a.y + a.z * a.z
