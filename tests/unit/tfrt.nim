## Unit tests of the Fortran runtime helpers (port/frt.nim) against values
## gfortran produces. Run: nim c -r --path:port tests/unit/tfrt.nim
import std/unittest
import frt

suite "frt":
  test "integer powers (libgcc __powidf2 / __powisf2)":
    let x = 10.0 / 10.0025
    check powi(x, 3) == x * x * x
    check powi(1.5, 4) == 5.0625
    check powi(2.0'f32, 3) == 8.0'f32
  test "NaN / Inf in F and D editing":
    let z = 0.0
    check fmtF(z / z, 8, 2) == "     NaN"
    check fmtF(1.0 / z, 8, 2) == "Infinity"
    check fmtF(-1.0 / z, 8, 2) == "    -Inf"
    check fmtF(-1.0 / z, 3, 0) == "***"
    check fmtE(-1.0 / z, 10, 3, 'D') == " -Infinity"
