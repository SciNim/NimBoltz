## Fortran-style arrays: column major, 1-based, fixed size.
##
## Used by the collision tables and the gas cross sections.

type
  Arr2*[T; N1, N2: static int] = object
    ## Fortran array A(N1, N2), column major, 1-based: `a[i, j]`
    data*: array[1..N2, array[1..N1, T]]
  Arr3*[T; N1, N2, N3: static int] = object
    ## Fortran array A(N1, N2, N3), column major, 1-based: `a[i, j, k]`
    data*: array[1..N3, array[1..N2, array[1..N1, T]]]

template `[]`*(a: Arr2, i, j: SomeInteger): untyped = a.data[j][i]
template `[]=`*(a: var Arr2, i, j: SomeInteger, v: untyped) = a.data[j][i] = v
template `[]`*(a: Arr3, i, j, k: SomeInteger): untyped = a.data[k][j][i]
template `[]=`*(a: var Arr3, i, j, k: SomeInteger, v: untyped) = a.data[k][j][i] = v
