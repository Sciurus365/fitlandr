# Align potential values

So that all path-potentials end up at same global min and then generate
potential surface with interpolation on a grid.

## Usage

``` r
align_pot_B(resultB, n = 200, digits = 2, linear = TRUE, ...)
```

## Arguments

- resultB:

  Result from
  [`path_integral_B()`](https://sciurus365.github.io/fitlandr/reference/path_integral_B.md).

- n:

  The number of equally spaced points in each axis, at which the
  landscape is to be estimated.

- digits:

  Currently, the raw sample points in some regions are too dense that
  may crashes interpolation. To avoid this problem, only one point of
  all with the same first several digits. is kept. Use this parameter to
  indicate how many digits are considered. Note that this is a temporary
  solution and might be changed in the near future.

- linear:

  logical – indicating whether linear or spline interpolation should be
  used.

- ...:

  Other parameters passed to
  [`akima::interp()`](https://rdrr.io/pkg/akima/man/interp.html)

## Value

list with 3 components:

- x,y:

  vectors of x- and y- coordinates of output grid, the same as the input
  argument `xo`, or `yo`, if present. Otherwise, their default, a vector
  40 points evenly spaced over the range of the input `x`.

- z:

  matrix of fitted z-values. The value `z[i,j]` is computed at the x,y
  point `xo[i], yo[j]`. `z` has dimensions `length(xo)` times
  `length(yo)`.

If input is a `SpatialPointsDataFrame` a `SpatialPixelssDataFrame` is
returned.
