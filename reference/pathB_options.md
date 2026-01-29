# Options controlling the path-integral algorithm

See
[`path_integral_B()`](https://sciurus365.github.io/fitlandr/reference/path_integral_B.md),
[`align_pot_B()`](https://sciurus365.github.io/fitlandr/reference/align_pot_B.md)
for details.

## Usage

``` r
pathB_options(
  vf,
  lims = rlang::expr(vf$lims),
  n_path_int = 20,
  stepsize = 0.01,
  tol = 0.01,
  numTimeSteps = 1400,
  n = 200,
  digits = 2,
  linear = TRUE,
  ...
)
```

## Arguments

- vf:

  A `vectorfield` object estimated by
  [`fit_2d_vf()`](https://sciurus365.github.io/fitlandr/reference/fit_2d_vf.md).

- lims:

  The limits of the range for the estimation as
  `c(<xl>, <xu>, <yl>, <yu>)`.

- n_path_int:

  The number of equally spaced points in each axis, at which the path
  integrals is to be calculated.

- stepsize:

  The stepsize for Euler–Maruyama simulation of the system.

- tol:

  The tolerance to test convergence.

- numTimeSteps:

  Number of time steps for integrating along each path (to ensure
  uniform arrays). Choose high-enough number for convergence with given
  stepsize.

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

  Not in use.

## Value

A list containing the parameters of the corresponding function. Only
intended to be used within
[`fit_3d_vfld()`](https://sciurus365.github.io/fitlandr/reference/fit_3d_vfld.md)
