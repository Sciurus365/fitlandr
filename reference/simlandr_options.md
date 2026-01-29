# Options controlling the landscape construction

To control the behavior of
[`simlandr::make_3d_static()`](https://sciurus365.github.io/simlandr/reference/make_3d_static.html),
but with default values accommodated for `fitlandr`. See
[`simlandr::make_3d_static()`](https://sciurus365.github.io/simlandr/reference/make_3d_static.html)
for details.

## Usage

``` r
simlandr_options(
  vf,
  x = rlang::expr(vf$x),
  y = rlang::expr(vf$y),
  lims = rlang::expr(vf$lims),
  kde_fun = c("ks", "MASS"),
  n = 200,
  adjust = 1,
  h,
  Umax = 5
)
```

## Arguments

- vf:

  A `vectorfield` object estimated by
  [`fit_2d_vf()`](https://sciurus365.github.io/fitlandr/reference/fit_2d_vf.md).

- x, y:

  The names of the target variables.

- lims:

  The limits of the range for the density estimator as `c(xl, xu)` for
  2D landscapes, `c(xl, xu, yl, yu)` for 3D landscapes,
  `c(xl, xu, yl, yu, zl, zu)` for 4D landscapes. If missing, the range
  of the data extended by 10% for both sides will be used. For
  landscapes based on multiple simulations, the largest range of all
  simulations (which means the lowest lower limit and the highest upper
  limit) will be used by default.

- kde_fun:

  Which kernel estimator to use? Choices: "ks"
  [`ks::kde()`](https://rdrr.io/pkg/ks/man/kde.html) (default; faster
  and using less memory); "base" `base::density()` (only for 2D
  landscapes); "MASS"
  [`MASS::kde2d()`](https://rdrr.io/pkg/MASS/man/kde2d.html) (only for
  3D landscapes).

- n:

  The number of equally spaced points in each axis, at which the density
  is to be estimated.

- adjust:

  The multiplier to the bandwidth. The bandwidth used is actually
  `adjust * h`. This makes it easy to specify values like "half the
  default" bandwidth.

- h:

  A number, or possibly a vector for 3D and 4D landscapes, specifying
  the smoothing bandwidth to be used. If missing, the default value of
  the kernel estimator will be used (but `bw = "SJ"` for
  `base::density()`). Note that the definition of bandwidth might be
  different for different kernel estimators. For landscapes based on
  multiple simulations, the largest `h` of all simulations will be used
  by default.

- Umax:

  The maximum displayed value of potential.

## Value

A list containing the parameters of the corresponding function. Only
intended to be used within
[`fit_3d_vfld()`](https://sciurus365.github.io/fitlandr/reference/fit_3d_vfld.md)
