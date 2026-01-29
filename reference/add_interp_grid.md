# Add a grid to a `vectorfield` object to enable linear interpolation

Add a grid to a `vectorfield` object to enable linear interpolation

## Usage

``` r
add_interp_grid(vf, lims = vf$lims, n = vf$n)
```

## Arguments

- vf:

  A `vectorfield` object estimated by
  [`fit_2d_vf()`](https://sciurus365.github.io/fitlandr/reference/fit_2d_vf.md).

- lims:

  The limits of the range for the vector field estimation as
  `c(<xl>, <xu>, <yl>, <yu>)`. If missing, the range of the data
  extended by 10% for both sides will be used.

- n:

  The number of equally spaced points in each axis, at which the vectors
  are to be estimated.

## Value

A `vectorfield` project with an `interp_grid` field.
