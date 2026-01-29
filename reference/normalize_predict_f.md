# Return a normalized prediction function

Return a normalized prediction function

## Usage

``` r
normalize_predict_f(vf)
```

## Arguments

- vf:

  A `vectorfield` object estimated by
  [`fit_2d_vf()`](https://sciurus365.github.io/fitlandr/reference/fit_2d_vf.md).

## Value

A function that takes a vector `x` and returns a list of `v`, the drift
part, and `a`, the diffusion part.
