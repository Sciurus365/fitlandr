# Multivariate vector field kernel estimator

See references for details.

## Usage

``` r
MVKE(d, v, h = 0.2, kernel = c("Gaussian", "exp"))
```

## Arguments

- d:

  The dataset. Should be a matrix or a data frame, with each row
  representing a random vector.

- v:

  The vectors corresponding to the dataset. Should be a matrix or a data
  frame with the same shape as `d`. If missing, then the vectors will be
  calculated from the dataset.

- h:

  The bandwidth for the kernel estimator.

- kernel:

  The type of kernel estimator used. "Gaussian" by default.

## Value

A function(x), which then returns the \\\mu\\ and \\a\\ estimators at
the position \\x\\.

## References

Bandi, F. M., & Moloche, G. (2018). On the functional estimation of
multivariate diffusion processes. Econometric Theory, 34(4), 896-946.
https://doi.org/10.1017/S0266466617000305
