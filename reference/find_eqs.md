# Find equilibrium points for a vector field

Find equilibrium points for a vector field

## Usage

``` r
find_eqs(vf, starts, jacobian_params = list(), ...)
```

## Arguments

- vf:

  A `vectorfield` object estimated by
  [`fit_2d_vf()`](https://sciurus365.github.io/fitlandr/reference/fit_2d_vf.md).

- starts:

  A vector indicating the starting value for solving the equilibrium
  point, or a list of vectors providing multiple starting values
  together.

- jacobian_params:

  Parameters passed to
  [`numDeriv::jacobian()`](https://rdrr.io/pkg/numDeriv/man/jacobian.html).

- ...:

  Parameters passed to
  [`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html).

## Value

A list of equilibrium points and their details. Use
`print.vectorfield_eqs()` to inspect it.
