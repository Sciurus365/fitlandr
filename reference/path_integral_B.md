# Bhattacharya method for path integration

A method to construct potential landscapes using path integration. See
references for details.

## Usage

``` r
path_integral_B(
  f,
  lims,
  n_path_int = 20,
  stepsize = 0.01,
  tol = 0.01,
  numTimeSteps = 1400,
  ...
)
```

## Arguments

- f:

  The vector field function. It should return `c(<dx/dt>, <dy/dt>)`

- lims:

  The limits of the range for the estimation as
  `c(<xl>, <xu>, <yl>, <yu>)`.

- n_path_int:

  The number of equally spaced points in each axis, at which the path
  integrals is to be calculated.

- stepsize:

  The time step used in each iteration.

- tol:

  The tolerance to test convergence.

- numTimeSteps:

  Number of time steps for integrating along each path (to ensure
  uniform arrays). Choose high-enough number for convergence with given
  stepsize.

- ...:

  Not in use.

## Value

A list with the following elements:

- `numPaths` Integer. Total Number of paths for defined grid spacing.

- `pot_path` Matrix. Potential along the paths.

- `path_tag` Vector. Tag for given paths.

- `attractors_pot` Vector. Potential value of each identified attractor
  by the path integral approach.

- `x_path` Vector. x-coord. along path.

- `y_path` Vector. y-coord. along path.

## References

Bhattacharya, S., Zhang, Q., & Andersen, M. E. (2011). A deterministic
map of Waddington’s epigenetic landscape for cell fate specification.
BMC Systems Biology, 5(1), 85. https://doi.org/10.1186/1752-0509-5-85.
The functions in this file were translated from the Matlab code provided
with the reference above, and its Python translation at
https://dynamo-release.readthedocs.io/en/v0.95.2/\_modules/dynamo/vectorfield/Bhattacharya.html
