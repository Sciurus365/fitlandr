# Simulation from vector fields

Parallel computing based on `future` is supported. Use
`future::plan("multisession")` to enable this.

## Usage

``` r
sim_vf(
  vf,
  noise = 1,
  noise_warmup = noise,
  chains = 10,
  length = 10000,
  discard = 0.3,
  stepsize = 0.01,
  sparse = 1,
  forbid_overflow = FALSE,
  linear_interp = FALSE,
  inits = matrix(c(stats::runif(chains, min = vf$lims[1], max = vf$lims[2]),
    stats::runif(chains, min = vf$lims[3], max = vf$lims[4])), ncol = 2)
)
```

## Arguments

- vf:

  A `vectorfield` object estimated by
  [`fit_2d_vf()`](https://sciurus365.github.io/fitlandr/reference/fit_2d_vf.md).

- noise:

  Relative noise of the simulation. Set this smaller when the simulation
  is unstable (e.g., when the elements in the diffusion matrix are not
  finite), and set this larger when the simulation converges too slowly.

- noise_warmup:

  The noise used for the warming-up period.

- chains:

  How many chains simulations should be performed?

- length:

  The simulation length for each chain.

- discard:

  How much of the starting part of each chain should be discarded?
  (Warming-up period.)

- stepsize:

  The stepsize for Euler–Maruyama simulation of the system.

- sparse:

  A number. How much do you want to sparse the output? When the noise is
  small, sparse the output may make the density estimation more
  efficient.

- forbid_overflow:

  If `TRUE`, when the simulated system runs out of the margins specified
  in `vf`, the system will be moved back to the previous value. This can
  help to stabilize the simulation. `FALSE` by default.

- linear_interp:

  Use linear interpolation method to estimate the drift vector (and the
  diffusion matrix). This can speed up the calculation. If `TRUE`, be
  sure that a linear grid was calculated for the vector field using
  `<vf> <- add_interp_grid(<vf>)`.

- inits:

  The initial values of each chain.

## Value

A matrix of the simulated data.
