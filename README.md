
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `fitlandr`: Fit Vector Fields and Potential Landscapes from Intensive Longitudinal Data <img src='man/figures/logo.png' align="right" height="138" />

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/fitlandr)](https://cran.r-project.org/package=fitlandr)
![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
[![](https://cranlogs.r-pkg.org/badges/fitlandr)](https://cran.r-project.org/package=fitlandr)
<!-- badges: end -->

A toolbox for estimating vector fields from intensive longitudinal data,
and construct potential landscapes thereafter. The vector fields can be
estimated with two nonparametric methods: the Sparse Vector Field
Consensus (SparseVFC) algorithm by Ma et al. (2013)
<doi:10.1016/j.patcog.2013.05.017> and the Multivariate Vector Field
Kernel Estimator (MVKE) by Bandi & Moloche (2018)
<doi:10.1017/S0266466617000305>. The potential landscapes are
constructed with a simulation-based approach with the `simlandr`
package.

## Installation

You can install the development version of `fitlandr` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Sciurus365/fitlandr")
```

## Example

We use the following bistable dynamic system to illustrate the use of
`fitlandr`. The test data set is created as follows.

``` r
single_output_grad <- simlandr::sim_fun_grad(length = 200, seed = 1614)

library(tidyverse)
ggplot(data = single_output_grad %>% as_tibble()) +
    geom_path(aes(x = 1:200, y = x), color = "blue") +
    geom_path(aes(x = 1:200, y = y), color = "red") +
    theme_bw()
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

Fit the vector field of the system with VFC (see `?SparseVFC`
for the explanations of parameters):

``` r
library(fitlandr)
v1 <- fit_2d_vf(single_output_grad, x = "x", y = "y", method = "VFC", silent = FALSE, beta = 2)
#> Start mismatch removal...
#> iterate: 1th, gamma: 0.900000, the energy change rate: 1.003978, sigma2=0.131801
#> iterate: 2th, gamma: 0.929648, the energy change rate: 0.390302, sigma2=0.047712
#> iterate: 3th, gamma: 0.929648, the energy change rate: 0.057607, sigma2=0.036509
#> iterate: 4th, gamma: 0.924623, the energy change rate: 0.010792, sigma2=0.033713
#> iterate: 5th, gamma: 0.924623, the energy change rate: 0.001748, sigma2=0.032990
#> iterate: 6th, gamma: 0.924623, the energy change rate: 0.000375, sigma2=0.032837
#> iterate: 7th, gamma: 0.924623, the energy change rate: 0.000084, sigma2=0.032802
#> iterate: 8th, gamma: 0.924623, the energy change rate: 0.000020, sigma2=0.032794
#> iterate: 9th, gamma: 0.924623, the energy change rate: 0.000005, sigma2=0.032792
#> Removing outliers succesfully completed.
plot(v1)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Fit the potential landscape:

``` r
future::plan("multisession")
set.seed(1614)
l1 <- fit_3d_vfld(v1, .sim_vf_options = sim_vf_options(chains = 16, stepsize = 1, forbid_overflow = TRUE), .simlandr_options = simlandr_options(adjust = 5, Umax = 5))
#> ℹ Simulating the model✔ Simulating the model [19s]
#> ℹ Constructing the landscape✔ Constructing the landscape [2s]
plot(l1, 2)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r
# equivalent:
# s1 <- sim_vf(v1, chains = 16, stepsize = 1, forbid_overflow = TRUE)
# l1 <- simlandr::make_3d_static(s1, x = "x", y = "y", lims = v1$lims, adjust = 5, Umax = 5)
```

Fit the vector field with MVKE (see `?MVKE` for the explanations of
parameters):

``` r
v2 <- fit_2d_vf(single_output_grad, x = "x", y = "y", method = "MVKE")
plot(v2)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

Fit the potential landscape:

``` r
set.seed(1614)
l2 <- fit_3d_vfld(v2, .sim_vf_options = sim_vf_options(chains = 16, stepsize = 1, forbid_overflow = TRUE), .simlandr_options = simlandr_options(adjust = 5, Umax = 4))
#> ℹ Simulating the model✔ Simulating the model [19.4s]
#> ℹ Constructing the landscape✔ Constructing the landscape [1s]
plot(l2, 2)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

``` r
# equivalent:
# s2 <- sim_vf(v2, chains = 16, stepsize = 1, forbid_overflow = TRUE)
# l2 <- simlandr::make_3d_static(s2, x = "x", y = "y", lims = v2$lims, adjust = 5, Umax = 4)
```
