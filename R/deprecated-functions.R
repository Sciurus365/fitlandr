# Deprecated fitlandr functions
#
# This file centralizes legacy functions that are retained only for backward
# compatibility. New development should use the current `make_*` / `autoplot()`
# interfaces instead.

#' @noRd
eval_pass_missing <- function(expr, ...) {
  if (rlang::is_missing(expr)) {
    return(rlang::missing_arg())
  }
  rlang::eval_tidy(expr, ...)
}

sim_vf_single <- function(init, f, length, noise, noise_warmup, stepsize, discard, sparse, lims, forbid_overflow) {
  dim <- length(init)
  prev <- init
  for (i in 1:(length * discard)) {
    dyn <- f(prev)
    next_point <- prev + stepsize * dyn$v + sqrt(stepsize) * MASS::mvrnorm(mu = rep(0, dim), Sigma = dyn$a * noise_warmup)
    if (forbid_overflow) {
      if (!(dplyr::between(next_point[1], lims[1], lims[2]) & dplyr::between(next_point[2], lims[3], lims[4]))) {
        next_point <- prev
      }
    }
    prev <- next_point
  }

  result <- matrix(NA_real_, nrow = length * (1 - discard), ncol = dim)
  result[1, ] <- prev
  for (i in 2:(length * (1 - discard))) {
    prev <- result[i - 1, ]
    dyn <- f(prev)
    next_point <- prev + stepsize * dyn$v + sqrt(stepsize) * MASS::mvrnorm(mu = rep(0, dim), Sigma = dyn$a * noise)
    if (forbid_overflow) {
      if (!(dplyr::between(next_point[1], lims[1], lims[2]) & dplyr::between(next_point[2], lims[3], lims[4]))) {
        next_point <- prev
      }
    }
    result[i, ] <- next_point
  }

  result[seq(1, nrow(result), by = sparse), ]
}

#' Estimate a 2D potential landscape from data with the MVKE method
#'
#' This function is a wrapper of the MVKE method (see [MVKE()]) that produces a 2D potential landscape from 1D data. The landscape is constructed by estimating the gradient of the data and then integrating it. The MVKE method is a non-parametric method that estimates the gradient of the data by using a kernel density estimator. The potential landscape is then constructed by integrating the gradient.
#'
#' @param data A data frame or matrix containing the data. The data frame should contain at least a column, with the column name indicated by `x`, that represents the dimension for landscape construction.
#' @param x The column name of the data frame that represents the dimension for landscape construction.
#' @param lims The limits of the range for the landscape calculation as `c(xl, xu)`.
#' @param n The number of equally spaced points in the axis, at which the landscape is to be estimated.
#' @param method The method used to estimate the gradient. Currently only "MVKE" is supported.
#' @param ... Additional arguments passed to [MVKE()]. (Not used for the `summary()` function).
#' @inheritParams stats::integrate
#' @inheritParams fit_2d_vf
#' @return A `2d_MVKE_landscape` object, which contains the following components:
#' \itemize{
#'   \item `dist`: A data frame containing the estimated potential landscape. The data frame has two columns: `x` and `U`, where `x` is the position and `U` is the potential.
#'   \item `plot`: A ggplot object containing the plot of the potential landscape.
#' }
#' @export
#'
#' @examples
#' # generate data
#' single_output_grad <- simlandr::sim_fun_grad(length = 200, seed = 1614)
#' # fit the landscape
#' l <- fit_2d_ld(single_output_grad, "x")
#'
#' summary(l)
#' autoplot(l)
#'
#' # different behaviors for different `na_action` choices
#'
#' l1 <- fit_2d_ld(data.frame(x = c(1, 2, 1, 2, NA, NA, NA, 10, 11, 10, 11)), "x")
#' autoplot(l1)
#'
#' l2 <- fit_2d_ld(data.frame(x = c(1, 2, 1, 2, NA, NA, NA, 10, 11, 10, 11)), "x",
#'   na_action = "omit_vectors"
#' )
#' autoplot(l2)
fit_2d_ld <- function(data, x, lims, n = 200L, vector_position = "start", na_action = "omit_data_points",
                      dayvar = NULL, beepvar = NULL, method = c("MVKE"), subdivisions = 100L, rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol, stop.on.error = TRUE, keep.xy = FALSE, aux = NULL, ...) {
  lifecycle::deprecate_warn(
    "0.2.0",
    "fit_2d_ld()",
    I("fit_1d_vf() + make_1d_ld()")
  )

  if (!is.null(dayvar) && !dayvar %in% colnames(data)) {
    cli::cli_abort("{.arg dayvar} must name a column in {.arg data}.")
  }
  if (!is.null(beepvar) && !beepvar %in% colnames(data)) {
    cli::cli_abort("{.arg beepvar} must name a column in {.arg data}.")
  }

  d <- insert_time_separators(data, x, dayvar = dayvar, beepvar = beepvar)
  warn_if_timevars_ineffective(dayvar = dayvar, beepvar = beepvar, na_action = na_action)

  if (is.data.frame(d)) {
    d_raw <- d[, c(x), drop = FALSE] %>% as.matrix()
  } else if (is.matrix(d)) {
    d_raw <- d[, c(x), drop = FALSE]
  } else {
    cli::cli_abort("{.arg data} must be a data frame or a matrix.")
  }

  if (na_action != "omit_data_points" & na_action != "omit_vectors") {
    cli::cli_abort('{.arg na_action} must be either "omit_data_points" or "omit_vectors".')
  }

  if (na_action == "omit_data_points" & any(is.na(d_raw))) {
    d_raw <- stats::na.omit(d_raw)
    cli::cli_inform("NA(s) found in the data. Those data points were omitted.")
  }

  v_mat <- diff(d_raw)

  if (vector_position == "start") {
    x_mat <- d_raw[1:(nrow(d_raw) - 1), , drop = FALSE]
  } else if (vector_position == "middle") {
    x_mat <- d_raw[1:(nrow(d_raw) - 1), , drop = FALSE] + 0.5 * v_mat
  } else if (vector_position == "end") {
    x_mat <- d_raw[2:nrow(d_raw), , drop = FALSE]
  } else {
    cli::cli_abort('{.arg vector_position} must be one of "start", "middle", or "end".')
  }

  data_vectors <- cbind(x_mat, v_mat) %>%
    `colnames<-`(c("x", "vx"))

  if (any(is.na(data_vectors)) && na_action == "omit_vectors") {
    data_vectors <- stats::na.omit(data_vectors)
    cli::cli_inform("NA(s) found in the data. Those vectors were omitted.")
  }

  lims <- simlandr::determine_lims(data, x, lims)
  MVKEresult <- MVKE(data_vectors[, 1, drop = FALSE], data_vectors[, 2, drop = FALSE], ...)

  xseq <- seq(lims[1], lims[2], length.out = n)
  Useq <- vector("numeric", length = n)

  Useq[1] <- 0
  for (i in 2:n) {
    Useq[i] <- Useq[i - 1] - stats::integrate(function(x) purrr::map_dbl(x, function(xx) MVKEresult(xx)$mu), xseq[i - 1], xseq[i], subdivisions = subdivisions, rel.tol = rel.tol, abs.tol = abs.tol, stop.on.error = stop.on.error, keep.xy = keep.xy, aux = aux)$value
  }

  dist <- data.frame(x = xseq, U = Useq)
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = dist, ggplot2::aes(x = x, y = U)) +
    ggplot2::theme_bw()

  structure(list(dist = dist, plot = p, MVKEresult = MVKEresult), class = c("2d_MVKE_landscape", "landscape"))
}

#' @export
#' @describeIn fit_2d_ld Find the local minima of the 2D potential landscape
#' @param object An object of class `2d_MVKE_landscape` returned by [fit_2d_ld()].
#' @method summary 2d_MVKE_landscape
summary.2d_MVKE_landscape <- function(object, ...) {
  local_minima <- which(diff(sign(diff(object$dist$U))) == 2) + 1
  cli::cli_inform("{length(local_minima)} local minima were found.")
  data.frame(x = object$dist$x[local_minima], U = object$dist$U[local_minima])
}

#' Simulation from vector fields
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated in favor of the more efficient and accurate `make_2d_ld` approach.
#'
#' Parallel computing based on `future` is supported. Use `future::plan("multisession")` to enable this.
#'
#' @inheritParams fit_3d_vfld
#' @param noise Relative noise of the simulation. Set this smaller when the simulation is unstable (e.g., when the elements in the diffusion matrix are not finite), and set this larger when the simulation converges too slowly.
#' @param noise_warmup The noise used for the warming-up period.
#' @param stepsize The stepsize for Euler–Maruyama simulation of the system.
#' @param sparse A number. How much do you want to sparse the output? When the noise is small, sparse the output may make the density estimation more efficient.
#' @param chains How many chains simulations should be performed?
#' @param length The simulation length for each chain.
#' @param discard How much of the starting part of each chain should be discarded? (Warming-up period.)
#' @param forbid_overflow If `TRUE`, when the simulated system runs out of the margins specified in `vf`, the system will be moved back to the previous value. This can help to stabilize the simulation. `FALSE` by default.
#' @param inits The initial values of each chain.
#' @inheritParams predict.vectorfield
#' @return A matrix of the simulated data.
#'
#' @keywords internal
#' @export
sim_vf <- function(vf, noise = 1, noise_warmup = noise, chains = 10, length = 1e4, discard = 0.3, stepsize = 0.01, sparse = 1, forbid_overflow = FALSE, linear_interp = FALSE, inits = matrix(c(
                     stats::runif(chains, min = vf$lims[1], max = vf$lims[2]),
                     stats::runif(chains, min = vf$lims[3], max = vf$lims[4])
                   ), ncol = 2)) {
  lifecycle::deprecate_warn(
    "0.2.0",
    "sim_vf()",
    "make_2d_ld()"
  )
  rlang::check_installed(
    "future.apply",
    reason = "for deprecated simulation via {.fn sim_vf()}. Install it with {.code install.packages(\"future.apply\") }."
  )
  rlang::check_installed(
    "MASS",
    reason = "for deprecated simulation via {.fn sim_vf()}. Install it with {.code install.packages(\"MASS\") }."
  )

  f <- function(x) {
    stats::predict(object = vf, pos = x, linear_interp = linear_interp, calculate_a = TRUE)
  }

  force(inits)
  result <- future.apply::future_apply(inits, MARGIN = 1, FUN = sim_vf_single, f = f, length = length, noise = noise, noise_warmup = noise_warmup, lims = vf$lims, forbid_overflow = forbid_overflow, stepsize = stepsize, sparse = sparse, discard = discard, simplify = FALSE, future.seed = TRUE, future.packages = "SparseVFC")
  result <- do.call(rbind, result)
  colnames(result) <- colnames(vf$data)
  result
}

#' Options controlling the vector field simulation
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated in favor of the more efficient and accurate `make_2d_ld` approach.
#'
#' See [sim_vf()] for details.
#' @inheritParams fit_3d_vfld
#' @inheritParams sim_vf
#'
#' @return A list containing the parameters of the corresponding function. Only intended to be used within [fit_3d_vfld()]
#' @export
#' @keywords internal
sim_vf_options <- function(vf, noise = 1, noise_warmup = noise, chains = 10, length = 1e4, discard = 0.3, stepsize = 0.01, sparse = 1, forbid_overflow = FALSE, linear_interp = FALSE, inits = rlang::expr(matrix(c(
                             stats::runif(chains, min = vf$lims[1], max = vf$lims[2]),
                             stats::runif(chains, min = vf$lims[3], max = vf$lims[4])
                           ), ncol = 2))) {
  lifecycle::deprecate_warn(
    "0.2.0",
    "sim_vf_options()",
    "make_2d_ld()"
  )

  if (!missing(vf)) {
    return(list(vf = vf, noise = noise, chains = chains, length = length, discard = discard, stepsize = stepsize, sparse = sparse, forbid_overflow = forbid_overflow, inits = eval(inits)))
  } else {
    return(list(vf = rlang::expr(vf), noise = noise, chains = chains, length = length, discard = discard, stepsize = stepsize, sparse = sparse, forbid_overflow = forbid_overflow, inits = inits))
  }
}

#' Options controlling the landscape construction
#'
#' To control the behavior of [simlandr::make_3d_static()], but with default values accommodated for `fitlandr`. See [simlandr::make_3d_static()] for details.
#' @inheritParams fit_3d_vfld
#' @inheritParams simlandr::make_3d_static
#'
#' @inherit sim_vf_options return
#' @export
simlandr_options <- function(vf, x = rlang::expr(vf$x), y = rlang::expr(vf$y), lims = rlang::expr(vf$lims), kde_fun = c("ks", "MASS"), n = 200, adjust = 1, h, Umax = 5) {
  if (!missing(vf)) {
    return(list(x = eval(x), y = eval(y), lims = eval(lims), kde_fun = kde_fun, n = n, adjust = adjust, h = rlang::maybe_missing(h), Umax = Umax))
  } else {
    return(list(x = x, y = y, lims = lims, kde_fun = kde_fun, n = n, adjust = adjust, h = rlang::maybe_missing(h), Umax = Umax))
  }
}

#' Reorder a simulation output in time order
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated in favor of the more efficient and accurate [make_2d_ld()] approach.
#'
#' Then [simlandr::check_conv()] can be used meaningfully.
#'
#' @param s A simulation output, possibly generated by [sim_vf()]
#' @inheritParams sim_vf
#'
#' @return A reordered matrix of the simulation output.
#' @keywords internal
#' @export
reorder_output <- function(s, chains) {
  lifecycle::deprecate_warn(
    "0.2.0",
    "reorder_output()",
    "make_2d_ld()"
  )
  reorder_index <- vector("integer", nrow(s))
  current_pos <- 1
  for (i in 1:(nrow(s) / chains)) {
    for (j in 0:(chains - 1)) {
      reorder_index[current_pos] <- (nrow(s) / chains) * j + i
      current_pos <- current_pos + 1
    }
  }
  s[reorder_index, ]
}

#' Estimate a 3D potential landscape from a vector field
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated in favor of the more efficient and accurate `make_2d_ld` approach.
#'
#' Two methods are available: `method = "pathB"` and `method = "simlandr"`. See *Details* section.
#'
#' @details
#' For `method = "simlandr"`, the landscape is constructed based on the generalized potential landscape by Wang et al. (2008), implemented by the `simlandr` package. This function is a wrapper of [sim_vf()] and [simlandr::make_3d_static()]. Use those two functions separately for more customization.
#'
#' For `method = "pathB"`, the landscape is constructed based on the deterministic path-integral quasi-potential defined by Bhattacharya et al. (2011).
#'
#' We recommend the `simlandr` method for psychological data because it is more stable.
#'
#' Parallel computing based on `future` is supported for both methods. Use `future::plan("multisession")` to enable this and speed up computation.
#'
#' @param vf A `vectorfield` object estimated by [fit_2d_vf()].
#' @param method The method used for landscape construction. Can be `pathB` or `simlandr`.
#' @param .pathB_options Only for `method = "pathB"`. Options controlling the path-integral algorithm. Should be generated by [sim_vf_options()].
#' @param .sim_vf_options Only for `method = "simlandr"`. Options controlling the vector field simulation. Should be generated by [sim_vf_options()].
#' @param .simlandr_options Only for `method = "simlandr"`. Options controlling the landscape construction. Should be generated by [simlandr_options()].
#' @inheritParams predict.vectorfield
#'
#' @return A `landscape` object as described in [simlandr::make_3d_static()], or a `3d_static_landscape_B` object, which inherits from the `landscape` class and contains the following elements: `dist`, the distribution estimation for landscapes; `plot`, a 3D plot using `plotly`; plot_2, a 2D plot using `ggplot2`; x, y, from `vf`.
#' @examplesIf interactive()
#' # generate data
#' single_output_grad <- simlandr::sim_fun_grad(length = 200, seed = 1614)
#' # fit the vector field
#' v2 <- fit_2d_vf(single_output_grad, x = "x", y = "y", method = "MVKE")
#' autoplot(v2)
#' # fit the landscape
#' future::plan("multisession")
#' set.seed(1614)
#' l2 <- fit_3d_vfld(v2,
#'   .sim_vf_options = sim_vf_options(chains = 16, stepsize = 1, forbid_overflow = TRUE),
#'   .simlandr_options = simlandr_options(adjust = 5, Umax = 4)
#' )
#' autoplot(l2)
#' future::plan("sequential")
#' @export
#' @keywords internal
fit_3d_vfld <- function(vf, method = c("simlandr", "pathB"), .pathB_options = pathB_options(vf), .sim_vf_options = sim_vf_options(vf), .simlandr_options = simlandr_options(vf), linear_interp = FALSE) {
  lifecycle::deprecate_warn(
    "0.2.0",
    "fit_3d_vfld()",
    "make_2d_ld()"
  )

  method <- match.arg(method[1], c("pathB", "simlandr"))
  if (method == "pathB") {
    all_pars <- .pathB_options %>% lapply(eval_pass_missing, list(vf = vf))
    all_pars$f <- function(x) {
      stats::predict(all_pars$vf, pos = x, linear_interp = linear_interp, calculate_a = FALSE)$v
    }
    cli::cli_progress_step("Calculating path integrals")
    resultB <- do.call(path_integral_B, all_pars)
    cli::cli_progress_step("Aligning potentials")
    out_B <- do.call(align_pot_B, c(list(resultB = resultB), all_pars))
    out_B$d <- out_B$z

    p <- plotly::plot_ly(x = out_B$x, y = out_B$y, z = out_B$z, type = "surface")
    p <- plotly::layout(p, scene = list(
      xaxis = list(title = vf$x),
      yaxis = list(title = vf$y), zaxis = list(title = "U")
    )) %>%
      plotly::colorbar(title = "U")
    p2 <- ggplot2::ggplot(simlandr::make_2d_tidy_dist(out_B), ggplot2::aes(
      x = x,
      y = y
    )) +
      ggplot2::geom_raster(ggplot2::aes(fill = d)) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(x = vf$x, y = vf$y, fill = "U") +
      ggplot2::theme_bw()
    result <- c(list(
      dist = out_B, plot = p, plot_2 = p2, x = vf$x,
      y = vf$y
    ), all_pars)
    class(result) <- c(
      "3d_static_landscape_B", "3d_static_landscape", "3d_landscape",
      "landscape"
    )
    return(result)
  } else if (method == "simlandr") {
    cli::cli_progress_step("Simulating the model")
    simulation_output <- do.call(sim_vf, c(.sim_vf_options, list(linear_interp = linear_interp)) %>% lapply(eval_pass_missing, list(vf = vf, chains = .$chains)))
    cli::cli_progress_step("Constructing the landscape")
    return(do.call(simlandr::make_3d_static, c(list(output = simulation_output), .simlandr_options %>% lapply(eval_pass_missing, list(vf = vf)))))
  }
}

utils::globalVariables("U")
