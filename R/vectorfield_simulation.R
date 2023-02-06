eval_pass_missing <- function(expr, ...) {
  if (rlang::is_missing(expr)) {
    return(rlang::missing_arg())
  }
  return(rlang::eval_tidy(expr, ...))
}


#' Simulation from vector fields
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
#' @export
sim_vf <- function(vf, noise = 1, noise_warmup = noise, chains = 10, length = 1e4, discard = 0.3, stepsize = 0.01, sparse = 1, forbid_overflow = FALSE, linear_interp = FALSE, inits = matrix(c(
                     stats::runif(chains, min = vf$lims[1], max = vf$lims[2]),
                     stats::runif(chains, min = vf$lims[3], max = vf$lims[4])
                   ), ncol = 2)) {
  f <- function(x) {
  	stats::predict(object = vf, pos = x, linear_interp = linear_interp, calculate_a = TRUE)
  }

  force(inits)
  result <- future.apply::future_apply(inits, MARGIN = 1, FUN = sim_vf_single, f = f, length = length, noise = noise, noise_warmup = noise_warmup, lims = vf$lims, forbid_overflow = forbid_overflow, stepsize = stepsize, sparse = sparse, discard = discard, simplify = FALSE, future.seed = TRUE, future.packages = "SparseVFC")
  result <- do.call(rbind, result)
  colnames(result) <- colnames(vf$data)
  return(result)
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

  result <- result[seq(1, nrow(result), by = sparse), ]

  return(result)
}



#' Options controlling the vector field simulation
#' See [sim_vf()] for details.
#' @inheritParams fit_3d_vfld
#' @inheritParams sim_vf
#' @export
sim_vf_options <- function(vf, noise = 1,noise_warmup = noise, chains = 10, length = 1e4, discard = 0.3, stepsize = 0.01, sparse = 1, forbid_overflow = FALSE, inits = rlang::expr(matrix(c(
                             stats::runif(chains, min = vf$lims[1], max = vf$lims[2]),
                             stats::runif(chains, min = vf$lims[3], max = vf$lims[4])
                           ), ncol = 2))) {
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
#' Then [simlandr::check_conv()] can be used meaningfully.
#'
#' @param s A simulation output, possibly generated by [sim_vf()]
#' @inheritParams sim_vf
#'
#' @return A reordered matrix of the simulation output.
#' @export
reorder_output <- function(s, chains) {
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
