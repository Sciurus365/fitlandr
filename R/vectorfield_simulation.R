#' Estimate a 2D potential landscape from a vector field
#'
#' This is done based on the generalized potential landscape by Wang et al (2008), implemented by the `simlandr` package.
fit_vfld_2d <- function(vf, method = "simlandr", .sim_vf_options = sim_vf_options(vf), .simlandr_options = simlandr_options(vf)) {

  simulation_output <- do.call(sim_vf, .sim_vf_options)
  return(do.call(simlandr::make_3d_static, c(list(output = simulation_output), .simlandr_options)))
}


#' Simulation from vector fields
#'
#' Parallel computing based on `future` is supported. Use `future::plan("multisession")` to enable this.
#'
#' @param chains How many chains simulations should be performed?
#' @param length The simulation length for each chain.
#' @param discard How much of the starting part of each chain should be discarded? (Warming-up period.)
#' @param inits The initial values of each chain.
sim_vf <- function(vf, chains = 10, length = 1e4, discard = 0.3, inits = matrix(c(
	runif(chains, min = vf$x_start, max = vf$x_end),
	runif(chains, min = vf$y_start, max = vf$y_end)
), ncol = 2)) {
	if (vf$method == "VFC") {
		f <- function(x) {
			v <- stats::predict(vf$VFCresult, x %>% normalize_v(vf$data_normalized)) %>% scale_up(vf$data_normalized)
			variance <- vf$VFCresult$sigma2 %>% scale_up2(vf$data_normalized)
			a <- matrix(c(variance, 0, 0, variance), nrow = 2)
			return(list(v = v, a = a))
		}
	} else if (vf$method == "MVKE") {
		f <- function(x) {
			result <- vf$MVKEresult(x %>% normalize_v(vf$data_normalized))
			v <- result$mu %>% scale_up(vf$data_normalized)
			a <- result$a %>% scale_up2(vf$data_normalized)
			return(list(v = v, a = a))
		}
	}

	force(inits)
  result <- future.apply::future_apply(inits, MARGIN = 1, FUN = sim_vf_single, f = f, length = length, simplify = FALSE, future.seed = TRUE, future.packages = "SparseVFC")
  # result <- furrr::future_map(inits %>% asplit(MARGIN = 1), .f = sim_vf_single, f = f, length = length, .options = furrr::furrr_options(seed = TRUE))
  result <- lapply(result, function(x) x[1:(nrow(x)*(1-discard)),]) %>% do.call(rbind, .)
  colnames(result) <- colnames(vf$data)
  return(result)
}

sim_vf_single <- function(init, f, length) {
  dim <- length(init)
  result <- matrix(NA_real_, nrow = length, ncol = dim)
  result[1, ] <- init
  for (i in 2:length) {
    prev <- result[i - 1, ]
    dyn <- f(prev)
    next_point <- prev + dyn$v + MASS::mvrnorm(mu = rep(0, dim), Sigma = dyn$a)
    result[i, ] <- next_point
  }
  return(result)
}




sim_vf_options <- function(vf, chains = 10, length = 1e4, discard = 0.3, inits = matrix(c(
	runif(chains, min = vf$x_start, max = vf$x_end),
	runif(chains, min = vf$y_start, max = vf$y_end)
), ncol = 2)) {
	list(vf = vf, chains = chains, length = length, discard = discard, inits = inits)
}

simlandr_options <- function(vf, x = vf$x, y = vf$y, lims = c(vf$x_start, vf$x_end, vf$y_start, vf$y_end), Umax = 8, h = 0.1) {
	list(x = x, y = y, lims = lims, Umax = Umax, h = h)
}
