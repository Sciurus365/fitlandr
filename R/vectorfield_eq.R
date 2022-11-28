#' Find equilibrium points for a vector field
#'
#' @inheritParams sim_vf
#' @param starts A vector indicating the starting value for solving the equilibrium point, or a list of vectors providing multiple starting values together.
#' @param jacobian_params Parameters passed to [numDeriv::jacobian()].
#' @param ... Parameters passed to [rootSolve::multiroot()].
#'
#' @return A list of equilibrium points and their details. Use [print.vectorfield_eqs()] to inspect it.
#' @export
find_eqs <- function(vf, starts, jacobian_params = list(), ...) {
  if (!is.list(starts)) {
    starts <- list(starts)
  }
  root_results <- lapply(starts, function(a) {
    if (vf$method == "VFC") {
      temp <- rootSolve::multiroot(function(x) predict(vf$VFCresult, normalize_v(x, ref = vf$data_normalized)), start = a, ...)
      temp$jacobian <- rlang::exec(numDeriv::jacobian, func = function(x) predict(vf$VFCresult, normalize_v(x, ref = vf$data_normalized)), x = temp$root, !!!jacobian_params)
    } else if (vf$method == "MVKE") {
      temp <- rootSolve::multiroot(function(x) vf$MVKEresult$mu(normalize_v(x, ref = vf$data_normalized)), start = a, ...)
      temp$jacobian <- rlang::exec(numDeriv::jacobian, function(x) vf$MVKEresult$mu(normalize_v(x, ref = vf$data_normalized)), x = temp$root, !!!jacobian_params)
    }

    if (temp$root[1] < vf$lims[1] | temp$root[1] > vf$lims[2] | temp$root[2] < vf$lims[3] | temp$root[2] > vf$lims[4]) {
      temp$check <- FALSE
    } else {
      temp$check <- TRUE
    }

  	temp$dorm_eigen_value <- sort(Re(eigen(temp$jacobian)$values) ,decreasing = TRUE)[1]
  	if(temp$dorm_eigen_value < 0) {
  		temp$stability <- "stable"
  	} else if(temp$dorm_eigen_value > 0) {
  		temp$stability <- "unstable"
  	} else {
  		temp$stability <- "unknown"
  	}
    temp
  })

  return(structure(root_results, class = "vectorfield_eqs"))
}

#' @export
print.vectorfield_eqs <- function(x, ...) {
  cat("Root(s):")
  lapply(x, function(x) cat("\n", x$root, "\t",
  													ifelse(x$stability == "stable", x$stability, cli::col_red(x$stability)), "\t",
  													ifelse(x$check, "[within range]", cli::col_red("[out of range]"))))
}
