#' Estimate a 2D potential landscape from data with the MVKE method
#'
#' This function is a wrapper of the MVKE method (see [MVKE()]) that produces a 2D potential landscape from 1D data. The landscape is constructed by estimating the gradient of the data and then integrating it. The MVKE method is a non-parametric method that estimates the gradient of the data by using a kernel density estimator. The potential landscape is then constructed by integrating the gradient.
#'
#' @param data A data frame or matrix containing the data. The data frame should contain at least a column, with the column name indicated by `x`, that represents the dimension for landscape construction.
#' @param x The column name of the data frame that represents the dimension for landscape construction.
#' @param lims The limits of the range for the landscape calculation as `c(xl, xu)`.
#' @param n The number of equally spaced points in the axis, at which the landscape is to be estimated.
#' @param method The method used to estimate the gradient. Currently only "MVKE" is supported.
#' @param ... Additional arguments passed to [MVKE()].
#' @inheritParams stats::integrate
#' @inheritParams fit_2d_vf
#' @return A `2d_MVKE_landscape` object, which contains the following components:
#' \itemize{
#'   \item `dist`: A data frame containing the estimated potential landscape. The data frame has two columns: `x` and `U`, where `x` is the position and `U` is the potential.
#'   \item `p`: A ggplot object containing the plot of the potential landscape.
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
#' plot(l)
#'
#' # different behaviors for different `na_action` choices
#'
#' l1 <- fit_2d_ld(data.frame(x = c(1,2,1,2,NA,NA,NA,10,11,10,11)), "x")
#' plot(l1)
#'
#' l2 <- fit_2d_ld(data.frame(x = c(1,2,1,2,NA,NA,NA,10,11,10,11)), "x", na_action = "omit_vectors")
#' plot(l2)
#'
#'
fit_2d_ld <- function(data, x, lims, n = 200L, vector_position = "start", na_action = "omit_data_points",
											method = c("MVKE"), subdivisions = 100L, rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol, stop.on.error = TRUE, keep.xy = FALSE, aux = NULL, ...) {
	d <- data
	# extract useful data for construction
	if (is.data.frame(d)) {
		d_raw <- d[, c(x), drop = FALSE] %>% as.matrix()
	} else if (is.matrix(d)) {
		d_raw <- d[, c(x), drop = FALSE]
	} else {
		rlang::abort("`d` must be a data frame or a matrix.")
	}

	if (na_action != "omit_data_points" & na_action != "omit_vectors") {
		rlang::abort('`na_action` must be either "omit_data_points" or "omit_vectors".')
	}

	if (na_action == "omit_data_points" & any(is.na(d_raw))) {
		d_raw <- stats::na.omit(d_raw)
		rlang::inform("NA(s) found in the data. Those data points were omitted.")
	}

	if (vector_position == "start") {
		x_mat <- d_raw[1:(nrow(d_raw) - 1), , drop = FALSE]
	} else if (vector_position == "middle") {
		x_mat <- d_raw[1:(nrow(d_raw) - 1), , drop = FALSE] + 0.5 * v_mat
	} else if (vector_position == "end") {
		x_mat <- d_raw[2:nrow(d_raw), , drop = FALSE]
	} else {
		rlang::abort('`vector_position` must be one of "start", "middle", or "end".')
	}

	v_mat <- diff(d_raw)

	data_vectors <- cbind(x_mat, v_mat) %>%
		`colnames<-`(c("x", "vx"))

	if (any(is.na(data_vectors))) {{ if (na_action == "omit_vectors") {
		data_vectors <- stats::na.omit(data_vectors)
		rlang::inform("NA(s) found in the data. Those vectors were omitted.")
	} }}

	lims <- determine_lims(data, x, lims)
  MVKEresult <- MVKE(data_vectors[,1, drop = FALSE], data_vectors[,2, drop = FALSE], ...)

  xseq <- seq(lims[1], lims[2], length.out = n)
  Useq <- vector("numeric", length = n)

  Useq[1] <- 0
  for (i in 2:n) {
    Useq[i] <- Useq[i - 1] - stats::integrate(function(x) purrr::map_dbl(x, function(xx) MVKEresult(xx)$mu), xseq[i - 1], xseq[i], subdivisions = subdivisions, rel.tol = rel.tol, abs.tol = abs.tol, stop.on.error = stop.on.error, keep.xy = keep.xy, aux = aux, ...)$value
  }

  dist <- data.frame(x = xseq, U = Useq)
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = dist, ggplot2::aes(x = x, y = U)) +
    ggplot2::theme_bw()


  return(structure(list(dist = dist, plot = p, MVKEresult = MVKEresult), class = c("2d_MVKE_landscape", "landscape")))
}

#' @export
#' @describeIn fit_2d_ld Find the local minima of the 2D potential landscape
#' @param object An object of class `2d_MVKE_landscape` returned by [fit_2d_ld()].
#' @param ... Not used.
summary.2d_MVKE_landscape <- function(object, ...) {
  # find the local minimum values in object$dist$U
	# return a data frame with the x and U values of the local minima

	local_minima <- which(diff(sign(diff(object$dist$U))) == 2) + 1
	cli::cli_inform("{length(local_minima)} local minima were found.")
  return(data.frame(x = object$dist$x[local_minima], U = object$dist$U[local_minima]))
}

# This function is taken from the `simlandr` package.
determine_lims <- function (output, var_names, lims)
{
	if (!rlang::is_missing(lims)) {
		return(lims)
	}
	if (is.list(output))
		output <- output[[1]]
	if (rlang::is_missing(lims)) {
		return(c(sapply(var_names, function(v) grDevices::extendrange(output[,
																																				 v], f = 0.1))))
	}
	if (any(is.infinite(lims)))
		stop("Non-infinite values found in `lims`.")
}

utils::globalVariables("U")
