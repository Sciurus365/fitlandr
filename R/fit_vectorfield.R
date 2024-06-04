#' Estimate a 2D vector field
#'
#' Estimate a 2D vector field from intensive longitudinal data. Two methods can be used: Multivariate Vector Field Kernel Estimator (MVKE, using [MVKE()]), or Sparse Vector Field Consensus (SparseVFC, using [SparseVFC::SparseVFC()]). Note that the input data are automatically normalized before being sent to the estimation engines to make sure the default parameter settings are close to the optimal. Therefore, you do not need to scale up or down the parameters of [MVKE()] or [SparseVFC::SparseVFC()]. We suggest the MVKE method to be used for psychological data because it has more realistic assumptions and produces more reasonable output.
#'
#' @param data The data set used for estimating the vector field.
#' Should be a data frame or a matrix.
#' @param x,y Characters to indicate the name of the two variables.
#' @param lims The limits of the range for the vector field estimation as `c(<xl>, <xu>, <yl>, <yu>)`. If missing, the range of the data extended by 10% for both sides will be used.
#' @param n The number of equally spaced points in each axis, at which the vectors are to be estimated.
#' @param vector_position Only useful if `method == "VFC"`. One of "start", "middle", or "end", representing the position of the vectors. If "start", for example, the starting point of a vector is regarded as the position of the vector.
#' @param na_action One of "omit_data_points" or "omit_vectors". If using "omit_data_points", then only the `NA` points are omitted, and the points before and after an `NA` will form a vector. If using "omit_vectors", then the vectors will be omitted if either of its points is `NA`.
#' @param method One of "MVKE" or "VFC".
#' @param ... Other parameters to be passed to [MVKE()] or [SparseVFC::SparseVFC()].
#'
#' @return A `vectorfield` object.
#' @seealso [plot.vectorfield()]
#'
#' @examples
#' # generate data
#' single_output_grad <- simlandr::sim_fun_grad(length = 200, seed = 1614)
#' # fit the vector field
#' v2 <- fit_2d_vf(single_output_grad, x = "x", y = "y", method = "MVKE")
#' plot(v2)
#' @export
fit_2d_vf <- function(data, x, y,
                      lims,
                      n = 20,
                      vector_position = "start",
                      na_action = "omit_data_points",
                      method = c("MVKE", "VFC"), ...) {
  d <- data
  # extract useful data for construction
  if (is.data.frame(d)) {
    d_raw <- d[, c(x, y)] %>% as.matrix()
  } else if (is.matrix(d)) {
    d_raw <- d[, c(x, y)]
  } else {
    rlang::abort("`d` must be a data frame or a matrix.")
  }

  if (na_action != "omit_data_points" & na_action != "omit_vectors") {
    rlang::abort('`na_action` must be either "omit_data_points" or "omit_vectors".')
  }

  dv <- normalize_vecs(d_raw)

  if (any(is.na(dv))) {{ if (na_action == "omit_data_points") {
    dv <- attr(dv, "x_noNA")
    rlang::inform("NA(s) found in the data. Those data points were omitted.")
  } }}

  v_mat <- diff(dv)

  if (vector_position == "start") {
    x_mat <- dv[1:(nrow(dv) - 1), ]
  } else if (vector_position == "middle") {
    x_mat <- dv[1:(nrow(dv) - 1), ] + 0.5 * v_mat
  } else if (vector_position == "end") {
    x_mat <- dv[2:nrow(dv), ]
  } else {
    rlang::abort('`vector_position` must be one of "start", "middle", or "end".')
  }

  original_vectors_normalized <- cbind(x_mat, v_mat) %>%
    `colnames<-`(c("x", "y", "vx", "vy"))

  if (any(is.na(original_vectors_normalized))) {{ if (na_action == "omit_vectors") {
    original_vectors_normalized <- stats::na.omit(original_vectors_normalized)
    rlang::inform("NA(s) found in the data. Those vectors were omitted.")
  } }}

  original_vectors <- original_vectors_normalized
  original_vectors[, "x"] <- original_vectors[, "x"] %>% denormalize_x(dv)
  original_vectors[, "y"] <- original_vectors[, "y"] %>% denormalize_y(dv)
  original_vectors[, "vx"] <- original_vectors[, "vx"] %>% scale_up(dv)
  original_vectors[, "vy"] <- original_vectors[, "vy"] %>% scale_up(dv)

  VFCresult <- MVKEresult <- NULL
  method <- toupper(method[1])
  if (method == "VFC") {
    VFCresult <- SparseVFC::SparseVFC(original_vectors_normalized[, 1:2], original_vectors_normalized[, 3:4], ...)
  } else if (method == "MVKE") {
    MVKEresult <- MVKE(original_vectors_normalized[, 1:2], ...)
  }

  lims <- determine_lims(d_raw, c(x, y), lims)

  vec <- tidyr::expand_grid(x = seq(lims[1], lims[2], length.out = n), y = seq(lims[3], lims[4], length.out = n))
  vec <- vec %>% dplyr::rowwise()

  if (method == "VFC") {
    vec <- vec %>% dplyr::mutate(v = list(stats::predict(VFCresult, c(x, y) %>% normalize_v(dv)) %>% scale_up(dv)))
  } else if (method == "MVKE") {
    vec <- vec %>% dplyr::mutate(v = list(MVKEresult(c(x, y) %>% normalize_v(dv))$mu %>% scale_up(dv)))
  }

  vec <- vec %>%
    dplyr::mutate(
      vx = v[1],
      vy = v[2],
      v_norm = (sum(v^2))^(1 / 2)
    ) %>%
    dplyr::select(-v) %>%
    dplyr::ungroup()

  result <- list(
    vec_grid = vec,
    VFCresult = VFCresult,
    MVKEresult = MVKEresult,
    data = d_raw,
    data_normalized = dv,
    original_vectors = original_vectors,
    original_vectors_normalized = original_vectors_normalized,
    x = x,
    y = y,
    lims = lims,
    n = n,
    method = method
  )

  class(result) <- "vectorfield"

  return(result)
}

# borrowed from simlandr 0.3.0
determine_lims <- function(output, var_names, lims) {
  if (!rlang::is_missing(lims)) {
    return(lims)
  }
  # if (is.list(output)) output <- output[[1]]
  if (rlang::is_missing(lims)) {
    return(c(sapply(var_names, function(v) grDevices::extendrange(output[, v], f = 0.1))))
  }
  if (any(is.infinite(lims))) stop("Non-infinite values found in `lims`.")
}

#' Calculate the vector value at a given position
#'
#' @param object A `vectorfield` project generated by [fit_2d_vf()].
#' @param pos A vector, the position of the vector.
#' @param linear_interp Use linear interpolation method to estimate the drift vector (and the diffusion matrix). This can speed up the calculation. If `TRUE`, be sure that a linear grid was calculated for the vector field using `<vf> <- add_interp_grid(<vf>)`.
#' @param calculate_a Effective when `linear_interp == TRUE`. Do you want to calculate the diffusion matrix? Use `FALSE` can save some time.
#' @param ... Not in use.
#'
#' @return A list of `v`, the drift part that is used for vector fields, and `a` (when `calculate_a == TRUE`), the diffusion part at a given position.
#'
#' @seealso [add_interp_grid()]
#' @export
predict.vectorfield <- function(object, pos, linear_interp = FALSE, calculate_a = TRUE, ...) {
	if(!linear_interp) {
		return(
			normalize_predict_f(object)(pos)
		)
	} else {
		if(is.null(object$interp_grid)) {
			stop("The grid for interpolation does not exist. Either use `<vf> <- add_interp_grid(<vf>)` to generate it, or use `linear_vec = TRUE`.")
		}
		result <- list()
		temp <- object$interp_grid
		result$v <- c(
			fast_bilinear(x = temp$x, y = temp$y, z = temp$z_vector_1, x0 = pos[1], y0 = pos[2])$z,
			fast_bilinear(x = temp$x, y = temp$y, z = temp$z_vector_2, x0 = pos[1], y0 = pos[2])$z
			)
		if (calculate_a) {
			result$a <- matrix(
				c(
					fast_bilinear(x = temp$x, y = temp$y, z = temp$z_diffusion_1, x0 = pos[1], y0 = pos[2])$z,
					fast_bilinear(x = temp$x, y = temp$y, z = temp$z_diffusion_2, x0 = pos[1], y0 = pos[2])$z,
					fast_bilinear(x = temp$x, y = temp$y, z = temp$z_diffusion_3, x0 = pos[1], y0 = pos[2])$z,
					fast_bilinear(x = temp$x, y = temp$y, z = temp$z_diffusion_4, x0 = pos[1], y0 = pos[2])$z
				), byrow = TRUE, nrow = 2, ncol = 2
			)
		}
		return(result)
	}
}


#' Add a grid to a `vectorfield` object to enable linear interpolation
#'
#' @inheritParams fit_3d_vfld
#' @inheritParams fit_2d_vf
#'
#' @return A `vectorfield` project with an `interp_grid` field.
#' @export
add_interp_grid <- function(vf, lims = vf$lims, n = vf$n) {
	if(!is.null(vf$interp_grid)) {
		message("There is already an `interp_grid` field in the `vectorfield`. It will be overwritten.")
	}
	x <- seq(lims[1], lims[2], length.out = n)
	y <- seq(lims[3], lims[4], length.out = n)
	z_vector_1 <- outer(x, y, Vectorize(function(x, y) normalize_predict_f(vf)(c(x, y))$v[1]))
	z_vector_2 <- outer(x, y, Vectorize(function(x, y) normalize_predict_f(vf)(c(x, y))$v[2]))
	z_diffusion_1 <- outer(x, y, Vectorize(function(x, y) normalize_predict_f(vf)(c(x, y))$a[1]))
	z_diffusion_2 <- outer(x, y, Vectorize(function(x, y) normalize_predict_f(vf)(c(x, y))$a[2]))
	z_diffusion_3 <- outer(x, y, Vectorize(function(x, y) normalize_predict_f(vf)(c(x, y))$a[3]))
	z_diffusion_4 <- outer(x, y, Vectorize(function(x, y) normalize_predict_f(vf)(c(x, y))$a[4]))
	vf$interp_grid <- list(x = x, y = y, z_vector_1 = z_vector_1, z_vector_2 = z_vector_2,
												 z_diffusion_1 = z_diffusion_1,
												 z_diffusion_2 = z_diffusion_2,
												 z_diffusion_3 = z_diffusion_3,
												 z_diffusion_4 = z_diffusion_4)
	return(vf)
}

#' Return a normalized prediction function
#'
#' @inheritParams fit_3d_vfld
#' @return A function that takes a vector `x` and returns a list of `v`, the drift part, and `a`, the diffusion part.
normalize_predict_f <- function(vf) {
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
	return(f)
}

#' A fast bilinear interpolation function
#'
#' It assumes equal grid intervals, thus can find the correct position in \eqn{O(1)} time.
#'
#' The following is from the documentation of [akima::bilinear()].
#'
#' This is an implementation of a bilinear interpolating function.
#' For a point (x0,y0) contained in a rectangle (x1,y1),(x2,y1), (x2,y2),(x1,y2) and x1<x2, y1<y2, the first step is to get z() at locations (x0,y1) and (x0,y2) as convex linear combinations z(x0,y*)=a*z(x1,y*)+(1-a)*z(x2,y*) where a=(x2-x1)/(x0-x1) for y*=y1,y2. In a second step z(x0,y0) is calculated as convex linear combination between z(x0,y1) and z(x0,y2) as z(x0,y1)=b*z(x0,y1)+(1-b)*z(x0,y2) where b=(y2-y1)/(y0-y1).
#'
#' Finally, z(x0,y0) is a convex linear combination of the z values at the corners of the containing rectangle with weights according to the distance from (x0,y0) to these corners.
#'
#' @inheritParams akima::bilinear
#' @param x0 `x` value used to interpolate at.
#' @param y0 `y` value used to interpolate at.
#'
#' @return A list which contains only one element, `z`.
#' @export
#' @keywords internal
#'
fast_bilinear <- function(x, y, z, x0, y0) {
  x_ind <- floor((x0 - x[1]) / (x[2] - x[1])) + 1
  y_ind <- floor((y0 - y[1]) / (y[2] - y[1])) + 1
  if (x_ind < 1 | x_ind > length(x) - 1 | y_ind < 1 | y_ind > length(y) - 1) {
    return(list(z = 0))
  }

  z11 <- z[x_ind, y_ind]
  z12 <- z[x_ind, y_ind + 1]
  z21 <- z[x_ind + 1, y_ind]
  z22 <- z[x_ind + 1, y_ind + 1]
  ex <- (x0 - x[x_ind]) / (x[x_ind + 1] - x[x_ind])
  ey <- (y0 - y[y_ind]) / (y[y_ind + 1] - y[y_ind])
  result <- (1 - ex) * (1 - ey) * z11 + (1 - ex) * (ey) * z12 + (ex) * (1 - ey) * z21 + (ex) * (ey) * z22

  return(list(z = result))
}
