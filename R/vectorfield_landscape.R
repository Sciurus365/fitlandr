#' Estimate a 2D vector field
#'
#' Estimate a 2D vector field from intensive longitudinal data. Two methods can be used: Sparse Vector Field Consensus (SparseVFC, using [SparseVFC::SparseVFC()]), or Multivariate Vector Field Kernel Estimator (MVKE, using [MVKE()]).
#'
#' @param d The data set used for estimating the vector field.
#' Should be a data frame or a matrix.
#' @param x,y Characters to indicate the name of the two variables.
#' @param x_start,x_end,x_by,y_start,y_end,y_by Specify the sample grid for the output: the starting, ending, and increment of the sequences for both axes.
#' @param vector_position One of "start", "middle", or "end". If "start",
#' @param na_action One of "omit_data_points" or "omit_vectors". If using "omit_data_points", then the vectors will be retained even if there are `NA` points between them. If using "omit_vectors", then the vectors will be omitted if either of its points is `NA`.
#' @param method One of "VFC" or "MVKE".
#' @param ... Other parameters to be passed to [SparseVFC::SparseVFC()] or [MVKE()].
#'
#' @return A `vectorfield` object.
#' @seealso [plot.vectorfield()]
#' @export
fit_vf_2d <- function(d, x, y,
											x_start, x_end, x_by, y_start, y_end, y_by,
											vector_position = "middle", na_action = "omit_data_points",
											method = "VFC", ...) {
  # extract useful data for construction
  if (is.data.frame(d)) {
    dv <- d[, c(x, y)] %>% as.matrix()
  } else if (is.matrix(d)) {
    dv <- d[, c(x, y)]
  } else {
    rlang::abort("`d` must be a data frame or a matrix.")
  }

	if(na_action != "omit_data_points" & na_action != "omit_vectors") {
		rlang::abort('`na_action` must be either "omit_data_points" or "omit_vectors".')
	}

	if(any(is.na(dv))){
		{
			if(na_action == "omit_data_points") {
				dv <- na.omit(dv)
				rlang::inform("NA(s) found in the data. Those data points were omitted.")
			}
		}
	}

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

  original_vectors <- cbind(x_mat, v_mat) %>%
  	`colnames<-`(c("x", "y", "vx", "vy"))

  if(any(is.na(original_vectors))){
  	{
  		if(na_action == "omit_vectors") {
  			original_vectors <- na.omit(original_vectors)
  			rlang::inform("NA(s) found in the data. Those vectors were omitted.")
  		}
  	}
  }

  VFCresult <- MVKEresult <- NULL
  method <- toupper(method)
  if(method == "VFC"){
  	VFCresult <- SparseVFC::SparseVFC(original_vectors[,1:2], original_vectors[,3:4], ...)
  } else if(method == "MVKE") {
  	MVKEresult <- MVKE(original_vectors[,1:2], ...)
  }

  vec <- expand.grid(x = seq(x_start, x_end, x_by), y = seq(y_start, y_end, y_by))
  vec <- vec %>% dplyr::rowwise()

  if(method == "VFC") {
  	vec <- vec %>% dplyr::mutate(v = list(predict(VFCresult, c(x, y))))
  } else if(method == "MVKE") {
  	vec <- vec %>% dplyr::mutate(v = list(MVKEresult(c(x, y))$mu))
  }

  vec <- vec %>%
  	dplyr::mutate(
  		vx = v[1],
  		vy = v[2]
  	) %>%
  	dplyr::select(-v) %>%
  	dplyr::ungroup()

  result <- list(
  	vec_grid = vec,
  	VFCresult = VFCresult,
  	MVKEresult = MVKEresult,
  	data = dv,
  	original_vectors = original_vectors,
  	x_start = x_start,
  	x = x,
  	y = y,
  	x_end = x_end,
  	x_by = x_by,
  	y_start = y_start,
  	y_end = y_end,
  	y_by = y_by,
  	method = method
  )

  class(result) <- "vectorfield"

	return(result)
}

#' Estimate a 2D potential landscape from a vector field
#'
#' This is done by a decomposition method described in the reference, using a modified version of [waydown::approxPot2D()].
#'
#' @param vf A `vectorfield` object estimated by [fit_vf_2d()]. The other parameters will be inherited from the `vectorfield` object if not otherwise specified.
#' @param ... Other parameters passed to [approxPot2D()].
#'
#' @export
#' @inheritParams fit_vf_2d
#' @references Rodríguez-Sánchez, P., Nes, E. H. van, & Scheffer, M. (2020). Climbing Escher’s stairs: A way to approximate stability landscapes in multidimensional systems. PLOS Computational Biology, 16(4), e1007788. https://doi.org/10.1371/journal.pcbi.1007788
#' @export
#' @seealso [plot.2d_vf_landscape()]
fit_vfld_2d <- function(vf,
												x_start = vf$x_start, x_end = vf$x_end, x_by = vf$x_by,
												y_start = vf$y_start, y_end = vf$y_end, y_by = vf$y_by, ...){
	if(!inherits(vf, "vectorfield")) rlang::abort("`vf` must be a `vectorfield` object.")
	xs <- seq(x_start, x_end, x_by)
	ys <- seq(y_start, y_end, y_by)

	if(vf$method == "VFC") f <- function(x) predict(vf$VFCresult, x)
	else if(vf$method == "MVKE") f <- function(x) vf$MVKEresult(x)$mu
	wdresult <- approxPot2D(f, xs = xs, ys = ys, ...)
	wdresult$xs <- xs
	wdresult$ys <- ys
	dist <- simlandr::make_2d_tidy_dist(list(x = xs, y = ys, z = wdresult$V))
	dist_error <- simlandr::make_2d_tidy_dist(list(x = xs, y = ys, z = wdresult$err))

	result <- list(
		dist = dist,
		dist_error = dist_error,
		dist_raw = wdresult,
		vf = vf
	)

	class(result) <- c("2d_vf_landscape", "landscape")

	return(result)
}
