#' Estimate a 2D vector field
#'
#' Estimate a 2D vector field from intensive longitudinal data. Two methods can be used: Sparse Vector Field Consensus (SparseVFC, using [SparseVFC::SparseVFC()]), or Multivariate Vector Field Kernel Estimator (MVKE, using [MVKE()]). Note that the input data are automatically normalized before being sent to the estimation engines to make sure the default parameter settings are close to the optimal. Therefore, you do not need to scale up or down the parameters of [SparseVFC::SparseVFC()] or [MVKE()].
#'
#' @param d The data set used for estimating the vector field.
#' Should be a data frame or a matrix.
#' @param x,y Characters to indicate the name of the two variables.
#' @param x_start,x_end,x_by,y_start,y_end,y_by Specify the sample grid for the output: the starting, ending, and increment of the sequences for both axes.
#' @param vector_position Only useful if `method == "VFC"`. One of "start", "middle", or "end", representing the position of the vectors. If "start", for example, the starting point of a vector is regarded as the position of the vector.
#' @param na_action One of "omit_data_points" or "omit_vectors". If using "omit_data_points", then only the `NA` points are omitted, and the points before and after an `NA` will form a vector. If using "omit_vectors", then the vectors will be omitted if either of its points is `NA`.
#' @param method One of "VFC" or "MVKE".
#' @param ... Other parameters to be passed to [SparseVFC::SparseVFC()] or [MVKE()].
#'
#' @return A `vectorfield` object.
#' @seealso [plot.vectorfield()]
#' @export
fit_vf_2d <- function(d, x, y,
											x_start = min(d[,x], na.rm = TRUE) - 0.1*(max(d[,x], na.rm = TRUE)-min(d[,x], na.rm = TRUE)),
											x_end = max(d[,x], na.rm = TRUE) + 0.1*(max(d[,x], na.rm = TRUE)-min(d[,x], na.rm = TRUE)),
											x_by = 0.05*(max(d[,x], na.rm = TRUE)-min(d[,x], na.rm = TRUE)),
											y_start = min(d[,y], na.rm = TRUE) - 0.1*(max(d[,y], na.rm = TRUE)-min(d[,y], na.rm = TRUE)),
											y_end = max(d[,y], na.rm = TRUE) + 0.1*(max(d[,y], na.rm = TRUE)-min(d[,y], na.rm = TRUE)),
											y_by = 0.05*(max(d[,x], na.rm = TRUE)-min(d[,x], na.rm = TRUE)),
											vector_position = "start",
											na_action = "omit_data_points",
											method = c("VFC", "MVKE"), ...) {
  # extract useful data for construction
  if (is.data.frame(d)) {
    d_raw <- d[, c(x, y)] %>% as.matrix()
  } else if (is.matrix(d)) {
    d_raw <- d[, c(x, y)]
  } else {
    rlang::abort("`d` must be a data frame or a matrix.")
  }

	if(na_action != "omit_data_points" & na_action != "omit_vectors") {
		rlang::abort('`na_action` must be either "omit_data_points" or "omit_vectors".')
	}

	dv <- normalize_vecs(d_raw)

	if(any(is.na(dv))){
		{
			if(na_action == "omit_data_points") {
				dv <- attr(dv, "x_noNA")
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

  original_vectors_normalized <- cbind(x_mat, v_mat) %>%
  	`colnames<-`(c("x", "y", "vx", "vy"))

  if(any(is.na(original_vectors_normalized))){
  	{
  		if(na_action == "omit_vectors") {
  			original_vectors_normalized <- stats::na.omit(original_vectors_normalized)
  			rlang::inform("NA(s) found in the data. Those vectors were omitted.")
  		}
  	}
  }

  original_vectors <- original_vectors_normalized
  original_vectors[,"x"] <- original_vectors[,"x"] %>% denormalize_x(dv)
  original_vectors[,"y"] <- original_vectors[,"y"] %>% denormalize_y(dv)
  original_vectors[,"vx"] <- original_vectors[,"vx"] %>% scale_up(dv)
  original_vectors[,"vy"] <- original_vectors[,"vy"] %>% scale_up(dv)

  VFCresult <- MVKEresult <- NULL
  method <- toupper(method[1])
  if(method == "VFC"){
  	VFCresult <- SparseVFC::SparseVFC(original_vectors_normalized[,1:2], original_vectors_normalized[,3:4], ...)
  } else if(method == "MVKE") {
  	MVKEresult <- MVKE(original_vectors_normalized[,1:2], ...)
  }

  vec <- expand.grid(x = seq(x_start, x_end, x_by), y = seq(y_start, y_end, y_by))
  vec <- vec %>% dplyr::rowwise()

  if(method == "VFC") {
  	vec <- vec %>% dplyr::mutate(v = list(stats::predict(VFCresult, c(x, y)%>% normalize_v(dv)) %>% scale_up(dv)))
  } else if(method == "MVKE") {
  	vec <- vec %>% dplyr::mutate(v = list(MVKEresult(c(x, y)%>% normalize_v(dv))$mu%>% scale_up(dv)))
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
  	data = d_raw,
  	data_normalized = dv,
  	original_vectors = original_vectors,
  	original_vectors_normalized = original_vectors_normalized,
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
#' @param x_sparse,y_sparse A number. For calculating the non-gradient part of the vector field, how much should the sample points be sparser than the sample points for the potential landscape?
#' @param ... Other parameters passed to [approxPot2D()].
#'
#' @export
#' @inheritParams fit_vf_2d
#' @inheritParams vectorfield_nongradient_2D
#' @references Rodríguez-Sánchez, P., Nes, E. H. van, & Scheffer, M. (2020). Climbing Escher’s stairs: A way to approximate stability landscapes in multidimensional systems. PLOS Computational Biology, 16(4), e1007788. https://doi.org/10.1371/journal.pcbi.1007788
#' @export
#' @seealso [plot.2d_vf_landscape()]
fit_vfld_2d <- function(vf,
												x_start = vf$x_start, x_end = vf$x_end, x_by = vf$x_by/2,
												y_start = vf$y_start, y_end = vf$y_end, y_by = vf$y_by/2, x_sparse = 2, y_sparse = 2, ...){
	if(!inherits(vf, "vectorfield")) rlang::abort("`vf` must be a `vectorfield` object.")
	xs <- seq(x_start, x_end, x_by)
	ys <- seq(y_start, y_end, y_by)

	if(vf$method == "VFC") f <- function(x) stats::predict(vf$VFCresult, x %>% normalize_v(vf$data_normalized)) %>% scale_up(vf$data_normalized)
	else if(vf$method == "MVKE") f <- function(x) vf$MVKEresult(x %>% normalize_v(vf$data_normalized))$mu %>% scale_up(vf$data_normalized)
	wdresult <- approxPot2D(f, xs = xs, ys = ys, ...)
	wdresult$xs <- xs
	wdresult$ys <- ys
	dist <- simlandr::make_2d_tidy_dist(list(x = xs, y = ys, z = wdresult$V))
	dist_error <- simlandr::make_2d_tidy_dist(list(x = xs, y = ys, z = wdresult$err))

	dist_nongrad <- vectorfield_nongradient_2D(wdresult, x_sparse, y_sparse)

	result <- list(
		dist = dist,
		dist_error = dist_error,
		dist_raw = wdresult,
		dist_nongrad = dist_nongrad,
		vf = vf
	)

	class(result) <- c("2d_vf_landscape", "landscape")

	return(result)
}
