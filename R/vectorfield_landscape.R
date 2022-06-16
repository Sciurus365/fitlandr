#' Estimate a 2D vector field
#'
#' @param d The data set used for estimating the vector field.
#' Should be a data frame or a matrix.
#' @param x,y Characters to indicate the name of the two variables.
#' @param vector_position One of `"start"`, `"middle"`, or `"end"`.
#' @param na_action One of `"omit_data_points"` or `"omit_vectors"`. If using `"omit_data_points"`, then the vectors will be retained even if there are `NA` points between them. If using `"omit_vectors"`, then the vectors will be omitted if either of its points is `NA`.
#' @param ... Other parameters to be passed to [SparseVFC::SparseVFC()].
#'
#' @return A `vectorfield` object.
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

  VFCresult <- mvkeresult <- NULL
  if(method == "VFC"){
  	VFCresult <- SparseVFC::SparseVFC(original_vectors[,1:2], original_vectors[,3:4], ...)
  } else if(method == "mvke") {
  	mvkeresult <- mvke(original_vectors[,1:2], ...)
  }

  vec <- expand.grid(x = seq(x_start, x_end, x_by), y = seq(y_start, y_end, y_by))
  vec <- vec %>% dplyr::rowwise()

  if(method == "VFC") {
  	vec <- vec %>% dplyr::mutate(v = list(predict(VFCresult, c(x, y))))
  } else if(method == "mvke") {
  	vec <- vec %>% dplyr::mutate(v = list(mvkeresult(c(x, y))$mu))
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
  	mvkeresult = mvkeresult,
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

#' @export
fit_vfld_2d <- function(...){
	UseMethod("fit_vfld_2d")
}

#' @export
fit_vfld_2d.data.frame <- fit_vfld_2d.matrix <- function(...){
	vf <- fit_vf_2d(...)
	fit_vfld_2d.vectorfield(vf)
}

#' @export
#' @inheritParams fit_vf_2d
fit_vfld_2d.vectorfield <- function(vf,
												x_start = vf$x_start, x_end = vf$x_end, x_by = vf$x_by,
												y_start = vf$y_start, y_end = vf$y_end, y_by = vf$y_by){
	if(!inherits(vf, "vectorfield")) rlang::abort("`vf` must be a `vectorfield` object.")
	xs <- seq(x_start, x_end, x_by)
	ys <- seq(y_start, y_end, y_by)

	if(vf$method == "VFC") f <- function(x) predict(vf$VFCresult, x)
	else if(vf$method == "mvke") f <- function(x) vf$mvkeresult(x)$mu
	wdresult <- waydown::approxPot2D(f, xs = xs, ys = ys)
	wdresult$xs <- xs
	wdresult$ys <- ys
	dist <- simlandr::make_2d_tidy_dist(list(x = xs, y = ys, z = wdresult$V))

	result <- list(
		dist = dist,
		dist_raw = wdresult,
		vf = vf
	)

	class(result) <- c("2d_vf_landscape", "landscape")

	return(result)
}
