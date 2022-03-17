#' Estimate a 2D vector field
#'
#' @param d The data set used for estimating the vector field.
#' Should be a data frame or a matrix.
#' @param x,y Characters to indicate the name of the two variables.
#' @param vector_position One of "start", "middle", or "end".
#' @param ... Other parameters to be passed to [SparseVFC].
#'
#' @return A `vectorfield` object.
#' @export
fit_vf_2d <- function(d, x, y, vector_position = "middle",
											x_start, x_end, x_by, y_start, y_end, y_by,
											...) {
  # extract useful data for construction
  if (is.data.frame(d)) {
    dv <- d[, c(x, y)] %>% as.matrix()
  } else if (is.matrix(d)) {
    dv <- d[, c(x, y)]
  } else {
    rlang::abort("`d` must be a data frame or a matrix.")
  }

	if(any(is.na(dv))){
		dv <- na.omit(dv)
		rlang::inform("NA(s) found in the data. Those cases were omitted.")
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

  VFCresult <- SparseVFC(x_mat, v_mat, ...)

  vec <- expand.grid(x = seq(x_start, x_end, x_by), y = seq(y_start, y_end, y_by))
  vec <- vec %>%
  	dplyr::mutate(v = purrr::map2(x, y, function(x, y){
  		output <- c(0,0)
  		for(i in 1:nrow(VFCresult$X)){
  			output <- output + con_K(c(x, y) %>% matrix(nrow = 1), VFCresult$X[i,] %>% matrix(nrow = 1), VFCresult$beta) %>%
  				as.numeric() * VFCresult$C[i,]
  		}
  		return(list(output))
  	}))

  vec <- vec %>%
  	dplyr::rowwise() %>%
  	dplyr::mutate(vx = v[1],
  				 vy = vx[2],
  				 vx = vx[1]) %>%
  	dplyr::select(-v)

  result <- list(
  	vec_grid = vec,
  	VFCresult = VFCresult,
  	x_start = x_start,
  	x = x,
  	y = y,
  	x_end = x_end,
  	x_by = x_by,
  	y_start = y_start,
  	y_end = y_end,
  	y_by = y_by
  )

  class(result) <- "vectorfield"

	return(result)
}

fit_vfld_2d <- function(...){
	UseMethod("fit_vfld_2d")
}

fit_vfld_2d.data.frame <- fit_vfld_2d.matrix <- function(...){
	vf <- fit_vf_2d(...)
	fit_vfld_2d.vectorfield(vf)
}

fit_vfld_2d.vectorfield <- function(vf,
												x_start = vf$x_start, x_end = vf$x_end, x_by = vf$x_by,
												y_start = vf$y_start, y_end = vf$y_end, y_by = vf$y_by){
	if(!inherits(vf, "vectorfield")) rlang::abort("`vf` must be a `vectorfield` object.")
	xs <- seq(x_start, x_end, x_by)
	ys <- seq(y_start, y_end, y_by)

	wdresult <- waydown::approxPot2D(VFCf(vf$VFCresult), xs = xs, ys = ys)
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

#' Generate a function to calculate the estimated vector field
#'
#' This will generate a helper function to estimate the vector field
#' for any position, based on the VFC estimation.
VFCf <- function(VFCresult){
	if(!inherits(VFCresult, "VFC")) rlang::abort("`VFCresult` must be a VFC object.")
	return(function(x){
		output <- c(0,0)
		for(i in 1:nrow(VFCresult$X)){
			output <- output + con_K(c(x[1], x[2]) %>% matrix(nrow = 1), VFCresult$X[i,] %>% matrix(nrow = 1), VFCresult$beta) %>%
				as.numeric() * VFCresult$C[i,]
		}
		return(output)
	})
}
