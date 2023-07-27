#' @inherit fit_2d_nf
#' @export
fit_nd_vf <- function(data, vars,
											lims,
											n = 20,
											vector_position = "start",
											na_action = "omit_data_points",
											method = c("MVKE", "VFC"), ...) {
	d <- data
	# extract useful data for construction
	if (is.data.frame(d)) {
		d_raw <- d[, vars] %>% as.matrix()
	} else if (is.matrix(d)) {
		d_raw <- d[, vars]
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
		`colnames<-`(vars, paste0("v", vars))

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
