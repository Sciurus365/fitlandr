#' Performs bootstrap resampling to estimate uncertainty in 2D vector field estimation.
#'
#' Moving block bootstrap (MBB) is used to account for temporal dependencies in the data.
#'
#' @param vf A `vectorfield` or `cv_vectorfield` object representing the fitted vector field.
#' @param block_length Length of each block for MBB. If NULL, it will be set to n^(1/3).
#' @param n_boot Number of bootstrap samples to generate (default 200).
#' @param seed Random seed for reproducibility (default 1614).
#' @param add_linear_interp Logical indicating whether to add linear interpolation in predictions.
#' @param ... Additional arguments passed to the vector field fitting function.
#' @return An object of class `bootstrap_2d_vf` containing:
#'         - `bootstrap_models`: A list of fitted vector field models from each bootstrap sample.
#'         - `original_vf`: The original fitted vector field model.
#'         - `block_length`: The block length used for MBB.
#'         - `n_boot`: The number of bootstrap samples.
#' @export
bootstrap_2d_vf <- function(vf, block_length = NULL, n_boot = 200, seed = 1614, add_linear_interp = TRUE, ...) {
	# Check the class of vf. If vf is from cv_vectorfield, extract the final_model

	if (inherits(vf, "cv_vectorfield")) {
		vf <- vf$final_model
	}

	if (!inherits(vf, "vectorfield")) {
		cli::cli_abort("Input 'vf' must be a 'vectorfield' or 'cv_vectorfield' object.")
	}

	# Extract the vectors (instead of the data points) as the basic unit for resampling.
	## This data frame is a part of the standard vectorfield object.
	## It contains 4 columns: x, y, vx, vy
	original_vectors <- vf$original_vectors
	original_vectors_normalized <- vf$original_vectors_normalized
	n_vec <- nrow(original_vectors_normalized)

	# Determine block length for MBB
	if(is.null(block_length)) {
		block_length <- ceiling(n_vec^(1/3))
		cli::cli_inform("Block length not provided. Using default block length = {block_length}.")
	}

	# Make blocks (moving blocks of consecutive indices)

	blocks <- lapply(1:(n_vec - block_length + 1), function(start_idx) {
		return(start_idx:(start_idx + block_length - 1))
	})

	n_blocks <- length(blocks)
	n_blocks_needed <- ceiling(n_vec / block_length)

	# Prepare for bootstrap

	bootstrap_models <- vector("list", n_boot)

	# retrieve the parameters from vf
	h <- environment(vf[["MVKEresult"]])[["h"]]
	kernel <- environment(vf[["MVKEresult"]])[["kernel"]]
	dv <- vf$data_normalized
	lims <- vf$lims
	vec <- vf$vec[,1:2] %>% dplyr::rowwise()
	x <- vf$x
	y <- vf$y
	n <- vf$n
	method <- vf$method
	d_raw <- vf$data
	normalize_v <- fitlandr:::normalize_v
	scale_up <- fitlandr:::scale_up

	p <- progressr::progressor(steps = n_boot)
	## Use future here for parallel processing. User can set up plan() outside this function.
	## If not set up, it will run sequentially.
	bootstrap_models <- furrr::future_map(
		1:n_boot,
		function(b) {
			# Sample blocks with replacement
			sampled_block_indices <- sample(1:n_blocks, n_blocks_needed, replace = TRUE)
			sampled_indices <- unlist(blocks[sampled_block_indices])
			sampled_indices <- sampled_indices[sampled_indices <= n_vec] # Trim to original size
			sampled_vectors <- original_vectors[sampled_indices, ]
			sampled_vectors_normalized <- original_vectors_normalized[sampled_indices, ]
			# Fit vector field to the sampled vectors
			MVKEresult <- fitlandr::MVKE(d = sampled_vectors_normalized[, 1:2], v = sampled_vectors_normalized[, 3:4],
																	 h = h,
																	 kernel = kernel)
			vec_temp <- vec %>% dplyr::mutate(v = list(MVKEresult(c(x, y) %>% normalize_v(dv))$mu %>% scale_up(dv)))
			vec_temp <- vec_temp %>% dplyr::mutate(vx = v[1], vy = v[2], v_norm = (sum(v^2))^(1/2)) %>%
				dplyr::select(-v) %>% dplyr::ungroup()
			result <- list(vec_grid = vec_temp, VFCresult = NULL, MVKEresult = MVKEresult,
										 data = d_raw, data_normalized = dv, original_vectors = sampled_vectors,
										 original_vectors_normalized = sampled_vectors_normalized,
										 x = x, y = y, lims = lims, n = n, method = method)
			class(result) <- "vectorfield"

			if(add_linear_interp) {
				result <- add_interp_grid(result)
			}
			p()
			return(result)
		},
		.options = furrr::furrr_options(seed = seed)
	)

	return(structure(list(
		bootstrap_models = bootstrap_models,
		original_vf = vf,
		block_length = block_length,
		n_boot = n_boot
	), class = "bootstrap_2d_vf"))
}

#' Generates bootstrap landscapes from bootstrap vector fields.
#'
#' @param boot_vf A `bootstrap_2d_vf` object containing bootstrap vector fields.
#' @param ... Additional arguments passed to `make_2d_ld`.
#'
#' @return An object of class `bootstrap_2d_ld` containing:
#'        - `bootstrap_lds`: A list of 2D landscapes from each bootstrap
#'        - `original_ld`: The original 2D landscape from the original vector field.
#'        - `n_boot`: The number of bootstrap samples.
#'
#' @export
bootstrap_2d_ld <- function(boot_vf, ...) {
	# check the class of boot_vf
	if (!inherits(boot_vf, "bootstrap_2d_vf")) {
		cli::cli_abort("Input 'boot_vf' must be a 'bootstrap_2d_vf' object.")
	}

	p <- progressr::progressor(steps = boot_vf$n_boot)

	# make landscapes from each bootstrap vector field
	# use future.apply for parallel processing

	boot_lds <- future.apply::future_lapply(boot_vf$bootstrap_models, function(vf){
		p()
		purrr::quietly(make_2d_ld)(vf, ...)$result})

	return(structure(list(
		bootstrap_lds = boot_lds,
		original_ld = make_2d_ld(boot_vf$original_vf),
		n_boot = boot_vf$n_boot
	), class = "bootstrap_2d_ld"))
}

#' @export
#' @rdname bootstrap_2d_ld
#' @param object A `bootstrap_2d_ld` object.
#' @param exclude_minor Logical indicating whether to exclude minor local minima. Default TRUE.
#' @param ... Additional arguments (not used).
summary.bootstrap_2d_ld <- function(object, exclude_minor = TRUE, ...) {
	# Find local minima for each bootstrap landscape as well as their U values
	# and positions.
	p <- progressr::progressor(steps = length(object$bootstrap_lds))

	lds <- object$bootstrap_lds

	boot_mins <- furrr::future_map(
		lds,
		function(ld) {
			p()
			find_loc_min(ld, exclude_minor = exclude_minor)
		},
		.progress = .progress
	)

	# reduce it into a data frame, with columns: boot_index, min_index, x, y, U
	boot_min_df <- do.call(rbind, lapply(1:length(boot_mins), function(i){
		mins <- boot_mins[[i]]$mins
		if(exclude_minor) {
			mins <- mins %>% dplyr::filter(!is_minor)
		}

		if(nrow(mins) == 0){
			return(NULL)
		}
		data.frame(
			boot_index = i,
			min_index = 1:nrow(mins),
			x = mins$x,
			y = mins$y,
			U = mins$U
		)
	}))

	# perform clustering based on all local minima of all bootstrapping samples
	# using DBSCAN

	dbscan_result <- dbscan::dbscan(
		boot_min_df[, c("x", "y")],
		eps = 0.05,
		minPts = 5
	)

	return(structure(list(
		boot_min_df = boot_min_df,
		original_ld = object$original_ld,
		n_boot = object$n_boot
	), class = "summary_bootstrap_2d_ld"))
}

# process_single_ld <- function(ld, p_func, find_func, ...) {
#   p_func()
#   find_func(ld, ...)
# }

#' Plots all the bootstrap minima positions
#'
#' Color represents the U value at that minimum.
#'
#' @param x A `summary_bootstrap_2d_ld` object.
plot.summary_bootstrap_2d_ld <- function(x, ...) {
	# jitter for a small value because some minima might overlap
	# how much jitter depends on the size of the landscape grid
	# available from the x_coords and y_coords attributes from original_ld$ss

	x_range <- range(attr(x$original_ld$ss, "x_coords"))
	y_range <- range(attr(x$original_ld$ss, "y_coords"))
	x_jitter_amount <- (x_range[2] - x_range[1]) * 0.02
	y_jitter_amount <- (y_range[2] - y_range[1]) * 0.02

	# the shape of the dot represents how many local minima are found in that bootstrap sample
	# so first calculate the counts

	boot_counts <- table(x$boot_min_df$boot_index)
	boot_count_df <- data.frame(
		boot_index = as.integer(names(boot_counts)),
		count = as.integer(boot_counts)
	)

	plot_df <- x$boot_min_df %>%
		dplyr::left_join(boot_count_df, by = "boot_index")

	# now make the plot
	# make sure that the shapes used are all solid shapes (i.e., don't use hollow shapes)
	# for better visibility
	# if there are more than 4 shapes needed, combine the shapes representing higher counts
	max_count <- max(plot_df$count)
	if(max_count > 4) {
		plot_df <- plot_df %>%
			dplyr::mutate(count = ifelse(count >= 4, 4, count))
		cli::cli_inform("More than 4 minima found in some bootstrap samples.
                  Combining counts >=4 into a single category for plotting.")
	}

	ggplot2::ggplot(
		plot_df,
		ggplot2::aes(x = x + runif(nrow(plot_df), -x_jitter_amount, x_jitter_amount),
								 y = y + runif(nrow(plot_df), -y_jitter_amount, y_jitter_amount),
								 color = U,
								 shape = as.factor(count))
	) +
		ggplot2::geom_point(size = 2, alpha = 0.7) +
		ggplot2::scale_color_viridis_c() +
		ggplot2::labs(
			x = x$original_ld$vf$x,
			y = x$original_ld$vf$y,
			color = "U value",
			shape = "Number of minima\nin bootstrap sample"
		) +
		ggplot2::theme_bw() +
		ggplot2::scale_shape_manual(values = c(16, 17, 15, 18, 3))
}
