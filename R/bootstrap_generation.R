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
bootstrap_2d_vf <- function(vf,
                            block_length = NULL,
                            n_boot = 200,
                            seed = 1614,
                            add_linear_interp = TRUE,
                            ...) {
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
  if (is.null(block_length)) {
    block_length <- ceiling(n_vec^(1 / 3))
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
  vec <- vf$vec_grid[, c("x", "y"), drop = FALSE]
  vec_xy <- as.matrix(vec)
  n_grid <- nrow(vec_xy)
  x <- vf$x
  y <- vf$y
  n <- vf$n
  method <- vf$method
  d_raw <- vf$data

  p <- progressr::progressor(steps = n_boot)
  bootstrap_models <- lapply(
    1:n_boot,
    function(b) {
      # Sample blocks with replacement
      sampled_block_indices <- sample(1:n_blocks, n_blocks_needed, replace = TRUE)
      sampled_indices <- unlist(blocks[sampled_block_indices])
      sampled_indices <- sampled_indices[sampled_indices <= n_vec] # Trim to original size
      sampled_vectors <- original_vectors[sampled_indices, ]
      sampled_vectors_normalized <- original_vectors_normalized[sampled_indices, ]
      # Fit vector field to the sampled vectors
      MVKEresult <- fitlandr::MVKE(
        d = sampled_vectors_normalized[, 1:2], v = sampled_vectors_normalized[, 3:4],
        h = h,
        kernel = kernel
      )
      v_mat <- matrix(NA_real_, nrow = n_grid, ncol = 2)
      for (idx in seq_len(n_grid)) {
        v_mat[idx, ] <- MVKEresult(normalize_v(vec_xy[idx, ], dv))$mu %>% scale_up(dv)
      }
      vec_temp <- vec
      vec_temp$vx <- v_mat[, 1]
      vec_temp$vy <- v_mat[, 2]
      vec_temp$v_norm <- sqrt(vec_temp$vx^2 + vec_temp$vy^2)

      result <- list(
        vec_grid = vec_temp, VFCresult = NULL, MVKEresult = MVKEresult,
        data = d_raw, data_normalized = dv, original_vectors = sampled_vectors,
        original_vectors_normalized = sampled_vectors_normalized,
        x = x, y = y, lims = lims, n = n, method = method
      )
      class(result) <- "vectorfield"

      if (add_linear_interp) {
        result <- add_interp_grid(result)
      }
      p()
      return(result)
    }
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
#'        - `bootstrap_lds`: A list of 2D landscapes from each bootstrap vector field. Note that the plots are removed to save space.
#'        - `original_ld`: The original 2D landscape from the original vector field.
#'        - `n_boot`: The number of bootstrap samples.
#'
#' @export
bootstrap_2d_ld <- function(boot_vf, ...) {
  # check the class of boot_vf
  if (!inherits(boot_vf, "bootstrap_2d_vf")) {
    cli::cli_abort("Input 'boot_vf' must be a 'bootstrap_2d_vf' object.")
  }

  original_ld <- purrr::quietly(make_2d_ld)(boot_vf$original_vf, ...)$result
  ref_x <- sort(unique(original_ld$dist$x))
  ref_y <- sort(unique(original_ld$dist$y))

  p <- progressr::progressor(steps = boot_vf$n_boot)

  boot_lds <- lapply(
    boot_vf$bootstrap_models,
    function(vf) {
      p()
      result <- purrr::quietly(make_2d_ld)(vf, ...)$result
      # to save space:
      result$plot <- NULL
      result$plot_2 <- NULL

      ux <- sort(unique(result$dist$x))
      uy <- sort(unique(result$dist$y))
      if (!identical(ux, ref_x) || !identical(uy, ref_y)) {
        cli::cli_abort("Bootstrap landscapes must use exactly the same grid as the original landscape.")
      }

      return(result)
    }
  )

  return(structure(list(
    bootstrap_lds = boot_lds,
    original_ld = original_ld,
    n_boot = boot_vf$n_boot
  ), class = "bootstrap_2d_ld"))
}


