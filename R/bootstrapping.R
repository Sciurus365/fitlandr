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


build_mean_potential_ld <- function(object) {
  ref_ld <- object$original_ld
  x_vals <- sort(unique(ref_ld$dist$x))
  y_vals <- sort(unique(ref_ld$dist$y))

  if (!length(object$bootstrap_lds)) {
    cli::cli_abort("bootstrap_2d_ld object contains no bootstrap landscapes.")
  }

  u_mats <- lapply(object$bootstrap_lds, function(ld) {
    ux <- sort(unique(ld$dist$x))
    uy <- sort(unique(ld$dist$y))

    if (!identical(ux, x_vals) || !identical(uy, y_vals)) {
      cli::cli_abort("All bootstrap landscapes must use the same grid to compute a mean potential surface.")
    }

    matrix(ld$dist$U, nrow = length(x_vals), ncol = length(y_vals), byrow = FALSE)
  })

  mean_u <- Reduce(`+`, u_mats) / length(u_mats)
  mean_d <- exp(-(mean_u - min(mean_u, na.rm = TRUE)))
  mean_d <- mean_d / sum(mean_d, na.rm = TRUE)

  mean_dist <- expand.grid(x = x_vals, y = y_vals)
  mean_dist$d <- as.vector(mean_d)
  mean_dist$U <- as.vector(mean_u)

  mean_ss <- mean_d
  attr(mean_ss, "x_coords") <- x_vals
  attr(mean_ss, "y_coords") <- y_vals

  ref_ld$dist <- mean_dist
  ref_ld$ss <- mean_ss
  ref_ld$plot <- NULL
  ref_ld$plot_2 <- NULL
  ref_ld
}


get_bootstrap_distance_scaler <- function(object) {
  dist <- object$original_ld$dist
  x_vals <- sort(unique(dist$x))
  y_vals <- sort(unique(dist$y))

  dx <- if (length(x_vals) > 1L) stats::median(diff(x_vals)) else 1
  dy <- if (length(y_vals) > 1L) stats::median(diff(y_vals)) else 1

  if (!is.finite(dx) || dx <= 0) {
    dx <- 1
  }
  if (!is.finite(dy) || dy <= 0) {
    dy <- 1
  }

  list(dx = dx, dy = dy)
}


standardize_xy <- function(x, y, scaler) {
  cbind(
    x = x / scaler$dx,
    y = y / scaler$dy
  )
}


compute_landscape_hessian <- function(ld, x, y) {
  dist <- ld$dist
  x_vals <- sort(unique(dist$x))
  y_vals <- sort(unique(dist$y))
  u_mat <- matrix(dist$U, nrow = length(x_vals), ncol = length(y_vals), byrow = FALSE)

  i <- which.min(abs(x_vals - x))
  j <- which.min(abs(y_vals - y))

  if (i <= 1 || i >= length(x_vals) || j <= 1 || j >= length(y_vals)) {
    return(matrix(NA_real_, nrow = 2, ncol = 2))
  }

  hx <- stats::median(diff(x_vals))
  hy <- stats::median(diff(y_vals))

  if (!is.finite(hx) || !is.finite(hy) || hx <= 0 || hy <= 0) {
    return(matrix(NA_real_, nrow = 2, ncol = 2))
  }

  dxx <- (u_mat[i + 1, j] - 2 * u_mat[i, j] + u_mat[i - 1, j]) / (hx^2)
  dyy <- (u_mat[i, j + 1] - 2 * u_mat[i, j] + u_mat[i, j - 1]) / (hy^2)
  dxy <- (u_mat[i + 1, j + 1] - u_mat[i + 1, j - 1] - u_mat[i - 1, j + 1] + u_mat[i - 1, j - 1]) / (4 * hx * hy)

  matrix(c(dxx, dxy, dxy, dyy), nrow = 2, byrow = TRUE)
}


ellipse_from_hessian <- function(hessian, level, n_runs) {
  out <- list(
    s_xx = NA_real_,
    s_yy = NA_real_,
    s_xy = NA_real_,
    angle = NA_real_,
    a_pred = NA_real_,
    b_pred = NA_real_,
    a_conf = NA_real_,
    b_conf = NA_real_,
    hessian_ok = FALSE
  )

  if (length(hessian) != 4 || any(!is.finite(hessian))) {
    return(out)
  }

  hessian <- (hessian + t(hessian)) / 2
  cov_mat <- tryCatch(solve(hessian), error = function(e) NULL)
  if (is.null(cov_mat) || any(!is.finite(cov_mat))) {
    return(out)
  }

  eig <- eigen((cov_mat + t(cov_mat)) / 2, symmetric = TRUE)
  lambda <- eig$values
  if (any(!is.finite(lambda)) || any(lambda <= 0)) {
    return(out)
  }

  c2 <- stats::qchisq(level, df = 2)
  n_runs <- as.integer(n_runs)
  if (!is.finite(n_runs) || n_runs < 1L) {
    n_runs <- 1L
  }

  out$s_xx <- cov_mat[1, 1]
  out$s_yy <- cov_mat[2, 2]
  out$s_xy <- cov_mat[1, 2]
  out$angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
  out$a_pred <- sqrt(lambda[1] * c2)
  out$b_pred <- sqrt(lambda[2] * c2)
  out$a_conf <- sqrt(lambda[1] * c2 / n_runs)
  out$b_conf <- sqrt(lambda[2] * c2 / n_runs)
  out$hessian_ok <- TRUE
  out
}


summarize_mean_potential_hessian <- function(object,
                                             boot_min_df,
                                             exclude_minor,
                                             min_barrier,
                                             level,
                                             one_per_run) {
  scaler <- get_bootstrap_distance_scaler(object)

  mean_ld <- build_mean_potential_ld(object)
  ref_mins <- find_loc_min(
    mean_ld,
    exclude_minor = exclude_minor,
    min_barrier = min_barrier
  )$mins

  if (exclude_minor && nrow(ref_mins) > 0) {
    ref_mins <- ref_mins[!ref_mins$is_minor, , drop = FALSE]
  }

  per_boot <- tibble::tibble(boot_index = seq_len(object$n_boot)) |>
    dplyr::left_join(
      boot_min_df |>
        dplyr::count(boot_index, name = "n_mins"),
      by = "boot_index"
    ) |>
    dplyr::mutate(n_mins = dplyr::coalesce(.data$n_mins, 0L))

  if (!nrow(ref_mins)) {
    out <- list(
      params = list(
        exclude_minor = exclude_minor,
        min_barrier = min_barrier,
        clustering_method = "mean_potential_hessian",
        minPts = NA_integer_,
        level = level,
        one_per_run = one_per_run
      ),
      per_boot = per_boot,
      per_point = boot_min_df |>
        dplyr::mutate(cluster = 0L, is_noise = TRUE),
      per_cluster = NULL,
      diagnostics = list(
        clustering_method = "mean_potential_hessian",
        n_reference = 0L,
        distance_scale_dx = scaler$dx,
        distance_scale_dy = scaler$dy,
        message = "No reference minima found on mean potential surface"
      ),
      original_ld = object$original_ld,
      n_boot = object$n_boot
    )
    return(structure(out, class = "summary_bootstrap_2d_ld"))
  }

  ref_xy <- standardize_xy(ref_mins$x, ref_mins$y, scaler = scaler)
  per_point <- dplyr::bind_rows(lapply(seq_len(object$n_boot), function(i) {
    df_run <- boot_min_df[boot_min_df$boot_index == i, , drop = FALSE]
    if (!nrow(df_run)) {
      return(NULL)
    }

    run_xy <- standardize_xy(df_run$x, df_run$y, scaler = scaler)
    df_run$cluster <- assign_run_to_reference(run_xy, ref_xy)
    df_run$is_noise <- df_run$cluster == 0L
    df_run
  }))

  if (is.null(per_point) || !nrow(per_point)) {
    per_point <- boot_min_df |>
      dplyr::mutate(cluster = 0L, is_noise = TRUE)
  }

  hessian_rows <- lapply(seq_len(nrow(per_point)), function(i) {
    row_i <- per_point[i, , drop = FALSE]
    if (isTRUE(row_i$is_noise[[1]])) {
      row_i$h_xx <- NA_real_
      row_i$h_xy <- NA_real_
      row_i$h_yy <- NA_real_
      row_i$hessian_ok <- FALSE
      return(row_i)
    }

    hessian <- compute_landscape_hessian(
      object$bootstrap_lds[[row_i$boot_index[[1]]]],
      x = row_i$x[[1]],
      y = row_i$y[[1]]
    )

    row_i$h_xx <- hessian[1, 1]
    row_i$h_xy <- hessian[1, 2]
    row_i$h_yy <- hessian[2, 2]
    row_i$hessian_ok <- all(is.finite(hessian))
    row_i
  })
  per_point <- dplyr::bind_rows(hessian_rows)

  per_cluster <- dplyr::bind_rows(lapply(seq_len(nrow(ref_mins)), function(cluster_id) {
    matched <- per_point |>
      dplyr::filter(cluster == cluster_id, !is_noise)
    matched_h <- matched |>
      dplyr::filter(.data$hessian_ok)

    h_bar <- if (nrow(matched_h)) {
      matrix(
        c(
          mean(matched_h$h_xx), mean(matched_h$h_xy),
          mean(matched_h$h_xy), mean(matched_h$h_yy)
        ),
        nrow = 2,
        byrow = TRUE
      )
    } else {
      matrix(NA_real_, nrow = 2, ncol = 2)
    }

    ellipse <- ellipse_from_hessian(
      hessian = h_bar,
      level = level,
      n_runs = dplyr::n_distinct(matched$boot_index)
    )

    tibble::tibble(
      cluster = cluster_id,
      n_points = nrow(matched),
      n_runs = dplyr::n_distinct(matched$boot_index),
      n_hessian = nrow(matched_h),
      mean_x = ref_mins$x[[cluster_id]],
      mean_y = ref_mins$y[[cluster_id]],
      mean_U = ref_mins$U[[cluster_id]],
      sd_x = if (nrow(matched) > 1) stats::sd(matched$x) else NA_real_,
      sd_y = if (nrow(matched) > 1) stats::sd(matched$y) else NA_real_,
      sd_U = if (nrow(matched) > 1) stats::sd(matched$U) else NA_real_,
      CI_U_lower = if (nrow(matched)) as.numeric(stats::quantile(matched$U, probs = (1 - level) / 2)) else NA_real_,
      CI_U_upper = if (nrow(matched)) as.numeric(stats::quantile(matched$U, probs = 1 - (1 - level) / 2)) else NA_real_,
      s_xx = ellipse$s_xx,
      s_yy = ellipse$s_yy,
      s_xy = ellipse$s_xy,
      stability = dplyr::n_distinct(matched$boot_index) / object$n_boot,
      angle = ellipse$angle,
      a_pred = ellipse$a_pred,
      b_pred = ellipse$b_pred,
      a_conf = ellipse$a_conf,
      b_conf = ellipse$b_conf,
      h_xx = h_bar[1, 1],
      h_xy = h_bar[1, 2],
      h_yy = h_bar[2, 2],
      hessian_ok = ellipse$hessian_ok
    )
  }))

  out <- list(
    params = list(
      exclude_minor = exclude_minor,
      min_barrier = min_barrier,
      clustering_method = "mean_potential_hessian",
      minPts = NA_integer_,
      level = level,
      one_per_run = one_per_run
    ),
    per_boot = per_boot,
    per_point = per_point,
    per_cluster = per_cluster,
    diagnostics = list(
      clustering_method = "mean_potential_hessian",
      n_reference = nrow(ref_mins),
      distance_scale_dx = scaler$dx,
      distance_scale_dy = scaler$dy,
      matched_runs_per_reference = per_cluster$n_runs,
      valid_hessian_per_reference = per_cluster$n_hessian,
      invalid_hessian_points = sum(!per_point$is_noise & !per_point$hessian_ok)
    ),
    original_ld = object$original_ld,
    n_boot = object$n_boot
  )

  structure(out, class = "summary_bootstrap_2d_ld")
}


summarize_mean_potential <- function(object,
                                     boot_min_df,
                                     exclude_minor,
                                     min_barrier,
                                     level,
                                     one_per_run) {
  scaler <- get_bootstrap_distance_scaler(object)

  mean_ld <- build_mean_potential_ld(object)
  ref_mins <- find_loc_min(
    mean_ld,
    exclude_minor = exclude_minor,
    min_barrier = min_barrier
  )$mins

  if (exclude_minor && nrow(ref_mins) > 0) {
    ref_mins <- ref_mins[!ref_mins$is_minor, , drop = FALSE]
  }

  per_boot <- tibble::tibble(boot_index = seq_len(object$n_boot)) |>
    dplyr::left_join(
      boot_min_df |>
        dplyr::count(boot_index, name = "n_mins"),
      by = "boot_index"
    ) |>
    dplyr::mutate(n_mins = dplyr::coalesce(.data$n_mins, 0L))

  if (!nrow(ref_mins)) {
    out <- list(
      params = list(
        exclude_minor = exclude_minor,
        min_barrier = min_barrier,
        clustering_method = "mean_potential",
        minPts = NA_integer_,
        level = level,
        one_per_run = one_per_run
      ),
      per_boot = per_boot,
      per_point = boot_min_df |>
        dplyr::mutate(cluster = 0L, is_noise = TRUE),
      per_cluster = NULL,
      diagnostics = list(
        clustering_method = "mean_potential",
        n_reference = 0L,
        distance_scale_dx = scaler$dx,
        distance_scale_dy = scaler$dy,
        message = "No reference minima found on mean potential surface"
      ),
      original_ld = object$original_ld,
      n_boot = object$n_boot
    )
    return(structure(out, class = "summary_bootstrap_2d_ld"))
  }

  ref_xy <- standardize_xy(ref_mins$x, ref_mins$y, scaler = scaler)
  per_point <- dplyr::bind_rows(lapply(seq_len(object$n_boot), function(i) {
    df_run <- boot_min_df[boot_min_df$boot_index == i, , drop = FALSE]
    if (!nrow(df_run)) {
      return(NULL)
    }

    run_xy <- standardize_xy(df_run$x, df_run$y, scaler = scaler)
    df_run$cluster <- assign_run_to_reference(run_xy, ref_xy)
    df_run$is_noise <- df_run$cluster == 0L
    df_run
  }))

  if (is.null(per_point) || !nrow(per_point)) {
    per_point <- boot_min_df |>
      dplyr::mutate(cluster = 0L, is_noise = TRUE)
  }

  per_cluster <- dplyr::bind_rows(lapply(seq_len(nrow(ref_mins)), function(cluster_id) {
    matched <- per_point |>
      dplyr::filter(cluster == cluster_id, !is_noise)

    s_xx <- if (nrow(matched) > 1L) stats::var(matched$x) else NA_real_
    s_yy <- if (nrow(matched) > 1L) stats::var(matched$y) else NA_real_
    s_xy <- if (nrow(matched) > 1L) stats::cov(matched$x, matched$y) else NA_real_

    sigma <- matrix(c(s_xx, s_xy, s_xy, s_yy), nrow = 2, byrow = TRUE)
    eig <- tryCatch(eigen(sigma, symmetric = TRUE), error = function(e) NULL)

    c2 <- stats::qchisq(level, df = 2)
    n_runs_i <- dplyr::n_distinct(matched$boot_index)

    if (is.null(eig)) {
      angle <- NA_real_
      a_pred <- NA_real_
      b_pred <- NA_real_
      a_conf <- NA_real_
      b_conf <- NA_real_
    } else {
      lambda <- pmax(eig$values, 0)
      angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
      a_pred <- sqrt(lambda[1] * c2)
      b_pred <- sqrt(lambda[2] * c2)
      a_conf <- sqrt(lambda[1] * c2 / max(n_runs_i, 1L))
      b_conf <- sqrt(lambda[2] * c2 / max(n_runs_i, 1L))
    }

    tibble::tibble(
      cluster = cluster_id,
      n_points = nrow(matched),
      n_runs = n_runs_i,
      mean_x = ref_mins$x[[cluster_id]],
      mean_y = ref_mins$y[[cluster_id]],
      mean_U = ref_mins$U[[cluster_id]],
      sd_x = if (nrow(matched) > 1) stats::sd(matched$x) else NA_real_,
      sd_y = if (nrow(matched) > 1) stats::sd(matched$y) else NA_real_,
      sd_U = if (nrow(matched) > 1) stats::sd(matched$U) else NA_real_,
      CI_U_lower = if (nrow(matched)) as.numeric(stats::quantile(matched$U, probs = (1 - level) / 2)) else NA_real_,
      CI_U_upper = if (nrow(matched)) as.numeric(stats::quantile(matched$U, probs = 1 - (1 - level) / 2)) else NA_real_,
      s_xx = s_xx,
      s_yy = s_yy,
      s_xy = s_xy,
      stability = n_runs_i / object$n_boot,
      angle = angle,
      a_pred = a_pred,
      b_pred = b_pred,
      a_conf = a_conf,
      b_conf = b_conf
    )
  }))

  out <- list(
    params = list(
      exclude_minor = exclude_minor,
      min_barrier = min_barrier,
      clustering_method = "mean_potential",
      minPts = NA_integer_,
      level = level,
      one_per_run = one_per_run
    ),
    per_boot = per_boot,
    per_point = per_point,
    per_cluster = per_cluster,
    diagnostics = list(
      clustering_method = "mean_potential",
      n_reference = nrow(ref_mins),
      distance_scale_dx = scaler$dx,
      distance_scale_dy = scaler$dy,
      matched_runs_per_reference = per_cluster$n_runs
    ),
    original_ld = object$original_ld,
    n_boot = object$n_boot
  )

  structure(out, class = "summary_bootstrap_2d_ld")
}


#' @rdname bootstrap_2d_ld
#'
#' @param object A `bootstrap_2d_ld` object with fields:
#'   - `bootstrap_lds`: list of landscapes,
#'   - `original_ld`: the original landscape,
#'   - `n_boot`: number of bootstrap runs.
#' @param exclude_minor Logical; exclude minor local minima. Default TRUE.
#' @param min_barrier Minimum barrier height fraction used by `find_loc_min()` to classify minor minima. Default 0.1.
#' @param clustering_method Clustering backend for pooled minima. One of
#'   `"hdbscan"` (default), `"hungarian"`, `"gmm_bic"`, `"mean_potential"`, `"mean_potential_hessian"`,
#'   `"pairwise_hungarian_graph"`.
#' @param minPts Integer; HDBSCAN minPts. Default 5.
#' @param pairwise_leiden_gamma Numeric Leiden CPM resolution parameter used only
#'   for `"pairwise_hungarian_graph"`. Default 0.01.
#' @param level Confidence level for ellipses (e.g., 0.95). Default 0.95.
#' @param one_per_run Logical; at most one point per run per cluster for summary stats. Default TRUE.
#' @param ... Unused.
#'
#' @return An object of class `"summary_bootstrap_2d_ld"` with components:
#'   - `params`, `per_boot`, `per_point`,
#'   - `per_cluster`: now includes `a_pred`, `b_pred`, `a_conf`, `b_conf`, `angle`, `c2`,
#'   - `diagnostics`, `original_ld`, `n_boot`.
#'
#' @details
#' For clustering-based methods, ellipse parameters are derived from the eigen-decomposition
#' of the 2×2 covariance matrix per cluster. For `"mean_potential"`, minima are first
#' detected on the mean potential surface across bootstrap landscapes; each bootstrap minimum is
#' then matched to the closest reference minimum, and ellipse parameters are derived from the
#' covariance of matched minima locations. For `"mean_potential_hessian"`, minima are first
#' detected on the mean potential surface across bootstrap landscapes; each bootstrap minimum is
#' then matched to the closest reference minimum, and ellipse parameters are derived from the
#' inverse of the mean matched Hessian.
#'
#' Semi-axes for the **prediction** ellipse are
#' \eqn{\sqrt{\lambda_i\,\chi^2_{2,\alpha}}}; for the **confidence** ellipse they are
#' \eqn{\sqrt{\lambda_i\,\chi^2_{2,\alpha}/n_{\text{runs}}}} where \eqn{n_{\text{runs}}}
#' is the number of bootstrap runs contributing to that minimum.
#' Ellipses are drawn with `ggforce::geom_ellipse()` which expects aesthetics
#' `x0`, `y0`, `a`, `b`, `angle`.  \[See ggforce docs.\]  # (geom_ellipse API)  [1](https://rstudio-pubs-static.s3.amazonaws.com/1236382_7016d680936d411e8fd45fc0b8b62b0c.html)
#' The \eqn{\chi^2_{2,\alpha}} quantile comes from the chi-square distribution
#' (e.g., 0.95 → 5.991).  \[See NIST/ITL table.\]  # (chi-square critical)  [2](https://stackoverflow.com/questions/70010774/dbscan-choice-of-epsilon-through-elbow-method)
#'
#' @export
#' @method summary bootstrap_2d_ld
summary.bootstrap_2d_ld <- function(object,
                                    exclude_minor = TRUE,
                                    min_barrier = 0.1,
                                    clustering_method = c("hdbscan", "hungarian", "gmm_bic", "mean_potential", "mean_potential_hessian", "pairwise_hungarian_graph"),
                                    minPts = 5,
                                    pairwise_leiden_gamma = 0.01,
                                    level = 0.95,
                                    one_per_run = TRUE,
                                    ...) {
  if (!is.list(object) || is.null(object$bootstrap_lds) || is.null(object$n_boot)) {
    cli::cli_abort("{.arg object} must be a {.cls bootstrap_2d_ld} object with {.field bootstrap_lds} and {.field n_boot}.")
  }
  clustering_method <- rlang::arg_match0(
    clustering_method,
    c("hdbscan", "hungarian", "gmm_bic", "mean_potential", "mean_potential_hessian", "pairwise_hungarian_graph")
  )

  # ---------- 1) Collect minima across bootstrap runs ----------
  p <- progressr::progressor(steps = length(object$bootstrap_lds))
  lds <- object$bootstrap_lds

  boot_mins <- lapply(lds, function(ld) {
    p()
    find_loc_min(ld, exclude_minor = exclude_minor, min_barrier = min_barrier)
  })

  boot_min_df <- do.call(rbind, lapply(seq_along(boot_mins), function(i) {
    mins <- boot_mins[[i]]$mins
    if (exclude_minor && nrow(mins) > 0) mins <- mins[!mins$is_minor, , drop = FALSE]
    if (!nrow(mins)) {
      return(NULL)
    }
    data.frame(
      boot_index = i,
      x = mins$x, y = mins$y, U = mins$U,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(boot_min_df) || !nrow(boot_min_df)) {
    out <- list(
      params = list(
        exclude_minor = exclude_minor, min_barrier = min_barrier, clustering_method = clustering_method,
        minPts = minPts, pairwise_leiden_gamma = pairwise_leiden_gamma, level = level,
        one_per_run = one_per_run
      ),
      per_boot = data.frame(boot_index = seq_len(object$n_boot), n_mins = 0L),
      per_point = NULL,
      per_cluster = NULL,
      diagnostics = list(message = "No minima found"),
      original_ld = object$original_ld,
      n_boot = object$n_boot
    )
    return(structure(out, class = "summary_bootstrap_2d_ld"))
  }

  if (identical(clustering_method, "mean_potential")) {
    return(summarize_mean_potential(
      object = object,
      boot_min_df = boot_min_df,
      exclude_minor = exclude_minor,
      min_barrier = min_barrier,
      level = level,
      one_per_run = one_per_run
    ))
  }

  if (identical(clustering_method, "mean_potential_hessian")) {
    return(summarize_mean_potential_hessian(
      object = object,
      boot_min_df = boot_min_df,
      exclude_minor = exclude_minor,
      min_barrier = min_barrier,
      level = level,
      one_per_run = one_per_run
    ))
  }

  # ---------- 2) Cluster pooled minima ----------
  cl_out <- cluster_bootstrap_minima(
    boot_min_df = boot_min_df,
    object = object,
    exclude_minor = exclude_minor,
    min_barrier = min_barrier,
    clustering_method = clustering_method,
    minPts = minPts,
    pairwise_leiden_gamma = pairwise_leiden_gamma
  )
  boot_min_df <- cl_out$boot_min_df

  # per-bootstrap counts
  per_boot <- boot_min_df |>
    dplyr::group_by(boot_index) |>
    dplyr::summarise(n_mins = dplyr::n(), .groups = "drop")

  # Exclude noise for per-cluster stats
  df_c <- dplyr::filter(boot_min_df, !is_noise)
  if (!nrow(df_c)) {
    out <- list(
      params = list(
        exclude_minor = exclude_minor, min_barrier = min_barrier, clustering_method = clustering_method,
        minPts = minPts, pairwise_leiden_gamma = pairwise_leiden_gamma, level = level,
        one_per_run = one_per_run
      ),
      per_boot = per_boot,
      per_point = boot_min_df,
      per_cluster = NULL,
      diagnostics = c(list(noise_frac = mean(boot_min_df$is_noise)), cl_out$diagnostics),
      original_ld = object$original_ld,
      n_boot = object$n_boot
    )
    return(structure(out, class = "summary_bootstrap_2d_ld"))
  }

  # Optionally: at most one point per run per cluster (closest to provisional center)
  if (one_per_run) {
    centers0 <- df_c |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(
        cx = mean(x),
        cy = mean(y),
        .groups = "drop"
      )
    df_c <- df_c |>
      dplyr::left_join(centers0, by = "cluster") |>
      dplyr::mutate(d2 = (x - cx)^2 + (y - cy)^2) |>
      dplyr::group_by(cluster, boot_index) |>
      dplyr::slice_min(order_by = d2, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::select(-dplyr::any_of(c("cx", "cy", "d2")))
  }

  # ---------- 3) Per-cluster summaries (variances/covariance) ----------
  per_cluster <- df_c |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      n_points = dplyr::n(),
      n_runs = dplyr::n_distinct(boot_index),
      mean_x = mean(x),
      mean_y = mean(y),
      mean_U = mean(U),
      sd_x = stats::sd(x),
      sd_y = stats::sd(y),
      sd_U = stats::sd(U),
      CI_U_lower = stats::quantile(U, probs = (1 - level) / 2),
      CI_U_upper = stats::quantile(U, probs = 1 - (1 - level) / 2),
      s_xx = stats::var(x),
      s_yy = stats::var(y),
      s_xy = stats::cov(x, y),
      .groups = "drop"
    ) |>
    dplyr::mutate(stability = n_runs / object$n_boot)

  # ---------- 4) Add ellipse parameters (prediction & confidence) ----------
  c2 <- stats::qchisq(level, df = 2) # e.g., 0.95 -> 5.991  [2](https://stackoverflow.com/questions/70010774/dbscan-choice-of-epsilon-through-elbow-method)

  ellipse_df <- lapply(seq_len(nrow(per_cluster)), function(i) {
    sigma <- matrix(c(per_cluster$s_xx[[i]], per_cluster$s_xy[[i]], per_cluster$s_xy[[i]], per_cluster$s_yy[[i]]), 2, 2)
    eig <- eigen(sigma, symmetric = TRUE)
    lambda <- pmax(eig$values, 0)
    angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
    n_runs_i <- per_cluster$n_runs[[i]]
    data.frame(
      angle = angle,
      a_pred = sqrt(lambda[1] * c2),
      b_pred = sqrt(lambda[2] * c2),
      a_conf = sqrt(lambda[1] * c2 / n_runs_i),
      b_conf = sqrt(lambda[2] * c2 / n_runs_i)
    )
  })
  per_cluster <- dplyr::bind_cols(per_cluster, dplyr::bind_rows(ellipse_df))

  out <- list(
    params = list(
      exclude_minor = exclude_minor, min_barrier = min_barrier, clustering_method = clustering_method,
      minPts = minPts, pairwise_leiden_gamma = pairwise_leiden_gamma, level = level,
      one_per_run = one_per_run
    ),
    per_boot = per_boot,
    per_point = boot_min_df,
    per_cluster = per_cluster,
    diagnostics = c(list(noise_frac = mean(boot_min_df$is_noise)), cl_out$diagnostics),
    original_ld = object$original_ld,
    n_boot = object$n_boot
  )
  structure(out, class = "summary_bootstrap_2d_ld")
}


solve_lsap_with_dummies <- function(cost, dummy_cost) {
  if (!is.matrix(cost) || !all(is.finite(cost))) {
    return(data.frame(i = integer(0), j = integer(0), d = numeric(0)))
  }

  n_a <- nrow(cost)
  n_b <- ncol(cost)
  if (n_a < 1L || n_b < 1L) {
    return(data.frame(i = integer(0), j = integer(0), d = numeric(0)))
  }

  # Square assignment with dummy rows/cols allows explicit unmatched minima.
  n_tot <- n_a + n_b
  big <- matrix(0, nrow = n_tot, ncol = n_tot)
  big[seq_len(n_a), seq_len(n_b)] <- cost
  big[seq_len(n_a), n_b + seq_len(n_a)] <- dummy_cost
  big[n_a + seq_len(n_b), seq_len(n_b)] <- dummy_cost

  assignment <- as.integer(clue::solve_LSAP(big))
  rows <- seq_len(n_a)
  cols <- assignment[rows]
  keep <- cols <= n_b
  if (!any(keep)) {
    return(data.frame(i = integer(0), j = integer(0), d = numeric(0)))
  }

  rows <- rows[keep]
  cols <- cols[keep]
  d <- cost[cbind(rows, cols)]
  keep2 <- is.finite(d) & d <= dummy_cost

  data.frame(
    i = rows[keep2],
    j = cols[keep2],
    d = d[keep2]
  )
}


connected_components_from_edges <- function(n_nodes, edge_u, edge_v) {
  cluster <- integer(n_nodes)
  if (n_nodes < 1L || length(edge_u) < 1L) {
    return(cluster)
  }

  adj <- vector("list", n_nodes)
  for (k in seq_along(edge_u)) {
    u <- as.integer(edge_u[[k]])
    v <- as.integer(edge_v[[k]])
    if (u < 1L || v < 1L || u > n_nodes || v > n_nodes || u == v) {
      next
    }
    adj[[u]] <- c(adj[[u]], v)
    adj[[v]] <- c(adj[[v]], u)
  }

  visited <- rep(FALSE, n_nodes)
  cid <- 0L
  for (start in seq_len(n_nodes)) {
    if (visited[[start]] || length(adj[[start]]) == 0L) {
      next
    }
    cid <- cid + 1L
    queue <- start
    visited[[start]] <- TRUE
    cluster[[start]] <- cid

    while (length(queue) > 0L) {
      u <- queue[[1]]
      queue <- queue[-1]
      nbrs <- unique(adj[[u]])
      for (v in nbrs) {
        if (!visited[[v]]) {
          visited[[v]] <- TRUE
          cluster[[v]] <- cid
          queue <- c(queue, v)
        }
      }
    }
  }

  cluster
}


cluster_pairwise_hungarian_graph <- function(boot_min_df, pairwise_leiden_gamma = 0.01) {
  n_points <- nrow(boot_min_df)
  if (n_points < 2L) {
    boot_min_df$cluster <- 0L
    boot_min_df$is_noise <- TRUE
    return(list(
      boot_min_df = boot_min_df,
      diagnostics = list(clustering_method = "pairwise_hungarian_graph", n_edges = 0L)
    ))
  }

  by_run <- split(boot_min_df, boot_min_df$boot_index)
  run_ids <- names(by_run)
  if (length(run_ids) < 2L) {
    boot_min_df$cluster <- 0L
    boot_min_df$is_noise <- TRUE
    return(list(
      boot_min_df = boot_min_df,
      diagnostics = list(clustering_method = "pairwise_hungarian_graph", n_edges = 0L)
    ))
  }

  edge_rows <- list()
  edge_counter <- 0L

  for (a in seq_len(length(run_ids) - 1L)) {
    for (b in (a + 1L):length(run_ids)) {
      df_a <- by_run[[a]]
      df_b <- by_run[[b]]
      xy_a <- as.matrix(df_a[, c("x_scaled", "y_scaled"), drop = FALSE])
      xy_b <- as.matrix(df_b[, c("x_scaled", "y_scaled"), drop = FALSE])

      cost <- outer(
        seq_len(nrow(xy_a)), seq_len(nrow(xy_b)),
        FUN = function(i, j) {
          dx <- xy_a[i, 1] - xy_b[j, 1]
          dy <- xy_a[i, 2] - xy_b[j, 2]
          sqrt(dx^2 + dy^2)
        }
      )

      # Adaptive dummy penalty so poor cross-run matches can stay unmatched.
      dvals <- as.numeric(cost)
      dvals <- dvals[is.finite(dvals)]
      if (!length(dvals)) {
        next
      }
      dummy_cost <- max(0.1, as.numeric(stats::quantile(dvals, probs = 0.75, na.rm = TRUE)))

      matches <- solve_lsap_with_dummies(cost = cost, dummy_cost = dummy_cost)
      if (!nrow(matches)) {
        next
      }

      edge_counter <- edge_counter + 1L
      edge_rows[[edge_counter]] <- data.frame(
        node_u = df_a$node_id[matches$i],
        node_v = df_b$node_id[matches$j],
        d = matches$d
      )
    }
  }

  if (!length(edge_rows)) {
    boot_min_df$cluster <- 0L
    boot_min_df$is_noise <- TRUE
    return(list(
      boot_min_df = boot_min_df,
      diagnostics = list(clustering_method = "pairwise_hungarian_graph", n_edges = 0L)
    ))
  }

  edges <- dplyr::bind_rows(edge_rows)
  sigma <- stats::median(edges$d, na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) {
    positive_d <- edges$d[is.finite(edges$d) & edges$d > 0]
    sigma <- if (length(positive_d)) {
      stats::median(positive_d, na.rm = TRUE)
    } else {
      1e-9
    }
  }
  edges$weight <- exp(-(edges$d^2) / (sigma^2))

  d_thr <- as.numeric(stats::quantile(edges$d, probs = 0.75, na.rm = TRUE))
  strong_edges <- edges[edges$d <= d_thr, , drop = FALSE]
  if (!nrow(strong_edges)) {
    strong_edges <- edges
    d_thr <- max(edges$d, na.rm = TRUE)
  }

  leiden_used <- FALSE
  leiden_error <- NULL
  cluster_vec <- rep(0L, n_points)

  if (requireNamespace("igraph", quietly = TRUE) && "cluster_leiden" %in% getNamespaceExports("igraph")) {
    g <- igraph::graph_from_data_frame(
      d = strong_edges[, c("node_u", "node_v", "weight")],
      directed = FALSE,
      vertices = data.frame(name = as.character(seq_len(n_points)))
    )

    lc <- tryCatch(
      igraph::cluster_leiden(
        g,
        weights = igraph::E(g)$weight,
        objective_function = "CPM",
        resolution = pairwise_leiden_gamma
      ),
      error = function(e) {
        leiden_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(lc)) {
      memb <- as.integer(igraph::membership(lc))
      if (length(memb) == n_points) {
        cluster_vec <- memb
        leiden_used <- TRUE
      }
    }
  }

  if (!leiden_used) {
    cluster_vec <- connected_components_from_edges(
      n_nodes = n_points,
      edge_u = strong_edges$node_u,
      edge_v = strong_edges$node_v
    )
  }

  boot_min_df$cluster <- as.integer(cluster_vec[boot_min_df$node_id])

  # Drop tiny components observed in only one bootstrap run.
  valid <- boot_min_df |>
    dplyr::filter(cluster > 0L) |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(n_runs = dplyr::n_distinct(boot_index), .groups = "drop") |>
    dplyr::filter(n_runs >= 2L)

  valid_clusters <- valid$cluster
  boot_min_df$cluster[!(boot_min_df$cluster %in% valid_clusters)] <- 0L

  uniq <- sort(unique(boot_min_df$cluster[boot_min_df$cluster > 0L]))
  if (length(uniq)) {
    remap <- setNames(seq_along(uniq), uniq)
    boot_min_df$cluster[boot_min_df$cluster > 0L] <- as.integer(remap[as.character(boot_min_df$cluster[boot_min_df$cluster > 0L])])
  }

  boot_min_df$is_noise <- boot_min_df$cluster == 0L

  list(
    boot_min_df = boot_min_df,
    diagnostics = list(
      clustering_method = "pairwise_hungarian_graph",
      graph_clustering = if (leiden_used) "leiden" else "connected_components_fallback",
      leiden_error = leiden_error,
      n_run_pairs = choose(length(run_ids), 2),
      n_edges = nrow(edges),
      n_edges_kept = nrow(strong_edges),
      pairwise_leiden_gamma = pairwise_leiden_gamma,
      edge_weight_sigma = sigma,
      edge_distance_threshold = d_thr
    )
  )
}

cluster_bootstrap_minima <- function(boot_min_df, object, exclude_minor, min_barrier, clustering_method, minPts, pairwise_leiden_gamma = 0.01) {
  if (is.null(boot_min_df$node_id)) {
    boot_min_df$node_id <- seq_len(nrow(boot_min_df))
  }

  scaler <- get_bootstrap_distance_scaler(object)
  boot_min_df$x_scaled <- boot_min_df$x / scaler$dx
  boot_min_df$y_scaled <- boot_min_df$y / scaler$dy

  if (clustering_method == "pairwise_hungarian_graph") {
    out <- cluster_pairwise_hungarian_graph(
      boot_min_df = boot_min_df,
      pairwise_leiden_gamma = pairwise_leiden_gamma
    )
    out$diagnostics$distance_scale_dx <- scaler$dx
    out$diagnostics$distance_scale_dy <- scaler$dy
    return(out)
  }

  if (clustering_method == "hdbscan") {
    db <- dbscan::hdbscan(boot_min_df[, c("x_scaled", "y_scaled")], minPts = minPts)
    boot_min_df$cluster <- as.integer(db$cluster)
    boot_min_df$is_noise <- boot_min_df$cluster == 0L
    return(list(
      boot_min_df = boot_min_df,
      diagnostics = list(
        clustering_method = "hdbscan",
        distance_scale_dx = scaler$dx,
        distance_scale_dy = scaler$dy
      )
    ))
  }

  if (clustering_method == "hungarian") {
    ref_mins <- find_loc_min(object$original_ld, exclude_minor = exclude_minor, min_barrier = min_barrier)$mins
    if (exclude_minor && nrow(ref_mins) > 0) {
      ref_mins <- ref_mins[!ref_mins$is_minor, , drop = FALSE]
    }

    if (!nrow(ref_mins)) {
      boot_min_df$cluster <- 0L
      boot_min_df$is_noise <- TRUE
      return(list(
        boot_min_df = boot_min_df,
        diagnostics = list(clustering_method = "hungarian", n_reference = 0L)
      ))
    }

    ref_xy <- standardize_xy(ref_mins$x, ref_mins$y, scaler = scaler)
    assignments <- lapply(split(boot_min_df, boot_min_df$boot_index), function(df_run) {
      run_xy <- standardize_xy(df_run$x, df_run$y, scaler = scaler)
      assign_run_to_reference(run_xy, ref_xy)
    })
    boot_min_df$cluster <- as.integer(unlist(assignments, use.names = FALSE))
    boot_min_df$is_noise <- boot_min_df$cluster == 0L
    return(list(
      boot_min_df = boot_min_df,
      diagnostics = list(
        clustering_method = "hungarian",
        n_reference = nrow(ref_xy),
        distance_scale_dx = scaler$dx,
        distance_scale_dy = scaler$dy
      )
    ))
  }

  # gmm_bic
  max_g <- min(10L, nrow(boot_min_df))
  gmm_error <- NULL
  boot_xy <- as.matrix(boot_min_df[, c("x_scaled", "y_scaled")])
  g_range <- 1:max_g

  # mclust::Mclust() resolves mclustBIC via eval(..., parent.frame()),
  # so ensure the symbol is available in this caller frame.
  mclustBIC <- get("mclustBIC", envir = asNamespace("mclust"))

  fit <- tryCatch(
    mclust::Mclust(boot_xy, G = g_range, verbose = FALSE),
    error = function(e) {
      gmm_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(fit) || is.null(fit$classification)) {
    boot_min_df$cluster <- 0L
    boot_min_df$is_noise <- TRUE
    return(list(
      boot_min_df = boot_min_df,
      diagnostics = list(
        clustering_method = "gmm_bic",
        distance_scale_dx = scaler$dx,
        distance_scale_dy = scaler$dy,
        gmm_fit_failed = TRUE,
        gmm_error = gmm_error
      )
    ))
  }

  boot_min_df$cluster <- as.integer(fit$classification)
  boot_min_df$is_noise <- FALSE
  list(
    boot_min_df = boot_min_df,
    diagnostics = list(
      clustering_method = "gmm_bic",
      distance_scale_dx = scaler$dx,
      distance_scale_dy = scaler$dy,
      gmm_G = fit$G,
      gmm_model = fit$modelName,
      gmm_bic = fit$bic
    )
  )
}

assign_run_to_reference <- function(run_xy, ref_xy) {
  run_mat <- as.matrix(run_xy)
  n_run <- nrow(run_mat)
  n_ref <- nrow(ref_xy)

  if (!n_run || !n_ref) {
    return(integer(n_run))
  }

  cost <- outer(
    seq_len(n_run), seq_len(n_ref),
    FUN = function(i, j) {
      dx <- run_mat[i, 1] - ref_xy[j, 1]
      dy <- run_mat[i, 2] - ref_xy[j, 2]
      dx^2 + dy^2
    }
  )

  assigned <- integer(n_run)
  if (n_run <= n_ref) {
    match_cols <- clue::solve_LSAP(cost)
    assigned[] <- as.integer(match_cols)
    return(assigned)
  }

  # More run minima than references: match references to best unique run points.
  match_rows <- clue::solve_LSAP(t(cost))
  assigned[as.integer(match_rows)] <- seq_len(n_ref)
  assigned
}


#' Quick autoplot: points + 95% prediction & confidence ellipses
#' @export
#' @rdname summary.bootstrap_2d_ld
#' @param object A `"summary_bootstrap_2d_ld"` object produced by `summary()`.
#' @param mode Plot style. `"clusters"` shows cluster-colored points and ellipses;
#'   `"minima"` reproduces the minima-U jitter plot.
#' @param show_ellipses Logical; draw ellipses (spread of points). Default TRUE.
#' @param show_original_major_minima Logical; overlay major minima estimated
#'   from the original (non-bootstrap) sample. Default TRUE.
#' @param point_alpha Numeric; alpha for points. Default 0.35.
#' @param ellipse_alpha Numeric; alpha for filled prediction ellipse. Default 0.15.
#' @param ... Unused.
#' @export
autoplot.summary_bootstrap_2d_ld <- function(object,
                                             mode = c("clusters", "minima"),
                                             show_ellipses = TRUE,
                                             show_original_major_minima = TRUE,
                                             point_alpha = 0.35,
                                             ellipse_alpha = 0.15,
                                             ...) {
  if (!inherits(object, "summary_bootstrap_2d_ld")) {
    cli::cli_abort("{.arg object} must inherit from {.cls summary_bootstrap_2d_ld}.")
  }
  if (is.null(object$per_point)) {
    return(ggplot2::ggplot())
  }
  mode <- rlang::arg_match0(mode, c("clusters", "minima"))
  df_points <- object$per_point
  df_cl <- object$per_cluster

  if (mode == "minima") {
    x_range <- range(attr(object$original_ld$ss, "x_coords"))
    y_range <- range(attr(object$original_ld$ss, "y_coords"))
    x_jitter_amount <- (x_range[2] - x_range[1]) * 0.02
    y_jitter_amount <- (y_range[2] - y_range[1]) * 0.02

    boot_counts <- table(df_points$boot_index)
    boot_count_df <- data.frame(
      boot_index = as.integer(names(boot_counts)),
      count = as.integer(boot_counts)
    )

    plot_df <- df_points %>%
      dplyr::left_join(boot_count_df, by = "boot_index")

    max_count <- max(plot_df$count)
    if (max_count > 4) {
      plot_df <- plot_df %>%
        dplyr::mutate(count = ifelse(count >= 4, 4, count))
      cli::cli_inform("More than 4 minima found in some bootstrap samples.
                  Combining counts >=4 into a single category for plotting.")
    }

    return(
      ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = x + stats::runif(nrow(plot_df), -x_jitter_amount, x_jitter_amount),
          y = y + stats::runif(nrow(plot_df), -y_jitter_amount, y_jitter_amount),
          color = U,
          shape = as.factor(count)
        )
      ) +
        ggplot2::geom_point(size = 2, alpha = 0.7) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::labs(
          x = object$original_ld$vf$x,
          y = object$original_ld$vf$y,
          color = "U value",
          shape = "Number of minima\nin bootstrap sample"
        ) +
        ggplot2::theme_bw() +
        ggplot2::scale_shape_manual(values = c(16, 17, 15, 18, 3))
    )
  }

  p <- ggplot2::ggplot(df_points) +
    ggplot2::geom_point(
      ggplot2::aes(x = x, y = y, color = factor(cluster)),
      alpha = point_alpha,
      size = 1
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(color = "cluster")

  if (!is.null(df_cl) && nrow(df_cl)) {
    if (show_ellipses) {
      p <- p + ggforce::geom_ellipse(
        data = df_cl,
        ggplot2::aes(
          x0 = mean_x, y0 = mean_y,
          a = a_pred, b = b_pred, angle = angle
        ),
        color = "firebrick", fill = "firebrick", alpha = ellipse_alpha
      )
    }
  }

  if (isTRUE(show_original_major_minima) && !is.null(object$original_ld)) {
    min_barrier <- object$params$min_barrier
    if (is.null(min_barrier) || !is.finite(min_barrier)) {
      min_barrier <- 0.1
    }

    orig_major <- tryCatch(
      {
        mins <- find_loc_min(
          object$original_ld,
          exclude_minor = TRUE,
          min_barrier = min_barrier
        )$mins
        if (is.null(mins) || !nrow(mins)) {
          NULL
        } else {
          mins
        }
      },
      error = function(e) NULL
    )

    if (!is.null(orig_major)) {
      p <- p +
        ggplot2::geom_point(
          data = orig_major,
          ggplot2::aes(x = x, y = y),
          inherit.aes = FALSE,
          shape = 4,
          stroke = 1.1,
          size = 3,
          color = "black"
        )
    }
  }
  p + ggplot2::theme_bw()
}


# process_single_ld <- function(ld, p_func, find_func, ...) {
#   p_func()
#   find_func(ld, ...)
# }

