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
                                             min_barrier_fraction,
                                             min_convex_hull_range_fraction,
                                             level,
                                             one_per_run) {
  scaler <- get_bootstrap_distance_scaler(object)

  mean_ld <- build_mean_potential_ld(object)
  ref_mins <- find_loc_min(
    mean_ld,
    exclude_minor = exclude_minor,
    min_barrier_fraction = min_barrier_fraction,
    min_convex_hull_range_fraction = min_convex_hull_range_fraction
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
        min_barrier_fraction = min_barrier_fraction,
        min_convex_hull_range_fraction = min_convex_hull_range_fraction,
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
      min_barrier_fraction = min_barrier_fraction,
      min_convex_hull_range_fraction = min_convex_hull_range_fraction,
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
                                     min_barrier_fraction,
                                     min_convex_hull_range_fraction,
                                     level,
                                     one_per_run) {
  scaler <- get_bootstrap_distance_scaler(object)

  mean_ld <- build_mean_potential_ld(object)
  ref_mins <- find_loc_min(
    mean_ld,
    exclude_minor = exclude_minor,
    min_barrier_fraction = min_barrier_fraction,
    min_convex_hull_range_fraction = min_convex_hull_range_fraction
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
        min_barrier_fraction = min_barrier_fraction,
        min_convex_hull_range_fraction = min_convex_hull_range_fraction,
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
      min_barrier_fraction = min_barrier_fraction,
      min_convex_hull_range_fraction = min_convex_hull_range_fraction,
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
#' @param min_barrier_fraction Minimum barrier height fraction, relative to the
#'   highest barrier, used by `find_loc_min()` to classify minor minima.
#'   Default 0.1.
#' @param min_convex_hull_range_fraction Minimum barrier height fraction,
#'   relative to the potential range inside the observed-data convex hull,
#'   used by `find_loc_min()` to classify minor minima. Default 0.01.
#' @param clustering_method Clustering backend for pooled minima. One of
#'   `"hungarian"` (default), `"mean_potential"`, `"hdbscan"`,
#'   `"pairwise_hungarian_graph"`, `"mean_potential_hessian"`, or `"gmm_bic"`.
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
