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

  p <- progressr::progressor(steps = boot_vf$n_boot)

  boot_lds <- lapply(
    boot_vf$bootstrap_models,
    function(vf) {
      p()
      result <- purrr::quietly(make_2d_ld)(vf, ...)$result
      # to save space:
      result$plot <- NULL
      result$plot_2 <- NULL

      return(result)
    }
  )

  return(structure(list(
    bootstrap_lds = boot_lds,
    original_ld = make_2d_ld(boot_vf$original_vf),
    n_boot = boot_vf$n_boot
  ), class = "bootstrap_2d_ld"))
}


#' @rdname bootstrap_2d_ld
#'
#' @param object A `bootstrap_2d_ld` object with fields:
#'   - `bootstrap_lds`: list of landscapes,
#'   - `original_ld`: the original landscape,
#'   - `n_boot`: number of bootstrap runs.
#' @param exclude_minor Logical; exclude minor local minima. Default TRUE.
#' @param clustering_method Clustering backend for pooled minima. One of
#'   `"hdbscan"` (default), `"hungarian"`, `"gmm_bic"`.
#' @param minPts Integer; HDBSCAN minPts. Default 5.
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
#' Ellipse parameters are derived from the eigen-decomposition of the 2×2 covariance
#' matrix per cluster. Semi-axes for the **prediction** ellipse are
#' \eqn{\sqrt{\lambda_i\,\chi^2_{2,\alpha}}}; for the **confidence** ellipse they are
#' \eqn{\sqrt{\lambda_i\,\chi^2_{2,\alpha}/n_{\text{runs}}}} where \eqn{n_{\text{runs}}}
#' is the number of bootstrap runs contributing to that cluster (uncertainty of the mean).
#' Ellipses are drawn with `ggforce::geom_ellipse()` which expects aesthetics
#' `x0`, `y0`, `a`, `b`, `angle`.  \[See ggforce docs.\]  # (geom_ellipse API)  [1](https://rstudio-pubs-static.s3.amazonaws.com/1236382_7016d680936d411e8fd45fc0b8b62b0c.html)
#' The \eqn{\chi^2_{2,\alpha}} quantile comes from the chi-square distribution
#' (e.g., 0.95 → 5.991).  \[See NIST/ITL table.\]  # (chi-square critical)  [2](https://stackoverflow.com/questions/70010774/dbscan-choice-of-epsilon-through-elbow-method)
#'
#' @export
#' @method summary bootstrap_2d_ld
summary.bootstrap_2d_ld <- function(object,
                                    exclude_minor = TRUE,
                                    clustering_method = c("hdbscan", "hungarian", "gmm_bic"),
                                    minPts = 5,
                                    level = 0.95,
                                    one_per_run = TRUE,
                                    ...) {
  if (!is.list(object) || is.null(object$bootstrap_lds) || is.null(object$n_boot)) {
    cli::cli_abort("{.arg object} must be a {.cls bootstrap_2d_ld} object with {.field bootstrap_lds} and {.field n_boot}.")
  }
  clustering_method <- rlang::arg_match0(clustering_method, c("hdbscan", "hungarian", "gmm_bic"))

  # ---------- 1) Collect minima across bootstrap runs ----------
  p <- progressr::progressor(steps = length(object$bootstrap_lds))
  lds <- object$bootstrap_lds

  boot_mins <- lapply(lds, function(ld) {
    p()
    find_loc_min(ld, exclude_minor = exclude_minor)
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
        exclude_minor = exclude_minor, clustering_method = clustering_method,
        minPts = minPts, level = level,
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

  # ---------- 2) Cluster pooled minima ----------
  cl_out <- cluster_bootstrap_minima(
    boot_min_df = boot_min_df,
    object = object,
    exclude_minor = exclude_minor,
    clustering_method = clustering_method,
    minPts = minPts
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
        exclude_minor = exclude_minor, clustering_method = clustering_method,
        minPts = minPts, level = level,
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
      exclude_minor = exclude_minor, clustering_method = clustering_method,
      minPts = minPts, level = level,
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

cluster_bootstrap_minima <- function(boot_min_df, object, exclude_minor, clustering_method, minPts) {
  if (clustering_method == "hdbscan") {
    db <- dbscan::hdbscan(boot_min_df[, c("x", "y")], minPts = minPts)
    boot_min_df$cluster <- as.integer(db$cluster)
    boot_min_df$is_noise <- boot_min_df$cluster == 0L
    return(list(
      boot_min_df = boot_min_df,
      diagnostics = list(clustering_method = "hdbscan")
    ))
  }

  if (clustering_method == "hungarian") {
    ref_mins <- find_loc_min(object$original_ld, exclude_minor = exclude_minor)$mins
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

    ref_xy <- as.matrix(ref_mins[, c("x", "y"), drop = FALSE])
    assignments <- lapply(split(boot_min_df, boot_min_df$boot_index), function(df_run) {
      assign_run_to_reference(df_run[, c("x", "y"), drop = FALSE], ref_xy)
    })
    boot_min_df$cluster <- as.integer(unlist(assignments, use.names = FALSE))
    boot_min_df$is_noise <- boot_min_df$cluster == 0L
    return(list(
      boot_min_df = boot_min_df,
      diagnostics = list(clustering_method = "hungarian", n_reference = nrow(ref_xy))
    ))
  }

  # gmm_bic
  max_g <- min(10L, nrow(boot_min_df))
  gmm_error <- NULL
  boot_xy <- as.matrix(boot_min_df[, c("x", "y")])
  g_range <- 1:max_g
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
#' @param point_alpha Numeric; alpha for points. Default 0.35.
#' @param ellipse_alpha Numeric; alpha for filled prediction ellipse. Default 0.15.
#' @param ... Unused.
#' @export
autoplot.summary_bootstrap_2d_ld <- function(object,
                                             mode = c("clusters", "minima"),
                                             show_ellipses = TRUE,
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
  p + ggplot2::theme_bw()
}


# process_single_ld <- function(ld, p_func, find_func, ...) {
#   p_func()
#   find_func(ld, ...)
# }

