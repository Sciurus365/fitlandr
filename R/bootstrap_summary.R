summary.bootstrap_2d_ld <- function(object,
                                    exclude_minor = TRUE,
                                    min_barrier = 0.1,
                                    clustering_method = c("hungarian", "mean_potential", "hdbscan", "pairwise_hungarian_graph", "mean_potential_hessian", "gmm_bic"),
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
    c("hungarian", "mean_potential", "hdbscan", "pairwise_hungarian_graph", "mean_potential_hessian", "gmm_bic")
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


