#' Bootstrap 1D vector fields
#'
#' @param vf A 1d_vectorfield or cv_vectorfield object.
#' @param block_length Length of each moving block.
#' @param n_boot Number of bootstrap samples.
#' @param seed Random seed.
#' @param ... Unused.
#'
#' @return A ootstrap_1d_vf object.
#' @export
bootstrap_1d_vf <- function(vf, block_length = NULL, n_boot = 200, seed = 1614, ...) {
  if (inherits(vf, "cv_vectorfield")) {
    vf <- vf$final_model
  }
  if (!inherits(vf, "1d_vectorfield")) {
    cli::cli_abort("Input {.arg vf} must be a {.cls 1d_vectorfield} or compatible {.cls cv_vectorfield} object.")
  }
  if (!identical(vf$method, "MVKE")) {
    cli::cli_abort("bootstrap_1d_vf() currently supports only {.code method = \"MVKE\"}.")
  }

  original_vectors <- vf$original_vectors
  original_vectors_normalized <- vf$original_vectors_normalized
  n_vec <- nrow(original_vectors_normalized)

  if (is.null(block_length)) {
    block_length <- ceiling(n_vec^(1 / 3))
  }

  blocks <- lapply(seq_len(n_vec - block_length + 1L), function(start_idx) start_idx:(start_idx + block_length - 1L))
  n_blocks <- length(blocks)
  n_blocks_needed <- ceiling(n_vec / block_length)
  h <- environment(vf[["MVKEresult"]])[["h"]]
  kernel <- environment(vf[["MVKEresult"]])[["kernel"]]
  dv <- vf$data_normalized
  lims <- vf$lims
  x_grid <- vf$vec_grid$x
  x_name <- vf$x
  n <- vf$n
  d_raw <- vf$data

  set.seed(seed)
  bootstrap_models <- lapply(seq_len(n_boot), function(b) {
    sampled_block_indices <- sample(seq_len(n_blocks), n_blocks_needed, replace = TRUE)
    sampled_indices <- unlist(blocks[sampled_block_indices])
    sampled_indices <- sampled_indices[sampled_indices <= n_vec]
    sampled_vectors <- original_vectors[sampled_indices, , drop = FALSE]
    sampled_vectors_normalized <- original_vectors_normalized[sampled_indices, , drop = FALSE]

    MVKEresult <- fitlandr::MVKE(
      d = as.matrix(sampled_vectors_normalized[, "x", drop = FALSE]),
      v = as.matrix(sampled_vectors_normalized[, "vx", drop = FALSE]),
      h = h,
      kernel = kernel
    )
    vx_grid <- vapply(x_grid, function(xx) {
      as.numeric(MVKEresult(normalize_x(xx, dv))$mu) %>% scale_up(dv)
    }, numeric(1))

    result <- list(
      vec_grid = data.frame(x = x_grid, vx = vx_grid, v_norm = abs(vx_grid)),
      VFCresult = NULL,
      MVKEresult = MVKEresult,
      data = d_raw,
      data_normalized = dv,
      original_vectors = sampled_vectors,
      original_vectors_normalized = sampled_vectors_normalized,
      x = x_name,
      y = NULL,
      lims = lims,
      n = n,
      method = "MVKE"
    )
    class(result) <- c("1d_vectorfield", "vectorfield")
    result
  })

  structure(list(
    bootstrap_models = bootstrap_models,
    original_vf = vf,
    block_length = block_length,
    n_boot = n_boot
  ), class = "bootstrap_1d_vf")
}

#' Bootstrap 1D landscapes
#'
#' @param boot_vf A `bootstrap_1d_vf` object.
#' @param ... Additional arguments passed to [make_1d_ld()].
#'
#' @return A `bootstrap_1d_ld` object.
#' @export
bootstrap_1d_ld <- function(boot_vf, ...) {
  if (!inherits(boot_vf, "bootstrap_1d_vf")) {
    cli::cli_abort("Input {.arg boot_vf} must be a {.cls bootstrap_1d_vf} object.")
  }

  original_ld <- purrr::quietly(make_1d_ld)(boot_vf$original_vf, ...)$result
  ref_x <- original_ld$dist$x
  boot_lds <- lapply(boot_vf$bootstrap_models, function(vf) {
    result <- purrr::quietly(make_1d_ld)(vf, ...)$result
    result$plot <- NULL
    result$plot_2 <- NULL
    if (!identical(result$dist$x, ref_x)) {
      cli::cli_abort("Bootstrap 1D landscapes must use exactly the same grid as the original landscape.")
    }
    result
  })

  structure(list(
    bootstrap_lds = boot_lds,
    original_ld = original_ld,
    n_boot = boot_vf$n_boot
  ), class = "bootstrap_1d_ld")
}

#' Summarize bootstrap minima for 1D landscapes
#'
#' @param object A `bootstrap_1d_ld` object.
#' @param exclude_minor Logical; exclude minor minima before clustering.
#' @param min_barrier Minimum barrier threshold used in minor-minimum
#'   detection.
#' @param clustering_method One of `"hungarian"` (default),
#'   `"mean_potential"`, `"pairwise_hungarian_graph"`, or `"gmm_bic"`.
#' @param level Coverage level for intervals.
#' @param one_per_run Logical; retain at most one point per run per cluster.
#' @param ... Unused.
#'
#' @return A `summary_bootstrap_1d_ld` object.
#' @method summary bootstrap_1d_ld
#' @export
summary.bootstrap_1d_ld <- function(object,
                                    exclude_minor = TRUE,
                                    min_barrier = 0.1,
                                    clustering_method = c("hungarian", "mean_potential", "pairwise_hungarian_graph", "gmm_bic"),
                                    level = 0.95,
                                    one_per_run = TRUE,
                                    ...) {
  if (!inherits(object, "bootstrap_1d_ld")) {
    cli::cli_abort("{.arg object} must be a {.cls bootstrap_1d_ld} object.")
  }

  clustering_method <- match.arg(clustering_method)

  boot_min_list <- lapply(seq_along(object$bootstrap_lds), function(i) {
    mins <- find_loc_min(object$bootstrap_lds[[i]], exclude_minor = exclude_minor, min_barrier = min_barrier)$mins
    if (exclude_minor && nrow(mins)) {
      mins <- mins[!mins$is_minor, , drop = FALSE]
    }
    if (!nrow(mins)) {
      return(NULL)
    }
    mins$boot_index <- i
    mins
  })
  boot_min_df <- dplyr::bind_rows(boot_min_list)
  per_boot <- tibble::tibble(boot_index = seq_len(object$n_boot)) |>
    dplyr::left_join(boot_min_df |>
      dplyr::count(boot_index, name = "n_mins"), by = "boot_index") |>
    dplyr::mutate(n_mins = dplyr::coalesce(.data$n_mins, 0L))

  if (!nrow(boot_min_df)) {
    out <- list(
      params = list(exclude_minor = exclude_minor, min_barrier = min_barrier, clustering_method = clustering_method, level = level, one_per_run = one_per_run),
      per_boot = per_boot,
      per_point = boot_min_df,
      per_cluster = NULL,
      diagnostics = list(clustering_method = clustering_method, n_reference = 0L),
      original_ld = object$original_ld,
      n_boot = object$n_boot
    )
    return(structure(out, class = "summary_bootstrap_1d_ld"))
  }

  x_vals <- object$original_ld$dist$x
  dx <- if (length(x_vals) > 1L) stats::median(diff(x_vals)) else 1
  if (!is.finite(dx) || dx <= 0) {
    dx <- 1
  }
  boot_min_df$node_id <- seq_len(nrow(boot_min_df))
  boot_min_df$x_scaled <- boot_min_df$x / dx

  if (clustering_method == "hungarian") {
    ref_mins <- find_loc_min(object$original_ld, exclude_minor = exclude_minor, min_barrier = min_barrier)$mins
    if (exclude_minor && nrow(ref_mins)) {
      ref_mins <- ref_mins[!ref_mins$is_minor, , drop = FALSE]
    }
    ref_x <- ref_mins$x / dx
    assignments <- lapply(split(boot_min_df, boot_min_df$boot_index), function(df_run) {
      assign_run_to_reference_1d(df_run$x / dx, ref_x)
    })
    boot_min_df$cluster <- as.integer(unlist(assignments, use.names = FALSE))
    boot_min_df$is_noise <- boot_min_df$cluster == 0L
    diagnostics <- list(clustering_method = "hungarian", n_reference = nrow(ref_mins), distance_scale_dx = dx)
  } else if (clustering_method == "mean_potential") {
    mean_ld <- build_mean_potential_ld_1d(object)
    ref_mins <- find_loc_min(mean_ld, exclude_minor = exclude_minor, min_barrier = min_barrier)$mins
    if (exclude_minor && nrow(ref_mins)) {
      ref_mins <- ref_mins[!ref_mins$is_minor, , drop = FALSE]
    }
    ref_x <- ref_mins$x / dx
    assignments <- lapply(split(boot_min_df, boot_min_df$boot_index), function(df_run) {
      assign_run_to_reference_1d(df_run$x / dx, ref_x)
    })
    boot_min_df$cluster <- as.integer(unlist(assignments, use.names = FALSE))
    boot_min_df$is_noise <- boot_min_df$cluster == 0L
    diagnostics <- list(clustering_method = "mean_potential", n_reference = nrow(ref_mins), distance_scale_dx = dx)
  } else if (clustering_method == "gmm_bic") {
    rlang::check_installed(
      "mclust",
      reason = "for {.code clustering_method = \"gmm_bic\"}. Install it with {.code install.packages(\"mclust\") }."
    )
    max_g <- min(10L, nrow(boot_min_df))
    fit <- tryCatch(mclust::Mclust(as.matrix(boot_min_df[, "x_scaled", drop = FALSE]), G = 1:max_g, verbose = FALSE), error = function(e) NULL)
    if (is.null(fit) || is.null(fit$classification)) {
      boot_min_df$cluster <- 0L
      boot_min_df$is_noise <- TRUE
      diagnostics <- list(clustering_method = "gmm_bic", distance_scale_dx = dx, gmm_fit_failed = TRUE)
    } else {
      boot_min_df$cluster <- as.integer(fit$classification)
      boot_min_df$is_noise <- FALSE
      diagnostics <- list(clustering_method = "gmm_bic", distance_scale_dx = dx, gmm_G = fit$G, gmm_model = fit$modelName)
    }
  } else {
    out <- cluster_pairwise_hungarian_graph_1d(boot_min_df, x_scale = dx)
    boot_min_df <- out$boot_min_df
    diagnostics <- out$diagnostics
  }

  df_c <- boot_min_df |>
    dplyr::filter(!is_noise, cluster > 0L)
  if (one_per_run && nrow(df_c)) {
    centers0 <- df_c |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(cx = mean(x), .groups = "drop")
    df_c <- df_c |>
      dplyr::left_join(centers0, by = "cluster") |>
      dplyr::mutate(d2 = (x - cx)^2) |>
      dplyr::group_by(cluster, boot_index) |>
      dplyr::slice_min(order_by = d2, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::select(-cx, -d2)
  }

  if (!nrow(df_c)) {
    out <- list(
      params = list(exclude_minor = exclude_minor, min_barrier = min_barrier, clustering_method = clustering_method, level = level, one_per_run = one_per_run),
      per_boot = per_boot,
      per_point = boot_min_df,
      per_cluster = NULL,
      diagnostics = diagnostics,
      original_ld = object$original_ld,
      n_boot = object$n_boot
    )
    return(structure(out, class = "summary_bootstrap_1d_ld"))
  }

  z_val <- stats::qnorm((1 + level) / 2)
  per_cluster <- df_c |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      n_points = dplyr::n(),
      n_runs = dplyr::n_distinct(boot_index),
      mean_x = mean(x),
      mean_U = mean(U),
      sd_x = stats::sd(x),
      sd_U = stats::sd(U),
      CI_U_lower = stats::quantile(U, probs = (1 - level) / 2),
      CI_U_upper = stats::quantile(U, probs = 1 - (1 - level) / 2),
      var_x = stats::var(x),
      stability = n_runs / object$n_boot,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      x_pred_lower = mean_x - z_val * sd_x,
      x_pred_upper = mean_x + z_val * sd_x,
      x_conf_lower = mean_x - z_val * sd_x / sqrt(pmax(n_runs, 1L)),
      x_conf_upper = mean_x + z_val * sd_x / sqrt(pmax(n_runs, 1L))
    )

  out <- list(
    params = list(exclude_minor = exclude_minor, min_barrier = min_barrier, clustering_method = clustering_method, level = level, one_per_run = one_per_run),
    per_boot = per_boot,
    per_point = boot_min_df,
    per_cluster = per_cluster,
    diagnostics = diagnostics,
    original_ld = object$original_ld,
    n_boot = object$n_boot
  )
  structure(out, class = "summary_bootstrap_1d_ld")
}

#' Autoplot a 1D bootstrap summary
#'
#' @param object A `summary_bootstrap_1d_ld` object.
#' @param mode One of `"clusters"` or `"minima"`.
#' @param show_intervals Logical; draw prediction/confidence intervals.
#' @param show_original_major_minima Logical; overlay original major minima.
#' @param point_alpha Numeric alpha for points.
#' @param ... Unused.
#'
#' @return A ggplot object.
#' @export
autoplot.summary_bootstrap_1d_ld <- function(object,
                                             mode = c("clusters", "minima"),
                                             show_intervals = TRUE,
                                             show_original_major_minima = TRUE,
                                             point_alpha = 0.35,
                                             ...) {
  if (!inherits(object, "summary_bootstrap_1d_ld")) {
    cli::cli_abort("{.arg object} must inherit from {.cls summary_bootstrap_1d_ld}.")
  }
  mode <- rlang::arg_match0(mode, c("clusters", "minima"))
  df_points <- object$per_point
  df_cl <- object$per_cluster

  if (mode == "minima") {
    boot_counts <- table(df_points$boot_index)
    boot_count_df <- data.frame(boot_index = as.integer(names(boot_counts)), count = as.integer(boot_counts))
    plot_df <- df_points |>
      dplyr::left_join(boot_count_df, by = "boot_index")

    return(
      ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = U, color = U, shape = as.factor(count))) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::theme_bw() +
        ggplot2::labs(x = object$original_ld$vf$x, y = "U", color = "U value", shape = "N minima")
    )
  }

  p <- ggplot2::ggplot(df_points) +
    ggplot2::geom_point(
      ggplot2::aes(x = x, y = U, color = factor(cluster)),
      alpha = point_alpha,
      size = 1
    ) +
    ggplot2::labs(x = object$original_ld$vf$x, y = "U", color = "cluster") +
    ggplot2::theme_bw()

  if (show_intervals && !is.null(df_cl) && nrow(df_cl)) {
    p <- p +
      ggplot2::geom_segment(
        data = df_cl,
        ggplot2::aes(x = x_pred_lower, xend = x_pred_upper, y = mean_U, yend = mean_U),
        inherit.aes = FALSE,
        alpha = 0.2,
        linewidth = 2,
        color = "firebrick"
      ) +
      ggplot2::geom_segment(
        data = df_cl,
        ggplot2::aes(x = x_conf_lower, xend = x_conf_upper, y = mean_U, yend = mean_U),
        inherit.aes = FALSE,
        linewidth = 1,
        color = "black"
      )
  }

  if (isTRUE(show_original_major_minima) && !is.null(object$original_ld)) {
    orig_major <- tryCatch({
      mins <- find_loc_min(object$original_ld, exclude_minor = TRUE, min_barrier = object$params$min_barrier)$mins
      if (!is.null(mins) && nrow(mins) && "is_minor" %in% names(mins)) {
        mins <- mins[!mins$is_minor, , drop = FALSE]
      }
      mins
    }, error = function(e) NULL)

    if (!is.null(orig_major) && nrow(orig_major)) {
      p <- p + ggplot2::geom_point(
        data = orig_major,
        ggplot2::aes(x = x, y = U),
        inherit.aes = FALSE,
        shape = 4,
        stroke = 1.1,
        size = 3,
        color = "black"
      )
    }
  }

  p
}
