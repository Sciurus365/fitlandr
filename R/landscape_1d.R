ss_fp_1d <- function(vf, linear_interp = TRUE, n_grid = 200L) {
  if (inherits(vf, "cv_vectorfield")) {
    vf <- vf$final_model
  } else if (!inherits(vf, "vectorfield")) {
    cli::cli_abort("Input {.arg vf} must be a {.cls vectorfield} or {.cls cv_vectorfield} object.")
  }

  if (!inherits(vf, "1d_vectorfield")) {
    cli::cli_abort("{.arg vf} must be a {.cls 1d_vectorfield} object.")
  }

  x_range <- vf$lims[1:2]
  x_coords <- seq(x_range[1], x_range[2], length.out = n_grid)
  mu <- numeric(n_grid)
  a <- numeric(n_grid)

  for (i in seq_along(x_coords)) {
    pred <- stats::predict(vf, x_coords[i], linear_interp = linear_interp)
    mu[i] <- as.numeric(pred$v)
    a[i] <- as.numeric(pred$a[1, 1])
  }

  bad_a <- !is.finite(a) | a <= 0
  if (all(bad_a)) {
    cli::cli_abort("All diffusion estimates are non-finite or non-positive on the 1D grid.")
  }
  if (any(bad_a)) {
    a[bad_a] <- min(a[!bad_a], na.rm = TRUE) * 0.1
    cli::cli_warn("Some 1D diffusion estimates were non-finite or non-positive and were replaced by a small positive fallback.")
  }

  integrand <- 2 * mu / a
  cum_int <- numeric(n_grid)
  if (n_grid > 1L) {
    dx <- diff(x_coords)
    for (i in 2:n_grid) {
      cum_int[i] <- cum_int[i - 1L] + 0.5 * (integrand[i - 1L] + integrand[i]) * dx[i - 1L]
    }
  }

  log_rho <- cum_int - log(a)
  log_rho <- log_rho - max(log_rho, na.rm = TRUE)
  rho <- exp(log_rho)
  rho <- rho / sum(rho, na.rm = TRUE)

  attr(rho, "x_coords") <- x_coords
  rho
}

#' Create a 1D potential landscape from a 1D vector field
#'
#' This is the 1D analogue of [make_2d_ld()], using the drift and diffusion
#' estimates from a fitted 1D vector field.
#'
#' @param vf A `1d_vectorfield` or `cv_vectorfield` object.
#' @param linear_interp Logical indicating whether to use interpolation in
#'   predictions.
#' @param n_grid Number of grid points.
#'
#' @return A `1d_static_ld` object.
#' @export
make_1d_ld <- function(vf, linear_interp = TRUE, n_grid = 200L) {
  ss <- ss_fp_1d(vf, linear_interp = linear_interp, n_grid = n_grid)
  x_coords <- attr(ss, "x_coords")

  U <- rep(NA_real_, length(ss))
  U[ss > 0] <- -log(ss[ss > 0])
  U_plot <- U
  U_plot[!is.finite(U_plot)] <- NA_real_

  dist <- data.frame(
    x = x_coords,
    d = as.numeric(ss),
    U = as.numeric(U),
    U_plot = as.numeric(U_plot)
  )

  plot <- plotly::plot_ly(dist, x = ~x, y = ~U_plot, type = "scatter", mode = "lines") %>%
    plotly::layout(xaxis = list(title = if (inherits(vf, "cv_vectorfield")) vf$final_model$x else vf$x),
                   yaxis = list(title = "U"))

  plot_2 <- ggplot2::ggplot(dist, ggplot2::aes(x = x, y = U_plot)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = if (inherits(vf, "cv_vectorfield")) vf$final_model$x else vf$x, y = "U")

  structure(list(
    dist = dist,
    plot = plot,
    plot_2 = plot_2,
    vf = vf,
    ss = ss
  ), class = c("1d_static_ld", "1d_ld", "landscape"))
}

find_loc_min_1d <- function(ld, exclude_minor = TRUE, min_barrier = 0.1) {
  dist <- ld$dist
  local_minima <- which(diff(sign(diff(dist$U))) == 2) + 1L
  local_maxima <- which(diff(sign(diff(dist$U))) == -2) + 1L

  mins <- data.frame(
    x = dist$x[local_minima],
    U = dist$U[local_minima]
  )

  n_mins <- nrow(mins)
  all_barriers <- matrix(NA_real_, nrow = n_mins, ncol = n_mins)
  minor_mins <- integer(0)

  if (n_mins > 1L && exclude_minor) {
    for (i in seq_len(n_mins - 1L)) {
      for (j in (i + 1L):n_mins) {
        idx_i <- local_minima[i]
        idx_j <- local_minima[j]
        between <- seq.int(idx_i, idx_j)
        saddle_u <- max(dist$U[between], na.rm = TRUE)
        all_barriers[i, j] <- saddle_u - mins$U[i]
        all_barriers[j, i] <- saddle_u - mins$U[j]
      }
    }

    finite_barriers <- all_barriers[is.finite(all_barriers)]
    max_barrier <- if (length(finite_barriers)) max(finite_barriers) else Inf
    all_barriers_copy <- all_barriers

    repeat {
      current_barriers <- all_barriers_copy[is.finite(all_barriers_copy)]
      if (!length(current_barriers)) {
        break
      }
      current_min_barrier <- min(current_barriers)
      if (is.infinite(current_min_barrier) || current_min_barrier >= min_barrier * max_barrier) {
        break
      }
      locs <- which(all_barriers_copy == current_min_barrier, arr.ind = TRUE)
      min_to_remove <- locs[1, 1]
      minor_mins <- unique(c(minor_mins, min_to_remove))
      all_barriers_copy[min_to_remove, ] <- NA_real_
      all_barriers_copy[, min_to_remove] <- NA_real_
    }
  }

  mins$is_minor <- seq_len(n_mins) %in% minor_mins
  structure(list(mins = mins, barriers = all_barriers, maxima_index = local_maxima), class = "ld_min")
}

build_mean_potential_ld_1d <- function(object) {
  ref_ld <- object$original_ld
  x_vals <- ref_ld$dist$x
  if (!length(object$bootstrap_lds)) {
    cli::cli_abort("bootstrap_1d_ld object contains no bootstrap landscapes.")
  }

  u_mats <- lapply(object$bootstrap_lds, function(ld) {
    if (!identical(ld$dist$x, x_vals)) {
      cli::cli_abort("All 1D bootstrap landscapes must use the same grid.")
    }
    ld$dist$U
  })

  mean_u <- Reduce(`+`, u_mats) / length(u_mats)
  mean_d <- exp(-(mean_u - min(mean_u, na.rm = TRUE)))
  mean_d <- mean_d / sum(mean_d, na.rm = TRUE)

  ref_ld$dist <- data.frame(x = x_vals, d = mean_d, U = mean_u, U_plot = mean_u)
  ref_ld$ss <- mean_d
  attr(ref_ld$ss, "x_coords") <- x_vals
  ref_ld$plot <- NULL
  ref_ld$plot_2 <- NULL
  ref_ld
}

assign_run_to_reference_1d <- function(run_x, ref_x) {
  n_run <- length(run_x)
  n_ref <- length(ref_x)
  if (!n_run || !n_ref) {
    return(integer(n_run))
  }

  cost <- outer(seq_len(n_run), seq_len(n_ref), FUN = function(i, j) (run_x[i] - ref_x[j])^2)
  assigned <- integer(n_run)
  if (n_run <= n_ref) {
    assigned[] <- as.integer(clue::solve_LSAP(cost))
    return(assigned)
  }

  match_rows <- clue::solve_LSAP(t(cost))
  assigned[as.integer(match_rows)] <- seq_len(n_ref)
  assigned
}

cluster_pairwise_hungarian_graph_1d <- function(boot_min_df, x_scale) {
  n_points <- nrow(boot_min_df)
  if (n_points < 2L) {
    boot_min_df$cluster <- 0L
    boot_min_df$is_noise <- TRUE
    return(list(boot_min_df = boot_min_df, diagnostics = list(clustering_method = "pairwise_hungarian_graph", n_edges = 0L)))
  }

  by_run <- split(boot_min_df, boot_min_df$boot_index)
  run_ids <- names(by_run)
  edge_rows <- list()
  edge_counter <- 0L

  for (a in seq_len(length(run_ids) - 1L)) {
    for (b in (a + 1L):length(run_ids)) {
      df_a <- by_run[[a]]
      df_b <- by_run[[b]]
      cost <- outer(df_a$x_scaled, df_b$x_scaled, FUN = function(i, j) abs(i - j))
      dvals <- as.numeric(cost)
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
    return(list(boot_min_df = boot_min_df, diagnostics = list(clustering_method = "pairwise_hungarian_graph", n_edges = 0L)))
  }

  edges <- dplyr::bind_rows(edge_rows)
  d_thr <- as.numeric(stats::quantile(edges$d, probs = 0.75, na.rm = TRUE))
  strong_edges <- edges[edges$d <= d_thr, , drop = FALSE]
  if (!nrow(strong_edges)) {
    strong_edges <- edges
  }

  cluster_vec <- connected_components_from_edges(
    n_nodes = n_points,
    edge_u = strong_edges$node_u,
    edge_v = strong_edges$node_v
  )

  boot_min_df$cluster <- as.integer(cluster_vec[boot_min_df$node_id])
  valid <- boot_min_df |>
    dplyr::filter(cluster > 0L) |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(n_runs = dplyr::n_distinct(boot_index), .groups = "drop") |>
    dplyr::filter(n_runs >= 2L)
  valid_clusters <- valid$cluster
  boot_min_df$cluster[!(boot_min_df$cluster %in% valid_clusters)] <- 0L
  uniq <- sort(unique(boot_min_df$cluster[boot_min_df$cluster > 0L]))
  if (length(uniq)) {
    remap <- stats::setNames(seq_along(uniq), uniq)
    boot_min_df$cluster[boot_min_df$cluster > 0L] <- as.integer(remap[as.character(boot_min_df$cluster[boot_min_df$cluster > 0L])])
  }
  boot_min_df$is_noise <- boot_min_df$cluster == 0L

  list(
    boot_min_df = boot_min_df,
    diagnostics = list(
      clustering_method = "pairwise_hungarian_graph",
      distance_scale_dx = x_scale,
      n_edges = nrow(edges),
      n_edges_kept = nrow(strong_edges)
    )
  )
}

