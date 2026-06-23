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

  g <- igraph::graph_from_data_frame(
    d = strong_edges[, c("node_u", "node_v", "weight")],
    directed = FALSE,
    vertices = data.frame(name = as.character(seq_len(n_points)))
  )

  lc <- igraph::cluster_leiden(
    g,
    weights = igraph::E(g)$weight,
    objective_function = "CPM",
    resolution = pairwise_leiden_gamma
  )
  cluster_vec <- as.integer(igraph::membership(lc))
  if (length(cluster_vec) != n_points) {
    cli::cli_abort("Leiden clustering returned an unexpected number of memberships.")
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
    remap <- stats::setNames(seq_along(uniq), uniq)
    boot_min_df$cluster[boot_min_df$cluster > 0L] <- as.integer(remap[as.character(boot_min_df$cluster[boot_min_df$cluster > 0L])])
  }

  boot_min_df$is_noise <- boot_min_df$cluster == 0L

  list(
    boot_min_df = boot_min_df,
    diagnostics = list(
      clustering_method = "pairwise_hungarian_graph",
      graph_clustering = "leiden",
      n_run_pairs = choose(length(run_ids), 2),
      n_edges = nrow(edges),
      n_edges_kept = nrow(strong_edges),
      pairwise_leiden_gamma = pairwise_leiden_gamma,
      edge_weight_sigma = sigma,
      edge_distance_threshold = d_thr
    )
  )
}

cluster_bootstrap_minima <- function(boot_min_df, object, exclude_minor, min_barrier_fraction, min_convex_hull_range_fraction, clustering_method, minPts, pairwise_leiden_gamma = 0.01) {
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
    rlang::check_installed(
      "dbscan",
      reason = "for {.code clustering_method = \"hdbscan\"}. Install it with {.code install.packages(\"dbscan\") }."
    )
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
    ref_mins <- find_loc_min(
      object$original_ld,
      exclude_minor = exclude_minor,
      min_barrier_fraction = min_barrier_fraction,
      min_convex_hull_range_fraction = min_convex_hull_range_fraction
    )$mins
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
  rlang::check_installed(
    "mclust",
    reason = "for {.code clustering_method = \"gmm_bic\"}. Install it with {.code install.packages(\"mclust\") }."
  )
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
