# Shared machinery for clustering fitted objects. Future stream-function or
# vector-field clustering can reuse these helpers with another extractor.
prepare_object_clustering_input <- function(objects, extractor, object_type) {
  if (!is.list(objects) || !length(objects)) {
    cli::cli_abort("{.arg objects} must be a non-empty list of {object_type} objects.")
  }

  extracted <- lapply(objects, extractor)
  reference <- extracted[[1L]]
  for (i in seq_along(extracted)[-1L]) {
    current <- extracted[[i]]
    if (!identical(current$dimension, reference$dimension) ||
        !identical(dim(current$values), dim(reference$values)) ||
        !isTRUE(all.equal(current$grid, reference$grid, tolerance = sqrt(.Machine$double.eps)))) {
      cli::cli_abort("All {object_type} objects must use the same dimension and grid.")
    }
  }

  features <- do.call(rbind, lapply(extracted, function(x) as.numeric(x$values)))
  if (any(!is.finite(features))) {
    cli::cli_abort("The clustering representation must contain only finite values.")
  }

  list(
    features = features,
    extracted = extracted,
    reference = reference
  )
}

with_clustering_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  if (length(seed) != 1L || !is.finite(seed)) {
    cli::cli_abort("{.arg seed} must be NULL or one finite number.")
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed))
  force(code)
}

fit_kmeans_features <- function(features, k, nstart, iter.max) {
  n_objects <- nrow(features)
  if (length(k) != 1L || !is.finite(k) || k != as.integer(k) ||
      k < 1L || k > n_objects || (k > 1L && k >= n_objects)) {
    cli::cli_abort("{.arg k} must be 1 or an integer smaller than the number of objects ({n_objects}).")
  }
  if (length(nstart) != 1L || !is.finite(nstart) || nstart < 1L) {
    cli::cli_abort("{.arg nstart} must be a positive integer.")
  }
  if (length(iter.max) != 1L || !is.finite(iter.max) || iter.max < 1L) {
    cli::cli_abort("{.arg iter.max} must be a positive integer.")
  }

  tryCatch(
    stats::kmeans(
      features,
      centers = as.integer(k),
      nstart = as.integer(nstart),
      iter.max = as.integer(iter.max)
    ),
    error = function(e) {
      cli::cli_abort(c(
        "K-means clustering failed for {.code k = {k}}.",
        "i" = conditionMessage(e)
      ))
    }
  )
}

extract_landscape_clustering_values <- function(landscape) {
  if (!inherits(landscape, "landscape") || is.null(landscape$ss)) {
    cli::cli_abort("Every element of {.arg landscapes} must be a landscape object with a steady-state distribution.")
  }

  density <- landscape$ss
  if (any(!is.finite(density)) || any(density < 0) || !any(density > 0)) {
    cli::cli_abort("Every landscape steady-state distribution must be finite, non-negative, and contain positive mass.")
  }

  if (inherits(landscape, "1d_ld")) {
    grid <- landscape$dist$x
    dimension <- "1d"
  } else if (inherits(landscape, "2d_ld")) {
    grid <- list(
      x = sort(unique(landscape$dist$x)),
      y = sort(unique(landscape$dist$y))
    )
    dimension <- "2d"
  } else {
    cli::cli_abort("Landscape clustering currently supports {.cls 1d_ld} and {.cls 2d_ld} objects.")
  }

  list(values = density, grid = grid, dimension = dimension)
}

prepare_landscape_clustering_input <- function(landscapes) {
  prepared <- prepare_object_clustering_input(
    objects = landscapes,
    extractor = extract_landscape_clustering_values,
    object_type = "landscape"
  )
  prepared$landscapes <- landscapes
  prepared
}

density_to_potential <- function(density) {
  potential <- rep(NA_real_, length(density))
  potential[density > 0] <- -log(density[density > 0])
  dim(potential) <- dim(density)
  potential
}

build_landscape_cluster_center <- function(template, density, cluster_id) {
  potential <- density_to_potential(density)
  center <- template
  center$vf <- NULL
  center$cluster_id <- cluster_id

  if (inherits(template, "1d_ld")) {
    x_coords <- template$dist$x
    attr(density, "x_coords") <- x_coords
    center$ss <- density
    center$dist <- data.frame(
      x = x_coords,
      d = as.numeric(density),
      U = as.numeric(potential),
      U_plot = as.numeric(potential)
    )
    center$plot_2 <- ggplot2::ggplot(
      center$dist,
      ggplot2::aes(x = x, y = U_plot)
    ) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "x", y = "U", title = paste("Cluster", cluster_id))
    center$plot <- plotly::plot_ly(
      center$dist,
      x = ~x,
      y = ~U_plot,
      type = "scatter",
      mode = "lines"
    )
  } else {
    x_coords <- sort(unique(template$dist$x))
    y_coords <- sort(unique(template$dist$y))
    attr(density, "x_coords") <- x_coords
    attr(density, "y_coords") <- y_coords
    center$ss <- density
    center$dist <- expand.grid(x = x_coords, y = y_coords)
    center$dist$d <- as.numeric(density)
    center$dist$U <- as.numeric(potential)
    center$dist$U_plot <- as.numeric(potential)
    center$plot_2 <- ggplot2::ggplot(
      center$dist,
      ggplot2::aes(x = x, y = y)
    ) +
      ggplot2::geom_raster(ggplot2::aes(fill = U_plot)) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::coord_equal() +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "x", y = "y", fill = "U", title = paste("Cluster", cluster_id))
    center$plot <- plotly::plot_ly(
      data = center$dist,
      x = x_coords,
      y = y_coords,
      z = matrix(center$dist$U_plot, nrow = length(x_coords), ncol = length(y_coords)),
      type = "surface"
    )
  }

  class(center) <- unique(c("landscape_cluster_center", class(template)))
  center
}

#' Evaluate candidate numbers of landscape clusters
#'
#' Runs K-means clustering over several candidate values of `k`. Landscapes
#' are represented by their flattened steady-state distributions, rather than
#' their negative-log potential values. Use [autoplot()] on the result to draw
#' an elbow plot of within-cluster variance against `k`.
#'
#' @param landscapes A non-empty list of 1D or 2D landscape objects. All
#'   landscapes must have the same dimension and grid.
#' @param k_values Integer candidate numbers of clusters. `NULL` defaults to 1
#'   through the smaller of 10, the number of distinct landscapes, and one less
#'   than the number of landscapes. The one-landscape case evaluates only 1.
#' @param method Clustering method. Currently only `"kmeans"` is available.
#' @param nstart Number of random K-means initializations.
#' @param iter.max Maximum number of K-means iterations.
#' @param seed Optional random seed. The caller's random-number state is
#'   restored afterward.
#'
#' @return A `landscape_cluster_evaluation` object containing the candidate
#'   metrics and fitted K-means models.
#' @export
evaluate_landscape_clusters <- function(landscapes,
                                        k_values = NULL,
                                        method = c("kmeans"),
                                        nstart = 25L,
                                        iter.max = 100L,
                                        seed = NULL) {
  method <- match.arg(method)
  prepared <- prepare_landscape_clustering_input(landscapes)
  n_objects <- nrow(prepared$features)
  n_distinct <- nrow(unique(prepared$features))
  max_k <- if (n_objects == 1L) 1L else min(n_distinct, n_objects - 1L)
  if (is.null(k_values)) {
    k_values <- seq_len(min(10L, max_k))
  }
  if (!length(k_values) || any(!is.finite(k_values)) ||
      any(k_values != as.integer(k_values)) || any(k_values < 1L) ||
      any(k_values > max_k)) {
    cli::cli_abort(c(
      "{.arg k_values} must contain integers between 1 and {max_k}.",
      "i" = "The maximum is constrained by {n_distinct} distinct landscape{?s} and {n_objects} total landscape{?s}."
    ))
  }
  k_values <- sort(unique(as.integer(k_values)))

  models <- with_clustering_seed(seed, lapply(k_values, function(k) {
    fit_kmeans_features(prepared$features, k, nstart, iter.max)
  }))
  names(models) <- as.character(k_values)
  metrics <- data.frame(
    k = k_values,
    within_variance = vapply(models, function(x) x$tot.withinss, numeric(1)) / n_objects,
    between_variance = vapply(models, function(x) x$betweenss, numeric(1)) / n_objects,
    total_variance = vapply(models, function(x) x$totss, numeric(1)) / n_objects
  )
  metrics$explained_variance <- ifelse(
    metrics$total_variance > 0,
    metrics$between_variance / metrics$total_variance,
    0
  )

  structure(
    list(
      metrics = metrics,
      models = models,
      method = method,
      n_objects = n_objects,
      object_names = names(landscapes)
    ),
    class = "landscape_cluster_evaluation"
  )
}

#' Cluster potential landscapes
#'
#' Performs K-means clustering on steady-state distributions. Cluster centers
#' are arithmetic means in density space and are transformed back to potential
#' landscapes using `U = -log(density)` for the returned output.
#'
#' @inheritParams evaluate_landscape_clusters
#' @param k Number of clusters.
#'
#' @return A `landscape_clusters` object containing cluster assignments,
#'   cluster-center landscapes, density-space centers, and the underlying
#'   K-means result.
#' @export
cluster_landscapes <- function(landscapes,
                               k,
                               method = c("kmeans"),
                               nstart = 25L,
                               iter.max = 100L,
                               seed = NULL) {
  method <- match.arg(method)
  prepared <- prepare_landscape_clustering_input(landscapes)
  n_objects <- nrow(prepared$features)
  n_distinct <- nrow(unique(prepared$features))
  max_k <- if (n_objects == 1L) 1L else min(n_distinct, n_objects - 1L)
  if (length(k) == 1L && is.finite(k) && k > max_k) {
    cli::cli_abort(c(
      "{.arg k} cannot exceed {max_k} for these landscapes.",
      "i" = "The maximum is constrained by {n_distinct} distinct landscape{?s} and {n_objects} total landscape{?s}."
    ))
  }
  model <- with_clustering_seed(
    seed,
    fit_kmeans_features(prepared$features, k, nstart, iter.max)
  )

  center_densities <- lapply(seq_len(model$centers |> nrow()), function(i) {
    density <- model$centers[i, ]
    dim(density) <- dim(prepared$reference$values)
    density
  })
  centers <- lapply(seq_along(center_densities), function(i) {
    build_landscape_cluster_center(
      template = landscapes[[1L]],
      density = center_densities[[i]],
      cluster_id = i
    )
  })
  object_names <- names(landscapes)
  if (is.null(object_names)) {
    object_names <- rep.int(NA_character_, length(landscapes))
  }
  assignments <- data.frame(
    landscape_index = seq_along(landscapes),
    landscape_name = object_names,
    cluster = as.integer(model$cluster),
    distance = sqrt(rowSums(
      (prepared$features - model$centers[model$cluster, , drop = FALSE])^2
    ))
  )

  structure(
    list(
      assignments = assignments,
      cluster = as.integer(model$cluster),
      sizes = as.integer(model$size),
      centers = centers,
      center_densities = center_densities,
      kmeans = model,
      method = method,
      k = as.integer(k),
      landscapes = landscapes
    ),
    class = "landscape_clusters"
  )
}

#' Autoplot landscape-cluster evaluation
#'
#' @param object A `landscape_cluster_evaluation` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return A ggplot elbow plot.
#' @export
autoplot.landscape_cluster_evaluation <- function(object, ...) {
  ggplot2::ggplot(
    object$metrics,
    ggplot2::aes(x = k, y = within_variance)
  ) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = object$metrics$k) +
    ggplot2::labs(
      x = "Number of clusters (K)",
      y = "Within-cluster variance (mean squared distance)",
      title = "Landscape clustering elbow plot"
    ) +
    ggplot2::theme_bw()
}
