extract_stream_clustering_values <- function(stream) {
  if (!inherits(stream, "2d_stream") || is.null(stream$grid)) {
    cli::cli_abort("Every element of {.arg streams} must be a {.cls 2d_stream} object with a stream-function grid.")
  }
  grid <- stream$grid
  required <- c("x", "y", "A")
  if (!all(required %in% names(grid)) || any(!is.finite(as.matrix(grid[required])))) {
    cli::cli_abort("Every stream grid must contain finite {.field x}, {.field y}, and {.field A} values.")
  }

  x_coords <- sort(unique(grid$x))
  y_coords <- sort(unique(grid$y))
  nx <- length(x_coords)
  ny <- length(y_coords)
  if (nrow(grid) != nx * ny || anyDuplicated(grid[c("x", "y")])) {
    cli::cli_abort("Every stream function must use a complete rectangular grid with one value per point.")
  }

  grid_index <- match(grid$x, x_coords) + (match(grid$y, y_coords) - 1L) * nx
  values <- numeric(nx * ny)
  values[grid_index] <- grid$A
  values <- values - mean(values)
  dim(values) <- c(nx, ny)

  list(
    values = values,
    grid = list(x = x_coords, y = y_coords),
    dimension = "2d_stream"
  )
}

prepare_stream_clustering_input <- function(streams) {
  if (inherits(streams, "group_dynamics")) {
    streams <- streams$streams
  }
  prepared <- prepare_object_clustering_input(
    objects = streams,
    extractor = extract_stream_clustering_values,
    object_type = "stream-function"
  )
  prepared$streams <- streams
  prepared
}

feature_cluster_limits <- function(features) {
  n_objects <- nrow(features)
  n_distinct <- nrow(unique(features))
  max_k <- if (n_objects == 1L) 1L else min(n_distinct, n_objects - 1L)
  list(n_objects = n_objects, n_distinct = n_distinct, max_k = max_k)
}

build_stream_cluster_center <- function(template, values, cluster_id) {
  x_coords <- sort(unique(template$grid$x))
  y_coords <- sort(unique(template$grid$y))
  grid <- expand.grid(x = x_coords, y = y_coords)
  grid$A <- as.numeric(values)
  x_label <- if (!is.null(template$pf$x)) template$pf$x else "x"
  y_label <- if (!is.null(template$pf$y)) template$pf$y else "y"

  structure(
    list(
      grid = grid,
      fitted_grid = NULL,
      residual_grid = NULL,
      rmse = NA_real_,
      relative_error = NA_real_,
      pf = list(x = x_label, y = y_label),
      cluster_id = cluster_id
    ),
    class = c("stream_cluster_center", "2d_stream")
  )
}

#' Evaluate candidate numbers of stream-function clusters
#'
#' Runs K-means over several candidate cluster counts. Each stream function is
#' represented by its common-grid `A` values after subtracting its grid mean,
#' because a stream function is identifiable only up to an additive constant.
#'
#' @param streams A non-empty list of `2d_stream` objects, or a
#'   `group_dynamics` object. All streams must use the same grid.
#' @param k_values Integer candidate cluster counts. `NULL` chooses the valid
#'   range up to 10 based on distinct streams and sample size.
#' @param method Clustering method. Currently only `"kmeans"` is available.
#' @param nstart Number of K-means random initializations.
#' @param iter.max Maximum number of K-means iterations.
#' @param seed Optional random seed. The caller's random-number state is
#'   restored afterward.
#'
#' @return A `stream_cluster_evaluation` object containing candidate metrics
#'   and fitted K-means models.
#' @export
evaluate_stream_clusters <- function(streams,
                                     k_values = NULL,
                                     method = c("kmeans"),
                                     nstart = 25L,
                                     iter.max = 100L,
                                     seed = NULL) {
  method <- match.arg(method)
  prepared <- prepare_stream_clustering_input(streams)
  streams <- prepared$streams
  limits <- feature_cluster_limits(prepared$features)
  if (is.null(k_values)) {
    k_values <- seq_len(min(10L, limits$max_k))
  }
  if (!length(k_values) || any(!is.finite(k_values)) ||
      any(k_values != as.integer(k_values)) || any(k_values < 1L) ||
      any(k_values > limits$max_k)) {
    cli::cli_abort(c(
      "{.arg k_values} must contain integers between 1 and {limits$max_k}.",
      "i" = "The maximum is constrained by {limits$n_distinct} distinct stream function{?s} and {limits$n_objects} total stream function{?s}."
    ))
  }
  k_values <- sort(unique(as.integer(k_values)))

  models <- with_clustering_seed(seed, lapply(k_values, function(k) {
    fit_kmeans_features(prepared$features, k, nstart, iter.max)
  }))
  names(models) <- as.character(k_values)
  metrics <- data.frame(
    k = k_values,
    within_variance = vapply(models, function(x) x$tot.withinss, numeric(1)) / limits$n_objects,
    between_variance = vapply(models, function(x) x$betweenss, numeric(1)) / limits$n_objects,
    total_variance = vapply(models, function(x) x$totss, numeric(1)) / limits$n_objects
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
      n_objects = limits$n_objects,
      object_names = names(streams)
    ),
    class = "stream_cluster_evaluation"
  )
}

#' Cluster stream functions
#'
#' Performs K-means clustering on mean-centered common-grid stream-function
#' values and returns arithmetic mean stream functions for each cluster.
#'
#' @inheritParams evaluate_stream_clusters
#' @param k Selected number of clusters.
#'
#' @return A `stream_clusters` object containing assignments, mean stream
#'   functions, and the underlying K-means result.
#' @export
cluster_streams <- function(streams,
                            k,
                            method = c("kmeans"),
                            nstart = 25L,
                            iter.max = 100L,
                            seed = NULL) {
  method <- match.arg(method)
  prepared <- prepare_stream_clustering_input(streams)
  streams <- prepared$streams
  limits <- feature_cluster_limits(prepared$features)
  if (length(k) == 1L && is.finite(k) && k > limits$max_k) {
    cli::cli_abort(c(
      "{.arg k} cannot exceed {limits$max_k} for these stream functions.",
      "i" = "The maximum is constrained by {limits$n_distinct} distinct stream function{?s} and {limits$n_objects} total stream function{?s}."
    ))
  }
  model <- with_clustering_seed(
    seed,
    fit_kmeans_features(prepared$features, k, nstart, iter.max)
  )

  center_values <- lapply(seq_len(nrow(model$centers)), function(i) {
    values <- model$centers[i, ]
    dim(values) <- dim(prepared$reference$values)
    values
  })
  centers <- lapply(seq_along(center_values), function(i) {
    build_stream_cluster_center(streams[[1L]], center_values[[i]], i)
  })
  object_names <- names(streams)
  if (is.null(object_names)) {
    object_names <- rep.int(NA_character_, length(streams))
  }
  assignments <- data.frame(
    stream_index = seq_along(streams),
    stream_name = object_names,
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
      center_values = center_values,
      kmeans = model,
      method = method,
      k = as.integer(k),
      streams = streams
    ),
    class = "stream_clusters"
  )
}

#' Autoplot stream-cluster evaluation
#'
#' @param object A `stream_cluster_evaluation` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return A ggplot elbow plot.
#' @export
autoplot.stream_cluster_evaluation <- function(object, ...) {
  ggplot2::ggplot(
    object$metrics,
    ggplot2::aes(x = .data$k, y = .data$within_variance)
  ) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = object$metrics$k) +
    ggplot2::labs(
      x = "Number of clusters (K)",
      y = "Within-cluster variance (mean squared distance)",
      title = "Stream-function clustering elbow plot"
    ) +
    ggplot2::theme_bw()
}

#' Autoplot stream-clustering centers
#'
#' @param object A `stream_clusters` object returned by [cluster_streams()].
#' @param contour Logical indicating whether to overlay contour lines.
#' @param ... Additional arguments, currently unused.
#'
#' @return A faceted ggplot of cluster-mean stream functions.
#' @export
autoplot.stream_clusters <- function(object, contour = TRUE, ...) {
  center_data <- do.call(rbind, lapply(
    seq_along(object$centers),
    function(i) {
      data <- object$centers[[i]]$grid
      data$cluster <- factor(i, levels = seq_along(object$centers))
      data
    }
  ))
  center_data <- add_stream_grid_cell_bounds(center_data)
  plot <- ggplot2::ggplot(center_data) +
    ggplot2::geom_rect(ggplot2::aes(
      fill = .data$A,
      xmin = .data$cell_xmin,
      xmax = .data$cell_xmax,
      ymin = .data$cell_ymin,
      ymax = .data$cell_ymax
    )) +
    ggplot2::facet_wrap(ggplot2::vars(cluster)) +
    ggplot2::scale_fill_viridis_c(name = "A") +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(x = "x", y = "y", title = "Cluster-mean stream functions") +
    ggplot2::theme_bw()
  if (isTRUE(contour)) {
    plot <- plot + ggplot2::geom_contour(
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        z = .data$A,
        group = .data$cluster
      ),
      color = "white",
      alpha = 0.6,
      show.legend = FALSE,
      inherit.aes = FALSE
    )
  }
  plot
}

resolve_joint_clustering_objects <- function(landscapes, streams) {
  if (inherits(landscapes, "group_dynamics")) {
    group <- landscapes
    landscapes <- group$landscapes
    if (is.null(streams)) {
      streams <- group$streams
    }
  }
  if (is.null(streams)) {
    cli::cli_abort("{.arg streams} must be supplied unless {.arg landscapes} is a {.cls group_dynamics} object.")
  }
  if (!is.list(landscapes) || !is.list(streams) || length(landscapes) != length(streams)) {
    cli::cli_abort("{.arg landscapes} and {.arg streams} must be lists of equal length.")
  }

  landscape_names <- names(landscapes)
  stream_names <- names(streams)
  if (!is.null(landscape_names) && !is.null(stream_names)) {
    if (anyDuplicated(landscape_names) || anyDuplicated(stream_names) ||
        !setequal(landscape_names, stream_names)) {
      cli::cli_abort("Named landscapes and streams must have the same unique names.")
    }
    streams <- streams[landscape_names]
  }

  list(landscapes = landscapes, streams = streams)
}

mean_pairwise_squared_distance <- function(features) {
  if (nrow(features) < 2L) {
    return(0)
  }
  centered <- sweep(features, 2L, colMeans(features), "-")
  2 * sum(centered^2) / (nrow(features) - 1L)
}

prepare_joint_clustering_input <- function(landscapes,
                                           streams,
                                           landscape_weight,
                                           stream_weight) {
  objects <- resolve_joint_clustering_objects(landscapes, streams)
  landscape_input <- prepare_landscape_clustering_input(objects$landscapes)
  stream_input <- prepare_stream_clustering_input(objects$streams)
  if (!identical(landscape_input$reference$dimension, "2d")) {
    cli::cli_abort("Joint landscape-stream clustering requires two-dimensional landscapes.")
  }
  if (!is.numeric(landscape_weight) || length(landscape_weight) != 1L ||
      !is.finite(landscape_weight) || landscape_weight < 0 ||
      !is.numeric(stream_weight) || length(stream_weight) != 1L ||
      !is.finite(stream_weight) || stream_weight < 0 ||
      landscape_weight + stream_weight <= 0) {
    cli::cli_abort("{.arg landscape_weight} and {.arg stream_weight} must be finite non-negative numbers with at least one positive value.")
  }

  landscape_mean_sq <- mean_pairwise_squared_distance(landscape_input$features)
  stream_mean_sq <- mean_pairwise_squared_distance(stream_input$features)
  tolerance <- sqrt(.Machine$double.eps)
  landscape_scale <- if (landscape_mean_sq > tolerance) sqrt(landscape_mean_sq) else 1
  stream_scale <- if (stream_mean_sq > tolerance) sqrt(stream_mean_sq) else 1
  landscape_scaled <- sqrt(landscape_weight) * landscape_input$features / landscape_scale
  stream_scaled <- sqrt(stream_weight) * stream_input$features / stream_scale
  features <- cbind(landscape_scaled, stream_scaled)

  list(
    features = features,
    landscape_input = landscape_input,
    stream_input = stream_input,
    landscapes = objects$landscapes,
    streams = objects$streams,
    normalization = list(
      landscape_mean_squared_distance = landscape_mean_sq,
      stream_mean_squared_distance = stream_mean_sq,
      landscape_scale = landscape_scale,
      stream_scale = stream_scale,
      landscape_weight = landscape_weight,
      stream_weight = stream_weight,
      landscape_active = landscape_mean_sq > tolerance && landscape_weight > 0,
      stream_active = stream_mean_sq > tolerance && stream_weight > 0,
      normalized_landscape_mean_squared_distance =
        mean_pairwise_squared_distance(landscape_scaled),
      normalized_stream_mean_squared_distance =
        mean_pairwise_squared_distance(stream_scaled)
    )
  )
}

#' Evaluate joint landscape-stream clusters
#'
#' Combines steady-state-distribution and stream-function representations for
#' joint K-means clustering. Each feature block is divided by the square root
#' of its mean pairwise squared distance before concatenation. Consequently,
#' with the default equal weights, landscape and stream variation make equal
#' average contributions to squared Euclidean distance. Stream functions are
#' mean-centered first to remove their arbitrary additive constants.
#'
#' @param landscapes A list of two-dimensional landscapes or a
#'   `group_dynamics` object.
#' @param streams Optional corresponding list of stream functions. It is taken
#'   from `landscapes$streams` when a `group_dynamics` object is supplied.
#' @param k_values Candidate cluster counts. `NULL` chooses the valid range up
#'   to 10 based on distinct joint representations and sample size.
#' @param landscape_weight,stream_weight Non-negative weights for the two
#'   normalized squared-distance contributions. Both default to 1.
#' @param method Clustering method. Currently only `"kmeans"` is available.
#' @param nstart Number of K-means random initializations.
#' @param iter.max Maximum number of K-means iterations.
#' @param seed Optional random seed. The caller's random-number state is
#'   restored afterward.
#'
#' @return A `landscape_stream_cluster_evaluation` object containing candidate
#'   metrics, fitted models, and modality-normalization diagnostics.
#' @export
evaluate_landscape_stream_clusters <- function(landscapes,
                                                streams = NULL,
                                                k_values = NULL,
                                                landscape_weight = 1,
                                                stream_weight = 1,
                                                method = c("kmeans"),
                                                nstart = 25L,
                                                iter.max = 100L,
                                                seed = NULL) {
  method <- match.arg(method)
  prepared <- prepare_joint_clustering_input(
    landscapes,
    streams,
    landscape_weight,
    stream_weight
  )
  limits <- feature_cluster_limits(prepared$features)
  if (is.null(k_values)) {
    k_values <- seq_len(min(10L, limits$max_k))
  }
  if (!length(k_values) || any(!is.finite(k_values)) ||
      any(k_values != as.integer(k_values)) || any(k_values < 1L) ||
      any(k_values > limits$max_k)) {
    cli::cli_abort(c(
      "{.arg k_values} must contain integers between 1 and {limits$max_k}.",
      "i" = "The maximum is constrained by {limits$n_distinct} distinct joint representation{?s} and {limits$n_objects} total observation{?s}."
    ))
  }
  k_values <- sort(unique(as.integer(k_values)))
  models <- with_clustering_seed(seed, lapply(k_values, function(k) {
    fit_kmeans_features(prepared$features, k, nstart, iter.max)
  }))
  names(models) <- as.character(k_values)
  metrics <- data.frame(
    k = k_values,
    within_variance = vapply(models, function(x) x$tot.withinss, numeric(1)) / limits$n_objects,
    between_variance = vapply(models, function(x) x$betweenss, numeric(1)) / limits$n_objects,
    total_variance = vapply(models, function(x) x$totss, numeric(1)) / limits$n_objects
  )
  metrics$explained_variance <- ifelse(
    metrics$total_variance > 0,
    metrics$between_variance / metrics$total_variance,
    0
  )
  object_names <- names(prepared$landscapes)
  if (is.null(object_names)) {
    object_names <- names(prepared$streams)
  }

  structure(
    list(
      metrics = metrics,
      models = models,
      normalization = prepared$normalization,
      method = method,
      n_objects = limits$n_objects,
      object_names = object_names
    ),
    class = "landscape_stream_cluster_evaluation"
  )
}

#' Cluster landscapes and streams jointly
#'
#' Fits a fixed-K solution using the normalized joint representation described
#' in [evaluate_landscape_stream_clusters()]. Cluster centers are returned
#' separately as mean steady-state landscapes and mean-centered stream
#' functions on their original scales.
#'
#' @inheritParams evaluate_landscape_stream_clusters
#' @param k Selected number of clusters.
#'
#' @return A `landscape_stream_clusters` object containing assignments,
#'   landscape and stream centers, normalization diagnostics, and K-means fit.
#' @export
cluster_landscape_streams <- function(landscapes,
                                      streams = NULL,
                                      k,
                                      landscape_weight = 1,
                                      stream_weight = 1,
                                      method = c("kmeans"),
                                      nstart = 25L,
                                      iter.max = 100L,
                                      seed = NULL) {
  method <- match.arg(method)
  prepared <- prepare_joint_clustering_input(
    landscapes,
    streams,
    landscape_weight,
    stream_weight
  )
  limits <- feature_cluster_limits(prepared$features)
  if (length(k) == 1L && is.finite(k) && k > limits$max_k) {
    cli::cli_abort(c(
      "{.arg k} cannot exceed {limits$max_k} for these joint representations.",
      "i" = "The maximum is constrained by {limits$n_distinct} distinct joint representation{?s} and {limits$n_objects} total observation{?s}."
    ))
  }
  model <- with_clustering_seed(
    seed,
    fit_kmeans_features(prepared$features, k, nstart, iter.max)
  )

  fitted_k <- nrow(model$centers)
  landscape_centers <- lapply(seq_len(fitted_k), function(i) {
    density <- colMeans(
      prepared$landscape_input$features[model$cluster == i, , drop = FALSE]
    )
    dim(density) <- dim(prepared$landscape_input$reference$values)
    build_landscape_cluster_center(prepared$landscapes[[1L]], density, i)
  })
  stream_centers <- lapply(seq_len(fitted_k), function(i) {
    values <- colMeans(
      prepared$stream_input$features[model$cluster == i, , drop = FALSE]
    )
    dim(values) <- dim(prepared$stream_input$reference$values)
    build_stream_cluster_center(prepared$streams[[1L]], values, i)
  })
  object_names <- names(prepared$landscapes)
  if (is.null(object_names)) {
    object_names <- names(prepared$streams)
  }
  if (is.null(object_names)) {
    object_names <- rep.int(NA_character_, length(prepared$landscapes))
  }
  assignments <- data.frame(
    object_index = seq_along(prepared$landscapes),
    object_name = object_names,
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
      centers = list(landscapes = landscape_centers, streams = stream_centers),
      landscape_centers = landscape_centers,
      stream_centers = stream_centers,
      normalization = prepared$normalization,
      kmeans = model,
      method = method,
      k = fitted_k,
      landscapes = prepared$landscapes,
      streams = prepared$streams
    ),
    class = "landscape_stream_clusters"
  )
}

#' Autoplot joint landscape-stream cluster evaluation
#'
#' @param object A `landscape_stream_cluster_evaluation` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return A ggplot elbow plot.
#' @export
autoplot.landscape_stream_cluster_evaluation <- function(object, ...) {
  ggplot2::ggplot(
    object$metrics,
    ggplot2::aes(x = .data$k, y = .data$within_variance)
  ) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = object$metrics$k) +
    ggplot2::labs(
      x = "Number of clusters (K)",
      y = "Within-cluster variance (mean squared distance)",
      title = "Joint landscape-stream clustering elbow plot"
    ) +
    ggplot2::theme_bw()
}

#' Autoplot joint landscape-stream cluster centers
#'
#' @param object A `landscape_stream_clusters` object.
#' @param type Which center representation to draw, `"landscape"` or
#'   `"stream"`.
#' @param ... Additional arguments passed to the selected center autoplot.
#'
#' @return A faceted ggplot of landscape or stream centers.
#' @export
autoplot.landscape_stream_clusters <- function(
    object,
    type = c("landscape", "stream"),
    ...) {
  type <- match.arg(type)
  if (type == "landscape") {
    centers <- structure(
      list(centers = object$landscape_centers),
      class = "landscape_clusters"
    )
  } else {
    centers <- structure(
      list(centers = object$stream_centers),
      class = "stream_clusters"
    )
  }
  autoplot(centers, ...)
}
