make_test_landscape <- function(density, x = 1:2, y = 1:2) {
  density <- matrix(density, nrow = length(x), ncol = length(y))
  potential <- -log(density)
  dist <- expand.grid(x = x, y = y)
  dist$d <- as.numeric(density)
  dist$U <- as.numeric(potential)
  dist$U_plot <- as.numeric(potential)
  structure(
    list(ss = density, dist = dist, plot = NULL, plot_2 = NULL, vf = NULL),
    class = c("2d_static_ld", "2d_ld", "landscape")
  )
}

test_that("landscape K-means uses density values and returns landscape centers", {
  landscapes <- list(
    a1 = make_test_landscape(c(0.70, 0.10, 0.10, 0.10)),
    a2 = make_test_landscape(c(0.65, 0.15, 0.10, 0.10)),
    b1 = make_test_landscape(c(0.10, 0.10, 0.10, 0.70)),
    b2 = make_test_landscape(c(0.10, 0.10, 0.15, 0.65))
  )

  result <- cluster_landscapes(landscapes, k = 2, seed = 1)

  expect_s3_class(result, "landscape_clusters")
  expect_equal(result$cluster[1], result$cluster[2])
  expect_equal(result$cluster[3], result$cluster[4])
  expect_false(result$cluster[1] == result$cluster[3])
  expect_length(result$centers, 2L)
  expect_true(all(vapply(result$centers, inherits, logical(1), "landscape")))
  expect_equal(
    result$centers[[1]]$dist$U,
    -log(result$centers[[1]]$dist$d)
  )
  evaluation <- evaluate_landscape_clusters(landscapes, k_values = 1:3, seed = 1)
  expect_equal(evaluation$metrics$k, 1:3)
  expect_true(all(diff(evaluation$metrics$within_variance) <= 0))
  expect_error(cluster_landscapes(rep(landscapes[1], 3), k = 2), "distinct landscape")
  expect_error(
    cluster_landscapes(list(landscapes[[1]], make_test_landscape(
      c(0.7, 0.1, 0.1, 0.1), x = c(1, 3)
    )), k = 2),
    "same dimension and grid"
  )
})

test_that("fit_group_dynamics reuses individual workflow with pooled limits", {
  calls <- list()

  local_mocked_bindings(
    fit_individual_dynamics = function(data, x, y, lims, ...) {
      calls[[length(calls) + 1L]] <<- list(data = data, lims = lims)
      index <- length(calls)
      landscape <- structure(
        list(ss = matrix(c(index, 1, 1, 1) / (index + 3), 2, 2)),
        class = c("2d_ld", "landscape")
      )
      structure(
        list(
          data = data,
          vectorfield = structure(list(), class = "vectorfield"),
          landscape = landscape,
          probability_flow = structure(list(), class = "2d_pf"),
          stream = structure(list(), class = "2d_stream")
        ),
        class = "individual_dynamics"
      )
    },
    .package = "fitlandr"
  )

  data <- list(
    first = data.frame(x = 0:1, y = 0:1),
    second = data.frame(x = 9:10, y = 19:20)
  )
  result <- fit_group_dynamics(
    data,
    x = "x",
    y = "y",
    h_values = 0.2,
    cv_folds = 2,
    n_grid = 5,
    flow_n = 3
  )

  expect_s3_class(result, "group_dynamics")
  expect_equal(names(result$vectorfields), c("first", "second"))
  expect_equal(names(result$landscapes), c("first", "second"))
  expect_equal(names(result$streams), c("first", "second"))
  expect_equal(result$lims, c(-1, 11, -2, 22))
  expect_equal(calls[[1]]$lims, result$lims)
  expect_equal(calls[[2]]$lims, result$lims)
})

test_that("stream clustering centers A and accepts group dynamics", {
  make_stream <- function(A, offset = 0) {
    grid <- expand.grid(x = c(0, 1, 3), y = c(0, 2, 5))
    grid$A <- A + offset
    structure(
      list(grid = grid, pf = list(x = "x", y = "y")),
      class = "2d_stream"
    )
  }
  streams <- list(
    clockwise_1 = make_stream(seq(-1, 1, length.out = 9), offset = 20),
    clockwise_2 = make_stream(seq(-0.9, 0.9, length.out = 9), offset = -7),
    counter_1 = make_stream(seq(1, -1, length.out = 9), offset = 12),
    counter_2 = make_stream(seq(0.9, -0.9, length.out = 9), offset = -3)
  )
  group <- structure(list(streams = streams), class = "group_dynamics")

  evaluation <- evaluate_stream_clusters(group, k_values = 1:2, seed = 1)
  clusters <- cluster_streams(group, k = 2, seed = 1)

  expect_s3_class(evaluation, "stream_cluster_evaluation")
  expect_s3_class(clusters, "stream_clusters")
  expect_equal(clusters$assignments$stream_name, names(streams))
  expect_equal(clusters$cluster[1], clusters$cluster[2])
  expect_equal(clusters$cluster[3], clusters$cluster[4])
  expect_false(clusters$cluster[1] == clusters$cluster[3])
  expect_true(all(vapply(
    clusters$centers,
    function(center) abs(mean(center$grid$A)) < 1e-12,
    logical(1)
  )))
  expect_s3_class(autoplot(evaluation), "ggplot")
  center_plot <- autoplot(clusters)
  expect_s3_class(center_plot, "ggplot")
})

test_that("joint landscape-stream clustering balances both modalities", {
  make_stream <- function(A, multiplier = 1, offset = 0) {
    grid <- expand.grid(x = 1:2, y = 1:2)
    grid$A <- multiplier * A + offset
    structure(
      list(grid = grid, pf = list(x = "x", y = "y")),
      class = "2d_stream"
    )
  }
  landscapes <- list(
    a = make_test_landscape(c(0.45, 0.25, 0.2, 0.1)),
    b = make_test_landscape(c(0.4, 0.3, 0.2, 0.1)),
    c = make_test_landscape(c(0.1, 0.2, 0.25, 0.45)),
    d = make_test_landscape(c(0.1, 0.2, 0.3, 0.4))
  )
  stream_patterns <- list(
    a = c(-1, 0, 0, 1),
    b = c(-0.8, 0, 0, 0.8),
    c = c(1, 0, 0, -1),
    d = c(0.8, 0, 0, -0.8)
  )
  streams <- Map(make_stream, stream_patterns, offset = c(10, -5, 20, -8))
  group <- structure(
    list(landscapes = landscapes, streams = streams),
    class = "group_dynamics"
  )

  evaluation <- evaluate_landscape_stream_clusters(
    group,
    k_values = 1:2,
    seed = 1
  )
  clusters <- cluster_landscape_streams(group, k = 2, seed = 1)

  expect_s3_class(evaluation, "landscape_stream_cluster_evaluation")
  expect_s3_class(clusters, "landscape_stream_clusters")
  expect_equal(
    evaluation$normalization$normalized_landscape_mean_squared_distance,
    1
  )
  expect_equal(
    evaluation$normalization$normalized_stream_mean_squared_distance,
    1
  )
  expect_length(clusters$landscape_centers, 2)
  expect_length(clusters$stream_centers, 2)
  scaled_streams <- Map(
    make_stream,
    stream_patterns,
    multiplier = rep(1000, 4),
    offset = c(100, -50, 200, -80)
  )
  scaled_evaluation <- evaluate_landscape_stream_clusters(
    landscapes,
    scaled_streams,
    k_values = 1:2,
    seed = 1
  )
  expect_equal(scaled_evaluation$metrics, evaluation$metrics, tolerance = 1e-10)
})
