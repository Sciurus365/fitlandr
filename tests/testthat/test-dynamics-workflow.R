test_that("fit_individual_dynamics composes the complete workflow", {
  calls <- new.env(parent = emptyenv())

  local_mocked_bindings(
    cv_fit_2d_vf = function(data, x, y, h_values, k, lims, n, ...) {
      calls$cv <- list(lims = lims, n = n, k = k)
      list(final_model = structure(list(x = x, y = y, lims = lims), class = "vectorfield"))
    },
    add_interp_grid = function(vf, lims, n) {
      calls$interp <- list(lims = lims, n = n)
      vf$interp_grid <- TRUE
      vf
    },
    make_2d_ld = function(vf, n_grid, ...) {
      calls$landscape <- n_grid
      structure(list(vf = vf, ss = matrix(1 / 4, 2, 2)), class = c("2d_ld", "landscape"))
    },
    make_2d_pf = function(vf, ld, n, divided_by_rho, ...) {
      calls$flow <- list(n = n, divided_by_rho = divided_by_rho)
      structure(list(vf = vf, ld = ld), class = "2d_pf")
    },
    make_2d_stream = function(pf) {
      calls$stream <- TRUE
      structure(list(pf = pf), class = "2d_stream")
    },
    .package = "fitlandr"
  )

  result <- fit_individual_dynamics(
    data.frame(x = 0:5, y = 5:0),
    x = "x",
    y = "y",
    lims = c(-1, 6, -1, 6),
    h_values = c(0.1, 0.2),
    cv_folds = 3,
    n_grid = 12,
    flow_n = 7
  )

  expect_s3_class(result, "individual_dynamics")
  expect_equal(result$lims, c(-1, 6, -1, 6))
  expect_equal(calls$cv, list(lims = c(-1, 6, -1, 6), n = 12L, k = 3L))
  expect_equal(calls$interp$n, 12L)
  expect_equal(calls$landscape, 12L)
  expect_equal(calls$flow, list(n = 7L, divided_by_rho = FALSE))
  expect_true(calls$stream)
  expect_s3_class(result$stream, "2d_stream")
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
  expect_equal(
    names(result),
    c("vectorfields", "landscapes", "streams", "x", "y", "lims", "settings")
  )
  expect_equal(result$lims, c(-1, 11, -2, 22))
  expect_equal(calls[[1]]$lims, result$lims)
  expect_equal(calls[[2]]$lims, result$lims)
})

test_that("group fitting suppresses routine messages but preserves warnings", {
  expect_warning(
    expect_message(
      value <- fitlandr:::fit_group_member_quietly({
        cli::cli_inform("routine internal output")
        warning("important fitting warning")
        42
      }),
      NA
    ),
    "important fitting warning"
  )
  expect_equal(value, 42)
})

test_that("group dynamics autoplot facets landscapes and streams", {
  make_landscape <- function(offset) {
    dist <- expand.grid(x = 1:3, y = 1:3)
    dist$U <- seq(0, 1, length.out = 9) + offset
    dist$U_plot <- dist$U
    structure(
      list(dist = dist, ss = matrix(exp(-dist$U), 3, 3)),
      class = c("2d_static_ld", "2d_ld", "landscape")
    )
  }
  make_stream <- function(offset) {
    grid <- expand.grid(x = c(0, 1, 3), y = c(0, 2, 5))
    grid$A <- seq(-1, 1, length.out = 9) + offset
    structure(
      list(grid = grid, pf = list(x = "x", y = "y")),
      class = "2d_stream"
    )
  }
  group <- structure(
    list(
      landscapes = list(first = make_landscape(10), second = make_landscape(-4)),
      streams = list(first = make_stream(20), second = make_stream(-7)),
      x = "state_x",
      y = "state_y"
    ),
    class = "group_dynamics"
  )

  landscape_plot <- autoplot(group, type = "landscape", ncol = 2)
  stream_plot <- autoplot(group, type = "stream", ncol = 2)

  expect_s3_class(landscape_plot, "ggplot")
  expect_s3_class(stream_plot, "ggplot")
  expect_equal(levels(landscape_plot$data$individual), c("first", "second"))
  expect_equal(levels(stream_plot$data$individual), c("first", "second"))
  expect_s3_class(stream_plot$layers[[1L]]$geom, "GeomRect")
  expect_warning(ggplot2::ggplot_build(landscape_plot), NA)
  expect_warning(ggplot2::ggplot_build(stream_plot), NA)
})

test_that("landscape clustering defaults respect the number of distinct inputs", {
  make_landscape <- function(density) {
    structure(
      list(
        ss = matrix(density, 2, 2),
        dist = expand.grid(x = 1:2, y = 1:2)
      ),
      class = c("2d_ld", "landscape")
    )
  }
  identical_landscapes <- list(
    make_landscape(rep(0.25, 4)),
    make_landscape(rep(0.25, 4))
  )

  evaluation <- evaluate_landscape_clusters(identical_landscapes)
  expect_equal(evaluation$metrics$k, 1L)
  expect_error(
    evaluate_landscape_clusters(identical_landscapes, k_values = 2),
    "maximum is constrained by 1 distinct landscape"
  )
  expect_error(
    cluster_landscapes(identical_landscapes, k = 2),
    "maximum is constrained by 1 distinct landscape"
  )
})

test_that("K-means cluster count must be smaller than the sample size", {
  make_landscape <- function(density) {
    structure(
      list(
        ss = matrix(density, 2, 2),
        dist = expand.grid(x = 1:2, y = 1:2)
      ),
      class = c("2d_ld", "landscape")
    )
  }
  two_landscapes <- list(
    make_landscape(c(0.4, 0.3, 0.2, 0.1)),
    make_landscape(c(0.1, 0.2, 0.3, 0.4))
  )

  expect_equal(evaluate_landscape_clusters(two_landscapes)$metrics$k, 1L)
  expect_error(
    cluster_landscapes(two_landscapes, k = 2),
    "2 total landscapes"
  )
})

test_that("landscape clustering accepts a group_dynamics object", {
  make_landscape <- function(density) {
    structure(
      list(
        ss = matrix(density, 2, 2),
        dist = expand.grid(x = 1:2, y = 1:2)
      ),
      class = c("2d_ld", "landscape")
    )
  }
  landscapes <- list(
    a = make_landscape(c(0.4, 0.3, 0.2, 0.1)),
    b = make_landscape(c(0.35, 0.35, 0.2, 0.1)),
    c = make_landscape(c(0.1, 0.2, 0.3, 0.4)),
    d = make_landscape(c(0.1, 0.2, 0.35, 0.35))
  )
  group <- structure(list(landscapes = landscapes), class = "group_dynamics")

  evaluation <- evaluate_landscape_clusters(group, k_values = 1:2, seed = 1)
  clusters <- cluster_landscapes(group, k = 2, seed = 1)

  expect_s3_class(evaluation, "landscape_cluster_evaluation")
  expect_s3_class(clusters, "landscape_clusters")
  expect_equal(clusters$assignments$landscape_name, names(landscapes))
  expect_s3_class(autoplot(clusters), "ggplot")
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
  expect_s3_class(center_plot$layers[[1L]]$geom, "GeomRect")
  expect_warning(ggplot2::ggplot_build(center_plot), NA)
})

test_that("joint landscape-stream clustering balances both modalities", {
  make_landscape <- function(density) {
    structure(
      list(
        ss = matrix(density, 2, 2),
        dist = expand.grid(x = 1:2, y = 1:2)
      ),
      class = c("2d_ld", "landscape")
    )
  }
  make_stream <- function(A, multiplier = 1, offset = 0) {
    grid <- expand.grid(x = 1:2, y = 1:2)
    grid$A <- multiplier * A + offset
    structure(
      list(grid = grid, pf = list(x = "x", y = "y")),
      class = "2d_stream"
    )
  }
  landscapes <- list(
    a = make_landscape(c(0.45, 0.25, 0.2, 0.1)),
    b = make_landscape(c(0.4, 0.3, 0.2, 0.1)),
    c = make_landscape(c(0.1, 0.2, 0.25, 0.45)),
    d = make_landscape(c(0.1, 0.2, 0.3, 0.4))
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
  expect_equal(evaluation$object_names, names(landscapes))
  expect_equal(clusters$assignments$object_name, names(landscapes))
  expect_length(clusters$landscape_centers, 2)
  expect_length(clusters$stream_centers, 2)
  expect_s3_class(autoplot(evaluation), "ggplot")
  expect_s3_class(autoplot(clusters, type = "landscape"), "ggplot")
  stream_plot <- autoplot(clusters, type = "stream")
  expect_s3_class(stream_plot, "ggplot")
  expect_s3_class(stream_plot$layers[[1L]]$geom, "GeomRect")
  expect_warning(ggplot2::ggplot_build(stream_plot), NA)

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
