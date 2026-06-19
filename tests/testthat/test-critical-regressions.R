test_that("autoplot.cv_vectorfield returns ggplot and plot is deprecated", {
  cv_obj <- structure(
    list(
      cv_results = data.frame(h = c(0.1, 0.2), cv_mse = c(1, 0.5)),
      h_optimal = 0.2
    ),
    class = "cv_vectorfield"
  )

  p <- autoplot(cv_obj)
  expect_s3_class(p, "ggplot")
  expect_warning(plot(cv_obj), "deprecated")
})

test_that("autoplot and plotly_ld expose landscape plots", {
  gg <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  plotly <- plotly::plot_ly(x = 1, y = 1, z = 1)
  ld <- structure(
    list(plot = plotly, plot_2 = gg),
    class = c("2d_static_ld", "2d_ld", "landscape")
  )

  expect_s3_class(autoplot(ld), "ggplot")
  expect_s3_class(plotly_ld(ld), "plotly")
  expect_warning(expect_s3_class(plot(ld, 2), "ggplot"), "deprecated")
  expect_warning(expect_s3_class(plot(ld), "plotly"), "deprecated")
})

test_that("autoplot.vectorfield returns ggplot and plot is deprecated", {
  vf <- structure(
    list(
      vec_grid = data.frame(x = 0, y = 0, vx = 1, vy = 1, v_norm = sqrt(2)),
      data = data.frame(state_x = 0, state_y = 0),
      x = "state_x",
      y = "state_y"
    ),
    class = "vectorfield"
  )

  expect_s3_class(autoplot(vf), "ggplot")
  expect_warning(expect_s3_class(plot(vf), "ggplot"), "deprecated")
})

test_that("autoplot.2d_pf plots probability-flow vectors", {
  pf <- structure(
    list(
      vec_grid = data.frame(x = 0, y = 0, vx = 1, vy = 1),
      x = "state_x",
      y = "state_y"
    ),
    class = c("2d_pf", "probabilityflow", "vectorfield")
  )

  expect_s3_class(autoplot(pf), "ggplot")
  expect_warning(expect_s3_class(plot(pf), "ggplot"), "deprecated")
})

test_that("fit_2d_ld supports vector_position = middle", {
  set.seed(1)
  d <- data.frame(x = cumsum(stats::rnorm(50)))

  ld <- fit_2d_ld(d, x = "x", vector_position = "middle", n = 20)
  expect_s3_class(ld, "2d_MVKE_landscape")
  expect_true(all(c("x", "U") %in% names(ld$dist)))
})

test_that("autoplot.summary_bootstrap_2d_ld minima mode consumes per_point field", {
  x_coords <- seq(0, 1, length.out = 5)
  y_coords <- seq(0, 1, length.out = 5)

  summary_obj <- structure(
    list(
      per_point = data.frame(
        boot_index = c(1L, 2L),
        x = c(0.2, 0.8),
        y = c(0.3, 0.7),
        U = c(1.1, 0.9)
      ),
      original_ld = list(
        ss = structure(matrix(0, nrow = 2, ncol = 2),
          x_coords = x_coords,
          y_coords = y_coords
        ),
        vf = list(x = "x", y = "y")
      )
    ),
    class = "summary_bootstrap_2d_ld"
  )

  p <- autoplot(summary_obj, mode = "minima")
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.summary_bootstrap_2d_ld overlays only original major minima", {
  original_ld <- structure(
    list(),
    class = c("2d_static_ld", "2d_ld", "landscape")
  )
  summary_obj <- structure(
    list(
      per_point = data.frame(
        boot_index = 1L,
        x = 0,
        y = 0,
        U = 0,
        cluster = 1L
      ),
      per_cluster = NULL,
      original_ld = original_ld,
      params = list(min_barrier = 0.1)
    ),
    class = "summary_bootstrap_2d_ld"
  )

  local_mocked_bindings(
    find_loc_min = function(...) {
      list(mins = data.frame(
        x = c(1, 2),
        y = c(3, 4),
        U = c(0.1, 0.2),
        is_minor = c(FALSE, TRUE)
      ))
    },
    .package = "fitlandr"
  )

  p <- autoplot(summary_obj, mode = "clusters")

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$layers[[2]]$data), 1L)
  expect_false(p$layers[[2]]$data$is_minor)
})

test_that("summary/bootstrap validation errors are explicit", {
  expect_error(
    summary.bootstrap_2d_ld(list()),
    "must be a .*bootstrap_2d_ld",
    perl = TRUE
  )

  expect_error(
    autoplot.summary_bootstrap_2d_ld(list()),
    "must inherit from .*summary_bootstrap_2d_ld",
    perl = TRUE
  )
})
