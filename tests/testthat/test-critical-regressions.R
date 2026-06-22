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
  x_coords <- seq(0, 1, length.out = 5)
  y_coords <- seq(0, 1, length.out = 5)
  original_ld <- structure(
    list(
      ss = structure(
        matrix(0, nrow = 5, ncol = 5),
        x_coords = x_coords,
        y_coords = y_coords
      )
    ),
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

test_that("fit_2d_vf excludes cross-day and non-consecutive transitions", {
  d <- data.frame(
    x = c(0, 1, 10, 11, 20),
    y = c(0, 1, 10, 11, 20),
    day = c(1, 1, 2, 2, 2),
    beep = c(1, 2, 1, 3, 4)
  )

  local_mocked_bindings(
    MVKE = function(d, v, ...) {
      force(d)
      force(v)
      function(pos) list(mu = c(0, 0), a = diag(2))
    },
    .package = "fitlandr"
  )

  vf <- fit_2d_vf(
    d,
    x = "x",
    y = "y",
    dayvar = "day",
    beepvar = "beep",
    na_action = "omit_vectors",
    method = "MVKE",
    n = 3
  )

  expect_equal(nrow(vf$original_vectors), 2L)
  expect_equal(unname(vf$original_vectors[, c("x", "y")]), matrix(c(0, 0, 11, 11), ncol = 2, byrow = TRUE))
})

test_that("fit_2d_vf warns when time-boundary arguments are ineffective", {
  d <- data.frame(
    x = c(0, 1, 10),
    y = c(0, 1, 10),
    day = c(1, 1, 2)
  )

  local_mocked_bindings(
    MVKE = function(d, v, ...) {
      force(d)
      force(v)
      function(pos) list(mu = c(0, 0), a = diag(2))
    },
    .package = "fitlandr"
  )

  expect_warning(
    fit_2d_vf(
      d,
      x = "x",
      y = "y",
      dayvar = "day",
      na_action = "omit_data_points",
      method = "MVKE",
      n = 3
    ),
    "vectors can still connect observations across day boundaries or non-consecutive beeps"
  )
})

test_that("fit_1d_vf and make_1d_ld return the new 1D classes", {
  set.seed(1)
  d <- data.frame(x = cumsum(stats::rnorm(50)))

  vf <- fit_1d_vf(d, x = "x", method = "MVKE", n = 25)
  vf <- add_interp_grid(vf)
  ld <- make_1d_ld(vf, n_grid = 50)

  expect_s3_class(vf, "1d_vectorfield")
  expect_s3_class(ld, "1d_static_ld")
  expect_true(all(c("x", "vx", "v_norm") %in% names(vf$vec_grid)))
  expect_true(all(c("x", "d", "U") %in% names(ld$dist)))
})

test_that("summary.bootstrap_1d_ld returns a structured summary", {
  set.seed(1)
  d <- data.frame(x = cumsum(stats::rnorm(60)))
  vf <- fit_1d_vf(d, x = "x", method = "MVKE", n = 25)
  vf <- add_interp_grid(vf)
  boot_vf <- bootstrap_1d_vf(vf, n_boot = 5, block_length = 4)
  boot_vf$bootstrap_models <- lapply(boot_vf$bootstrap_models, add_interp_grid)
  boot_ld <- bootstrap_1d_ld(boot_vf, n_grid = 40)
  smry <- summary(boot_ld, clustering_method = "mean_potential")

  expect_s3_class(smry, "summary_bootstrap_1d_ld")
  expect_true(all(c("boot_index", "n_mins") %in% names(smry$per_boot)))
  expect_true(!is.null(smry$original_ld))
})
