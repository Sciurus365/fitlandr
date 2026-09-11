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

test_that("2D FVM uses half the infinitesimal covariance as diffusivity", {
  vf <- structure(list(lims = c(0, 1, 0, 1)), class = "vectorfield")
  local_mocked_bindings(
    predict.vectorfield = function(object, pos, ...) {
      list(v = c(0, 0), a = 2 * diag(2))
    },
    .package = "fitlandr"
  )

  rho <- ss_fp_2d(vf, n_grid = 2)
  generator <- as.matrix(attr(rho, "M"))
  faces <- attr(rho, "fvm_faces")

  expect_equal(generator[2, 1], 4, tolerance = 1e-12)
  expect_equal(generator[1, 1], -8, tolerance = 1e-12)
  expect_equal(max(abs(faces$divergence)), 0, tolerance = 1e-12)
})

test_that("make_2d_stream recovers a compatible staggered-grid stream function", {
  x_faces <- seq(0, 1, length.out = 4)
  y_faces <- seq(-1, 1, length.out = 4)
  x_centers <- (x_faces[-1L] + x_faces[-length(x_faces)]) / 2
  y_centers <- (y_faces[-1L] + y_faces[-length(y_faces)]) / 2
  hx <- diff(x_faces)[1L]
  hy <- diff(y_faces)[1L]
  A <- outer(
    x_faces,
    y_faces,
    function(x, y) sin(pi * x) * cos(pi * y / 2)
  )
  Jx <- (A[, -ncol(A), drop = FALSE] - A[, -1L, drop = FALSE]) / hy
  Jy <- (A[-1L, , drop = FALSE] - A[-nrow(A), , drop = FALSE]) / hx
  ld <- structure(
    list(
      ss = matrix(1, nrow = 3, ncol = 3),
      fvm_faces = list(
        Jx = Jx,
        Jy = Jy,
        x_faces = x_faces,
        y_faces = y_faces,
        x_centers = x_centers,
        y_centers = y_centers,
        hx = hx,
        hy = hy
      ),
      fvm_compatible = TRUE,
      cross_diffusion_mode = "drop"
    ),
    class = c("2d_static_ld", "2d_ld", "landscape")
  )
  vf <- structure(
    list(lims = c(0, 1, -1, 1), x = "state_x", y = "state_y"),
    class = "vectorfield"
  )
  pf <- make_2d_pf(vf, ld, n = 3)

  stream <- make_2d_stream(pf)

  expect_identical(pf$method, "fvm")
  expect_equal(nrow(pf$vec_grid), 9L)
  expect_equal(nrow(pf$face_grid), 24L)
  expect_s3_class(stream, "2d_stream")
  expect_identical(stream$method, "fvm")
  expect_equal(stream$corner_grid$A, as.numeric(A), tolerance = 1e-10)
  expect_equal(nrow(stream$grid), 9L)
  expect_equal(nrow(stream$corner_grid), 16L)
  expect_equal(nrow(stream$face_grid), 24L)
  expect_lt(stream$rmse, 1e-10)
  expect_lt(stream$relative_error, 1e-10)
})

test_that("make_2d_stream validates its probability-flow input", {
  divided_pf <- structure(
    list(
      vec_grid = expand.grid(x = 1:2, y = 1:2) |>
        transform(vx = 0, vy = 0),
      divided_by_rho = TRUE
    ),
    class = "2d_pf"
  )

  expect_error(make_2d_stream(divided_pf), "divided by density")
  divided_pf$divided_by_rho <- FALSE
  expect_error(make_2d_stream(divided_pf), "requires an FVM probability-flow object")
})

test_that("compatible stream method rejects unsupported landscape discretizations", {
  pf <- structure(
    list(
      method = "fvm",
      fvm_faces = list(invalid = TRUE),
      cross_diffusion_mode = "full",
      divided_by_rho = FALSE
    ),
    class = "2d_pf"
  )

  expect_error(make_2d_stream(pf), "cross_diffusion_mode")
})

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
  expect_true(all(result$assignments$distance >= 0))
  expect_length(result$centers, 2L)
  expect_true(all(vapply(result$centers, inherits, logical(1), "landscape")))
  expect_equal(
    result$centers[[1]]$dist$U,
    -log(result$centers[[1]]$dist$d)
  )
})

test_that("landscape cluster evaluation returns elbow metrics and plot", {
  landscapes <- list(
    make_test_landscape(c(0.70, 0.10, 0.10, 0.10)),
    make_test_landscape(c(0.65, 0.15, 0.10, 0.10)),
    make_test_landscape(c(0.10, 0.10, 0.10, 0.70)),
    make_test_landscape(c(0.10, 0.10, 0.15, 0.65))
  )

  evaluation <- evaluate_landscape_clusters(
    landscapes,
    k_values = 1:3,
    seed = 1
  )

  expect_s3_class(evaluation, "landscape_cluster_evaluation")
  expect_equal(evaluation$metrics$k, 1:3)
  expect_true(all(diff(evaluation$metrics$within_variance) <= 0))
  expect_s3_class(autoplot(evaluation), "ggplot")
})

test_that("landscape clustering requires a common grid", {
  landscapes <- list(
    make_test_landscape(c(0.70, 0.10, 0.10, 0.10)),
    make_test_landscape(c(0.70, 0.10, 0.10, 0.10), x = c(1, 3))
  )

  expect_error(cluster_landscapes(landscapes, k = 2), "same dimension and grid")
})

test_that("landscape clustering supports 1D landscapes", {
  make_1d_test_landscape <- function(density) {
    potential <- -log(density)
    structure(
      list(
        ss = density,
        dist = data.frame(
          x = seq_along(density),
          d = density,
          U = potential,
          U_plot = potential
        ),
        plot = NULL,
        plot_2 = NULL,
        vf = NULL
      ),
      class = c("1d_static_ld", "1d_ld", "landscape")
    )
  }
  landscapes <- list(
    make_1d_test_landscape(c(0.7, 0.2, 0.1)),
    make_1d_test_landscape(c(0.65, 0.25, 0.1)),
    make_1d_test_landscape(c(0.1, 0.2, 0.7))
  )

  result <- cluster_landscapes(landscapes, k = 2, seed = 1)

  expect_s3_class(result$centers[[1]], "1d_ld")
  expect_equal(length(result$centers[[1]]$ss), 3L)
  expect_s3_class(autoplot(result$centers[[1]]), "ggplot")
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

test_that("time separators tolerate separator rows already present", {
  d <- data.frame(
    x = c(0, 1, NA, 10, 11),
    y = c(0, 1, NA, 10, 11),
    day = c(1, 1, NA, 2, 2),
    beep = c(1, 2, NA, 1, 2)
  )

  separated <- insert_time_separators(
    d,
    columns = c("x", "y"),
    dayvar = "day",
    beepvar = "beep"
  )

  expect_equal(nrow(separated), nrow(d))
  expect_equal(which(is.na(separated$x)), 3L)
})

test_that("blocked CV works with day and beep separators in interior folds", {
  d <- data.frame(
    x = seq_len(30),
    y = seq_len(30) / 2,
    day = rep(seq_len(6), each = 5),
    beep = rep(seq_len(5), 6)
  )

  local_mocked_bindings(
    MVKE = function(d, v, ...) {
      force(d)
      force(v)
      function(pos) list(mu = c(0, 0), a = diag(2))
    },
    .package = "fitlandr"
  )

  cv <- cv_fit_2d_vf(
    d,
    x = "x",
    y = "y",
    dayvar = "day",
    beepvar = "beep",
    h_values = 0.2,
    k = 3,
    n = 3,
    lims = c(0, 31, 0, 16)
  )

  expect_s3_class(cv, "cv_vectorfield")
  expect_true(is.finite(cv$cv_results$cv_mse))
  expect_equal(cv$final_model$n, 3)
  expect_equal(nrow(cv$final_model$vec_grid), 9L)
})

test_that("blocked CV reports complete candidate failure explicitly", {
  d <- data.frame(x = seq_len(10), y = seq_len(10))

  local_mocked_bindings(
    fit_2d_vf = function(...) stop("synthetic fit failure"),
    .package = "fitlandr"
  )

  expect_error(
    cv_fit_2d_vf(d, x = "x", y = "y", h_values = 0.2, k = 2),
    "failed for every candidate bandwidth"
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
