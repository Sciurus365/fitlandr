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

test_that("fit_2d_vf excludes cross-day and non-consecutive transitions", {
  d <- data.frame(
    x = c(0, 1, 10, 11, 20),
    y = c(0, 1, 10, 11, 20),
    day = c(1, 1, 2, 2, 2),
    beep = c(1, 2, 1, 3, 4)
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

test_that("bandwidth units agree and every CV fold uses the same actual width", {
  d <- data.frame(
    x = c(seq_len(10), 30 + seq_len(10), 100 + seq_len(10)),
    y = sin(seq_len(30)) * 4
  )
  original_fit_1d <- fit_1d_vf
  original_fit_2d <- fit_2d_vf
  widths <- numeric()
  record_width <- function(model) {
    widths <<- c(widths, environment(model$MVKEresult)$h *
                   attr(model$data_normalized, "scale"))
    model
  }
  local_mocked_bindings(
    fit_1d_vf = function(...) record_width(original_fit_1d(...)),
    fit_2d_vf = function(...) record_width(original_fit_2d(...)),
    .package = "fitlandr"
  )

  for (dimension in 1:2) {
    fit <- if (dimension == 1) fit_1d_vf else fit_2d_vf
    cv_fit <- if (dimension == 1) cv_fit_1d_vf else cv_fit_2d_vf
    args <- list(data = d, x = "x", n = 3)
    if (dimension == 2) args$y <- "y"

    normalized_fit <- do.call(fit, c(args, list(h = 0.2)))
    original_fit <- do.call(fit, c(args, list(h = normalized_fit$h, h_scale = "original")))
    expect_equal(normalized_fit$vec_grid, original_fit$vec_grid)
    expect_equal(original_fit$h_normalized, 0.2)

    candidates <- c(0.2, 0.6)
    widths <- numeric()
    normalized_cv <- do.call(cv_fit, c(args, list(h_values = candidates, k = 3)))
    actual_candidates <- candidates * attr(normalized_fit$data_normalized, "scale")
    expect_equal(widths, c(rep(actual_candidates, each = 3), normalized_cv$h_optimal))
    expect_equal(normalized_cv$cv_results$h, actual_candidates)
    expect_equal(normalized_cv$cv_results$h_normalized, candidates)
    expect_equal(normalized_cv$final_model$h, normalized_cv$h_optimal)
    expect_equal(autoplot(normalized_cv)$data$h, actual_candidates)

    original_cv <- do.call(cv_fit, c(args, list(h_values = actual_candidates,
                                               h_scale = "original", k = 3)))
    expect_equal(normalized_cv$cv_results, original_cv$cv_results)
    expect_equal(normalized_cv$final_model$vec_grid, original_cv$final_model$vec_grid)
  }
})

test_that("summary.bootstrap_1d_ld returns a structured summary", {
  set.seed(1)
  d <- data.frame(x = cumsum(stats::rnorm(60)))
  vf <- fit_1d_vf(d, x = "x", method = "MVKE", n = 25)
  vf <- add_interp_grid(vf)
  boot_vf <- bootstrap_1d_vf(vf, n_boot = 5, block_length = 4)
  boot_vf$bootstrap_models <- lapply(boot_vf$bootstrap_models, add_interp_grid)
  boot_ld <- bootstrap_1d_ld(boot_vf, n_grid = 40)
  smry <- summary(boot_ld, minima_method = "mean_potential")

  expect_s3_class(smry, "summary_bootstrap_1d_ld")
  expect_true(all(c("bootstrap_minima", "per_minimum") %in% names(smry)))
  expect_identical(smry$params$minima_method, "mean_potential")
  expect_true(all(c("boot_index", "n_mins") %in% names(smry$per_boot)))
  expect_s3_class(smry$original_ld, "1d_static_ld")
})
