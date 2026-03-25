test_that("plot.cv_vectorfield uses x object and returns ggplot", {
  cv_obj <- structure(
    list(
      cv_results = data.frame(h = c(0.1, 0.2), cv_mse = c(1, 0.5)),
      h_optimal = 0.2
    ),
    class = "cv_vectorfield"
  )

  p <- plot(cv_obj)
  expect_s3_class(p, "ggplot")
})

test_that("fit_2d_ld supports vector_position = middle", {
  set.seed(1)
  d <- data.frame(x = cumsum(stats::rnorm(50)))

  ld <- fit_2d_ld(d, x = "x", vector_position = "middle", n = 20)
  expect_s3_class(ld, "2d_MVKE_landscape")
  expect_true(all(c("x", "U") %in% names(ld$dist)))
})

test_that("plot.summary_bootstrap_2d_ld consumes per_point field", {
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

  p <- plot(summary_obj)
  expect_s3_class(p, "ggplot")
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
