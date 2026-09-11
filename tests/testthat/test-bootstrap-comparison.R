test_that("compare_minima_depths uses paired bootstrap contrasts", {
  object <- structure(
    list(
      params = list(level = 0.95),
      per_cluster = data.frame(
        cluster = c(1L, 2L),
        mean_U = c(2, 0.5)
      ),
      per_point = data.frame(
        boot_index = rep(1:3, each = 2),
        cluster = rep(1:2, 3),
        U = c(2, 1, 4, 2, 5, 2),
        is_noise = FALSE
      )
    ),
    class = "summary_bootstrap_2d_ld"
  )

  out <- compare_minima_depths(object, c(1, 2), level = 0.8)

  expect_s3_class(out, "minima_depth_comparison")
  expect_equal(out$result$delta_U, 1.5)
  expect_equal(out$result$bootstrap_mean, 2)
  expect_equal(out$per_boot$delta_U, c(1, 2, 3))
  expect_equal(out$result$n_paired, 3)
})

test_that("compare_minima_depths validates selected minima", {
  object <- structure(
    list(
      params = list(level = 0.95),
      per_cluster = data.frame(cluster = 1L, mean_U = 0),
      per_point = data.frame(boot_index = 1L, cluster = 1L, U = 0)
    ),
    class = "summary_bootstrap_2d_ld"
  )

  expect_error(compare_minima_depths(object, c(1, 2)), "Missing cluster")
  expect_error(compare_minima_depths(object, c(1, 1)), "two distinct")
})
