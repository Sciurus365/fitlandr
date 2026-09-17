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

test_that("HDBSCAN returns no clusters for fewer pooled minima than minPts", {
  boot_min_df <- data.frame(
    x = c(-1, 0, 1, 2),
    y = c(0, 1, 0, 1),
    U = c(0.1, 0.2, 0.3, 0.4),
    boot_index = seq_len(4)
  )
  object <- list(
    original_ld = list(
      dist = expand.grid(x = 0:1, y = 0:1)
    )
  )

  out <- fitlandr:::cluster_bootstrap_minima(
    boot_min_df = boot_min_df,
    object = object,
    exclude_minor = TRUE,
    min_barrier_fraction = 0.1,
    min_convex_hull_range_fraction = 0.01,
    clustering_method = "hdbscan",
    minPts = 5
  )

  expect_true(all(out$boot_min_df$is_noise))
  expect_equal(out$boot_min_df$cluster, rep(0L, 4))
  expect_identical(
    out$diagnostics$clustering_status,
    "fewer_pooled_minima_than_minPts"
  )
})
