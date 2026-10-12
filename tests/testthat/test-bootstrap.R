make_test_ld <- function(U_mat) {
  x_coords <- seq(-1, 1, length.out = nrow(U_mat))
  y_coords <- seq(-1, 1, length.out = ncol(U_mat))

  d_mat <- exp(-(U_mat - min(U_mat)))
  d_mat <- d_mat / sum(d_mat)

  dist <- expand.grid(x = x_coords, y = y_coords)
  dist$d <- as.vector(d_mat)
  dist$U <- as.vector(U_mat)

  ss <- d_mat
  attr(ss, "x_coords") <- x_coords
  attr(ss, "y_coords") <- y_coords

  structure(
    list(
      dist = dist,
      plot = NULL,
      plot_2 = NULL,
      vf = list(x = "x", y = "y"),
      ss = ss
    ),
    class = c("2d_static_ld", "2d_ld", "landscape")
  )
}

test_that("moving-block bootstraps retain the original number of vectors", {
  data <- data.frame(
    x = sin(seq_len(30) / 3),
    y = cos(seq_len(30) / 5)
  )
  fields <- list(
    fit_1d_vf(data, x = "x", n = 5, h = 0.5),
    fit_2d_vf(data, x = "x", y = "y", n = 3, h = 0.5)
  )
  bootstrappers <- list(bootstrap_1d_vf, bootstrap_2d_vf)

  for (dimension in seq_along(fields)) {
    vf <- fields[[dimension]]
    n_vectors <- nrow(vf$original_vectors)
    for (block_length in 4L) {
      boot <- bootstrappers[[dimension]](
        vf, block_length = block_length, n_boot = 2,
        add_linear_interp = FALSE
      )
      for (model in boot$bootstrap_models) {
        expect_equal(nrow(model$original_vectors), n_vectors)
        expect_equal(nrow(model$original_vectors_normalized), n_vectors)
      }
    }
  }
})

test_that("Hungarian summaries use matched bootstrap means and covariance", {
  x_coords <- y_coords <- seq(-1, 1, length.out = 9)
  original_U <- outer(x_coords, y_coords, function(x, y) (x^2 - 0.25)^2 + y^2)
  centers <- data.frame(x = c(0.5, 0.75, 0.75, 1), y = c(0, 0.25, 0, -0.25))
  boot_obj <- structure(
    list(
      bootstrap_lds = lapply(seq_len(nrow(centers)), function(i) {
        make_test_ld(outer(x_coords, y_coords, function(x, y) {
          (x - centers$x[i])^2 + (y - centers$y[i])^2
        }))
      }),
      original_ld = make_test_ld(original_U),
      n_boot = nrow(centers)
    ),
    class = "bootstrap_2d_ld"
  )

  out <- summary(boot_obj, minima_method = "hungarian", exclude_minor = FALSE)
  minimum <- out$per_minimum
  expect_equal(nrow(minimum), 1L)
  expect_equal(minimum$minimum, 2L)
  expect_equal(minimum$mean_x, 0.75)
  expect_equal(minimum$mean_y, 0)
  expect_equal(mean(out$bootstrap_minima$x), 0.75)
  expect_equal(minimum$s_xx, var(centers$x))
  expect_equal(minimum$s_yy, var(centers$y))
  expect_equal(minimum$s_xy, cov(centers$x, centers$y))
  eigenvalues <- eigen(cov(centers), symmetric = TRUE)$values
  expect_equal(c(minimum$a_pred, minimum$b_pred),
               sqrt(eigenvalues * qchisq(0.95, df = 2)))
  expect_equal(minimum$n_runs, 4L)

  mean_potential <- summary(boot_obj, minima_method = "mean_potential", exclude_minor = FALSE)
  expect_equal(mean_potential$per_minimum$mean_x, 0.75)
  expect_equal(mean_potential$per_minimum$mean_y, 0)
  expect_equal(mean_potential$per_minimum$s_xx, minimum$s_xx)
  expect_equal(mean_potential$per_minimum$s_yy, minimum$s_yy)
  expect_equal(mean_potential$per_minimum$s_xy, minimum$s_xy)
})

test_that("compare_minima_depths uses paired bootstrap contrasts", {
  object <- structure(
    list(
      params = list(level = 0.95),
      per_minimum = data.frame(
        minimum = c(1L, 2L),
        mean_U = c(2, 0.5)
      ),
      bootstrap_minima = data.frame(
        boot_index = rep(1:3, each = 2),
        minimum = rep(1:2, 3),
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
  expect_error(compare_minima_depths(object, c(1, 3)), "Missing minimum")
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
