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


test_that("mean_potential_hessian summary returns compatible ellipse output", {
  x_coords <- seq(-1, 1, length.out = 5)
  y_coords <- seq(-1, 1, length.out = 5)
  U_base <- outer(x_coords, y_coords, function(x, y) x^2 + y^2)

  boot_obj <- structure(
    list(
      bootstrap_lds = list(
        make_test_ld(U_base),
        make_test_ld(U_base + 0.05),
        make_test_ld(U_base * 1.10)
      ),
      original_ld = make_test_ld(U_base),
      n_boot = 3L
    ),
    class = "bootstrap_2d_ld"
  )

  out <- summary(
    boot_obj,
    clustering_method = "mean_potential_hessian",
    exclude_minor = TRUE,
    min_barrier = 0.05,
    level = 0.95
  )

  expect_s3_class(out, "summary_bootstrap_2d_ld")
  expect_equal(out$params$clustering_method, "mean_potential_hessian")
  expect_equal(nrow(out$per_boot), 3L)
  expect_equal(out$per_boot$n_mins, c(1L, 1L, 1L))
  expect_true(nrow(out$per_cluster) >= 1L)
  expect_true(all(c("a_pred", "b_pred", "a_conf", "b_conf", "h_xx", "h_xy", "h_yy") %in% names(out$per_cluster)))
  expect_true(any(is.finite(out$per_cluster$a_pred)))
  expect_true(any(out$per_cluster$hessian_ok))
})


test_that("mean_potential_hessian excludes bootstrap runs without matched minima", {
  x_coords <- seq(-1, 1, length.out = 5)
  y_coords <- seq(-1, 1, length.out = 5)
  U_base <- outer(x_coords, y_coords, function(x, y) x^2 + y^2)
  U_flat <- matrix(1, nrow = 5, ncol = 5)

  boot_obj <- structure(
    list(
      bootstrap_lds = list(
        make_test_ld(U_base),
        make_test_ld(U_flat),
        make_test_ld(U_base)
      ),
      original_ld = make_test_ld(U_base),
      n_boot = 3L
    ),
    class = "bootstrap_2d_ld"
  )

  out <- summary(
    boot_obj,
    clustering_method = "mean_potential_hessian",
    exclude_minor = TRUE,
    min_barrier = 0.05,
    level = 0.95
  )

  expect_equal(out$per_boot$n_mins, c(1L, 0L, 1L))
  expect_equal(out$per_cluster$n_runs[[1]], 2L)
  expect_equal(out$per_cluster$n_hessian[[1]], 2L)
  expect_true(out$per_cluster$stability[[1]] < 1)
})