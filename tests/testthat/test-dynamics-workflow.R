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
    evaluate_landscape_clusters = function(landscapes, ...) {
      structure(list(metrics = data.frame(k = 1:2)), class = "landscape_cluster_evaluation")
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
  expect_equal(names(result$members), c("first", "second"))
  expect_true(all(vapply(result$members, inherits, logical(1), "individual_dynamics")))
  expect_equal(result$lims, c(-1, 11, -2, 22))
  expect_equal(calls[[1]]$lims, result$lims)
  expect_equal(calls[[2]]$lims, result$lims)
  expect_null(result$clustering)
})

test_that("add_group_clusters attaches a selected solution", {
  object <- structure(
    list(
      landscapes = list(a = 1, b = 2),
      settings = list(cluster_method = "kmeans", nstart = 25L, iter.max = 100L, seed = 4L),
      clustering = NULL
    ),
    class = "group_dynamics"
  )

  local_mocked_bindings(
    cluster_landscapes = function(landscapes, k, method, nstart, iter.max, seed) {
      structure(list(k = k, cluster = c(1L, 2L)), class = "landscape_clusters")
    },
    .package = "fitlandr"
  )

  result <- add_group_clusters(object, k = 2)
  expect_s3_class(result$clustering, "landscape_clusters")
  expect_equal(result$clustering$k, 2)
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
