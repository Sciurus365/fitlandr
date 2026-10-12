make_minima_test_landscape <- function(minima, grid = 1:7, observed = c(3, 5)) {
  dist <- expand.grid(x = grid, y = grid)
  dist$U <- 10
  for (i in seq_len(nrow(minima))) {
    dist$U[dist$x == minima[i, 1] & dist$y == minima[i, 2]] <- minima[i, 3]
  }
  dist$d <- exp(-dist$U)
  structure(list(
    dist = dist,
    vf = list(
      data = expand.grid(x = observed, y = observed),
      x = "x", y = "y"
    )
  ), class = c("2d_static_ld", "2d_ld", "landscape"))
}

test_that("convex-hull exclusion precedes barrier comparisons", {
  ld <- make_minima_test_landscape(
    rbind(c(3, 3, 0), c(5, 5, 1), c(8, 8, -10)),
    grid = 1:9, observed = c(2, 6)
  )
  result <- find_loc_min(ld, min_barrier_fraction = 0.6,
                         min_convex_hull_range_fraction = 0)
  expect_equal(result$mins$x, c(3, 5, 8))
  expect_equal(result$mins$is_minor, c(FALSE, FALSE, TRUE))
  expect_true(all(is.na(result$barriers[3, ])))
  expect_true(all(is.na(result$barriers[, 3])))
})

test_that("exclude_minor = FALSE disables all minor-minimum calculations", {
  ld <- make_minima_test_landscape(rbind(c(4, 4, 0), c(2, 2, 0.5)))
  result <- find_loc_min(ld, exclude_minor = FALSE)
  expect_true(all(!result$mins$is_minor))
  expect_true(all(result$mins$exclusion_reason == "retained_major"))
  expect_true(all(is.na(result$barriers)))
})

test_that("a sole outside-hull minimum is excluded unless explicitly retained", {
  ld <- make_minima_test_landscape(rbind(c(2, 2, 0)))
  excluded <- find_loc_min(ld)$mins
  retained <- find_loc_min(ld, retain_single_outside_hull = TRUE)$mins
  expect_true(excluded$is_minor[[1]])
  expect_identical(excluded$exclusion_reason[[1]], "outside_observed_data_convex_hull")
  expect_false(retained$is_minor[[1]])
  expect_identical(retained$exclusion_reason[[1]],
                   "retained_major_outside_observed_data_convex_hull")
})
