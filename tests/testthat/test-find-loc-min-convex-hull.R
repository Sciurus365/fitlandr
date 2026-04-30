test_that("find_loc_min marks outside-hull minima as minor", {
  x_vals <- seq(1, 7)
  y_vals <- seq(1, 7)
  dist <- expand.grid(x = x_vals, y = y_vals)
  dist$U <- 10
  dist$d <- exp(-dist$U)

  dist$U[dist$x == 4 & dist$y == 4] <- 0
  dist$U[dist$x == 2 & dist$y == 2] <- 0.5

  ld <- structure(
    list(
      dist = dist,
      vf = list(
        data = rbind(
          c(3, 3),
          c(5, 3),
          c(5, 5),
          c(3, 5)
        )
      )
    ),
    class = c("2d_static_ld", "2d_ld", "landscape")
  )

  mins <- find_loc_min(ld, exclude_minor = FALSE)$mins

  inside_idx <- which(mins$x == 4 & mins$y == 4)
  outside_idx <- which(mins$x == 2 & mins$y == 2)

  expect_length(inside_idx, 1)
  expect_length(outside_idx, 1)
  expect_false(mins$is_minor[inside_idx])
  expect_true(mins$is_minor[outside_idx])
})
