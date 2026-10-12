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
    for (block_length in c(4L, 5L)) {
      boot <- bootstrappers[[dimension]](
        vf, block_length = block_length, n_boot = 3,
        add_linear_interp = FALSE
      )
      for (model in boot$bootstrap_models) {
        expect_equal(nrow(model$original_vectors), n_vectors)
        expect_equal(nrow(model$original_vectors_normalized), n_vectors)
      }
    }
  }
})
