#' Estimate a 2D vector field
#'
#' Estimate a 2D vector field from intensive longitudinal data. Two methods can be used: Sparse Vector Field Consensus (SparseVFC, using [SparseVFC::SparseVFC()]), or Multivariate Vector Field Kernel Estimator (MVKE, using [MVKE()]). Note that the input data are automatically normalized before being sent to the estimation engines to make sure the default parameter settings are close to the optimal. Therefore, you do not need to scale up or down the parameters of [SparseVFC::SparseVFC()] or [MVKE()].
#'
#' @param data The data set used for estimating the vector field.
#' Should be a data frame or a matrix.
#' @param x,y Characters to indicate the name of the two variables.
#' @param lims The limits of the range for the vector field estimation as `c(<xl>, <xu>, <yl>, <yu>)`. If missing, the range of the data extended by 10% for both sides will be used.
#' @param n The number of equally spaced points in each axis, at which the vectors are to be estimated.
#' @param vector_position Only useful if `method == "VFC"`. One of "start", "middle", or "end", representing the position of the vectors. If "start", for example, the starting point of a vector is regarded as the position of the vector.
#' @param na_action One of "omit_data_points" or "omit_vectors". If using "omit_data_points", then only the `NA` points are omitted, and the points before and after an `NA` will form a vector. If using "omit_vectors", then the vectors will be omitted if either of its points is `NA`.
#' @param method One of "VFC" or "MVKE".
#' @param ... Other parameters to be passed to [SparseVFC::SparseVFC()] or [MVKE()].
#'
#' @return A `vectorfield` object.
#' @seealso [plot.vectorfield()]
#' @export
fit_2d_vf <- function(data, x, y,
                      lims,
                      n = 20,
                      vector_position = "start",
                      na_action = "omit_data_points",
                      method = c("VFC", "MVKE"), ...) {
  d <- data
  # extract useful data for construction
  if (is.data.frame(d)) {
    d_raw <- d[, c(x, y)] %>% as.matrix()
  } else if (is.matrix(d)) {
    d_raw <- d[, c(x, y)]
  } else {
    rlang::abort("`d` must be a data frame or a matrix.")
  }

  if (na_action != "omit_data_points" & na_action != "omit_vectors") {
    rlang::abort('`na_action` must be either "omit_data_points" or "omit_vectors".')
  }

  dv <- normalize_vecs(d_raw)

  if (any(is.na(dv))) {{ if (na_action == "omit_data_points") {
    dv <- attr(dv, "x_noNA")
    rlang::inform("NA(s) found in the data. Those data points were omitted.")
  } }}

  v_mat <- diff(dv)

  if (vector_position == "start") {
    x_mat <- dv[1:(nrow(dv) - 1), ]
  } else if (vector_position == "middle") {
    x_mat <- dv[1:(nrow(dv) - 1), ] + 0.5 * v_mat
  } else if (vector_position == "end") {
    x_mat <- dv[2:nrow(dv), ]
  } else {
    rlang::abort('`vector_position` must be one of "start", "middle", or "end".')
  }

  original_vectors_normalized <- cbind(x_mat, v_mat) %>%
    `colnames<-`(c("x", "y", "vx", "vy"))

  if (any(is.na(original_vectors_normalized))) {{ if (na_action == "omit_vectors") {
    original_vectors_normalized <- stats::na.omit(original_vectors_normalized)
    rlang::inform("NA(s) found in the data. Those vectors were omitted.")
  } }}

  original_vectors <- original_vectors_normalized
  original_vectors[, "x"] <- original_vectors[, "x"] %>% denormalize_x(dv)
  original_vectors[, "y"] <- original_vectors[, "y"] %>% denormalize_y(dv)
  original_vectors[, "vx"] <- original_vectors[, "vx"] %>% scale_up(dv)
  original_vectors[, "vy"] <- original_vectors[, "vy"] %>% scale_up(dv)

  VFCresult <- MVKEresult <- NULL
  method <- toupper(method[1])
  if (method == "VFC") {
    VFCresult <- SparseVFC::SparseVFC(original_vectors_normalized[, 1:2], original_vectors_normalized[, 3:4], ...)
  } else if (method == "MVKE") {
    MVKEresult <- MVKE(original_vectors_normalized[, 1:2], ...)
  }

  lims <- determine_lims(d_raw, c(x, y), lims)

  vec <- tidyr::expand_grid(x = seq(lims[1], lims[2], length.out = n), y = seq(lims[3], lims[4], length.out = n))
  vec <- vec %>% dplyr::rowwise()

  if (method == "VFC") {
    vec <- vec %>% dplyr::mutate(v = list(stats::predict(VFCresult, c(x, y) %>% normalize_v(dv)) %>% scale_up(dv)))
  } else if (method == "MVKE") {
    vec <- vec %>% dplyr::mutate(v = list(MVKEresult(c(x, y) %>% normalize_v(dv))$mu %>% scale_up(dv)))
  }

  vec <- vec %>%
    dplyr::mutate(
      vx = v[1],
      vy = v[2],
      v_norm = (sum(v^2))^(1 / 2)
    ) %>%
    dplyr::select(-v) %>%
    dplyr::ungroup()

  result <- list(
    vec_grid = vec,
    VFCresult = VFCresult,
    MVKEresult = MVKEresult,
    data = d_raw,
    data_normalized = dv,
    original_vectors = original_vectors,
    original_vectors_normalized = original_vectors_normalized,
    x = x,
    y = y,
    lims = lims,
    n = n,
    method = method
  )

  class(result) <- "vectorfield"

  return(result)
}

# borrowed from simlandr 0.3.0
determine_lims <- function(output, var_names, lims) {
  if (!rlang::is_missing(lims)) {
    return(lims)
  }
  if (is.list(output)) output <- output[[1]]
  if (rlang::is_missing(lims)) {
    return(c(sapply(var_names, function(v) grDevices::extendrange(output[, v], f = 0.1))))
  }
  if (any(is.infinite(lims))) stop("Non-infinite values found in `lims`.")
}
