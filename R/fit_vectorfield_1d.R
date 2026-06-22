#' Estimate a 1D vector field
#'
#' Estimate a one-dimensional vector field from intensive longitudinal data.
#' The interface mirrors [fit_2d_vf()], but for a single state variable.
#'
#' @param data The data set used for estimating the vector field.
#' @param x Character naming the state variable.
#' @param lims The limits of the range for vector-field estimation as
#'   `c(xl, xu)`. If missing, the data range extended by 10% is used.
#' @param n Number of equally spaced grid points at which the drift is
#'   estimated.
#' @param vector_position One of `"start"`, `"middle"`, or `"end"`.
#' @param na_action One of `"omit_data_points"` or `"omit_vectors"`.
#' @param dayvar Optional character scalar naming the day variable.
#' @param beepvar Optional character scalar naming the within-day assessment
#'   order variable.
#' @param method One of `"MVKE"` or `"VFC"`.
#' @param ... Additional arguments passed to [MVKE()] or
#'   [SparseVFC::SparseVFC()].
#'
#' @return A `1d_vectorfield` object.
#' @export
fit_1d_vf <- function(data, x,
                      lims,
                      n = 200L,
                      vector_position = "start",
                      na_action = "omit_data_points",
                      dayvar = NULL,
                      beepvar = NULL,
                      method = c("MVKE", "VFC"), ...) {
  d <- data
  if (!is.null(dayvar) && !dayvar %in% colnames(d)) {
    cli::cli_abort("{.arg dayvar} must name a column in {.arg data}.")
  }
  if (!is.null(beepvar) && !beepvar %in% colnames(d)) {
    cli::cli_abort("{.arg beepvar} must name a column in {.arg data}.")
  }

  d_sep <- insert_time_separators(d, x, dayvar = dayvar, beepvar = beepvar)
  warn_if_timevars_ineffective(dayvar = dayvar, beepvar = beepvar, na_action = na_action)

  if (is.data.frame(d)) {
    d_input_raw <- d[, x, drop = FALSE] %>% as.matrix()
  } else if (is.matrix(d)) {
    d_input_raw <- d[, x, drop = FALSE]
  } else {
    cli::cli_abort("{.arg data} must be a data frame or a matrix.")
  }

  if (is.data.frame(d_sep)) {
    d_raw <- d_sep[, x, drop = FALSE] %>% as.matrix()
  } else if (is.matrix(d_sep)) {
    d_raw <- d_sep[, x, drop = FALSE]
  } else {
    cli::cli_abort("{.arg data} must be a data frame or a matrix.")
  }

  if (na_action != "omit_data_points" && na_action != "omit_vectors") {
    cli::cli_abort('{.arg na_action} must be either "omit_data_points" or "omit_vectors".')
  }

  dv <- normalize_vecs(d_raw)
  if (any(is.na(dv)) && na_action == "omit_data_points") {
    dv <- attr(dv, "x_noNA")
    cli::cli_inform("NA(s) found in the data. Those data points were omitted.")
  }

  v_mat <- diff(dv[, 1])
  if (vector_position == "start") {
    x_mat <- dv[1:(nrow(dv) - 1), 1]
  } else if (vector_position == "middle") {
    x_mat <- dv[1:(nrow(dv) - 1), 1] + 0.5 * v_mat
  } else if (vector_position == "end") {
    x_mat <- dv[2:nrow(dv), 1]
  } else {
    cli::cli_abort('{.arg vector_position} must be one of "start", "middle", or "end".')
  }

  original_vectors_normalized <- cbind(x = x_mat, vx = v_mat)
  if (any(is.na(original_vectors_normalized)) && na_action == "omit_vectors") {
    original_vectors_normalized <- stats::na.omit(original_vectors_normalized)
    cli::cli_inform("NA(s) found in the data. Those vectors were omitted.")
  }

  original_vectors <- original_vectors_normalized
  original_vectors[, "x"] <- denormalize_x(original_vectors[, "x"], dv)
  original_vectors[, "vx"] <- scale_up(original_vectors[, "vx"], dv)

  VFCresult <- MVKEresult <- NULL
  method <- toupper(method[1])
  if (method == "VFC") {
    VFCresult <- SparseVFC::SparseVFC(
      as.matrix(original_vectors_normalized[, "x", drop = FALSE]),
      as.matrix(original_vectors_normalized[, "vx", drop = FALSE]),
      ...
    )
  } else if (method == "MVKE") {
    MVKEresult <- MVKE(
      as.matrix(original_vectors_normalized[, "x", drop = FALSE]),
      as.matrix(original_vectors_normalized[, "vx", drop = FALSE]),
      ...
    )
  } else {
    cli::cli_abort('{.arg method} must be one of "MVKE" or "VFC".')
  }

  lims <- simlandr::determine_lims(d_input_raw, x, lims)
  x_grid <- seq(lims[1], lims[2], length.out = n)
  vx_grid <- numeric(length(x_grid))

  for (i in seq_along(x_grid)) {
    pos_norm <- normalize_x(x_grid[i], dv)
    if (method == "VFC") {
      vx_grid[i] <- as.numeric(stats::predict(VFCresult, matrix(pos_norm, ncol = 1))) %>% scale_up(dv)
    } else {
      vx_grid[i] <- as.numeric(MVKEresult(pos_norm)$mu) %>% scale_up(dv)
    }
  }

  vec_grid <- data.frame(x = x_grid, vx = vx_grid, v_norm = abs(vx_grid))

  result <- list(
    vec_grid = vec_grid,
    VFCresult = VFCresult,
    MVKEresult = MVKEresult,
    data = d_input_raw,
    data_normalized = dv,
    original_vectors = original_vectors,
    original_vectors_normalized = original_vectors_normalized,
    x = x,
    y = NULL,
    lims = lims,
    n = n,
    method = method
  )
  class(result) <- c("1d_vectorfield", "vectorfield")
  result
}

#' 10-fold cross-validation for 1D vector fields
#'
#' @inheritParams fit_1d_vf
#' @param h_values Numeric vector of candidate bandwidths.
#' @param k Number of folds.
#'
#' @return A `cv_vectorfield` object.
#' @export
cv_fit_1d_vf <- function(data, x, dayvar = NULL, beepvar = NULL,
                         h_values = exp(seq(log(0.01), log(2), length.out = 20)),
                         k = 10, ...) {
  if (!x %in% colnames(data)) {
    cli::cli_abort("Data must contain a column named '{x}'.")
  }
  if (!is.null(dayvar) && !dayvar %in% colnames(data)) {
    cli::cli_abort("Data must contain the day variable column named '{dayvar}'.")
  }
  if (!is.null(beepvar) && !beepvar %in% colnames(data)) {
    cli::cli_abort("Data must contain the beep variable column named '{beepvar}'.")
  }

  n <- nrow(data)
  data <- as.data.frame(data)
  fold_sizes <- rep.int(floor(n / k), k)
  if (n %% k > 0) {
    fold_sizes[seq_len(n %% k)] <- fold_sizes[seq_len(n %% k)] + 1L
  }
  folds <- rep.int(seq_len(k), times = fold_sizes)
  verbose <- isTRUE(getOption("fitlandr.verbose", TRUE))
  mse_by_h <- numeric(length(h_values))

  prepare_train_data <- function(indices, full_data, folds) {
    train_indices <- which(folds %in% indices)
    train_indices_sorted <- sort(train_indices)
    breaks <- c(0, diff(train_indices_sorted)) > 1

    train_data_list <- list()
    start_idx <- 1L
    for (i in seq_along(train_indices_sorted)) {
      if (breaks[i]) {
        train_data_list[[length(train_data_list) + 1L]] <- full_data[train_indices_sorted[start_idx:(i - 1L)], , drop = FALSE]
        na_row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(full_data)))
        colnames(na_row) <- colnames(full_data)
        train_data_list[[length(train_data_list) + 1L]] <- na_row
        start_idx <- i
      }
    }
    train_data_list[[length(train_data_list) + 1L]] <- full_data[train_indices_sorted[start_idx:length(train_indices_sorted)], , drop = FALSE]
    do.call(rbind, train_data_list)
  }

  is_valid_transition <- function(full_data, idx) {
    valid <- idx < nrow(full_data)
    if (!is.null(dayvar)) {
      valid <- valid & (full_data[[dayvar]][idx] == full_data[[dayvar]][idx + 1L])
    }
    if (!is.null(beepvar)) {
      valid <- valid & ((full_data[[beepvar]][idx + 1L] - full_data[[beepvar]][idx]) == 1)
    }
    valid
  }

  for (j in seq_along(h_values)) {
    h_curr <- h_values[j]
    fold_mse <- numeric(k)

    for (i in seq_len(k)) {
      test_indices <- which(folds == i)
      train_data <- prepare_train_data(setdiff(seq_len(k), i), data, folds)
      test_set_start_indices <- test_indices[is_valid_transition(data, test_indices)]
      if (!length(test_set_start_indices)) {
        fold_mse[i] <- NA_real_
        next
      }

      X_t <- data[test_set_start_indices, x, drop = TRUE]
      X_t_plus_1 <- data[test_set_start_indices + 1L, x, drop = TRUE]
      true_empirical_drift <- X_t_plus_1 - X_t

      model_train <- tryCatch(
        suppressMessages(fit_1d_vf(train_data, x = x, h = h_curr, dayvar = dayvar, beepvar = beepvar, ...)),
        error = function(e) {
          if (verbose) {
            cli::cli_alert("Error in fit_1d_vf for fold {i}, h = {h_curr}: {e$message}")
          }
          NULL
        }
      )

      if (is.null(model_train)) {
        fold_mse[i] <- NA_real_
        next
      }

      predicted_drift <- vapply(X_t, function(xx) as.numeric(stats::predict(model_train, xx)$v), numeric(1))
      sq_error <- (true_empirical_drift - predicted_drift)^2
      fold_mse[i] <- mean(sq_error, na.rm = TRUE)
    }

    mse_by_h[j] <- mean(fold_mse, na.rm = TRUE)
  }

  cv_results <- data.frame(h = h_values, cv_mse = mse_by_h)
  h_optimal <- cv_results$h[which.min(cv_results$cv_mse)]
  final_model <- fit_1d_vf(data, x = x, h = h_optimal, dayvar = dayvar, beepvar = beepvar, ...)

  structure(list(
    final_model = final_model,
    cv_results = cv_results,
    h_optimal = h_optimal,
    data = data,
    folds = folds,
    x = x,
    y = NULL,
    dayvar = dayvar,
    beepvar = beepvar
  ), class = c("cv_1d_vectorfield", "cv_vectorfield"))
}

