#' @title 10-Fold Cross-Validation for fit_2d_vf Bandwidth Selection
#'
#' @description Performs 10-fold cross-validation to find the optimal bandwidth 'h'
#'              by minimizing the MSE of the predicted drift (v) against the
#'              empirical drift of the test set (X_t+1 - X_t).
#'              Folds are assigned as contiguous blocks to respect temporal order.
#'              Handles single-input prediction and NA insertion.
#'
#' @param data A matrix/data.frame with columns x and y representing the STATE coordinates (X_t, Y_t).
#' @param x The name of the column representing the X coordinate.
#' @param y The name of the column representing the Y coordinate.
#' @param dayvar Optional character scalar naming the day variable. When
#'   supplied, overnight transitions are excluded from both training and
#'   validation.
#' @param beepvar Optional character scalar naming the within-day assessment
#'   order variable. When supplied, only consecutive beeps are used in
#'   training and validation.
#' @param h_values A numeric vector of candidate bandwidths to test.
#' @param k The number of folds for cross-validation (default is 10).
#' @param n Number of vector-field grid points per axis, passed to every
#'   [fit_2d_vf()] call.
#' @param na_action Missing-data handling for the final full-data fit. Fold
#'   fits always use `"omit_vectors"` so vectors cannot bridge across held-out
#'   blocks. Defaults to `"omit_vectors"`, which is also required for
#'   `dayvar` and `beepvar` separators to remain effective in the final fit.
#' @param ... Additional arguments passed to [fit_2d_vf()] for every
#'   fold-specific fit and the final full-data fit. These include `lims`,
#'   `vector_position`, and `method`.
#' @return A list containing:
#'         - 'cv_results': A data frame of 'h' values and their corresponding 'cv_mse'.
#'         - 'h_optimal': The bandwidth that yielded the minimum 'cv_mse'.
#'
#' @export
cv_fit_2d_vf <- function(data, x, y, dayvar = NULL, beepvar = NULL, h_values = exp(seq(log(0.01), log(2), length.out = 20)), k = 10, n = 20, na_action = "omit_vectors", ...) {
  # 1. Setup and Initialization
  if (!all(c(x, y) %in% colnames(data))) {
    cli::cli_abort("Data must contain columns named '{x}' and '{y}'.")
  }
  if (!is.null(dayvar) && !dayvar %in% colnames(data)) {
    cli::cli_abort("Data must contain the day variable column named '{dayvar}'.")
  }
  if (!is.null(beepvar) && !beepvar %in% colnames(data)) {
    cli::cli_abort("Data must contain the beep variable column named '{beepvar}'.")
  }

  n_obs <- nrow(data)
  data <- as.data.frame(data)
  fold_sizes <- rep.int(floor(n_obs / k), k)
  if (n_obs %% k > 0) {
    fold_sizes[seq_len(n_obs %% k)] <- fold_sizes[seq_len(n_obs %% k)] + 1L
  }
  folds <- rep.int(seq_len(k), times = fold_sizes)
  verbose <- isTRUE(getOption("fitlandr.verbose", TRUE))

  is_valid_transition <- function(full_data, idx) {
    valid <- idx < nrow(full_data)
    if (!is.null(dayvar)) {
      valid <- valid & (full_data[[dayvar]][idx] == full_data[[dayvar]][idx + 1L])
    }
    if (!is.null(beepvar)) {
      valid <- valid & ((full_data[[beepvar]][idx + 1L] - full_data[[beepvar]][idx]) == 1)
    }
    valid[is.na(valid)] <- FALSE
    valid
  }

  mse_by_h <- numeric(length(h_values))

  if (verbose) {
    cli::cli_inform("Starting {k}-Fold CV for {length(h_values)} candidate bandwidths.")
  }

  # --- Helper function to prepare training data (NA insertion) ---
  # Inserts NA rows between non-consecutive segments to prevent artifact vectors
  prepare_train_data <- function(indices, full_data, folds) {
    train_indices <- which(folds %in% indices)

    # Sort indices to identify contiguous segments
    train_indices_sorted <- sort(train_indices)

    # Find breaks (non-consecutive indices)
    breaks <- c(0, diff(train_indices_sorted)) > 1

    train_data_list <- list()
    start_idx <- 1

    for (i in seq_along(train_indices_sorted)) {
      if (breaks[i]) {
        # If a break is found, add the previous segment
        segment <- full_data[train_indices_sorted[start_idx:(i - 1)], ]
        train_data_list[[length(train_data_list) + 1]] <- segment

        # Add an NA row (the separator)
        na_row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(full_data)))
        colnames(na_row) <- colnames(full_data)
        train_data_list[[length(train_data_list) + 1]] <- na_row

        start_idx <- i
      }
    }
    # Add the final segment (the remaining data)
    train_data_list[[length(train_data_list) + 1]] <- full_data[train_indices_sorted[start_idx:length(train_indices_sorted)], ]

    # Combine into one data frame
    return(do.call(rbind, train_data_list))
  }

  # 2. Loop over Candidate Bandwidths
  for (j in seq_along(h_values)) {
    h_curr <- h_values[j]
    fold_mse <- numeric(k)

    if (verbose) {
      cli::cli_inform("Testing h = {h_curr}")
    }

    # 3. Loop over Folds
    for (i in 1:k) {
      test_indices <- which(folds == i)

      # --- A. Prepare Training Data ---
      train_indices_to_use <- setdiff(1:k, i)
      train_data <- prepare_train_data(train_indices_to_use, data, folds)

      # --- B. Define Test Set and Calculate True Empirical Drift ---

      # The test set consists of intervals (X_t, X_t+1) where both X_t and X_t+1 belong to the fold.
      # A better approach here is to define the test set based on the start point X_t.

      # We test on all points X_t in the fold *except* the very last point of the series (n)
      # if the last point of the series happens to be in the current fold.
      test_set_start_indices <- test_indices[is_valid_transition(data, test_indices)]

      if (length(test_set_start_indices) == 0) next

      # X_t is the state where prediction is made
      X_t <- data[test_set_start_indices, c(x, y)]
      # X_t_plus_1 is the state used for the true target
      X_t_plus_1 <- data[test_set_start_indices + 1, c(x, y)]

      # TRUE TARGET: The empirical drift (Y_t) = X_t+1 - X_t
      # This is the vector the predicted drift 'v' should match.
      true_empirical_drift <- X_t_plus_1 - X_t


      # --- C. Train Model and Predict ---
      model_train <- tryCatch(
        suppressMessages(fit_2d_vf(
          train_data,
          x = x,
          y = y,
          h = h_curr,
          dayvar = dayvar,
          beepvar = beepvar,
          n = n,
          na_action = "omit_vectors",
          ...
        )),
        error = function(e) {
          if (verbose) {
            cli::cli_alert("Error in fit_2d_vf for fold {i}, h = {h_curr}: {e$message}")
          }
          return(NULL)
        }
      )

      if (is.null(model_train)) {
        fold_mse[i] <- NA
        next
      }

      # Predict the drift 'v' at the test states X_t individually
      predicted_drift_matrix <- matrix(NA_real_, nrow = nrow(X_t), ncol = 2)
      for (idx in seq_len(nrow(X_t))) {
        pred <- stats::predict(model_train, c(X_t[idx, 1], X_t[idx, 2]))
        predicted_drift_matrix[idx, ] <- pred$v
      }

      # --- D. Calculate Mean Squared Error (MSE) ---

      # The error is calculated between the true empirical drift and the predicted drift
      drift_error <- true_empirical_drift - predicted_drift_matrix

      # MSE: sum of squares across the two dimensions (x and y component of drift)
      sq_error <- rowSums(drift_error^2)

      fold_mse[i] <- mean(sq_error, na.rm = TRUE)
    }

    # 4. Store Average MSE for the current bandwidth h
    finite_fold_mse <- fold_mse[is.finite(fold_mse)]
    mse_by_h[j] <- if (length(finite_fold_mse)) {
      mean(finite_fold_mse)
    } else {
      NA_real_
    }
  }

  # 5. Determine Optimal Bandwidth
  cv_results <- data.frame(h = h_values, cv_mse = mse_by_h)
  if (!any(is.finite(cv_results$cv_mse))) {
    cli::cli_abort(c(
      "Cross-validation failed for every candidate bandwidth.",
      "i" = "Check missing values, time-order variables, and whether each fold contains enough valid transitions."
    ))
  }
  h_optimal <- cv_results$h[which.min(cv_results$cv_mse)]

  if (verbose) {
    cli::cli_inform("Cross-Validation complete.")
    cli::cli_inform("Optimal Bandwidth (h) selected: {round(h_optimal, 4)}")
    cli::cli_inform("Fitting the optimal model on the full dataset...")
  }
  final_model <- fit_2d_vf(
    data,
    x = x,
    y = y,
    h = h_optimal,
    dayvar = dayvar,
    beepvar = beepvar,
    n = n,
    na_action = na_action,
    ...
  )
  if (verbose) {
    cli::cli_inform("Final model fitted.")
  }

  return(structure(list(
    final_model = final_model,
    cv_results = cv_results,
    h_optimal = h_optimal,
    data = data,
    folds = folds,
    x = x,
    y = y,
    dayvar = dayvar,
    beepvar = beepvar
  ), class = "cv_vectorfield"))
}

#' Autoplot 2D vector-field cross-validation results
#'
#' Draw the cross-validation error over the candidate bandwidth values and
#' mark the selected bandwidth.
#'
#' @param object An object of class 'cv_vectorfield' returned by cv_fit_2d_vf.
#' @param ... Additional arguments (not used).
#'
#' @return A ggplot object.
#' @export
autoplot.cv_vectorfield <- function(object, ...) {
  cv_data <- object$cv_results
  ggplot2::ggplot(cv_data, ggplot2::aes(x = h, y = cv_mse)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::annotate(
      x = object$h_optimal, y = min(cv_data$cv_mse),
      geom = "point", color = "red"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "Bandwidth (h, log scale)",
      y = "Cross-Validation MSE"
    ) +
    ggplot2::theme_bw()
}

#' Plot cross-validation results
#'
#' `r lifecycle::badge("deprecated")`
#'
#' `plot.cv_vectorfield()` is deprecated in favor of [autoplot()].
#'
#' @param x An object of class `cv_vectorfield` returned by [cv_fit_2d_vf()].
#' @param ... Arguments passed to [autoplot()].
#'
#' @export
plot.cv_vectorfield <- function(x, ...) {
  lifecycle::deprecate_warn(
    "0.2.0",
    "plot.cv_vectorfield()",
    "autoplot()"
  )
  autoplot(x, ...)
}
