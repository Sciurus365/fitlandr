#' @title 10-Fold Cross-Validation for fit_2d_vf Bandwidth Selection
#'
#' @description Performs 10-fold cross-validation to find the optimal bandwidth 'h'
#'              by minimizing the MSE of the predicted drift (v) against the
#'              empirical drift of the test set (X_t+1 - X_t).
#'              Handles single-input prediction and NA insertion.
#'
#' @param data A matrix/data.frame with columns x and y representing the STATE coordinates (X_t, Y_t).
#' @param x The name of the column representing the X coordinate.
#' @param y The name of the column representing the Y coordinate.
#' @param h_values A numeric vector of candidate bandwidths to test.
#' @param k The number of folds for cross-validation (default is 10).
#' @param ... Additional arguments passed to fit_2d_vf (e.g., method, lims).
#' @return A list containing:
#'         - 'cv_results': A data frame of 'h' values and their corresponding 'cv_mse'.
#'         - 'h_optimal': The bandwidth that yielded the minimum 'cv_mse'.
#'
#' @export
cv_fit_2d_vf <- function(data, x, y, h_values = exp(seq(log(0.01), log(2), length.out = 20)), k = 10, ...) {
  # 1. Setup and Initialization
  if (!all(c(x, y) %in% colnames(data))) {
    cli::cli_abort("Data must contain columns named '{x}' and '{y}'.")
  }

  n <- nrow(data)
  data <- as.data.frame(data)
  folds <- sample(rep(1:k, length.out = n)) # Assign data points to folds

  mse_by_h <- numeric(length(h_values))

  cli::cli_inform("Starting {k}-Fold CV for {length(h_values)} candidate bandwidths.")

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

    cli::cli_inform("Testing h = {h_curr}")

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
      test_set_start_indices <- test_indices[test_indices < n]

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
        suppressMessages(fit_2d_vf(train_data, x = x, y = y, h = h_curr, ...)),
        error = function(e) {
          cli::cli_alert("Error in fit_2d_vf for fold {i}, h = {h_curr}: {e$message}")
          return(NULL)
        }
      )

      if (is.null(model_train)) {
        fold_mse[i] <- NA
        next
      }

      # Predict the drift 'v' at the test states X_t individually
      predicted_v_list <- apply(X_t, 1, function(row) {
        # 'row' is a two-element vector c(x, y), which is passed to predict()
        pred <- stats::predict(model_train, c(row[1], row[2]))
        return(pred$v)
      })

      # Transpose the list result to an N x 2 matrix for comparison
      predicted_drift_matrix <- t(predicted_v_list)

      # --- D. Calculate Mean Squared Error (MSE) ---

      # The error is calculated between the true empirical drift and the predicted drift
      drift_error <- true_empirical_drift - predicted_drift_matrix

      # MSE: sum of squares across the two dimensions (x and y component of drift)
      sq_error <- rowSums(drift_error^2)

      fold_mse[i] <- mean(sq_error, na.rm = TRUE)
    }

    # 4. Store Average MSE for the current bandwidth h
    mse_by_h[j] <- mean(fold_mse, na.rm = TRUE)
  }

  # 5. Determine Optimal Bandwidth
  cv_results <- data.frame(h = h_values, cv_mse = mse_by_h)
  h_optimal <- cv_results$h[which.min(cv_results$cv_mse)]

  cli::cli_inform("Cross-Validation complete.")
  cli::cli_inform("Optimal Bandwidth (h) selected: {round(h_optimal, 4)}")

  cli::cli_inform("Fitting the optimal model on the full dataset...")
  final_model <- fit_2d_vf(data, x = x, y = y, h = h_optimal, ...)
  cli::cli_inform("Final model fitted.")

  return(structure(list(
    final_model = final_model,
    cv_results = cv_results,
    h_optimal = h_optimal,
    data = data,
    folds = folds,
    x = x,
    y = y
  ), class = "cv_vectorfield"))
}

#' @rdname cv_fit_2d_vf
#' @export
#'
#' @param x An object of class 'cv_vectorfield' returned by cv_fit_2d_vf.
#' @param ... Additional arguments (not used).
plot.cv_vectorfield <- function(x, ...) {
  cv_data <- x$cv_results
  ggplot2::ggplot(cv_data, ggplot2::aes(x = h, y = cv_mse)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::annotate(
      x = x$h_optimal, y = min(cv_data$cv_mse),
      geom = "point", color = "red"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "Bandwidth (h, log scale)",
      y = "Cross-Validation MSE"
    ) +
    ggplot2::theme_bw()
}
