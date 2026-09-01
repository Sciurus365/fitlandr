#' Multivariate vector field kernel estimator
#'
#' See references for details.
#'
#' @param d The dataset. Should be a matrix or a data frame, with each row representing a random vector.
#' @param v The vectors corresponding to the dataset. Should be a matrix or a data frame with the same shape as `d`. If missing, then the vectors will be calculated from the dataset.
#' @param h The bandwidth for the kernel estimator.
#' @param kernel The type of kernel estimator used. "Gaussian" by default.
#' @param prior_mean_drift Logical. If `TRUE`, apply a weak prior drift that points
#' toward the mean position of `d`. The prior is only active in very low-support
#' regions and is negligible close to the data.
#' @param prior_strength Non-negative numeric scalar controlling the maximum strength
#' of the prior in extremely low-support regions.
#' @param prior_n_eff_scale Positive numeric scalar controlling how quickly the prior
#' decays as local effective sample size increases.
#' @param prior_distance_linear Logical. If `TRUE`, use a linear distance-scaled
#' prior value: the prior drift magnitude is proportional to
#' `distance_to_mean / median_distance_to_mean_data`.
#'
#' @return A function(x), which then returns the drift estimator \eqn{\mu} and
#'   the infinitesimal covariance estimator
#'   \eqn{a = \sigma \sigma^{\mathsf{T}}} at position \eqn{x}. The estimator
#'   assumes unit time intervals.
#' @references Bandi, F. M., & Moloche, G. (2018). On the functional estimation of multivariate diffusion processes. Econometric Theory, 34(4), 896-946. https://doi.org/10.1017/S0266466617000305
#' @export
MVKE <- function(d, v, h = 0.2, kernel = c("Gaussian", "exp"),
                 prior_mean_drift = FALSE,
                 prior_strength = 0.2,
                 prior_n_eff_scale = 2,
                 prior_distance_linear = FALSE) {
  if (is.data.frame(d)) d <- as.matrix(d)
  if (!is.matrix(d)) cli::cli_abort("{.arg d} should be a data frame or a matrix.")
  if (any(is.na(d))) cli::cli_abort("There are missing values in {.arg d}.")
  if (missing(v)) {
    v <- diff(d)
    d <- d[1:(nrow(d) - 1), , drop = FALSE]
  } else {
    if (is.data.frame(v)) v <- as.matrix(v)
    if (!is.matrix(v)) cli::cli_abort("{.arg v} should be a data frame or a matrix.")
    if (any(is.na(v))) cli::cli_abort("There are missing values in {.arg v}.")
    if (!all(dim(v) == dim(d))) cli::cli_abort("{.arg v} should have the same shape as {.arg d}.")
  }


  # d <- stats::na.omit(d)
  dim <- ncol(d)
  d_mean <- colMeans(d)
  d_centered <- sweep(d, 2, d_mean, "-")
  prior_dist_ref <- stats::median(sqrt(Rfast::rowsums(d_centered^2)))
  if (!is.finite(prior_dist_ref) || prior_dist_ref <= 0) prior_dist_ref <- 1

  temp_diff <- v
  prior_scale <- stats::median(sqrt(Rfast::rowsums(temp_diff^2)))
  if (!is.finite(prior_scale) || prior_scale <= 0) prior_scale <- 1
  kernel <- kernel[1]
  if (kernel == "Gaussian") {
    log_K <- log_K_gaussian_mat
  } else if (kernel == "exp") {
    log_K <- log_K_exp_mat
  } else {
    cli::cli_abort('{.arg kernel} must be one of "Gaussian" or "exp".')
  }
  if (!is.logical(prior_mean_drift) || length(prior_mean_drift) != 1) {
    cli::cli_abort("{.arg prior_mean_drift} must be TRUE or FALSE.")
  }
  if (!is.numeric(prior_strength) || length(prior_strength) != 1 || !is.finite(prior_strength) || prior_strength < 0) {
    cli::cli_abort("{.arg prior_strength} must be a non-negative finite number.")
  }
  if (!is.numeric(prior_n_eff_scale) || length(prior_n_eff_scale) != 1 || !is.finite(prior_n_eff_scale) || prior_n_eff_scale <= 0) {
    cli::cli_abort("{.arg prior_n_eff_scale} must be a positive finite number.")
  }
  if (!is.logical(prior_distance_linear) || length(prior_distance_linear) != 1) {
    cli::cli_abort("{.arg prior_distance_linear} must be TRUE or FALSE.")
  }

  force(h)
  function(x) {
    if (length(x) != dim) cli::cli_abort("Input has wrong dimension.")

    # Compute log-kernel weights once. The same weights are used in both
    # drift and diffusion estimators.
    log_w <- log_K(d, x, h = h)

    # Find a common constant to shift by (usually the max of the denominator weights)
    max_log <- max(log_w)

    # Shift and exponentiate: exp(log_w - max_log)
    # This brings the largest value to 1, others will be relative to it
    w_shifted <- exp(log_w - max_log)

    # The constant (exp(max_log)) cancels out in the numerator and denominator
    denom_sum <- sum(w_shifted)
    mu_data <- colSums(w_shifted * temp_diff) / denom_sum

    if (prior_mean_drift) {
      # Effective sample size is small in low-support regions.
      n_eff <- denom_sum^2 / sum(w_shifted^2)
      if (!is.finite(n_eff)) n_eff <- 0

      pull_vec <- d_mean - x
      pull_norm <- sqrt(sum(pull_vec^2))
      dist_factor <- if (prior_distance_linear) pull_norm / prior_dist_ref else 1

      # Prior weight depends only on local effective support.
      prior_weight <- prior_strength * exp(-n_eff / prior_n_eff_scale)

      if (pull_norm > 0) {
        # Prior value itself can scale linearly with distance to mean.
        prior_mu <- prior_scale * dist_factor * pull_vec / pull_norm
      } else {
        prior_mu <- rep(0, dim)
      }
      mu_est <- (mu_data + prior_weight * prior_mu) / (1 + prior_weight)
    } else {
      mu_est <- mu_data
    }

    # Sum_i w_i * (v_i %*% t(v_i)) using a weighted crossproduct.
    weighted_diff <- temp_diff * sqrt(w_shifted)
    a_est <- crossprod(weighted_diff) / denom_sum

    return(list(
      mu = mu_est,
      a = a_est
    ))
  }
}

log_K_gaussian_mat <- function(mat, x, h) {
  dim <- length(x)
  # Calculate the squared distances scaled by h
  # Using sweep or scale-like logic for better efficiency than matrix(rep...)
  z <- sweep(mat, 2, x, "-") / h

  # Log of the Gaussian product
  # log(dnorm(u)) is -0.5 * u^2 - log(sqrt(2*pi))
  log_probs <- -0.5 * z^2 - log(sqrt(2 * pi))

  # Sum across dimensions for each row, then adjust for h^dim
  log_values <- Rfast::rowsums(log_probs) - (dim * log(h))
  return(log_values)
}

log_K_exp_mat <- function(mat, x, h) {
  dim <- length(x)

  # Calculate absolute differences scaled by h
  # sweep() is efficient for row-wise or column-wise operations
  z <- abs(sweep(mat, 2, x, "-")) / h

  # The log of the exponential part is just -z.
  # We sum these logs across the dimensions (rowSums)
  # and subtract the normalization constant for d dimensions.
  log_values <- Rfast::rowsums(-z) - (dim * log(2 * h))

  return(log_values)
}
