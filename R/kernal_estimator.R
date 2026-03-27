#' Multivariate vector field kernel estimator
#'
#' See references for details.
#'
#' @param d The dataset. Should be a matrix or a data frame, with each row representing a random vector.
#' @param v The vectors corresponding to the dataset. Should be a matrix or a data frame with the same shape as `d`. If missing, then the vectors will be calculated from the dataset.
#' @param h The bandwidth for the kernel estimator.
#' @param kernel The type of kernel estimator used. "Gaussian" by default.
#'
#' @return A function(x), which then returns the \eqn{\mu} and \eqn{a} estimators at the position \eqn{x}.
#' @references Bandi, F. M., & Moloche, G. (2018). On the functional estimation of multivariate diffusion processes. Econometric Theory, 34(4), 896-946. https://doi.org/10.1017/S0266466617000305
#' @export
MVKE <- function(d, v, h = 0.2, kernel = c("Gaussian", "exp")) {
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

  temp_diff <- v
  kernel <- kernel[1]
  if (kernel == "Gaussian") {
    log_K <- log_K_gaussian_mat
  } else if (kernel == "exp") {
    log_K <- log_K_exp_mat
  } else {
    cli::cli_abort('{.arg kernel} must be one of "Gaussian" or "exp".')
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

    # Sum_i w_i * (v_i %*% t(v_i)) using a weighted crossproduct.
    weighted_diff <- temp_diff * sqrt(w_shifted)
    a_est <- crossprod(weighted_diff) / denom_sum

    return(list(
      mu = colSums(w_shifted * temp_diff) / denom_sum,
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
