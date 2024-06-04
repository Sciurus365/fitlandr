#' Multivariate vector field kernel estimator
#'
#' See references for details.
#'
#' @param d The dataset. Should be a matrix or a data frame, with each row representing a random vector.
#' @param h The bandwidth for the kernel estimator.
#' @param kernel The type of kernel estimator used. "exp" by default ([exp()]), and if "Gaussian" then [stats::dnorm()] will be used.
#'
#' @return A function(x), which then returns the \eqn{\mu} and \eqn{a} estimators at the position \eqn{x}.
#' @references Bandi, F. M., & Moloche, G. (2018). On the functional estimation of multivariate diffusion processes. Econometric Theory, 34(4), 896-946. https://doi.org/10.1017/S0266466617000305
#' @export
MVKE <- function(d, h = 0.2, kernel = c("exp", "Gaussian")) {
  if (is.data.frame(d)) d <- as.matrix(d)
  if (!is.matrix(d)) stop("`d` should be a data.frame or a matrix.")

  d <- stats::na.omit(d)
  dim <- ncol(d)

  temp_d <- d[1:(nrow(d) - 1), , drop = FALSE]
  temp_diff <- diff(d)
  temp_norm <- apply(temp_diff, MARGIN = 1, FUN = function(x) norm(x, "2"))
  temp_diff_tcrossprod <- apply(temp_diff,
    MARGIN = 1,
    FUN = function(x) {
      tcrossprod(x, x)
    }, simplify = FALSE
  )
  kernel <- kernel[1]
  if (kernel == "Gaussian") {
    K <- K_gaussian_mat
  } else if (kernel == "exp") {
    K <- K_exp_mat
  } else {
    stop('`kernel` must be one of "Gaussian" or "exp".')
  }

  force(h)
  function(x) {
    if (length(x) != dim) stop("Input of wrong dimension.")
    temp_kernel_term_upper <- K_gaussian_mat(temp_d, x, h = h)
    temp_kernel_term_lower <- K_gaussian_mat(d, x, h = h)
    return(list(
      mu = colSums(temp_kernel_term_upper * temp_diff) / sum(temp_kernel_term_lower),
      a = mapply(`*`, temp_kernel_term_upper, temp_diff_tcrossprod, SIMPLIFY = FALSE) %>% Reduce(`+`, .) / sum(temp_kernel_term_lower)
    ))
  }
}

# K_gaussian <- function(x, h) {
#   dim <- length(x)
#   1 / (h^dim) * prod(stats::dnorm(x / h))
# }

K_gaussian_mat <- function(mat, x, h) {
  dim <- length(x)
  mat <- mat - matrix(rep(x, nrow(mat)), ncol = dim, byrow = TRUE)
  mat <- stats::dnorm(mat / h)
  values <- 1 / (h^dim) * Rfast::rowprods(mat)
  return(values)
}

K_exp_mat <- function(mat, x, h) {
  dim <- length(x)
  mat <- mat - matrix(rep(x, nrow(mat)), ncol = dim, byrow = TRUE)
  mat <- exp(mat / h)
  values <- 1 / (h^dim) * Rfast::rowprods(mat)
  return(values)
}
