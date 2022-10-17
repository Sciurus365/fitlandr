#' Multivariate vector field kernel estimator
#'
#' Gaussian product kernel is used.
#'
#' @param d The dataset. Should be a matrix or a data frame, with each row representing a random vector.
#' @param h The bandwidth for the kernel estimator.
#'
#' @return A function(x), which then returns the \eqn{\mu} and \eqn{a} estimators at the position \eqn{x}.
#' @references Bandi, F. M., & Moloche, G. (2018). On the functional estimation of multivariate diffusion processes. Econometric Theory, 34(4), 896-946. https://doi.org/10.1017/S0266466617000305
#' @export
MVKE <- function(d, h = 1) {
  if (is.data.frame(d)) d <- as.matrix(d)
  if (!is.matrix(d)) stop("`d` should be a data.frame or a matrix.")

  d <- na.omit(d)
  dim <- ncol(d)

  temp_d <- d[1:(nrow(d) - 1), ]
  temp_diff <- diff(d)
  temp_norm <- apply(temp_diff, MARGIN = 1, FUN = function(x) norm(x, "2"))

  force(h)
  function(x) {
  	if(length(x) != dim) stop("Input of wrong dimension.")
    temp_kernel_term_upper <- apply(temp_d, MARGIN = 1, FUN = function(y) {
      K_gaussian(y - x, h = h)
    })
    temp_kernel_term_lower <- apply(d, MARGIN = 1, FUN = function(y) {
      K_gaussian(y - x, h = h)
    })
    return(list(
      mu = colSums(temp_kernel_term_upper * temp_diff) / sum(temp_kernel_term_lower),
      a = mapply(`*`, temp_kernel_term_upper, apply(temp_diff,
        MARGIN = 1,
        FUN = function(x) {
          x %*% t(x)
        }, simplify = FALSE
      ), SIMPLIFY = FALSE) %>% Reduce(`+`, .) / sum(temp_kernel_term_lower)
    ))
  }
}

K_gaussian <- function(x, h) {
  dim <- length(x)
  1 / (h^dim) * prod(dnorm(x / h))
}
