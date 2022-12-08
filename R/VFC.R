#' Sparse Vector Field Consensus
#'
#' The main function for the SparseVFC algorithm.
#' See `References` for more information. This is a modified version of the original implementation ([SparseVFC::SparseVFC()])
#'
#' @param X The position of the vectors.
#' @param Y The value of the vectors.
#' @param ctrl_pts The control points for constructing kernels. Either a character string "random" (then the control points will be randomly picked from the dataset), or a number (then a grid of control points will be made, with this number of points in each axis.)
#' @param M The number of the basis functions used for sparse approximation. Default value is 16. Only effective when setting `ctrl_points = "random"`.
#' @param MaxIter Maximum iteration times. Default value is 500.
#' @param gamma Percentage of inliers in the samples. This is an initial value for EM iteration, and it is not important. Default value is 0.9.
#' @param beta Parameter of Gaussian Kernel, \eqn{k(x, y) = exp(-beta*||x-y||^2)}. Default value is 0.1.
#' @param lambda Represents the trade-off between the goodness of data fit and smoothness of the field. Default value is 3.
#' @param theta If the posterior probability of a sample being an inlier is larger than theta, then it will be regarded as an inlier. Default value is 0.75.
#' @param a Parameter of the uniform distribution. We assume that the outliers obey a uniform distribution \eqn{1/a}. Default Value is 10.
#' @param ecr The minimum limitation of the energy change rate in the iteration process. Default value is 1e-5.
#' @param minP The posterior probability Matrix P may be singular for matrix inversion. We set the minimum value of P as `minP`. Default value is 1e-5.
#' @param silent Should the messages be suppressed? Default value is `TRUE`.
#'
#' @references
#' The algorithm is described in Ma et al. (2013) \doi{10.1016/j.patcog.2013.05.017}.
#' This function is translated with permission from Jiayi Ma's Matlab function at \url{https://github.com/jiayi-ma/VFC}.
#' Also see Zhao et al. (2011) \doi{10.1109/CVPR.2011.5995336} for the earlier VFC algorithm.
#'
#' @return A `VFC` object, which is a list containing the following elements:
#' \describe{
#' \item{X}{A matrix of the positions of kernels.}
#' \item{Y}{A matrix of the input vectors.}
#' \item{beta}{The input value of `beta`.}
#' \item{V}{A matrix of the estimated vectors.}
#' \item{C}{A matrix of the coefficients of each kernel.}
#' \item{P}{A vector of the posterior probability of the input vectors (`Y`) being an inlier.}
#' \item{VFCIndex}{A vector of indices of the inliers.}
#' \item{sigma2}{The \eqn{\sigma^2} of the estimations weighted by `P`.}
#' }
#'
#' @export
SparseVFC <- function(X, Y, M = 16, MaxIter = 500,
                      gamma = 0.9, beta = 0.1,
                      lambda = 3, theta = 0.75, a = 10, ecr = 1e-5,
                      minP = 1e-5, ctrl_pts = "random", silent = TRUE) {
  if (!silent) message("Start mismatch removal...\n")
  N <- nrow(Y)
  D <- ncol(Y)

  # Construct kernel matrix K
  if (ctrl_pts == "random") {
    tmp_X <- unique(X)
    idx <- sample.int(nrow(tmp_X))
    idx <- idx[1:min(M, nrow(tmp_X))]
    ctrl_pts <- tmp_X[idx, ]
  } else if (is.numeric(ctrl_pts) & length(ctrl_pts) == 1) {
    max1 <- max(X[, 1])
    min1 <- min(X[, 1])
    max2 <- max(X[, 2])
    min2 <- min(X[, 2])
    ctrl_pts <- as.matrix(expand.grid(
      seq(min1, max1, length.out = ctrl_pts),
      seq(min2, max2, length.out = ctrl_pts)
    ))
  } else {
    stop("`ctrl_pts` must be 'random', or a number.")
  }

  K <- con_K(ctrl_pts, ctrl_pts, beta)
  U <- con_K(X, ctrl_pts, beta)
  M <- nrow(ctrl_pts)

  # Initialization
  V <- matrix(0, nrow = N, ncol = D)
  iter <- 1
  tecr <- 1
  C <- matrix(0, nrow = M, ncol = D)
  E <- 1
  sigma2 <- sum((Y - V)^2) / (N * D)

  while (iter < MaxIter & tecr > ecr) {
    # E-step.
    E_old <- E
    temp_PE <- get_P(Y, V, sigma2, gamma, a)
    P <- temp_PE$P
    E <- temp_PE$E
    tecr <- abs((E - E_old) / E)
    if (!silent) {
      message(
        sprintf("iterate: %dth, gamma: %f, the energy change rate: %f, sigma2=%f\n", iter, gamma, tecr, sigma2)
      )
    }

    # M-step. Solve linear system for C.
    P <- pmax(P, minP)
    C <- solve((t(U) * matrix(rep(P, M), nrow = M, byrow = TRUE)) %*% U + lambda * sigma2 * K) %*%
      ((t(U) * matrix(rep(P, M), nrow = M, byrow = TRUE)) %*% Y)

    # Update V and sigma^2
    V <- U %*% C
    Sp <- sum(P)
    sigma2 <- sum(P * rowSums((Y - V)^2)) / (Sp * D)

    # Update gamma
    numcorr <- sum(P > theta)
    gamma <- numcorr / nrow(X)
    if (gamma > 0.95) gamma <- 0.95
    if (gamma < 0.05) gamma <- 0.05

    iter <- iter + 1
  }

  if (!silent) message("Removing outliers succesfully completed.")

  result <- list(
    X = ctrl_pts,
    Y = Y,
    beta = beta,
    V = V,
    C = C,
    P = P,
    VFCIndex = which(P > theta),
    sigma2 = sigma2
  )

  class(result) <- c("VFC")
  return(
    result
  )
}

#' Construct the kernel K
#'
#' \deqn{K[i, j] = k(x[i,], y[j,]) = exp(-beta*||x[i,]-y[j,]||^2)}
#'
#' @param x,y,beta The variable and parameter values for the funcion
#' described above.
#' @seealso [SparseVFC()]
#' @return The kernel estimation result.
#' @noRd
con_K <- function(x, y, beta) {
  n <- nrow(x)
  m <- nrow(y)
  if (ncol(x) != ncol(y)) stop("ncol(x) != ncol(y)")
  d <- ncol(x)

  if (n == 1) {
    K <- matrix(sum((x - y)^2)^(1 / 2))
  } else {
    K <- as.matrix(purrr::quietly(pdist::pdist)(x, y)$result)
  }
  K <- exp(-beta * (K^2))

  return(K)
}

#' Estimate the posterior probability and part of the energy
#' @noRd
get_P <- function(Y, V, sigma2, gamma, a) {
  D <- ncol(Y)
  temp1 <- exp(-rowSums((Y - V)^2) / (2 * sigma2))
  temp2 <- (2 * pi * sigma2)^(D / 2) * (1 - gamma) / (gamma * a)
  P <- temp1 / (temp1 + temp2)
  E <- sum(P * rowSums((Y - V)^2) / (2 * sigma2)) + sum(P) * log(sigma2) * D / 2
  return(list(P = P, E = E))
}

#' Normalize (a Matrix of) Vectors
#'
#' Normalize the data so that the mean of the vectors is **0** and the variance of the vectors is 1. Here the variance of vectors is calculated by interpreting the deviation as the Euclidean distance, which means the trace of the (population) covariance matrix is 1.
#'
#' @param x The matrix to be normalized. Each row of `x` represent a vector.
#' @return The normalized matrix with an attribution `scale`, which is the scale factor used for normalization.
#'
#' @export
#' @examples
#' norm_vecs(matrix(rep(1, 100), nrow = 2))
norm_vecs <- function(x) {
  n <- nrow(x)
  xm <- colMeans(x)
  x <- x - matrix(rep(xm, n), nrow = n, byrow = TRUE)
  xscale <- sqrt(sum(x^2) / n)

  X <- x / xscale

  return(
    structure(
      X,
      scale = xscale
    )
  )
}

#' Predict Method for VFC Fits.
#'
#' Predicted values based on `VFC` objects.
#'
#' @param object A `VFC` object generated by [SparseVFC()].
#' @param newdata A vector specifying the position.
#' @param ... Not in use.
#' @return A vector.
#' @seealso [SparseVFC()]
#' @export
predict.VFC <- function(object, newdata, ...) {
  if (ncol(object$X) != length(newdata)) {
    stop("The dimension of the `newdata` should be the same as in the model.")
  }
  output <- apply(cbind(object$X, object$C), 1, function(x, beta) {
    con_K(newdata |> matrix(nrow = 1), x[1:(length(x) / 2)] |> matrix(nrow = 1), beta) |>
      as.numeric() * x[(length(x) / 2 + 1):length(x)]
  }, beta = object$beta)

  output <- rowSums(simplify2array(output))
  return(output)
}
