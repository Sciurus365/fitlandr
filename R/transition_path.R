#' Calculate the Minimum Action Path using gMAM
#'
#' @param vf A `vectorfield` or `cv_vectorfield` object.
#' @param linear_interp Logical indicating whether to use linear interpolation in predictions.
#' @param start A numeric vector c(x, y) for the starting point.
#' @param end A numeric vector c(x, y) for the destination point.
#' @param n_beads Number of points to represent the path (default 50).
#' @param max_iter Number of optimization-reparameterization cycles.
#' @param progress Logical indicating whether to show a progress bar during optimization.
#' @return A data frame containing the x and y coordinates of the optimal path.
find_gmam_path <- function(vf, linear_interp = TRUE, start, end, n_beads = 10, max_iter = 10, progress = interactive()) {
  # 1. Initialization: Create a straight line between start and end
  path_x <- seq(start[1], end[1], length.out = n_beads)
  path_y <- seq(start[2], end[2], length.out = n_beads)

  get_phys <- function(pos) {
    pred <- stats::predict(vf, pos, linear_interp = linear_interp)
    # Inverse of diffusion matrix a
    a_inv <- solve(pred$a + diag(1e-8, 2)) # Small epsilon for stability
    return(list(b = pred$v, a_inv = a_inv))
  }

  # 2. Define the Geometric Action Function
  # par: numeric vector of length 2*(n_beads-2) representing internal points
  calculate_action <- function(internal_points) {
    # Reconstruct full path
    coords <- matrix(0, nrow = n_beads, ncol = 2)
    coords[1, ] <- start
    coords[n_beads, ] <- end
    coords[2:(n_beads - 1), ] <- matrix(internal_points, ncol = 2)

    S <- 0
    # Numerical integration over segments
    for (i in 1:(n_beads - 1)) {
      # Midpoint for evaluating fields
      mid <- (coords[i, ] + coords[i + 1, ]) / 2
      phys <- get_phys(mid)

      # Tangent vector (delta x)
      dx <- coords[i + 1, ] - coords[i, ]

      # Term 1: ||phi'||_a * ||b||_a
      # ||v||_a = sqrt(v^T * a_inv * v)
      norm_dx <- sqrt(as.numeric(t(dx) %*% phys$a_inv %*% dx))
      norm_b <- sqrt(as.numeric(t(phys$b) %*% phys$a_inv %*% phys$b))

      # Term 2: <phi', b>_a = dx^T * a_inv * b
      work_term <- as.numeric(t(dx) %*% phys$a_inv %*% phys$b)

      S <- S + (norm_dx * norm_b - work_term)
    }
    return(S)
  }

  # 3. Iterative Optimization and Reparameterization
  current_internal <- as.vector(cbind(path_x[2:(n_beads - 1)], path_y[2:(n_beads - 1)]))

  if (isTRUE(progress)) {
    cli::cli_progress_bar("Optimizing Path", total = max_iter)
  }


  for (iter in 1:max_iter) {
    # A. Optimize the positions (L-BFGS-B is efficient for high dimensions)
    res <- stats::optim(
      par = current_internal,
      fn = calculate_action,
      method = "L-BFGS-B",
      control = list(maxit = 100)
    )
    current_internal <- res$par

    # B. Re-parameterize (Keep beads equidistant along arc length)
    # This prevents beads from bunching up in low-drift areas
    full_coords <- rbind(start, matrix(current_internal, ncol = 2), end)

    # Calculate cumulative arc length
    diffs <- apply(full_coords, 2, diff)
    dists <- sqrt(rowSums(diffs^2))
    cum_dist <- c(0, cumsum(dists))
    total_len <- max(cum_dist)

    # Interpolate new points at equal intervals
    new_s <- seq(0, total_len, length.out = n_beads)
    new_x <- stats::approx(cum_dist, full_coords[, 1], xout = new_s)$y
    new_y <- stats::approx(cum_dist, full_coords[, 2], xout = new_s)$y

    current_internal <- as.vector(cbind(new_x[2:(n_beads - 1)], new_y[2:(n_beads - 1)]))

    if (isTRUE(progress)) {
      cli::cli_progress_update()
    }
  }

  final_path <- data.frame(
    x = c(start[1], new_x[2:(n_beads - 1)], end[1]),
    y = c(start[2], new_y[2:(n_beads - 1)], end[2])
  )

  return(structure(list(path = final_path), class = "gmam_path"))
}

#' @export
autolayer.gmam_path <- function(object, ...) {
  list(ggplot2::geom_path(
    data = object$path,
    ggplot2::aes(x = x, y = y),
    color = "red",
    ...
  ), ggplot2::geom_point(
    data = object$path[c(1), ],
    ggplot2::aes(x = x, y = y),
    color = "red", size = 3, shape = 16
  ), ggplot2::geom_point(
    data = object$path[nrow(object$path), ],
    ggplot2::aes(x = x, y = y),
    color = "red", size = 3, shape = 17
  ))
}

#' @importFrom ggplot2 autolayer autoplot
NULL
