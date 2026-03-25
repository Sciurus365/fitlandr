# The 2nd-generation functions for calculating landscapes from vector fields

#' @import Matrix
#' @importClassesFrom Matrix dgCMatrix
NULL


#' Estimates the steady-state distribution using Central Differencing for an
#' arbitrary rectangular domain defined by x_range and y_range.
#'
#' @param vf A `vectorfield` or `cv_vectorfield` object representing the vector field.
#' @param linear_interp Logical indicating whether to use linear interpolation in predictions.
#' @param n_grid The number of grid points along one dimension (e.g., 100).
#' @return matrix of probability density.
ss_fp_2d <- function(vf, linear_interp = TRUE, n_grid = 100) {
  # Extract drift and diffusion functions from the vector field object

  if (inherits(vf, "cv_vectorfield")) {
    vf <- vf$final_model
  } else if (!inherits(vf, "vectorfield")) {
    cli::cli_abort("Input {.arg vf} must be a {.cls vectorfield} or {.cls cv_vectorfield} object.")
  }

  drift_func <- function(x, y) {
    pred <- stats::predict(vf, c(x, y), linear_interp = linear_interp)
    return(pred$v)
  }
  diffusion_func <- function(x, y) {
    pred <- stats::predict(vf, c(x, y), linear_interp = linear_interp)
    return(pred$a)
  }

  x_range <- vf$lims[1:2]
  y_range <- vf$lims[3:4]

  # --- 1. Setup Grid & Pre-Calculate Fields ---

  cli::cli_progress_step("Setting up grid and pre-calculating fields...")

  N <- n_grid * n_grid
  hx <- (x_range[2] - x_range[1]) / n_grid
  hy <- (y_range[2] - y_range[1]) / n_grid
  x_coords <- seq(x_range[1] + hx / 2, x_range[2] - hx / 2, length.out = n_grid)
  y_coords <- seq(y_range[1] + hy / 2, y_range[2] - hy / 2, length.out = n_grid)

  # Pre-calculate fields
  Ax <- matrix(0, n_grid, n_grid)
  Ay <- matrix(0, n_grid, n_grid)
  Dxx <- matrix(0, n_grid, n_grid)
  Dyy <- matrix(0, n_grid, n_grid)
  Dxy <- matrix(0, n_grid, n_grid)

  for (i in 1:n_grid) {
    for (j in 1:n_grid) {
      drt <- drift_func(x_coords[i], y_coords[j])
      diff <- diffusion_func(x_coords[i], y_coords[j])
      Ax[i, j] <- drt[1]
      Ay[i, j] <- drt[2]
      Dxx[i, j] <- diff[1, 1]
      Dyy[i, j] <- diff[2, 2]
      Dxy[i, j] <- diff[1, 2]
    }
  }

  cli::cli_progress_step("Building sparse matrix representation...")

  # Triplet accumulation
  i_vec <- integer(N * 25)
  j_vec <- integer(N * 25)
  x_vec <- numeric(N * 25)
  cnt <- 0

  add_entry <- function(row, col, val) {
    if (val == 0) {
      return()
    }
    cnt <<- cnt + 1
    i_vec[cnt] <<- row
    j_vec[cnt] <<- col
    x_vec[cnt] <<- val
  }

  k_idx <- function(i, j) i + (j - 1) * n_grid

  # 1. X-Interfaces (Flow Jx between (i,j) and (i+1, j))
  for (j in 1:n_grid) {
    for (i in 1:(n_grid - 1)) {
      k1 <- k_idx(i, j)
      k2 <- k_idx(i + 1, j)

      # A. Standard Drift & Dxx contribution
      A_mid <- (Ax[i, j] + Ax[i + 1, j]) / 2
      coeff_k1 <- (A_mid / 2) + (Dxx[i, j] / hx)
      coeff_k2 <- (A_mid / 2) - (Dxx[i + 1, j] / hx)

      # B. Mixed Dxy contribution: -d(Dxy * rho)/dy
      # We use a 4-point central difference for the y-gradient at the interface
      jp <- if (j == n_grid) j else j + 1
      jm <- if (j == 1) j else j - 1
      denom_y <- if (j == 1 || j == n_grid) hy else 2 * hy

      # Contribution from neighbors to the flux Jx
      # Jx_mixed = - [ (Dxy*rho)_{i+1/2, j+1} - (Dxy*rho)_{i+1/2, j-1} ] / 2hy
      coeff_mixed <- 1 / (2 * denom_y)

      # Add to matrix (Rate of change = -div(J))
      # Normal Flow
      add_entry(k1, k1, -coeff_k1 / hx)
      add_entry(k1, k2, -coeff_k2 / hx)
      add_entry(k2, k1, coeff_k1 / hx)
      add_entry(k2, k2, coeff_k2 / hx)

      # Mixed Flow (affects k1 and k2 by drawing from surrounding y-cells)
      for (curr_i in c(i, i + 1)) {
        mult <- if (curr_i == i) -1 else 1 # Sign change for k1 vs k2
        add_entry(if (mult < 0) k1 else k2, k_idx(curr_i, jp), -mult * Dxy[curr_i, jp] * coeff_mixed / hx)
        add_entry(if (mult < 0) k1 else k2, k_idx(curr_i, jm), mult * Dxy[curr_i, jm] * coeff_mixed / hx)
      }
    }
  }

  # 2. Y-Interfaces (Flow Jy between (i,j) and (i, j+1))
  for (i in 1:n_grid) {
    for (j in 1:(n_grid - 1)) {
      k1 <- k_idx(i, j)
      k2 <- k_idx(i, j + 1)

      A_mid <- (Ay[i, j] + Ay[i, j + 1]) / 2
      coeff_k1 <- (A_mid / 2) + (Dyy[i, j] / hy)
      coeff_k2 <- (A_mid / 2) - (Dyy[i, j + 1] / hy)

      # Mixed Dxy contribution: -d(Dxy * rho)/dx
      ip <- if (i == n_grid) i else i + 1
      im <- if (i == 1) i else i - 1
      denom_x <- if (i == 1 || i == n_grid) hx else 2 * hx
      coeff_mixed <- 1 / (2 * denom_x)

      add_entry(k1, k1, -coeff_k1 / hy)
      add_entry(k1, k2, -coeff_k2 / hy)
      add_entry(k2, k1, coeff_k1 / hy)
      add_entry(k2, k2, coeff_k2 / hy)

      for (curr_j in c(j, j + 1)) {
        mult <- if (curr_j == j) -1 else 1
        add_entry(if (mult < 0) k1 else k2, k_idx(ip, curr_j), -mult * Dxy[ip, curr_j] * coeff_mixed / hy)
        add_entry(if (mult < 0) k1 else k2, k_idx(im, curr_j), mult * Dxy[im, curr_j] * coeff_mixed / hy)
      }
    }
  }


  M <- Matrix::sparseMatrix(i = i_vec[1:cnt], j = j_vec[1:cnt], x = x_vec[1:cnt], dims = c(N, N))

  # 3. Final Solve
  cli::cli_progress_step("Solving for steady-state distribution...")
  # Constraint: Integral of rho = 1
  ones <- rep(1, N) * (hx * hy)
  M_aug <- rbind(cbind(M, rep(1, N)), c(ones, 0))
  b_aug <- c(rep(0, N), 1)

  sol <- solve(M_aug, b_aug)
  rho_ss <- matrix(sol[1:N], n_grid, n_grid)
  attr(rho_ss, "M") <- M
  attr(rho_ss, "x_coords") <- x_coords
  attr(rho_ss, "y_coords") <- y_coords
  cli::cli_progress_done()
  return(rho_ss)
}

#' Creates a 2D Potential Landscape object from a vector field.
#' @param vf A `vectorfield` or `cv_vectorfield` object representing the vector field.
#' @param linear_interp Logical indicating whether to use linear interpolation in predictions.
#' @param n_grid The number of grid points along one dimension (e.g., 100).
#' @return An object of class `2d_static_ld` containing:
#'        - `dist`: A data frame with columns x, y, d (steady-state distribution), U (potential).
#'        - `plot`: A plotly surface plot of the potential landscape U.
#'        - `plot_2`: A ggplot2 raster plot of the potential landscape U.
#'        - `vf`: The input vector field object.
#'        - `ss`: The steady-state distribution matrix.
#'
#' @export
make_2d_ld <- function(vf, linear_interp = TRUE, n_grid = 100) {
  ss <- ss_fp_2d(vf, linear_interp = linear_interp)
  if (min(ss) <= 0) {
    cli::cli_warn("Steady-state distribution contains non-positive values,
                  cannot compute potential landscape directly.
                  The smallest value is {min(ss)}.
                  I will try adding a small constant to ss to fix this.")
    ss <- ss + abs(min(ss)) + 1e-10
  }
  U <- -log(ss)

  # Make a regular data frame for plotting. It contains x, y, d, U. (d: steady-state distribution, ss)

  x_coords <- attr(ss, "x_coords")
  y_coords <- attr(ss, "y_coords")
  dist <- expand.grid(x = x_coords, y = y_coords)
  dist$d <- as.vector(ss)
  dist$U <- as.vector(U)

  plot <- plotly::plot_ly(
    data = dist,
    x = x_coords, y = y_coords, z = U,
    type = "surface"
  ) %>%
    plotly::layout(scene = list(
      xaxis = list(title = vf$x),
      yaxis = list(title = vf$y), zaxis = list(title = "U")
    )) %>%
    plotly::colorbar(title = "U")

  plot_2 <- ggplot2::ggplot(
    dist,
    ggplot2::aes(x = x, y = y)
  ) +
    ggplot2::geom_raster(ggplot2::aes(fill = U)) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(
      x = vf$x,
      y = vf$y, fill = "U"
    ) +
    ggplot2::theme_bw()

  return(structure(list(
    dist = dist,
    plot = plot,
    plot_2 = plot_2,
    vf = vf,
    ss = ss
  ), class = c("2d_static_ld", "2d_ld", "landscape")))
}


#' Finds local minima in a 2D landscape object.
#'
#' @param ld A `2d_ld` or `2d_static_ld` object representing the landscape.
#' @param min_barrier When the barrier height between a local minimum and any of the
#' saddle points connecting it to other local minima is less than `min_barrier` times
#' the highest barrier height value, the local minimum will be considered minor.
#' Default is 0.1.
#' @param exclude_minor Logical indicating whether to mark minor local minima
#' based on the barrier height criterion, so that they can be easily excluded
#' from subsequent calculations. Default is TRUE.
#' @return A data frame with columns x, y, U for each local minimum found.
#' @export
find_loc_min <- function(ld, exclude_minor = TRUE, min_barrier = 0.05) {
  if (!inherits(ld, "2d_ld") && !inherits(ld, "2d_static_ld")) {
    cli::cli_abort("Input {.arg ld} must be a {.cls 2d_ld} or {.cls 2d_static_ld} object.")
  }

  dist <- ld$dist
  U_matrix <- matrix(dist$U, nrow = length(unique(dist$x)), ncol = length(unique(dist$y)))

  n_x <- nrow(U_matrix)
  n_y <- ncol(U_matrix)

  local_mins <- data.frame(x = numeric(0), y = numeric(0), U = numeric(0))

  for (i in 2:(n_x - 1)) {
    for (j in 2:(n_y - 1)) {
      current_U <- U_matrix[i, j]
      neighbors <- c(
        U_matrix[i - 1, j], U_matrix[i + 1, j],
        U_matrix[i, j - 1], U_matrix[i, j + 1],
        U_matrix[i - 1, j - 1], U_matrix[i - 1, j + 1],
        U_matrix[i + 1, j - 1], U_matrix[i + 1, j + 1]
      )

      if (all(current_U < neighbors)) {
        local_mins <- rbind(local_mins, data.frame(
          x = unique(dist$x)[i],
          y = unique(dist$y)[j],
          U = current_U
        ))
      }
    }
  }

  # calculate the barrier height for each pair of local minima
  # this is not needed if only one local minimum is found
  # store the results in a matrix. Each row is a minimum, each column is also a minimum
  # and each cell represents the barrier height between the two minima
  # for each pair, there are two barrier heights (from min1 to the saddle point, and from min2 to the saddle point)
  # therefore, from i to the barrier and from j to the barrier will be stored at (i, j) and (j, i), respectively
  # if i == j, the barrier height is NA

  # but the problem with this approach is that if there are two "twin" local minima close to each other
  # they are "disqualified" together
  # stepwise? then it's ... complicated
  # Oh no. it's not. no recalculation needed.

  n_mins <- nrow(local_mins)
  if (n_mins <= 1 || exclude_minor == FALSE) {
    all_barriers <- matrix(NA, nrow = n_mins, ncol = n_mins)
    minor_mins <- c()
  } else {
    all_barriers <- matrix(NA, nrow = n_mins, ncol = n_mins)
    ld_reformated <- ld
    ld_reformated$dist <- transform_to_grid(ld$dist)
    for (i in 1:(n_mins - 1)) {
      for (j in (i + 1):n_mins) {
        barrier_info <- calculate_barrier.2d_ld(
          l = ld_reformated,
          start_location_value = c(local_mins$x[i], local_mins$y[i]),
          start_r = 1e-5,
          end_location_value = c(local_mins$x[j], local_mins$y[j]),
          end_r = 1e-5,
          expand = FALSE,
          Umax = Inf
        ) # something wrong here. why do we have 0 barrier heights?
        all_barriers[i, j] <- summary(barrier_info)[1] %>% as.numeric()
        all_barriers[j, i] <- summary(barrier_info)[2] %>% as.numeric()
      }
    }

    # find the highest barrier height
    max_barrier <- max(all_barriers, na.rm = TRUE)


    # do the following until no remaining barriers are lower than min_barrier * max_barrier
    # first, find the lowest barrier
    # label this minimum as minor
    # remove this minimum from the matrix (set the corresponding row and column to NA)
    # repeat

    minor_mins <- c()
    all_barriers_copy <- all_barriers
    repeat {
      current_min_barrier <- min(all_barriers_copy, na.rm = TRUE)
      if (is.infinite(current_min_barrier) || current_min_barrier >= min_barrier * max_barrier) {
        break
      }
      locs <- which(all_barriers_copy == current_min_barrier, arr.ind = TRUE)
      min_to_remove <- locs[1, 1] # arbitrarily choose the first
      minor_mins <- c(minor_mins, min_to_remove)
      all_barriers_copy[min_to_remove, ] <- NA
      all_barriers_copy[, min_to_remove] <- NA
    }
  }

  local_mins <- local_mins %>%
    dplyr::mutate(is_minor = ifelse(dplyr::row_number() %in% minor_mins, TRUE, FALSE))

  return(structure(list(mins = local_mins, barriers = all_barriers), class = "ld_min"))
}


transform_to_grid <- function(dist) {
  x_vals <- sort(unique(dist$x))
  y_vals <- sort(unique(dist$y))

  # Determine if data is ordered by y first or x first
  # Here we assume each y row has all x's in order
  d_mat <- matrix(dist$d, nrow = length(y_vals), ncol = length(x_vals), byrow = FALSE)

  list(
    x = x_vals,
    y = y_vals,
    d = d_mat
  )
}


get_min_pos <- function(ld_min, index) {
  if (!inherits(ld_min, "ld_min")) {
    cli::cli_abort("Input {.arg ld_min} must be an {.cls ld_min} object.")
  }

  if (index < 1 || index > nrow(ld_min$mins)) {
    cli::cli_abort("{.arg index} is out of bounds.")
  }

  return(c(ld_min$mins$x[index], ld_min$mins$y[index]))
}
