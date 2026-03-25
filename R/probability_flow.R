#' Generates Probability Flow vectors and a corresponding plot
#' for a given vector field and landscape.
#'
#' @param vf A `vectorfield` or `cv_vectorfield` object representing the vector field.
#' @param ld A `2d_ld` or `2d_static_ld` object representing the landscape.
#' @param n Number of flow vectors to generate along each dimension (default 20).
#' @param divided_by_rho Logical indicating whether to divide flow vectors by steady-state density.
#' @return An object of class `2d_pf` containing:
#'         - `pf_data`: Data frame with columns x, y, Jx, Jy.
#'         - `plot`: ggplot2 object visualizing the probability flow vectors.
#'         - `vf`: The input vector field object.
#'         - `ld`: The input landscape object.
make_2d_pf <- function(vf, ld, n = 20, divided_by_rho = FALSE) {
  if (inherits(vf, "cv_vectorfield")) {
    vf <- vf$final_model
  } else if (!inherits(vf, "vectorfield")) {
    cli::cli_abort("Input {.arg vf} must be a {.cls vectorfield} or {.cls cv_vectorfield} object.")
  }

  rho <- ld$ss
  drift_func <- function(x, y) {
    pred <- stats::predict(vf, c(x, y))
    return(pred$v)
  }
  diffusion_func <- function(x, y) {
    pred <- stats::predict(vf, c(x, y))
    return(pred$a)
  }
  pf_data <- calculate_probability_flow(
    rho = rho,
    drift_func = drift_func,
    diffusion_func = diffusion_func,
    x_range = vf$lims[1:2],
    y_range = vf$lims[3:4],
    n_flow = n,
    devided_by_rho = divided_by_rho
  )
  plot <- ggplot2::ggplot(
    pf_data,
    ggplot2::aes(x = x, y = y)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        xend = x + Jx,
        yend = y + Jy
      ),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm")),
      alpha = 0.7
    ) +
    ggplot2::labs(
      x = vf$x,
      y = vf$y,
      title = "Probability Flow Vectors"
    ) +
    ggplot2::theme_bw()

  return(structure(list(
    vec_grid = pf_data
  ), class = c("2d_pf", "probabilityflow", "vectorfield")))
}


#' Calculates Probability Flow (J) vectors at a subset of grid points
#'
#' @param rho The steady-state density matrix (n_grid x n_grid).
#' @param drift_func Function returning c(Ax, Ay).
#' @param diffusion_func Function returning 2x2 matrix.
#' @param x_range Vector c(xmin, xmax).
#' @param y_range Vector c(ymin, ymax).
#' @param n_flow Number of points along one dimension to sample (default 20).
#' @param devided_by_rho Logical indicating whether to divide flow vectors by steady-state density.
#' This may be useful when performing the force decomposition.
#' @return Data frame with columns: x, y, Jx, Jy.
calculate_probability_flow <- function(rho, drift_func, diffusion_func,
                                       x_range = c(0, 1), y_range = c(0, 1),
                                       n_flow = 20, devided_by_rho = FALSE) {
  # 1. Recover Grid Parameters from Input Rho
  n_grid <- nrow(rho) # Assuming square grid based on previous code
  Lx <- x_range[2] - x_range[1]
  Ly <- y_range[2] - y_range[1]
  hx <- Lx / n_grid
  hy <- Ly / n_grid

  # Full grid coordinates (centers)
  x_coords <- seq(x_range[1] + hx / 2, x_range[2] - hx / 2, length.out = n_grid)
  y_coords <- seq(y_range[1] + hy / 2, y_range[2] - hy / 2, length.out = n_grid)

  # 2. Select Subset Indices (Evenly spaced)
  # We use round() to pick the closest integer indices in the original grid
  idx_sub_x <- unique(round(seq(1, n_grid, length.out = n_flow)))
  idx_sub_y <- unique(round(seq(1, n_grid, length.out = n_flow)))

  # Initialize output vectors
  n_total <- length(idx_sub_x) * length(idx_sub_y)
  out_x <- numeric(n_total)
  out_y <- numeric(n_total)
  out_Jx <- numeric(n_total)
  out_Jy <- numeric(n_total)
  out_Jx_rho <- numeric(n_total)
  out_Jy_rho <- numeric(n_total)

  cnt <- 0

  # 3. Helper: Get D * rho for a specific grid index (i, j)
  # We calculate D locally to avoid storing huge matrices for the full grid
  get_D_rho <- function(i, j) {
    # Clamp indices to boundaries (consistent with your solver's logic)
    i_c <- if (i < 1) 1 else if (i > n_grid) n_grid else i
    j_c <- if (j < 1) 1 else if (j > n_grid) n_grid else j

    val_rho <- rho[i_c, j_c]
    val_diff <- diffusion_func(x_coords[i_c], y_coords[j_c])

    # Return list of D_ab * rho
    list(
      xx = val_diff[1, 1] * val_rho,
      yy = val_diff[2, 2] * val_rho,
      xy = val_diff[1, 2] * val_rho
    )
  }

  # 4. Loop over subset to calculate Flux
  for (i in idx_sub_x) {
    for (j in idx_sub_y) {
      cnt <- cnt + 1

      # Coordinate of current point
      x <- x_coords[i]
      y <- y_coords[j]
      r <- rho[i, j]

      # --- A. Advection Part (Drift * rho) ---
      drift <- drift_func(x, y)
      Adv_x <- drift[1] * r
      Adv_y <- drift[2] * r

      # --- B. Diffusion Gradient Part ---
      # We need derivatives of (D * rho).
      # Using indices i+1/i-1 and j+1/j-1 creates the Central Difference.
      # Note: We clamp indices inside get_D_rho just like the solver did.

      # Neighbors for X-derivative
      val_p_x <- get_D_rho(i + 1, j)
      val_m_x <- get_D_rho(i - 1, j)

      # Neighbors for Y-derivative
      val_p_y <- get_D_rho(i, j + 1)
      val_m_y <- get_D_rho(i, j - 1)

      # Gradients (Central Difference: (f(x+h) - f(x-h)) / 2h)
      # d(Dxx * rho) / dx
      d_Dxx_dx <- (val_p_x$xx - val_m_x$xx) / (2 * hx)
      # d(Dxy * rho) / dx
      d_Dxy_dx <- (val_p_x$xy - val_m_x$xy) / (2 * hx)

      # d(Dyy * rho) / dy
      d_Dyy_dy <- (val_p_y$yy - val_m_y$yy) / (2 * hy)
      # d(Dxy * rho) / dy
      d_Dxy_dy <- (val_p_y$xy - val_m_y$xy) / (2 * hy)

      # --- C. Combine ---
      # Jx = Ax*rho - d(Dxx*rho)/dx - d(Dxy*rho)/dy
      Jx <- Adv_x - d_Dxx_dx - d_Dxy_dy

      # Jy = Ay*rho - d(Dyy*rho)/dy - d(Dxy*rho)/dx
      Jy <- Adv_y - d_Dyy_dy - d_Dxy_dx

      # Store
      out_x[cnt] <- x
      out_y[cnt] <- y
      out_Jx[cnt] <- Jx
      out_Jy[cnt] <- Jy
      out_Jx_rho[cnt] <- Jx / r
      out_Jy_rho[cnt] <- Jy / r
    }
  }
  if (devided_by_rho) {
    return(data.frame(x = out_x, y = out_y, vx = out_Jx_rho, vy = out_Jy_rho))
  } else {
    return(data.frame(x = out_x, y = out_y, vx = out_Jx, vy = out_Jy))
  }
}
