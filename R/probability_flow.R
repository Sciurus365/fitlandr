#' Generates probability-flow vectors
#' for a given vector field and landscape.
#'
#' @param vf A `vectorfield` or `cv_vectorfield` object representing the vector field.
#' @param ld A `2d_ld` or `2d_static_ld` object representing the landscape.
#' @param n Number of flow vectors to generate along each dimension (default 20).
#' @param divided_by_rho Logical indicating whether to divide flow vectors by steady-state density.
#' @return An object of class `2d_pf` containing:
#'         - `vec_grid`: Data frame with columns x, y, vx, vy.
#'         - `vf`: The input vector field object.
#'         - `ld`: The input landscape object.
#' @export
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
  return(structure(list(
    vec_grid = pf_data,
    x = vf$x,
    y = vf$y,
    vf = vf,
    ld = ld,
    divided_by_rho = divided_by_rho
  ), class = c("2d_pf", "probabilityflow", "vectorfield")))
}


#' Estimate a stream function from a two-dimensional probability flow
#'
#' Estimates a scalar stream function `A` whose perpendicular gradient
#' `(-dA/dy, dA/dx)` is closest in least-squares distance to the estimated
#' probability flow. Derivatives are represented by sparse finite-difference
#' matrices, using centered differences in the grid interior and one-sided
#' differences at its boundary.
#'
#' A stream function is identifiable only up to an additive constant. The
#' returned solution is anchored by setting its value at the first grid point
#' to zero. If the estimated probability flow is not exactly divergence-free,
#' the fitted flow is its least-squares stream-function component and the
#' remaining component is available in `residual_grid`.
#'
#' @param pf A `2d_pf` object returned by [make_2d_pf()].
#'
#' @return An object of class `2d_stream` containing:
#'   - `grid`: grid coordinates and estimated stream-function values.
#'   - `fitted_grid`: observed and fitted probability-flow vectors.
#'   - `residual_grid`: residual probability-flow vectors.
#'   - `rmse`: root mean squared residual across both flow components.
#'   - `relative_error`: residual norm divided by the observed-flow norm.
#'   - `pf`: the input probability-flow object.
#'
#' @export
make_2d_stream <- function(pf) {
  if (!inherits(pf, "2d_pf")) {
    cli::cli_abort("Input {.arg pf} must be a {.cls 2d_pf} object.")
  }
  if (isTRUE(pf$divided_by_rho)) {
    cli::cli_abort(c(
      "Cannot estimate a probability-current stream function from flow divided by density.",
      "i" = "Recreate {.arg pf} with {.code divided_by_rho = FALSE}."
    ))
  }

  flow <- pf$vec_grid
  required_columns <- c("x", "y", "vx", "vy")
  if (!all(required_columns %in% names(flow))) {
    cli::cli_abort("{.arg pf$vec_grid} must contain columns {.field x}, {.field y}, {.field vx}, and {.field vy}.")
  }
  if (any(!is.finite(as.matrix(flow[required_columns])))) {
    cli::cli_abort("{.arg pf$vec_grid} must contain only finite coordinates and flow values.")
  }

  x_coords <- sort(unique(flow$x))
  y_coords <- sort(unique(flow$y))
  nx <- length(x_coords)
  ny <- length(y_coords)
  if (nx < 2L || ny < 2L) {
    cli::cli_abort("The probability-flow grid must contain at least two points along each axis.")
  }
  if (nrow(flow) != nx * ny || anyDuplicated(flow[c("x", "y")])) {
    cli::cli_abort("{.arg pf$vec_grid} must be a complete rectangular grid with one flow vector per point.")
  }

  derivative_matrix <- function(coords) {
    n <- length(coords)
    if (n == 2L) {
      h <- coords[2L] - coords[1L]
      return(Matrix::sparseMatrix(
        i = c(1L, 1L, 2L, 2L),
        j = c(1L, 2L, 1L, 2L),
        x = c(-1, 1, -1, 1) / h,
        dims = c(n, n)
      ))
    }

    interior <- 2L:(n - 1L)
    h_left <- coords[interior] - coords[interior - 1L]
    h_right <- coords[interior + 1L] - coords[interior]
    lower <- -h_right / (h_left * (h_left + h_right))
    center <- (h_right - h_left) / (h_left * h_right)
    upper <- h_left / (h_right * (h_left + h_right))

    Matrix::sparseMatrix(
      i = c(1L, 1L, rep(interior, each = 3L), n, n),
      j = c(
        1L,
        2L,
        as.vector(rbind(interior - 1L, interior, interior + 1L)),
        n - 1L,
        n
      ),
      x = c(
        -1 / (coords[2L] - coords[1L]),
        1 / (coords[2L] - coords[1L]),
        as.vector(rbind(lower, center, upper)),
        -1 / (coords[n] - coords[n - 1L]),
        1 / (coords[n] - coords[n - 1L])
      ),
      dims = c(n, n)
    )
  }

  dx_1d <- derivative_matrix(x_coords)
  dy_1d <- derivative_matrix(y_coords)
  dx_2d <- Matrix::kronecker(Matrix::Diagonal(ny), dx_1d)
  dy_2d <- Matrix::kronecker(dy_1d, Matrix::Diagonal(nx))
  perpendicular_gradient <- rbind(-dy_2d, dx_2d)

  grid_index <- match(flow$x, x_coords) + (match(flow$y, y_coords) - 1L) * nx
  u <- numeric(nx * ny)
  v <- numeric(nx * ny)
  u[grid_index] <- flow$vx
  v[grid_index] <- flow$vy
  observed <- c(u, v)

  anchor <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 1, dims = c(1L, nx * ny))
  system_matrix <- rbind(perpendicular_gradient, anchor)
  stream_values <- as.numeric(Matrix::solve(Matrix::qr(system_matrix), c(observed, 0)))
  fitted <- as.numeric(perpendicular_gradient %*% stream_values)
  residual <- observed - fitted

  ordered_grid <- expand.grid(x = x_coords, y = y_coords)
  fitted_grid <- transform(
    ordered_grid,
    vx = u,
    vy = v,
    fitted_vx = fitted[seq_len(nx * ny)],
    fitted_vy = fitted[nx * ny + seq_len(nx * ny)]
  )
  residual_grid <- transform(
    ordered_grid,
    vx = residual[seq_len(nx * ny)],
    vy = residual[nx * ny + seq_len(nx * ny)]
  )
  observed_norm <- sqrt(sum(observed^2))

  structure(
    list(
      grid = transform(ordered_grid, A = stream_values),
      fitted_grid = fitted_grid,
      residual_grid = residual_grid,
      rmse = sqrt(mean(residual^2)),
      relative_error = if (observed_norm == 0) 0 else sqrt(sum(residual^2)) / observed_norm,
      pf = pf
    ),
    class = "2d_stream"
  )
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
