# The 2nd-generation functions for calculating landscapes from vector fields

#' @import Matrix
#' @importClassesFrom Matrix dgCMatrix
NULL


#' Estimates the steady-state distribution from a 2D drift-diffusion field.
#'
#' @param vf A `vectorfield` or `cv_vectorfield` object representing the vector field.
#' @param linear_interp Logical indicating whether to use linear interpolation in predictions.
#' @param n_grid The number of grid points along one dimension (e.g., 100).
#' @param drift_scheme Drift discretization scheme. One of `"upwind"` or `"central"`.
#' @param cross_diffusion_mode Cross-diffusion handling mode. One of
#' `"drop"` (default, ignore `Dxy`) or `"full"` (use full diffusion tensor).
#' @param boundary_mode Boundary handling mode. One of `"reflective"` (explicit
#' zero-normal-flux handling) or `"legacy_implicit"` (index-clamping behavior).
#' @return matrix of probability density.
ss_fp_2d <- function(vf,
                     linear_interp = TRUE,
                     n_grid = 100,
                     drift_scheme = c("upwind", "central"),
                     cross_diffusion_mode = c("drop", "full"),
                     boundary_mode = c("reflective", "legacy_implicit")) {
  # Extract drift and diffusion functions from the vector field object

  drift_scheme <- match.arg(drift_scheme)
  cross_diffusion_mode <- match.arg(cross_diffusion_mode)
  boundary_mode <- match.arg(boundary_mode)

  if (inherits(vf, "cv_vectorfield")) {
    vf <- vf$final_model
  } else if (!inherits(vf, "vectorfield")) {
    cli::cli_abort("Input {.arg vf} must be a {.cls vectorfield} or {.cls cv_vectorfield} object.")
  }

  predict_both <- function(x, y) {
    stats::predict(vf, c(x, y), linear_interp = linear_interp)
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
      pred <- predict_both(x_coords[i], y_coords[j])
      drt <- pred$v
      diff <- pred$a
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

  add_entries <- function(rows, cols, vals) {
    keep <- vals != 0
    if (!any(keep)) {
      return()
    }
    rows <- rows[keep]
    cols <- cols[keep]
    vals <- vals[keep]
    n_add <- length(vals)
    idx <- (cnt + 1):(cnt + n_add)
    i_vec[idx] <<- rows
    j_vec[idx] <<- cols
    x_vec[idx] <<- vals
    cnt <<- cnt + n_add
  }

  reflect_index <- function(idx, n) {
    if (idx < 1L) {
      return(2L - idx)
    }
    if (idx > n) {
      return(2L * n - idx)
    }
    idx
  }

  k_idx <- function(i, j) i + (j - 1) * n_grid

  # 1. X-Interfaces (Flow Jx between (i,j) and (i+1, j))
  for (j in 1:n_grid) {
    for (i in 1:(n_grid - 1)) {
      k1 <- k_idx(i, j)
      k2 <- k_idx(i + 1, j)

      # A. Drift + Dxx contribution
      A_mid <- (Ax[i, j] + Ax[i + 1, j]) / 2
      if (drift_scheme == "upwind") {
        adv_k1 <- if (A_mid >= 0) A_mid else 0
        adv_k2 <- if (A_mid < 0) A_mid else 0
        coeff_k1 <- adv_k1 + (Dxx[i, j] / hx)
        coeff_k2 <- adv_k2 - (Dxx[i + 1, j] / hx)
      } else {
        coeff_k1 <- (A_mid / 2) + (Dxx[i, j] / hx)
        coeff_k2 <- (A_mid / 2) - (Dxx[i + 1, j] / hx)
      }

      # Add to matrix (Rate of change = -div(J))
      # Normal Flow
      add_entries(
        rows = c(k1, k1, k2, k2),
        cols = c(k1, k2, k1, k2),
        vals = c(-coeff_k1 / hx, -coeff_k2 / hx, coeff_k1 / hx, coeff_k2 / hx)
      )

      if (cross_diffusion_mode != "drop") {
        # B. Mixed Dxy contribution: -d(Dxy * rho)/dy
        if (boundary_mode == "reflective") {
          jp <- reflect_index(j + 1L, n_grid)
          jm <- reflect_index(j - 1L, n_grid)
          denom_y <- 2 * hy
        } else {
          jp <- if (j == n_grid) j else j + 1
          jm <- if (j == 1) j else j - 1
          denom_y <- if (j == 1 || j == n_grid) hy else 2 * hy
        }

        coeff_mixed <- 1 / denom_y

        # Mixed Flow (affects k1 and k2 by drawing from surrounding y-cells)
        for (curr_i in c(i, i + 1)) {
          mult <- if (curr_i == i) -1 else 1
          row_target <- if (mult < 0) k1 else k2

          dxy_jp <- Dxy[curr_i, jp]
          dxy_jm <- Dxy[curr_i, jm]

          add_entries(
            rows = c(row_target, row_target),
            cols = c(k_idx(curr_i, jp), k_idx(curr_i, jm)),
            vals = c(-mult * dxy_jp * coeff_mixed / hx, mult * dxy_jm * coeff_mixed / hx)
          )
        }
      }
    }
  }

  # 2. Y-Interfaces (Flow Jy between (i,j) and (i, j+1))
  for (i in 1:n_grid) {
    for (j in 1:(n_grid - 1)) {
      k1 <- k_idx(i, j)
      k2 <- k_idx(i, j + 1)

      A_mid <- (Ay[i, j] + Ay[i, j + 1]) / 2
      if (drift_scheme == "upwind") {
        adv_k1 <- if (A_mid >= 0) A_mid else 0
        adv_k2 <- if (A_mid < 0) A_mid else 0
        coeff_k1 <- adv_k1 + (Dyy[i, j] / hy)
        coeff_k2 <- adv_k2 - (Dyy[i, j + 1] / hy)
      } else {
        coeff_k1 <- (A_mid / 2) + (Dyy[i, j] / hy)
        coeff_k2 <- (A_mid / 2) - (Dyy[i, j + 1] / hy)
      }

      add_entries(
        rows = c(k1, k1, k2, k2),
        cols = c(k1, k2, k1, k2),
        vals = c(-coeff_k1 / hy, -coeff_k2 / hy, coeff_k1 / hy, coeff_k2 / hy)
      )

      if (cross_diffusion_mode != "drop") {
        # Mixed Dxy contribution: -d(Dxy * rho)/dx
        if (boundary_mode == "reflective") {
          ip <- reflect_index(i + 1L, n_grid)
          im <- reflect_index(i - 1L, n_grid)
          denom_x <- 2 * hx
        } else {
          ip <- if (i == n_grid) i else i + 1
          im <- if (i == 1) i else i - 1
          denom_x <- if (i == 1 || i == n_grid) hx else 2 * hx
        }
        coeff_mixed <- 1 / denom_x

        for (curr_j in c(j, j + 1)) {
          mult <- if (curr_j == j) -1 else 1
          row_target <- if (mult < 0) k1 else k2

          dxy_ip <- Dxy[ip, curr_j]
          dxy_im <- Dxy[im, curr_j]

          add_entries(
            rows = c(row_target, row_target),
            cols = c(k_idx(ip, curr_j), k_idx(im, curr_j)),
            vals = c(-mult * dxy_ip * coeff_mixed / hy, mult * dxy_im * coeff_mixed / hy)
          )
        }
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

  if (cross_diffusion_mode == "full") {
    n_negative <- sum(rho_ss < 0, na.rm = TRUE)
    neg_warn_threshold <- max(5L, ceiling(0.01 * N))
    if (n_negative >= neg_warn_threshold) {
      cli::cli_warn(paste0(
        "Using cross_diffusion_mode='full' produced ", n_negative,
        " negative density cell(s). Consider cross_diffusion_mode='drop' (default) if this causes instability."
      ))
    }
  }

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
#' @param drift_scheme Drift discretization scheme. One of `"upwind"` or `"central"`.
#' @param cross_diffusion_mode Cross-diffusion handling mode. One of
#' `"drop"` (default, ignore `Dxy`) or `"full"` (use full diffusion tensor).
#' @param boundary_mode Boundary handling mode. One of `"reflective"` (explicit
#' zero-normal-flux handling) or `"legacy_implicit"` (index-clamping behavior).
#' @return An object of class `2d_static_ld` containing:
#'        - `dist`: A data frame with columns x, y, d (steady-state distribution), U (potential).
#'        - `plot`: A plotly surface plot of the potential landscape U, retained
#'          for compatibility. Use [plotly_ld()] to access it.
#'        - `plot_2`: A ggplot2 raster plot of the potential landscape U,
#'          retained for compatibility. Use [autoplot()] to access it.
#'        - `vf`: The input vector field object.
#'        - `ss`: The steady-state distribution matrix.
#'
#' @export
make_2d_ld <- function(vf,
                       linear_interp = TRUE,
                       n_grid = 100,
                       drift_scheme = c("upwind", "central"),
                       cross_diffusion_mode = c("drop", "full"),
                       boundary_mode = c("reflective", "legacy_implicit")) {
  drift_scheme <- match.arg(drift_scheme)
  cross_diffusion_mode <- match.arg(cross_diffusion_mode)
  boundary_mode <- match.arg(boundary_mode)

  ss_raw <- ss_fp_2d(
    vf,
    linear_interp = linear_interp,
    n_grid = n_grid,
    drift_scheme = drift_scheme,
    cross_diffusion_mode = cross_diffusion_mode,
    boundary_mode = boundary_mode
  )

  # Preserve solver attributes before any post-processing.
  x_coords <- attr(ss_raw, "x_coords")
  y_coords <- attr(ss_raw, "y_coords")
  M <- attr(ss_raw, "M")

  nonfinite_mask <- !is.finite(ss_raw)
  nonpos_mask <- is.finite(ss_raw) & ss_raw <= 0
  correction_mask <- nonfinite_mask | nonpos_mask

  # Backup correction: remove problematic cells and renormalize remaining mass.
  # This keeps corrected cells out of the density support instead of shifting
  # the full surface.
  ss <- ss_raw
  if (any(correction_mask)) {
    ss[nonfinite_mask] <- 0
    ss[is.finite(ss) & ss <= 0] <- 0

    remaining_mass <- sum(ss, na.rm = TRUE)
    if (remaining_mass > 0) {
      ss <- ss / remaining_mass
    } else {
      # Degenerate safeguard: if everything is invalid, use a uniform fallback.
      ss[] <- 1 / length(ss)
    }
  }

  ss[!is.finite(ss)] <- 0
  ss[ss < 0] <- 0

  # For plotting, blank cells that were numerically problematic in the raw
  # solver output (non-finite or non-positive).
  bad_mask <- correction_mask
  n_bad <- sum(bad_mask)

  n_nonfinite <- sum(nonfinite_mask)
  n_nonpos <- sum(nonpos_mask)
  n_problem <- n_nonfinite + n_nonpos

  U <- matrix(NA_real_, nrow = nrow(ss), ncol = ncol(ss))
  U[ss > 0] <- -log(ss[ss > 0])
  U_plot <- U
  U_plot[bad_mask] <- NA_real_

  if (n_problem > 0 || n_bad > 0) {
    total_cells <- length(ss_raw)
    pct_problem <- 100 * n_problem / total_cells
    pct_bad <- 100 * n_bad / total_cells
    strong_warn_threshold <- max(5L, ceiling(0.01 * total_cells))

    bad_idx <- which(bad_mask, arr.ind = TRUE)
    n_boundary_bad <- if (length(bad_idx) > 0) {
      sum(
        bad_idx[, 1] %in% c(1L, nrow(ss_raw)) |
          bad_idx[, 2] %in% c(1L, ncol(ss_raw))
      )
    } else {
      0L
    }

    if (n_problem < strong_warn_threshold) {
      cli::cli_warn(paste0(
        "Detected ", n_problem, " numerically problematic grid cell(s) (",
        sprintf("%.2f", pct_problem),
        "% of grid) with non-finite/non-positive steady-state density.",
        " Plotting blanks ", n_bad, " corrected/non-finite cell(s) (",
        sprintf("%.2f", pct_bad), "% of grid)",
        if (n_boundary_bad > 0) {
          paste0(" (", n_boundary_bad, " at the boundary)")
        } else {
          ""
        },
        ". Breakdown: ", n_nonpos, " non-positive, ", n_nonfinite, " non-finite",
        ". These pixels are shown as blank in landscape plots due to local numerical issues; ",
        "this usually does not affect the main results."
      ))
    } else {
      cli::cli_warn(paste0(
        "Detected ", n_problem, " numerically problematic grid cell(s) (",
        sprintf("%.2f", pct_problem),
        "% of grid) with non-finite/non-positive steady-state density.",
        " Plotting blanks ", n_bad, " corrected/non-finite cell(s) (",
        sprintf("%.2f", pct_bad), "% of grid)",
        if (n_boundary_bad > 0) {
          paste0(" (", n_boundary_bad, " at the boundary)")
        } else {
          ""
        },
        ". Breakdown: ", n_nonpos, " non-positive, ", n_nonfinite, " non-finite",
        ". These pixels are shown as blank in landscape plots, but the count is high and may indicate ",
        "broader numerical instability. Please inspect this fit more closely."
      ))
    }
  }

  # Make a regular data frame for plotting. It contains x, y, d, U. (d: steady-state distribution, ss)
  dist <- expand.grid(x = x_coords, y = y_coords)
  dist$d <- as.vector(ss)
  dist$U <- as.vector(U)
  dist$U_plot <- as.vector(U_plot)

  plot <- plotly::plot_ly(
    data = dist,
    x = x_coords, y = y_coords, z = U_plot,
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
    ggplot2::geom_raster(ggplot2::aes(fill = U_plot)) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(
      x = vf$x,
      y = vf$y, fill = "U"
    ) +
    ggplot2::theme_bw()

  attr(ss, "x_coords") <- x_coords
  attr(ss, "y_coords") <- y_coords
  attr(ss, "M") <- M

  return(structure(list(
    dist = dist,
    plot = plot,
    plot_2 = plot_2,
    vf = vf,
    ss = ss
  ), class = c("2d_static_ld", "2d_ld", "landscape")))
}


is_inside_convex_hull <- function(points, hull_vertices, tol = 1e-10) {
  if (nrow(points) == 0) {
    return(logical(0))
  }
  if (nrow(hull_vertices) < 3) {
    return(rep(TRUE, nrow(points)))
  }

  x1 <- hull_vertices[, 1]
  y1 <- hull_vertices[, 2]
  x2 <- c(x1[-1], x1[1])
  y2 <- c(y1[-1], y1[1])

  out <- logical(nrow(points))
  for (k in seq_len(nrow(points))) {
    px <- points[k, 1]
    py <- points[k, 2]
    cross <- (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1)
    non_zero <- cross[abs(cross) > tol]
    out[k] <- length(non_zero) == 0 || all(non_zero >= 0) || all(non_zero <= 0)
  }

  out
}


#' Finds local minima in a 2D landscape object.
#'
#' @param ld A `2d_ld` or `2d_static_ld` object representing the landscape.
#' @param min_barrier_fraction When the barrier height between a local minimum
#'   and any of the saddle points connecting it to other local minima is less
#'   than `min_barrier_fraction` times the highest barrier height value, the
#'   local minimum will be considered minor. Default is 0.1.
#' @param min_convex_hull_range_fraction Additional lower bound for barrier-based
#'   retention. In 2D, this is expressed as a fraction of the potential range
#'   within the observed-data convex hull. In 1D, it is expressed as a fraction
#'   of the potential range over the observed data range. Default is 0.01.
#' @param exclude_minor Logical indicating whether to mark minor local minima
#' based on the barrier height criterion, so that they can be easily excluded
#' from subsequent calculations. Default is TRUE.
#' @param use_convex_hull Logical indicating whether to mark minima outside
#' the convex hull of observed data points as minor before barrier-based
#' minor-minimum detection. Default is TRUE.
#' @return An object of class `ld_min`. Its `mins` component contains one row
#'   per local minimum with coordinates, potential, `is_minor`, and
#'   `exclusion_reason`. The exclusion reason is one of
#'   `"retained_major"`, `"retained_major_outside_observed_data_convex_hull"`,
#'   `"outside_observed_data_convex_hull"`, `"insufficient_barrier_separation"`,
#'   or both minor criteria joined by `"; "`.
#' @export
find_loc_min <- function(ld,
                         exclude_minor = TRUE,
                         min_barrier_fraction = 0.1,
                         min_convex_hull_range_fraction = 0.01,
                         use_convex_hull = TRUE) {
  if (inherits(ld, "1d_ld") || inherits(ld, "1d_static_ld")) {
    return(find_loc_min_1d(
      ld,
      exclude_minor = exclude_minor,
      min_barrier_fraction = min_barrier_fraction,
      min_convex_hull_range_fraction = min_convex_hull_range_fraction
    ))
  }

  if (!inherits(ld, "2d_ld") && !inherits(ld, "2d_static_ld")) {
    cli::cli_abort("Input {.arg ld} must be a {.cls 1d_ld}, {.cls 1d_static_ld}, {.cls 2d_ld}, or {.cls 2d_static_ld} object.")
  }

  dist <- ld$dist
  U_matrix <- matrix(dist$U, nrow = length(unique(dist$x)), ncol = length(unique(dist$y)))

  n_x <- nrow(U_matrix)
  n_y <- ncol(U_matrix)

  local_mins <- data.frame(x = numeric(0), y = numeric(0), U = numeric(0))
  if (n_x >= 1 && n_y >= 1) {
    is_min <- matrix(FALSE, nrow = n_x, ncol = n_y)

    for (i in seq_len(n_x)) {
      for (j in seq_len(n_y)) {
        neighbor_i <- max(1L, i - 1L):min(n_x, i + 1L)
        neighbor_j <- max(1L, j - 1L):min(n_y, j + 1L)
        neighbor_coords <- expand.grid(ii = neighbor_i, jj = neighbor_j)
        neighbor_coords <- neighbor_coords[!(neighbor_coords$ii == i & neighbor_coords$jj == j), , drop = FALSE]

        if (nrow(neighbor_coords) == 0) {
          next
        }

        neighbor_u <- U_matrix[cbind(neighbor_coords$ii, neighbor_coords$jj)]
        is_min[i, j] <- all(U_matrix[i, j] < neighbor_u)
      }
    }

    idx <- which(is_min, arr.ind = TRUE)
    if (nrow(idx) > 0) {
      idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
      ux <- unique(dist$x)
      uy <- unique(dist$y)
      local_mins <- data.frame(
        x = ux[idx[, 1]],
        y = uy[idx[, 2]],
        U = U_matrix[idx],
        row.names = NULL
      )
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

  hull_vertices <- NULL
  hull_range_threshold <- 0
  if (exclude_minor && isTRUE(use_convex_hull) && n_mins > 0 && !is.null(ld$vf) && !is.null(ld$vf$data)) {
    data_xy <- as.matrix(ld$vf$data[, c(ld$vf$x, ld$vf$y), drop = FALSE])
    data_xy <- data_xy[stats::complete.cases(data_xy), , drop = FALSE]
      if (nrow(data_xy) >= 3) {
        hull_idx <- grDevices::chull(data_xy[, 1], data_xy[, 2])
        hull_vertices <- data_xy[hull_idx, , drop = FALSE]

        grid_xy <- as.matrix(dist[, c("x", "y"), drop = FALSE])
        grid_inside_hull <- is_inside_convex_hull(grid_xy, hull_vertices)
        hull_u <- dist$U[grid_inside_hull]
        hull_u <- hull_u[is.finite(hull_u)]
        if (length(hull_u) > 0) {
          hull_range_threshold <- min_convex_hull_range_fraction * (max(hull_u) - min(hull_u))
        }
      }
  }

  barrier_minor_mins <- integer(0)
  if (n_mins <= 1 || exclude_minor == FALSE) {
    all_barriers <- matrix(NA, nrow = n_mins, ncol = n_mins)
    minor_mins <- integer(0)
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
    minor_mins <- integer(0)
    all_barriers_copy <- all_barriers

    finite_barriers <- all_barriers_copy[is.finite(all_barriers_copy)]
    max_barrier <- if (length(finite_barriers) > 0) max(finite_barriers) else Inf


    # do the following until no remaining barriers are lower than the
    # relative highest-barrier threshold and the convex-hull-range threshold
    # first, find the lowest barrier
    # label this minimum as minor
    # remove this minimum from the matrix (set the corresponding row and column to NA)
    # repeat

    repeat {
      current_barriers <- all_barriers_copy[is.finite(all_barriers_copy)]
      if (length(current_barriers) == 0) {
        break
      }
      current_min_barrier <- min(current_barriers)
      barrier_threshold <- max(
        min_barrier_fraction * max_barrier,
        hull_range_threshold
      )
      if (is.infinite(current_min_barrier) || current_min_barrier >= barrier_threshold) {
        break
      }
      locs <- which(all_barriers_copy == current_min_barrier, arr.ind = TRUE)
      min_to_remove <- locs[1, 1] # arbitrarily choose the first
      barrier_minor_mins <- unique(c(barrier_minor_mins, min_to_remove))
      minor_mins <- unique(c(minor_mins, min_to_remove))
      all_barriers_copy[min_to_remove, ] <- NA
      all_barriers_copy[, min_to_remove] <- NA
    }
  }

  hull_minor_mins <- integer(0)
  hull_outside_but_retained <- integer(0)
  if (exclude_minor && !is.null(hull_vertices)) {
    mins_xy <- as.matrix(local_mins[, c("x", "y"), drop = FALSE])
    inside <- is_inside_convex_hull(mins_xy, hull_vertices)
    outside_idx <- which(!inside)
    if (length(outside_idx) > 0) {
      major_after_barrier <- setdiff(seq_len(n_mins), barrier_minor_mins)
      outside_major_after_barrier <- intersect(major_after_barrier, outside_idx)
      if (length(major_after_barrier) == 1L && major_after_barrier %in% outside_idx) {
        hull_outside_but_retained <- major_after_barrier
        hull_minor_mins <- setdiff(outside_idx, major_after_barrier)
      } else {
        hull_minor_mins <- outside_idx
      }
    }
  }

  minor_mins <- unique(c(barrier_minor_mins, hull_minor_mins))

  exclusion_reason <- rep("retained_major", n_mins)
  if (n_mins > 0) {
    reason_parts <- rep("", n_mins)
    if (length(hull_outside_but_retained) > 0) {
      reason_parts[hull_outside_but_retained] <- "retained_major_outside_observed_data_convex_hull"
    }
    if (length(hull_minor_mins) > 0) {
      reason_parts[hull_minor_mins] <- "outside_observed_data_convex_hull"
    }
    if (length(barrier_minor_mins) > 0) {
      reason_parts[barrier_minor_mins] <- ifelse(
        nzchar(reason_parts[barrier_minor_mins]),
        paste(reason_parts[barrier_minor_mins], "insufficient_barrier_separation", sep = "; "),
        "insufficient_barrier_separation"
      )
    }
    exclusion_reason[nzchar(reason_parts)] <- reason_parts[nzchar(reason_parts)]
  }

  local_mins <- local_mins %>%
    dplyr::mutate(
      is_minor = ifelse(dplyr::row_number() %in% minor_mins, TRUE, FALSE),
      exclusion_reason = exclusion_reason
    )

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
