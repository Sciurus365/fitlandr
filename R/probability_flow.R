#' Generates probability-flow vectors
#' for a given vector field and landscape.
#'
#' @param vf A `vectorfield` or `cv_vectorfield` object representing the vector field.
#' @param ld A `2d_ld` or `2d_static_ld` object representing the landscape.
#' @param n Number of flow vectors to generate along each dimension (default 20).
#' @param divided_by_rho Logical indicating whether to divide flow vectors by steady-state density.
#' @param cross_diffusion_mode Cross-diffusion handling mode. By default, this
#'   is inherited from `ld`. If supplied, it must match the mode used to
#'   estimate the landscape.
#' @return An object of class `2d_pf` containing:
#'         - `vec_grid`: Data frame with columns x, y, vx, vy.
#'         - `face_grid`: Staggered conservative face currents.
#'         - `vf`: The input vector field object.
#'         - `ld`: The input landscape object.
#' @export
make_2d_pf <- function(vf, ld, n = 20, divided_by_rho = FALSE,
                       cross_diffusion_mode = NULL) {
  if (inherits(vf, "cv_vectorfield")) {
    vf <- vf$final_model
  } else if (!inherits(vf, "vectorfield")) {
    cli::cli_abort("Input {.arg vf} must be a {.cls vectorfield} or {.cls cv_vectorfield} object.")
  }
  if (!inherits(ld, "2d_static_ld")) {
    cli::cli_abort("Input {.arg ld} must be a {.cls 2d_static_ld} object.")
  }
  if (is.null(ld$cross_diffusion_mode) && is.null(cross_diffusion_mode)) {
    cli::cli_abort(c(
      "The landscape does not record its cross-diffusion mode.",
      "i" = "Recreate {.arg ld} with the current version of {.fn make_2d_ld}, or supply {.arg cross_diffusion_mode} explicitly."
    ))
  }
  if (is.null(cross_diffusion_mode)) {
    cross_diffusion_mode <- ld$cross_diffusion_mode
  } else {
    cross_diffusion_mode <- match.arg(cross_diffusion_mode, c("drop", "full"))
    if (!is.null(ld$cross_diffusion_mode) &&
        !identical(cross_diffusion_mode, ld$cross_diffusion_mode)) {
      cli::cli_abort("{.arg cross_diffusion_mode} must match the mode used to estimate {.arg ld}.")
    }
  }

  make_2d_pf_fvm(
    vf = vf,
    ld = ld,
    n = n,
    divided_by_rho = divided_by_rho,
    cross_diffusion_mode = cross_diffusion_mode
  )
}
make_2d_pf_fvm <- function(vf, ld, n, divided_by_rho, cross_diffusion_mode) {
  if (!identical(cross_diffusion_mode, "drop")) {
    cli::cli_abort(c(
      "The compatible FVM probability flow currently requires {.code cross_diffusion_mode = \"drop\"}.",
      "i" = "The full cross-diffusion stencil does not yet expose one conservative current per shared face.",
      "i" = "Recreate the landscape with the default cross-diffusion mode."
    ))
  }
  if (!isTRUE(ld$fvm_compatible) || is.null(ld$fvm_faces)) {
    cli::cli_abort(c(
      "The landscape does not contain an unmodified compatible FVM face current.",
      "i" = "Recreate {.arg ld} with the current version of {.fn make_2d_ld}.",
      "i" = "If its stationary density required numerical correction, inspect that fit."
    ))
  }
  if (length(n) != 1L || !is.finite(n) || n != as.integer(n) || n < 2L) {
    cli::cli_abort("{.arg n} must be an integer of at least 2.")
  }

  faces <- ld$fvm_faces
  required <- c(
    "Jx", "Jy", "x_faces", "y_faces", "x_centers", "y_centers", "hx", "hy"
  )
  if (!all(required %in% names(faces))) {
    cli::cli_abort("The landscape contains incomplete finite-volume face-current data.")
  }
  nx <- length(faces$x_centers)
  ny <- length(faces$y_centers)
  if (!identical(dim(faces$Jx), c(nx + 1L, ny)) ||
      !identical(dim(faces$Jy), c(nx, ny + 1L)) ||
      any(!is.finite(c(faces$Jx, faces$Jy)))) {
    cli::cli_abort("The landscape contains invalid finite-volume face-current data.")
  }

  center_vx <- (faces$Jx[seq_len(nx), , drop = FALSE] +
    faces$Jx[seq_len(nx) + 1L, , drop = FALSE]) / 2
  center_vy <- (faces$Jy[, seq_len(ny), drop = FALSE] +
    faces$Jy[, seq_len(ny) + 1L, drop = FALSE]) / 2
  if (isTRUE(divided_by_rho)) {
    center_vx <- center_vx / ld$ss
    center_vy <- center_vy / ld$ss
  }

  idx_x <- unique(round(seq(1L, nx, length.out = n)))
  idx_y <- unique(round(seq(1L, ny, length.out = n)))
  grid_index <- expand.grid(i = idx_x, j = idx_y)
  vec_grid <- transform(
    grid_index,
    x = faces$x_centers[grid_index$i],
    y = faces$y_centers[grid_index$j],
    vx = center_vx[cbind(grid_index$i, grid_index$j)],
    vy = center_vy[cbind(grid_index$i, grid_index$j)]
  )[c("x", "y", "vx", "vy")]

  vertical_faces <- expand.grid(x = faces$x_faces, y = faces$y_centers)
  horizontal_faces <- expand.grid(x = faces$x_centers, y = faces$y_faces)
  face_grid <- rbind(
    transform(vertical_faces, component = "x", current = as.numeric(faces$Jx)),
    transform(horizontal_faces, component = "y", current = as.numeric(faces$Jy))
  )

  structure(list(
    vec_grid = vec_grid,
    face_grid = face_grid,
    fvm_faces = faces,
    x = vf$x,
    y = vf$y,
    vf = vf,
    ld = ld,
    divided_by_rho = divided_by_rho,
    cross_diffusion_mode = cross_diffusion_mode,
    method = "fvm"
  ), class = c("2d_pf", "probabilityflow", "vectorfield"))
}


#' Estimate a stream function from a two-dimensional probability flow
#'
#' Estimates a scalar stream function `A` whose perpendicular gradient
#' `(-dA/dy, dA/dx)` is closest in least-squares distance to the estimated
#' probability flow. By default, the function uses the conservative face
#' currents retained by [make_2d_ld()] and a compatible staggered-grid
#' perpendicular-gradient operator. This uses exactly the same faces and
#' reflective boundary conditions as the finite-volume stationary-density
#' solve.
#'
#' A stream function is identifiable only up to an additive constant. The
#' returned solution is anchored by setting its value at the first grid corner
#' to zero.
#'
#' @param pf A `2d_pf` object returned by [make_2d_pf()].
#'
#' @return An object of class `2d_stream` containing:
#'   - `grid`: cell-center coordinates and estimated stream-function values.
#'   - `corner_grid`: corner coordinates and staggered stream-function values
#'     on the finite-volume grid.
#'   - `face_grid`: observed, fitted, and residual face currents for
#'     the finite-volume grid.
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

  make_2d_stream_fvm(pf)
}


make_2d_stream_fvm <- function(pf) {
  if (!identical(pf$method, "fvm") || is.null(pf$fvm_faces)) {
    cli::cli_abort(c(
      "The compatible stream method requires an FVM probability-flow object.",
      "i" = "Recreate {.arg pf} with {.fn make_2d_pf}."
    ))
  }
  if (!identical(pf$cross_diffusion_mode, "drop")) {
    cli::cli_abort(c(
      "The compatible FVM stream function currently requires {.code cross_diffusion_mode = \"drop\"}.",
      "i" = "The full cross-diffusion stencil does not yet expose one conservative current per shared face.",
      "i" = "Recreate the landscape with the default cross-diffusion mode."
    ))
  }

  faces <- pf$fvm_faces
  required <- c(
    "Jx", "Jy", "x_faces", "y_faces", "x_centers", "y_centers", "hx", "hy"
  )
  if (!all(required %in% names(faces))) {
    cli::cli_abort("The landscape contains incomplete finite-volume face-current data.")
  }
  nx <- length(faces$x_centers)
  ny <- length(faces$y_centers)
  if (!identical(dim(faces$Jx), c(nx + 1L, ny)) ||
      !identical(dim(faces$Jy), c(nx, ny + 1L)) ||
      any(!is.finite(c(faces$Jx, faces$Jy)))) {
    cli::cli_abort("The landscape contains invalid finite-volume face-current data.")
  }

  n_x_faces <- (nx + 1L) * ny
  n_y_faces <- nx * (ny + 1L)
  n_faces <- n_x_faces + n_y_faces
  n_corners <- (nx + 1L) * (ny + 1L)
  x_face_idx <- function(i, j) i + (j - 1L) * (nx + 1L)
  y_face_idx <- function(i, j) n_x_faces + i + (j - 1L) * nx
  corner_idx <- function(i, j) i + (j - 1L) * (nx + 1L)

  row_idx <- integer(2L * n_faces)
  col_idx <- integer(2L * n_faces)
  values <- numeric(2L * n_faces)
  pos <- 0L
  for (j in seq_len(ny)) {
    for (i in seq_len(nx + 1L)) {
      idx <- pos + seq_len(2L)
      row_idx[idx] <- x_face_idx(i, j)
      col_idx[idx] <- c(corner_idx(i, j), corner_idx(i, j + 1L))
      values[idx] <- c(1 / faces$hy, -1 / faces$hy)
      pos <- pos + 2L
    }
  }
  for (j in seq_len(ny + 1L)) {
    for (i in seq_len(nx)) {
      idx <- pos + seq_len(2L)
      row_idx[idx] <- y_face_idx(i, j)
      col_idx[idx] <- c(corner_idx(i, j), corner_idx(i + 1L, j))
      values[idx] <- c(-1 / faces$hx, 1 / faces$hx)
      pos <- pos + 2L
    }
  }
  perpendicular_gradient <- Matrix::sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = values,
    dims = c(n_faces, n_corners)
  )

  observed <- c(as.numeric(faces$Jx), as.numeric(faces$Jy))
  anchor <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 1, dims = c(1L, n_corners))
  system_matrix <- rbind(perpendicular_gradient, anchor)
  stream_values <- as.numeric(Matrix::solve(Matrix::qr(system_matrix), c(observed, 0)))
  fitted <- as.numeric(perpendicular_gradient %*% stream_values)
  residual <- observed - fitted

  A_corner <- matrix(stream_values, nrow = nx + 1L, ncol = ny + 1L)
  A_center <- (
    A_corner[seq_len(nx), seq_len(ny), drop = FALSE] +
      A_corner[seq_len(nx) + 1L, seq_len(ny), drop = FALSE] +
      A_corner[seq_len(nx), seq_len(ny) + 1L, drop = FALSE] +
      A_corner[seq_len(nx) + 1L, seq_len(ny) + 1L, drop = FALSE]
  ) / 4

  fitted_Jx <- matrix(fitted[seq_len(n_x_faces)], nrow = nx + 1L, ncol = ny)
  fitted_Jy <- matrix(fitted[n_x_faces + seq_len(n_y_faces)], nrow = nx, ncol = ny + 1L)
  residual_Jx <- faces$Jx - fitted_Jx
  residual_Jy <- faces$Jy - fitted_Jy
  center_component <- function(x_faces, y_faces) {
    list(
      x = (x_faces[seq_len(nx), , drop = FALSE] +
        x_faces[seq_len(nx) + 1L, , drop = FALSE]) / 2,
      y = (y_faces[, seq_len(ny), drop = FALSE] +
        y_faces[, seq_len(ny) + 1L, drop = FALSE]) / 2
    )
  }
  observed_center <- center_component(faces$Jx, faces$Jy)
  fitted_center <- center_component(fitted_Jx, fitted_Jy)
  residual_center <- center_component(residual_Jx, residual_Jy)

  ordered_grid <- expand.grid(x = faces$x_centers, y = faces$y_centers)
  corner_grid <- expand.grid(x = faces$x_faces, y = faces$y_faces)
  vertical_faces <- expand.grid(x = faces$x_faces, y = faces$y_centers)
  horizontal_faces <- expand.grid(x = faces$x_centers, y = faces$y_faces)
  face_grid <- rbind(
    transform(
      vertical_faces,
      component = "x",
      current = as.numeric(faces$Jx),
      fitted = as.numeric(fitted_Jx),
      residual = as.numeric(residual_Jx)
    ),
    transform(
      horizontal_faces,
      component = "y",
      current = as.numeric(faces$Jy),
      fitted = as.numeric(fitted_Jy),
      residual = as.numeric(residual_Jy)
    )
  )
  observed_norm <- sqrt(sum(observed^2))

  structure(
    list(
      grid = transform(ordered_grid, A = as.numeric(A_center)),
      corner_grid = transform(corner_grid, A = stream_values),
      face_grid = face_grid,
      fitted_grid = transform(
        ordered_grid,
        vx = as.numeric(observed_center$x),
        vy = as.numeric(observed_center$y),
        fitted_vx = as.numeric(fitted_center$x),
        fitted_vy = as.numeric(fitted_center$y)
      ),
      residual_grid = transform(
        ordered_grid,
        vx = as.numeric(residual_center$x),
        vy = as.numeric(residual_center$y)
      ),
      rmse = sqrt(mean(residual^2)),
      relative_error = if (observed_norm == 0) 0 else sqrt(sum(residual^2)) / observed_norm,
      method = "fvm",
      pf = pf
    ),
    class = "2d_stream"
  )
}
