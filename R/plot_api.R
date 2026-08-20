#' Autoplot a landscape or probability flow
#'
#' These methods return the two-dimensional ggplot representation of a
#' landscape or probability-flow object.
#'
#' @param object A fitlandr object.
#' @param ... Additional arguments, currently unused.
#'
#' @return A ggplot object.
#' @importFrom ggplot2 autoplot
#' @export
autoplot.2d_static_ld <- function(object, ...) {
  object$plot_2
}

#' @rdname autoplot.2d_static_ld
#' @export
autoplot.1d_static_ld <- function(object, ...) {
  object$plot_2
}

#' @rdname autoplot.2d_static_ld
#' @export
autoplot.1d_ld <- function(object, ...) {
  object$plot_2
}

#' @rdname autoplot.2d_static_ld
#' @export
autoplot.2d_ld <- function(object, ...) {
  object$plot_2
}

#' @rdname autoplot.2d_static_ld
#' @export
autoplot.2d_MVKE_landscape <- function(object, ...) {
  object$plot
}

#' @rdname autoplot.2d_static_ld
#' @export
autoplot.2d_pf <- function(object, ...) {
  ggplot2::ggplot(object$vec_grid, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = x + vx, yend = y + vy),
      arrow = grid::arrow(length = grid::unit(0.1, "cm")),
      alpha = 0.7
    ) +
    ggplot2::labs(
      x = object$x,
      y = object$y,
      title = "Probability Flow Vectors"
    ) +
    ggplot2::theme_bw()
}

#' Autoplot a two-dimensional stream function
#'
#' @param object A `2d_stream` object returned by [make_2d_stream()].
#' @param contour Logical indicating whether to overlay contour lines.
#' @param ... Additional arguments, currently unused.
#'
#' @return A ggplot object.
#' @export
autoplot.2d_stream <- function(object, contour = TRUE, ...) {
  grid <- object$grid
  cell_bounds <- function(coords) {
    coords <- sort(unique(coords))
    edges <- c(
      coords[1L] - (coords[2L] - coords[1L]) / 2,
      (coords[-1L] + coords[-length(coords)]) / 2,
      coords[length(coords)] +
        (coords[length(coords)] - coords[length(coords) - 1L]) / 2
    )
    data.frame(
      coord = coords,
      lower = edges[-length(edges)],
      upper = edges[-1L]
    )
  }
  x_bounds <- cell_bounds(grid$x)
  y_bounds <- cell_bounds(grid$y)
  x_match <- match(grid$x, x_bounds$coord)
  y_match <- match(grid$y, y_bounds$coord)
  grid$cell_xmin <- x_bounds$lower[x_match]
  grid$cell_xmax <- x_bounds$upper[x_match]
  grid$cell_ymin <- y_bounds$lower[y_match]
  grid$cell_ymax <- y_bounds$upper[y_match]

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_rect(ggplot2::aes(
      fill = A,
      xmin = .data$cell_xmin,
      xmax = .data$cell_xmax,
      ymin = .data$cell_ymin,
      ymax = .data$cell_ymax
    ), inherit.aes = FALSE) +
    ggplot2::scale_fill_viridis_c(name = "A") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = object$pf$x,
      y = object$pf$y,
      title = "Probability-flow stream function"
    ) +
    ggplot2::theme_bw()

  if (isTRUE(contour)) {
    p <- p + ggplot2::geom_contour(
      ggplot2::aes(z = A),
      color = "white",
      alpha = 0.6,
      show.legend = FALSE
    )
  }

  p
}

#' @export
plotly_ld.2d_static_ld <- function(object, ...) {
  object$plot
}

#' @export
plotly_ld.1d_static_ld <- function(object, ...) {
  object$plot
}

#' @export
plotly_ld.1d_ld <- function(object, ...) {
  object$plot
}

#' @export
plotly_ld.2d_ld <- function(object, ...) {
  object$plot
}

deprecate_landscape_plot <- function(old, index) {
  replacement <- if (identical(index, 2) || identical(index, "2")) {
    "autoplot()"
  } else {
    "plotly_ld()"
  }
  lifecycle::deprecate_warn("0.2.0", old, replacement)
}

legacy_landscape_plot <- function(x, index, ...) {
  if (identical(index, 1) || identical(index, "1")) {
    plotly_ld(x, ...)
  } else if (identical(index, 2) || identical(index, "2")) {
    autoplot(x, ...)
  } else if (identical(index, 3) || identical(index, "3") ||
             identical(index, "mat_3d")) {
    graphics::plot(x$mat_3d, ...)
  } else {
    cli::cli_abort("{.arg index} must be 1, 2, 3, or {.val mat_3d}.")
  }
}

#' Deprecated landscape and probability-flow plot methods
#'
#' `r lifecycle::badge("deprecated")`
#'
#' Use [plotly_ld()] for interactive three-dimensional landscapes and
#' [autoplot()] for ggplot output.
#'
#' @param x A fitlandr object.
#' @param index The legacy landscape plot index.
#' @param ... Arguments passed to the replacement method.
#'
#' @export
plot.2d_static_ld <- function(x, index = 1, ...) {
  deprecate_landscape_plot("plot.2d_static_ld()", index)
  legacy_landscape_plot(x, index, ...)
}

#' @rdname plot.2d_static_ld
#' @export
plot.1d_static_ld <- function(x, index = 1, ...) {
  deprecate_landscape_plot("plot.1d_static_ld()", index)
  legacy_landscape_plot(x, index, ...)
}

#' @rdname plot.2d_static_ld
#' @export
plot.1d_ld <- function(x, index = 1, ...) {
  deprecate_landscape_plot("plot.1d_ld()", index)
  legacy_landscape_plot(x, index, ...)
}

#' @rdname plot.2d_static_ld
#' @export
plot.2d_ld <- function(x, index = 1, ...) {
  deprecate_landscape_plot("plot.2d_ld()", index)
  legacy_landscape_plot(x, index, ...)
}

#' @rdname plot.2d_static_ld
#' @export
plot.2d_MVKE_landscape <- function(x, index = 1, ...) {
  lifecycle::deprecate_warn(
    "0.2.0",
    "plot.2d_MVKE_landscape()",
    "autoplot()"
  )
  autoplot(x, ...)
}

#' @rdname plot.2d_static_ld
#' @export
plot.2d_pf <- function(x, ...) {
  lifecycle::deprecate_warn(
    "0.2.0",
    "plot.2d_pf()",
    "autoplot()"
  )
  autoplot(x, ...)
}

#' @importFrom simlandr plotly_ld
#' @export
simlandr::plotly_ld

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
