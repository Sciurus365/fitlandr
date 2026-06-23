#' @export
autoplot.summary_bootstrap_2d_ld <- function(object,
                                             mode = c("clusters", "minima"),
                                             show_ellipses = TRUE,
                                             show_original_major_minima = TRUE,
                                             point_alpha = 0.35,
                                             ellipse_alpha = 0.15,
                                             ...) {
  if (!inherits(object, "summary_bootstrap_2d_ld")) {
    cli::cli_abort("{.arg object} must inherit from {.cls summary_bootstrap_2d_ld}.")
  }
  if (is.null(object$per_point)) {
    return(ggplot2::ggplot())
  }
  mode <- rlang::arg_match0(mode, c("clusters", "minima"))
  df_points <- object$per_point
  df_cl <- object$per_cluster

  if (mode == "minima") {
    x_range <- range(attr(object$original_ld$ss, "x_coords"))
    y_range <- range(attr(object$original_ld$ss, "y_coords"))
    x_jitter_amount <- (x_range[2] - x_range[1]) * 0.02
    y_jitter_amount <- (y_range[2] - y_range[1]) * 0.02

    boot_counts <- table(df_points$boot_index)
    boot_count_df <- data.frame(
      boot_index = as.integer(names(boot_counts)),
      count = as.integer(boot_counts)
    )

    plot_df <- df_points %>%
      dplyr::left_join(boot_count_df, by = "boot_index")

    max_count <- max(plot_df$count)
    if (max_count > 4) {
      plot_df <- plot_df %>%
        dplyr::mutate(count = ifelse(count >= 4, 4, count))
      cli::cli_inform("More than 4 minima found in some bootstrap samples.
                  Combining counts >=4 into a single category for plotting.")
    }

    return(
      ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = x + stats::runif(nrow(plot_df), -x_jitter_amount, x_jitter_amount),
          y = y + stats::runif(nrow(plot_df), -y_jitter_amount, y_jitter_amount),
          color = U,
          shape = as.factor(count)
        )
      ) +
        ggplot2::geom_point(size = 2, alpha = 0.7) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::labs(
          x = object$original_ld$vf$x,
          y = object$original_ld$vf$y,
          color = "U value",
          shape = "Number of minima\nin bootstrap sample"
        ) +
        ggplot2::theme_bw() +
        ggplot2::scale_shape_manual(values = c(16, 17, 15, 18, 3)) +
        ggplot2::coord_fixed(xlim = x_range, ylim = y_range)
    )
  }

  x_range <- range(attr(object$original_ld$ss, "x_coords"))
  y_range <- range(attr(object$original_ld$ss, "y_coords"))

  p <- ggplot2::ggplot(df_points) +
    ggplot2::geom_point(
      ggplot2::aes(x = x, y = y, color = factor(cluster)),
      alpha = point_alpha,
      size = 1
    ) +
    ggplot2::coord_fixed(xlim = x_range, ylim = y_range) +
    ggplot2::labs(color = "cluster")

  if (!is.null(df_cl) && nrow(df_cl)) {
    if (show_ellipses) {
      p <- p + ggforce::geom_ellipse(
        data = df_cl,
        ggplot2::aes(
          x0 = mean_x, y0 = mean_y,
          a = a_pred, b = b_pred, angle = angle
        ),
        color = "firebrick", fill = "firebrick", alpha = ellipse_alpha
      )
    }
  }

  if (isTRUE(show_original_major_minima) && !is.null(object$original_ld)) {
    min_barrier_fraction <- object$params$min_barrier_fraction
    min_convex_hull_range_fraction <- object$params$min_convex_hull_range_fraction
    if (is.null(min_barrier_fraction) || !is.finite(min_barrier_fraction)) {
      min_barrier_fraction <- 0.1
    }
    if (is.null(min_convex_hull_range_fraction) || !is.finite(min_convex_hull_range_fraction)) {
      min_convex_hull_range_fraction <- 0.01
    }

    orig_major <- tryCatch(
      {
        mins <- find_loc_min(
          object$original_ld,
          exclude_minor = TRUE,
          min_barrier_fraction = min_barrier_fraction,
          min_convex_hull_range_fraction = min_convex_hull_range_fraction
        )$mins
        if (!is.null(mins) && nrow(mins) && "is_minor" %in% names(mins)) {
          mins <- mins[!mins$is_minor, , drop = FALSE]
        }
        if (is.null(mins) || !nrow(mins)) {
          NULL
        } else {
          mins
        }
      },
      error = function(e) NULL
    )

    if (!is.null(orig_major)) {
      p <- p +
        ggplot2::geom_point(
          data = orig_major,
          ggplot2::aes(x = x, y = y),
          inherit.aes = FALSE,
          shape = 4,
          stroke = 1.1,
          size = 3,
          color = "black"
        )
    }
  }
  p + ggplot2::theme_bw()
}
