#' Compare the potential values of two bootstrap minima
#'
#' Computes a paired bootstrap interval for the potential difference between
#' two selected minima in a summarized bootstrap landscape. Pairing by
#' bootstrap replication makes the comparison invariant to the arbitrary
#' additive constant of each bootstrap potential.
#'
#' @param object A `summary_bootstrap_2d_ld` object.
#' @param minima Integer vector of length two giving cluster identifiers. The
#'   reported contrast is the first minimum minus the second, `U[1] - U[2]`.
#' @param level Coverage level for the bootstrap interval. By default, uses the
#'   level stored in `object$params`, or 0.95 if it is unavailable.
#'
#' @return An object of class `minima_depth_comparison` with components
#'   `result`, a one-row data frame containing the potential contrast and its
#'   interval, and `per_boot`, the paired bootstrap contrasts.
#' @export
compare_minima_depths <- function(object, minima, level = NULL) {
  if (!inherits(object, "summary_bootstrap_2d_ld")) {
    cli::cli_abort("{.arg object} must inherit from {.cls summary_bootstrap_2d_ld}.")
  }
  if (length(minima) != 2L || anyNA(minima) || length(unique(minima)) != 2L) {
    cli::cli_abort("{.arg minima} must contain two distinct cluster identifiers.")
  }
  if (is.null(object$per_point) || !nrow(object$per_point)) {
    cli::cli_abort("{.arg object} does not contain bootstrap minima in {.field per_point}.")
  }
  if (is.null(object$per_cluster) || !nrow(object$per_cluster)) {
    cli::cli_abort("{.arg object} does not contain summarized minima in {.field per_cluster}.")
  }

  minima <- as.integer(minima)
  available <- as.integer(object$per_cluster$cluster)
  missing_minima <- setdiff(minima, available)
  if (length(missing_minima)) {
    cli::cli_abort(c(
      "Selected minima are not present in {.field per_cluster}.",
      "x" = "Missing cluster identifier{?s}: {paste(missing_minima, collapse = ', ')}."
    ))
  }

  if (is.null(level)) {
    level <- object$params$level
    if (is.null(level)) {
      level <- 0.95
    }
  }
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    cli::cli_abort("{.arg level} must be one number strictly between 0 and 1.")
  }

  points <- object$per_point
  if ("is_noise" %in% names(points)) {
    points <- points[!points$is_noise, , drop = FALSE]
  }
  points <- points[points$cluster %in% minima, , drop = FALSE]

  per_minimum <- points |>
    dplyr::group_by(.data$boot_index, .data$cluster) |>
    dplyr::summarise(U = mean(.data$U), .groups = "drop")

  first <- per_minimum[per_minimum$cluster == minima[[1]], c("boot_index", "U")]
  second <- per_minimum[per_minimum$cluster == minima[[2]], c("boot_index", "U")]
  names(first)[[2]] <- "U_first"
  names(second)[[2]] <- "U_second"

  per_boot <- dplyr::inner_join(first, second, by = "boot_index") |>
    dplyr::mutate(delta_U = .data$U_first - .data$U_second) |>
    dplyr::arrange(.data$boot_index)

  if (!nrow(per_boot)) {
    cli::cli_abort("The selected minima do not co-occur in any bootstrap replication.")
  }

  alpha <- (1 - level) / 2
  interval <- stats::quantile(
    per_boot$delta_U,
    probs = c(alpha, 1 - alpha),
    names = FALSE,
    na.rm = TRUE
  )

  cluster_result <- object$per_cluster[
    match(minima, object$per_cluster$cluster),
    ,
    drop = FALSE
  ]
  result <- data.frame(
    minimum_first = minima[[1]],
    minimum_second = minima[[2]],
    delta_U = cluster_result$mean_U[[1]] - cluster_result$mean_U[[2]],
    bootstrap_mean = mean(per_boot$delta_U),
    CI_lower = interval[[1]],
    CI_upper = interval[[2]],
    level = level,
    n_paired = nrow(per_boot)
  )

  structure(
    list(result = result, per_boot = per_boot),
    class = "minima_depth_comparison"
  )
}


#' @export
print.minima_depth_comparison <- function(x, ...) {
  print(x$result, row.names = FALSE)
  invisible(x)
}
