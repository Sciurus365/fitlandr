#' Insert separator rows between invalid time transitions
#'
#' Internal helper to prevent lag vectors from being formed across day
#' boundaries or non-consecutive beeps.
#'
#' @param data A data frame or matrix.
#' @param columns Character vector of columns to retain.
#' @param dayvar Optional character scalar naming the day variable.
#' @param beepvar Optional character scalar naming the beep variable.
#'
#' @return A data frame containing the retained columns plus optional
#'   `dayvar` and `beepvar`, with all-`NA` separator rows inserted where
#'   transitions should be broken.
#' @keywords internal
insert_time_separators <- function(data, columns, dayvar = NULL, beepvar = NULL) {
  keep <- unique(c(columns, dayvar, beepvar))
  out <- tibble::as_tibble(data[, keep, drop = FALSE])

  if (nrow(out) <= 1L || (is.null(dayvar) && is.null(beepvar))) {
    return(out)
  }

  break_after <- rep(FALSE, nrow(out) - 1L)
  state_missing <- !stats::complete.cases(out[, columns, drop = FALSE])
  transition_has_state_missing <- state_missing[-nrow(out)] | state_missing[-1L]

  add_time_breaks <- function(invalid_transition) {
    # Existing state-NA rows already separate vectors. A missing time value on
    # an otherwise observed pair is itself treated as an invalid transition.
    invalid_transition[is.na(invalid_transition)] <-
      !transition_has_state_missing[is.na(invalid_transition)]
    break_after <<- break_after | invalid_transition
  }

  if (!is.null(dayvar)) {
    add_time_breaks(out[[dayvar]][-nrow(out)] != out[[dayvar]][-1L])
  }

  if (!is.null(beepvar)) {
    add_time_breaks((out[[beepvar]][-1L] - out[[beepvar]][-nrow(out)]) != 1)
  }

  if (!any(break_after)) {
    return(out)
  }

  pieces <- vector("list", length = sum(break_after) * 2L + 1L)
  piece_i <- 1L
  start_i <- 1L

  for (i in which(break_after)) {
    pieces[[piece_i]] <- out[start_i:i, , drop = FALSE]
    piece_i <- piece_i + 1L
    pieces[[piece_i]] <- out[NA_integer_, , drop = FALSE]
    piece_i <- piece_i + 1L
    start_i <- i + 1L
  }

  pieces[[piece_i]] <- out[start_i:nrow(out), , drop = FALSE]
  dplyr::bind_rows(pieces)
}

#' Warn when time-boundary arguments cannot take effect
#'
#' @param dayvar Optional day variable name.
#' @param beepvar Optional beep variable name.
#' @param na_action The requested NA handling mode.
#'
#' @return Invisibly `NULL`.
#' @keywords internal
warn_if_timevars_ineffective <- function(dayvar = NULL, beepvar = NULL, na_action) {
  if ((is.null(dayvar) && is.null(beepvar)) || identical(na_action, "omit_vectors")) {
    return(invisible(NULL))
  }

  time_args <- c()
  if (!is.null(dayvar)) {
    time_args <- c(time_args, "dayvar")
  }
  if (!is.null(beepvar)) {
    time_args <- c(time_args, "beepvar")
  }

  cli::cli_warn(c(
    "{.arg {paste(time_args, collapse = ', ')}} specified with {.code na_action = \"{na_action}\"}.",
    "i" = "Time-boundary separators only prevent bridging vectors when {.code na_action = \"omit_vectors\"}.",
    "x" = "With the current setting, vectors can still connect observations across day boundaries or non-consecutive beeps."
  ))

  invisible(NULL)
}
