# borrowed from simlandr 0.3.0
determine_lims <- function(output, var_names, lims) {
  if (!rlang::is_missing(lims)) {
    return(lims)
  }
  # if (is.list(output)) output <- output[[1]]
  if (rlang::is_missing(lims)) {
    return(c(sapply(var_names, function(v) grDevices::extendrange(output[, v], f = 0.1))))
  }
  if (any(is.infinite(lims))) cli::cli_abort("Non-infinite values found in {.arg lims}.")
}
