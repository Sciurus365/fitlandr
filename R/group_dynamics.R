normalize_group_data <- function(data, id) {
  if (is.data.frame(data) || is.matrix(data)) {
    if (!is.null(id)) {
      if (!is.data.frame(data)) {
        cli::cli_abort("{.arg id} can only be used when {.arg data} is a data frame.")
      }
      if (length(id) != 1L || !is.character(id) || !id %in% names(data)) {
        cli::cli_abort("{.arg id} must name one column in {.arg data}.")
      }
      if (anyNA(data[[id]])) {
        cli::cli_abort("The grouping column {.field {id}} must not contain missing values.")
      }
      id_values <- unique(data[[id]])
      groups <- lapply(id_values, function(value) data[data[[id]] == value, , drop = FALSE])
      names(groups) <- as.character(id_values)
      return(groups)
    }
    return(list(group_1 = data))
  }

  if (!is.list(data) || !length(data)) {
    cli::cli_abort("{.arg data} must be a non-empty list of datasets, or a data frame with an optional {.arg id} column.")
  }
  if (!is.null(id)) {
    cli::cli_abort("{.arg id} must be {.code NULL} when {.arg data} is already a list.")
  }
  valid <- vapply(data, function(x) is.data.frame(x) || is.matrix(x), logical(1))
  if (!all(valid)) {
    cli::cli_abort("Every element of {.arg data} must be a data frame or matrix.")
  }

  group_names <- names(data)
  if (is.null(group_names)) {
    group_names <- rep.int("", length(data))
  }
  missing_names <- is.na(group_names) | !nzchar(group_names)
  group_names[missing_names] <- paste0("group_", which(missing_names))
  names(data) <- make.unique(group_names)
  data
}

validate_group_variables <- function(groups, x, y) {
  if (length(x) != 1L || length(y) != 1L || !is.character(x) || !is.character(y)) {
    cli::cli_abort("{.arg x} and {.arg y} must each be one column name.")
  }
  missing_variables <- vapply(
    groups,
    function(data) !all(c(x, y) %in% colnames(data)),
    logical(1)
  )
  if (any(missing_variables)) {
    cli::cli_abort(c(
      "Every group dataset must contain columns {.field {x}} and {.field {y}}.",
      "x" = "Missing from: {paste(names(groups)[missing_variables], collapse = ', ')}."
    ))
  }
}

determine_group_lims <- function(groups, x, y, lims) {
  if (!rlang::is_missing(lims)) {
    if (!is.numeric(lims) || length(lims) != 4L || any(!is.finite(lims)) ||
        lims[1L] >= lims[2L] || lims[3L] >= lims[4L]) {
      cli::cli_abort("{.arg lims} must be four finite increasing limits in the order c(x_min, x_max, y_min, y_max).")
    }
    return(as.numeric(lims))
  }

  pooled_x <- unlist(lapply(groups, function(data) data[, x]), use.names = FALSE)
  pooled_y <- unlist(lapply(groups, function(data) data[, y]), use.names = FALSE)
  pooled_x <- pooled_x[is.finite(pooled_x)]
  pooled_y <- pooled_y[is.finite(pooled_y)]
  if (!length(pooled_x) || !length(pooled_y)) {
    cli::cli_abort("The grouping variables must contain finite observations from which common limits can be determined.")
  }

  c(
    grDevices::extendrange(pooled_x, f = 0.1),
    grDevices::extendrange(pooled_y, f = 0.1)
  )
}

validate_stage_args <- function(args, reserved, argument) {
  if (!is.list(args) || (length(args) && (is.null(names(args)) || any(!nzchar(names(args)))))) {
    cli::cli_abort("{.arg {argument}} must be a named list.")
  }
  duplicated <- intersect(names(args), reserved)
  if (length(duplicated)) {
    cli::cli_abort("{.arg {argument}} must not redefine workflow arguments: {paste(duplicated, collapse = ', ')}.")
  }
  args
}

#' Fit vector field, landscape, and stream dynamics for one dataset
#'
#' Runs the recommended two-dimensional fitlandr workflow for one intensive
#' longitudinal dataset: cross-validated vector-field estimation, common-grid
#' interpolation, potential-landscape estimation, probability-flow
#' calculation, and stream-function estimation.
#'
#' @param data A data frame or matrix containing one time series.
#' @param x,y Column names containing the two state variables.
#' @param lims Limits as `c(x_min, x_max, y_min, y_max)`. If omitted, limits
#'   are calculated from the data and extended by 10 percent.
#' @param h_values Candidate bandwidths passed to [cv_fit_2d_vf()].
#' @param cv_folds Number of contiguous cross-validation folds.
#' @param n_grid Number of common grid points per axis for vector-field
#'   interpolation and landscape estimation.
#' @param flow_n Number of probability-flow vectors per axis.
#' @param cv_args Named list of additional arguments passed to
#'   [cv_fit_2d_vf()], such as `method`, `dayvar`, `beepvar`, or `na_action`.
#'   The cross-validation default is `na_action = "omit_vectors"`.
#' @param landscape_args Named list of additional arguments passed to
#'   [make_2d_ld()].
#' @param flow_args Named list of additional arguments passed to [make_2d_pf()].
#'   `divided_by_rho` is fixed to `FALSE` because a stream function is fitted.
#'
#' @return An `individual_dynamics` object containing the input data,
#'   cross-validation result, interpolated vector field, landscape,
#'   probability flow, stream function, and workflow settings.
#'
#' @export
fit_individual_dynamics <- function(
    data,
    x,
    y,
    lims,
    h_values = exp(seq(log(0.01), log(2), length.out = 20)),
    cv_folds = 10L,
    n_grid = 50L,
    flow_n = 20L,
    cv_args = list(),
    landscape_args = list(),
    flow_args = list()) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    cli::cli_abort("{.arg data} must be a data frame or matrix.")
  }
  validate_group_variables(list(individual = data), x, y)
  individual_lims <- determine_group_lims(list(individual = data), x, y, lims)
  if (length(cv_folds) != 1L || !is.finite(cv_folds) ||
      cv_folds != as.integer(cv_folds) || cv_folds < 2L) {
    cli::cli_abort("{.arg cv_folds} must be an integer of at least 2.")
  }
  if (length(n_grid) != 1L || !is.finite(n_grid) ||
      n_grid != as.integer(n_grid) || n_grid < 2L) {
    cli::cli_abort("{.arg n_grid} must be an integer of at least 2.")
  }
  if (length(flow_n) != 1L || !is.finite(flow_n) ||
      flow_n != as.integer(flow_n) || flow_n < 2L) {
    cli::cli_abort("{.arg flow_n} must be an integer of at least 2.")
  }
  cv_args <- validate_stage_args(
    cv_args,
    c("data", "x", "y", "h_values", "k", "lims", "n"),
    "cv_args"
  )
  landscape_args <- validate_stage_args(
    landscape_args,
    c("vf", "n_grid"),
    "landscape_args"
  )
  flow_args <- validate_stage_args(
    flow_args,
    c("vf", "ld", "n", "divided_by_rho"),
    "flow_args"
  )

  cv <- do.call(
    cv_fit_2d_vf,
    c(
      list(
        data = data,
        x = x,
        y = y,
        h_values = h_values,
        k = as.integer(cv_folds),
        lims = individual_lims,
        n = as.integer(n_grid)
      ),
      cv_args
    )
  )
  vectorfield <- add_interp_grid(
    cv$final_model,
    lims = individual_lims,
    n = as.integer(n_grid)
  )
  landscape <- do.call(
    make_2d_ld,
    c(list(vf = vectorfield, n_grid = as.integer(n_grid)), landscape_args)
  )
  probability_flow <- do.call(
    make_2d_pf,
    c(
      list(
        vf = vectorfield,
        ld = landscape,
        n = as.integer(flow_n),
        divided_by_rho = FALSE
      ),
      flow_args
    )
  )
  stream <- make_2d_stream(probability_flow)

  structure(
    list(
      data = data,
      cv = cv,
      vectorfield = vectorfield,
      landscape = landscape,
      probability_flow = probability_flow,
      stream = stream,
      x = x,
      y = y,
      lims = individual_lims,
      settings = list(
        h_values = h_values,
        cv_folds = as.integer(cv_folds),
        n_grid = as.integer(n_grid),
        flow_n = as.integer(flow_n),
        cv_args = cv_args,
        landscape_args = landscape_args,
        flow_args = flow_args
      )
    ),
    class = "individual_dynamics"
  )
}

#' Fit group-level vector fields, landscapes, and streams
#'
#' Runs the complete two-dimensional fitlandr workflow for multiple datasets.
#' Each group receives a cross-validated vector field, potential landscape,
#' probability flow, and stream function. Landscapes are estimated on one
#' pooled range and common grid before candidate K-means solutions are
#' evaluated in steady-state-density space.
#'
#' Supply either a list of data frames/matrices or one data frame together with
#' an `id` column. Use [autoplot()] on the returned object to inspect the elbow
#' plot, then use [add_group_clusters()] to attach a selected solution without
#' repeating the fitted dynamics.
#'
#' @param data A non-empty list of group datasets, or one data frame or matrix.
#' @param x,y Column names containing the two state variables.
#' @param id Optional grouping-column name when `data` is one data frame.
#' @param lims Common limits as `c(x_min, x_max, y_min, y_max)`. If omitted,
#'   limits are calculated from the pooled data and extended by 10 percent.
#' @param h_values Candidate bandwidths passed to [cv_fit_2d_vf()].
#' @param cv_folds Number of contiguous cross-validation folds.
#' @param n_grid Common number of grid points per axis for vector-field
#'   interpolation and landscape estimation.
#' @param flow_n Number of probability-flow vectors per axis.
#' @param k_values Candidate cluster counts. `NULL` respects the limits imposed
#'   by the number of distinct landscapes and the total sample size.
#' @param seed Optional K-means seed.
#' @param cv_args Named list of additional arguments passed to
#'   [cv_fit_2d_vf()], such as `method`, `dayvar`, `beepvar`, or `na_action`.
#'   The cross-validation default is `na_action = "omit_vectors"`.
#' @param landscape_args Named list of additional arguments passed to
#'   [make_2d_ld()].
#' @param flow_args Named list of additional arguments passed to [make_2d_pf()].
#'   `divided_by_rho` is fixed to `FALSE` because a stream function is fitted.
#' @param cluster_method Clustering method passed to
#'   [evaluate_landscape_clusters()]. Currently only `"kmeans"` is available.
#' @param nstart Number of K-means random initializations.
#' @param iter.max Maximum number of K-means iterations.
#'
#' @return A `group_dynamics` object containing member-level fits, convenient
#'   lists of each fitted object, the cluster-count evaluation, and workflow
#'   settings. Its `clustering` field is `NULL` until [add_group_clusters()] is
#'   called.
#'
#' @export
fit_group_dynamics <- function(data,
                               x,
                               y,
                               id = NULL,
                               lims,
                               h_values = exp(seq(log(0.01), log(2), length.out = 20)),
                               cv_folds = 10L,
                               n_grid = 50L,
                               flow_n = 20L,
                               k_values = NULL,
                               seed = NULL,
                               cv_args = list(),
                               landscape_args = list(),
                               flow_args = list(),
                               cluster_method = "kmeans",
                               nstart = 25L,
                               iter.max = 100L) {
  groups <- normalize_group_data(data, id)
  validate_group_variables(groups, x, y)
  common_lims <- determine_group_lims(groups, x, y, lims)

  progress_id <- cli::cli_progress_bar(
    "Fitting group dynamics",
    total = length(groups)
  )
  on.exit(cli::cli_progress_done(progress_id), add = TRUE)
  members <- lapply(seq_along(groups), function(i) {
    member <- fit_individual_dynamics(
      data = groups[[i]],
      x = x,
      y = y,
      lims = common_lims,
      h_values = h_values,
      cv_folds = cv_folds,
      n_grid = n_grid,
      flow_n = flow_n,
      cv_args = cv_args,
      landscape_args = landscape_args,
      flow_args = flow_args
    )
    cli::cli_progress_update(id = progress_id)
    member
  })
  names(members) <- names(groups)

  vectorfields <- lapply(members, `[[`, "vectorfield")
  landscapes <- lapply(members, `[[`, "landscape")
  probability_flows <- lapply(members, `[[`, "probability_flow")
  streams <- lapply(members, `[[`, "stream")
  evaluation <- evaluate_landscape_clusters(
    landscapes,
    k_values = k_values,
    method = cluster_method,
    nstart = nstart,
    iter.max = iter.max,
    seed = seed
  )

  structure(
    list(
      members = members,
      vectorfields = vectorfields,
      landscapes = landscapes,
      probability_flows = probability_flows,
      streams = streams,
      cluster_evaluation = evaluation,
      clustering = NULL,
      x = x,
      y = y,
      lims = common_lims,
      settings = list(
        h_values = h_values,
        cv_folds = cv_folds,
        n_grid = as.integer(n_grid),
        flow_n = as.integer(flow_n),
        cluster_method = cluster_method,
        nstart = as.integer(nstart),
        iter.max = as.integer(iter.max),
        seed = seed
      )
    ),
    class = "group_dynamics"
  )
}

#' Autoplot individual dynamics
#'
#' @param object An `individual_dynamics` object.
#' @param type The fitted component to draw: `"landscape"`, `"stream"`,
#'   `"vectorfield"`, or `"probability_flow"`.
#' @param ... Additional arguments passed to the component's [autoplot()]
#'   method.
#'
#' @return A ggplot object.
#' @export
autoplot.individual_dynamics <- function(
    object,
    type = c("landscape", "stream", "vectorfield", "probability_flow"),
    ...) {
  type <- match.arg(type)
  autoplot(object[[type]], ...)
}

#' Add a selected clustering solution to group dynamics
#'
#' @param object A `group_dynamics` object returned by
#'   [fit_group_dynamics()].
#' @param k Selected number of clusters.
#' @param method,nstart,iter.max,seed Clustering settings. By default, these
#'   reuse the settings from the cluster-count evaluation.
#'
#' @return The `group_dynamics` object with a `landscape_clusters` result in
#'   its `clustering` field.
#' @export
add_group_clusters <- function(object,
                               k,
                               method = object$settings$cluster_method,
                               nstart = object$settings$nstart,
                               iter.max = object$settings$iter.max,
                               seed = object$settings$seed) {
  if (!inherits(object, "group_dynamics")) {
    cli::cli_abort("{.arg object} must be a {.cls group_dynamics} object.")
  }
  object$clustering <- cluster_landscapes(
    object$landscapes,
    k = k,
    method = method,
    nstart = nstart,
    iter.max = iter.max,
    seed = seed
  )
  object
}

#' Autoplot group dynamics
#'
#' @param object A `group_dynamics` object.
#' @param type Either `"elbow"` for the cluster-count evaluation or
#'   `"centers"` for selected cluster-mean landscapes.
#' @param ... Additional arguments passed to the elbow autoplot. Currently
#'   unused for `type = "centers"`.
#'
#' @return A ggplot object.
#' @export
autoplot.group_dynamics <- function(object, type = c("elbow", "centers"), ...) {
  type <- match.arg(type)
  if (type == "elbow") {
    return(autoplot(object$cluster_evaluation, ...))
  }
  if (is.null(object$clustering)) {
    cli::cli_abort(c(
      "No selected clustering solution is attached to {.arg object}.",
      "i" = "Call {.fn add_group_clusters} first."
    ))
  }

  center_data <- do.call(rbind, lapply(
    seq_along(object$clustering$centers),
    function(i) {
      data <- object$clustering$centers[[i]]$dist
      data$U_relative <- data$U - min(data$U, na.rm = TRUE)
      data$cluster <- factor(i, levels = seq_along(object$clustering$centers))
      data
    }
  ))
  ggplot2::ggplot(
    center_data,
    ggplot2::aes(x = .data$x, y = .data$y, fill = .data$U_relative)
  ) +
    ggplot2::geom_raster() +
    ggplot2::facet_wrap(ggplot2::vars(cluster)) +
    ggplot2::scale_fill_viridis_c(name = "Relative U") +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(x = object$x, y = object$y, title = "Cluster-mean landscapes") +
    ggplot2::theme_bw()
}
