#' @importFrom simlandr calculate_barrier calculate_barrier_3d_core
NULL


#' @export
calculate_barrier.2d_ld <- function(l, start_location_value, start_r, end_location_value,
                                    end_r, Umax, expand = TRUE, omit_unstable = FALSE, base = exp(1),
                                    ...) {
  if (missing(Umax)) {
    Umax <- l$Umax
  }

  simlandr::calculate_barrier_3d_core(
    d = l$dist,
    x_label = l$x,
    y_label = l$y,
    start_location_value = start_location_value,
    start_r = start_r,
    end_location_value = end_location_value,
    end_r = end_r,
    Umax = Umax,
    expand = expand,
    omit_unstable = omit_unstable,
    base = base
  )
}
