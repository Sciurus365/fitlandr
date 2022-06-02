#' @export
plot.vectorfield <- function(x, arrow = grid::arrow(length = grid::unit(0.1, "cm")),
														 estimated_vector_enlarge = 1,
														 estimated_vector_options = list(),
														 show_point = TRUE,
														 point_options = list(size = 0.5),
														 show_original_vector = FALSE,
                             original_vector_enlarge = 1,
														 original_vector_options = list(),
														 show_used_vector = FALSE,
														 used_vector_options = list(color = "red"),
														 ...) {
  p <-
    ggplot2::ggplot(x$vec_grid, ggplot2::aes(x = x, y = y)) +
    do.call(
      ggplot2::geom_segment,
      c(
        list(
          mapping = ggplot2::aes(
          	xend = x + vx * estimated_vector_enlarge,
          	yend = y + vy * estimated_vector_enlarge),
          arrow = arrow
        ),
        estimated_vector_options
      )
    ) +
    ggplot2::theme_bw() +
  	ggplot2::labs(x = x$x, y = x$y)

  if (show_point) {
    p <- p +
      do.call(
        ggplot2::geom_point,
        c(
          list(
            data = x$data %>% as.data.frame(),
            mapping = ggplot2::aes(x = .data[[x$x]], y = .data[[x$y]])
          ),
          point_options
        )
      )
  }

  if (show_original_vector) {
  	p <- p +
  		do.call(
  			ggplot2::geom_segment,
  			c(
  				list(
  					data = x$original_vectors %>% as.data.frame(),
  					mapping = ggplot2::aes(
  						x = x,
  						y = y,
  						xend = x + vx * original_vector_enlarge,
  						yend = y + vy * original_vector_enlarge),
  					arrow = arrow
  				),
  				original_vector_options
  			)
  		)
  }

  if (show_used_vector) {
  	p <- p +
  		do.call(
  			ggplot2::geom_segment,
  			c(
  				list(
  					data = x$original_vectors[x$VFCresult$VFCIndex,] %>% as.data.frame(),
  					mapping = ggplot2::aes(
  						x = x,
  						y = y,
  						xend = x + vx * original_vector_enlarge,
  						yend = y + vy * original_vector_enlarge),
  					arrow = arrow
  				),
  				used_vector_options
  			)
  		)
  }

  return(p)
}

plot.2d_vf_landscape <- function(x, ...){

}
