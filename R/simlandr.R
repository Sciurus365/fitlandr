#' @importFrom simlandr calculate_barrier
NULL


#' @export
calculate_barrier.2d_ld <- function(l, start_location_value, start_r, end_location_value,
                                    end_r, Umax, expand = TRUE, omit_unstable = FALSE, base = exp(1),
                                    ...) {
  # from simlandr:::calculate_barrier.3d_landscape
  d <- l$dist
  if (missing(Umax)) {
    Umax <- l$Umax
  }
  local_min_start <- find_local_min_3d(d, start_location_value,
    start_r, Umax,
    expand = expand
  )
  local_min_end <- find_local_min_3d(d, end_location_value,
    end_r, Umax,
    expand = expand
  )
  if (is.na(local_min_start$U) | is.na(local_min_end$U)) {
    min_path <- NULL_path()
    saddle_point <- NULL_point()
  } else {
    min_path_index <- dijkstra(
      log(d$d, base = base), local_min_start$location[1:2],
      local_min_end$location[1:2]
    )
    min_path <- min_path_index %>%
      unlist() %>%
      matrix(
        ncol = 2,
        byrow = TRUE
      ) %>%
      as.data.frame() %>%
      {
        colnames(.) <- c("x_index", "y_index")
        .
      } %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        x_value = d$x[x_index],
        y_value = d$y[y_index]
      ) %>%
      dplyr::mutate(U = -log(d$d[
        x_index,
        y_index
      ], base = base)) %>%
      dplyr::ungroup()
    s_U <- max(min_path$U)
    s_location_path_index <- as.integer(stats::median(which(min_path$U ==
      s_U)))
    s_location <- c(s_location_path_index, as.numeric(min_path[
      s_location_path_index[1],
      1:4
    ]))
    names(s_location) <- c(
      "path_index", "x_index", "y_index",
      "x_value", "y_value"
    )
    saddle_point <- list(U = s_U, location = s_location)
    if (omit_unstable & local_min_start$location["x_index"] ==
      saddle_point$location["x_index"] & local_min_start$location["y_index"] ==
      saddle_point$location["y_index"]) {
      local_min_start <- NULL_point()
      saddle_point <- NULL_point()
      min_path <- NULL_path()
    }
    if (omit_unstable & local_min_end$location["x_index"] ==
      saddle_point$location["x_index"] & local_min_end$location["y_index"] ==
      saddle_point$location["y_index"]) {
      local_min_end <- NULL_point()
      saddle_point <- NULL_point()
      min_path <- NULL_path()
    }
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_path(
      data = min_path,
      mapping = ggplot2::aes(x = x_value, y = y_value)
    ) +
    ggplot2::geom_point(ggplot2::aes(
      x = local_min_start$location["x_value"],
      y = local_min_start$location["y_value"]
    ), color = "black") +
    ggplot2::geom_point(ggplot2::aes(
      x = local_min_end$location["x_value"],
      y = local_min_end$location["y_value"]
    ), color = "black") +
    ggplot2::geom_point(ggplot2::aes(
      x = saddle_point$location["x_value"],
      y = saddle_point$location["y_value"]
    ), color = "red") +
    ggplot2::labs(x = l$x, y = l$y)
  result <- list(
    local_min_start = local_min_start, local_min_end = local_min_end,
    saddle_point = saddle_point, min_path = min_path, delta_U_start = saddle_point$U -
      local_min_start$U, delta_U_end = saddle_point$U -
      local_min_end$U, plot = p, geom = list(
      ggplot2::geom_path(
        data = min_path,
        mapping = ggplot2::aes(x = x_value, y = y_value),
        color = "white"
      ), ggplot2::geom_point(ggplot2::aes(
        x = local_min_start$location["x_value"],
        y = local_min_start$location["y_value"]
      ), color = "white"),
      ggplot2::geom_point(ggplot2::aes(
        x = local_min_end$location["x_value"],
        y = local_min_end$location["y_value"]
      ), color = "white"),
      ggplot2::geom_point(ggplot2::aes(
        x = saddle_point$location["x_value"],
        y = saddle_point$location["y_value"]
      ), color = "red")
    ),
    x = l$x, y = l$y, Umax = l$Umax
  )
  class(result) <- c("3d_barrier", "barrier")
  return(result)
}

find_local_min_3d <- function(dist, localmin, r, Umax, expand = TRUE, first_called = TRUE) {
  # from simlandr:::find_local_min_3d
  if (!is.matrix(dist$d)) {
    cli::cli_abort("Wrong input. {.arg dist} should be a list with {.field x}, {.field y}, and {.field d}, and {.field d} should be a matrix.")
  }
  x1 <- localmin[1]
  y1 <- localmin[2]
  if (length(r) == 1) {
    r <- rep(r, 2)
  }
  effective_dist <- dist$d[dist$x > x1 - r[1] & dist$x < x1 +
    r[1], dist$y > y1 - r[2] & dist$y < y1 + r[2]]
  max_dist <- max(effective_dist)
  min_U <- -log(max_dist)
  if (min_U > Umax) {
    if (expand) {
      if (first_called) {
        cli::cli_inform("The U in this range is too high. Searching range expanded...")
      }
      return(find_local_min_3d(dist, localmin, c(r[1] +
        dist$x[2] - dist$x[1], r[2] + dist$y[2] - dist$y[1]),
      Umax,
      first_called = FALSE
      ))
    } else {
      return(NULL_point())
    }
  }
  location_index <- which(dist$d == max_dist, arr.ind = TRUE) %>%
    apply(2, function(x) as.integer(stats::median(x)))
  location_value <- c(dist$x[location_index[1]], dist$y[location_index[2]])
  location <- c(location_index, location_value)
  names(location) <- c("x_index", "y_index", "x_value", "y_value")
  if (!first_called) {
    cli::cli_inform("r = c({round(r[1], 6)}, {round(r[2], 6)})")
  }
  return(list(U = min_U, location = location))
}


# dijkstra <- function (d, s, e)
# {
# 	# from simlandr:::dijkstra
# 	nr <- nrow(d)
# 	nc <- ncol(d)
# 	D <- matrix(Inf, nr, nc)
# 	P <- array(NA, dim = c(nr, nc, 2))
# 	Q <- matrix(Inf, nr, nc)
# 	Q[s[1], s[2]] <- 0
# 	Q_size <- 1
# 	while (Q_size > 0) {
# 		v <- arrayInd(which.min(Q), dim(Q))
# 		v_value <- Q[v[1], v[2]]
# 		Q[v[1], v[2]] <- Inf
# 		Q_size <- Q_size - 1
# 		D[v[1], v[2]] <- v_value
# 		neighs <- get_neighbor_idx(v[1], v[2], nr, nc)
# 		grad_m <- get_gradient_magnitude(v[1], v[2], nr, nc,
# 																		 d)
# 		for (w in neighs) {
# 			vwLength <- D[v[1], v[2]] + sqrt((v[1] - w[1])^2 +
# 																			 	(v[2] - w[2])^2) * grad_m
# 			if (D[w[1], w[2]] != Inf) {
# 				if (vwLength < D[v[1], v[2]])
# 					stop("ValueError: vwLength < D[v[1], v[2]]")
# 			}
# 			else if (vwLength < Q[w[1], w[2]]) {
# 				if (Q[w[1], w[2]] == Inf)
# 					Q_size <- Q_size + 1
# 				Q[w[1], w[2]] <- vwLength
# 				P[w[1], w[2], ] <- c(v[1], v[2])
# 			}
# 		}
# 	}
# 	path <- list()
# 	while (1) {
# 		path <- append(path, list(c(e[1], e[2])))
# 		if (e[1] == s[1] & e[2] == s[2])
# 			break
# 		e <- P[e[1], e[2], ]
# 	}
# 	path <- rev(path)
# 	return(path)
# }
#
# get_neighbor_idx <- function (x, y, nr, nc)
# {
# 	res <- list()
# 	if (x == 1) {
# 		ilist <- c(0, 1)
# 	}
# 	else if (x == nr) {
# 		ilist <- c(-1, 0)
# 	}
# 	else {
# 		ilist <- c(-1, 0, 1)
# 	}
# 	if (y == 1) {
# 		jlist <- c(0, 1)
# 	}
# 	else if (y == nc) {
# 		jlist <- c(-1, 0)
# 	}
# 	else {
# 		jlist <- c(-1, 0, 1)
# 	}
# 	for (i in ilist) {
# 		for (j in jlist) {
# 			if (i == 0 & j == 0)
# 				next
# 			res <- append(res, list(c(x + i, y + j)))
# 		}
# 	}
# 	return(res)
# }
#
# get_gradient_magnitude <- function (x, y, nr, nc, d)
# {
# 	# from simlandr:::get_gradient_magnitude
# 	if (x == 1) {
# 		dx <- d[x + 1, y] - d[x, y]
# 	}
# 	else if (x == nr) {
# 		dx <- d[x, y] - d[x - 1, y]
# 	}
# 	else {
# 		dx <- (d[x + 1, y] - d[x - 1, y]) / 2
# 	}
# 	if (y == 1) {
# 		dy <- d[x, y + 1] - d[x, y]
# 	}
# 	else if (y == nc) {
# 		dy <- d[x, y] - d[x, y - 1]
# 	}
# 	else {
# 		dy <- (d[x, y + 1] - d[x, y - 1]) / 2
# 	}
# 	return(sqrt(dx^2 + dy^2))
# }
