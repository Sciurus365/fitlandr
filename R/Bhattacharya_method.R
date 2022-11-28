#' Title
#'
#' @return
#' @export
#'
#' @references Bhattacharya, S., Zhang, Q., & Andersen, M. E. (2011). A deterministic map of Waddington’s epigenetic landscape for cell fate specification. BMC Systems Biology, 5(1), 85. https://doi.org/10.1186/1752-0509-5-85
path_integral_B <- function(f, lims, n = 200, timestep = 1e-2, tol = 1e-2, maxIter = 1400) {
  numPaths <- n^2
  x_path <- matrix(0, nrow = numPaths, ncol = maxIter)
  y_path <- matrix(0, nrow = numPaths, ncol = maxIter)
  pot_path <- matrix(0, nrow = numPaths, ncol = maxIter)

  path_tag <- rep(1L, numPaths)

  path_counter <- 1
  num_attractors <- 0
  num_sepx <- 0

  attractors_num_X_Y <- NULL
  attractors_pot <- NULL
  numPaths_att <- NULL

  # Assign array to keep track of separatrices
  sepx_old_new_pathNum <- NULL

  for (i in seq(lims[1], lims[2], length.out = n)) {
    for (j in seq(lims[3], lims[4], length.out = n)) {
      x0 <- i
      y0 <- j
      p0 <- 0

      x_p <- x0
      y_p <- y0

      Pot <- p0
      Pot_old <- Inf

      x_path[path_counter, 1] <- x_p
      y_path[path_counter, 1] <- y_p
      pot_path[path_counter, 1] <- Pot

      for (n_steps in 2:maxIter) {
        Pot_old <- Pot
        temp_v <- f(c(x_p, y_p))
        dx_dt <- temp_v[1]
        dy_dt <- temp_v[2]
        dx <- dx_dt * timestep
        dy <- dy_dt * timestep

        x_p <- x_p + dx
        y_p <- y_p + dy

        x_path[path_counter, n_steps] <- x_p
        y_path[path_counter, n_steps] <- y_p
        dPot <- -dx_dt * dx - dy_dt * dy
        Pot <- Pot_old + dPot
        pot_path[path_counter, n_steps] <- Pot
      }

      if (path_counter == 1) {
        num_attractors <- num_attractors + 1
        current_att_num_X_Y <- c(num_attractors, x_p, y_p)
        attractors_num_X_Y <- rbind(attractors_num_X_Y, current_att_num_X_Y)
        attractors_pot <- c(attractors_pot, Pot)
        path_tag[path_counter] <- num_attractors
        numPaths_att <- c(numPaths_att, 1)
      } else {
        path_tag[path_counter] <- path_tag[path_counter - 1]

        x0_lastPath <- x_path[path_counter - 1, 1]
        y0_lastPath <- y_path[path_counter - 1, 1]

        xp_lastPath <- x_path[path_counter - 1, maxIter]
        yp_lastPath <- y_path[path_counter - 1, maxIter]

        pot_p_lastPath <- pot_path[path_counter - 1, maxIter]

        startPt_dist_sqr <- (x0 - x0_lastPath)^2 + (y0 - y0_lastPath)^2
        endPt_dist_sqr <- (x_p - xp_lastPath)^2 + (y_p - yp_lastPath)^2

        if (endPt_dist_sqr > 2 * (lims[2] - lims[1]) / n * (lims[4] - lims[3]) / n) {
          new_attr_found <- TRUE
          for (k in 1:num_attractors) {
            x_att <- attractors_num_X_Y[k, 2]
            y_att <- attractors_num_X_Y[k, 3]
            if (abs(x_p - x_att) < (lims[2] - lims[1]) / n & abs(y_p - y_att) < (lims[4] - lims[3]) / n) {
              new_attr_found <- FALSE
              path_tag[path_counter] <- k
              numPaths_att[k] <- numPaths_att[k] + 1
              break
            }
          }

          if (new_attr_found) {
            num_attractors <- num_attractors + 1
            current_att_num_X_Y <- c(num_attractors, x_p, y_p)
            attractors_num_X_Y <- rbind(attractors_num_X_Y, current_att_num_X_Y)
            path_tag[path_counter] <- num_attractors
            numPaths_att <- c(numPaths_att, 1)

            if (startPt_dist_sqr < 2 * (lims[2] - lims[1]) / n * (lims[4] - lims[3]) / n) {
              curr_sepx <- c(
                path_tag[path_counter - 1],
                path_tag[path_counter],
                path_counter - 1
              )
              sepx_old_new_pathNum <- rbind(sepx_old_new_pathNum, curr_sepx)
              attractors_pot <- c(attractors_pot, Pot)
              num_sepx <- num_sepx + 1
            }
          } else {
            prev_attr_new <- TRUE

            for (k in 1:num_sepx) {
              attr1 <- sepx_old_new_pathNum[k, 1]
              attr2 <- sepx_old_new_pathNum[k, 2]

              if (path_tag[path_counter - 1] == attr1 | path_tag[path_counter - 1] == attr2) {
                prev_attr_new <- FALSE
                break
              }
            }

            if (prev_attr_new) {
              if (startPt_dist_sqr < 2 * (lims[2] - lims[1]) / n * (lims[4] - lims[3]) / n) {
                curr_sepx <- c(path_tag[path_counter - 1], path_tag[path_counter], path_counter - 1)
                sepx_old_new_pathNum <- c(sepx_old_new_pathNum, curr_sepx)
                attractors_pot <- c(attractors_pot, pot_p_lastPath)
                num_sepx <- num_sepx + 1
              }
            }
          }
        } else {
          tag <- path_tag[path_counter]
          numPaths_att[tag] <- numPaths_att[tag] + 1
        }
      }
      path_counter <- path_counter + 1
    }
  }

  return(
  	list(
  		attractors_num_X_Y = attractors_num_X_Y,
  		sepx_old_new_pathNum = sepx_old_new_pathNum,
  		numPaths_att = numPaths_att,
  		num_attractors = num_attractors,
  		numPaths = numPaths,
  		maxIter = maxIter,
  		pot_path = pot_path,
  		path_tag = path_tag,
  		attractors_pot = attractors_pot,
  		x_path = x_path,
  		y_path = y_path
  	)
  )
}
