#' Bhattacharya method for path integration
#'
#' @param f The vector field function. It should return `c(<dx/dt>, <dy/dt>)`
#' @param lims The limits of the range for the estimation as `c(<xl>, <xu>, <yl>, <yu>)`.
#' @param n_path_int The number of equally spaced points in each axis, at which the path integrals is to be calculated.
#' @param stepsize The time step used in each iteration.
#' @param tol The tolerance to test convergence.
#' @param numTimeSteps Number of time steps for integrating along each path (to ensure uniform arrays). Choose high-enough number for convergence with given stepsize.
#' @param ... Not in use.
#'
#'
#' @return A list with the following elements:
#' - `numPaths` Integer. Total Number of paths for defined grid spacing.
#' - `pot_path` Matrix. Potential along the paths.
#' - `path_tag` Vector. Tag for given paths.
#' - `attractors_pot` Vector. Potential value of each identified attractors by the path integral approach.
#' - `x_path` Vector. x-coord. along path.
#' - `y_path` Vector. y-coord. along path.
#'
#' @export
#' @keywords internal
#' @references Bhattacharya, S., Zhang, Q., & Andersen, M. E. (2011). A deterministic map of Waddington’s epigenetic landscape for cell fate specification. BMC Systems Biology, 5(1), 85. https://doi.org/10.1186/1752-0509-5-85.
#' The functions in this file were translated from the Matlab code provided with the reference above, and its Python translation at https://dynamo-release.readthedocs.io/en/v0.95.2/_modules/dynamo/vectorfield/Bhattacharya.html
path_integral_B <- function(f, lims, n_path_int = 20, stepsize = 1e-2, tol = 1e-2, numTimeSteps = 1400, ...) {
  numPaths <- n_path_int^2

  # Initialize "path" variable matrices
  x_path <- matrix(0, nrow = numPaths, ncol = numTimeSteps)
  y_path <- matrix(0, nrow = numPaths, ncol = numTimeSteps)
  pot_path <- matrix(0, nrow = numPaths, ncol = numTimeSteps)

  path_tag <- rep(1L, numPaths) # tag for given path (to denote basin of attraction). initialized to 1 for all paths.

	# Initialize "Path counter" to 1
  path_counter <- 1

  # Initialize no. of attractors and separatrices (basin boundaries)
  num_attractors <- 0
  num_sepx <- 0

  # Assign array to keep track of attractors and their coordinates; and pot.
  attractors_num_X_Y <- NULL
  attractors_pot <- NULL

  # Assign array to keep track of no. of paths per attractor
  numPaths_att <- NULL

  # Assign array to keep track of separatrices
  sepx_old_new_pathNum <- NULL

  ## Set up a progress bar
  cli::cli_progress_bar("Looping over the x-y grid ...", total = n_path_int^2)

  # Loop over x-y grid
  for (i in seq(lims[1], lims[2], length.out = n_path_int)) {
    for (j in seq(lims[3], lims[4], length.out = n_path_int)) {
    	# Init conds for given (x,y)
    	# Initialize coords.
      x0 <- i
      y0 <- j

      # Set initial value of "potential" to 0 (to facilitate comparison of "potential drop")
      p0 <- 0

      # Initialize "path" variables
      x_p <- x0
      y_p <- y0

      # Initialize accumulators for "potential" along path
      Pot <- p0
      Pot_old <- Inf

      # Initialize global arrays (time t = 0 counts as "time step #1")
      x_path[path_counter, 1] <- x_p
      y_path[path_counter, 1] <- y_p
      pot_path[path_counter, 1] <- Pot

      # Evaluate potential (Integrate) over trajectory from init cond to stable steady state
      for (n_steps in 2:numTimeSteps) {
      	# record "old" values of variables
        Pot_old <- Pot

        # update dx/dt, dy/dt
        temp_v <- f(c(x_p, y_p))
        dx_dt <- temp_v[1]
        dy_dt <- temp_v[2]
        dx <- dx_dt * stepsize
        dy <- dy_dt * stepsize

        # update x, y
        x_p <- x_p + dx
        y_p <- y_p + dy

        x_path[path_counter, n_steps] <- x_p
        y_path[path_counter, n_steps] <- y_p

        # update "potential"
        dPot <- -dx_dt * dx - dy_dt * dy # signs ensure that "potential" decreases as "velocity" increases
        Pot <- Pot_old + dPot
        pot_path[path_counter, n_steps] <- Pot

        # end integration over path
      }

      # check for convergence
      if (abs(Pot - Pot_old) > tol) {
      	message(glue::glue("Warning: Not converged.
      										 \t Start point: {i}, {j}.
      										 \t End point: {x_p}, {y_p}. \n"))
      }

      # assign path tag (to track multiple basins of attraction)
      if (path_counter == 1) {
      	# record attractor of first path and its coords
        num_attractors <- num_attractors + 1
        current_att_num_X_Y <- c(num_attractors, x_p, y_p)
        attractors_num_X_Y <- rbind(attractors_num_X_Y, current_att_num_X_Y)
        attractors_pot <- c(attractors_pot, Pot)
        path_tag[path_counter] <- num_attractors
        numPaths_att <- c(numPaths_att, 1)
      } else {
      	# i.e., if path counter > 1
      	# set path tag to that of previous path (default)
        path_tag[path_counter] <- path_tag[path_counter - 1]

				# record info of previous path
        x0_lastPath <- x_path[path_counter - 1, 1]
        y0_lastPath <- y_path[path_counter - 1, 1]

        xp_lastPath <- x_path[path_counter - 1, numTimeSteps]
        yp_lastPath <- y_path[path_counter - 1, numTimeSteps]

        pot_p_lastPath <- pot_path[path_counter - 1, numTimeSteps]

        # calculate distance between "start points" of current and previous paths
        startPt_dist_sqr <- (x0 - x0_lastPath)^2 + (y0 - y0_lastPath)^2
        # calculate distance between "end points" of current and previous
        endPt_dist_sqr <- (x_p - xp_lastPath)^2 + (y_p - yp_lastPath)^2

        # check if the current path *ended* in a different point compared to previous path
        # (x-y grid spacing used as a "tolerance" for distance)
        if (endPt_dist_sqr > 2 * (lims[2] - lims[1]) / n_path_int * (lims[4] - lims[3]) / n_path_int) {
        	# check if this "different" attractor has been identified before
          new_attr_found <- TRUE
          for (k in 1:num_attractors) {
            x_att <- attractors_num_X_Y[k, 2]
            y_att <- attractors_num_X_Y[k, 3]
            if (abs(x_p - x_att) < (lims[2] - lims[1]) / n_path_int & abs(y_p - y_att) < (lims[4] - lims[3]) / n_path_int) {
            	# this attractor has been identified before
              new_attr_found <- FALSE
              path_tag[path_counter] <- k
              numPaths_att[k] <- numPaths_att[k] + 1
              break # exit for-loop
            }
          }

          if (new_attr_found) {
            num_attractors <- num_attractors + 1
            current_att_num_X_Y <- c(num_attractors, x_p, y_p)
            attractors_num_X_Y <- rbind(attractors_num_X_Y, current_att_num_X_Y)
            path_tag[path_counter] <- num_attractors
            numPaths_att <- c(numPaths_att, 1)

            # check if start points of current and previous paths are "adjacent"
            # if so, assign separatrix
            if (startPt_dist_sqr < 2 * (lims[2] - lims[1]) / n_path_int * (lims[4] - lims[3]) / n_path_int) {
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
          	# check if the attractor of the *previous* path as been encountered in a separatrix before
          	# (note that current path tag has already been set above)
            prev_attr_new <- TRUE

            for (k in 1:num_sepx) {
              attr1 <- sepx_old_new_pathNum[k, 1]
              attr2 <- sepx_old_new_pathNum[k, 2]

              if (path_tag[path_counter - 1] == attr1 | path_tag[path_counter - 1] == attr2) {
              	# this attractor has been identified before
                prev_attr_new <- FALSE
                break
              }
            }

            if (prev_attr_new) {
            	# check if start points of current and previous paths are "adjacent"
            	# if so, assign separatrix
              if (startPt_dist_sqr < 2 * (lims[2] - lims[1]) / n_path_int * (lims[4] - lims[3]) / n_path_int) {
                curr_sepx <- c(path_tag[path_counter - 1], path_tag[path_counter], path_counter - 1)
                sepx_old_new_pathNum <- rbind(sepx_old_new_pathNum, curr_sepx)
                attractors_pot <- c(attractors_pot, pot_p_lastPath)
                num_sepx <- num_sepx + 1
              }
            }
          }
        } else {
        	# i.e. current path converged at same pt. as previous path
          tag <- path_tag[path_counter]
          numPaths_att[tag] <- numPaths_att[tag] + 1
        }
      }
      path_counter <- path_counter + 1

      cli::cli_progress_update()
    }
  }
  cli::cli_process_done()

  return(
    list(
      attractors_num_X_Y = attractors_num_X_Y,
      sepx_old_new_pathNum = sepx_old_new_pathNum,
      numPaths_att = numPaths_att,
      num_attractors = num_attractors,
      numPaths = numPaths,
      numTimeSteps = numTimeSteps,
      pot_path = pot_path,
      path_tag = path_tag,
      attractors_pot = attractors_pot,
      x_path = x_path,
      y_path = y_path
    )
  )
}


#' Align potential values
#'
#' So all path-potentials end up at same global min and then generate potential surface with interpolation on a grid.
#'
#' @param resultB Result from [path_integral_B()].
#' @param n The number of equally spaced points in each axis, at which the landscape is to be estimated.
#' @param digits Currently, the raw sample points in some regions are too dense that may crashes interpolation. To avoid this problem, only one point of all with the  same first several digits. is kept. Use this parameter to indicate how many digits are considered. Note that this is a temporary solution and might be changed in the near future.
#' @inheritParams akima::interp
#' @param ... Other parameters passed to [akima::interp()]
#'
#' @keywords internal
#' @export
#'
#' @inherit akima::interp return
align_pot_B <- function(resultB,
											n = 200,
											digits = 2,
											linear = TRUE,
											...) {
	list2env(resultB, envir=environment())

	list_size <- numPaths * numTimeSteps
	x_p_list <- rep(0, list_size)
	y_p_list <- rep(0, list_size)
	pot_p_list <- rep(0, list_size)

	n_list <- 1

	cli::cli_progress_bar("Align potential values", total = numPaths)
	# "Align" potential values so all path-potentials end up at same global min.
	for(n_path in 1:numPaths) {
		tag <- path_tag[n_path]
		del_pot <- pot_path[n_path, numTimeSteps] - attractors_pot[tag]

		# align pot. at each time step along path
		for(n_steps in 1:numTimeSteps) {
			pot_old <- pot_path[n_path, n_steps]
			pot_path[n_path, n_steps] <- pot_old - del_pot

			# add data point to list
			x_p_list[n_list] <- x_path[n_path, n_steps]
			y_p_list[n_list] <- y_path[n_path, n_steps]
			pot_p_list[n_list] <- pot_path[n_path, n_steps]

			n_list <- n_list + 1
		}
		cli::cli_progress_update()
	}
	cli::cli_process_done()

	# Generate surface interpolation grid
	xlin <- seq(min(x_p_list), max(x_p_list), length.out = n)
	ylin <- seq(min(y_p_list), max(y_p_list), length.out = n)

	df <- data.frame(x = x_p_list, y = y_p_list, z = pot_p_list)
	df_sparse <- df[!duplicated(round(df[,c("x","y")], digits = digits)),]

	if(any(!is.finite(df_sparse$z))) {
		warning("`z` contains non-finite values. Removed automatically, but be careful with the result!")
		df_sparse <- df_sparse %>%
			filter(is.finite(z))
	}

	return(R.utils::doCall(akima::interp, x = df_sparse$x, y = df_sparse$y, z = df_sparse$z, xo = xlin, yo = ylin, linear = linear, args = list(...)))
}

#' Options controlling the path-integral algorithm
#' See [path_integral_B()], [align_pot_B()] for details.
#' @inheritParams sim_vf
#' @inheritParams path_integral_B
#' @inheritParams align_pot_B
#' @export
pathB_options <- function(vf, lims = rlang::expr(vf$lims), n_path_int = 20, stepsize = 1e-2, tol = 1e-2, numTimeSteps = 1400, n = 200, digits = 2, linear = TRUE, ...) {
	if(!missing(vf)) return(list(vf = vf, lims = eval(lims), n_path_int = n_path_int, stepsize = stepsize, tol = tol, numTimeSteps = numTimeSteps, n = n, digits = digits, linear = linear, ...))
	else return(list(vf = rlang::expr(vf), lims = lims, n_path_int = n_path_int, stepsize = stepsize, tol = tol, numTimeSteps = numTimeSteps, n = n, digits = digits, linear = linear, ...))
}
