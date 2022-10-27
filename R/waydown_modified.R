### The code in this file is *modified* from the waydown package, with the purpose (1) to calculate the error (or non-gradient) vector field. During the decomposition performed by waydown, the symmetric part is used to construct the potential landscape, but the skew part was also used to calculate the norm to represent the magnitude of the error. Here I also want to see the exact error field, for example, if there are any interesting fluxes in the vector field. To do so, I need the original form of the skew part, not only the norm of it. That is why I need to *modify* some of the code in waydown. And (2) to provide the possibility to use a faster method for numDeriv::jacobian(). (3) to add a tanh regularization.###
### The original code is licensed by the MIT license which is included at the end of this script. ###
### The link to the waydown package: https://CRAN.R-project.org/package=waydown; https://github.com/PabRod/waydown


#' Approximate potential difference between two points
#'
#' @param x Position where we want to know the approximate potential
#' @param x0 Reference position (center of the Taylor expansion)
#' @param f Flow equations (right hand side of differential equation)
#' @param normType (default: 'f') Matrix norm used to compute the error
#' @param jacobianMethod Passed to `method` of [numDeriv::jacobian()].
#' @param ... Other parameters passed to [numDeriv::jacobian()].
#'
#' @return A list containing the approximate potential difference between x and x0 and the estimated error
#' @export
#' @keywords internal
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io}); modified by Jingmeng Cui.
#' @references \url{https://doi.org/10.1371/journal.pcbi.1007788}
#'
#'
#'
#' @examples
#' # Two dimensional flow
#' f <- function(x) {
#'   c(
#'     -2 * x[1] * x[2],
#'     -x[1]^2 - 1
#'   )
#' }
#'
#' # Evaluation points
#' x0 <- matrix(c(1, 2), ncol = 1)
#' x1 <- matrix(c(0.98, 2.01), ncol = 1)
#'
#' dV <- deltaV(f, x1, x0)
deltaV <- function(f, x, x0, normType = "f", jacobianMethod = "Richardson", reg_tanh = FALSE, ...) {
  # Calculate the local Jacobian
  J0 <- numDeriv::jacobian(f, x0, method = jacobianMethod, ...)

  # Perform the skew/symmetric decomposition
  J_symm <- Matrix::symmpart(J0)
  J_skew <- Matrix::skewpart(J0)

  # Use J_symm to estimate the difference in potential as 2nd order Taylor expansion
  #
  # Detailed information available at https://doi.org/10.1371/journal.pcbi.1007788
  dV <- as.numeric(
    -f(x0) %*% (x - x0) + # Linear term
      -0.5 * t(x - x0) %*% J_symm %*% (x - x0) # Quadratic term
  )

  if(reg_tanh != FALSE) {
  	dV <- tanh(reg_tanh * dV)/reg_tanh
  }

  # Use J_skew to estimate the relative error
  #
  # Detailed information available at https://doi.org/10.1371/journal.pcbi.1007788
  rel_err <- norm(J_skew, type = normType) / (norm(J_skew, type = normType) + norm(J_symm, type = normType))

  # Return
  ls <- list(dV = dV, err = rel_err, J_skew_12 = J_skew[1, 2])
  return(ls)
}

#' Approximate potential in two dimensions
#'
#' @param f Two-dimensional representing the flow (right hand side of differential equation).
#' @param xs Vector xs positions to evaluate.
#' @param ys Vector of ys positions to evaluate.
#' @param V0 (Optional) Value of V at first element of (xs,ys). When default, the global minimum is assigned 0.
#' @param mode (Optional) Integration mode. Options are horizontal (default), vertical and mixed.
#' @param ... Other parameters passed to [deltaV()].
#' @inheritParams deltaV
#'
#' @return The potential estimated at each point (xs, ys)
#' @export
#' @keywords internal
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io}); modified by Jingmeng Cui.
#' @references \url{https://doi.org/10.1371/journal.pcbi.1007788}
#'
#'
#'
#' @examples
#' # Flow
#' f <- function(x) {
#'   c(-x[1] * (x[1]^2 - 1.1), -x[2] * (x[2]^2 - 1))
#' }
#'
#' # Sampling points
#' xs <- seq(-1.5, 1.5, length.out = 10)
#' ys <- seq(-1.5, 1.5, length.out = 15)
#'
#' # Approximated potential
#' Vs <- approxPot2D(f, xs, ys, mode = "horizontal")
approxPot2D <- function(f, xs, ys, V0 = "auto", mode = "mixed", jacobianMethod = "Richardson", reg_tanh = FALSE, ...) {
  # Initialize
  V <- matrix(0, nrow = length(xs), ncol = length(ys))
  err <- matrix(0, nrow = length(xs), ncol = length(ys))
  J_skew_12 <- matrix(0, nrow = length(xs), ncol = length(ys))

  # Assign initial value
  # The algorithm is a recursion relationship. It needs an initial potential at the first integration point
  if (V0 == "auto") {
    V[1, 1] <- 0 # Assign any, it will be overriden later
  } else {
    V[1, 1] <- V0 # Assign the desired reference potential
  }

  # Compute
  # We first compute along the first column...
  for (i in 2:length(xs)) {
    temp <- deltaV(f, c(xs[i], ys[1]), c(xs[i - 1], ys[1]), jacobianMethod = jacobianMethod, ...)
    V[i, 1] <- V[i - 1, 1] + temp$dV
    err[i, 1] <- temp$err
    J_skew_12[i - 1, 1] <- temp$J_skew_12
  }

  # ... and then along the first row...
  for (j in 2:length(ys)) {
    temp <- deltaV(f, c(xs[1], ys[j]), c(xs[1], ys[j - 1]), jacobianMethod = jacobianMethod, ...)
    V[1, j] <- V[1, j - 1] + temp$dV
    err[1, j] <- temp$err
    J_skew_12[1, j - 1] <- temp$J_skew_12
  }

  # ... and last but not least, we fill the inside gaps
  for (i in 2:length(xs)) {
    for (j in 2:length(ys)) {
      if (mode == "horizontal") {
        # Sweep horizontally

        temp <- deltaV(f, c(xs[i], ys[j]), c(xs[i - 1], ys[j]), jacobianMethod = jacobianMethod, reg_tanh = reg_tanh, ...)
        V[i, j] <- V[i - 1, j] + temp$dV
        err[i, j] <- temp$err
        J_skew_12[i - 1, j] <- temp$J_skew_12
      } else if (mode == "vertical") {
        # Sweep vertically

        temp <- deltaV(f, c(xs[i], ys[j]), c(xs[i], ys[j - 1]), jacobianMethod = jacobianMethod, reg_tanh = reg_tanh, ...)
        V[i, j] <- V[i, j - 1] + temp$dV
        err[i, j] <- temp$err
        J_skew_12[i, j - 1] <- temp$J_skew_12
      } else if (mode == "mixed") {
        # Sweep in both directions, then take the mean

        temp_hor <- deltaV(f, c(xs[i], ys[j]), c(xs[i - 1], ys[j]), jacobianMethod = jacobianMethod, reg_tanh = reg_tanh, ...)
        V_hor <- V[i - 1, j] + temp_hor$dV
        J_skew_12[i - 1, j] <- temp_hor$J_skew_12
        temp_ver <- deltaV(f, c(xs[i], ys[j]), c(xs[i], ys[j - 1]), jacobianMethod = jacobianMethod, reg_tanh = reg_tanh, ...)
        V_ver <- V[i, j - 1] + temp_ver$dV
        J_skew_12[i, j - 1] <- temp_ver$J_skew_12
        V[i, j] <- mean(c(V_hor, V_ver))
        err[i, j] <- mean(c(temp_hor$err, temp_ver$err))
      } else {
        stop("Error: supported modes are horizontal (default), vertical and mixed")
      }
    }
  }

  if (V0 == "auto") {
    V <- V - min(c(V)) # Make V_min = 0
  }

  for (i in 1:length(xs)) {
    for (j in 1:length(ys)) {
      if (is.na(J_skew_12[i, j])) {
        J_skew_12[i, j] <- numDeriv::jacobian(f, c(xs[i], ys[j]), method = jacobianMethod, ...) |> Matrix::skewpart()
      }
    }
  }

  return(list(xs = xs, ys = ys, V = V, err = err, J_skew_12 = J_skew_12))
}

#' Calculate the nongradient part of the vector field.
#'
#' @param wdresult Output from `approxPot2D()`.
#' @param x_sparse,y_sparse Should the sample points for calculating the non-gradient part of the field be sparser than the sample points for calculating the landscape? Integer numbers larger than 1.
#' @export
#' @keywords internal
#' @author Jingmeng Cui
vectorfield_nongradient_2D <- function(wdresult, x_sparse, y_sparse, ...) {
	if(x_sparse < 1 | y_sparse < 1) {
		warning("`x_sparse` and `y_sparse` should be larger than 1. Using 1 in this calculation.")
	} else if(x_sparse - floor(x_sparse) > 1e-10 | y_sparse - floor(y_sparse)) {
		warning("`x_sparse` and `y_sparse` should be integers. Rounding them up in this calculation.")
	}

	unit_x <- wdresult$xs[2] - wdresult$xs[1]
	unit_y <- wdresult$ys[2] - wdresult$ys[1]

	x_sparse <- max(as.integer(x_sparse), 1L)
	y_sparse <- max(as.integer(y_sparse), 1L)
	xs_index <- seq(1, length(wdresult$xs), x_sparse)
	ys_index <- seq(1, length(wdresult$ys), y_sparse)
	xs <- wdresult$xs[xs_index]
	ys <- wdresult$ys[ys_index]

	ng_vf <- expand.grid(x = xs, y = ys)
	ng_vf$vx <- ng_vf$vy <- NA
	row_index <- 1

	for (j in ys_index) {
		for (i in xs_index) {
			eval_points_h <- list()
			eval_points_v <- list()
			if(i > 1) {
				eval_points_h <- append(eval_points_h, list(c(x_eval_index = i - 1,
																						 y_eval_index = j,
																						 delta_x = unit_x,
																						 delta_y = 0)))
			}
			if(i < length(xs_index)) {
				eval_points_h <- append(eval_points_h, list(c(x_eval_index = i + 1,
																						 y_eval_index = j,
																						 delta_x = -unit_x,
																						 delta_y = 0)))
			}
			if(j > 1) {
				eval_points_v <- append(eval_points_v, list(c(x_eval_index = i,
																						 y_eval_index = j - 1,
																						 delta_x = 0,
																						 delta_y = unit_y)))
			}
			if(j < length(ys_index)) {
				eval_points_v <- append(eval_points_v, list(c(x_eval_index = i,
																						 y_eval_index = j + 1,
																						 delta_x = 0,
																						 delta_y = -unit_y)))
			}
			eval_FUN <- function(x) {
				J_skew_12 <- wdresult$J_skew_12[x["x_eval_index"], x["y_eval_index"]]
				J_skew <- matrix(c(0, J_skew_12, -J_skew_12, 0), ncol = 2, nrow = 2, byrow = TRUE)
				delta_vector <- c(x["delta_x"], x["delta_y"])
				return(as.vector(J_skew %*% delta_vector))
			}
			eval_values_h <- lapply(eval_points_h, FUN = eval_FUN)
			eval_values_v <- lapply(eval_points_v, FUN = eval_FUN)

			if(length(eval_values_h) == 1) eval_values_h <- eval_values_h[[1]]
			else if(length(eval_values_h) == 2) eval_values_h <- rowMeans(cbind(eval_values_h[[1]], eval_values_h[[2]]))
			else stop("Wrong length `eval_values_h`")

			if(length(eval_values_v) == 1) eval_values_v <- eval_values_v[[1]]
			else if(length(eval_values_v) == 2) eval_values_v <- rowMeans(cbind(eval_values_v[[1]], eval_values_v[[2]]))
			else stop("Wrong length `eval_values_v`")

			result <- rowMeans(cbind(eval_values_h, eval_values_v))
			ng_vf[row_index, c("vx", "vy")] <- result
			row_index <- row_index + 1
		}
	}

	return(ng_vf)
}

# The following is the license of the waydown package:
#
# Based on http://opensource.org/licenses/MIT
#
# This is a template. Complete and ship as file LICENSE the following 2
# lines (only)
#
# YEAR:
# 	COPYRIGHT HOLDER:
#
# 	and specify as
#
# License: MIT + file LICENSE
#
# Copyright (c) <YEAR>, <COPYRIGHT HOLDER>
#
# 	Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# 																														"Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# 	The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# YEAR: 2018
# COPYRIGHT HOLDER: Pablo Rodríguez-Sánchez
