### The code in this file is *modified* from the waydown package, with the purpose to calculate the error (or non-gradient) vector field. During the decomposition performed by waydown, the symmetric part is used to construct the potential landscape, but the skew part was also used to calculate the norm to represent the magnitude of the error. Here I also want to see the exact error field, for example, if there are any interesting fluxes in the vector field. To do so, I need the original form of the skew part, not only the norm of it. That is why I need to *modify* some of the code in waydown. ###
### Another purpose is to provide the possibility to use a faster method for numDeriv::jacobian(). ###
### The original code is licensed by the MIT license which is included at the end of this script. ###
### The link to the waydown package: https://CRAN.R-project.org/package=waydown; https://github.com/PabRod/waydown


#' Approximate potential difference between two points
#'
#' @param x Position where we want to know the approximate potential
#' @param x0 Reference position (center of the Taylor expansion)
#' @param f Flow equations (right hand side of differential equation)
#' @param normType (default: 'f') Matrix norm used to compute the error
#'
#' @return A list containing the approximate potential difference between x and x0 and the estimated error
#' @export
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#' @references \url{https://doi.org/10.1371/journal.pcbi.1007788}
#'
#'
#' @seealso \code{\link{approxPot1D}, \link{approxPot2D}, \link{norm}}
#'
#' @examples
#' # One dimensional flow
#' f <- function(x) { cos(x) }
#'
#' # Evaluation points
#' x0 <- 1
#' x1 <- 1.02
#'
#' dV <- deltaV(f, x1, x0)
#'
#'  # Two dimensional flow
#' f <- function(x) { c(
#'  -2*x[1]*x[2],
#'  -x[1]^2 - 1
#' )}
#'
#' # Evaluation points
#' x0 <- matrix(c(1,2), ncol = 1)
#' x1 <- matrix(c(0.98,2.01), ncol = 1)
#'
#' dV <- deltaV(f, x1, x0)
deltaV <- function(f, x, x0, normType='f') {

  # Calculate the local Jacobian
  J0 <- numDeriv::jacobian(f, x0)

  # Perform the skew/symmetric decomposition
  J_symm <- Matrix::symmpart(J0)
  J_skew <- Matrix::skewpart(J0)

  # Use J_symm to estimate the difference in potential as 2nd order Taylor expansion
  #
  # Detailed information available at https://doi.org/10.1371/journal.pcbi.1007788
  dV <- as.numeric(
        -f(x0) %*% (x - x0) +  # Linear term
        -0.5 * t(x-x0) %*% J_symm %*% (x - x0) # Quadratic term
  )

  # Use J_skew to estimate the relative error
  #
  # Detailed information available at https://doi.org/10.1371/journal.pcbi.1007788
  rel_err <- norm(J_skew, type = normType)/(norm(J_skew, type = normType) + norm(J_symm, type = normType))

  # Return
  ls <- list(dV = dV, err = rel_err)
  return(ls)
}


#' Approximate potential in two dimensions
#'
#' @param f Two-dimensional representing the flow (right hand side of differential equation)
#' @param xs Vector xs positions to evaluate
#' @param ys Vector of ys positions to evaluate
#' @param V0 (Optional) Value of V at first element of (xs,ys). When default, the global minimum is assigned 0
#' @param mode (Optional) Integration mode. Options are horizontal (default), vertical and mixed
#'
#' @return The potential estimated at each point (xs, ys)
#' @export
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#' @references \url{https://doi.org/10.1371/journal.pcbi.1007788}
#'
#'
#' @seealso \code{\link{approxPot1D}, \link{deltaV}}
#'
#' @examples
#' # Flow
#' f = function(x) {c(-x[1]*(x[1]^2 - 1.1), -x[2]*(x[2]^2 - 1))}
#'
#' # Sampling points
#' xs <- seq(-1.5, 1.5, length.out = 10)
#' ys <- seq(-1.5, 1.5, length.out = 15)
#'
#' # Approximated potential
#' Vs <- approxPot2D(f, xs, ys, mode = 'horizontal')
approxPot2D <- function(f, xs, ys, V0 = 'auto', mode = 'mixed') {
  # Initialize
  V <- matrix(0, nrow = length(xs), ncol = length(ys))
  err <- matrix(0, nrow = length(xs), ncol = length(ys))

  # Assign initial value
  # The algorithm is a recursion relationship. It needs an initial potential at the first integration point
  if (V0 == 'auto') {
    V[1,1] <- 0 # Assign any, it will be overriden later
  } else {
    V[1,1] <- V0 # Assign the desired reference potential
  }

  # Compute
  # We first compute along the first column...
  for(i in 2:length(xs)) {
    temp <- deltaV(f, c(xs[i], ys[1]), c(xs[i-1], ys[1]))
    V[i,1] <- V[i-1,1] + temp$dV
    err[i,1] <- temp$err
  }

  # ... and then along the first row...
  for(j in 2:length(ys)) {
    temp <- deltaV(f, c(xs[1], ys[j]), c(xs[1], ys[j-1]))
    V[1,j] <- V[1,j-1] + temp$dV
    err[1,j] <- temp$err
  }

  # ... and last but not least, we fill the inside gaps
  for(i in 2:length(xs)) {
    for(j in 2:length(ys)) {

      if(mode == 'horizontal') { # Sweep horizontally

        temp <- deltaV(f, c(xs[i], ys[j]), c(xs[i-1], ys[j]))
        V[i,j] <- V[i-1,j] + temp$dV
        err[i,j] <- temp$err

      } else if(mode == 'vertical') { # Sweep vertically

        temp <- deltaV(f, c(xs[i], ys[j]), c(xs[i], ys[j-1]))
        V[i,j] <- V[i,j-1] + temp$dV
        err[i,j] <- temp$err

      } else if(mode == 'mixed') { # Sweep in both directions, then take the mean

        temp_hor <- deltaV(f, c(xs[i], ys[j]), c(xs[i-1], ys[j]))
        V_hor <- V[i-1,j] + temp_hor$dV
        temp_ver <- deltaV(f, c(xs[i], ys[j]), c(xs[i], ys[j-1]))
        V_ver <- V[i,j-1] + temp_ver$dV
        V[i,j] <- mean(c(V_hor, V_ver))
        err[i,j] <- mean(c(temp_hor$err, temp_ver$err))

      } else {

        stop('Error: supported modes are horizontal (default), vertical and mixed')

      }
    }
  }

  if(V0 == 'auto') {
    V <- V - min(c(V)) # Make V_min = 0
  }

  return(list(V = V, err = err))
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
