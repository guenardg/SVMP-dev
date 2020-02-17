## **************************************************************************
##
##    (c) 2020 Guillaume Guénard
##        Université de Montréal, Montreal, Quebec, Canada
##
##    **Spatial covariance function wrapper**
##
##    This file is part of SVMP
##
##    SVMP is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    SVMP is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with SVMP.  If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##
#' Calculation of Spatial Covariance Models
#' 
#' Function \code{sp.cov} carries out the calculation of the distance-based
#' spatial covariance models.
#' 
#' @param d Pairwise distances among objects (see details).
#' @param type The type of covariance model (see details).
#' @param alpha Shape parameter of the covariance model.
#' @param beta Range parameter of the covariance model.
#' 
#' @return An array with the same dimensions as argument \code{d}.
#' 
#' @details Pairwise distances among objects can be an array with arbitrary
#' dimensions or a \code{\link{dist}-class} object, in which case a
#' one-dimensional array will be returned.
#' 
#' Argument \code{type} can be one of "spherical" (the default), "exponential",
#' "power", "hyperbolic", "superelliptic", or any unambiguous abbreviation; and
#' correspond to the particular function used to obtain the weights as a
#' function of the distances. The function does not directly return a
#' covariance matrix but weights that are used calculate the spatial
#' covariances.
#' 
#' Distance-weight relationships are monotonically decreasing functions
#' associating the level of ressemblance of two locations with the distance
#' separating them. For all covariance models (argument \code{type} above), the
#' weight is 1 at \code{d=0} and 0 beyond the range (i.e., when
#' \code{d >= beta}). The shape of the relationship is controled by the shape
#' parameter given as argument \code{alpha}. That parameter controls the
#' gentleness of the slope of the distance-weight relationship. It is bounded
#' between values of 0 (the steepest slope: a Kronecker's delta). and 1 (a
#' straight line).
#' 
#' @author Guillaume Guénard \email{guillaume.guenard@umontreal.ca}
#' 
#' @seealso The \code{\link{dist}-class} and associated methods.
#' 
#' @references
#' [To be included...]
#' 
#' @examples
#' ##
#' ### First example: locations spread throughout the world
#' ##
#' ## [Examples here...]
#' 
#' @useDynLib SVMP, .registration = TRUE
#' 
#' @export
sp.cov <- function(d,
                   type = c("spherical","exponential","power",
                            "hyperbolic","superelliptic"),
                   alpha = 1, beta = 1) {
  storage.mode(d) <- "double"
  type <- sprintf("scf_%s",substr(match.arg(type),1L,5L))
  storage.mode(alpha) <- "double"
  storage.mode(beta) <- "double"
  n <- length(d)
  if(any(c(length(alpha),length(beta)) > 1L)) {
    alpha <- rep(alpha, length.out = n)
    beta <- rep(beta, length.out = n)
    out <- .C(type,d,alpha,beta,n,0L,double(n),NAOK=TRUE)[[6L]]
  } else {
    out <- .C(type,d,alpha,beta,n,1L,double(n),NAOK=TRUE)[[6L]]
  }
  dim(out) <- dim(d)
  return(out)
}
