## **************************************************************************
##
##    (c) 2020 Guillaume Guénard
##        Université de Montréal, Montreal, Quebec, Canada
##
##    **Euclidean distance between two point sets A and B**
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
#' Calculation of the Euclidean Distance Among Two Sets of Points
#' 
#' Function \code{EuclidAB} calculates the Euclidean distance between two sets
#' of points A and B.
#' 
#' @param a First matrix of cartesian coordinates.
#' @param b Second matrix of cartesian coordinates.
#' @param squared Should the squared Euclidean distances be returned (default:
#' FALSE).
#' 
#' @return An matrix with as many rows as \code{a} and as many columns as there
#' are rows in \code{b}.
#' 
#' @details The standard \code{R} function used to calculate the Euclidean
#' distance, namely \code{\link{dist}}, does only allow one to calculate
#' pairwise distances between the rows of a single matrix of cartesian
#' coordinates and return a \code{\link{dist}}-class object, which is a
#' one-dimensional array that is meant to be interpreted as a triangular matrix.
#' Function \code{Euclid} allows the user to provide two data matrices \code{a}
#' and \code{b} and output a rectangular Euclidean distance matrix.
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
#' ### First example
#' ##
#' ## [Examples here...]
#' 
#' @useDynLib SVMP, .registration = TRUE
#' 
#' @export
EuclidAB <- function(a, b, squared=FALSE) {
  m <- NCOL(a)
  if(m!=NCOL(b))
    stop("'a' and 'b' must have the equal numbers of columns!")
  storage.mode(a) <- "double"
  storage.mode(b) <- "double"
  na <- NROW(a)
  nb <- NROW(b)
  res <- .C("dist_Euclid",a,b,na,nb,m,double(na*nb),
            as.integer(squared))[[6L]]
  dim(res) <- c(nb,na)
  dimnames(res) <- list(rownames(b),rownames(a))
  return(res)
}
