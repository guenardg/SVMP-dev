## **************************************************************************
##
##    (c) 2020 Guillaume Guénard
##        Université de Montréal, Montreal, Quebec, Canada
##
##    **Centering functions**
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
#' @name centering
#' 
#' Matrix Centering and Re-centering Functions
#'
#' The function allows one to calculate, in a single step, column means, row
#' means, and the overall mean (\code{center} and \code{get.center}) of a matrix
#' an return the centered matrix (\code{center}), or center a new matrix on
#' values calculated previouly (\code{recenter}).
#' 
#' @param x A data matrix.
#' @param row Whether the centering also on the rows.
#' @param object An object returned by any of the functions listed here.
#' @param newx A new data matrix.
#' 
#' @return A matrix with dimensions equal to \code{x} (\code{center} or
#' \code{get.center}) or \code{newx} (\code{recenter}) and column means
#' (optionally, row means and the overall mean) stored as (an) attribute(s).
#' 
#' @details Function \code{center} calculates the column means (optionally, the
#' row means and overall mean) of a data matrix and returns the latter with its
#' columns centered on means of 0 (optionally, with both its columns and rows
#' centered on means of 0). Function \code{get.center} proceeds similarly as
#' \code{center}, but returns the matrix unchanged. Function \code{recenter}
#' takes column means (optionally, row means and an overall mean) from an object
#' returned by \code{center} or \code{get.center} and centers the columns
#' (optionally the columns and rows) of a new matrix on these values.
#' 
#' The functions are generally intended to be used internally but are made
#' available to the user willing to experiment with new ideas. They are based
#' on C language implementations that can also be used internally in C language
#' computer code.
#' 
#' @author Guillaume Guénard \email{guillaume.guenard@umontreal.ca}
#' 
#' @seealso The \code{\link{dist}-class} and associated methods.
#' 
#' @references
#' [To be included...]
#' 
#' @examples
#' [Examples here...]
#' 
#' @useDynLib SVMP, .registration = TRUE
#' 
NULL

#' @rdname centering
#' @export
center <- function(x, row = FALSE) {
    if(!is.matrix(x))
        stop("'x' must be a matrix!")
    storage.mode(x) <- "double"
    dim <- dim(x)
    if(row) {
        res <- .C("mat_center",x,dim[1L],dim[2L],double(dim[1L]),
                  double(dim[2L]),double(1L),1L)
        out <- res[[1L]]
        attr(out,"rowMeans") <- res[[4L]]
        attr(out,"colMeans") <- res[[5L]]
        attr(out,"mean") <- res[[6L]]
    } else {
        res <- .C("mat_center",x,dim[1L],dim[2L],double(),
                  double(dim[2L]),double(),0L)
        out <- res[[1L]]
        attr(out,"colMeans") <- res[[5L]]
    }
    dim(out) <- dim
    return(out)
}

#' @rdname centering
#' @export
get.center <- function(x, row = FALSE) {
    if(!is.matrix(x))
        stop("'x' must be a matrix!")
    storage.mode(x) <- "double"
    dim <- dim(x)
    if(row) {
        res <- .C("get_center",x,dim[1L],dim[2L],double(dim[1L]),
                  double(dim[2L]),double(1L),1L)
        out <- res[[1L]]
        attr(out,"rowMeans") <- res[[4L]]
        attr(out,"colMeans") <- res[[5L]]
        attr(out,"mean") <- res[[6L]]
    } else {
        res <- .C("get_center",x,dim[1L],dim[2L],double(),
                  double(dim[2L]),double(),0L)
        out <- res[[1L]]
        attr(out,"colMeans") <- res[[5L]]
    }
    dim(out) <- dim
    return(out)
}

#' @rdname centering
#' @export
recenter <- function(object, newx, row) {
    if(!is.matrix(newx))
        stop("'x' must be a matrix!")
    dim <- dim(newx)
    if(length(attr(object,"colMeans"))!=dim[2L])
        stop("'newx' has ",ncol(newx)," columns, but 'object' involves ",
             length(attr(object,"colMeans")),"!")
    if(missing(row)) row <- !is.null(attr(object,"mean"))
    if(row && is.null(attr(object,"mean")))
        stop("No row centering has been performed on 'object' but 'row' is TRUE!")
    storage.mode(newx) <- "double"
    if(row) {
        res <- .C("mat_recenter",newx,dim[1L],dim[2L],double(dim[1L]),
                  attr(object,"colMeans"),attr(object,"mean"),1L)
        out <- res[[1L]]
        attr(out,"rowMeans") <- res[[4L]]
        attr(out,"colMeans") <- res[[5L]]
        attr(out,"mean") <- res[[6L]]
    } else {
        res <- .C("mat_recenter",newx,dim[1L],dim[2L],double(),
                  attr(object,"colMeans"),double(),0L)
        out <- res[[1L]]
        attr(out,"colMeans") <- res[[5L]]
    }
    dim(out) <- dim
    return(out)
}
##
