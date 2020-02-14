##
##    (c) 2019 Guillaume Guénard
##        Université de Montréal, Montreal, Quebec, Canada
##
##    This file is part of spectR
##
##    spectR is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    spectR is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
sp.cov <- function(d,
                   type = c("spherical","exponential","power",
                            "hyperbolic","superelliptic"),
                   alpha = 1, beta = 1) {
    storage.mode(d) <- "double"
    type <- substr(match.arg(type),1L,5L)
    storage.mode(alpha) <- "double"
    storage.mode(beta) <- "double"
    n <- length(d)
    if(any(c(length(alpha),length(beta)) > 1L)) {
        alpha <- rep(alpha, length.out = n)
        beta <- rep(beta, length.out = n)
        out <- .C(sprintf("scf_%s",type),d,alpha,beta,
                  n,0L,double(n),NAOK=TRUE)[[6L]]
    } else {
        out <- .C(sprintf("scf_%s",type),d,alpha,beta,
                  n,1L,double(n),NAOK=TRUE)[[6L]]
    }
    dim(out) <- dim(d)
    return(out)
}
##
### I was thinking about a new implementation with internal recycling of the
### parameters but I postponed that implementation to work on more pressing
### issues.
##
## sp.cov2 <- function(d,
##                     type = c("spherical","exponential","power",
##                              "hyperbolic","superelliptic"),
##                     alpha = 1, beta = 1) {
##     storage.mode(d) <- "double"
##     type <- substr(match.arg(type),1L,5L)
##     storage.mode(alpha) <- "double"
##     storage.mode(beta) <- "double"
##     n <- c(length(d),length(alpha),length(beta))
##     out <- .C(sprintf("scf2_%s",type),d,alpha,beta,
##               n,any(n[2L:3L]>1L),double(n),NAOK=TRUE)[[6L]]
##     dim(out) <- dim(d)
##     return(out)
## }
##
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
##
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
    return(out)
}
##
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
    return(out)
}
##
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
    return(out)
}
##
