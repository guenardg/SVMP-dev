## **************************************************************************
##
##    (c) 2020 Guillaume Guénard
##        Université de Montréal, Montreal, Quebec, Canada
##
##    **Moran's coefficient**
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
#' Calculation and Test of Moran's Coefficent
#' 
#' Carries out the calculation of Moran's autocorrelation coefficient and test
#' its significance using Monte-Carlo permutations.
#' 
#' @param x A vector of one-dimensional array.
#' @param w A matrix of connection weights.
#' @param test The type of permutation test (default: "two-tail", see details)
#' @param nperm The number of random permutations of \code{x} (default: 999).
#' @param saveIp Whether or not to save the permuted coefficient values
#' (default: FALSE).
#' 
#' @return A list with Moran's I values and testing information.
#' 
#' @details The Moran's autocorrelation coefficient allows one to test the bulk
#' autocorrelation in a data set. It can be positive, in which case observations
#' located closer tend to be more similar than random observations drawn
#' irrespective of the location, or negative, in which case observations located
#' closer tend to be more dissimilar (or less similar) than random observations
#' drawn irrespective of the location. In the (somewhat unlikely) situation
#' where signal components associated with both positive and negative
#' autocorrelation would be present in a data set, bulk autocorrelation as
#' quantified by Moran's coefficient may appear small and non-significant in
#' spite of the presence of a structured signal.
#' 
#' The matrice of connection weights can be a binary matrix of a weighted
#' matrix obtained from function \code{\link{sp.cov}}.
#' 
#' Options for the permutation test (argument \code{test}) are either "two-tail"
#' for a two-tail test, "lower-tail" to perform a single-tail test under the
#' alternate hypothesis that the coefficient is smaller (more negative) than
#' the value expected for the null hypothesis, "upper-tail" to perform a
#' single-tail test under the alternate hypothesis that the coefficient is
#' larger (more positive) than the value expected for the null hypothesis, or
#' "none" to perform no permutation test; only returning the value of the
#' coefficient without testing information.
#' 
#' Optionally, the permuted coefficient values can be returned, allowing one to
#' plot the frequency distribution of the coefficient under the null hypothesis
#' of the absence of autocorrelation.
#' 
#' @author Guillaume Guénard \email{guillaume.guenard@umontreal.ca}
#' 
#' @seealso Function \code{\link{sp.cov}}, which can be used to obtain a weight
#' matrix.
#' 
#' @references
#' Legendre, P. and L. Legendre. 2012. Numerical ecology, 3rd English edition.
#' Elsevier Science BV, Amsterdam.
#' 
#' [Moran's I ref here...]
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
MoranI <- function(x, w, test=c("two-tail","lower-tail","upper-tail","none"),
                   nperm=999L, saveIp=FALSE) {
  ## test="two-tail"
  test <- match.arg(test)
  storage.mode(x) <- "double"
  storage.mode(w) <- "double"
  if(test=="none")
    nperm <- 0L
  out <- .C("moran",x,w,length(x),double(1L),as.integer(nperm),
            if(nperm) integer(3L) else integer(),as.integer(saveIp),
            if(saveIp) double(nperm) else double())
  Ip <- out[[8L]]
  tval <- out[[6L]]
  out <- out[[4L]]
  if(nperm) {
    pval <- if(test == "upper-tail") {
      tval[3L]/(nperm+1L)
    } else if(test == "lower-tail") {
      tval[1L]/(nperm+1L)
    } else if(test == "two-tail")
      (tval[1L]+tval[3L])/(nperm+1L)
  }
  if(nperm) {
    attr(out,"tval") <- tval
    attr(out,"pval") <- pval
  }
  attr(out,"test") <- test
  attr(out,"call") <- match.call()
  if(saveIp)
    attr(out,"Permuted_I") <- Ip
  return(out)
}
