##
### From Guillaume Guénard - Université de Montréal
### August 2019
### License: GPL version 3
### R code
##
MoranI <- function(x, w,test=c("none","two-tail","lower-tail","upper-tail"),
                   nperm=999L,saveIp=FALSE) {
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
