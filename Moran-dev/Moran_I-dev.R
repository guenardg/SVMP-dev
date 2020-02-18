## rm(list=ls())
##
## compile <- function() {
##   try(dyn.unload("../SVMP/src/Moran_I.so"),silent=TRUE)
##   system("R CMD SHLIB ../SVMP/src/Moran_I.c")
##   dyn.load("../SVMP/src/Moran_I.so")
##   source("../SVMP/R/Moran_I.R")
## }
## compile()
##
library(SVMP)
##
x <- c(9,8,7,8,1,1,2,4,1,4,7,8,9,7,3,4,4,5,6,7)
plot(x,type="l")
w <- matrix(0,length(x),length(x))
w[cbind(2L:length(x),1L:(length(x)-1L))] <- 1
w[cbind(1L:(length(x)-1L),2L:length(x))] <- 1
##
if(FALSE) {
  xc <- scale(x,scale=FALSE)
  acc <- 0
  W <- 0
  ssx <- 0
  for(i in 1L:length(x)) {
    ssx <- ssx + xc[i]*xc[i]
    for(h in (1L:length(x))[-i]) {
      acc <- acc + w[h,i]*xc[h]*xc[i]
      W <- W + w[h,i]
    }
  }
  ssx/length(x) ; acc / W
  cat("acc =",acc,"; W =",W,"; ssx =",ssx,"; I =",length(x)*acc/(W*ssx),"\n")
  rm(xc,acc,W,ssx,h,i)
}
##
MoranPerm <- function(x, w, nperm=999,
                      test=c("two-tail","low-tail","up_tail")) {
  test <- match.arg(test)
  xc <- as.matrix(x - mean(x))
  cp <- length(x)/(sum(xc*xc)*sum(w))
  EI <- -1/(length(xc)-1)
  I <- t(xc)%*%w%*%xc*cp-EI
  permv <- integer(3L)
  permv[3L] <- 1L
  for(p in 1L:nperm) {
    i <- sample(1L:length(x),replace=FALSE)
    Ip <- t(xc[i,])%*%w%*%xc[i,]*cp-EI
    if(Ip <= -abs(I)) {
      permv[1L] <- permv[1L] + 1
    } else if(Ip >= abs(I)) {
      permv[3L] <- permv[3L] + 1
    } else permv[2L] <- permv[2L] + 1
  }
  pval <- if(test == "up_tail") {
    permv[3L]/(nperm+1L)
  } else if(test == "low-tail") {
    permv[1L]/(nperm+1L)
  } else (permv[1L]+permv[3L])/(nperm+1L)
  return(list(I,permv,pval,test))
}
##
res <- MoranPerm(x,w)
##
rm(res,MoranPerm)
##
compile()
resc <-
  .C(
    "moran",
    as.double(x),
    as.double(w),
    length(x),
    double(1L),
    999L,
    integer(3L),
    0L,
    double()
    )
resc
##
resc <-
  .C(
    "moran",
    as.double(x),
    as.double(w),
    length(x),
    double(1L),
    999L,
    integer(3L),
    1L,
    double(999L)
  )
resc
rm(resc)
##
res <- MoranI(x,w)
res <- MoranI(x,w,"u")
res <- MoranI(x,w,"l")
res <- MoranI(x,w,"t")
res <- MoranI(x,w,"t",99999L)
res <- MoranI(x,w,"t",99L,TRUE)
res <- MoranI(x,w,"t",9999999L)
##
library(parallel)
cl <- makeForkCluster(8L)
typeI <- parSapplyLB(
  cl,
  1L:100000L,
  function(X,n,nperm) {
    xt <- rnorm(n,0,1)
    tmp <- MoranI(xt,w,"t",nperm)
    return(attr(tmp,"pval"))
  },
  n=length(x),
  nperm=9999L
)
mean(typeI<0.05)
stopCluster(cl)
##
