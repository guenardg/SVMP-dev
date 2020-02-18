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
plot(x, type="l")
##
### Connections between immediate neighbors
w <- matrix(0, length(x), length(x))
w[cbind(2L:length(x),1L:(length(x)-1L))] <- 1
w[cbind(1L:(length(x)-1L),2L:length(x))] <- 1
##
res <- MoranI(x, w, "t")
##
N <- 1500L
x_big <- rnorm(N,0,1)
##
### Connections between immediate neighbors
w_big <- matrix(0,N,N)
w_big[cbind(2L:N,1L:(N-1L))] <- 1
w_big[cbind(1L:(N-1L),2L:N)] <- 1
##
tms_big <- system.time({res_big <- MoranI(x_big, w_big, "t")})
tms_big
res_big
## rm(x_big,w_big)
##
