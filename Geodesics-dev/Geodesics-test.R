##
### Testing script
##
## rm(list=ls())
##
compile <- function() {
  try(dyn.unload("../SVMP-dev/src/geodesics.so"),silent=TRUE)
  system("R CMD SHLIB ../SVMP-dev/src/geodesics.c")
  dyn.load("../SVMP-dev/src/geodesics.so")
  source("../SVMP-dev/R/geodesics.R")
}
compile()
##
### Parameters:
## radius <- 6.371e6         ## in meters
## a <- 6378137.0            ## Semi-major axis in metres (WGS-84)
## f <- 1/298.257223563      ## Flattening factor (WGS-84)
## maxiter <- 1024                   ## Maxiumum number of iteraction
## tol <- .Machine$double.eps^0.75   ## Calculation tolerance
##
N <- 1500L
coords <- cbind(runif(N,-90,90),runif(N,-180,180))
tms_hav <- system.time({res_hav <- geodesics(coords,method="h")})
tms_vif <- system.time({res_vif <- geodesics(coords,method="V")})
summary(res_hav)
summary(res_vif)
summary(res_vif-res_hav)
##
