##
### Examples of confounding spatial structure.
##
## rm(list=ls())
##
compile <- function() {
  try(dyn.unload("../src/spectR.so"),silent=TRUE)
  system("R CMD SHLIB ../src/spectR.c")
  dyn.load("../src/spectR.so")
  source("../R/spectR.R")
  source("../spectR-aux.R")
  library(magrittr)
  library(sf)
  library(sp)
}
compile()
##
pts <- cbind(x=c(0,1,5,8,11,12,15),y=c(1,5,2,4,3,1,2)) %>%
  SpatialPointsDataFrame(
    data=as.data.frame(matrix(NA,nrow(.),0L))
  )
plot(pts)
d <- EuclidAB(a = pts@coords, b = pts@coords)
p <- list(type="sup", alpha=0.85, beta=15)
DWmat <- c(list(d), p) %>% do.call(sp.cov,.)
##
if(FALSE) {
  ## Plotting an example of pairwise (non-zero) connections
  X11(width=7.25,height=7.25)
  par(mar=c(4.25,4.25,0.5,0.5))
  pts@coords %>%
    plot.scmat(.,.,DWmat, wh=c(3), xlim=c(-1,16), ylim=c(0,6),
               asp=1, xlab="Eastings (km)", ylab="Northings (km)")
}
##
DWmatRC <- DWmat %>% center(TRUE)
svdDWmatRC <- DWmatRC %>% svdTrimMP
## svdDWmatRC$logDet
## (svdDWmatRC$d - 1)*nrow(pts@coords)/sum(DWmat)
##
grid <- list(grain = 250L)
grid[["fineGrd"]] <-
  cbind(x=rep(seq(-3,18,length.out=grid[["grain"]]),
              each=grid[["grain"]]),
        y=rep(seq(-3,18,length.out=grid[["grain"]]),
              N=grid[["grain"]]))
grid[["d"]] <-
  EuclidAB(a=pts@coords, b=grid[["fineGrd"]])
grid[["DWmat"]] <- 
  c(grid["d"], p) %>% do.call(sp.cov,.)
grid[["UpRC"]] <-
  DWmatRC %>% recenter(grid[["DWmat"]]) %*%
  svdDWmatRC$v %*% diag(svdDWmatRC$d^-1)
##
if(FALSE) {
  X11(width=7.25,height=7.25)
  par(mar=c(4.25,4.25,0.5,0.5))
  plot.eigenfunctions(
    svdDWmat=svdDWmatRC,
    Up=grid[["UpRC"]],
    dstGrd=grid[["d"]],
    beta=p$beta,
    grain=grid$grain,
    pts=pts@coords,
    DWmat=DWmat,
    xlim=c(-3,18),
    ylim=c(-3,18))
}
if(FALSE) {
  plot.eigenfunctionsRGL(
    svdDWmat=svdDWmatRC,
    Up=grid[["UpRC"]],
    dstGrd=grid[["d"]],
    beta=p$beta,
    grain=grid$grain,
    pts=pts@coords,
    DWmat=DWmat,
    xlim=c(-3,18),
    ylim=c(-3,18),3)
}





