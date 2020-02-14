##
##    (c) 2020 Guillaume Guénard
##        Université de Montréal, Montreal, Quebec, Canada
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
##    R testing and development file
##
## rm(list=ls())     ## Clean the workspace
library(magrittr)
library(sf)
library(sp)
compile <- function() {
  try(dyn.unload("src/SVMP.so"),silent=TRUE)
  try(dyn.unload("src/matrix.so"),silent=TRUE)
  system("R CMD SHLIB src/SVMP.c")
  system("R CMD SHLIB src/matrix.c")
  dyn.load("src/SVMP.so")
  dyn.load("src/matrix.so")
  source("R/SVMP.R")
  source("SVMP-aux.R")
}
compile()
##
### Display the SCF functions
##
if(FALSE) {
  ## X11(width=7.0,height=6.0)
  dst <- seq(0,5,0.01)
  alpha <- c(1,0.75,0.5,0.25,0.1,0.0)
  beta <- 4
  par(mfrow=c(3L,2L))
  mar <- list(spherical=c(3,5,2,0),exponential=c(3,3,2,2),power=c(4,5,1,0),
              hyperbolic=c(4,3,1,2),superelliptic=c(5,5,0,0))
  xlab <- c(spherical="",exponential="",power="",
            hyperbolic="italic(d)",superelliptic="italic(d)")
  for(f in eval(formals(sp.cov)$type)) {
    par(mar=mar[[f]])
    plot(y=sp.cov(dst,f,alpha[1L],beta),x=dst,type="l",las=1L,ylim=c(0,1),
         xlab=parse(text=xlab[f]),ylab="",cex.lab=1.6)
    for(a in 2L:length(alpha))
      lines(y=sp.cov(dst,f,alpha[a],beta),x=dst,lty=a)
    text(5,0.9,toupper(f),1,font=2L,cex=1.5)
  }
  par(mar=c(2,3,1,2))
  plot.new()
  legend(x=0,y=1,lty=1L:6L,title=expression(beta==4),ncol=2,
         legend=parse(text=sprintf("alpha == %0.2f",alpha)),cex=1.5)
  mtext(expression(italic(w)[list(alpha,beta)](italic(d))),
        2L,-1.75,TRUE,0.53)
  rm(dst,alpha,beta,mar,f,a)
  ## dev.copy2eps(file="../Image/MEM_weighting_functions.eps") ; dev.off()
}
##
### Create a 2D set of sampling points:
if(FALSE) {
  X11(height=7, width=7)
  par(mar=c(4.25,4.25,0.5,0.5))
  plot(NA,xlim=c(-10,10),ylim=c(-10,10),asp=1,
       xlab="Eastings (km)",ylab="Northings (km)")
  sp <- 2.5 ; sigma <- 0.15
  dat <- list(
    x_even=seq(-10,10,sp),
    x_odd=c(rev(seq(-0.5*sp,-10,-sp)),seq(0.5*sp,10,sp)),
    y=c(rev(seq(0,-10,-0.5*sp*sqrt(3))),seq(0.5*sp*sqrt(3),10,0.5*sp*sqrt(3)))
  )
  pts <- matrix(NA,0,2L,dimnames=list(NULL,c("x","y")))
  for(i in 1L:length(dat$y))
    pts <- rbind(pts,cbind(dat[[i%%2+1]],dat$y[i]))
  reg <- pts
  points(reg)
  repeat {
    pts[] <- pts + rnorm(length(pts),0,sigma)
    pts[pts>10] <- 10 ; pts[pts<(-10)] <- -10
    plot(NA,xlim=c(-10,10),ylim=c(-10,10),asp=1,
         xlab="Eastings (km)",ylab="Northings (km)")
    points(pts)
    if(is.null(locator(1L))) break
  }
  plot(pts,xlim=c(-11,11),ylim=c(-11,11),asp=1,
       xlab="Eastings (km)",ylab="Northings (km)")
  km <- kmeans(pts,16L,nstart=1000L)
  points(km$centers,pch=3L)
  save(reg,pts,km,file="spdat.rda")
  rm(sp,sigma,dat,i)
  dev.copy2eps(file="../Image/Plot_example.eps") ; dev.off()
} else load(file="spdat.rda")
##
### Begin by having a framework the works good for individual observations
dat <- reg %>% ## pts %>% 
  SpatialPointsDataFrame(
    data=as.data.frame(matrix(NA,nrow(.),0L))
  ) %>% list(spdf=.)
## dat[["spdf"]]@coords
## dat[["spdf"]]@data
##
dat[["d"]] <- EuclidAB(a = dat[["spdf"]]@coords, b = dat[["spdf"]]@coords)
dat[["par"]] <- list(type="spherical",     alpha=0.25, beta=15)
## dat[["par"]] <- list(type="exponential",   alpha=0.5, beta=15)
## dat[["par"]] <- list(type="power",         alpha=0.5, beta=15)
## dat[["par"]] <- list(type="hyperbolic",    alpha=0.5, beta=15)
## dat[["par"]] <- list(type="superelliptic", alpha=0.5, beta=15)
##
dat[["DWmat"]] <- c(dat["d"], dat[["par"]]) %>% do.call(sp.cov,.)
##
if(FALSE) {
  ## Plotting an example of pairwise (non-zero) connections
  X11(width=7.25,height=7.25)
  par(mar=c(4.25,4.25,0.5,0.5))
  dat[["spdf"]]@coords %>%
    plot.scmat(.,.,dat[["DWmat"]], wh=c(15,45), xlim=c(-10,10), ylim=c(-10,10),
               asp=1, xlab="Eastings (km)", ylab="Northings (km)")
}
##
## Centred on the columns of the matrix
dat[["DWmatRC"]] <- dat[["DWmat"]] %>% center(TRUE)
dat[["svdDWmatRC"]] <- dat[["DWmatRC"]] %>% svdTrimMP
## dat[["svdDWmatRC"]]$logDet
##
## This one is the classical definition (double-centred):
## (dat[["svdDWmatRC"]]$d - 1)*nrow(dat[["spdf"]]@coords)/sum(dat[["DWmat"]])
## dat[["svdDWmatRC"]]$d - dat[["DWmatRC"]]%>%eigen%>%{.$values[-length(.$values)]}
## round(dat[["svdDWmatRC"]]$u - recenter(dat[["DWmatRC"]],dat[["DWmat"]])%*%
##         dat[["svdDWmatRC"]]$v%*%diag(dat[["svdDWmatRC"]]$d^-1),8)
##
### It is also possible to center only on the columns: the covariance factors
### matrix would remain positive semi-definite anyway.
### Here, I went with the classical double-centring scheme before transpose-
### multiplying, so the singular values remained linked to the Moran's I. This
### is not a necessity, column-only centring works too, but double centring does
### no harm with predictions, either from individual points or multiple points
### kernels.
##
## A square grid to calculate the eigenfunctions in-between the points of the
## spatial eigenvectors
dat[["grid"]] <- list(grain = 250L)
dat[["grid"]][["fineGrd"]] <-
  cbind(x=rep(seq(-12,12,length.out=dat[["grid"]][["grain"]]),
              each=dat[["grid"]][["grain"]]),
        y=rep(seq(-12,12,length.out=dat[["grid"]][["grain"]]),
              N=dat[["grid"]][["grain"]]))
dat[["grid"]][["d"]] <-
  EuclidAB(a=dat[["spdf"]]@coords, b=dat[["grid"]][["fineGrd"]])
dat[["grid"]][["DWmat"]] <- 
  c(dat[["grid"]]["d"], dat[["par"]]) %>% do.call(sp.cov,.)
dat[["grid"]][["UpRC"]] <-
  dat[["DWmatRC"]] %>% recenter(dat[["grid"]][["DWmat"]]) %*%
  dat[["svdDWmatRC"]]$v %*% diag(dat[["svdDWmatRC"]]$d^-1)
##
if(FALSE) {
  X11(width=7.25,height=7.25)
  par(mar=c(4.25,4.25,0.5,0.5))
  plot.eigenfunctions(
    svdDWmat=dat[["svdDWmatRC"]],
    Up=dat[["grid"]][["UpRC"]],
    dstGrd=dat[["grid"]][["d"]],
    beta=dat[["par"]]$beta,
    grain=dat[["grid"]]$grain,
    pts=dat[["spdf"]]@coords,
    DWmat=dat[["DWmat"]])
}
##
## dat[["par"]]
profile <- list(
  conditions=list(
    type  = "spherical",
    ## type  = "exponential",
    ## type  = "power",
    ## type  = "hyperbolic",
    ## type  = "superelliptic",
    alpha = 10^seq(log10(0.01),log10(1),length.out=100L),
    beta  = 10^seq(log10(1),log10(100),length.out=100L)
  )
)
##
dat[["spdf"]]@data$y <-
  MASS::mvrnorm(10L,rep(0,nrow(dat[["svdDWmatRC"]]$sigma)),
                Sigma=dat[["svdDWmatRC"]]$sigma) %>% t
##
profile[["matrix"]] <-
  matrix(NA,length(profile[["conditions"]][["alpha"]]),
         length(profile[["conditions"]][["beta"]]))
for(i in 1L:length(profile[["conditions"]][["alpha"]]))
  for(j in 1L:length(profile[["conditions"]][["beta"]])) {
    profile[["matrix"]][i,j] <-
      objf_test(
        y=dat[["spdf"]]@data$y,d=dat[["d"]],
        f=profile[["conditions"]][["type"]],
        alpha=profile[["conditions"]][["alpha"]][i],
        beta=profile[["conditions"]][["beta"]][j]
      )
  }
##
plot.profile(profile)
##
library(DEoptim)
library(GenSA)
##
parOpt <- list()
for(f in formals(sp.cov)$type %>% eval) {
  ## f="spherical"
  parOpt[[f]] <- list()
  parOpt[[f]]$opt <-
    optim(par=c(-1,1), fn=objfw_test, method="L-BFGS-B",
          lower=c(-2,0), upper=c(0,2), y=dat[["spdf"]]@data$y,
          d=dat[["d"]], f=f)
  parOpt[[f]]$DEopt <-
   DEoptim(fn = objfw_test, lower=c(-2,0), upper=c(0,2),
           y=dat[["spdf"]]@data$y, d=dat[["d"]], f=f)
  parOpt[[f]]$SA <-
    optim(par=c(-1,1), fn=objfw_test_bounded, method="SANN",
          lb=c(0.01,1), ub=c(1,100), y=dat[["spdf"]]@data$y,
          d=dat[["d"]], f=f)
} ; rm(f)
##
lapply(parOpt,function(x) {
  cat("Difference: ",x$DEopt$optim$bestval-x$opt$value,x$SA$value-x$opt$value,"\n")
  cat("Optim:      ",c(x$opt$value,10^x$opt$par),"\n")
  cat("DEoptim:    ",c(x$DEopt$optim$bestval,10^x$DEopt$optim$bestmem),"\n")
  cat("OptimSA:    ",c(x$SA$value,c(0.01,1)+(c(1,100)-c(0.01,1))/(exp(-x$SA$par)+1)),"\n\n")
  invisible(NULL)
})
##
### On my system, R's C headers are located in /usr/share/R/include/








##
### For the Oribatid data...
load("../Data/Oribates/mite.rds")
mite

















### 2. The non-classical framework where groupes of closely-packed observations
###    are replaced by their centroids, which effectively work as kernels
##
dstK <- EuclidAB(a = km$centers, b = pts)  ## d2
DWmatK <- dstK%>%sp.cov(scf,alpha,beta)    ## w2
##
if(FALSE) {
  X11(width=7.25,height=7.25)
  par(mar=c(4.25,4.25,0.5,0.5))
  plot.scmatK(km$centers,pts,DWmatK,wh=c(1,5),xlim=c(-10,10),ylim=c(-10,10),asp=1,
              xlab="Eastings (km)",ylab="Northings (km)")
}
##
DWmatKC <- DWmatK %>% center(FALSE)
svdDWmatKC <- DWmatKC %>% svdTrimMP    ## svdDWmatKRC$logDet
## (svdDWmatKC$d - 1)*length(pts)/sum(DWmatK)
##
DWmatKRC <- DWmatK %>% center(TRUE)
svdDWmatKRC <- DWmatKRC %>% svdTrimMP      ## svdDWmatKC$logDet
## (svdDWmatKRC$d - 1)*length(pts)/sum(DWmatK)
## round(svdDWmatKRC$u - recenter(DWmatKRC,DWmatK)%*%svdDWmatKRC$v%*%diag(svdDWmatKRC$d^-1),8)
##
dstGrdK <- EuclidAB(a=km$centers, b=fineGrd)
DWmatGrdK <- dstGrdK %>% sp.cov(scf,alpha,beta)
UpKC <- DWmatKC %>% recenter(DWmatGrdK) %*% svdDWmatKC$v %*% diag(svdDWmatKC$d^-1)
UpKRC <- DWmatKRC %>% recenter(DWmatGrdK) %*% svdDWmatKRC$v %*% diag(svdDWmatKRC$d^-1)
##
if(FALSE) {
  X11(width=7.25,height=7.25)
  par(mar=c(4.25,4.25,0.5,0.5))
  plot.eigenfunctionsK(svdDWmatKRC,UpKRC,dstGrdK,beta,grain,pts,DWmatK,km)
}
##
### Seems to work fine enough
##
### The covariance matrices associated with these problems are:
svdDWmatC$sigma   ## C1m <- U1%*%diag(lambda1^2)%*%t(U1)
svdDWmatC$inv     ## C1m_inv <- U1%*%diag(lambda1^-2)%*%t(U1)
svdDWmatC$logDet  ## logDetC1m <- sum(2*log(lambda1))
svdDWmatC$rank    ## rankC1m <- length(lambda1)
##
### or, upon double centring:
svdDWmatRC$sigma
svdDWmatRC$inv
svdDWmatRC$logDet
svdDWmatRC$rank
##
### When using kernels instead of the individual points
svdDWmatKC$sigma
svdDWmatKC$inv
svdDWmatKC$logDet
svdDWmatKC$rank
##
### or, upon double centring:
svdDWmatKRC$sigma
svdDWmatKRC$inv
svdDWmatKRC$logDet
svdDWmatKRC$rank
##
### Is there a smooth MLE for dmax and alpha given that framework?
### The profile is flat for any constant term because of the centering
### of the weight matrix. It won't be so for the other model coefficients.
### Verify that all is well.
##
























