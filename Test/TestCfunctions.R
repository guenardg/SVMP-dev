##
### Test C functions
##
rm(list=ls())     ## Clean the wworkspace
library(magrittr)
compile <- function() {
  try(dyn.unload("../Analysis/src/spectR.so"),silent=TRUE)
  system("R CMD SHLIB ../Analysis/src/spectR.c")
  dyn.load("../Analysis/src/spectR.so")
  source("../Analysis/R/spectR.R")
}
compile()
##
### Testing mat_center
Q <- function(n) diag(n)-matrix(1/n,n,n)
compile()
mat <- matrix(c(1,2,5,3,7,5,0,3,6,5,2,5),4L,3L)
res <- .C("mat_center",mat,nrow(mat),ncol(mat),double(nrow(mat)),
          double(ncol(mat)),double(1L),1L)
rbind(colMeans(mat),res[[5L]])
rbind(rowMeans(mat),res[[4L]])
c(mean(mat),res[[6L]])
Q(4) %*% mat %*% Q(3)
colMeans(res[[1L]])
rowMeans(res[[1L]])
res[[1L]]
resgc <- .C("get_center",mat,nrow(mat),ncol(mat),double(nrow(mat)),
            double(ncol(mat)),double(1L),1L)
## svd(resgc[[1L]])
(mat - resgc[[1L]]) - res[[1L]]
if(FALSE) {
  mat <- matrix(rnorm(1e+6*200,0,1),1e+6,200)
  res <- .C("mat_center",mat,nrow(mat),ncol(mat),double(nrow(mat)),
            double(ncol(mat)),double(1L),0L)
}
## Rather fast to execute: will no be a computational bottleneck
rm(mat, res, Q) ; gc()
##
### Testing EuclidAB and the sp.cov
x <- seq(-1,1,0.01)
d <- EuclidAB(x,x)
w <- d%>%sp.cov("spher",0.5,1)
wcc <- w%>%center(TRUE)
eig_wcc <- wcc%>%eigen
svd_wcc <- wcc%>%svd
## round(eig_wcc$values-svd_wcc$d,8)  ## Check for equality
## length(x)*(svd_wcc$d-1)/sum(w)     ## Associated Moran's I
##
### Moran's I generalizes for "external" references
xx <- x[seq(2L,length(x),2L)]
dd <- EuclidAB(xx,x)                  ## dim(dd)
ww <- dd%>%sp.cov("spher",0.5,1)
wwcc <- ww%>%center(TRUE)
svd_wwcc <- wwcc%>%svd
## svd_wwcc$d
## length(x)*(svd_wwcc$d-1)/sum(ww)   ## Associated Moran's I
##
rm(x,d,w,eig_wcc,svd_wcc,xx,dd,ww,wwcc,svd_wwcc)
##
### Testing the re-centering routine:
compile()
x <- matrix(c(1,2,5,3,7,5,0,3,6,5,2,5),4L,3L)
a <- center(x,FALSE)
b <- .C("mat_recenter",x,nrow(x),ncol(x),double(),attr(a,"colMeans"),double(),0L)
b[[1L]]
a
a <- center(x,TRUE)
b <- .C("mat_recenter",x,nrow(x),ncol(x),double(nrow(x)),
        attr(a,"colMeans"),attr(a,"mean"),1L)
b[[1L]]
a
recenter(center(x,FALSE),x)
recenter(center(x,TRUE),x)
recenter(center(x,TRUE),x,FALSE)
rm(x,a,b)
##
load(file="../Analysis/spdat.rda")
##
### Exemplifying how predictions are made
### 1. Classical framework where all observation points are kernels
##
d1 <- EuclidAB(a = pts, b = pts)
f <- "spher" ; alpha <- 0.75 ; beta <- 8
## f <- "expon" ; alpha <- 0.25 ; beta <- 8
## f <- "power" ; alpha <- 0.25 ; beta <- 8
## f <- "hyper" ; alpha <- 0.25 ; beta <- 8
## f <- "super" ; alpha <- 0.25 ; beta <- 8
w1 <- d1%>%sp.cov(f,alpha,beta)
##
X11(width=7.25,height=7.25)
par(mar=c(4.25,4.25,0.5,0.5))
plot(NA,xlim=c(-10,10),ylim=c(-10,10),asp=1,
     xlab="Eastings (km)",ylab="Northings (km)")
for(i in 2L:nrow(pts)) {
  ## i=6L
  jj <- !!w1[i,1L:i]
  if(any(jj))
    for(j in which(jj))
      segments(
        x0=pts[i,1L],x1=pts[j,1L],
        y0=pts[i,2L],y1=pts[j,2L],lty=3L,
        lwd=if(xor(i==45,j==45)) 1.5 else 0.5,
        col=if(xor(i==45,j==45)) "black" else "grey"
      )
}
points(pts,pch=21L,bg="black")
##
w1cc <- w1%>%center(TRUE)
svd_w1cc <- w1cc%>%svd
## svd_w1cc$d
U1 <- svd_w1cc$u[,svd_w1cc$d>.Machine$double.eps^0.5]
V1 <- svd_w1cc$v[,svd_w1cc$d>.Machine$double.eps^0.5]
lambda1 <- svd_w1cc$d[svd_w1cc$d>.Machine$double.eps^0.5]
## dim(svd_w1cc$u) ## dim(svd_w1cc$v)
## round(svd_w1cc$u%*%diag(svd_w1cc$d)%*%t(svd_w1cc$v) - w1cc,8)
## round(U1%*%diag(lambda1)%*%t(V1) - w1cc,8)
## 
## (lambda1 - 1)*length(pts)/sum(w1)  ## (l-1)*n/t(ones(n))%*%w%*%ones(n)
## round(U - recenter(wcc,w)%*%V%*%diag(lambda^-1),8)
grain <- 250L
fineGrd <- cbind(x=rep(seq(-12,12,length.out=grain),each=grain),
                 y=rep(seq(-12,12,length.out=grain),N=grain))
dd1 <- EuclidAB(a=pts, b=fineGrd)    ## dim(dd)
ww1 <- dd1%>%sp.cov(f,alpha,beta)
## dim(ww)
U1p <- w1cc%>%recenter(ww1)%*%V1%*%diag(lambda1^-1)
##
## X11(width=7.25,height=7.25)
k <- 1L
zrng <- range(U1[,k],U1p[,k]) ; cols <- rainbow(1200L)[1L:1000L]
par(mar=c(4.25,4.25,0.5,0.5))
tmp <- U1p[,k]
tmp[rowSums(dd1<=beta)==0L] <- NA
image(z=matrix(tmp,grain,grain,byrow=TRUE),zlim=zrng,
      x=seq(-12,12,length.out=grain),xlim=c(-12,12),xlab="Eastings (km)",
      y=seq(-12,12,length.out=grain),ylim=c(-12,12),ylab="Northings (km)",
      col=cols,asp=1,axes=FALSE)
axis(1L) ; axis(2L)
for(i in 2L:nrow(pts)) {
  ## i=6L
  jj <- !!w1[i,1L:i]
  if(any(jj))
    for(j in which(jj))
      segments(
        x0=pts[i,1L],x1=pts[j,1L],
        y0=pts[i,2L],y1=pts[j,2L],lty=3L,lwd=0.5,col="grey"
      )
}
points(pts,pch=21L,
       bg=cols[floor(999*(U1[,k] - zrng[1L])/(zrng[2L] - zrng[1L])) + 1L])
## dev.copy2eps(file="../Image/SpatialWeightingExample1.eps") ; dev.off()
##
### 2. The non-classical framework where groupes of closely-packed observations
###    are replaced by their centroids, which effectively work as kernels
##
d2 <- EuclidAB(a = km$centers, b = pts)
w2 <- d2%>%sp.cov(f,alpha,beta+4)
##
par(mar=c(4.25,4.25,0.5,0.5))
plot(NA,xlim=c(-10,10),ylim=c(-10,10),
     xlab="Eastings (km)",ylab="Northings (km)")
for(i in 1L:nrow(km$centers)) {
  jj <- !!w2[,i]
  if(any(jj))
    for(j in which(jj))
      segments(
        x0=km$centers[i,1L],x1=pts[j,1L],
        y0=km$centers[i,2L],y1=pts[j,2L],
        lty=3L,lwd=0.5,col="grey"
      )
}
points(pts,pch=21L,bg="black")
points(km$centers,pch=3L,cex=2)
##
w2cc <- w2%>%center(TRUE)
svd_w2cc <- w2cc%>%svd
## svd_w2cc$d
U2 <- svd_w2cc$u[,svd_w2cc$d>.Machine$double.eps^0.5]
V2 <- svd_w2cc$v[,svd_w2cc$d>.Machine$double.eps^0.5]
lambda2 <- svd_w2cc$d[svd_w2cc$d>.Machine$double.eps^0.5]
## dim(svd_w2cc$u) ## dim(svd_w2cc$v)
## round(svd_w2cc$u%*%diag(svd_w2cc$d)%*%t(svd_w2cc$v) - w2cc,8)
## round(U2%*%diag(lambda2)%*%t(V2) - w2cc,8)
## 
## (lambda2 - 1)*length(pts)/sum(w2)  ## (l-1)*n/t(ones(n))%*%w%*%ones(n)
dd2 <- EuclidAB(a=km$centers, b=fineGrd)  ## dim(dd)
ww2 <- dd2%>%sp.cov(f,alpha,beta+4)
## dim(ww2)
U2p <- w2cc%>%recenter(ww2)%*%V2%*%diag(lambda2^-1)
##
## X11(width=7.25,height=7.25)
k <- 1L
zrng <- range(U2[,k],U2p[,k]) ; cols <- rainbow(1200L)[1L:1000L]
par(mar=c(4.25,4.25,0.5,0.5))
tmp <- U2p[,k]
## tmp[rowSums(dd2<=beta)==0L] <- NA
image(z=matrix(tmp,grain,grain,byrow=TRUE),zlim=zrng,
      x=seq(-12,12,length.out=grain),xlim=c(-12,12),xlab="Eastings (km)",
      y=seq(-12,12,length.out=grain),ylim=c(-12,12),ylab="Northings (km)",
      col=cols,asp=1,axes=FALSE)
axis(1L) ; axis(2L)
for(i in 1L:nrow(km$centers)) {
  jj <- !!w2[,i]
  if(any(jj))
    for(j in which(jj))
      segments(
        x0=km$centers[i,1L],x1=pts[j,1L],
        y0=km$centers[i,2L],y1=pts[j,2L],
        lty=3L,lwd=0.5,col="grey"
      )
}
points(pts,pch=21L,
       bg=cols[floor(999*(U2[,k] - zrng[1L])/(zrng[2L] - zrng[1L])) + 1L])
points(km$centers,pch=3L,cex=2)
## dev.copy2eps(file="../Image/SpatialWeightingExample2.eps") ; dev.off()
##
### Seems to work fine enough
##
### The covariance matrices associated with these problems are:
C1m <- U1%*%diag(lambda1^2)%*%t(U1)
C1m_inv <- U1%*%diag(lambda1^-2)%*%t(U1)
logDetC1m <- sum(2*log(lambda1))     ## exp(logDetCm)   ## det(Cm)
rankC1m <- length(lambda1)
##
### and
C2m <- U2%*%diag(lambda2^2)%*%t(U2)
C2m_inv <- U2%*%diag(lambda2^-2)%*%t(U2)
logDetC2m <- sum(2*log(lambda2))   ## exp(logDetC2m)  ## det(C2m)
rankC2m <- length(lambda2)
##
### Is there a smooth MLE for dmax and alpha given that framework?
### The profile is flat for any constant term because of the centering
### of the weight matrix. It won't be so for the other model coefficients.
### Verify that all is well.
##
y <- MASS::mvrnorm(1L,rep(0,nrow(C1m)),Sigma=C1m)
##
objf <- function(y, x, d, f, alpha, beta, tol=.Machine$double.eps^0.5) {
  w <- sp.cov(d,f,alpha,beta)
  wcc <- center(w,TRUE)
  svd_wcc <- svd(wcc)
  U <- svd_wcc$u[,svd_wcc$d>tol]
  V <- svd_wcc$v[,svd_wcc$d>tol]
  L <- svd_wcc$d[svd_wcc$d>tol]
  C <- U%*%diag(L^2)%*%t(U)
  C_inv <- U%*%diag(L^-2)%*%t(U)
  logDetC <- sum(2*log(L))     ## exp(logDetCmm)   ## det(Cmm)
  rankC <- length(L)
  ##
  n <- length(y)
  ## Only useful for non-constant columns of x
  if(!missing(x)) {
    b <- solve(t(x)%*%C_inv%*%x)%*%t(x)%*%C_inv%*%y
    yc <- y - x%*%b
  } else {
    yc <- y - mean(y)  ## Not really consequential: no gradient for a constant
  }
  sigma <- t(yc)%*%C_inv%*%yc/n
  return(n + rankC*log(2*pi) + n*log(sigma) + logDetC)
  ##
}
##
## f == "spher"
alpha <- 10^seq(log10(0.01),log10(1),length.out=100L)
beta <- 10^seq(log10(1),log10(100),length.out=100L)
## f == "expon"
alpha <- 10^seq(log10(0.01),log10(1),length.out=100L)
beta <- 10^seq(log10(1),log10(100),length.out=100L)
## f == "power"
alpha <- 10^seq(log10(0.01),log10(1),length.out=100L)
beta <- 10^seq(log10(1),log10(100),length.out=100L)
## f == "hyper"
alpha <- 10^seq(log10(0.01),log10(1),length.out=100L)
beta <- 10^seq(log10(1),log10(100),length.out=100L)
## f == "super"
alpha <- 10^seq(log10(0.01),log10(1),length.out=100L)
beta <- 10^seq(log10(1),log10(100),length.out=100L)
##
profile <- matrix(NA,length(alpha),length(beta))
for(i in 1L:length(alpha))
  for(j in 1L:length(beta))
    profile[i,j] <- objf(y=y,d=EuclidAB(pts,pts),f=f,
                         alpha=alpha[i],beta=beta[j])
image(z=profile,x=alpha,y=beta,log="xy",col=rainbow(1200L)[1L:1000L])
wmp <- which.min(profile)
i <- (wmp-1L)%%nrow(profile)+1L
j <- (wmp-1L)%/%nrow(profile)+1L
profile[i,j]
## min(profile)
points(x=alpha[i],y=beta[j])
##
maxd <- EuclidAB(pts,pts)%>%max
plot(y=sp.cov(seq(0,maxd,length.out=1000L),f,alpha[i],beta[j]),
     x=seq(0,maxd,length.out=1000L),type="l",ylim=c(0,1))
points(x=EuclidAB(pts,pts),y=sp.cov(EuclidAB(pts,pts),f,alpha[i],beta[j]))
##
objfw <- function(par,y,d,f) {
  res <- objf(y=y, d=d, f=f, alpha=10^par[1L], beta=10^par[2L])
  cat(10^par[1L],10^par[2L],":",res,"\n")
  res
}
##
opt <- optim(par=c(-1,1),fn=objfw,method="L-BFGS-B",
             lower=c(-3,0),upper=c(0,2),
             y=y,d=EuclidAB(pts,pts),f="spher") ; f <- "spher"
opt <- optim(par=c(-1,20),fn=objfw,method="L-BFGS-B",
             lower=c(-3,0),upper=c(0,2),
             y=y,d=EuclidAB(pts,pts),f="expon") ; f <- "expon" ## Does not seem smooth
opt <- optim(par=c(-1,1),fn=objfw,method="L-BFGS-B",
             lower=c(-3,0),upper=c(0,2),
             y=y,d=EuclidAB(pts,pts),f="power") ; f <- "power"
opt <- optim(par=c(-1,1),fn=objfw,method="L-BFGS-B",
             lower=c(-3,0),upper=c(0,2),
             y=y,d=EuclidAB(pts,pts),f="hyper") ; f <- "hyper"
opt <- optim(par=c(-1,1),fn=objfw,method="L-BFGS-B",
             lower=c(-3,0),upper=c(0,2),
             y=y,d=EuclidAB(pts,pts),f="super") ; f <- "super"
##
plot(y=MEMweight(seq(0,maxd,length.out=1000L),f,10^opt$par[1L],10^opt$par[2L]),
     x=seq(0,maxd,length.out=1000L),type="l",ylim=c(0,1))
points(x=EuclidAB(pts,pts),
       y=MEMweight(EuclidAB(pts,pts),f,10^opt$par[1L],10^opt$par[2L]))
##
seq(0, 1, 0.01) -> x
plot(y = -0.5*x, x = x, type = "l")
plot(y = 1 - 0.5*x, x = x, type = "l")
plot(y = 1 - x, x = x, type = "l")






