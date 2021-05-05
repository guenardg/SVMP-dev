##
### Development of MEM and PCNM functions with prediction functionalities
##
## rm(list=ls())     ## Clean the workspace
##
### Would be good to derive a pdist function from package stat dist to help in
### calculating pairwise distances from a pair of matrices x and y rather than
### a single matrix as does the actual dist implementation.
##
library(magrittr)
data("mafragh",package="ade4")
## names(mafragh)
##
xy <- mafragh$xy %>% as.matrix
## plot(xy, asp=1L)
dd <- xy %>% dist
##
dc <- function(x) {
  rm <- rowMeans(x)
  cm <- colMeans(x)
  xm <- mean(x)
  list(xc=t(t(x - rm) - cm) + xm,rm=rm,cm=cm,xm=xm)
}
##
tr <- function(x) (x - min(x))/(max(x) - min(x))
##
sn <- function(x) c("#000000","#FFFFFF")[(sign(x) + 3L)/2]
##
wf <- list(
  PCNM = function(d,threshold,a) {
    w <- d
    w[d>threshold] <- 4*threshold
    w <- -2*w^2
    w
  },
  MEM1 = function(d,threshold,a) {
    w <- d
    w[d<=threshold] <- 1 - w[d<=threshold]/threshold
    w[d>threshold] <- 0
    w
  },
  MEM2 = function(d,threshold,a) {
    w <- d
    w[d<=threshold] <- 1 - (w[d<=threshold]/threshold)^a
    w[d>threshold] <- 0
    w
  },
  MEM3 = function(d,threshold,a) {
    w <- d
    w[d<=threshold] <- 1 / (w[d<=threshold]/threshold)^a
    w[(d==0)|(d>threshold)] <- 0
    w
  }
)
##
plotBubble <- function(xy,U,pch=21L,...)
  plot(xy, pch=pch, bg=sn(U), asp=1, cex=2*tr(U), ...)

##
threshold <- 50
a <- c(PCNM=NA,MEM1=NA,MEM2=0.5,MEM3=0.1)
ww <- list()
xc <- list()
U <- list()
lambda <- list()
##
### Calculations:
##
for(case in names(wf)) {
  ww[[case]] <- dd %>% wf[[case]](threshold,a[case])
  xc[[case]] <- ww[[case]] %>% as.matrix %>% dc
  eigxc <- eigen(xc[[case]]$xc)
  ##
  U[[case]] <- eigxc$vectors[,abs(eigxc$values)>.Machine$double.eps^0.5]
  lambda[[case]] <- eigxc$values[abs(eigxc$values)>.Machine$double.eps^0.5]
}
rm(case,eigxc)
##
par(mfrow=c(2L,2L))
for(case in names(wf)) lambda[[case]] %>% barplot(main=case)
##
par(mfrow=c(2L,2L)) ; i <- 5L
for(case in names(wf)) xy %>% plotBubble(U=U[[case]][,i],main=case)
##
rm(case,i)
##
rng <- xy %>%
  apply(2L, range) %>%
  {list(
    x=seq(floor(.[1L,1L]),ceiling(.[2L,1L])),
    y=seq(floor(.[1L,2L]),ceiling(.[2L,2L]))
  )}
##
da <- array(NA,c(length(rng$y),length(rng$x),nrow(xy)))
wwg <- list()
Ug <- list()
##
## i=1L ; j=1L
for(i in 1L:length(rng$y))
  for(j in 1L:length(rng$x))
    da[i,j,] <- sqrt((xy[,"y"] - rng$y[i])^2 + (xy[,"x"] - rng$x[j])^2)
##
for(case in names(wf)) {
  wwg[[case]] <- da %>% wf[[case]](threshold,a[case])
  Ug[[case]] <- array(NA,c(length(rng$y),length(rng$x),length(lambda[[case]])))
  for(i in 1L:length(rng$y))
    for(j in 1L:length(rng$x)) {
      w <- wwg[[case]][i,j,] %>%
        {. - mean(.) - xc[[case]]$cm + xc[[case]]$xm} %>%
        t
      Ug[[case]][i,j,] <- w %*% U[[case]] %*% diag(lambda[[case]]^(-1))
    }
}
rm(case,i,j,w)
##
i <- 1L
par(mfrow=c(2L,2L))
for(case in names(wf)) {
  s <- seq(0,1,length.out=1000L)
  cols <- c(rgb(1,s,s),rgb(1-s,1-s,1))
  image(z=t(Ug[[case]][,,i]),x=rng$x,y=rng$y,col=cols,xlab="",ylab="",
        zlim=max(abs(Ug[[case]][,,i]))*c(-1,1),main=case)
  points(xy)
}
rm(i,case)
##
### PCNM and MEM3 are not smooth spatial variance estimators
### MEM1-2 are smooth spatial variance estimators
##








