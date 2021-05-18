##
### Toward a Hilbert operator to represent asymmetric procecesses
##
### Temporal processes are inherently asymmetrical as time cannot be reversed.
### Also, some spatial problems involve an asymmetric distance metric. For
### instance, observations performed along a stream of running water, where
### dissemination is made harder when going upstream, and easier when going
### downstream, then in the absence of a water flow. Wind might also have a
### similar effect. An REML estimator may also be developed to test hypothesis
### about influencial processes using similar ideas.
##
### The main interest I could see is to amend the MCA to test codependence in
### terms of coherence and phase relationship.
##
## rm(list=ls())
N <- 100L
tol <- .Machine$double.eps^0.5
##
m <- matrix(0+0i,N,N)
m[cbind(1L:(N-1L),2L:N)] <- complex(mod=1,arg=0.5*pi/N)
m[cbind(2L:N,1L:(N-1L))] <- complex(mod=1,arg=-0.5*pi/N)
## m[1L,N] <- complex(mod=1,arg=-pi/50)
## m[N,1L] <- complex(mod=1,arg=pi/50)
mm <- list(row=rowMeans(m),col=colMeans(m),mean=mean(m))
##
mc <- t(t(m) - mm$col) - mm$row + mm$mean
## rowMeans(mc) ## colMeans(mc) ## mean(mc)
emc <- eigen(mc, symmetric = TRUE)
emc$vectors <- emc$vectors[,abs(emc$values) >= tol]
emc$values <- emc$values[abs(emc$values) >= tol]
## round(t(Conj(emc$vectors)) %*% emc$vectors,3)
##
i <- 1L
par(mfrow=c(2L,1L))
while(i <= 9) {  ## length(emc$values)
  par(mar=c(1,5,4,2))
  plot(emc$values,type="l",xaxt="n",main=sprintf("U[%d]",i))
  points(x=i,y=emc$values[i],pch=21L,bg="black")
  par(mar=c(5,5,1,2))
  plot(Re(emc$vectors[,i]),type="l",
       ylim=max(abs(range(Re(emc$vectors[,i]),Im(emc$vectors[,i]))))*c(-1,1))
  lines(Im(emc$vectors[,i]),col="red")
  ## if(is.null(locator(1L))) break
  dev.copy2eps(file=sprintf("./Image/U_%02d.eps",i))
  i <- i + 1L
}
##
U <- cbind(1,emc$vectors)
y <- sin(3.1*pi*(0:(N - 1L))/N) + sin(3.3*pi*(0:(N - 1L))/N) +
  sin(13.3*pi*(0:(N - 1L))/N) + sin(15.1*pi*(0:(N - 1L))/N) +
  rnorm(N,0,0.1)
b <- solve(t(Conj(U))%*%U) %*% t(Conj(U)) %*% y
##
par(mfrow=c(3L,1L), par(mar=c(4,5,2,2)))
plot(y, type="l", xlab="Location", ylab="Response", las=1L)
par(mar=c(4,5,1,2))
plot(Mod(b[-1L]), type="l", xlab="Index", ylab="Power", las=1L)
par(mar=c(4.5,5,1,2))
plot(180*Arg(b[-1L])/pi, ylim=c(-180,180), type="l", xlab="Index", ylab="Phase",
     las=1L)
dev.copy2eps(file="./Image/Signal_power.eps")
##
i <- 5L
par(mfrow=c(2L,1L))
for(a in seq(0, 2*pi, length.out=100L)) {
  par(mar=c(1,5,4,2))
  plot(emc$values, type="l", xaxt="n", main=sprintf("U[%d]",i), las=1L)
  points(x=i,y=emc$values[i], pch=21L, bg="black")
  par(mar=c(5,5,1,2))
  bi <- complex(mod=1, arg=a)
  yhat <- emc$vectors[,i,drop=FALSE] %*% as.matrix(bi)
  plot(Re(yhat), type="l", ylim=c(-0.15,0.15), las=1L)
  lines(Im(yhat), col="red", lty=3L)
  Sys.sleep(0.05)
}
##
par(mfrow=c(1L,1L))
plot(b)
##
### Spiral representation
rgl::plot3d(NA, xlim=c(0,N-1L), zlim=c(-0.05,0.05), ylim=c(-0.05,0.05),
            xlab="", zlab="", ylab="", box=FALSE, axes=FALSE)
rgl::axes3d(edges = c("x--", "y--", "z-+"))
rgl::title3d(xlab = "Time/Space", ylab = "Im(U[i])")
rgl::mtext3d("Re(U[i])", "z-+", line = 2)
rgl::segments3d(x=c(0,N-1L),z=c(0,0),y=c(0,0))
i <- 1L
rgl::lines3d(x=0:(N - 1L), z=Re(emc$vectors[,i]), y=Im(emc$vectors[,i]),
             color="red", lwd=3)
i <- 2L
rgl::lines3d(x=0:(N - 1L), z=Re(emc$vectors[,i]), y=Im(emc$vectors[,i]),
             color="orange", lwd=3)
i <- 5L
rgl::lines3d(x=0:(N - 1L), z=Re(emc$vectors[,i]), y=Im(emc$vectors[,i]),
             color="green", lwd=3)
i <- 15L
rgl::lines3d(x=0:(N - 1L), z=Re(emc$vectors[,i]), y=Im(emc$vectors[,i]),
             color="blue", lwd=3)
i <- 35L
rgl::lines3d(x=0:(N - 1L), z=Re(emc$vectors[,i]), y=Im(emc$vectors[,i]),
             color="purple", lwd=3)
##
### Works as expected.
##
### The set of coefficients b is not conjugate complementary from end to end, as
### is the case for FT's orthogonal bases.
###
##
### Two-demensional example
##
library(magrittr)
##
N <- c(10,10)
##
cbind(
  rep(scale(1:N[1L], center=TRUE, scale=FALSE), N[2L]),
  rep(scale(1:N[2L], center=TRUE, scale=FALSE), each=N[1L])
) -> coords
##
coords %>%
  dist %>%
  as.matrix -> m
##
mn <- m %>% apply(1L,function(x) min(x[x!=0])) %>% max
##
mag <- 1.5
a <- pi*45/180
th <- sqrt(2)
w <- 1
##
plot(coords, asp=1)
arrows(x0=0, x1=mag*cos(a), y0=0, y1=mag*sin(a))
##
m <- matrix(0,nrow(coords),nrow(coords))
## i=2 ; j=1
for(i in 2:nrow(coords))
  for(j in 1:(i-1)) {
    ## points(x=coords[i,1],coords[i,2],pch=21,bg="red")
    ## points(x=coords[j,1],coords[j,2],pch=21,bg="green")
    d <- sum((coords[i,]-coords[j,])^2)^0.5
    aa1 <- atan2(coords[i,2]-coords[j,2],coords[i,1]-coords[j,1])
    aa2 <- atan2(coords[j,2]-coords[i,2],coords[j,1]-coords[i,1])
    if(d <= th) {
      m[i,j] <- complex(
        real=1-(d/th)^w,
        imaginary=mag*cos(a-aa1)/nrow(coords)
      )
      m[j,i] <- complex(
        real=1-(d/th)^w,
        imaginary=mag*cos(a-aa2)/nrow(coords)
      )
    }
  }
## m[1:5,1:5]
##
mm <- list(row=rowMeans(m),col=colMeans(m),mean=mean(m))
##
mc <- t(t(m) - mm$col) - mm$row + mm$mean
## rowMeans(mc) ## colMeans(mc) ## mean(mc)
emc <- eigen(mc, symmetric = TRUE)
emc$vectors <- emc$vectors[,abs(emc$values) >= tol]
emc$values <- emc$values[abs(emc$values) >= tol]
## round(t(Conj(emc$vectors)) %*% emc$vectors,3)
##
i <- 1L
par(mar=c(5,5,2,2))
plot(x=coords[,1],y=coords[,2],asp=1,pch=21,
     bg=c("black","white")[0.5*(sign(Re(emc$vectors[,i]))+1)+1],
     cex=10*abs(Re(emc$vectors[,i])))
##
arrows(x0=0, x1=mag*cos(a), y0=0, y1=mag*sin(a))
##
points(x=coords[,1]+0.5*cos(a),y=coords[,2]+0.5*sin(a),asp=1,pch=21,
       bg=c("red","green")[0.5*(sign(Im(emc$vectors[,i]))+1)+1],
       cex=10*abs(Im(emc$vectors[,i])))

emc$vectors[,i]





while(i <= 9) {  ## length(emc$values)
  par(mar=c(1,5,4,2))
  plot(emc$values,type="l",xaxt="n",main=sprintf("U[%d]",i))
  points(x=i,y=emc$values[i],pch=21L,bg="black")
  par(mar=c(5,5,1,2))
  plot(Re(emc$vectors[,i]),type="l",
       ylim=max(abs(range(Re(emc$vectors[,i]),Im(emc$vectors[,i]))))*c(-1,1))
  lines(Im(emc$vectors[,i]),col="red")
  ## if(is.null(locator(1L))) break
  dev.copy2eps(file=sprintf("./Image/U_%02d.eps",i))
  i <- i + 1L
}
##
U <- cbind(1,emc$vectors)
y <- sin(3.1*pi*(0:(N - 1L))/N) + sin(3.3*pi*(0:(N - 1L))/N) +
  sin(13.3*pi*(0:(N - 1L))/N) + sin(15.1*pi*(0:(N - 1L))/N) +
  rnorm(N,0,0.1)
b <- solve(t(Conj(U))%*%U) %*% t(Conj(U)) %*% y
##
par(mfrow=c(3L,1L), par(mar=c(4,5,2,2)))
plot(y, type="l", xlab="Location", ylab="Response", las=1L)
par(mar=c(4,5,1,2))
plot(Mod(b[-1L]), type="l", xlab="Index", ylab="Power", las=1L)
par(mar=c(4.5,5,1,2))
plot(180*Arg(b[-1L])/pi, ylim=c(-180,180), type="l", xlab="Index", ylab="Phase",
     las=1L)
dev.copy2eps(file="./Image/Signal_power.eps")
##
i <- 5L
par(mfrow=c(2L,1L))
for(a in seq(0, 2*pi, length.out=100L)) {
  par(mar=c(1,5,4,2))
  plot(emc$values, type="l", xaxt="n", main=sprintf("U[%d]",i), las=1L)
  points(x=i,y=emc$values[i], pch=21L, bg="black")
  par(mar=c(5,5,1,2))
  bi <- complex(mod=1, arg=a)
  yhat <- emc$vectors[,i,drop=FALSE] %*% as.matrix(bi)
  plot(Re(yhat), type="l", ylim=c(-0.15,0.15), las=1L)
  lines(Im(yhat), col="red", lty=3L)
  Sys.sleep(0.05)
}
##
par(mfrow=c(1L,1L))
plot(b)
##
### Spiral representation
rgl::plot3d(NA, xlim=c(0,N-1L), zlim=c(-0.05,0.05), ylim=c(-0.05,0.05),
            xlab="", zlab="", ylab="", box=FALSE, axes=FALSE)
rgl::axes3d(edges = c("x--", "y--", "z-+"))
rgl::title3d(xlab = "Time/Space", ylab = "Im(U[i])")
rgl::mtext3d("Re(U[i])", "z-+", line = 2)
rgl::segments3d(x=c(0,N-1L),z=c(0,0),y=c(0,0))
i <- 1L
rgl::lines3d(x=0:(N - 1L), z=Re(emc$vectors[,i]), y=Im(emc$vectors[,i]),
             color="red", lwd=3)
i <- 2L
rgl::lines3d(x=0:(N - 1L), z=Re(emc$vectors[,i]), y=Im(emc$vectors[,i]),
             color="orange", lwd=3)
i <- 5L
rgl::lines3d(x=0:(N - 1L), z=Re(emc$vectors[,i]), y=Im(emc$vectors[,i]),
             color="green", lwd=3)
i <- 15L
rgl::lines3d(x=0:(N - 1L), z=Re(emc$vectors[,i]), y=Im(emc$vectors[,i]),
             color="blue", lwd=3)
i <- 35L
rgl::lines3d(x=0:(N - 1L), z=Re(emc$vectors[,i]), y=Im(emc$vectors[,i]),
             color="purple", lwd=3)
##



