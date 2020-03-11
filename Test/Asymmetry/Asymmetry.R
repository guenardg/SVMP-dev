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
m[cbind(1L:(N-1L),2L:N)] <- complex(mod=1,arg=pi/50)
m[cbind(2L:N,1L:(N-1L))] <- complex(mod=1,arg=-pi/50)
## m[1L,N] <- complex(mod=1,arg=-pi/50)
## m[N,1L] <- complex(mod=1,arg=pi/50)
mm <- list(row=rowMeans(m),col=colMeans(m),mean=mean(m))
##
mc <- t(t(m) - mm$col) - mm$row + mm$mean
## rowMeans(mc) ## colMeans(mc) ## mean(mc)
emc <- eigen(mc)
emc$vectors <- emc$vectors[,abs(emc$values) >= tol]
emc$values <- emc$values[abs(emc$values) >= tol]
## round(t(Conj(emc$vectors)) %*% emc$vectors,3)
##
i <- 1L
par(mfrow=c(2L,1L))
while(i <= length(emc$values)) {
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
par(mfrow=c(2L,1L), par(mar=c(4,5,2,2)))
plot(y, type="l", xlab="Location", ylab="Response")
par(mar=c(4.5,5,1,2))
plot(Mod(b[-1L]), type="l", xlab="Index", ylab="Power")
dev.copy2eps(file="./Image/Signal_power.eps")
##
i <- 1L
par(mfrow=c(2L,1L))
for(a in seq(0, 2*pi, length.out=100L)) {
  par(mar=c(1,5,4,2))
  plot(emc$values, type="l", xaxt="n", main=sprintf("U[%d]",i))
  points(x=i,y=emc$values[i], pch=21L, bg="black")
  par(mar=c(5,5,1,2))
  bi <- complex(mod=1, arg=a)
  yhat <- emc$vectors[,i,drop=FALSE] %*% as.matrix(bi)
  plot(Re(yhat), type="l", ylim=c(-0.25,0.25))
  lines(Im(yhat), col="red", lty=3L)
  Sys.sleep(0.1)
}
##
par(mfrow=c(1L,1L))
plot(b)
##
### Works as expected.
##
### The set of coefficients b is not conjugate complementary from end to end, as
### is the case for FT's orthogonal bases.
###
##
