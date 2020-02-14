##
### Development of covariance functions
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
## Spherical: based on an order three polynomial whose derivative increases
## from a value of -1.5 at a distance of 0 to a value of 0 at d=beta (range)
## in a purely quadratic manner: d(w)/dd = -1.5 + 1.5x²
##
spher = function(d,alpha,beta) {
  dd <- d/beta
  n <- alpha / (alpha - 1)
  out <- 1 - alpha*dd + (alpha - 1)*dd^n  ## alpha/n == alpha - 1
  out[d>beta] <- 0
  out
}
plot(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),1,3),type="l")
lines(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),1.2,3),col="red3")
lines(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),1.5,3),col="red2",lty=2L)
lines(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),2.3,3),col="red")
lines(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),10,3),col="orange")
lines(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),-1,3),col="green")
lines(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),-2,3),col="yellow")
lines(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),-10,3),col="turquoise")
lines(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),100,3),lty=3L)
lines(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),-100,3),lty=3L)
lines(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),-0.5,3),col="blue")
lines(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),-0.1,3),col="purple")
lines(x=seq(0,5,0.01),y=spher(seq(0,5,0.01),-0.001,3),col="grey",lty=2L)
1 - (10000*0.5) + (10000-1)*0.5^(10000/(10000-1))
1 - (-10000*0.5) + (-10000-1)*0.5^(-10000/(-10000-1))
0.5*log(0.5) - 0.5 + 1
##
## Limit: alpha toward Inf (or -Inf) : dd * log(dd) - dd + 1
##
x = seq(-1,2,0.01)
## y = 0.5-0.25/(x-0.5) ; y[abs(y)>4] <- NA
## y = (1/2)-1/(4*(x-(1/2))) ; y[abs(y)>4] <- NA
y = (x-1)/(2*x-1) ; y[abs(y)>4] <- NA
plot(x,y,type="l",lty=2L) ; abline(h=0) ; abline(v=0) ; abline(v=0.5,lty=3L)
## z = (y-1)/(2*y-1)
## plot(z,y,type="l",lty=2L) ; abline(h=0) ; abline(v=0) ; abline(v=0.5,lty=3L)
rect(xleft=-1,xright=2,ybottom=0,ytop=1,col=rgb(0.5,0.5,0.5,0.5))
abline(h=c(-2,2))
0.5-0.25/(0.25-0.5)
##
hyper = function(d,alpha,beta) {  ## -1/alpha changed for -alpha
  a <- alpha
  if(is.finite(alpha)) {
    out <- ((1+d)^(1/a) - (1+beta)^(1/a))/(1 - (1+beta)^(1/a))
  } else {
    out <- 1 - log(d + 1)/log(beta + 1)
  }
  out[d>beta] <- 0
  out
}
plot(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),1,3),type="l")
lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),1.2,3),col="red3")
lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),1.5,3),col="red2",lty=2L)
lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),2.3,3),col="red")
lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),10,3),col="orange")
## lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),100,3),lty=3L)
lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),Inf,3),lty=3L,lwd=3)
## lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),-100,3),lty=3L)
lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),-10,3),col="turquoise")
lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),-2.3,3),col="yellow")
lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),-1,3),col="green")
lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),-0.5,3),col="blue")
lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),-0.1,3),col="purple")
lines(x=seq(0,5,0.01),y=hyper(seq(0,5,0.01),-0.001,3),col="grey",lty=2L)
##
### Re-thinking the exponential
##
expon = function(d,alpha,beta) {
  out <- exp(-(d/beta))
  out[d>beta] <- 0
  out
}
plot(x=seq(0,5,0.01),y=expon(seq(0,5,0.01),1,4),type="l")
## Derivative: f'(d) = -exp(-(d/beta))/beta
## The point where f'(d) = -1/beta is that where exp(-(d/beta)) = 1
## d = -beta*log(1) = 0, obviously!
## When alpha = 1: we must take the exponential at d = 0 and stretch it
## between d = 0 and d = beta.
## When alpha = 0: we must take the exponential between d = 0 and d = +Inf and
## shrink it between d = 0 and d = beta.
tr <- function(x) x/(1 - x)
plot(x=seq(0,1,0.01),y=tr(seq(0,1,0.01)),type="l")
##
expon2 = function(d,alpha,beta) {
  a <- (alpha-1)/alpha
  out <- (exp(a*(d/beta)) - exp(a))/(1 - exp(a))
  out[d>beta] <- 0
  out
}
plot(x=seq(0,5,0.01),y=expon2(seq(0,5,0.01),0.1,4),type="l",ylim=c(0,1))
##
### Testing the weighting functions
##
dst <- seq(0,5,0.01)
wflist <- list(
  ## Inherited from variogram
  spher = function(d,alpha,beta) {
    dd <- d/beta
    if(alpha==0.5) {
      out <- dd * log(dd) - dd + 1
      out[d==0] <- 1
    } else if(alpha > 0) {
      a <- 0.5 - 0.25/(0.5 - alpha)
      n <- a / (a - 1)
      out <- 1 - a*dd + (a - 1)*dd^n
    } else {
      out <- numeric(length(d))
      out[d==0] <- 1
    }
    out[d>beta] <- 0
    out
  },
  ## Inherited from variogram
  ## Parametrization allows to cover the Gaussian case as well
  expon = function(d,alpha,beta) {
    if(alpha<1) {
      if(alpha>0) {
        a <- (alpha-1)/alpha
        out <- (exp(a*(d/beta)) - exp(a))/(1 - exp(a))
      } else {
        out <- numeric(length(d))
        out[d==0] <- 1
      }
    } else {
      out <- 1-d/beta
    }
    out[d>beta] <- 0
    out
  },
  ## Works generally fine
  ## Erratic behaviour (non-smooth likelihood surface) when alpha > 0.7
  power = function(d,alpha,beta) {
    out <- (1-d/beta)^(1/alpha)
    out[d>beta] <- 0
    out
  },
  ## This one behaves rather strangely but might perform well at time...
  hyper = function(d,alpha,beta) {  ## -1/alpha changed for -alpha
    a <- alpha
    if(alpha == 0.5) {
      out <- 1 - log(d + 1)/log(beta + 1)
    } else if (alpha > 0) {
      a <- 0.5 - 0.25/(0.5 - alpha)
      out <- ((1+d)^(1/a) - (1+beta)^(1/a))/(1 - (1+beta)^(1/a))
    } else {
      out <- numeric(length(d))
      out[d==0] <- 1
    }
    out[d>beta] <- 0
    out
  },
  ## Works generally fine
  super = function(d,alpha,beta) {  ## Works fine
    if(alpha > 0) {
      out <- 1 - (1/beta)*(beta^(1/alpha) - (beta-d)^(1/alpha))^alpha
      out[d>beta] <- 0
    } else {
      out <- numeric(length(d))
      out[d==0] <- 1
    }
    out
  }
)
##
## plot(x=seq(0,5,0.01),y=wflist[["spher"]](seq(0,5,0.01),1,3),type="l")
## plot(x=seq(0,5,0.01),y=wflist[["expon"]](seq(0,5,0.01),1,3),type="l")
## plot(x=seq(0,5,0.01),y=wflist[["hyper"]](seq(0,5,0.01),0.25,3),type="l")
## plot(x=seq(0,5,0.01),y=MEMweight(seq(0,5,0.01),"hyper",0.25,3),type="l")
## plot(x=seq(0,5,0.01),y=MEMweight(seq(0,5,0.01),"expon",0.25,3),type="l")
alpha <- c(1,0.75,0.5,0.25,0.1,0.0)
beta <- 4
par(mar=c(5,5,3,2))
for(f in names(wflist)) {
  ## f=names(wflist)[4L]
  plot(y=wflist[[f]](dst,alpha[1L],beta),x=dst,type="l",
       ylim=c(0,1),xlab=expression(italic(d)),main=sprintf("Function: %s",f),
       ylab=expression(italic(w)[list(alpha,beta)](italic(d))))
  for(a in 2L:length(alpha))
    lines(y=wflist[[f]](dst,alpha[a],beta),x=dst,lty=a)
  if(is.null(locator(1L))) break
  for(a in 1L:length(alpha))
    lines(y=sp.cov(dst,f,alpha[a],beta),x=dst,lty=a,col=2,lwd=2)
  if(is.null(locator(1L))) break
  for(a in 1L:length(alpha))
    lines(y=sp.cov(dst,f,rep(alpha[a],2L),beta),x=dst,lty=a,col=3,lwd=2)
  if(is.null(locator(1L))) break
}
rm(dst,wflist,alpha,beta,f,a)
##

