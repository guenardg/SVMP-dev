##
### Auxiliary code for spectR-dev.R
##
plot.scmat <- function(x, y, w, wh = c(), line.col=c("black","grey"), ...) {
  plot(NA, ...)
  for(i in 2L:nrow(x)) {
    ## i=6L
    jj <- !!w[i,1L:i]
    if(any(jj))
      for(j in which(jj))
        segments(
          x0=x[i,1L],x1=y[j,1L],
          y0=x[i,2L],y1=y[j,2L],lty=3L,
          lwd=if(xor(i%in%wh,j%in%wh)) 1.5 else 0.5,
          col=if(xor(i%in%wh,j%in%wh)) line.col[1L] else line.col[2L]
        )
  }
  points(x,pch=21L,bg=line.col[1L])
  invisible(NULL)
}
##
### Trim the too small singular values from an SVD
svdTrimMP <- function(x, tol=.Machine$double.eps^0.5) {
  out <- svd(x)
  discard <- which(out$d<tol)
  if(length(discard)) {
    out$u <- out$u[,-discard]
    out$v <- out$v[,-discard]
    out$d <- out$d[-discard]
  }
  out$sigma <- out$u%*%diag(out$d^2)%*%t(out$u)
  out$inv <- out$u%*%diag(out$d^-2)%*%t(out$u)
  out$logDet <- 2*sum(log(out$d))
  out$rank <- length(out$d)
  out
}
##
plot.eigenfunctions <- function(svdDWmat,Up,dstGrd,beta,grain,pts,DWmat,sl=1,
                                cols=rainbow(1200L)[1L:1000L],
                                xlim=c(-12,12),ylim=c(-12,12)) {
  for(k in 1L:length(svdDWmat$d)) {
    zrng <- range(svdDWmat$u[,k],Up[,k])
    tmp <- Up[,k]
    tmp[rowSums(dstGrd<=beta)==0L] <- NA
    image(z=matrix(tmp,grain,grain,byrow=TRUE),zlim=zrng,
          x=seq(xlim[1L],xlim[2L],length.out=grain),xlim=xlim,xlab="Eastings (km)",
          y=seq(ylim[1L],ylim[2L],length.out=grain),ylim=ylim,ylab="Northings (km)",
          col=cols,asp=1,axes=FALSE)
    axis(1L) ; axis(2L)
    for(i in 2L:nrow(pts)) {
      ## i=6L
      jj <- !!DWmat[i,1L:i]
      if(any(jj))
        for(j in which(jj))
          segments(
            x0=pts[i,1L],x1=pts[j,1L],
            y0=pts[i,2L],y1=pts[j,2L],lty=3L,lwd=0.5,col="grey"
          )
    }
    points(pts,pch=21L,
           bg=cols[floor(999*(svdDWmat$u[,k] - zrng[1L])/(zrng[2L] - zrng[1L])) + 1L])
    ## if(!!scan(what=character(), quiet=TRUE)%>%length) break
    Sys.sleep(sl)
  }
  return(invisible(NULL))
}
##
plot.eigenfunctionsRGL <- function(svdDWmat,Up,dstGrd,beta,grain,pts,DWmat,
                                   xlim=c(-12,12),ylim=c(-12,12),ef=1) {
  rgl.open()
  for(k in 1L:length(svdDWmat$d)) {
    tmp <- Up[,k]
    tmp[rowSums(dstGrd<=beta)==0L] <- NA
    rgl.surface(y=ef*matrix(tmp,grain,grain,byrow=TRUE),
                x=seq(xlim[1L],xlim[2L],length.out=grain),
                z=seq(ylim[1L],ylim[2L],length.out=grain))
    if(length(scan(what=character()))) break
    if(k < length(svdDWmat$d)) rgl.clear()
  }
  return(invisible(NULL))
}
##
objf_test <- function(y, x, d, f, alpha, beta, tol=.Machine$double.eps^0.5) {
  w <- sp.cov(d,f,alpha,beta)
  wcc <- center(w,TRUE)
  svd_wcc <- svd(wcc)
  U <- svd_wcc$u[,svd_wcc$d>tol]
  V <- svd_wcc$v[,svd_wcc$d>tol]
  L <- svd_wcc$d[svd_wcc$d>tol]
  C <- U%*%diag(L^2)%*%t(U)
  C_inv <- U%*%diag(L^-2)%*%t(U)
  logDetC <- 2*sum(log(L))     ## exp(logDetCmm)   ## det(Cmm)
  rankC <- length(L)
  ##
  n <- NROW(y)
  m <- NCOL(y)
  ## Only useful for non-constant columns of x
  if(!missing(x)) {
    b <- solve(t(x)%*%C_inv%*%x)%*%t(x)%*%C_inv%*%y
    yc <- y - x%*%b
  } else {
    yc <- t(t(y) - colMeans(y))
  }
  ## sigma <- (t(yc)%*%C_inv%*%yc)/rankC
  ## return(rankC*(1 + log(2*pi) + log(sigma)) + logDetC)
  res <- 0
  for(i in 1L:m) {
    sigma <- (t(yc[,i])%*%C_inv%*%yc[,i])/n
    res <- res + (n + rankC*log(2*pi) + n*log(sigma) + logDetC)
  }
  return(res)
  ##
}
##
plot.profile <- function(x, ...) {
  image(z=x[["matrix"]],x=x[["conditions"]][["alpha"]],
        y=x[["conditions"]][["beta"]],log="xy",
        col=rainbow(1200L)[1L:1000L])
  wmp <- which.min(x[["matrix"]])
  i <- (wmp-1L)%%nrow(x[["matrix"]])+1L
  j <- (wmp-1L)%/%nrow(x[["matrix"]])+1L
  ## x[["matrix"]][i,j]
  ## min(profile)
  points(x=x[["conditions"]][["alpha"]][i],
         y=x[["conditions"]][["beta"]][j])
}
##
objfw_test <- function(par,y,d,f) {
  res <- objf_test(y=y, d=d, f=f, alpha=10^par[1L], beta=10^par[2L])
  cat(10^par[1L],10^par[2L],":",res,"\n")
  res
}
##
objfw_test_bounded <- function(par,y,d,f,lb,ub) {
  logpar <- lb+(ub-lb)/(exp(-par)+1)
  res <- objf_test(y=y, d=d, f=f, alpha=logpar[1L], beta=logpar[2L])
  cat(alpha,beta,":",res,"\n")
  res
}
##
plot.scmatK <- function(x, y, w, wh = c(), line.col=c("black","grey"), ...) {
  plot(NA, ...)
  for(i in 1L:nrow(x)) {
    jj <- !!w[,i]
    if(any(jj))
      for(j in which(jj))
        segments(
          x0=x[i,1L],x1=y[j,1L],
          y0=x[i,2L],y1=y[j,2L],lty=3L,
          lwd=if(i%in%wh) 1.5 else 0.5,
          col=if(i%in%wh) line.col[1L] else line.col[2L]
        )
  }
  points(y,pch=21L,bg=line.col[1L])
  points(x,pch=21L, bg="red")
  invisible(NULL)
}
##
plot.eigenfunctionsK <- function(svdDWmatK,UpK,dstGrdK,beta,grain,pts,DWmatK,
                                km,sl=1,cols=rainbow(1200L)[1L:1000L],
                                xlim=c(-12,12),ylim=c(-12,12)) {
  for(k in 1L:length(svdDWmatK$d)) {
    ## k=1L
    zrng <- range(svdDWmatK$u[,k],UpK[,k])
    tmp <- UpK[,k]
    tmp[rowSums(dstGrdK<=beta)==0L] <- NA
    image(z=matrix(tmp,grain,grain,byrow=TRUE),zlim=zrng,
          x=seq(xlim[1L],xlim[2L],length.out=grain),xlim=xlim,xlab="Eastings (km)",
          y=seq(ylim[1L],ylim[2L],length.out=grain),ylim=ylim,ylab="Northings (km)",
          col=cols,asp=1,axes=FALSE)
    axis(1L) ; axis(2L)
    for(i in 1L:nrow(km$centers)) {
      ## i=6L
      jj <- !!DWmatK[i,1L:i]
      if(any(jj))
        for(j in which(jj))
          segments(
            x0=km$centers[i,1L],x1=pts[j,1L],
            y0=km$centers[i,2L],y1=pts[j,2L],lty=3L,lwd=0.5,col="grey"
          )
    }
    points(pts,pch=21L,
           bg=cols[floor(999*(svdDWmatK$u[,k] - zrng[1L])/
                           (zrng[2L] - zrng[1L])) + 1L])
    points(km$centers,pch=21L, bg="red")
    Sys.sleep(sl)
  }
  return(invisible(NULL))
}
##




