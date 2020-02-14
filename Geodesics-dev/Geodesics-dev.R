##
## rm(list=ls())
##
compile <- function() {
  try(dyn.unload("src/geodesics.so"),silent=TRUE)
  system("R CMD SHLIB src/geodesics.c")
  dyn.load("src/geodesics.so")
  source("R/geodesics.R")
}
compile()
##
### Harversine
##
### Parameters:
radius <- 6.371e6         ## in meters
##
### Geographic coordinates:
Lat <- c(44.50, 45.50)    ## Latitude of the points
Lon <- c(73.40, 73.40)    ## Longitude of the points
##
### Nearly antipodal points:
## Lat <- c(0, 0.5) ; Lon <- c(0,179.5)   ## correct result: 19936288.579m
## Lat <- c(0, 0.5) ; Lon <- c(0,179.7)   ## correct result: 19944127.421m
##
### Calculations:
Lat <- pi*Lat/180
Lon <- pi*Lon/180
##
delta.lat <- Lat[2L] - Lat[1L]
delta.lon <- Lon[2L] - Lon[1L]
a <- sin(delta.lat/2)^2 +
  cos(Lat[1L])*cos(Lat[2L])*sin(delta.lon/2)^2
d <- 2 * asin(min(1,sqrt(a))) * radius
## 111194.9 ## d-19936288.579 == 172.0293 ## d-19944127.421 == 6122.366
##
coords <- cbind(seq(42,44.5,0.5),rep(c(72,73),3))
##
compile()
res <- .C(
  "dist_geo_hvs",
  as.double(coords),
  as.double(coords),
  c(6L,6L),
  1L,
  double(15L),
  6.371e6
)
res
##
res <- .C(
  "dist_geo_hvs",
  coords[1L:3L,],
  coords[4L:5L,],
  c(3L,2L),
  0L,
  double(6L),
  6.371e6
)
res
##
## rm(list=ls())
##
### Vincenty
##
### Parameters:
a <- 6378137.0            ## Semi-major axis in metres (WGS-84)
f <- 1/298.257223563      ## Flattening factor (WGS-84)
maxiter <- 1024                   ## Maxiumum number of iteraction
tol <- .Machine$double.eps^0.75   ## Calculation tolerance
##
### Geographic coordinates:
Lat <- c(44.50, 45.50)    ## Latitude of the points
Lon <- c(73.40, 73.40)    ## Longitude of the points
##
### Nearly antipodal points:
## Lat <- c(0, 0.5) ; Lon <- c(0,179.5)   ## correct result: 19936288.579m
## Lat <- c(0, 0.5) ; Lon <- c(0,179.7)
##
### Calculations:
b <- (1 - f) * a           ## 6356752.314245 meters in WGS-84
Lat <- pi*Lat/180
Lon <- pi*Lon/180
L <- diff(Lon)            ## Difference in longitude of two points
U <- atan((1-f)*tan(Lat)) ## Reduced latitude (latitude on the auxiliary sphere)
##
### Set lambda to L as a starting point
lambda <- L
##
i <- 1L
### Iteractions begin here
repeat {
  sin_sigma <- ((cos(U[2L])*sin(lambda))^2 +
    (cos(U[1L])*sin(U[2L]) - sin(U[1L])*cos(U[2L])*cos(lambda))^2)^0.5
  cos_sigma <- sin(U[1L])*sin(U[2L]) + cos(U[1L])*cos(U[2L])*cos(lambda)
  sigma <- atan2(sin_sigma,cos_sigma)
  sin_alpha <- cos(U[1L])*cos(U[2L])*sin(lambda)/sin_sigma
  cos_2sigma_m <- cos_sigma - 2*sin(U[1L])*sin(U[2L])/(1-sin_alpha^2)
  alpha <- asin(sin_alpha)
  C <- f*cos(alpha)^2*(4+f*(4-3*cos(alpha)^2))/16
  lambda_prev <- lambda
  lambda <- L + (1-C)*f*sin_alpha*
    (sigma + C*sin(sigma)*(cos_2sigma_m + C*cos_sigma*(-1+2*cos_2sigma_m^2)))
  if((i>maxiter) || (abs(lambda-lambda_prev)<=tol)) break else i <- i + 1L
}
##
u2 <- cos(alpha)^2*((a^2-b^2)/b^2)
A <- 1 + u2/16384*(4096+u2*(-768+u2*(320-175*u2)))
B <- u2/1024*(256+u2*(-128+u2*(74-47*u2)))
deltaSigma <- B*sin_sigma*
  (cos_2sigma_m + 1/4*B*
     (cos_sigma*(-1+2*cos_2sigma_m^2) -
        B/6*cos_2sigma_m*(-3+4*sin_sigma^2)*(-3+4*cos_2sigma_m^2)))
s <- b*A*(sigma-deltaSigma)
az <- c(
  atan2(cos(U[2L])*sin(lambda),
        cos(U[1L])*sin(U[2L]) - sin(U[1L])*cos(U[2L])*cos(lambda)),
  atan2(cos(U[1L])*sin(lambda),
        -sin(U[1L])*cos(U[2L]) + cos(U[1L])*sin(U[2L])*cos(lambda))
)
## 180*az/pi
## Difference: 19936288.579-s  ## Gives the correct result!!!
##
### Alternative calculations:
k1 <- ((1 + u2)^0.5 - 1)/((1 + u2)^0.5 + 1)  ## ()^0.5 == sqrt()
A <- (1 + 0.25*k1^2)/(1-k1)                  ## 1/2 : 0.5
B <- k1*(1-0.375*k1^2)                       ## 3/8 : 0.375
deltaSigma <- B*sin_sigma*
  (cos_2sigma_m + 0.25*B*                    ## 1/4 : 0.25
     (cos_sigma*(-1+2*cos_2sigma_m^2) -
        B/6*cos_2sigma_m*(-3+4*sin_sigma^2)*(-3+4*cos_2sigma_m^2)))
s <- b*A*(sigma-deltaSigma)
az <- c(
  atan2(cos(U[2L])*sin(lambda),
        cos(U[1L])*sin(U[2L]) - sin(U[1L])*cos(U[2L])*cos(lambda)),
  atan2(cos(U[1L])*sin(lambda),
        -sin(U[1L])*cos(U[2L]) + cos(U[1L])*sin(U[2L])*cos(lambda))
)
## Difference: 19936288.579-s  ## Gives the correct result!!!
##
coords <- cbind(c(44.50,45.50,0,0.5,0.5),
                c(73.40,73.40,0,179.5,179.7))
##
compile()
res <- .C(
  "dist_geo_vif",
  as.double(coords),
  as.double(coords),
  c(5L,5L),
  1L,
  double(10L),
  integer(10L),
  a,
  f,
  1024L,
  .Machine$double.eps^0.75
)
res
##
res <- .C(
  "dist_geo_vif",
  coords[c(1L,3L,5L),],
  coords[c(2L,4L),],
  c(3L,2L),
  0L,
  double(6L),
  integer(6L),
  a,
  f,
  1024L,
  .Machine$double.eps^0.75
)
res
## 19936288.579-res[[5L]][4L]
##
N <- 15000L
coords <- cbind(runif(N,-90,90),runif(N,-180,180))
tms <- system.time(
  {
    res <- .C(
      "dist_geo_vif",
      as.double(coords),
      as.double(coords),
      c(N,N),
      1L,
      double(N*(N-1)/2),
      integer(N*(N-1)/2),
      a,
      f,
      1024L,
      .Machine$double.eps^0.75
    )
  }
)
tms
res[[6L]]   ## summary(res[[6L]]) ## sum(res[[9L]]>1024) ## mean(res[[6L]]>1024)
res[[5L]]   ## summary(res[[5L]]) ## summary(res[[5L]][res[[6L]]<=1024])
## rm(res) ; gc()
##
compile()
coords <- cbind(c(44.50,45.50,0,0.5,0.5),
                c(73.40,73.40,0,179.5,179.7))
res <- geodesics(coords)
## str(res)
res <- geodesics(coords[c(1L,3L,5L),],coords[c(2L,4L),])
## str(res)
res <- geodesics(coords,method="V")
## str(res)
res <- geodesics(coords[c(1L,3L,5L),],coords[c(2L,4L),],method="V")
## str(res)
19936288.579-res[2L,2L]
##
N <- 15000L
coords <- cbind(runif(N,-90,90),runif(N,-180,180))
tms <- system.time({res <- geodesics(coords,method="V")})
rm(res)
##
### Package examples
##
### First example: locations spread throughout the world
coords <- cbind(c(43,22,9,12,-40,72,-86,-22),
                c(-135,22,0,1,-45,12,27,-139))
res_hav <- geodesics(coords)  ## Default: the haversine formula
res_hav
res_vif <- geodesics(coords,method = "Vincenty")
res_vif
attr(res_vif,"niter") ## The numbers of iterations
res_vif-res_hav       ## Absolute difference
200*(res_vif-res_hav)/(res_vif+res_hav) ## Large relative difference
##
### Second example: locations nearer from one another
coords <- cbind(c(45.01,44.82,45.23,44.74),
                c(72.03,72.34,71.89,72.45))
res_hav <- geodesics(coords)
res_vif <- geodesics(coords,method = "Vincenty")
res_vif-res_hav       ## Absolute difference
200*(res_vif-res_hav)/(res_vif+res_hav) ## Relative difference are smaller
##
## Implement Vincenty's solution to nearly antipodal location pairs whereby
## we can't achieve convergence using the standard calculations.
## Implement the simple solution for points along the equator or
## Implement the simple solution for points at the same latitude.
##
## rm(list=ls())
##
a <- 6378137.0            ## Semi-major axis in metres (WGS-84)
f <- 1/298.257223563      ## Flattening factor (WGS-84)
##
b <- a*(1 - f)
##
f <- function(x,a,b) sqrt(b*b - b*b*x*x/(a*a))
x <- seq(0,a,length.out=1000L)
plot(x=x, y=f(x,a,b),type="l")
##
g <- function(phi,a,b) cbind(x=a*cos(phi),y=b*sin(phi))
phi <- seq(0,pi/2,length.out=10L)
xy <- g(phi,a,b)
points(xy,pch=21L,bg="red")
##
2*pi*a  ## 40,008 km ok
##
### Calculation issue when both (Lat[1L]==0) && (Lat[2L]==0)
Lon <- c(-90,90)
Lon <- pi*Lon/180
Lat <- c(0,0)
Lat <- pi*Lat/180
deltaLon <- diff(Lon)
s <- a*cos(Lat[1L])*deltaLon
## s <- a*L in the c code
##
### Would be best to implement the Newton method proposed by Karney (2013) as
### it works for any pair of input points. I attempted to implement Vicenty
### modification for nearly antipodal points without success (1.25 work day
### wasted on that attempt.)
##







