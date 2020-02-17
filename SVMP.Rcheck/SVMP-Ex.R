pkgname <- "SVMP"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "SVMP-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('SVMP')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("EuclidAB")
### * EuclidAB

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: EuclidAB
### Title: Calculation of the Euclidean Distance Among Two Sets of Points
### Aliases: EuclidAB

### ** Examples

##
### First example: locations spread throughout the world
##
## [Examples here...]




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("EuclidAB", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("geodesics")
### * geodesics

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: geodesics
### Title: Calculation Of Geodesic Distances
### Aliases: geodesics

### ** Examples

##
### First example: locations spread throughout the world
##
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
##
coords <- cbind(c(45.01,44.82,45.23,44.74),
                c(72.03,72.34,71.89,72.45))
res_hav <- geodesics(coords)
res_vif <- geodesics(coords,method = "Vincenty")
res_vif-res_hav       ## Absolute difference
200*(res_vif-res_hav)/(res_vif+res_hav) ## Relative difference are smaller
##




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("geodesics", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sp.cov")
### * sp.cov

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sp.cov
### Title: Calculation of Spatial Covariance Models
### Aliases: sp.cov

### ** Examples

##
### First example: locations spread throughout the world
##
## [Examples here...]




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sp.cov", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
