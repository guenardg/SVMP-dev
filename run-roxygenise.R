#!/usr/bin/env Rscript
setwd("SVMP")
system("R CMD SHLIB -o src/SVMP.so src/*.c")
roxygen2::roxygenise()
system("rm -f src/*.so")
system("rm -f src/*.o")
system("rm -rf *~")
system("find -type f \\( -not -name \"MD5\" \\) -exec md5sum '{}' \\; > MD5")
