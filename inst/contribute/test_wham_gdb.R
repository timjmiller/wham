# use TMB::gdbsource locally if you've changed the .cpp and now it doesn't compile
# two differences:
#   1. compile using "-O0 -g" flag
#   2. can't use fit_wham, must use TMB::MakeADFun
#   3. gdbsource('filename.R')

# This file uses the 'here' package to work across computers, as long as you 
#   *open from the top-level wham folder* (where you git cloned wham into).
# If you do not have 'here' installed:
# install.packages("here")
library(here)
pkg.dir <- here()

# Compile model with "-O0 -g" flag, allows dbs debugger to be used
library(TMB)
compile(file.path(pkg.dir, "src", "wham_v0.cpp"),"-O0 -g")

# ---------------------------------------------------------
# here's an example, you can create a new file in a new location
gdbsource(file.path(pkg.dir,"inst","contribute","run_gdb.R"), TRUE)
