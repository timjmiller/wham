# Run all WHAM examples and save output for debugging
#   uses ex scripts modified to work from wham source directory

# source(file.path(here::here(),"inst","contribute","run_all_examples_local.R"))

# This file uses the 'here' package to work across computers, as long as you 
#   *open from the top-level wham folder* (where you git cloned wham into).
# If you do not have 'here' installed:
# install.packages("here")
library(here)
pkg.dir <- here()

# choose where to save output
# I have a .gitignore'd folder 'sandbox' in the wham folder
main.dir <- file.path(pkg.dir, "sandbox", paste0("runall-",format(Sys.Date(), "%Y%m%d")))
if(!dir.exists(main.dir)) dir.create(main.dir)

# load local wham, recompile if cpp changed
devtools::load_all(path = pkg.dir)

# Ex 1
write.dir <- file.path(main.dir,"ex1")
if(!dir.exists(write.dir)) dir.create(write.dir)
source(file.path(pkg.dir, "inst", "contribute", "ex1_basics.R"))

# Ex 2
write.dir <- file.path(main.dir,"ex2")
if(!dir.exists(write.dir)) dir.create(write.dir)
source(file.path(pkg.dir, "inst", "contribute", "ex2_CPI_recruitment.R"))

# Ex 3
write.dir <- file.path(main.dir,"ex3")
if(!dir.exists(write.dir)) dir.create(write.dir)
source(file.path(pkg.dir, "inst", "contribute", "ex3_projections.R"))

# Ex 4
write.dir <- file.path(main.dir,"ex4")
if(!dir.exists(write.dir)) dir.create(write.dir)
source(file.path(pkg.dir, "inst", "contribute", "ex4_selectivity.R"))

# Ex 5
write.dir <- file.path(main.dir,"ex5")
if(!dir.exists(write.dir)) dir.create(write.dir)
source(file.path(pkg.dir, "inst", "contribute", "ex5_M_GSI.R"))

# Ex 6
write.dir <- file.path(main.dir,"ex6")
if(!dir.exists(write.dir)) dir.create(write.dir)
source(file.path(pkg.dir, "inst", "contribute", "ex6_NAA.R"))

# Ex 8
write.dir <- file.path(main.dir,"ex8")
if(!dir.exists(write.dir)) dir.create(write.dir)
source(file.path(pkg.dir, "inst", "contribute", "ex8_compare_asap.R"))

# Ex 9
write.dir <- file.path(main.dir,"ex9")
if(!dir.exists(write.dir)) dir.create(write.dir)
source(file.path(pkg.dir, "inst", "contribute", "ex9_retro_pred.R"))

# Ex 10
write.dir <- file.path(main.dir,"ex10")
if(!dir.exists(write.dir)) dir.create(write.dir)
source(file.path(pkg.dir, "inst", "contribute", "ex10_om_mse.R"))

