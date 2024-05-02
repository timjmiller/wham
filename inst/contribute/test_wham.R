# Run WHAM tests on your local repo copy (ie not yet committed and pushed to github)
# If all of your tests pass:
#   Great!
# If your tests do not pass:
#   Bummer... the tests do not save objects for debugging, but the tests now provide some information on where failures are occuring within each test. 
# To run all example scripts, see inst/contribute/run_all_examples_local.R.
# If you've changed the .cpp and now it doesn't compile:
#   Check out test_wham_gdb.R and run_gdb.R to use TMB::gdbsource

# This file uses the 'here' package to work across computers, as long as you 
#   *open from the top-level wham folder* (where you git cloned wham into).
# If you do not have 'here' installed:
# install.packages("here")
library(here)
pkg.dir <- here()

# confirm that the .cpp compiles and wham can be loaded
devtools::load_all(path = pkg.dir)

# run all tests
devtools::test(pkg = pkg.dir)

# to run 1 test at a time
# testthat::test_file(file.path(pkg.dir,"tests","testthat","test_ex1_SNEMAYT.R"))
# testthat::test_file(file.path(pkg.dir,"tests","testthat","test_ex2_CPI_SNEMAYT.R"))
# testthat::test_file(file.path(pkg.dir,"tests","testthat","test_ex3_projections.R"))
# testthat::test_file(file.path(pkg.dir,"tests","testthat","test_ex4_selectivity.R"))
# testthat::test_file(file.path(pkg.dir,"tests","testthat","test_ex5_M.R"))
# testthat::test_file(file.path(pkg.dir,"tests","testthat","test_ex6_NAA.R"))
# testthat::test_file(file.path(pkg.dir,"tests","testthat","test_ex8_compare.R"))
# testthat::test_file(file.path(pkg.dir,"tests","testthat","test_ex9_retro_pred.R"))
# testthat::test_file(file.path(pkg.dir,"tests","testthat","test_ex11_catchability.R"))
