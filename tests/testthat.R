library(testthat)
library(wham)

# Check that test results are available locally
path_to_examples <- system.file("extdata", package="wham")
file.exists(file.path(path_to_examples,"ex1_test_results.rds"))

testthat::test_check("wham")
