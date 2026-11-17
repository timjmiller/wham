library(wham)
# Check that test results are available locally
file.exists("fixtures/ex01_test_results.RDS")

testthat::test_check("wham")
