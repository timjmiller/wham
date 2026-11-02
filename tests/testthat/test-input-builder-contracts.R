small_input <- function(catch_info = NULL, index_info = NULL, F = NULL) {
  suppressMessages(prepare_wham_input(
    basic_info = list(years = 2000:2005, ages = 1:4),
    NAA_re = list(sigma = "rec"),
    catch_info = catch_info,
    index_info = index_info,
    F = F
  ))
}

test_that("prepared inputs have consistent observation and parameter dimensions", {
  input <- small_input()

  expect_identical(input$years, 2000:2005)
  expect_identical(input$data$n_years_model, 6L)
  expect_identical(input$data$n_ages, 4L)
  expect_equal(input$data$n_fleets, 1L)
  expect_equal(input$data$n_indices, 1L)
  expect_identical(dim(input$data$agg_catch), c(6L, 1L))
  expect_identical(dim(input$data$catch_paa), c(1L, 6L, 4L))
  expect_identical(dim(input$data$agg_indices), c(6L, 1L))
  expect_identical(dim(input$data$index_paa), c(1L, 6L, 4L))
  expect_identical(dim(input$par$F_pars), c(6L, 1L))
  expect_equal(
    apply(input$data$catch_paa, c(1L, 2L), sum),
    matrix(1, 1L, 6L)
  )
  expect_equal(
    apply(input$data$index_paa, c(1L, 2L), sum),
    matrix(1, 1L, 6L)
  )
})

test_that("catch and index observations retain supplied values and scales", {
  catch <- matrix(seq_len(6L), 6L, 1L)
  catch_paa <- array(rep(c(0.1, 0.2, 0.3, 0.4), each = 6L), c(1L, 6L, 4L))
  indices <- matrix(seq(10, 60, by = 10), 6L, 1L)
  index_paa <- array(rep(c(0.4, 0.3, 0.2, 0.1), each = 6L), c(1L, 6L, 4L))

  input <- small_input(
    catch_info = list(
      fleet_names = "test fleet",
      agg_catch = catch,
      agg_catch_cv = matrix(0.2, 6L, 1L),
      catch_paa = catch_paa,
      catch_Neff = matrix(25, 6L, 1L)
    ),
    index_info = list(
      index_names = "test index",
      agg_indices = indices,
      agg_index_cv = matrix(0.3, 6L, 1L),
      index_paa = index_paa,
      index_Neff = matrix(30, 6L, 1L)
    )
  )

  expect_identical(input$fleet_names, "test fleet")
  expect_identical(input$index_names, "test index")
  expect_equal(input$data$agg_catch, catch)
  expect_equal(input$data$catch_paa, catch_paa)
  expect_equal(input$data$catch_Neff, matrix(25, 6L, 1L))
  expect_equal(input$data$agg_catch_sigma, matrix(sqrt(log(1 + 0.2^2)), 6L, 1L))
  expect_equal(input$data$agg_indices, indices)
  expect_equal(input$data$index_paa, index_paa)
  expect_equal(input$data$index_Neff, matrix(30, 6L, 1L))
  expect_equal(input$data$agg_index_sigma, matrix(sqrt(log(1 + 0.3^2)), 6L, 1L))
})

test_that("fully selected fishing mortality supports both parameterizations", {
  input <- small_input()
  annual_F <- matrix(c(0.1, 0.2, 0.4, 0.8, 0.6, 0.3), 6L, 1L)
  map_F <- matrix(c(NA_integer_, seq_len(5L)), 6L, 1L)

  deviations <- set_F(input, list(F = annual_F, F_config = 1L, map_F = map_F))
  annual_values <- set_F(input, list(F = annual_F, F_config = 2L))

  expect_identical(deviations$data$F_config, 1L)
  expect_equal(cumsum(deviations$par$F_pars[, 1L]), log(annual_F[, 1L]))
  expect_identical(annual_values$data$F_config, 2L)
  expect_equal(annual_values$par$F_pars, log(annual_F))
  expect_true(is.factor(deviations$map$F_pars))
  expect_equal(as.integer(deviations$map$F_pars), as.integer(factor(map_F)))
})

test_that("basic years and ages must be integer vectors", {
  expect_error(
    suppressMessages(prepare_wham_input(
      basic_info = list(years = as.numeric(2000:2005), ages = 1:4)
    )),
    "basic_info\\$years"
  )
  expect_error(
    suppressMessages(prepare_wham_input(
      basic_info = list(years = 2000:2005, ages = as.numeric(1:4))
    )),
    "basic_info\\$ages"
  )
})
