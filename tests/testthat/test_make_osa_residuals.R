# btime <- Sys.time(); devtools::test(filter = "make_osa"); etime <- Sys.time(); runtime = etime - btime; runtime;

context("One-step-ahead residuals")
test_that("OSA validates methods, sdreport, and unsupported age compositions", {
  testthat::local_mocked_bindings(verify_version = function(model) NULL,
                                 .package = "wham")
  model <- list(
    env = list(data = list(do_SPR_BRPs = 1, do_MSY_BRPs = 1)),
    is_sdrep = TRUE,
    input = list(data = list(age_comp_model_fleets = 1,
                             age_comp_model_indices = 1))
  )
  expect_error(make_osa_residuals(model, list(method = "invalid")),
               "Only osa methods allowed", fixed = TRUE)
  model$is_sdrep <- FALSE
  expect_error(make_osa_residuals(model), "TMB::sdreport completed", fixed = TRUE)
  model$is_sdrep <- TRUE
  model$input$data$age_comp_model_fleets <- 10
  expect_error(make_osa_residuals(model), "mvtweedie", fixed = TRUE)
  model$input$data$age_comp_model_fleets <- 1
  model$input$data$age_comp_model_indices <- 10
  expect_error(make_osa_residuals(model), "mvtweedie", fixed = TRUE)
})
test_that("OSA conditions on earlier groups and masks terminal multinomial ages", {
  testthat::local_mocked_bindings(verify_version = function(model) NULL,
                                 .package = "wham")
  obs <- data.frame(
    type = c("catchpaa", "logcatch", "Ecov", "indexpaa", "logindex",
             "catchpaa", "indexpaa", "other"),
    age = c(1, NA, NA, 1, NA, 2, 2, NA),
    fleet = c("fleet_1", "fleet_1", NA, "index_1", "index_1",
              "fleet_1", "index_1", NA),
    val = seq_len(8)
  )
  model <- list(
    env = list(data = list(do_SPR_BRPs = 1, do_MSY_BRPs = 2)),
    is_sdrep = FALSE,
    input = list(data = list(
      obs = as.list(obs), age_comp_model_fleets = 1,
      age_comp_model_indices = 11, condition_no_osa = 8L
    ))
  )
  calls <- list()
  testthat::local_mocked_bindings(
    oneStepPredict = function(obj, method, subset, conditional, observation.name,
                              data.term.indicator, discrete, parallel) {
      calls[[length(calls) + 1L]] <<- list(
        subset = subset, conditional = conditional, method = method,
        parallel = parallel
      )
      expect_equal(obj$env$data$do_SPR_BRPs, 0)
      expect_equal(obj$env$data$do_MSY_BRPs, 0)
      expect_identical(observation.name, "obsvec")
      expect_identical(data.term.indicator, "keep")
      expect_false(discrete)
      data.frame(residual = subset / 10)
    },
    .package = "TMB"
  )
  for (method in c("oneStepGaussianOffMode", "oneStepGaussian", "oneStepGeneric")) {
    calls <- list()
    opts <- if (method == "oneStepGaussianOffMode") list() else
      list(method = method, parallel = FALSE)
    result <- suppressMessages(make_osa_residuals(model, opts, sdrep_required = FALSE))
    expect_length(calls, 3)
    expect_identical(calls[[1]]$subset, 3L)
    expect_null(calls[[1]]$conditional)
    expect_identical(calls[[2]]$subset, c(2L, 5L))
    expect_identical(calls[[2]]$conditional, 3L)
    expect_identical(calls[[3]]$subset, c(1L, 4L, 6L, 7L))
    expect_identical(calls[[3]]$conditional, c(3L, 2L, 5L, 8L))
    for (call in calls) {
      expect_identical(call$method, method)
      expect_identical(call$parallel, method == "oneStepGaussianOffMode")
    }
    expect_equal(result$osa$residual, c(0.1, 0.2, 0.3, 0.4, 0.5, NA, NA, NA))
    expect_equal(result$osa[names(obs)], obs)
    expect_equal(result$OSA.Ecov$residual, 0.3)
    expect_equal(result$OSA.aggregate$residual, c(0.2, 0.5))
    expect_equal(result$OSA.agecomp$residual, c(0.1, 0.4, 0.6, 0.7))
    expect_equal(result$env$data, model$env$data)
  }
})
