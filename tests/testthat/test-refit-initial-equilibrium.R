test_that("RT refits preserve historical process scales and likelihood", {
  path <- test_path("..", "..", "oldfitobjects", "BSB_2023RTfit.RDS")
  skip_if_not(file.exists(path), "Local BSB RT fit is unavailable")
  fit <- readRDS(path)
  original_input <- fit$input
  fit$call <- quote(fit_wham(
    input = input,
    do.retro = FALSE,
    do.osa = FALSE,
    do.sdrep = FALSE,
    do.brps = FALSE,
    MakeADFun.silent = TRUE
  ))
  capture.output(result <- suppressMessages(refit_wham(fit)))
  expect_equal(
    result$opt$objective,
    fit$opt$objective,
    tolerance = 1e-6,
    scale = 1
  )
  expect_equal(result$input$data$N1_model, c(1, 1))
  expect_equal(result$rep$N1, fit$rep$N1, tolerance = 1e-6)
  expect_equal(result$rep$NAA, fit$rep$NAA, tolerance = 1e-6)
  expect_equal(
    result$parList$sel_repars[, 2:3] / 2,
    fit$parList$sel_repars[, 2:3],
    tolerance = 1e-6
  )
  rho <- tanh(result$parList$Ecov_process_pars[3, ] / 2)
  expect_equal(
    result$parList$Ecov_process_pars[2, ] - 0.5 * log1p(-rho^2),
    fit$parList$Ecov_process_pars[2, ],
    tolerance = 1e-6
  )
  expect_identical(fit$input, original_input)
  scales <- refit_process_scales(result$input, result$rep)
  expect_identical(scales$input, result$input)
  expect_identical(scales$selectivity, FALSE)
  expect_length(scales$ecov, 0)
})

test_that("equilibrium refits with selectivity random effects use the correction", {
  path <- test_path(
    "..",
    "..",
    "oldfitobjects",
    "sel_re_logistic_normal_fit_5.rds"
  )
  skip_if_not(file.exists(path), "Local selectivity fit is unavailable")
  fit <- readRDS(path)
  fit$call <- quote(fit_wham(
    input = input,
    do.retro = FALSE,
    do.osa = FALSE,
    do.sdrep = FALSE,
    do.brps = FALSE,
    MakeADFun.silent = TRUE
  ))
  warnings <- character()
  capture.output(
    result <- withCallingHandlers(
      suppressMessages(refit_wham(fit)),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    )
  )
  expect_match(paste(warnings, collapse = "\n"), "historical equilibrium bug")
  expect_equal(result$input$data$N1_model, 1L)
  expect_equal(result$opt$convergence, 0)
  input <- result$input
  input$par <- result$parList
  expect_identical(refit_initial_equilibrium(input, result$rep$N1), input)
  local_edition(3)
  expect_snapshot({
    expect_identical(refit_initial_equilibrium(input, result$rep$N1 * 2), input)
    expect_identical(
      refit_initial_equilibrium(input, array(0, c(1L, 1L, 1L))),
      input
    )
  })
})

test_that("Run4 refits warn and use corrected equilibrium", {
  path <- test_path("..", "..", "oldfitobjects", "BSB_2023Run4fit.RDS")
  skip_if_not(file.exists(path), "Local BSB Run4 fit is unavailable")
  fit <- readRDS(path)
  original_input <- fit$input
  fit$call <- quote(fit_wham(
    input = input,
    do.retro = FALSE,
    do.osa = FALSE,
    do.sdrep = FALSE,
    do.brps = FALSE,
    MakeADFun.silent = TRUE
  ))
  warnings <- character()
  capture.output(
    result <- withCallingHandlers(
      suppressMessages(refit_wham(fit)),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    )
  )
  expect_match(paste(warnings, collapse = "\n"), "historical equilibrium bug")
  expect_equal(result$input$data$N1_model, 1L)
  expect_equal(result$opt$convergence, 0)
  expect_equal(is.finite(result$opt$objective), TRUE)
  expect_identical(fit$input, original_input)

  input <- result$input
  input$par <- result$parList
  expect_identical(refit_initial_equilibrium(input, result$rep$N1), input)
  model <- TMB::MakeADFun(
    input$data,
    input$par,
    map = input$map,
    DLL = "wham",
    silent = TRUE
  )
  on.exit(TMB::FreeADFun(model))
  model$fn(model$par)
  current_N1 <- model$report()$N1
  expect_equal(result$rep$N1, current_N1, tolerance = 1e-8)
  expect_identical(refit_initial_equilibrium(input, current_N1), input)
})

test_that("equilibrium warnings handle saved compatibility modes and missing reports", {
  local_edition(3)
  input <- list(data = list(N1_model = c(0L, -1L, 2L)))
  corrected <- input
  corrected$data$N1_model[2L] <- 1L
  expect_snapshot({
    expect_identical(refit_initial_equilibrium(input, NULL), corrected)
    expect_identical(refit_initial_equilibrium(corrected, NULL), corrected)
  })
  input$data$N1_model <- c(0L, 2L)
  expect_identical(refit_initial_equilibrium(input, NULL), input)
})
