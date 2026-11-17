# btime <- Sys.time(); devtools::test(filter = "self_test"); etime <- Sys.time(); runtime = etime - btime; runtime;

context("Self-test simulations")
test_that("self_test requires a model and valid seed indices", {
  expect_error(self_test(), "Provide fit_RDS", fixed = TRUE)
  fit_path <- tempfile(fileext = ".rds")
  on.exit(unlink(fit_path), add = TRUE)
  saveRDS(list(), fit_path)
  expect_error(
    suppressWarnings(suppressMessages(
      self_test(fit_path, seeds = c(11, 22), which_seeds = 3,
                do_parallel = FALSE)
    )),
    "some of which_seeds are outside 1 to n.", fixed = TRUE
  )
})
test_that("serial self tests simulate, refit, and save conditional and unconditional inputs", {
  fit_path <- tempfile(fileext = ".rds")
  result_dir <- tempfile("wham-self-test-")
  dir.create(result_dir)
  on.exit(unlink(c(fit_path, result_dir), recursive = TRUE), add = TRUE)
  fit <- list(
    input = list(
      data = list(do_SPR_BRPs = 1, do_simulate_N_re = 1,
                  do_simulate_M_re = 1, do_simulate_data = c(1, 1, 1)),
      par = list(beta = -1),
      map = list(beta = factor(1)),
      random = "N_re"
    ),
    parList = list(beta = 2)
  )
  saveRDS(fit, fit_path)
  report <- list(SSB = matrix(c(10, 20), 2), log_F_tot = log(c(0.1, 0.2)),
                 NAA = array(1:4, c(1, 1, 2, 2)))
  calls <- list()
  testthat::local_mocked_bindings(
    fit_wham = function(input, do.fit = TRUE, ...) {
      calls[[length(calls) + 1L]] <<- list(input = input, do.fit = do.fit,
                                         options = list(...))
      if (!do.fit) {
        return(list(
          par = c(beta = 2), rep = report,
          simulate = function(complete) {
            expect_true(complete)
            data <- input$data
            data$simulated <- runif(1)
            data
          }
        ))
      }
      list(opt = list(obj = 12, par = c(beta = 3)), final_gradient = 0.001,
           rep = report)
    },
    .package = "wham"
  )
  new_map <- list(beta = factor(NA))
  for (conditional in c(TRUE, FALSE)) {
    calls <- list()
    result <- suppressWarnings(suppressMessages(
      self_test(fit_path, n = 2, seeds = c(11, 22), conditional = conditional,
                map_change = new_map, do_parallel = FALSE,
                res_dir = result_dir, save_inputs = TRUE)
    ))
    expect_identical(result$seeds, c(11, 22))
    expect_length(result$self_test_results, 2)
    expect_length(calls, 4)
    for (i in 1:2) {
      simulation <- calls[[2 * i - 1]]
      estimation <- calls[[2 * i]]
      expect_false(simulation$do.fit)
      expect_null(simulation$input$random)
      expect_equal(simulation$input$par, fit$parList)
      expect_equal(simulation$input$map, new_map)
      expect_equal(simulation$input$data$do_SPR_BRPs, 0)
      expect_equal(simulation$input$data$do_simulate_N_re, as.numeric(!conditional))
      expect_equal(simulation$input$data$do_simulate_M_re, as.numeric(!conditional))
      expect_equal(simulation$input$data$do_simulate_data, c(1, 1, 1))
      expect_true(estimation$do.fit)
      expect_identical(estimation$input$random, "N_re")
      expect_false(estimation$options$do.osa)
      expect_false(estimation$options$do.retro)
      expect_false(estimation$options$do.sdrep)
      out <- result$self_test_results[[i]]
      expect_equal(out$obj, 12)
      expect_equal(out$par, c(beta = 3))
      expect_equal(out$grad, 0.001)
      expect_equal(out$SSB, report$SSB)
      expect_equal(out$F, c(0.1, 0.2))
      expect_equal(out$NAA, report$NAA)
      expect_equal(out$seed, result$seeds[i])
      prefix <- if (conditional) "cond_sim" else "sim"
      expect_equal(readRDS(file.path(result_dir, paste0(prefix, "_", i, ".RDS"))),
                   out)
      saved_input <- readRDS(file.path(result_dir, paste0(prefix, "_input_", i, ".RDS")))
      expect_equal(saved_input, estimation$input)
    }
    expect_false(identical(calls[[2]]$input$data$simulated,
                           calls[[4]]$input$data$simulated))
  }
})
