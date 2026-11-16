# btime <- Sys.time(); devtools::test(filter = "jitter_wham"); etime <- Sys.time(); runtime = etime - btime; runtime;

context("Jittered model fits")

test_that("jitter_wham validates the model, covariance, and requested rows", {
  expect_error(jitter_wham(), "Provide fit_RDS", fixed = TRUE)
  fit_path <- tempfile(fileext = ".rds")
  on.exit(unlink(fit_path), add = TRUE)
  saveRDS(list(), fit_path)
  expect_error(
    suppressWarnings(suppressMessages(jitter_wham(fit_path, do_parallel = FALSE))),
    "must have do_sdreport() completed", fixed = TRUE
  )
  expect_error(
    suppressWarnings(suppressMessages(
      jitter_wham(fit_path, initial_vals = matrix(1:4, 2),
                  which_rows = 3, do_parallel = FALSE)
    )),
    "some of which_rows are outside the rows of initial_vals.", fixed = TRUE
  )
})
test_that("serial jitters refit supplied starting values and persist results", {
  fit_path <- tempfile(fileext = ".rds")
  result_dir <- tempfile("wham-jitter-")
  dir.create(result_dir)
  on.exit(unlink(c(fit_path, result_dir), recursive = TRUE), add = TRUE)
  saveRDS(list(par = c(a = 0, b = 0),
               env = list(data = list(do_SPR_BRPs = 1),
                          last.par.best = c(0, 0, 0))), fit_path)
  initial <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  calls <- list()
  testthat::local_mocked_bindings(
    fit_tmb = function(model, n.newton, do.sdrep) {
      calls[[length(calls) + 1L]] <<- model
      expect_equal(n.newton, 0)
      expect_false(do.sdrep)
      model$opt <- list(obj = sum(model$par^2), par = model$par + 1)
      model$gr <- function(par) 2 * par
      model$rep <- list(SSB = 10)
      model$parList <- list(beta = model$opt$par)
      model$env$last.par.best <- c(model$opt$par, 5)
      model
    },
    .package = "wham"
  )
  result <- suppressWarnings(suppressMessages(
    jitter_wham(fit_path, n_jitter = 2, initial_vals = initial,
                do_parallel = FALSE, res_dir = result_dir)
  ))
  expect_equal(result$initial_values, initial)
  expect_equal(readRDS(file.path(result_dir, "initial_values.RDS")), initial)
  expect_length(calls, 2)
  for (i in 1:2) {
    expect_equal(unname(calls[[i]]$par), initial[i, ])
    expect_equal(calls[[i]]$env$data$do_SPR_BRPs, 0)
    out <- result$jitter_results[[i]]
    expect_equal(out$initial_vals_row, i)
    expect_equal(out$initial_vals, initial[i, ])
    expect_equal(out$obj, sum(initial[i, ]^2))
    expect_equal(unname(out$par), initial[i, ] + 1)
    expect_equal(out$grad, 2 * out$par)
    expect_equal(out$rep, list(SSB = 10))
    expect_equal(out$parList, list(beta = out$par))
    expect_equal(out$last.par.best, c(out$par, 5))
    expect_equal(readRDS(file.path(result_dir, paste0("jitter_sim_", i, ".RDS"))),
                 out)
  }
})
