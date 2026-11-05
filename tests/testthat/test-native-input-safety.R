native_safety_input <- function(process = "ar1", sigma = "est_re") {
  suppressMessages(prepare_wham_input(
    basic_info = list(years = 2000:2005, ages = 1:4),
    NAA_re = list(sigma = "rec"),
    ecov = list(
      label = "temperature",
      year = 2001:2003,
      mean = matrix(c(-0.2, 0.1, 0.3), 3, 1),
      use_obs = matrix(1, 3, 1),
      logsigma = sigma,
      process_model = process,
      recruitment_how = matrix("controlling-lag-1-linear", 1, 1),
      M_how = array("lag-0-linear", c(1, 1, 4, 1)),
      q_how = matrix("lag-2-linear", 1, 1)
    )
  ))
}

native_safety_objective <- function(input, mapped = TRUE) {
  if (!mapped) input$map <- list()
  validate_input_safety(input)
  TMB::MakeADFun(
    data = input$data,
    parameters = input$par,
    map = if (mapped) input$map else NULL,
    random = if (mapped) input$random else NULL,
    DLL = "wham",
    silent = TRUE
  )
}

test_that("R validation rejects historical NAA matrix ranks safely", {
  for (field in c("log_NAA_sigma", "trans_NAA_rho")) {
    input <- native_safety_input()
    input$par[[field]] <- matrix(input$par[[field]], nrow = 1L)
    objective <- NULL
    error <- tryCatch(
      {
        objective <- native_safety_objective(input, mapped = FALSE)
        objective$fn(objective$par)
      },
      error = identity,
      finally = {
        if (!is.null(objective)) TMB::FreeADFun(objective)
      }
    )
    expect_s3_class(error, "error")
    expect_match(
      conditionMessage(error),
      paste0(field, " must have dimensions"),
      fixed = TRUE
    )
  }
})

test_that("R validation guards malformed Ecov dimensions and endpoints", {
  input <- native_safety_input()
  malformed <- list(
    function(x) {
      x$par$Ecov_re <- matrix(0, nrow(x$par$Ecov_re), 0)
      x
    },
    function(x) {
      x$par$Ecov_process_pars <- matrix(0, 2, 1)
      x
    },
    function(x) {
      x$par$Ecov_obs_logsigma_re <- matrix(0, 1, 1)
      x
    },
    function(x) {
      x$par$Ecov_obs_sigma_par <- matrix(0, 1, 1)
      x
    },
    function(x) {
      x$data$ind_Ecov_out_end_R <- x$data$ind_Ecov_out_end_R - 1L
      x
    },
    function(x) {
      x$data$years_use_Ecov[1] <- x$data$n_years_Ecov
      x
    }
  )
  patterns <- c(
    "Ecov_re must span",
    "Ecov_process_pars must have dimensions",
    "Ecov_obs_logsigma_re must have dimensions",
    "Ecov_obs_sigma_par must have dimensions",
    "start/end indices",
    "years_use_Ecov contains invalid"
  )
  for (i in seq_along(malformed)) {
    objective <- NULL
    error <- tryCatch(
      {
        objective <- native_safety_objective(
          malformed[[i]](input),
          mapped = FALSE
        )
        objective$fn(objective$par)
      },
      error = identity,
      finally = {
        if (!is.null(objective)) TMB::FreeADFun(objective)
      }
    )
    expect_s3_class(error, "error")
    expect_match(conditionMessage(error), patterns[i], fixed = TRUE)
  }
})

test_that("native objectives and gradients remain finite after safe peels", {
  for (process in c("ar1", "rw")) {
    for (sigma in list(-1.7, "est_re")) {
      input <- native_safety_input(process, sigma)
      input$data$selblock_pointer_fleets[, 1] <- c(1, 1, 1, 1, 3, 3)
      input <- set_selectivity(
        input,
        list(
          n_selblocks = 3,
          model = rep("logistic", 3),
          re = c("none", "none", "ar1_y")
        )
      )
      peeled <- reduce_input(input, 2004:2005)
      expect_equal(peeled$data$n_years_selblocks[3], 0)
      expect_equal(peeled$data$selblock_models_re[3], 0)
      for (candidate in list(input, peeled)) {
        objective <- native_safety_objective(candidate)
        value <- objective$fn(objective$par)
        gradient <- objective$gr(objective$par)
        expect_equal(all(is.finite(value)), TRUE)
        expect_equal(all(is.finite(gradient)), TRUE)
        TMB::FreeADFun(objective)
      }
    }
  }
})
