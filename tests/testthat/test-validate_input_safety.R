validation_safety_input <- function() {
  suppressMessages(prepare_wham_input(
    NAA_re = list(sigma = "rec"),
    basic_info = list(years = 2000:2005, ages = 1:4),
    ecov = list(
      label = "temperature",
      year = 2000:2005,
      mean = matrix(1:6, 6, 1),
      use_obs = matrix(1, 6, 1),
      logsigma = "est_re",
      process_model = "ar1",
      recruitment_how = matrix("controlling-lag-1-linear", 1, 1)
    )
  ))
}

test_that("prepared inputs pass safety validation without modification", {
  input <- validation_safety_input()
  result <- withVisible(validate_input_safety(input))
  expect_identical(result$value, input)
  expect_identical(result$visible, FALSE)
  no_ecov <- suppressMessages(prepare_wham_input(
    basic_info = list(years = 2000:2005, ages = 1:4)
  ))
  expect_identical(validate_input_safety(no_ecov), no_ecov)
})

test_that("R guards cover vector lengths, random-walk predecessors and projections", {
  input <- validation_safety_input()
  broken <- input
  broken$data$NAA_re_model <- integer()
  local_edition(3)
  expect_snapshot(error = TRUE, validate_input_safety(broken))
  broken <- input
  broken$data$selblock_models <- integer()
  expect_snapshot(error = TRUE, validate_selpars_re(broken))
  broken <- input
  broken$data$Ecov_model[] <- 1L
  broken$data$years_use_Ecov <- c(0L, 1L, 0L)
  expect_snapshot(error = TRUE, validate_input_safety(broken))
  broken <- input
  broken$data$Ecov_model[] <- 3L
  expect_snapshot(error = TRUE, validate_input_safety(broken))
  input$data$n_years_proj <- 2L
  input$data$proj_Ecov_opt <- 3L
  expect_identical(validate_input_safety(input), input)
  broken <- input
  broken$data$n_years_proj <- -1L
  expect_snapshot(error = TRUE, validate_input_safety(broken))
  broken <- input
  broken$data$proj_Ecov_opt <- 0L
  expect_snapshot(error = TRUE, validate_input_safety(broken))
  broken$data$proj_Ecov_opt <- 1L
  expect_snapshot(error = TRUE, validate_input_safety(broken))
  broken$data$proj_Ecov_opt <- 2L
  broken$data$avg_years_Ecov <- integer()
  expect_snapshot(error = TRUE, validate_input_safety(broken))
  broken$data$proj_Ecov_opt <- 4L
  broken$data$Ecov_use_proj <- matrix(0, 1L, 1L)
  expect_snapshot(error = TRUE, validate_input_safety(broken))
  broken$data$Ecov_use_proj <- matrix(0, 2L, 1L)
  expect_identical(validate_input_safety(broken), broken)
})

test_that("unsafe environmental dimensions, mappings and values identify their fields", {
  input <- validation_safety_input()
  cases <- list(
    list("n_years_model", function(x) {
      x$data$n_years_model <- 0
      x
    }),
    list("n_Ecov", function(x) {
      x$data$n_Ecov <- NA_real_
      x
    }),
    list("n_years_Ecov", function(x) {
      x$data$n_years_Ecov <- 1.5
      x
    }),
    list("at least two model years", function(x) {
      x$data$n_years_model <- 1
      x
    }),
    list("map$Ecov_re", function(x) {
      x$map$Ecov_re <- 1:7
      x
    }),
    list("map$Ecov_re", function(x) {
      x$map$Ecov_re <- factor(1)
      x
    }),
    list("random contains unknown", function(x) {
      x$random <- "not_a_parameter"
      x
    }),
    list("Ecov_model", function(x) {
      x$data$Ecov_model <- NA_real_
      x
    }),
    list("Ecov_use_re", function(x) {
      x$data$Ecov_use_re <- numeric()
      x
    }),
    list("Ecov_obs_sigma_opt", function(x) {
      x$data$Ecov_obs_sigma_opt <- 1.5
      x
    }),
    list("Ecov_process_pars", function(x) {
      x$par$Ecov_process_pars <- matrix(0, 1, 3)
      x$map$Ecov_process_pars <- NULL
      x
    }),
    list("Ecov_re must span", function(x) {
      x$par$Ecov_re <- matrix(0, 1, 1)
      x$map$Ecov_re <- NULL
      x
    }),
    list("Ecov_obs must have dimensions", function(x) {
      x$data$Ecov_obs <- matrix(0, 1, 1)
      x
    }),
    list("Ecov_use_obs must have dimensions", function(x) {
      x$data$Ecov_use_obs <- matrix(0, 1, 1)
      x
    }),
    list("Ecov_obs_logsigma must have dimensions", function(x) {
      x$par$Ecov_obs_logsigma <- matrix(0, 1, 1)
      x$map$Ecov_obs_logsigma <- NULL
      x
    }),
    list("Ecov_obs_logsigma_re must have dimensions", function(x) {
      x$par$Ecov_obs_logsigma_re <- matrix(0, 1, 1)
      x$map$Ecov_obs_logsigma_re <- NULL
      x
    }),
    list("Ecov_obs_sigma_par must have dimensions", function(x) {
      x$par$Ecov_obs_sigma_par <- matrix(0, 1, 2)
      x$map$Ecov_obs_sigma_par <- NULL
      x
    }),
    list("finite starting values", function(x) {
      x$par$Ecov_obs_logsigma_re[1, 1] <- Inf
      x
    }),
    list("finite starting values", function(x) {
      x$par$Ecov_obs_sigma_par[1, 1] <- NA_real_
      x
    }),
    list("used environmental observations must be finite", function(x) {
      x$data$Ecov_obs[2, 1] <- NA_real_
      x
    }),
    list("finite log standard deviations", function(x) {
      x$data$Ecov_obs_sigma_opt <- 1
      x$par$Ecov_obs_logsigma[2, 1] <- Inf
      x
    }),
    list("years_use_Ecov", function(x) {
      x$data$years_use_Ecov <- -1
      x
    }),
    list("years_use_Ecov", function(x) {
      x$data$years_use_Ecov <- nrow(x$par$Ecov_re)
      x
    }),
    list("years_use_Ecov", function(x) {
      x$data$years_use_Ecov <- numeric()
      x
    })
  )
  for (case in cases) {
    error <- tryCatch(validate_input_safety(case[[2]](input)), error = identity)
    expect_s3_class(error, "error")
    expect_match(conditionMessage(error), "Unsafe WHAM input:", fixed = TRUE)
    expect_match(conditionMessage(error), case[[1]], fixed = TRUE)
  }
  input$data$Ecov_obs[2, 1] <- NA_real_
  input$data$Ecov_use_obs[2, 1] <- 0
  expect_identical(validate_input_safety(input), input)
})

test_that("active environmental families validate every lag endpoint and degree", {
  for (effect in c("R", "M", "q", "mu")) {
    input <- validation_safety_input()
    how <- paste0("Ecov_how_", effect)
    start <- paste0("ind_Ecov_out_start_", effect)
    end <- paste0("ind_Ecov_out_end_", effect)
    degree <- paste0("n_poly_Ecov_", effect)
    # Give movement a real destination dimension, without invoking native code.
    if (effect == "mu") {
      for (field in c(how, start, end, degree)) {
        input$data[[field]] <- array(1, c(1, 1, 4, 1, 2, 1))
      }
    }
    input$data[[how]][] <- 1
    input$data[[start]][] <- 1
    input$data[[end]][] <- 6
    input$data[[degree]][] <- 1
    expect_identical(validate_input_safety(input), input)
    for (field in c(start, end, degree)) {
      broken <- input
      broken$data[[field]] <- as.vector(broken$data[[field]])
      error <- tryCatch(validate_input_safety(broken), error = identity)
      expect_match(
        conditionMessage(error),
        paste0("invalid Ecov ", effect),
        fixed = TRUE
      )
    }
    for (value in c(-1, 2)) {
      broken <- input
      broken$data[[start]][1] <- value
      error <- tryCatch(validate_input_safety(broken), error = identity)
      expect_match(conditionMessage(error), "start/end indices", fixed = TRUE)
    }
    for (value in c(0, 6)) {
      broken <- input
      broken$data[[degree]][1] <- value
      error <- tryCatch(validate_input_safety(broken), error = identity)
      expect_match(conditionMessage(error), "polynomial degree", fixed = TRUE)
    }
  }
})

test_that("native environmental option codes and observation flags are bounded", {
  input <- validation_safety_input()
  for (field in c(
    "Ecov_model",
    "Ecov_use_re",
    "Ecov_obs_sigma_opt",
    "Ecov_use_obs"
  )) {
    broken <- input
    broken$data[[field]][1] <- 99
    error <- tryCatch(validate_input_safety(broken), error = identity)
    expect_s3_class(error, "error")
    if (inherits(error, "error")) {
      expect_match(conditionMessage(error), "Unsafe WHAM input:", fixed = TRUE)
      expect_match(conditionMessage(error), field, fixed = TRUE)
    }
  }
  input$data$Ecov_how_R[1] <- NA_real_
  error <- tryCatch(validate_input_safety(input), error = identity)
  expect_s3_class(error, "error")
  expect_match(conditionMessage(error), "Unsafe WHAM input:", fixed = TRUE)
  expect_match(conditionMessage(error), "Ecov_how_R", fixed = TRUE)
})

test_that("selectivity guards reject inconsistent counts and empty active blocks", {
  input <- validation_safety_input()
  cases <- list(
    list("n_selblocks", function(x) {
      x$data$n_selblocks <- 0
      x
    }),
    list("selblock_years", function(x) {
      x$data$selblock_years <- matrix(0, 1, 1)
      x
    }),
    list("column sums", function(x) {
      x$data$n_years_selblocks[1] <- 0
      x
    }),
    list("selpars_re", function(x) {
      x$par$selpars_re <- array(0, c(1, 1, 1))
      x$map$selpars_re <- NULL
      x
    }),
    list("sel_repars", function(x) {
      x$par$sel_repars <- matrix(0, 1, 1)
      x$map$sel_repars <- NULL
      x
    }),
    list("random-effect counts", function(x) {
      x$data$selblock_models_re[1] <- 1
      x$data$n_selpars_re[1] <- 0
      x
    }),
    list("no retained years", function(x) {
      x$data$selblock_years[, 1] <- 0
      x$data$n_years_selblocks[1] <- 0
      x$data$selblock_models_re[1] <- 1
      x$data$n_selpars_re[1] <- 1
      x
    })
  )
  for (case in cases) {
    error <- tryCatch(validate_input_safety(case[[2]](input)), error = identity)
    expect_s3_class(error, "error")
    expect_match(conditionMessage(error), case[[1]], fixed = TRUE)
  }
})

test_that("fit_peel guards reject invalid requests before native model construction", {
  input <- validation_safety_input()
  for (peel in list(
    NULL,
    numeric(),
    c(1, 2),
    NA_real_,
    Inf,
    "1",
    0,
    -1,
    1.5,
    5,
    6
  )) {
    error <- tryCatch(
      fit_peel(peel, input, retro.silent = TRUE),
      error = identity
    )
    expect_s3_class(error, "error")
    expect_match(
      conditionMessage(error),
      "peel must be a positive integer retaining at least two model years",
      fixed = TRUE
    )
  }
  input$data$Ecov_obs[2, 1] <- Inf
  error <- tryCatch(fit_peel(1, input, retro.silent = TRUE), error = identity)
  expect_match(
    conditionMessage(error),
    "Unsafe WHAM input: used environmental observations",
    fixed = TRUE
  )
  input <- validation_safety_input()
  input$data$n_poly_Ecov_R[] <- 4
  error <- tryCatch(fit_peel(2, input, retro.silent = TRUE), error = identity)
  expect_match(conditionMessage(error), "polynomial degree", fixed = TRUE)
})

test_that("retro catches each failed peel and continues without native fitting", {
  input <- validation_safety_input()
  input$data$Ecov_obs[2, 1] <- Inf
  model <- list(input = input, parList = input$par)
  warnings <- character()
  peels <- withCallingHandlers(
    retro(model, n.peels = 3, check.version = FALSE, retro.silent = TRUE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(peels, 3)
  expect_length(warnings, 3)
  for (i in 1:3) {
    expect_equal(peels[[i]]$peel, i)
    expect_match(
      peels[[i]]$err,
      paste0("Retrospective peel ", i, ": Unsafe WHAM input:"),
      fixed = TRUE
    )
    expect_identical(warnings[i], peels[[i]]$err)
  }
})

test_that("retro forwards MLEs and random overrides and continues after a middle failure", {
  input <- validation_safety_input()
  mle <- input$par
  mle$Ecov_re[] <- 4.2
  model <- list(input = input, parList = mle)
  # Only the cloned closure sees the R-only fitter; package/global bindings stay intact.
  local_retro <- retro
  scope <- new.env(parent = environment(retro))
  scope$fit_peel <- function(peel, input, ...) {
    if (peel == 2) {
      stop("deliberate R-only peel failure")
    }
    list(peel = peel, input = reduce_input(input, tail(input$years, peel)))
  }
  environment(local_retro) <- scope
  for (use_mle in c(TRUE, FALSE)) {
    warnings <- character()
    result <- withCallingHandlers(
      local_retro(
        model,
        3,
        ran = "F_pars",
        use.mle = use_mle,
        check.version = FALSE
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_length(result, 3)
    expect_identical(
      warnings,
      "Retrospective peel 2: deliberate R-only peel failure"
    )
    expect_equal(result[[2]]$peel, 2)
    expect_equal(result[[3]]$peel, 3)
    for (i in c(1, 3)) {
      expect_identical(result[[i]]$input$random, "F_pars")
      source_par <- if (use_mle) mle else input$par
      expect_equal(
        result[[i]]$input$par$Ecov_re,
        source_par$Ecov_re[
          seq_len(nrow(source_par$Ecov_re) - i),
          ,
          drop = FALSE
        ]
      )
    }
  }
})

test_that("projection padding keeps fixed sigma inactive and gives sigma REs new maps", {
  years <- 2000:2005
  fixed_sigma <- matrix(seq(-2.2, -1.1, length.out = length(years)), ncol = 1L)
  input <- suppressMessages(prepare_wham_input(
    NAA_re = list(sigma = "rec"),
    basic_info = list(years = years, ages = 1:4),
    ecov = list(
      label = c("fixed", "estimated"),
      year = years,
      mean = cbind(seq_along(years), rev(seq_along(years))),
      use_obs = matrix(1, length(years), 2L),
      logsigma = list(cbind(fixed_sigma, NA_real_), c(NA, "est_re")),
      process_model = c("ar1", "ar1"),
      recruitment_how = matrix("controlling-lag-1-linear", 2L, 1L)
    )
  ))
  model <- list(
    input = input,
    parList = input$par,
    years = years,
    rep = list(FAA = array(1, c(1L, length(years), 4L)))
  )

  projected <- suppressMessages(prepare_projection(
    model,
    list(
      n.yrs = 3L,
      use.last.F = TRUE,
      use.avg.F = FALSE,
      use.FXSPR = FALSE,
      use.FMSY = FALSE,
      cont.ecov = TRUE,
      use.last.ecov = FALSE
    ),
    check.version = FALSE
  ))

  old_map <- matrix(as.integer(input$map$Ecov_obs_logsigma_re), ncol = 2L)
  new_map <- matrix(as.integer(projected$map$Ecov_obs_logsigma_re), ncol = 2L)
  padded_rows <- nrow(old_map) + 1:3
  expect_equal(
    dim(projected$par$Ecov_obs_logsigma_re),
    c(nrow(old_map) + 3L, 2L)
  )
  expect_identical(new_map[seq_len(nrow(old_map)), , drop = FALSE], old_map)
  expect_true(all(is.na(new_map[, 1L])))
  expect_true(all(is.na(new_map[padded_rows, 1L])))
  expect_equal(new_map[padded_rows, 2L], max(old_map[, 2L], na.rm = TRUE) + 1:3)
  expect_equal(length(unique(new_map[padded_rows, 2L])), 3L)
  expect_equal(
    projected$par$Ecov_obs_logsigma_re[padded_rows, 2L],
    rep(input$par$Ecov_obs_sigma_par[1L, 2L], 3L)
  )
})
