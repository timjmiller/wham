context("Prepare WHAM input from fitted model")

test_that("input is rebuilt from fitted model options", {
  asap3 <- read_asap3_dat(
    file.path(system.file("extdata", package = "wham"), "ex1_SNEMAYT.dat")
  )
  input <- prepare_wham_input(asap3)
  fit <- list(input = input)

  expect_equal(
    prepare_wham_input_from_fit(fit)$model_name,
    "WHAM for unnamed stock"
  )
})

test_that("missing recruitment decoupling option defaults to coupled", {
  asap3 <- read_asap3_dat(
    file.path(system.file("extdata", package = "wham"), "ex1_SNEMAYT.dat")
  )
  input <- prepare_wham_input(asap3)
  input$options$NAA_re$decouple_recruitment <- NULL
  input$data$decouple_recruitment <- NULL

  rebuilt <- prepare_wham_input_from_fit(list(input = input))

  expect_identical(rebuilt$data$decouple_recruitment, 0L)
})

test_that("invalid fitted model files produce an informative error", {
  expect_error(
    prepare_wham_input_from_fit(list()),
    "does not contain an \\$input component"
  )
})

test_that("recruitment decoupling uses stored flags and active correlations", {
  fit <- list(input = list(data = list(), par = list(trans_NAA_rho = numeric(3))))
  expect_false(refit_decouple_recruitment(fit))
  fit$opt$par <- setNames(numeric(3), rep("trans_NAA_rho", 3))
  expect_true(refit_decouple_recruitment(fit))
  fit$input$map$trans_NAA_rho <- factor(c(1L, 2L, NA_integer_))
  expect_false(refit_decouple_recruitment(fit))
  fit$input$map$trans_NAA_rho <- factor(1:3)
  expect_true(refit_decouple_recruitment(fit))
  fit$input$data$decouple_rec <- 0L
  expect_false(refit_decouple_recruitment(fit))
  fit$input$data$decouple_rec <- 1L
  expect_true(refit_decouple_recruitment(fit))
  fit$input$options$NAA_re$decouple_recruitment <- FALSE
  expect_false(refit_decouple_recruitment(fit))
})

test_that("rebuilding preserves decoupling when preparation options omit it", {
  asap3 <- read_asap3_dat(
    file.path(system.file("extdata", package = "wham"), "ex1_SNEMAYT.dat")
  )
  input <- prepare_wham_input(asap3)
  input$options$NAA_re$decouple_recruitment <- NULL
  input$data$decouple_recruitment <- 1L
  expect_identical(prepare_wham_input_from_fit(list(input = input))$data$decouple_recruitment, 1L)
})

test_that("models without options return their stored input with a warning", {
  input <- list(data = list(), model_name = "old")

  expect_warning(
    expect_identical(prepare_wham_input_from_fit(list(input = input)), input),
    "does not contain preparation options"
  )
})

test_that("refitting supplies missing year indices and preserves stored indices", {
  fit <- list(
    input = list(data = list(n_years_model = 3L)),
    opt = list(par = c(parameter = 1), objective = 0, convergence = 0)
  )
  refit <- refit_wham
  environment(refit) <- list2env(list(
    fit_wham = function(input, ...) {
      result <- fit
      result$input <- input
      result
    }
  ), parent = environment(refit_wham))

  capture.output(result <- refit(fit))
  expect_identical(result$input$data$years_use, 0:2)
  expect_identical(result$input$data$decouple_recruitment, 0L)
  expect_null(fit$input$data$years_use)

  fit$input$data$years_use <- c(0L, 2L)
  fit$input$data$decouple_recruitment <- 1L
  capture.output(result <- refit(fit))
  expect_identical(result$input, fit$input)
})

test_that("refitting saved code 3 uses current AR1 without changing the saved input", {
  local_edition(3)
  input <- suppressMessages(prepare_wham_input(
    NAA_re = list(sigma = "rec"),
    basic_info = list(years = 2000:2005, ages = 1:4),
    ecov = list(year = 2000:2005, mean = matrix(1:6, 6L, 1L),
      use_obs = matrix(1L, 6L, 1L),
      logsigma = -1, process_model = "ar1",
      recruitment_how = matrix("controlling-lag-1-linear", 1L, 1L))
  ))
  input$data$Ecov_model[] <- 3L
  fit <- list(input = input,
    opt = list(par = c(parameter = 1), objective = 0, convergence = 0))
  refit <- refit_wham
  environment(refit) <- list2env(list(fit_wham = function(input, ...) {
    validate_input_safety(input)
    model <- TMB::MakeADFun(input$data, input$par, map = input$map,
      DLL = "wham", silent = TRUE)
    on.exit(TMB::FreeADFun(model))
    expect_equal(is.finite(model$fn(model$par)), TRUE)
    expect_equal(all(is.finite(model$gr(model$par))), TRUE)
    list(input = input, opt = fit$opt)
  }), parent = environment(refit_wham))
  expect_snapshot(capture.output(result <- suppressMessages(refit(fit))))
  expect_equal(as.integer(result$input$data$Ecov_model), 2L)
  expect_identical(fit$input, input)
})

test_that("refitting preserves stored mappings and uses fitted starting values", {
  fit <- list(
    input = list(
      data = list(n_years_model = 3L),
      par = list(parameter = c(0, 7, 0)),
      map = list(parameter = factor(c(1, NA, 1))),
      random = character(),
      options = list()
    ),
    parList = list(parameter = c(2, 7, 2)),
    opt = list(par = c(parameter = 2), objective = 0, convergence = 0)
  )
  original_input <- fit$input
  refit <- refit_wham
  environment(refit) <- list2env(list(
    prepare_wham_input_from_fit = function(...) stop("Input must not be rebuilt"),
    fit_wham = function(input, ...) {
      result <- fit
      result$input <- input
      result
    }
  ), parent = environment(refit_wham))
  capture.output(result <- refit(fit))
  expect_identical(result$input$map, original_input$map)
  expect_identical(result$input$random, original_input$random)
  expect_identical(result$input$par, fit$parList)
  expect_identical(fit$input, original_input)
})

test_that("projection calls refit stored inputs without evaluating the old model", {
  for (call in list(
    quote(project_wham(unavailable_model, unavailable_options,
      do.sdrep = FALSE, n.newton = 0L, check.version = FALSE)),
    quote(wham::project_wham(model = unavailable_model,
      proj.opts = unavailable_options, do.sdrep = FALSE, n.newton = 0L))
  )) {
    fit <- list(
      input = list(data = list(n_years_model = 3L, do_proj = 1L,
        n_years_proj = 2L, proj_F_opt = c(5L, 3L))),
      opt = list(par = c(parameter = 1), objective = 0, convergence = 0),
      call = call
    )
    refit <- refit_wham
    environment(refit) <- list2env(list(
      fit_wham = function(input, ...) {
        arguments <- list(...)
        expect_null(arguments$model)
        expect_null(arguments$proj.opts)
        expect_null(arguments$check.version)
        expect_false(arguments$do.proj)
        expect_false(arguments$do.osa)
        expect_false(arguments$do.post.samp)
        expect_false(arguments$do.sdrep)
        expect_identical(arguments$n.newton, 0L)
        expect_false(arguments$do.retro)
        list(input = input, opt = fit$opt)
      }
    ), parent = environment(refit_wham))

    capture.output(result <- refit(fit))

    expect_identical(result$input$data$do_proj, 1L)
    expect_identical(result$input$data$n_years_proj, 2L)
    expect_identical(result$input$data$proj_F_opt, c(5L, 3L))
  }
})

test_that("stored projections retain controls and separate averaging windows", {
  input <- list(
    data = list(n_years_model = 4L, n_years_proj = 2L, do_proj = 1L,
      n_stocks = 1L, n_regions = 1L, n_fleets = 1L, n_ages = 2L,
      n_selblocks = 1L, selblock_models = 1L, selblock_models_re = 1L,
      selpars_est = matrix(0L, 1L, 6L), proj_F_opt = c(5L, 3L),
      proj_Fcatch = c(100, 0), proj_M_opt = 2L,
      avg_years_ind = c(2L, 3L), avg_years_ind_static = c(0L, 1L),
      FXSPR_init = rep(0.5, 6L), FMSY_init = rep(0.5, 6L)),
    par = list(selpars_re = 0), map = list(selpars_re = factor(NA)),
    years = 2000:2003, years_full = 2000:2005
  )
  fit <- list(input = input,
    opt = list(par = c(parameter = 1), objective = 0, convergence = 0))
  refit <- refit_wham
  environment(refit) <- list2env(list(
    fit_wham = function(input, ...) list(input = input, opt = fit$opt)
  ), parent = environment(refit_wham))

  capture.output(result <- refit(fit))
  data <- result$input$data
  expect_identical(data$do_proj, 1L)
  expect_identical(data$n_years_proj, 2L)
  expect_identical(data$proj_F_opt, input$data$proj_F_opt)
  expect_identical(as.vector(data$proj_Fcatch), input$data$proj_Fcatch)
  expect_identical(dim(data$proj_Fcatch), c(2L, 1L))
  expect_identical(data$proj_M_opt, 2L)
  expect_equal(data$avg_years_ind_sel[, 1L], c(2, 2, 3, 0, 0))
  expect_equal(data$avg_years_ind_static_sel[, 1L], c(2, 0, 1, 0, 0))
  expect_equal(as.vector(data$avg_years_ind_M), c(2, 2, 3, 0, 0))
  expect_equal(as.vector(data$avg_years_ind_static_M), c(2, 0, 1, 0, 0))
  expect_identical(result$input$years_full, input$years_full)
  expect_identical(result$input$par$logR_proj, NULL)
  expect_identical(fit$input, input)

  fit$input$data$avg_years_ind <- c(3L, 4L)
  expect_error(refit(fit), "invalid stored averaging years")
})

test_that("stored historical projections preserve likelihoods and forecast reports", {
  path <- test_path("..", "..", "oldfitobjects", "BSB_2023Run34projfit.RDS")
  skip_if_not(file.exists(path), "Local historical projection fit is unavailable")
  fit <- readRDS(path)
  original_input <- fit$input
  fit$call$do.sdrep <- FALSE
  refit <- refit_wham
  environment(refit) <- list2env(list(
    fit_wham = function(input, ...) {
      expect_identical(input$data$do_proj, original_input$data$do_proj)
      expect_identical(input$data$n_years_proj, original_input$data$n_years_proj)
      expect_identical(input$data$proj_F_opt, original_input$data$proj_F_opt)
      expect_identical(as.vector(input$data$proj_Fcatch), as.vector(original_input$data$proj_Fcatch))
      expect_identical(input$map$log_NAA, original_input$map$log_NAA)
      expect_identical(input$random, original_input$random)
      expect_identical(input$years_full, original_input$years_full)
      for (field in names(fit$parList)) {
        if (!field %in% c("selpars_re", "sel_repars", "Ecov_process_pars")) {
          expect_identical(input$par[[field]], fit$parList[[field]], info = field)
        }
      }
      expected_selectivity <- fit$parList$sel_repars
      expected_selectivity[, 2:3] <- 2 * expected_selectivity[, 2:3]
      expect_equal(input$par$sel_repars, expected_selectivity, tolerance = 1e-12)
      expected_environment <- fit$parList$Ecov_process_pars
      columns <- which(input$data$Ecov_model == 2L & input$data$Ecov_use_re == 1L)
      correlation <- tanh(expected_environment[3L, columns] / 2)
      expected_environment[2L, columns] <- expected_environment[2L, columns] +
        0.5 * log1p(-correlation^2)
      expect_equal(input$par$Ecov_process_pars, expected_environment, tolerance = 1e-12)
      objective <- TMB::MakeADFun(input$data, input$par, random = input$random,
        map = input$map, DLL = "wham", silent = TRUE)
      expect_equal(as.numeric(objective$fn(objective$par)), as.numeric(fit$opt$objective),
        tolerance = 1e-6, scale = 1)
      TMB::FreeADFun(objective)
      objective <- TMB::MakeADFun(input$data, input$par,
        map = input$map, DLL = "wham", silent = TRUE)
      objective$fn(objective$par)
      report <- objective$report()
      for (field in c("MAA", "NAA", "SSB", "pred_catch", "pred_indices")) {
        expect_equal(dim(report[[field]]), dim(fit$rep[[field]]), info = field)
        expect_equal(as.vector(report[[field]]), as.vector(fit$rep[[field]]),
          tolerance = 1e-8, info = field)
      }
      TMB::FreeADFun(objective)
      list(input = input, opt = fit$opt, parList = input$par)
    }
  ), parent = environment(refit_wham))

  capture.output(result <- suppressMessages(refit(fit)))

  expect_identical(fit$input, original_input)
  expect_identical(result$input$data$proj_F_opt, c(5, 3, 3))
})

test_that("packed and array selectivity migrations retain means and shared effects", {
  for (layout in c("packed", "array")) for (model in c(4L, 5L)) {
    shared <- model == 4L
    values <- if (shared) c(0.1, 0.2, 0.1, 0.2) else seq_len(4) / 10
    mapping <- if (shared) c(2, 1, 2, 1) else c(4, 3, 2, 1)
    deviations <- values
    if (layout == "array") {
      deviations <- array(0, c(1L, 2L, 4L))
      deviations[1, , 1:2] <- values
      array_mapping <- array(NA_integer_, dim(deviations))
      array_mapping[1, , 1:2] <- mapping
      mapping <- array_mapping
    }
    input <- list(
      data = list(
        n_years_model = 2L, n_selblocks = 1L, n_ages = 2L,
        selblock_models = 1L, selblock_models_re = model,
        selblock_years = matrix(1L, 2L, 1L), n_selpars_est = 2L,
        selpars_est = matrix(c(1L, 1L, rep(0L, 6)), 1L),
        selpars_lower = matrix(0, 1L, 8L), selpars_upper = matrix(1, 1L, 8L),
        do_proj = 0L, proj_Fcatch = matrix(0)
      ),
      par = list(logit_selpars = matrix(seq_len(8), 1L), selpars_re = deviations),
      map = list(logit_selpars = factor(c(2, 1, rep(NA, 6))), selpars_re = factor(mapping))
    )
    fit <- list(input = input, opt = list(par = c(parameter = 1), objective = 0, convergence = 0))
    refit <- refit_wham
    environment(refit) <- list2env(list(
      set_proj = identity,
      set_brp_static = function(input, ...) input,
      fit_wham = function(input, ...) list(input = input, opt = fit$opt)
    ), parent = environment(refit_wham))
    capture.output(result <- refit(fit))
    migrated <- result$input
    expect_identical(migrated$par$logit_selpars, input$par$logit_selpars)
    expect_identical(migrated$map$logit_selpars, input$map$logit_selpars)
    expect_identical(migrated$data$selpars_lower, input$data$selpars_lower)
    expect_identical(migrated$data$selpars_upper, input$data$selpars_upper)
    expect_equal(migrated$data$n_selpars_re, if (shared) 1L else 2L)
    expect_equal(migrated$data$selpars_re_index[1, 1:2], if (shared) c(1L, 1L) else 1:2)
    expect_equal(as.vector(migrated$par$selpars_re[1, , seq_len(if (shared) 1L else 2L)]),
                 if (shared) values[1:2] else values)
    expect_identical(levels(migrated$map$selpars_re), levels(input$map$selpars_re))
  }
})

test_that("stored BSB selectivity arrays preserve parameters and likelihood", {
  fit_file <- test_path("..", "..", "oldfitobjects", "BSB_2024MTfit.RDS")
  skip_if_not(file.exists(fit_file), "Local BSB fit is unavailable")
  fit <- readRDS(fit_file)
  original_input <- fit$input
  refit <- refit_wham
  environment(refit) <- list2env(list(
    fit_wham = function(input, ...) {
      expect_identical(input$par$logit_selpars, fit$parList$logit_selpars)
      expect_identical(input$map$logit_selpars, fit$input$map$logit_selpars)
      expect_identical(input$map$log_index_sig_scale, fit$input$map$log_index_sig_scale)
      expect_identical(input$data$obsvec, fit$input$data$obsvec)
      model <- TMB::MakeADFun(
        input$data, input$par, DLL = "wham", random = input$random,
        map = input$map, silent = TRUE
      )
      expect_equal(model$par, fit$opt$par)
      expect_equal(as.numeric(model$fn(model$par)), as.numeric(fit$opt$objective), tolerance = 1e-8)
      expect_lt(max(abs(model$gr(model$par))), 1e-6)
      report <- model$report()
      for (field in c("SSB", "NAA", "pred_catch", "pred_indices", "nll_sel")) {
        expect_equal(as.vector(report[[field]]), as.vector(fit$rep[[field]]), tolerance = 1e-8)
      }
      list(input = input, opt = fit$opt, parList = input$par)
    }
  ), parent = environment(refit_wham))
  capture.output(result <- suppressMessages(refit(fit)))
  expect_length(result$opt$par, 211L)
  expect_identical(fit$input, original_input)
})

test_that("legacy migration defaults missing recruitment information before rebuilding", {
  input <- list(
    asap3 = list(n_ages = 5L),
    data = list(
      Ecov_model = 0, selblock_models_re = 1, M_re_model = 1,
      M_model = 2, use_steepness = 0, do_proj = 0,
      decouple_rec_unused = 1L
    )
  )
  migrate <- migrate_legacy_refit_input
  environment(migrate) <- list2env(list(
    prepare_wham_input = function(..., NAA_re) {
      expect_identical(NAA_re$decouple_recruitment, FALSE)
      stop("Reached input rebuilding")
    }
  ), parent = environment(migrate_legacy_refit_input))

  expect_error(migrate(list(input = input)), "Reached input rebuilding")
  expect_null(input$data[["decouple_rec", exact = TRUE]])
})

test_that("unsupported legacy configurations are rejected before rebuilding", {
  input <- list(
    asap3 = list(n_ages = 5L),
    data = list(
      Ecov_model = 0, selblock_models_re = 1, M_re_model = 1,
      M_model = 2, decouple_rec = 0, use_steepness = 0, do_proj = 0
    )
  )
  unsupported <- list(
    Ecov_model = 3, selblock_models_re = 6, M_re_model = 6,
    M_model = 4, use_steepness = 1, do_proj = 1
  )
  for (field in names(unsupported)) {
    legacy_input <- input
    legacy_input$data[[field]] <- unsupported[[field]]
    expect_error(
      migrate_legacy_refit_input(list(input = legacy_input)),
      "Cannot safely migrate this legacy WHAM input"
    )
  }
})

test_that("packed legacy age-composition positions preserve series order", {
  data <- list(
    n_age_comp_pars_fleets = c(0L, 1L, 2L, 3L),
    n_age_comp_pars_indices = c(1L, 0L)
  )
  expect_equal(
    refit_legacy_age_comp_indices(data, "catch_paa_pars"),
    c(2L, 3L, 7L, 4L, 8L, 12L)
  )
  expect_equal(refit_legacy_age_comp_indices(data, "index_paa_pars"), 1L)
  data$n_age_comp_pars_fleets <- c(0L, 4L)
  expect_error(refit_legacy_age_comp_indices(data, "catch_paa_pars"), "Cannot safely migrate")
})

test_that("v1.0.0 fits migrate without ASAP3 data", {
  release_dir <- test_path("..", "..", "oldfitobjects", "release_refits", "v1.0.0")
  paths <- list.files(release_dir, "^m[1-4]_original\\.rds$", full.names = TRUE)
  skip_if_not(length(paths) > 0L, "Generate fixtures with tests/manual/refit-releases.R")
  for (path in paths) {
    fit <- readRDS(path)
    original_input <- fit$input
    expect_null(fit$input$asap3)
    input <- suppressMessages(migrate_legacy_refit_input(fit))
    expect_identical(fit$input, original_input)
    expect_identical(input$years, fit$input$years)
    expect_identical(input$ages.lab, fit$input$ages.lab)
    expect_equal(input$data$which_F_age, rep(fit$input$data$which_F_age, length(fit$years)))
    expect_equal(as.vector(input$par$Mpars), fit$parList$M0 + fit$parList$M_a)
    expect_equal(as.vector(input$par$trans_NAA_rho), c(2 * fit$parList$trans_NAA_rho, 0))
    expect_true(is.na(tail(input$map$trans_NAA_rho, 1)))
    for (field in c("agg_catch", "catch_paa", "agg_indices", "index_paa",
                    "use_agg_catch", "use_catch_paa", "use_indices", "use_index_paa")) {
      expect_identical(input$data[[field]], fit$input$data[[field]])
    }
    expect_equal(dim(input$data$keep_Cpaa), c(input$data$n_fleets, input$data$n_years_model, 2L))
    expect_equal(length(input$data$agesvec), length(input$data$obsvec))
    if (any(fit$input$data$age_comp_model_fleets == 7L)) {
      expect_true(all(input$data$age_comp_model_fleets == 5L))
      expect_true(all(input$data$age_comp_model_indices == 5L))
      expect_true(all(input$data$catch_Neff == 1))
      expect_true(all(input$data$index_Neff == 1))
      expect_identical(input$legacy_age_comp_Neff$catch, fit$input$data$catch_Neff)
      expect_identical(input$legacy_age_comp_Neff$index, fit$input$data$index_Neff)
    }
    model <- TMB::MakeADFun(
      input$data, input$par, DLL = "wham", random = input$random,
      map = input$map, silent = TRUE
    )
    expect_equal(as.numeric(model$fn(model$par)), as.numeric(fit$opt$objective), tolerance = 1e-6, scale = 1)
    report <- model$report()
    for (field in c("MAA", "NAA", "SSB", "pred_catch", "pred_indices")) {
      expect_equal(as.vector(report[[field]]), as.vector(fit$rep[[field]]), tolerance = 1e-8)
    }
    TMB::FreeADFun(model)
    fit$call <- quote(fit_wham(input = input, do.retro = FALSE, do.osa = FALSE,
      do.sdrep = FALSE, do.brps = FALSE, do.post.samp = FALSE, MakeADFun.silent = TRUE))
    capture.output(result <- suppressWarnings(suppressMessages(refit_wham(fit))))
    expect_equal(result$opt$objective, fit$opt$objective, tolerance = 1e-6, scale = 1)
    expect_equal(result$opt$convergence, 0)
    expect_identical(fit$input, original_input)
    TMB::FreeADFun(result)
  }
  fit <- readRDS(paths[1])
  fit$input$map$M0 <- factor(1L)
  expect_error(suppressMessages(migrate_legacy_refit_input(fit)), "estimated legacy M0")
  fit$input$years <- NULL
  expect_error(suppressMessages(migrate_legacy_refit_input(fit)), "stored model years")
  fit <- readRDS(paths[1])
  fit$input$data$age_comp_model_fleets[] <- 5L
  expect_error(suppressMessages(migrate_legacy_refit_input(fit)), "unsupported early age-composition likelihood")
})

test_that("tagged v1.0.1 through v1.0.6 fits preserve stored likelihoods", {
  root <- test_path("..", "..", "oldfitobjects", "release_refits")
  paths <- unlist(lapply(paste0("v1.0.", 1:6), function(tag) {
    list.files(file.path(root, tag), "_original\\.rds$", full.names = TRUE)
  }), use.names = FALSE)
  skip_if_not(length(paths) > 0L, "Generate fixtures with tests/manual/refit-releases.R")
  for (path in paths) {
    fit <- readRDS(path)
    reference <- readRDS(sub("_original\\.rds$", "_reference.rds", path))
    original_input <- fit$input
    input <- suppressMessages(migrate_legacy_refit_input(fit))
    expect_identical(fit$input, original_input, info = path)
    expect_equal(input$years, original_input$years, info = path)
    expect_length(input$data$which_F_age, input$data$n_years_model)
    expect_equal(input$data$bias_correct_pe, original_input$data$bias_correct_pe)
    expect_equal(input$data$bias_correct_oe, original_input$data$bias_correct_oe)
    if (is.matrix(original_input$par$catch_paa_pars)) {
      expect_identical(input$data$age_comp_model_fleets, original_input$data$age_comp_model_fleets)
      expect_identical(input$data$age_comp_model_indices, original_input$data$age_comp_model_indices)
      expect_identical(input$data$catch_Neff, original_input$data$catch_Neff)
      expect_identical(input$data$index_Neff, original_input$data$index_Neff)
      expect_null(input$legacy_age_comp_Neff)
    }
    model <- TMB::MakeADFun(
      input$data, input$par, DLL = "wham", random = input$random,
      map = input$map, silent = TRUE
    )
    expect_equal(as.numeric(model$fn(model$par)), as.numeric(fit$opt$objective),
      tolerance = 1e-6, scale = 1, info = path)
    TMB::FreeADFun(model)
    model <- TMB::MakeADFun(
      input$data, input$par, DLL = "wham", map = input$map, silent = TRUE
    )
    model$fn(model$par)
    report <- model$report()
    for (field in c("MAA", "NAA", "SSB", "pred_catch", "pred_indices")) {
      expect_equal(as.vector(report[[field]]), as.vector(reference[[field]]),
        tolerance = 1e-7, info = paste(path, field))
    }
    TMB::FreeADFun(model)
  }
})

test_that("early logistic-normal weighting requires a matching stored report", {
  root <- test_path("..", "..", "oldfitobjects", "release_refits")
  paths <- file.path(root, c("v1.0.2", "v1.0.3"), "m2_original.rds")
  skip_if_not(all(file.exists(paths)), "Generate v1.0.2 and v1.0.3 fixtures")
  for (position in seq_along(paths)) {
    fit <- readRDS(paths[position])
    input <- suppressMessages(migrate_legacy_refit_input(fit))
    expected <- fit$input$data$catch_Neff
    if (position == 1L) expected[] <- 1
    expect_identical(input$data$catch_Neff, expected)
    expect_identical(input$legacy_age_comp_Neff$catch, fit$input$data$catch_Neff)
  }
  fit <- readRDS(paths[2L])
  fit$rep$nll_catch_acomp <- fit$rep$nll_catch_acomp + 1
  expect_error(suppressMessages(migrate_legacy_refit_input(fit)),
    "Cannot safely infer historical age-composition weighting")
  fit$rep$nll_catch_acomp <- NULL
  expect_error(suppressMessages(migrate_legacy_refit_input(fit)),
    "without stored nll_catch_acomp")
})

test_that("legacy process and environmental features preserve stored likelihoods", {
  root <- test_path("..", "..", "oldfitobjects", "release_refits")
  paths <- unlist(lapply(c("v1.0.0", "v1.0.5", "v1.0.6"), function(tag) {
    list.files(file.path(root, tag, "features"), "_original\\.rds$", full.names = TRUE)
  }), use.names = FALSE)
  skip_if_not(length(paths) > 0L, "Generate feature fixtures with tests/manual/refit-releases.R")
  for (path in paths) {
    fit <- readRDS(path)
    reference <- readRDS(sub("_original\\.rds$", "_reference.rds", path))
    original_input <- fit$input
    if (isTRUE(original_input$data$use_steepness == 1L)) next
    warnings <- character()
    input <- withCallingHandlers(suppressMessages(migrate_legacy_refit_input(fit)),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    if (any(grepl("Historical environmental AR1", warnings, fixed = TRUE))) {
      expect_equal(as.integer(input$data$Ecov_model), 2L)
      next
    }
    expect_identical(fit$input, original_input, info = path)
    expect_equal(input$data$M_re_model[1], c(1L, 4L, 2L, 3L, 4L)[original_input$data$M_re_model])
    expect_equal(input$data$selblock_models_re, original_input$data$selblock_models_re - 1L)
    expect_equal(input$data$M_model, if (original_input$data$M_model == 3L) 2L else 1L)
    model <- TMB::MakeADFun(input$data, input$par, random = input$random,
      map = input$map, DLL = "wham", silent = TRUE)
    expect_length(model$par, length(fit$opt$par))
    expect_equal(as.numeric(model$fn(model$par)), as.numeric(fit$opt$objective),
      tolerance = 1e-6, scale = 1, info = path)
    TMB::FreeADFun(model)
    model <- TMB::MakeADFun(input$data, input$par, map = input$map, DLL = "wham", silent = TRUE)
    model$fn(model$par)
    report <- model$report()
    for (field in names(reference)) {
      expect_equal(as.vector(report[[field]]), as.vector(reference[[field]]),
        tolerance = 1e-7, info = paste(path, field))
    }
    TMB::FreeADFun(model)
  }
})

test_that("historical environmental AR1 migrates to the current process", {
  path <- test_path("..", "..", "oldfitobjects", "release_refits", "v1.0.0",
    "features", "environment_ar1_recruit_original.rds")
  skip_if_not(file.exists(path), "Generate the earliest environmental fixture")
  fit <- readRDS(path)
  local_edition(3)
  expect_snapshot(input <- suppressMessages(migrate_legacy_refit_input(fit)))
  expect_equal(as.integer(input$data$Ecov_model), 2L)
  for (field in grep("^do_simulate_", names(input$data), value = TRUE)) input$data[[field]][] <- 0L
  input$data$do_simulate_Ecov_re <- 1L
  input$data$do_simulate_period <- c(1L, 0L)
  model <- TMB::MakeADFun(input$data, input$par, map = input$map, DLL = "wham", silent = TRUE)
  on.exit(TMB::FreeADFun(model))
  simulated <- model$simulate(model$par)
  expect_equal(all(is.finite(simulated$Ecov_re)), TRUE)
  expect_equal(as.vector(simulated$Ecov_x),
    as.vector(simulated$Ecov_re + input$par$Ecov_process_pars[1L, 1L]))
})

test_that("legacy steepness is not silently treated as alpha and beta", {
  local_edition(3)
  root <- test_path("..", "..", "oldfitobjects", "release_refits", "v1.0.6", "features")
  paths <- file.path(root, paste0(c("steepness_bh", "steepness_ricker"), "_original.rds"))
  skip_if_not(all(file.exists(paths)), "Generate the steepness fixtures")
  for (path in paths) {
    fit <- readRDS(path)
    expect_snapshot(error = TRUE, migrate_legacy_refit_input(fit))
  }
})

test_that("legacy Butterfish migration preserves the likelihood and estimates", {
  fit_file <- test_path("..", "..", "oldfitobjects", "Butterfishfit.rds")
  skip_if_not(file.exists(fit_file), "Local legacy fit is unavailable")
  fit <- readRDS(fit_file)
  original_input <- fit$input
  input <- suppressMessages(migrate_legacy_refit_input(fit))
  expect_identical(fit$input, original_input)
  expect_equal(input$data$years_use, 0:34)
  expect_equal(input$data$bias_correct_pe, fit$input$data$bias_correct_pe)
  expect_equal(input$data$bias_correct_oe, fit$input$data$bias_correct_oe)
  expect_identical(input$data$obsvec, fit$input$data$obsvec)
  expect_identical(input$map$logit_q, fit$input$map$logit_q)
  expect_equal(as.vector(input$par$Mpars), fit$parList$M_a)
  expect_equal(as.vector(input$par$trans_NAA_rho), 2 * fit$parList$trans_NAA_rho)

  model <- TMB::MakeADFun(
    input$data, input$par, DLL = "wham", random = input$random,
    map = input$map, silent = TRUE
  )
  expect_equal(as.numeric(model$fn(model$par)), as.numeric(fit$opt$objective), tolerance = 1e-8)
  report <- model$report()
  for (field in c("nll_NAA", "nll_agg_catch", "nll_agg_indices", "nll_catch_acomp",
                  "nll_index_acomp", "SSB", "NAA", "pred_catch", "pred_indices")) {
    expect_equal(as.vector(report[[field]]), as.vector(fit$rep[[field]]), tolerance = 1e-8)
  }

  fit$call <- quote(fit_wham(input = input, do.retro = FALSE, do.osa = FALSE,
    do.sdrep = FALSE, do.brps = FALSE, do.post.samp = FALSE, MakeADFun.silent = TRUE))
  output <- capture.output(result <- suppressMessages(refit_wham(fit)))
  expect_equal(result$opt$convergence, 0)
  expect_equal(as.numeric(result$opt$objective), as.numeric(fit$opt$objective), tolerance = 1e-8)
  expect_true(any(grepl("max_abs_par_diff", output)))
})

test_that("legacy decoupled Cod migration preserves the likelihood", {
  fit_file <- test_path("..", "..", "oldfitobjects", "GBK_Cod_2024MTfit.rds")
  skip_if_not(file.exists(fit_file), "Local decoupled fit is unavailable")
  fit <- readRDS(fit_file)
  original_input <- fit$input
  input <- suppressMessages(migrate_legacy_refit_input(fit))
  expect_identical(fit$input, original_input)
  expect_identical(input$data$decouple_recruitment, 1L)
  expect_equal(as.vector(input$par$trans_NAA_rho), 2 * fit$parList$trans_NAA_rho[c(1L, 3L, 2L)])

  model <- TMB::MakeADFun(
    input$data, input$par, DLL = "wham", random = input$random,
    map = input$map, silent = TRUE
  )
  expect_equal(as.numeric(model$fn(model$par)), as.numeric(fit$opt$objective), tolerance = 1e-8)

  refit <- refit_wham
  environment(refit) <- list2env(list(
    fit_wham = function(input, ...) list(input = input, opt = fit$opt, parList = input$par)
  ), parent = environment(refit_wham))
  expect_message(
    output <- capture.output(result <- suppressMessages(refit(fit))),
    NA
  )
  expect_identical(result$input$data$decouple_recruitment, 1L)
  expect_true(any(grepl("max_abs_par_diff", output)))

  fit$input$data$decouple_rec <- NULL
  inferred <- suppressMessages(migrate_legacy_refit_input(fit))
  expect_identical(inferred$data$decouple_recruitment, 1L)
  expect_equal(inferred$par, input$par)
  expect_identical(inferred$map, input$map)
})
