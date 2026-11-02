args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) args[1] else "generate"
tag <- if (length(args) >= 2L) args[2] else "v1.0.0"
output_root <- Sys.getenv(
  "WHAM_RELEASE_RESULTS",
  file.path("oldfitobjects", "release_refits")
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, winslash = "/")
release_dir <- file.path(output_root, tag)
if (stage %in% c("generate-features", "refit-features")) {
  release_dir <- file.path(release_dir, "features")
}
dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)

git_output <- function(arguments) {
  output <- system2("git", arguments, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(output, "status"))) stop(paste(output, collapse = "\n"))
  output
}

save_reference_report <- function(fit, path) {
  pars <- fit$input$par
  pars[names(fit$parList)] <- fit$parList
  objective <- TMB::MakeADFun(fit$input$data, pars, map = fit$input$map,
    DLL = "wham", silent = TRUE)
  on.exit(TMB::FreeADFun(objective))
  objective$fn(objective$par)
  report <- objective$report()
  saveRDS(report[c("MAA", "NAA", "SSB", "pred_catch", "pred_indices")], path)
}

if (stage == "summarize-features") {
  tags <- if (length(args) > 1L) args[-1L] else c("v1.0.0", "v1.0.5", "v1.0.6")
  expected <- c("selectivity_iid", "selectivity_ar1", "selectivity_shared_year",
    "mortality_ar1_y", "mortality_constant_ar1", "mortality_2dar1_fixed",
    "environment_ar1_recruit", "environment_rw_recruit", "environment_ar1_mortality",
    "mortality_weight", "mortality_weight_ar1", "steepness_bh", "steepness_ricker")
  summary <- do.call(rbind, lapply(tags, function(release) {
    target <- file.path(output_root, release, "features")
    fits <- list.files(target, "_original\\.rds$")
    refitted <- read.csv(file.path(target, "refit-features_summary.csv"),
      colClasses = c(warnings = "character", error = "character"))
    complete <- nrow(refitted) == length(expected) &&
      setequal(refitted$model, expected) &&
      setequal(sub("_original\\.rds$", "", fits), expected)
    data.frame(tag = release, cases = nrow(refitted),
      stored_nll_matches = sum(refitted$stored_nll_matches),
      refit_nll_matches = sum(refitted$nll_matches),
      max_stored_nll_difference = max(abs(refitted$stored_nll_difference)),
      max_refit_nll_difference = max(abs(refitted$nll_difference)),
      original_convergence_warnings = sum(refitted$convergence_original != 0L),
      cases_with_warnings = sum(nzchar(refitted$warnings)),
      errors = sum(nzchar(refitted$error)),
      passed = complete && all(refitted$stored_nll_matches) &&
        all(refitted$nll_matches) && all(refitted$convergence == 0L) &&
        all(refitted$error == ""))
  }))
  write.csv(summary, file.path(output_root, "feature_summary.csv"), row.names = FALSE)
  print(summary)
  stopifnot(all(summary$passed))
} else if (stage == "reference") {
  library(TMB)
  dyn.load(TMB::dynlib(file.path(release_dir, "source", "src", "wham")))
  paths <- list.files(release_dir, "_original\\.rds$", full.names = TRUE)
  if (!length(paths)) stop("Generate historical fits first.")
  for (path in paths) {
    save_reference_report(readRDS(path), sub("_original\\.rds$", "_reference.rds", path))
  }
  message("Saved ", length(paths), " historical parameter-level reports for ", tag)
} else if (stage == "summarize") {
  tags <- if (length(args) > 1L) args[-1L] else paste0("v1.0.", 1:6)
  expected <- c(paste0("m", 1:4), "scaa", "rec_ar1_y", "naa_ar1_a", "naa_ar1_y", "naa_2dar1")
  summary <- do.call(rbind, lapply(tags, function(release) {
    target <- file.path(output_root, release)
    generated <- read.csv(file.path(target, "generation_summary.csv"))
    refitted <- read.csv(file.path(target, "refit_summary.csv"),
      colClasses = c(warnings = "character", error = "character"))
    archived <- !is.na(generated$archived_nll)
    complete <- nrow(generated) == length(expected) && nrow(refitted) == length(expected) &&
      setequal(generated$model, expected) && setequal(refitted$model, expected)
    data.frame(
      tag = release,
      cases = nrow(refitted),
      archived_matches = sum(archived & is.finite(generated$difference) & abs(generated$difference) <= 1e-6),
      stored_nll_matches = sum(refitted$stored_nll_matches),
      refit_nll_matches = sum(refitted$nll_matches),
      max_stored_nll_difference = max(abs(refitted$stored_nll_difference)),
      max_refit_nll_difference = max(abs(refitted$nll_difference)),
      cases_with_warnings = sum(nzchar(refitted$warnings)),
      errors = sum(nzchar(refitted$error)),
      passed = complete && sum(archived) == 4L &&
        all(abs(generated$difference[archived]) <= 1e-6) &&
        all(refitted$stored_nll_matches) && all(refitted$nll_matches) &&
        all(refitted$convergence == 0L) && all(refitted$error == "")
    )
  }))
  write.csv(summary, file.path(output_root, "release_summary.csv"), row.names = FALSE)
  print(summary)
  stopifnot(all(summary$passed))
} else if (stage == "run") {
  tags <- if (length(args) > 1L) args[-1L] else paste0("v1.0.", 1:6)
  script <- normalizePath("tests/manual/refit-releases.R", winslash = "/")
  executable <- file.path(R.home("bin"), "Rscript.exe")
  if (.Platform$OS.type != "windows") executable <- file.path(R.home("bin"), "Rscript")
  status <- data.frame()
  for (release in tags) {
    target <- file.path(output_root, release)
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
    for (action in c("generate", "refit")) {
      message(release, " ", action)
      log_file <- file.path(target, paste0(action, ".log"))
      exit_code <- system2(executable, c(shQuote(script), action, release),
        stdout = log_file, stderr = log_file)
      status <- rbind(status, data.frame(tag = release, stage = action,
        exit_code = exit_code, log = log_file))
      write.csv(status, file.path(output_root, "release_status.csv"), row.names = FALSE)
      if (exit_code != 0L) stop("Failed ", release, " ", action, "; see ", log_file)
    }
  }
  print(status)
} else if (stage %in% c("generate", "generate-features")) {
  library(TMB)
  library(dplyr)
  library(tidyr)
  if (tag == "current") {
    current_library <- Sys.getenv("WHAM_CURRENT_LIBRARY")
    if (nzchar(current_library)) library(wham, lib.loc = current_library) else library(wham)
    source_dir <- normalizePath(".", winslash = "/")
    commit <- git_output(c("rev-parse", "HEAD"))
    description <- read.dcf(file.path(source_dir, "DESCRIPTION"))
    stopifnot(as.character(packageVersion("wham")) == description[1, "Version"])
    historical <- asNamespace("wham")
  } else {
    commit <- git_output(c("rev-parse", paste0(tag, "^{commit}")))
    source_dir <- file.path(output_root, tag, "source")
    archive <- file.path(release_dir, "source.zip")
    git_output(c("archive", "--format=zip", paste0("--output=", shQuote(archive)), tag))
    unzip(archive, exdir = source_dir)
    description <- read.dcf(file.path(source_dir, "DESCRIPTION"))
    stopifnot(paste0("v", description[1, "Version"]) == tag)
    original_dir <- getwd()
    setwd(file.path(source_dir, "src"))
    if (!file.exists("wham.cpp")) {
      template <- if (file.exists("wham_v0.cpp")) "wham_v0.cpp" else "multi_wham.cpp"
      stopifnot(file.copy(template, "wham.cpp"))
    }
    makevars <- if (.Platform$OS.type == "windows" && file.exists("Makevars.win")) {
      "Makevars.win"
    } else "Makevars"
    framework <- if (any(grepl("^[^#]*-DTMBAD_FRAMEWORK", readLines(makevars, warn = FALSE)))) {
      "TMBad"
    } else "CppAD"
    if (!file.exists(dynlib("wham"))) {
      compile("wham.cpp", flags = "-O1 -g0", framework = framework)
    }
    dyn.load(dynlib("wham"))
    setwd(original_dir)
    historical <- new.env(parent = globalenv())
    for (path in list.files(file.path(source_dir, "R"), "\\.[Rr]$", full.names = TRUE)) {
      sys.source(path, envir = historical)
    }
  }
  provenance <- list(
    tag = tag,
    commit = commit,
    version = description[1, "Version"],
    R = R.version.string,
    TMB = as.character(packageVersion("TMB")),
    AD_framework = if (tag == "current") "installed" else framework,
    workspace_dirty = length(git_output(c("status", "--porcelain"))) > 0L,
    generated = Sys.time()
  )
  saveRDS(provenance, file.path(release_dir, "provenance.rds"))
  asap3 <- historical$read_asap3_dat(
    file.path(source_dir, "inst", "extdata", "ex1_SNEMAYT.dat")
  )
  example <- new.env(parent = historical)
  example$asap3 <- asap3
  example_paths <- file.path(source_dir, "tests", "testthat",
    c("test_ex1_SNEMAYT.R", "test_ex01_SNEMAYT.R"))
  expressions <- parse(example_paths[file.exists(example_paths)][1L])
  test <- Filter(function(expression) {
    is.call(expression) && identical(expression[[1]], as.name("test_that"))
  }, as.list(expressions))[[1L]]
  baseline_call <- NULL
  setup_expressions <- function(expression) {
    if (!is.call(expression)) return(list())
    operator <- as.character(expression[[1L]])[1L]
    if (operator %in% c("<-", "=")) return(list(expression))
    if (!operator %in% c("{", "suppressWarnings", "suppressMessages")) return(list())
    unlist(lapply(as.list(expression)[-1L], setup_expressions), recursive = FALSE)
  }
  for (expression in setup_expressions(test[[3L]])) {
    if (!is.call(expression) || !as.character(expression[[1L]]) %in% c("<-", "=")) next
    target <- expression[[2L]]
    while (is.call(target)) target <- target[[2L]]
    if (identical(target, as.name("basic_info"))) eval(expression, example)
    if (is.symbol(target) && as.character(target) %in% paste0("input", 1:4)) {
      eval(expression, example)
      if (identical(expression[[2L]], as.name("input1"))) baseline_call <- expression[[3L]]
    }
  }
  stopifnot(is.call(baseline_call))
  cases <- list(
    scaa = NULL,
    rec_ar1_y = list(sigma = "rec", cor = "ar1_y"),
    naa_ar1_a = list(sigma = "rec+1", cor = "ar1_a"),
    naa_ar1_y = list(sigma = "rec+1", cor = "ar1_y"),
    naa_2dar1 = list(sigma = "rec+1", cor = "2dar1")
  )
  feature_cases <- list()
  if (stage == "generate-features") {
    environment_data <- read.csv(file.path(source_dir, "inst", "extdata", "CPI.csv"))
    selection <- list(model = rep("logistic", 3), re = c("iid", "none", "none"),
      initial_pars = rep(list(c(2, 0.3)), 3), fix_pars = rep(list(NULL), 3))
    feature_cases <- list(
      selectivity_iid = list(selectivity = selection),
      mortality_ar1_y = list(M = list(model = "age-specific", re = "ar1_y")),
      environment_ar1_recruit = list(ecov = list(
        label = "CPI", mean = as.matrix(environment_data$CPI),
        logsigma = as.matrix(log(environment_data$CPI_sigma)),
        year = environment_data$Year,
        use_obs = matrix(1, nrow(environment_data), 1),
        lag = 1, process_model = "ar1", where = "recruit", how = 1
      ))
    )
    if ("age_comp" %in% names(formals(historical$prepare_wham_input))) {
      feature_cases$selectivity_iid$age_comp <- "logistic-normal-miss0"
    }
    feature_cases$selectivity_ar1 <- feature_cases$selectivity_iid
    feature_cases$selectivity_ar1$selectivity$re[1] <- "ar1"
    feature_cases$selectivity_shared_year <- feature_cases$selectivity_iid
    feature_cases$selectivity_shared_year$selectivity <- list(
      model = rep("age-specific", 3), re = c("ar1_y", "none", "none"),
      initial_pars = list(c(0.1, 0.5, 0.5, 1, 1, 1), c(0.5, 0.5, 0.5, 1, 0.5, 0.5), c(0.5, 0.5, 1, 1, 1, 1)),
      fix_pars = list(4:6, 4, 3:6)
    )
    feature_cases$mortality_constant_ar1 <- list(M = list(
      model = "constant", initial_means = 0.3, re = "ar1_y"
    ))
    feature_cases$mortality_2dar1_fixed <- list(M = list(model = "age-specific", re = "2dar1"))
    feature_cases$environment_rw_recruit <- feature_cases$environment_ar1_recruit
    feature_cases$environment_rw_recruit$ecov$process_model <- "rw"
    feature_cases$environment_ar1_mortality <- feature_cases$environment_ar1_recruit
    feature_cases$environment_ar1_mortality$ecov$where <- "M"
    feature_cases$environment_ar1_mortality$ecov$lag <- 0
    feature_cases$mortality_weight <- list(M = list(
      model = "weight-at-age", re = "none", est_ages = 1L, initial_means = 0.3
    ))
    feature_cases$mortality_weight_ar1 <- list(M = list(
      model = "weight-at-age", re = "ar1_y", est_ages = 1L, initial_means = 0.3
    ))
    feature_cases$steepness_bh <- list(recruit_model = 3L,
      NAA_re = list(sigma = "rec+1", cor = "iid", use_steepness = TRUE,
        recruit_pars = c(0.7, 100000)))
    feature_cases$steepness_ricker <- list(recruit_model = 4L,
      NAA_re = list(sigma = "rec+1", cor = "iid", use_steepness = TRUE,
        recruit_pars = c(0.7, 100000)))
  }
  reference <- readRDS(file.path(source_dir, "inst", "extdata", "ex1_test_results.rds"))
  summary <- data.frame()
  model_ids <- if (stage == "generate-features") names(feature_cases) else c(paste0("m", seq_len(4L)), names(cases))
  if (stage == "generate-features" && length(args) > 2L) model_ids <- intersect(model_ids, args[-c(1L, 2L)])
  summary_path <- file.path(release_dir, "generation_summary.csv")
  if (stage == "generate-features" && file.exists(summary_path)) {
    summary <- read.csv(summary_path)
    summary <- summary[!summary$model %in% model_ids, , drop = FALSE]
  }
  generation_errors <- data.frame(model = character(), error = character())
  for (model_id in model_ids) {
    tryCatch({
    if (stage == "generate-features") {
      prepare_call <- baseline_call
      prepare_call$NAA_re <- list(sigma = "rec+1", cor = "iid")
      for (option in names(feature_cases[[model_id]])) {
        prepare_call[[option]] <- feature_cases[[model_id]][[option]]
      }
      input <- eval(prepare_call, example)
      if (startsWith(model_id, "steepness_")) {
        years <- input$data$n_years_model
        if (is.matrix(input$data$mature)) {
          input$data$mature[,] <- rep(input$data$mature[1L, ], each = years)
        } else {
          for (year in seq_len(years)[-1L]) {
            input$data$mature[, year, ] <- input$data$mature[, 1L, ]
          }
        }
        pointer <- input$data$waa_pointer_ssb[1L]
        weights <- input$data$waa[pointer, 1L, ]
        for (year in seq_len(input$data$n_years_model)) {
          input$data$waa[pointer, year, ] <- weights
        }
        input$data$fracyr_SSB[] <- input$data$fracyr_SSB[1L]
        input$data$M_re_model[] <- 1L
        input$data$use_steepness <- 1L
        settings <- feature_cases[[model_id]]$NAA_re$recruit_pars
        input$par$mean_rec_pars[] <- c(
          if (model_id == "steepness_bh") qlogis((settings[1L] - 0.2) / 0.8) else
            log(settings[1L] - 0.2),
          log(settings[2L])
        )
      }
      if (startsWith(model_id, "selectivity_") && !"age_comp" %in% names(formals(historical$prepare_wham_input))) {
        input$data$age_comp_model_fleets <- rep(7L, input$data$n_fleets)
        input$data$age_comp_model_indices <- rep(7L, input$data$n_indices)
        input$data$n_age_comp_pars_fleets <- rep(1L, input$data$n_fleets)
        input$data$n_age_comp_pars_indices <- rep(1L, input$data$n_indices)
        input$par$catch_paa_pars <- rep(0, input$data$n_fleets)
        input$par$index_paa_pars <- rep(0, input$data$n_indices)
        input$map[c("catch_paa_pars", "index_paa_pars")] <- NULL
      }
      if (model_id == "mortality_2dar1_fixed") {
        input$par$M_repars <- c(log(0.1), atanh(0.2), atanh(0.3))
        input$map$M_repars <- factor(rep(NA_integer_, 3L))
      }
      if (startsWith(model_id, "mortality_weight")) {
        weights <- input$data$waa
        count <- dim(weights)[1L]
        input$data$waa <- array(0, c(count + 1L, dim(weights)[-1L]))
        input$data$waa[seq_len(count), , ] <- weights
        input$data$waa[count + 1L, , ] <- matrix(
          seq(0.1, 1, length.out = input$data$n_ages),
          input$data$n_years_model, input$data$n_ages, byrow = TRUE)
        input$data$waa_pointer_jan1 <- count + 1L
        if (model_id == "mortality_weight_ar1") {
          input$par$M_repars <- c(log(0.1), 0, atanh(0.3))
          input$map$M_repars <- factor(rep(NA_integer_, 3L))
        }
      }
    } else if (model_id %in% names(cases)) {
      prepare_call <- baseline_call
      prepare_call$NAA_re <- cases[[model_id]]
      input <- eval(prepare_call, example)
    } else {
      input <- example[[sub("m", "input", model_id)]]
    }
    message("Generating ", tag, " ", model_id)
    fit <- historical$fit_wham(
      input,
      do.retro = FALSE,
      do.osa = FALSE,
      do.sdrep = FALSE,
      MakeADFun.silent = TRUE
    )
    fit$release_test <- list(tag = tag, commit = commit, case = model_id)
    saveRDS(fit, file.path(release_dir, paste0(model_id, "_original.rds")))
    save_reference_report(fit, file.path(release_dir, paste0(model_id, "_reference.rds")))
    archived_nll <- if (model_id %in% names(reference$nll)) reference$nll[[model_id]] else NA_real_
    summary <- rbind(summary, data.frame(
      tag = tag, commit = commit, model = model_id,
      nll = as.numeric(fit$opt$objective),
      archived_nll = archived_nll,
      difference = as.numeric(fit$opt$objective) - archived_nll,
      convergence = fit$opt$convergence,
      n_NAA_sigma = if (is.null(input$data$n_NAA_sigma)) NA_integer_ else input$data$n_NAA_sigma,
      max_gradient = max(abs(fit$final_gradient))
    ))
    write.csv(summary, file.path(release_dir, "generation_summary.csv"), row.names = FALSE)
    TMB::FreeADFun(fit)
    }, error = function(error) {
      message("Generation failed for ", model_id, ": ", conditionMessage(error))
      generation_errors <<- rbind(generation_errors,
        data.frame(model = model_id, error = conditionMessage(error)))
    })
  }
  write.csv(generation_errors, file.path(release_dir, "generation_errors.csv"), row.names = FALSE)
  print(summary)
  checked <- !is.na(summary$archived_nll)
  stopifnot(nrow(generation_errors) == 0L, all(is.finite(summary$nll)), sum(checked) == if (stage == "generate") 4L else 0L,
    all(abs(summary$difference[checked]) <= 1e-6))
} else if (stage %in% c("refit", "refit-defaults", "refit-features")) {
  current_library <- Sys.getenv("WHAM_CURRENT_LIBRARY")
  if (nzchar(current_library)) library(wham, lib.loc = current_library) else library(wham)
  updated <- new.env(parent = asNamespace("wham"))
  sys.source("R/refit_legacy_compatibility.R", updated)
  sys.source("R/migrate_wham_input.R", updated)
  sys.source("R/refit_wham.R", updated)
  sys.source("R/fit_tmb.R", updated)
  sys.source("R/fit_wham.R", updated)
  paths <- list.files(release_dir, "_original\\.rds$", full.names = TRUE)
  if (!length(paths)) stop("Generate historical fits first.")
  summary <- data.frame()
  for (path in paths) {
    model_id <- sub("_original\\.rds$", "", basename(path))
    message("Refitting ", tag, " ", model_id)
    warnings <- character()
    result <- data.frame(
      tag = tag, model = model_id, current_version = as.character(packageVersion("wham")),
      nll_original = NA_real_, nll_refit = NA_real_, nll_difference = NA_real_,
      nll_at_stored = NA_real_, stored_nll_difference = NA_real_,
      convergence_original = NA_integer_, convergence = NA_integer_,
      max_gradient_original = NA_real_, max_gradient_refit = NA_real_,
      stored_nll_matches = FALSE, nll_matches = FALSE, warnings = "", error = ""
    )
    tryCatch(withCallingHandlers({
      fit <- readRDS(path)
      if (startsWith(model_id, "steepness_") &&
          !isTRUE(fit$input$data$use_steepness == 1L)) {
        stop("Historical fixture does not have steepness enabled; regenerate it.")
      }
      result$nll_original <- as.numeric(fit$opt$objective)
      result$convergence_original <- fit$opt$convergence
      result$max_gradient_original <- max(abs(fit$final_gradient))
      legacy <- is.null(fit$input$data$n_stocks) &&
        !is.null(fit$input$data$n_NAA_sigma)
      if (legacy) {
        input <- suppressMessages(updated$migrate_legacy_refit_input(fit))
      } else {
        prepare <- updated$refit_wham
        environment(prepare) <- list2env(list(
          fit_wham = function(input, ...) stop(structure(
            list(message = "Input prepared", input = input),
            class = c("wham_prepared_input", "error", "condition")))
        ), parent = updated)
        input <- tryCatch(suppressMessages(prepare(fit)),
          wham_prepared_input = function(condition) condition$input)
      }
      objective <- TMB::MakeADFun(input$data, input$par, random = input$random,
        map = input$map, DLL = "wham", silent = TRUE)
      result$nll_at_stored <- as.numeric(objective$fn(objective$par))
      TMB::FreeADFun(objective)
      result$stored_nll_difference <- result$nll_at_stored - result$nll_original
      result$stored_nll_matches <- is.finite(result$stored_nll_difference) &&
        abs(result$stored_nll_difference) <= 1e-6
      if (stage != "refit-defaults") {
        fit$call <- quote(fit_wham(input = input, do.retro = FALSE, do.osa = FALSE,
          do.sdrep = FALSE, do.brps = FALSE, do.post.samp = FALSE, MakeADFun.silent = TRUE))
      }
      output <- capture.output(refit <- suppressMessages(updated$refit_wham(fit)))
      writeLines(output, file.path(release_dir, paste0(model_id, "_", stage, ".txt")))
      result$nll_refit <- as.numeric(refit$opt$objective)
      result$nll_difference <- result$nll_refit - result$nll_original
      result$convergence <- refit$opt$convergence
      result$max_gradient_refit <- max(abs(refit$final_gradient))
      result$nll_matches <- is.finite(result$nll_difference) && abs(result$nll_difference) <= 1e-6
      saveRDS(refit, file.path(release_dir, paste0(model_id, "_", stage, ".rds")))
      TMB::FreeADFun(refit)
    }, warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }), error = function(error) {
      result$error <<- conditionMessage(error)
    })
    result$warnings <- paste(unique(warnings), collapse = "\n")
    summary <- rbind(summary, result)
    write.csv(summary, file.path(release_dir, paste0(stage, "_summary.csv")), row.names = FALSE)
  }
  print(summary[, setdiff(names(summary), "warnings")])
  stopifnot(all(summary$stored_nll_matches), all(summary$nll_matches),
    all(summary$convergence == 0L), all(summary$error == ""))
} else {
  stop("Unknown stage: ", stage)
}