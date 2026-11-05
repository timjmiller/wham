refit_correction_warning <- function(issue, confirmed = TRUE) {
  warning(
    "WHAM refit correction: ", issue, ". ",
    if (confirmed) "The original fit used an affected historical calculation. " else
      "The stored report differs from the current calculation; the cause is not uniquely identified. ",
    "Refitting uses the current implementation; the corrected NLL and parameter estimates may differ from the original fit.",
    call. = FALSE
  )
}

# Recover historical process scales only when stored likelihood components
# identify a reproducible transform. Return input plus comparison-scale flags;
# reject ambiguous fixed/shared AR1 scales rather than change the model.
refit_process_scales <- function(input, stored_report) {
  if (is.null(stored_report$nll_sel) && is.null(stored_report$nll_Ecov)) {
    return(list(input = input, selectivity = FALSE, ecov = integer()))
  }
  report_at_parameters <- function(input) {
    validate_input_safety(input)
    model <- TMB::MakeADFun(
      input$data,
      input$par,
      map = input$map,
      DLL = "wham",
      silent = TRUE
    )
    on.exit(TMB::FreeADFun(model))
    model$fn(model$par)
    model$report()
  }
  matches <- function(x, y) {
    if (is.null(x)) {
      return(all(is.finite(y) & y == 0))
    }
    isTRUE(all.equal(x, y, tolerance = 1e-8))
  }
  current <- report_at_parameters(input)
  selectivity <- FALSE
  ecov <- integer()
  if (
    !is.null(stored_report$nll_sel) &&
      !matches(current$nll_sel, stored_report$nll_sel)
  ) {
    candidate <- input
    candidate$par$sel_repars[, 2:3] <- candidate$par$sel_repars[, 2:3] * 2
    if (
      !matches(report_at_parameters(candidate)$nll_sel, stored_report$nll_sel)
    ) {
      stop("Cannot safely reproduce stored selectivity process likelihood.")
    }
    input <- candidate
    selectivity <- TRUE
  }
  if (
    !is.null(stored_report$nll_Ecov) &&
      !matches(current$nll_Ecov, stored_report$nll_Ecov)
  ) {
    ecov <- which(input$data$Ecov_model == 2L & input$data$Ecov_use_re == 1L)
    candidate <- input
    rho <- tanh(input$par$Ecov_process_pars[3, ecov] / 2)
    candidate$par$Ecov_process_pars[2, ecov] <- input$par$Ecov_process_pars[
      2,
      ecov
    ] +
      0.5 * log1p(-rho^2)
    if (
      !length(ecov) ||
        !matches(
          report_at_parameters(candidate)$nll_Ecov,
          stored_report$nll_Ecov
        )
    ) {
      stop("Cannot safely reproduce stored environmental process likelihood.")
    }
    mapping <- input$map$Ecov_process_pars
    if (!is.null(mapping)) {
      mapping <- matrix(
        as.character(mapping),
        nrow(input$par$Ecov_process_pars)
      )
      sigma <- mapping[2, ecov]
      varying_rho <- !is.na(mapping[3, ecov])
      if (
        any(
          varying_rho &
            (is.na(sigma) |
              duplicated(sigma) |
              duplicated(sigma, fromLast = TRUE))
        )
      ) {
        stop(
          "Cannot safely migrate fixed or shared historical environmental AR1 scales."
        )
      }
    }
    input <- candidate
  }
  list(input = input, selectivity = selectivity, ecov = ecov)
}

# Compare equilibrium abundance at stored latent values, not optimized modes.
# Historical bug replication is no longer supported.
refit_initial_equilibrium <- function(input, stored_N1) {
  historical <- which(input$data$N1_model == -1L)
  if (length(historical)) {
    input$data$N1_model[historical] <- 1L
  }
  stocks <- which(input$data$N1_model == 1L)
  if (!length(stocks)) {
    return(input)
  }
  warn_equilibrium <- function(stocks, reason) {
    warning(
      "Initial NAA equilibrium for stock(s) ",
      paste(stocks, collapse = ", "),
      ": ",
      reason,
      " The model will be refit using the corrected equilibrium ",
      "calculation; parameter estimates and likelihoods may change.",
      call. = FALSE
    )
  }
  if (length(historical)) {
    warn_equilibrium(
      historical,
      "the saved input used the historical equilibrium bug."
    )
  }
  stocks <- setdiff(stocks, historical)
  if (!length(stocks)) {
    return(input)
  }
  if (is.null(stored_N1)) {
    warn_equilibrium(
      stocks,
      "stored initial abundance is unavailable, so exposure to the historical equilibrium bug cannot be checked."
    )
    return(input)
  }
  validate_input_safety(input)
  initial_numbers <- function(input) {
    # Compare reports at stored latent values, not newly optimized modes.
    model <- TMB::MakeADFun(
      input$data,
      input$par,
      map = input$map,
      DLL = "wham",
      silent = TRUE
    )
    on.exit(TMB::FreeADFun(model))
    model$fn(model$par)
    model$report()$N1
  }
  current <- initial_numbers(input)
  if (!identical(dim(current), dim(stored_N1))) {
    warn_equilibrium(
      stocks,
      "stored initial abundance has incompatible dimensions, so exposure to the historical equilibrium bug cannot be checked."
    )
    return(input)
  }
  affected <- integer()
  for (stock in stocks) {
    target <- stored_N1[stock, , , drop = FALSE]
    if (
      isTRUE(all.equal(
        current[stock, , , drop = FALSE],
        target,
        tolerance = 1e-8
      ))
    ) {
      next
    }
    affected <- c(affected, stock)
  }
  if (length(affected)) {
    warn_equilibrium(
      affected,
      "stored initial abundance differs from the corrected calculation and may reflect the historical equilibrium bug."
    )
  }
  input
}

# Recover recruitment decoupling from explicit options or flags first, then
# active correlation mappings. Default to coupled when evidence is absent so
# migration does not introduce an unsupported recruitment correlation.
refit_decouple_recruitment <- function(fit) {
  input <- fit$input
  options <- input$options
  NAA_re <- options[["NAA_re"]]
  if (is.null(NAA_re)) {
    NAA_re <- options[["NAA"]]
  }
  flags <- list(
    NAA_re$decouple_recruitment,
    input$data[["decouple_recruitment", exact = TRUE]],
    input$data[["decouple_rec", exact = TRUE]]
  )
  for (flag in flags) {
    if (length(flag) == 1L && !is.na(flag)) return(isTRUE(as.logical(flag)))
  }

  mapping <- input$map$trans_NAA_rho
  if (!is.null(mapping)) {
    dimensions <- dim(input$par$trans_NAA_rho)
    if (is.null(dimensions)) {
      dimensions <- c(1L, 1L, length(mapping))
    }
    if (tail(dimensions, 1L) == 3L) {
      active <- array(!is.na(mapping), dim = c(prod(head(dimensions, -1L)), 3L))
      return(any(active[, 3L]))
    }
    return(FALSE)
  }
  n_stocks <- input$data$n_stocks
  if (is.null(n_stocks)) {
    n_stocks <- 1L
  }
  n_stocks == 1L && sum(names(fit$opt$par) == "trans_NAA_rho") == 3L
}

# Migrate supported single-stock pre-2.0 fits, retaining fitted values, maps,
# observations, and correlation conventions. Reject unsupported structures or
# dimension mismatches instead of silently rebuilding a different model.
refit_legacy_age_comp_indices <- function(data, parameter) {
  counts <- if (parameter == "catch_paa_pars") {
    data$n_age_comp_pars_fleets
  } else {
    data$n_age_comp_pars_indices
  }
  if (is.null(counts) || anyNA(counts) || any(!counts %in% 0:3)) {
    stop("Cannot safely migrate packed age-composition parameters: ", parameter)
  }
  unlist(lapply(seq_along(counts), function(series) {
    series + length(counts) * (seq_len(counts[series]) - 1L)
  }), use.names = FALSE)
}

refit_legacy_age_comp_weights <- function(input, stored_report) {
  if (!is.list(stored_report)) return(input)
  report_at_parameters <- function(candidate) {
    validate_input_safety(candidate)
    model <- TMB::MakeADFun(candidate$data, candidate$par, map = candidate$map,
      DLL = "wham", silent = TRUE)
    on.exit(TMB::FreeADFun(model))
    model$fn(model$par)
    model$report()
  }
  weighted_report <- NULL
  for (kind in c("catch", "index")) {
    group <- if (kind == "catch") "fleets" else "indices"
    models <- input$data[[paste0("age_comp_model_", group)]]
    use <- input$data[[paste0("use_", kind, "_paa")]]
    likelihood <- paste0("nll_", kind, "_acomp")
    target <- stored_report[[likelihood]]
    if (!is.matrix(target) || !identical(dim(target), dim(use))) next
    for (series in which(models == 5L)) {
      active <- which(use[, series] == 1L)
      if (!length(active)) next
      if (is.null(weighted_report)) {
        weighted_report <- tryCatch(report_at_parameters(input), error = function(e) NULL)
        if (is.null(weighted_report)) return(input)
      }
      if (!is.matrix(weighted_report[[likelihood]]) ||
          !identical(dim(weighted_report[[likelihood]]), dim(use))) next
      matches <- isTRUE(all.equal(
        as.numeric(weighted_report[[likelihood]][active, series]),
        as.numeric(target[active, series]), tolerance = 1e-8, scale = 1
      ))
      if (!matches) {
        refit_correction_warning(
          "historical logistic-normal effective-sample-size weighting",
          confirmed = FALSE
        )
        return(input)
      }
    }
  }
  input
}

migrate_refit_selectivity <- function(input) {
  data <- input$data
  pars <- input$par
  width <- max(data$n_ages, 4L)
  shape <- c(data$n_selblocks, data$n_years_model, width)
  deviations <- array(0, shape)
  deviation_map <- array(NA_character_, shape)
  indices <- matrix(0L, data$n_selblocks, width)
  counts <- rep(1L, data$n_selblocks)
  old_map <- input$map$selpars_re
  if (is.null(old_map)) old_map <- factor(seq_along(pars$selpars_re))
  if (length(old_map) != length(pars$selpars_re)) {
    stop("Invalid stored selectivity random-effect dimensions.")
  }
  if (all(data$selblock_models_re == 1L) && length(pars$selpars_re) == 1L && all(is.na(old_map))) {
    pars$selpars_re <- numeric()
    old_map <- factor(character(), levels = levels(old_map))
  }
  array_layout <- !is.null(dim(pars$selpars_re))
  if (array_layout && !identical(as.integer(dim(pars$selpars_re)), as.integer(shape))) {
    stop("Invalid stored selectivity random-effect dimensions.")
  }
  array_map <- if (array_layout) array(as.character(old_map), shape) else NULL
  used <- array(FALSE, shape)
  offset <- 0L
  for (block in seq_len(data$n_selblocks)) {
    model <- data$selblock_models[block]
    columns <- if (model == 1L) seq_len(data$n_ages) else if (model == 3L) {
      data$n_ages + 2L + seq_len(4L)
    } else data$n_ages + seq_len(2L)
    if (data$selblock_models_re[block] > 1L) {
      years <- which(data$selblock_years[, block] == 1L)
      active <- which(data$selpars_est[block, columns] == 1L)
      count <- data$n_selpars_est[block]
      if (length(active) != count || !length(years) || count < 1L || count > width) {
        stop("Cannot safely migrate stored selectivity random effects for block ", block, ".")
      }
      if (array_layout) {
        values <- matrix(pars$selpars_re[block, years, seq_len(count)], nrow = length(years))
        mapping <- matrix(array_map[block, years, seq_len(count)], nrow = length(years))
        used[block, years, seq_len(count)] <- TRUE
      } else {
        positions <- offset + seq_len(length(years) * count)
        if (max(positions) > length(pars$selpars_re)) stop("Invalid stored selectivity random-effect dimensions.")
        values <- matrix(pars$selpars_re[positions], nrow = length(years))
        mapping <- matrix(as.character(old_map[positions]), nrow = length(years))
        offset <- offset + length(positions)
      }
      if (data$selblock_models_re[block] == 4L) {
        if (!identical(as.vector(mapping), rep(mapping[, 1L], count)) ||
            !isTRUE(all.equal(as.vector(values), rep(values[, 1L], count)))) {
          stop("Cannot safely migrate unshared annual selectivity random effects.")
        }
        counts[block] <- 1L
        indices[block, active] <- 1L
        deviations[block, years, 1L] <- values[, 1L]
        deviation_map[block, years, 1L] <- mapping[, 1L]
      } else {
        counts[block] <- count
        indices[block, active] <- seq_len(count)
        deviations[block, years, seq_len(count)] <- values
        deviation_map[block, years, seq_len(count)] <- mapping
      }
    }
  }
  if ((!array_layout && offset != length(pars$selpars_re)) ||
      (array_layout && any(!is.na(array_map[!used])))) {
    stop("Invalid stored selectivity random-effect dimensions.")
  }
  input$par <- pars
  input$par$selpars_re <- deviations
  input$map$selpars_re <- factor(deviation_map, levels = levels(old_map))
  input$data$selpars_re_index <- indices
  input$data$n_selpars_re <- counts
  input$data$selblock_models_re <- as.integer(data$selblock_models_re - 1L)
  input
}

refit_legacy_ecov_layout <- function(data, beta) {
  dimensions <- dim(beta)
  expanded <- length(dimensions) == 4L
  age_specific <- length(dimensions) == 3L
  if (!expanded && !age_specific && length(dimensions) != 2L) stop("Invalid legacy environmental coefficient dimensions.")
  degree <- dimensions[if (expanded) 2L else 1L]
  covariates <- data$n_Ecov
  effects <- 2L + data$n_indices
  if ((expanded && !identical(as.integer(dimensions), as.integer(c(effects, degree, covariates, data$n_ages)))) ||
      (!expanded && dimensions[2L] != covariates) ||
      (age_specific && dimensions[3L] != data$n_ages)) stop("Invalid legacy environmental coefficient dimensions.")
  source <- array(seq_along(beta), dimensions)
  where <- if (expanded) data$Ecov_where else {
    result <- matrix(0L, covariates, effects)
    result[, 1L] <- data$Ecov_where == 1L
    result[, 2L] <- data$Ecov_where == 2L
    result
  }
  if (!identical(as.integer(dim(where)), as.integer(c(covariates, effects)))) {
    stop("Invalid legacy environmental effect destinations.")
  }
  if (expanded) where[data$Ecov_how == 0L, 1L] <- 0L else where[data$Ecov_how == 0L, ] <- 0L
  indices <- list(
    Ecov_beta_R = array(0L, c(1L, covariates, degree)),
    Ecov_beta_M = array(0L, c(1L, data$n_ages, 1L, covariates, degree)),
    Ecov_beta_q = array(0L, c(data$n_indices, covariates, degree))
  )
  for (covariate in seq_len(covariates)) for (power in seq_len(degree)) {
    if (where[covariate, 1L] != 0L) indices$Ecov_beta_R[1L, covariate, power] <-
      if (expanded) source[1L, power, covariate, 1L] else if (age_specific) source[power, covariate, 1L] else source[power, covariate]
    if (where[covariate, 2L] != 0L) indices$Ecov_beta_M[1L, , 1L, covariate, power] <-
      if (expanded) source[2L, power, covariate, ] else if (age_specific) source[power, covariate, ] else source[power, covariate]
    for (index in seq_len(data$n_indices)) if (where[covariate, 2L + index] != 0L) {
      indices$Ecov_beta_q[index, covariate, power] <- source[2L + index, power, covariate, 1L]
    }
  }
  start <- if (expanded) data$ind_Ecov_out_start else matrix(data$ind_Ecov_out_start, covariates, effects)
  end <- if (expanded) data$ind_Ecov_out_end else matrix(data$ind_Ecov_out_end, covariates, effects)
  if (anyNA(start) || anyNA(end) || any(start < 0L) || any(end >= data$n_years_Ecov) ||
      any(end - start + 1L != data$n_years_model)) stop("Invalid stored environmental lag indices.")
  list(indices = indices, where = where, degree = degree, start = start, end = end)
}

refit_legacy_ecov_process <- function(input, stored_report) {
  parameters <- input$par$Ecov_process_pars
  mapping <- matrix(as.character(input$map$Ecov_process_pars), nrow(parameters))
  likelihood <- function(parameters, deviations, historical) {
    correlation <- tanh(parameters[3L] / 2)
    sigma <- exp(parameters[2L])
    conditional_mean <- correlation * (head(deviations, -1L) - if (historical) parameters[1L] else 0)
    -dnorm(deviations[1L], 0, sigma / sqrt(1 - correlation^2), log = TRUE) -
      sum(dnorm(deviations[-1L], conditional_mean, sigma, log = TRUE))
  }
  legacy <- integer()
  for (covariate in which(input$data$Ecov_model == 2L & input$data$Ecov_use_re == 1L)) {
    deviations <- input$par$Ecov_re[, covariate]
    mean <- parameters[1L, covariate]
    correlation <- tanh(parameters[3L, covariate] / 2)
    if ((mean == 0 && is.na(mapping[1L, covariate])) ||
        (correlation == 0 && is.na(mapping[3L, covariate]))) next
    target <- stored_report$nll_Ecov
    if (!is.matrix(target) || ncol(target) < covariate) {
      stop("Cannot infer historical environmental AR1 without stored process likelihoods.")
    }
    target <- sum(target[, covariate])
    reported_parameters <- parameters[, covariate]
    if (identical(dim(stored_report$Ecov_process_pars), dim(parameters)) &&
        identical(dim(stored_report$Ecov_x), dim(input$par$Ecov_re))) {
      reported_parameters <- stored_report$Ecov_process_pars[, covariate]
      deviations <- stored_report$Ecov_x[, covariate] - reported_parameters[1L]
    }
    standard <- likelihood(reported_parameters, deviations, FALSE)
    historical <- likelihood(reported_parameters, deviations, TRUE)
    matches <- vapply(c(standard, historical), function(value) {
      is.finite(value) && is.finite(target) && abs(value - target) <= 1e-7
    }, logical(1))
    if (sum(matches) != 1L) stop("Cannot safely infer historical environmental AR1 convention for covariate ", covariate, ".")
    if (matches[2L]) legacy <- c(legacy, covariate)
  }
  if (length(legacy)) {
    warning("Historical environmental AR1 is no longer supported; refitting with the current AR1 model. Likelihoods and estimates may change.", call. = FALSE)
  }
  input
}

migrate_legacy_refit_input <- function(fit) {
  old <- fit$input
  data <- old$data
  pars <- old$par
  if (!is.null(fit$parList)) pars[names(fit$parList)] <- fit$parList
  environmental <- any(data$Ecov_model != 0L)
  packed_age_comp <- !is.matrix(old$par$catch_paa_pars) &&
    !is.matrix(old$par$index_paa_pars)
  data$decouple_rec <- as.integer(refit_decouple_recruitment(fit))
  unsupported <- character()
  if (any(!data$Ecov_model %in% 0:2)) {
    unsupported <- c(unsupported, "unknown environmental process")
  }
  if (any(!data$selblock_models_re %in% 1:5)) {
    unsupported <- c(unsupported, "unknown selectivity random-effect model")
  }
  if (length(data$M_re_model) != 1L || !data$M_re_model %in% 1:5) {
    unsupported <- c(unsupported, "unknown mortality random-effect model")
  }
  if (!data$M_model %in% 1:3) {
    unsupported <- c(unsupported, "unknown mortality mean model")
  }
  if (isTRUE(data$use_steepness == 1)) {
    unsupported <- c(unsupported, "steepness parameterization")
  }
  if (isTRUE(data$do_proj == 1)) {
    unsupported <- c(unsupported, "stored projections")
  }
  if (is.null(data$agesvec) && packed_age_comp && any(!c(
    data$age_comp_model_fleets, data$age_comp_model_indices
  ) %in% c(1L, 2L, 7L))) {
    unsupported <- c(unsupported, "unsupported early age-composition likelihood")
  }
  if (length(unsupported)) {
    stop(
      "Cannot safely migrate this legacy WHAM input: ",
      paste(unsupported, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (is.null(old$asap3$n_ages)) {
    years <- old[["years"]]
    if (length(years) != data$n_years_model || anyNA(years)) {
      stop("Cannot safely migrate legacy input without stored model years.")
    }
    input <- prepare_wham_input(
      basic_info = list(
        years = as.integer(years),
        ages = seq_len(data$n_ages)
      ),
      catch_info = list(n_fleets = data$n_fleets),
      index_info = list(n_indices = data$n_indices),
      selectivity = list(n_selblocks = data$n_selblocks),
      model_name = old$model_name,
      recruit_model = data$recruit_model,
      NAA_re = list(decouple_recruitment = as.logical(data$decouple_rec))
    )
    input$years <- input$years_full <- years
    input$ages.lab <- old$ages.lab
  } else {
    asap3 <- align_fitted_asap3_indices(old$asap3)
    input <- prepare_wham_input(
      asap3 = list(list(dat = asap3)),
      model_name = old$model_name,
      recruit_model = data$recruit_model,
      NAA_re = list(decouple_recruitment = as.logical(data$decouple_rec))
    )
  }
  if (environmental) {
    if (any(!data$Ecov_obs_sigma_opt %in% 1:4)) stop("Unsupported legacy environmental observation model.")
    if (is.null(pars$Ecov_obs_logsigma_re) && any(data$Ecov_obs_sigma_opt == 4L)) {
      stop("Cannot safely migrate early environmental observation-error random effects.")
    }
    years <- data$year1_Ecov + seq_len(data$n_years_Ecov) - 1L
    input$by_pwi <- TRUE
    input <- set_ecov(input, list(
      mean = data$Ecov_obs, use_obs = data$Ecov_use_obs, year = years,
      logsigma = pars$Ecov_obs_logsigma,
      label = as.character(unlist(data$Ecov_label)),
      process_model = c(NA_character_, "rw", "ar1")[data$Ecov_model + 1L]
    ))
    input$by_pwi <- NULL
    input$years_Ecov <- years
    input$data$n_years_Ecov <- data$n_years_Ecov
    input$data$Ecov_obs <- data$Ecov_obs
    input$data$Ecov_use_obs <- data$Ecov_use_obs
    for (parameter in c("Ecov_re", "Ecov_obs_logsigma", "Ecov_obs_logsigma_re")) {
      input$par[[parameter]] <- matrix(0, data$n_years_Ecov, data$n_Ecov)
      input$map[[parameter]] <- factor(rep(NA_integer_, data$n_years_Ecov * data$n_Ecov))
    }
  }
  for (field in intersect(names(data), names(input$data))) {
    value <- data[[field]]
    target <- input$data[[field]]
    if (is.null(dim(target)) || is.data.frame(target) || field == "waa") {
      input$data[[field]] <- value
    } else if (length(value) == length(target)) {
      target[] <- value
      input$data[[field]] <- target
    }
  }
  input$data$waa_pointer_indices <- data$waa_pointer_indices[seq_len(
    data$n_indices
  )]
  input$data$waa_pointer_M <- data$waa_pointer_jan1
  input <- migrate_legacy_refit_model_codes(input, data)
  input$data$M_re_model[] <- c(1L, 4L, 2L, 3L, 4L)[data$M_re_model]
  input$data$n_M_re[] <- data$n_M_a
  input$data$M_re_index[] <- if (data$n_M_a == 1L) 1L else seq_len(data$n_ages)
  input$data$selblock_models_re <- rep(0L, data$n_selblocks)
  input$data$NAA_re_model <- min(data$n_NAA_sigma, 2L)
  input$data$decouple_recruitment <- data$decouple_rec
  input$data$use_alt_AR1 <- 0L
  if (length(data$which_F_age) == 1L) {
    input$data$which_F_age <- rep(data$which_F_age, data$n_years_model)
  }
  if (!is.null(data$static_FXSPR_init)) {
    input$data$FXSPR_static_init <- data$static_FXSPR_init
  }
  if (!is.null(data$avg_years_ind)) {
    for (field in grep("^avg_years_ind_", names(input$data), value = TRUE)) {
      target <- input$data[[field]]
      if (length(dim(target)) == 2) {
        target[,] <- 0
        target[seq_len(length(data$avg_years_ind) + 1L), ] <- c(
          length(data$avg_years_ind),
          data$avg_years_ind
        )
      } else if (length(dim(target)) == 3) {
        target[] <- 0
        target[,, seq_len(length(data$avg_years_ind) + 1L)] <- rep(
          c(length(data$avg_years_ind), data$avg_years_ind),
          each = prod(dim(target)[1:2])
        )
      }
      input$data[[field]] <- target
    }
  }

  selection <- old
  selection$par <- pars
  selection <- migrate_refit_selectivity(selection)
  input$par$selpars_re <- selection$par$selpars_re
  input$map$selpars_re <- selection$map$selpars_re
  for (field in c("selpars_re_index", "n_selpars_re", "selblock_models_re")) {
    input$data[[field]] <- selection$data[[field]]
  }
  if (!is.null(pars$M0) && data$M_model == 2L &&
      (is.null(old$map$M0) || any(!is.na(old$map$M0)))) {
    stop("Cannot safely migrate an estimated legacy M0 with age-specific offsets.")
  }
  transfer <- function(
    target,
    source = target,
    indices = seq_along(pars[[source]]),
    scale = 1
  ) {
    values <- pars[[source]][indices] * scale
    if (length(values) != length(input$par[[target]])) {
      stop(
        "Cannot safely migrate legacy parameter dimensions: ",
        source,
        call. = FALSE
      )
    }
    input$par[[target]][] <<- values
    mapping <- old$map[[source]]
    if (is.null(mapping)) {
      mapping <- factor(seq_along(pars[[source]]))
    }
    input$map[target] <<- list(factor(
      as.character(mapping[indices]),
      levels = levels(mapping)
    ))
  }
  renamed <- c(
    "log_N1_pars",
    "log_NAA_sigma",
    "mean_rec_pars",
    "F_devs",
    "log_F1",
    "M_a",
    "M_re",
    "trans_NAA_rho",
    "catch_paa_pars",
    "index_paa_pars"
  )
  inactive <- c("selpars_re", "Ecov_beta", "Ecov_obs_sigma_par")
  for (parameter in setdiff(names(pars), c(renamed, inactive))) {
    if (parameter %in% names(input$par)) transfer(parameter)
  }
  if (data$M_model == 3L) {
    input <- migrate_legacy_weight_mortality(input, old, pars, data)
  } else {
    transfer("M_re", indices = if (data$n_M_a == 1L) {
      rep(seq_len(data$n_years_model), data$n_ages)
    } else seq_along(pars$M_re))
  }
  if (environmental) {
    use_process <- data$Ecov_use_re
    if (is.matrix(use_process)) {
      if (any(vapply(seq_len(ncol(use_process)), function(covariate) {
        length(unique(use_process[, covariate])) != 1L
      }, logical(1)))) stop("Cannot safely migrate year-specific environmental process masking.")
      use_process <- use_process[1L, ]
    }
    input$data$Ecov_use_re <- as.integer(use_process)
    input$data$years_use_Ecov <- seq_len(data$n_years_Ecov) - 1L
    if (length(pars$Ecov_obs_sigma_par) == length(input$par$Ecov_obs_sigma_par)) {
      transfer("Ecov_obs_sigma_par")
    }
    layout <- refit_legacy_ecov_layout(data, pars$Ecov_beta)
    source_map <- old$map$Ecov_beta
    if (is.null(source_map)) source_map <- factor(seq_along(pars$Ecov_beta))
    assigned <- character()
    for (parameter in names(layout$indices)) {
      positions <- layout$indices[[parameter]]
      active <- positions > 0L
      values <- array(0, dim(positions))
      mapping <- array(NA_character_, dim(positions))
      values[active] <- pars$Ecov_beta[positions[active]]
      mapping[active] <- as.character(source_map[positions[active]])
      labels <- unique(mapping[!is.na(mapping)])
      if (length(intersect(assigned, labels))) {
        stop("Cannot safely migrate environmental coefficients shared across effect types.")
      }
      assigned <- c(assigned, labels)
      input$par[[parameter]] <- values
      input$map[[parameter]] <- factor(mapping, levels = levels(source_map)[levels(source_map) %in% labels])
    }
    if (length(setdiff(levels(source_map), assigned))) stop("Cannot safely migrate unused estimated environmental coefficients.")
    input$data$Ecov_how_R[, 1L] <- layout$where[, 1L] * data$Ecov_how
    input$data$Ecov_how_q[,] <- layout$where[, 2L + seq_len(data$n_indices), drop = FALSE]
    input$data$Ecov_how_M[, 1L, , 1L] <- rep(layout$where[, 2L], data$n_ages)
    for (suffix in c("R", "M", "q")) {
      input$data[[paste0("n_poly_Ecov_", suffix)]][] <- layout$degree
      for (boundary in c("start", "end")) {
        pointer <- layout[[boundary]]
        target <- paste0("ind_Ecov_out_", boundary, "_", suffix)
        if (suffix == "R") input$data[[target]][, 1L] <- pointer[, 1L]
        if (suffix == "M") input$data[[target]][, 1L, , 1L] <- rep(pointer[, 2L], data$n_ages)
        if (suffix == "q") input$data[[target]][, ] <- pointer[, 2L + seq_len(data$n_indices), drop = FALSE]
      }
    }
  }
  for (parameter in c("catch_paa_pars", "index_paa_pars")) {
    if (is.matrix(pars[[parameter]])) {
      transfer(parameter)
    } else {
      positions <- refit_legacy_age_comp_indices(data, parameter)
      if (length(positions) != length(pars[[parameter]])) {
        stop("Cannot safely migrate legacy parameter dimensions: ", parameter)
      }
      input$par[[parameter]][] <- 0
      input$par[[parameter]][positions] <- pars[[parameter]]
      mapping <- old$map[[parameter]]
      if (is.null(mapping)) mapping <- factor(seq_along(pars[[parameter]]))
      if (length(mapping) != length(positions)) {
        stop("Cannot safely migrate legacy parameter mapping: ", parameter)
      }
      target_map <- rep(NA_character_, length(input$par[[parameter]]))
      target_map[positions] <- as.character(mapping)
      input$map[[parameter]] <- factor(target_map, levels = levels(mapping))
    }
  }
  for (parameter in c("log_N1", "mean_rec_pars")) {
    source <- if (parameter == "log_N1") "log_N1_pars" else parameter
    values <- pars[[source]]
    input$par[[parameter]][seq_along(values)] <- values
    mapping <- rep(NA_integer_, length(input$par[[parameter]]))
    source_map <- old$map[[source]]
    if (is.null(source_map)) {
      source_map <- factor(seq_along(values))
    }
    mapping[seq_along(values)] <- as.integer(source_map)
    input$map[[parameter]] <- factor(mapping)
  }
  pointers <- data$NAA_sigma_pointers
  if (data$n_NAA_sigma > 0) {
    transfer("log_NAA_sigma", indices = pointers)
  }
  if (!is.null(pars$M0) && data$M_model %in% c(1L, 3L)) {
    transfer("Mpars", "M0", rep(1L, data$n_ages))
  } else {
    transfer(
      "Mpars",
      "M_a",
      if (data$M_model %in% c(1L, 3L)) rep(1L, data$n_ages) else seq_len(data$n_ages)
    )
    if (!is.null(pars$M0)) input$par$Mpars[] <- input$par$Mpars + pars$M0
  }
  if (length(pars$trans_NAA_rho) == 2L && data$decouple_rec == 0L) {
    input$par$trans_NAA_rho[] <- c(2 * pars$trans_NAA_rho, 0)
    mapping <- old$map$trans_NAA_rho
    if (is.null(mapping)) mapping <- factor(seq_len(2L))
    if (length(mapping) != 2L) stop("Invalid legacy abundance correlation mapping.")
    input$map$trans_NAA_rho <- factor(
      c(as.character(mapping), NA_character_), levels = levels(mapping)
    )
  } else {
    transfer("trans_NAA_rho", scale = 2)
  }
  if (data$decouple_rec == 1) {
    input$par$trans_NAA_rho[] <- 2 * pars$trans_NAA_rho[c(1L, 3L, 2L)]
    mapping <- input$map$trans_NAA_rho
    input$map$trans_NAA_rho <- factor(
      as.character(mapping[c(1L, 3L, 2L)]),
      levels = levels(mapping)
    )
    if (data$n_NAA_sigma > 1) {
      rho <- tanh(pars$trans_NAA_rho)
      if (rho[1]^2 != 0 || rho[2]^2 != rho[3]^2) {
        refit_correction_warning("historical decoupled recruitment used the wrong correlations in NAA marginal variances")
      }
    }
  }
  if (!is.null(pars$q_repars)) {
    input$par$q_repars[, 2] <- 2 * pars$q_repars[, 2]
  }
  input$par$sel_repars[, 2:3] <- 2 * pars$sel_repars[, 2:3]
  input$par$M_repars[1, 1, 2:3] <- 2 * pars$M_repars[2:3]
  input$par$F_pars[1, ] <- pars$log_F1
  input$par$F_pars[-1, ] <- pars$F_devs
  map_F <- matrix(NA_character_, data$n_years_model, data$n_fleets)
  for (parameter in c("log_F1", "F_devs")) {
    mapping <- old$map[[parameter]]
    if (is.null(mapping)) {
      mapping <- factor(seq_along(pars[[parameter]]))
    }
    values <- ifelse(
      is.na(mapping),
      NA_character_,
      paste(parameter, as.character(mapping))
    )
    if (parameter == "log_F1") map_F[1, ] <- values else map_F[-1, ] <- values
  }
  input$map$F_pars <- factor(
    map_F,
    levels = unique(as.vector(map_F)[!is.na(map_F)])
  )
  input$random <- old$random
  input$options <- NULL
  if (is.null(data$agesvec)) {
    catch_logistic <- packed_age_comp & data$age_comp_model_fleets == 7L
    index_logistic <- packed_age_comp & data$age_comp_model_indices == 7L
    if (any(catch_logistic) || any(index_logistic)) {
      input$data$age_comp_model_fleets[catch_logistic] <- 5L
      input$data$age_comp_model_indices[index_logistic] <- 5L
    }
    input <- set_osa_obs(input)
    if (any(catch_logistic) || any(index_logistic)) {
      input <- refit_legacy_age_comp_weights(input, fit[["rep"]])
    }
  }
  if (environmental && !is.null(data$agesvec)) input <- set_osa_obs(input)
  if (environmental) input <- refit_legacy_ecov_process(input, fit$rep)
  input
}
