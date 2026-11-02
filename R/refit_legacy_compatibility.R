#' Translate legacy mortality and recruitment model codes
#'
#' Convert legacy model codes to the current codes used by the WHAM objective.
#'
#' @param input A partially migrated WHAM input.
#' @param data The original legacy input data.
#'
#' @return The migrated input with current mortality and recruitment model codes.
#' @keywords internal
migrate_legacy_refit_model_codes <- function(input, data) {
  input$data$M_model <- if (data$M_model == 3L) 2L else 1L
  input
}

#' Convert legacy steepness parameters to current stock-recruit coefficients
#'
#' A legacy steepness model is equivalent to its current Beverton-Holt or
#' Ricker form when unfished SPR is constant over time. Reject fits where that
#' equivalence would depend on time-varying mortality or spawning inputs.
#'
#' @param input The partially migrated WHAM input.
#' @param fit The original fitted legacy model.
#' @param data The original legacy input data.
#' @param pars The fitted legacy parameter values.
#'
#' @return The migrated input with current stock-recruit coefficients and
#'   comparison metadata.
#' @keywords internal
migrate_legacy_steepness <- function(input, fit, data, pars) {
  if (!isTRUE(data$use_steepness == 1L)) {
    return(input)
  }
  if (
    input$data$n_stocks != 1L ||
      input$data$n_regions != 1L ||
      input$data$n_seasons != 1L ||
      isTRUE(data$do_proj == 1L)
  ) {
    stop(
      "Cannot safely migrate legacy steepness outside a single-stock, single-region, single-season model without projections.",
      call. = FALSE
    )
  }
  if (length(data$M_re_model) != 1L || data$M_re_model != 1L) {
    stop(
      "Cannot safely migrate legacy steepness with time-varying mortality random effects.",
      call. = FALSE
    )
  }
  if (
    !is.null(input$data$Ecov_how_M) &&
      (anyNA(input$data$Ecov_how_M) || any(input$data$Ecov_how_M != 0L))
  ) {
    stop(
      "Cannot safely migrate legacy steepness with environmental effects on mortality.",
      call. = FALSE
    )
  }

  years <- data$n_years_model
  ages <- data$n_ages
  mortality <- fit$rep$log_M
  expected_dims <- c(1L, 1L, years, ages)
  if (is.null(mortality) ||
      !identical(as.integer(dim(mortality)), expected_dims)) {
    stop(
      "Cannot safely migrate legacy steepness without stored annual mortality estimates.",
      call. = FALSE
    )
  }
  if (
    !identical(as.integer(dim(data$mature)), c(1L, years, ages)) ||
      !identical(
        as.integer(dim(input$data$waa))[-1L],
        c(years, ages)
      ) ||
      length(input$data$waa_pointer_ssb) != 1L ||
      input$data$waa_pointer_ssb < 1L ||
      input$data$waa_pointer_ssb > dim(input$data$waa)[1L] ||
      !is.matrix(input$data$fracyr_SSB) ||
      !identical(as.integer(dim(input$data$fracyr_SSB)), c(years, 1L))
  ) {
    stop(
      "Cannot safely migrate legacy steepness with incomplete spawning inputs.",
      call. = FALSE
    )
  }
  if (
    data$M_model == 3L &&
      (
        length(input$data$waa_pointer_M) != 1L ||
          input$data$waa_pointer_M < 1L ||
          input$data$waa_pointer_M > dim(input$data$waa)[1L]
      )
  ) {
    stop(
      "Cannot safely migrate legacy steepness with an invalid mortality weight pointer.",
      call. = FALSE
    )
  }
  maturity <- matrix(data$mature, nrow = years, ncol = ages)
  spawning_weights <- matrix(
    input$data$waa[input$data$waa_pointer_ssb[1L], , ],
    nrow = years,
    ncol = ages
  )
  spawn_fraction <- as.numeric(input$data$fracyr_SSB[, 1L])
  is_year_invariant <- function(value) {
    values <- matrix(value, nrow = dim(value)[1L])
    all(vapply(seq_len(ncol(values)), function(column) {
      isTRUE(all.equal(
        values[, column],
        rep(values[1L, column], nrow(values)),
        tolerance = 1e-10
      ))
    }, logical(1)))
  }
  if (
    !is_year_invariant(maturity) ||
      !is_year_invariant(spawning_weights) ||
      !isTRUE(all.equal(
        spawn_fraction,
        rep(spawn_fraction[1L], years),
        tolerance = 1e-10
      )) ||
      (data$M_model == 3L &&
        !is_year_invariant(input$data$waa[input$data$waa_pointer_M[1L], , ]))
  ) {
    stop(
      "Cannot safely migrate legacy steepness when spawning inputs or mortality vary by year.",
      call. = FALSE
    )
  }

  annual_mortality <- matrix(mortality[1L, 1L, , ], nrow = years, ncol = ages)
  if (!is_year_invariant(annual_mortality)) {
    stop(
      "Cannot safely migrate legacy steepness with time-varying mortality.",
      call. = FALSE
    )
  }
  spr <- vapply(seq_len(years), function(year) {
    mortality_at_age <- exp(annual_mortality[year, ])
    survival <- cumprod(c(1, exp(-mortality_at_age[-ages])))
    survival[ages] <- survival[ages] / (1 - exp(-mortality_at_age[ages]))
    sum(
      survival * maturity[year, ] * spawning_weights[year, ] *
        exp(-spawn_fraction[year] * mortality_at_age)
    )
  }, numeric(1))
  if (!all(is.finite(spr)) || spr[1L] <= 0 ||
      !isTRUE(all.equal(spr, rep(spr[1L], years), tolerance = 1e-10))) {
    stop(
      "Cannot safely migrate legacy steepness because unfished SPR varies by year.",
      call. = FALSE
    )
  }

  model <- data$recruit_model[1L]
  steepness <- if (model == 3L) {
    0.2 + 0.8 * plogis(pars$mean_rec_pars[1L, 1L])
  } else {
    0.2 + exp(pars$mean_rec_pars[1L, 1L])
  }
  r0 <- exp(pars$mean_rec_pars[1L, 2L])
  if (model == 3L) {
    alpha <- 4 * steepness / (spr[1L] * (1 - steepness))
    beta <- (5 * steepness - 1) /
      ((1 - steepness) * spr[1L] * r0)
  } else {
    alpha <- 1.25 * log(5 * steepness) / spr[1L]
    beta <- alpha / r0
  }
  if (!all(is.finite(c(alpha, beta))) || any(c(alpha, beta) <= 0)) {
    stop("Cannot safely convert legacy steepness parameters.", call. = FALSE)
  }
  input$data$recruit_model <- as.integer(model)
  input$par$mean_rec_pars[1L, ] <- log(c(alpha, beta))
  input$legacy_steepness <- list(recruit_model = model, spr = spr[1L])
  input
}

#' Restore legacy steepness parameters for refit comparisons
#'
#' Convert current stock-recruit coefficients back to the original steepness
#' and unfished recruitment parameterization.
#'
#' @param pars Parameters from the refitted model.
#' @param metadata The steepness metadata returned with the migrated input.
#'
#' @return The refitted parameters with legacy steepness values restored.
#' @keywords internal
refit_legacy_steepness_parameters <- function(pars, metadata) {
  if (is.null(metadata)) {
    return(pars)
  }
  alpha <- exp(pars$mean_rec_pars[1L, 1L])
  beta <- exp(pars$mean_rec_pars[1L, 2L])
  spr <- metadata$spr
  if (metadata$recruit_model == 3L) {
    steepness <- alpha * spr / (4 + alpha * spr)
    r0 <- (5 * steepness - 1) /
      ((1 - steepness) * spr * beta)
    pars$mean_rec_pars[1L, ] <- c(
      qlogis((steepness - 0.2) / 0.8),
      log(r0)
    )
  } else {
    steepness <- exp(alpha * spr / 1.25) / 5
    r0 <- alpha / beta
    pars$mean_rec_pars[1L, ] <- c(log(steepness - 0.2), log(r0))
  }
  pars
}

#' Migrate legacy weight-based mortality deviations
#'
#' Preserve fitted annual mortality deviations when the legacy model represented
#' one shared intercept process across ages.
#'
#' @param input A partially migrated WHAM input.
#' @param old The original fitted input.
#' @param pars The fitted legacy parameter values.
#' @param data The original legacy input data.
#'
#' @return The input with migrated mortality parameters and mappings.
#' @keywords internal
migrate_legacy_weight_mortality <- function(input, old, pars, data) {
  if (
    data$n_ages < 2L ||
      !is.matrix(pars$M_re) ||
      any(vapply(seq_len(ncol(pars$M_re)), function(age) {
        !isTRUE(all.equal(pars$M_re[, age], pars$M_re[, 1L]))
      }, logical(1)))
  ) {
    stop(
      "Cannot safely migrate age-varying deviations in weight-based mortality.",
      call. = FALSE
    )
  }
  mapping <- old$map$M_re
  if (is.null(mapping)) {
    mapping <- factor(seq_along(pars$M_re))
  }
  old_mapping <- matrix(as.character(mapping), nrow(pars$M_re))
  if (any(vapply(seq_len(ncol(old_mapping)), function(age) {
    !identical(old_mapping[, age], old_mapping[, 1L])
  }, logical(1)))) {
    stop(
      "Cannot safely migrate unshared weight-based mortality deviations.",
      call. = FALSE
    )
  }

  input$par$M_re[] <- 0
  input$par$M_re[1L, 1L, , 1L] <- pars$M_re[, 1L]
  target_map <- array(NA_character_, dim(input$par$M_re))
  target_map[1L, 1L, , 1L] <- old_mapping[, 1L]
  input$map$M_re <- factor(target_map, levels = levels(mapping))
  input$data$n_M_re[] <- 1L
  input$data$log_b_model <- 1L
  input
}

#' Restore legacy mortality parameterization for refit comparisons
#'
#' Match the original mortality baseline and age-specific offsets when
#' comparing fitted parameter estimates after migration.
#'
#' @param pars Parameters from the refitted model.
#' @param fit The original fitted legacy model.
#'
#' @return The refitted parameters expressed in their legacy parameterization.
#' @keywords internal
refit_legacy_mortality_parameters <- function(pars, fit) {
  old_M0 <- fit$parList$M0
  if (is.null(old_M0)) {
    old_M0 <- fit$input$par$M0
  }
  if (!is.null(old_M0)) {
    pars$M0 <- if (fit$input$data$M_model %in% c(1L, 3L)) {
      pars$Mpars[1L]
    } else {
      old_M0
    }
    pars$M_a <- pars$M_a - pars$M0
  }
  pars
}
