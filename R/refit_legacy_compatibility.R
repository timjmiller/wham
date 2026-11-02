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

#' Preserve legacy steepness parameters in the native stock-recruit model
#'
#' Retain the original parameterization and calculate stock-recruit
#' coefficients from annual unfished SPR inside the native objective.
#'
#' @param input The partially migrated WHAM input.
#' @param fit The original fitted legacy model.
#' @param data The original legacy input data.
#' @param pars The fitted legacy parameter values.
#'
#' @return The migrated input with native legacy recruitment codes and
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
  model <- data$recruit_model[1L]
  if (!model %in% 3:4 || length(pars$mean_rec_pars) != 2L ||
      any(!is.finite(pars$mean_rec_pars))) {
    stop("Cannot safely migrate invalid legacy steepness parameters.", call. = FALSE)
  }
  input$data$recruit_model <- as.integer(model + 2L)
  input$legacy_steepness <- list(recruit_model = model, native = TRUE)
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
  if (is.null(metadata) || isTRUE(metadata$native)) {
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
