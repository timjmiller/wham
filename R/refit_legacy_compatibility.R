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
