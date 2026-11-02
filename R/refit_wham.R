#' Refit a WHAM model from a fitted model object
#'
#' Uses the stored input and fitted parameter values, preserving parameter
#' mappings and manual input changes. Refits use the \code{\link{fit_wham}}
#' options recorded in \code{fit$call} (e.g. \code{do.sdrep}, \code{do.retro},
#' \code{do.osa}, \code{do.brps}). Parameter estimates and negative
#' log-likelihoods are compared between the original and refitted models.
#'
#' For each parameter, the difference between the refit and original estimate
#' is divided by the original standard error (from \code{fit$sdrep}) to give a
#' z-like value. Parameters with \eqn{|z| > 1.96} are flagged when standard
#' errors are available. This is a descriptive check, not a formal test,
#' because both fits use the same data.
#'
#' Older inputs are migrated to current structures while retaining fitted
#' values, parameter mappings, observations, and supported likelihoods.
#' Migration covers historical selectivity, recruitment, mortality,
#' environmental, and age-composition layouts when their conventions can be
#' identified. Unsupported or ambiguous layouts are rejected rather than
#' guessed. Native compatibility modes require rebuilding WHAM.
#'
#' Parameter differences are compared in the original parameterization, even
#' when migration changes parameter names or order. Recruitment decoupling is
#' recovered from stored options, data flags, or active correlation parameters.
#'
#' @param fit a fitted WHAM model object containing \code{$input} and
#'   \code{$opt}.
#'   Older fits without a usable \code{$call} are refit using WHAM defaults.
#'
#' @return The refitted WHAM model object, returned invisibly. A comparison
#'   table and flagged parameters are printed to the console.
#'
#' @seealso \code{\link{prepare_wham_input_from_fit}}, \code{\link{fit_wham}}
#'
#' @export
refit_wham <- function(fit = NULL) {
  if (is.null(fit) || is.null(fit$opt$par)) {
    stop("Provide a fitted WHAM model object.")
  }

  legacy <- is.null(fit$input$data$n_stocks) &&
    !is.null(fit$input$data$n_NAA_sigma)
  input <- if (legacy) migrate_legacy_refit_input(fit) else fit$input
  if (!legacy && !is.null(fit$parList)) {
    input$par[names(fit$parList)] <- fit$parList
  }
  if (
    is.null(input$data$selpars_re_index) && !is.null(input$data$selpars_est)
  ) {
    data <- input$data
    input <- migrate_refit_selectivity(input)
    if (isTRUE(data$do_proj == 1L)) {
      stop("Cannot safely migrate stored projections in this WHAM input.")
    }
    defaults <- set_proj(set_brp_static(input, list()))$data
    missing <- setdiff(names(defaults), names(input$data))
    input$data[missing] <- defaults[missing]
    if (!is.matrix(input$data$proj_Fcatch)) {
      input$data$proj_Fcatch <- matrix(input$data$proj_Fcatch)
    }
    if (is.null(input$data$bias_correct_brps)) {
      input$data$bias_correct_brps <- 0L
    }
    if (is.null(input$data$do_simulate_period)) {
      input$data$do_simulate_period <- c(1L, 1L)
    }
    if (!is.null(data$avg_years_ind)) {
      for (field in grep("^avg_years_ind_", missing, value = TRUE)) {
        target <- input$data[[field]]
        year_axis <- length(dim(target))
        if (year_axis == 2L) {
          target[] <- 0
          target[seq_len(length(data$avg_years_ind) + 1L), ] <- c(
            length(data$avg_years_ind),
            data$avg_years_ind
          )
        } else {
          dim(target) <- c(
            prod(dim(target)[-year_axis]),
            dim(target)[year_axis]
          )
          target[] <- 0
          target[, seq_len(length(data$avg_years_ind) + 1L)] <- rep(
            c(length(data$avg_years_ind), data$avg_years_ind),
            each = nrow(target)
          )
          dim(target) <- dim(input$data[[field]])
        }
        input$data[[field]] <- target
      }
    }
  }
  if (is.null(input$data$decouple_recruitment)) {
    input$data$decouple_recruitment <- as.integer(refit_decouple_recruitment(
      fit
    ))
  }

  # Zero, fixed correlations are invariant to historical transform conventions.
  # Preserve stock-level sharing when adding the regional dimension.
  naa_matrices <- c(log_NAA_sigma = input$data$n_ages, trans_NAA_rho = 2L)
  old_naa <- vapply(
    names(naa_matrices),
    function(parameter) {
      value <- input$par[[parameter]]
      !is.null(value) &&
        is.matrix(value) &&
        identical(
          as.integer(dim(value)),
          as.integer(c(input$data$n_stocks, naa_matrices[[parameter]]))
        )
    },
    logical(1)
  )
  if (any(old_naa)) {
    valid_maps <- TRUE
    inactive <- length(input$data$NAA_re_model) > 0L &&
      all(!is.na(input$data$NAA_re_model) & input$data$NAA_re_model == 0L)
    for (parameter in names(naa_matrices)) {
      value <- input$par[[parameter]]
      if (is.null(value)) {
        next
      }
      mapping <- input$map[[parameter]]
      valid_maps <- valid_maps &&
        !is.null(mapping) &&
        length(mapping) == length(value)
      inactive <- inactive && !is.null(mapping) && all(is.na(mapping))
    }
    independent <-
      length(input$data$NAA_re_model) == input$data$n_stocks &&
      all(input$data$NAA_re_model %in% 0:2) &&
      !is.null(input$par$trans_NAA_rho) &&
      all(is.finite(input$par$trans_NAA_rho) & input$par$trans_NAA_rho == 0) &&
      !is.null(input$map$trans_NAA_rho) &&
      all(is.na(input$map$trans_NAA_rho))
    if (!valid_maps || !(inactive || independent)) {
      stop(
        "Cannot safely migrate active or unsupported stored NAA matrices; legacy correlation transforms are not inferred."
      )
    }
    expand_regions <- function(value) {
      array(value, c(input$data$n_stocks, 1L, ncol(value)))[,
        rep(1L, input$data$n_regions),
        ,
        drop = FALSE
      ]
    }
    if (old_naa[["log_NAA_sigma"]]) {
      mapping <- input$map$log_NAA_sigma
      mapping_array <- expand_regions(matrix(
        as.character(mapping),
        input$data$n_stocks
      ))
      input$par$log_NAA_sigma <- expand_regions(input$par$log_NAA_sigma)
      if (input$data$n_regions > 1L) {
        input$map$log_NAA_sigma <- factor(
          as.vector(mapping_array),
          levels = levels(mapping)
        )
      }
    }
    if (old_naa[["trans_NAA_rho"]]) {
      rho <- expand_regions(input$par$trans_NAA_rho)
      extra <- input$data$n_stocks * input$data$n_regions
      input$par$trans_NAA_rho <- array(
        c(as.vector(rho), rep(0, extra)),
        c(input$data$n_stocks, input$data$n_regions, 3L)
      )
      mapping <- input$map$trans_NAA_rho
      input$map$trans_NAA_rho <- factor(
        rep(NA_character_, length(input$par$trans_NAA_rho)),
        levels = levels(mapping)
      )
    }
  }

  # Early no-movement inputs used a scalar instead of the region matrix.
  if (
    is.null(dim(input$data$mu_model)) &&
      identical(as.numeric(input$data$mu_model), 1) &&
      length(input$data$can_move) > 0L &&
      all(is.finite(input$data$can_move) & input$data$can_move == 0)
  ) {
    input$data$mu_model <- matrix(
      1L,
      input$data$n_regions,
      input$data$n_regions - 1L
    )
  }

  if (is.null(input$data$years_use)) {
    input$data$years_use <- seq_len(input$data$n_years_model) - 1L
  }

  scales <- list(input = input, selectivity = FALSE, ecov = integer())
  if (!legacy && !is.null(fit$input$data$selpars_est)) {
    scales <- refit_process_scales(input, fit$rep)
    input <- scales$input
  }

  if (legacy || !is.null(fit$input$data$selpars_est)) {
    input <- refit_initial_equilibrium(input, fit$rep$N1)
  }

  if (is.call(fit$call)) {
    fit_args <- as.list(fit$call)[-1]
    fit_args$input <- NULL
  } else {
    # Old saved fits may not contain the call recorded by fit_wham().
    # Continue with fit_wham() defaults and infer retrospective settings below.
    fit_args <- list()
  }

  # A recorded call can contain symbols for options supplied by the caller.
  # Resolve those in the refit's calling environment before falling back to
  # base constants; evaluating only in baseenv() breaks otherwise valid calls.
  caller_env <- parent.frame()
  fit_args <- lapply(fit_args, function(arg) {
    if (!is.language(arg)) {
      return(arg)
    }
    tryCatch(
      eval(arg, envir = caller_env),
      error = function(error) eval(arg, envir = baseenv())
    )
  })

  # Older fits (including <= 1.0.9) may not retain the original fit options.
  # The presence of peels is the only reliable indication that a retrospective
  # was completed, but do not override an explicit value in the saved call.
  if (is.null(fit_args$do.retro)) {
    fit_args$do.retro <- !is.null(fit$peels)
  }
  if (
    isTRUE(fit_args$do.retro) &&
      is.null(fit_args$n.peels) &&
      !is.null(fit$peels)
  ) {
    fit_args$n.peels <- length(fit$peels)
  }
  legacy_steepness <- input$legacy_steepness
  input$legacy_steepness <- NULL
  refit <- do.call(fit_wham, c(list(input = input), fit_args))

  est_fit <- fit$opt$par
  est_refit <- refit$opt$par
  if (scales$selectivity || length(scales$ecov)) {
    pars <- refit$parList
    if (scales$selectivity) {
      pars$sel_repars[, 2:3] <- pars$sel_repars[, 2:3] / 2
    }
    if (length(scales$ecov)) {
      columns <- scales$ecov
      rho <- tanh(pars$Ecov_process_pars[3, columns] / 2)
      pars$Ecov_process_pars[2, columns] <- pars$Ecov_process_pars[2, columns] -
        0.5 * log1p(-rho^2)
    }
    for (parameter in c("sel_repars", "Ecov_process_pars")) {
      values <- as.vector(pars[[parameter]])
      mapping <- input$map[[parameter]]
      if (!is.null(mapping)) {
        values <- values[match(levels(mapping), as.character(mapping))]
      }
      est_refit[names(est_refit) == parameter] <- values
    }
  }
  if (length(est_fit) != length(est_refit)) {
    stop(
      "The original and refitted models have different numbers of parameters."
    )
  }
  if (legacy) {
    est_refit <- est_fit
    pars <- refit$parList
    if (any(fit$input$data$Ecov_model != 0L)) {
      beta <- fit$input$par$Ecov_beta
      layout <- refit_legacy_ecov_layout(fit$input$data, beta)
      for (parameter in names(layout$indices)) {
        positions <- layout$indices[[parameter]]
        active <- positions > 0L
        beta[positions[active]] <- pars[[parameter]][active]
      }
      pars$Ecov_beta <- beta
    }
    pars$log_N1_pars <- as.vector(pars$log_N1)[seq_along(
      fit$input$par$log_N1_pars
    )]
    pars$log_F1 <- pars$F_pars[1, ]
    pars$F_devs <- pars$F_pars[-1, , drop = FALSE]
    pars$M_a <- as.vector(pars$Mpars)[seq_along(fit$input$par$M_a)]
    pars <- refit_legacy_mortality_parameters(pars, fit)
    pars <- refit_legacy_steepness_parameters(pars, legacy_steepness)
    for (parameter in c("catch_paa_pars", "index_paa_pars")) {
      if (!is.matrix(fit$input$par[[parameter]])) {
        pars[[parameter]] <- pars[[parameter]][
          refit_legacy_age_comp_indices(fit$input$data, parameter)
        ]
      }
    }
    pointers <- fit$input$data$NAA_sigma_pointers
    if (refit_decouple_recruitment(fit)) {
      rho <- tanh(as.vector(pars$trans_NAA_rho) / 2)
      if (fit$input$data$n_NAA_sigma > 1) {
        pars$log_NAA_sigma[1] <- pars$log_NAA_sigma[1] + 0.5 * log1p(-rho[1]^2)
        pars$log_NAA_sigma[-1] <- pars$log_NAA_sigma[-1] -
          0.5 * log1p(-rho[2]^2) +
          0.5 * log1p(-rho[3]^2)
      }
      pars$trans_NAA_rho[] <- as.vector(pars$trans_NAA_rho)[c(1L, 3L, 2L)]
    }
    pars$log_NAA_sigma <- as.vector(pars$log_NAA_sigma)[
      match(seq_along(fit$input$par$log_NAA_sigma), pointers)
    ]
    pars$trans_NAA_rho <- as.vector(pars$trans_NAA_rho)[
      seq_along(fit$input$par$trans_NAA_rho)
    ] / 2
    pars$q_repars[, 2] <- pars$q_repars[, 2] / 2
    pars$sel_repars[, 2:3] <- pars$sel_repars[, 2:3] / 2
    pars$M_repars <- as.vector(pars$M_repars)
    pars$M_repars[2:3] <- pars$M_repars[2:3] / 2
    for (parameter in unique(names(est_fit))) {
      values <- as.vector(pars[[parameter]])
      mapping <- fit$input$map[[parameter]]
      if (!is.null(mapping)) {
        values <- values[match(levels(mapping), as.character(mapping))]
      }
      positions <- which(names(est_fit) == parameter)
      if (length(values) < length(positions)) {
        stop("Could not compare migrated parameter: ", parameter)
      }
      est_refit[positions] <- values[seq_along(positions)]
    }
  }

  se <- rep(NA_real_, length(est_fit))
  cov_fixed <- fit$sdrep$cov.fixed
  if (!is.null(cov_fixed) && nrow(cov_fixed) == length(est_fit)) {
    se <- sqrt(diag(cov_fixed))
  }

  fit_comparison <- data.frame(
    nll_original = fit$opt$objective,
    nll_refit = refit$opt$objective,
    nll_diff = refit$opt$objective - fit$opt$objective,
    max_abs_par_diff = max(abs(est_refit - est_fit)),
    conv_original = fit$opt$convergence,
    conv_refit = refit$opt$convergence
  )

  param_check <- data.frame(
    parameter = names(est_fit),
    estimate_fit = est_fit,
    estimate_refit = est_refit,
    diff = est_refit - est_fit,
    se_fit = se,
    z = (est_refit - est_fit) / se
  )
  param_check$significant <- !is.na(param_check$z) & abs(param_check$z) > 1.96

  print(fit_comparison)
  n_sig <- sum(param_check$significant)
  message(
    n_sig,
    " of ",
    nrow(param_check),
    " parameters flagged with |z| > 1.96",
    if (all(is.na(se))) {
      " (no standard errors available, so no parameters could be checked)"
    } else {
      ""
    }
  )
  if (n_sig > 0) {
    print(param_check[param_check$significant, ], row.names = FALSE)
  }

  invisible(refit)
}
