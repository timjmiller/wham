#' Add TMB sdreport object to WHAM model
#'
#' Runs \code{\link[TMB:sdreport]{TMB::sdreport}} and adds the object to the fitted (and projected) model list. E.g., fit$sdrep.
#'
#' @param model a fitted WHAM model object returned by fit_wham or project_wham.
#' @param save.sdrep  T/F, save the full \code{\link[TMB]{TMB::sdreport}} object? If \code{FALSE}, only save \code{\link[TMB:summary.sdreport]{summary.sdreport}} to reduce model object file size. Default = \code{TRUE}.
#' @param TMB.bias.correct T/F whether to use the bias.correct feature of TMB::sdreport. Default = \code{FALSE}.
#' @param TMB.jointPrecision T/F whether to return the joint precision matrix for the fixed and random effects from TMB::sdreport. Default = \code{FALSE}.
#'
#' @export
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{project_wham}}
#'
do_sdreport <- function(model, save.sdrep = TRUE, TMB.bias.correct = FALSE, TMB.jointPrecision = FALSE) {
  if(model$is_sdrep & !is.null(model$opt)){ #much faster than old way below
    model$rep <- model$report(model$env$last.par.best)
    model$sdrep <- try(TMB::sdreport(model, par.fixed = model$opt$par, hessian.fixed = solve(model$sdrep$cov.fixed), bias.correct = TMB.bias.correct, getJointPrecision = TMB.jointPrecision))
  } else {
    model$sdrep <- try(TMB::sdreport(model, bias.correct = TMB.bias.correct, getJointPrecision = TMB.jointPrecision))
  }
  # model$sdrep <- try(TMB::sdreport(model))
  model$is_sdrep <- !is.character(model$sdrep)
  model$na_sdrep <- ifelse(model$is_sdrep, any(is.na(summary(model$sdrep,"fixed")[,2])), NA)
  # if(model$is_sdrep) model$na_sdrep <- any(is.na(summary(model$sdrep,"fixed")[,2])) else model$na_sdrep = NA
  if(!save.sdrep) model$sdrep <- summary(model$sdrep) # only save summary to reduce model object size
  return(model)
}
