#' Add TMB sdreport object to WHAM model
#'
#' Runs \code{\link[TMB:sdreport]{TMB::sdreport}} and adds the object to the fitted (and projected) model list. E.g., fit$sdrep.
#'
#' @param model a fitted WHAM model object returned by fit_wham or project_wham.
#' @param save.sdrep  T/F, save the full \code{\link[TMB]{TMB::sdreport}} object? If \code{FALSE}, only save \code{\link[TMB:summary.sdreport]{summary.sdreport}} to reduce model object file size. Default = \code{TRUE}.
#'
#' @export
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{project_wham}}
#'
do_sdreport <- function(model, save.sdrep = TRUE) {
  model$sdrep <- try(TMB::sdreport(model))
  model$is_sdrep <- !is.character(model$sdrep)
  if(model$is_sdrep) model$na_sdrep <- any(is.na(summary(model$sdrep,"fixed")[,2])) else model$na_sdrep = NA
  if(!save.sdrep) model$sdrep <- summary(model$sdrep) # only save summary to reduce model object size
  return(model)
}
