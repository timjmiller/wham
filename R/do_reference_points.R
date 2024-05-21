#' Add reporting of biological reference points to WHAM model
#'
#' Changes internal flags to do the extra calculations and reporting for reference points.
#'
#' @param model a fitted WHAM model object returned by fit_wham or project_wham.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? See \code{\link[TMB]{sdreport}}. Default = \code{FALSE}.
#' @param save.sdrep  T/F, save the full \code{\link[TMB]{TMB::sdreport}} object? If \code{FALSE}, only save \code{\link[TMB:summary.sdreport]{summary.sdreport}} to reduce model object file size. Default = \code{TRUE}.
#'
#' @export
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{project_wham}}
#'

do_reference_points <- function(model, do.sdrep = FALSE, save.sdrep = TRUE){
  model$input$data$do_SPR_BRPs <- model$env$data$do_SPR_BRPs <- 1
  if(any(model$input$data$recruit_model %in% 3:4)) model$input$data$do_MSY_BRPs <- model$env$data$do_MSY_BRPs <- 1
  model <- check_which_F_age(model)
  # model$rep = model$report() #par values don't matter because function has not been evaluated
  model <- check_FXSPR(model)
  if(any(model$input$data$can_move==1) & any(model$input$data$mig_type == 1)){
    warning("Cannot currently calculate standard errors of biological reference points internally when survival and movement are simultaneous for any stock.")
  } else {
    if(do.sdrep) model <- do_sdreport(model, save.sdrep)
  }
  return(model)
}
