#' Fit model peeling off \emph{i} years of data
#'
#' Fits a single retrospective peel by removing the specified number of years
#' of data and calling \code{\link{fit_tmb}}. Called by
#' \code{\link{retro}} and available for fitting an individual peel directly.
#'
#' @param peel Integer, number of years of data to remove before model fitting.
#' @param input input with same structure as that provided by \code{\link{prepare_wham_input}}. May want to use input$par = model$parList to start at MLEs.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? Default = \code{FALSE}.
#' @param n.newton integer, number of additional Newton steps after optimizafit_tmbtion for each peel. Default = \code{3}.
#' @param MakeADFun.silent T/F, Passed to silent argument of \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}. Default = \code{FALSE}.
#' @param retro.silent T/F, suppress the message identifying the peel number?
#'   Default = \code{FALSE}.
#' @param save.input T/F, should modified input list be saved? Necessary to project from a peel but increases model object size. Default = \code{FALSE}.
#'
#' @return \code{out}, output of \code{\link{fit_tmb}} for peel \emph{i}
#'
#' @export
#' 
#' @seealso \code{\link{fit_wham}}, \code{\link{retro}}, \code{\link{fit_tmb}}
#'
fit_peel <- function(peel, input, do.sdrep = FALSE, n.newton = 3, MakeADFun.silent = FALSE, retro.silent = FALSE, save.input = FALSE)
{
  out <- list()
  validate_peel_count(peel, input$years)
  if(!retro.silent) message(paste0("Retro Peel: ", peel))
  validate_input_safety(input)
  temp <- reduce_input(input, tail(input$years,peel))
  validate_input_safety(temp)
  temp.mod <- TMB::MakeADFun(temp$data, temp$par, DLL="wham", random = temp$random, map = temp$map, silent = MakeADFun.silent)

  out <- fit_tmb(temp.mod, do.sdrep = do.sdrep, n.newton = n.newton, do.check=FALSE)
  out$peel <- peel
  out$ages.lab <- input$ages.lab
  out$model_name <- paste0(input$model_name, " peel ",peel)
  out$years <- temp$years
  if(save.input){
    out$input <- temp
    out$input$model_name <- paste0(input$model_name, " peel ",peel)
  }

  return(out)
}
