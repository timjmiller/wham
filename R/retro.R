#' Run retrospective analysis
#'
#' Internal function called by \code{\link{fit_wham}}. Calls \code{\link{fit_peel}}
#' to fit the model peeling off \code{1, 2, ..., n.peels} years of data.
#'
#' @param model Optimized TMB model, output from \code{\link{fit_tmb}}.
#' @param n.peels Integer, number of peels to use in retrospective analysis. Default = \code{7}.
#' @param ran Character, specifies which parameters to treat as random effects. Default = \code{"log_NAA"}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters for each peel? Default = \code{FALSE}.
#' @param n.newton integer, number of additional Newton steps after optimization for each peel. Default = \code{0}.
#' @param MakeADFun.silent T/F, Passed to silent argument of \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}. Default = \code{FALSE}.
#' @param retro.silent T/F, Passed to argument of internal fit_peel function. Determines whether peel number is printed to screen. Default = \code{FALSE}.
#' @param save.input T/F, should modified input list be saved for every peel? Necessary to project from a peel but increases model object size. Default = \code{FALSE}.
#' 
#' @return \code{peels}, a list of length \code{n.peels}, where entry \emph{i} is a model
#' fit by peeling off \emph{i} years of data.
#'
#' @export
#' 
#' @seealso \code{\link{fit_wham}}, \code{\link{fit_peel}}
#'
retro = function(model, n.peels = 7, ran = "log_NAA", do.sdrep = FALSE, n.newton = 0, MakeADFun.silent = FALSE, retro.silent = FALSE, save.input = FALSE)
{
  temp = list(data = model$env$data, par = model$parList, map = model$env$map, random = ran, years=model$years, years_full=model$years_full, ages.lab=model$ages.lab, model_name=model$model_name)
  if(n.peels>0) peels = list(fit_peel(1, input = temp, do.sdrep = do.sdrep, n.newton = n.newton, MakeADFun.silent = MakeADFun.silent, retro.silent = retro.silent, save.input = save.input))
  if(n.peels>1) for(i in 2:n.peels) peels[[i]] = fit_peel(i, input = temp, do.sdrep = do.sdrep, n.newton = n.newton, MakeADFun.silent = MakeADFun.silent, retro.silent = retro.silent, save.input = save.input)
  return(peels)
}
