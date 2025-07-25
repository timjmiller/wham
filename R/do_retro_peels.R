#' Fit retrospective peels and add them to the fitted model object
#'
#' Function to add fitted retrospective peels to fitted model returned by \code{\link{fit_wham}}. Calls \code{\link{retro}}.
#' to fit the model peeling off \code{1, 2, ..., n.peels} years of data.
#' This is just a wrapper for retro that instead of returning just the list of peels, returns the fitted model with the peels. 
#'
#' @param model Optimized TMB model, output from \code{\link{fit_tmb}}.
#' @param n.peels Integer, number of peels to use in retrospective analysis. Default = \code{7}.
#' @param ran Character, specifies which parameters to treat as random effects. Default = \code{"model$input$random"}.
#' @param use.mle T/F, use MLEs from full model fit as initial values for each peel? If not, the initial values from full model input are used. Default = \code{TRUE}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters for each peel? Default = \code{FALSE}.
#' @param n.newton integer, number of additional Newton steps after optimization for each peel. Default = \code{0}.
#' @param MakeADFun.silent T/F, Passed to silent argument of \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}. Default = \code{FALSE}.
#' @param retro.silent T/F, Passed to argument of internal fit_peel function. Determines whether peel number is printed to screen. Default = \code{FALSE}.
#' @param save.input T/F, should modified input list be saved for every peel? Necessary to project from a peel but increases model object size. Default = \code{FALSE}.
#' @param do.brps T/F, calculate and report biological reference points
#' @param check.version T/F, whether to verify the wham package commit and version for the fitted model are the same as the currently used package.
#' @param save.sdrep T/F, save the full \code{\link[TMB]{TMB::sdreport}} object? If \code{FALSE}, only save \code{\link[TMB:summary.sdreport]{summary.sdreport}} to reduce model object file size. Default = \code{FALSE}.
#' 
#' @return \code{peels}, a list of length \code{n.peels}, where entry \emph{i} is a model
#' fit by peeling off \emph{i} years of data.
#'
#' @export
#' 
#' @seealso \code{\link{fit_wham}}, \code{\link{retro}}, \code{\link{fit_peel}}
#'
do_retro_peels <- function(model, n.peels = 7, ran = NULL, use.mle = TRUE, do.sdrep = FALSE, n.newton = 0, MakeADFun.silent = FALSE, retro.silent = FALSE, save.input = FALSE, do.brps = FALSE, check.version = TRUE, save.sdrep = FALSE) {
  model$peels <- retro(model, n.peels, ran, use.mle, do.sdrep, n.newton, MakeADFun.silent, retro.silent, save.input, do.brps, check.version, save.sdrep)
  return(model)
}
