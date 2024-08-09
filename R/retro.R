#' Run retrospective analysis
#'
#' Internal function called by \code{\link{fit_wham}}. Calls \code{\link{fit_peel}}
#' to fit the model peeling off \code{1, 2, ..., n.peels} years of data.
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
#' 
#' @return \code{peels}, a list of length \code{n.peels}, where entry \emph{i} is a model
#' fit by peeling off \emph{i} years of data.
#'
#' @export
#' 
#' @seealso \code{\link{fit_wham}}, \code{\link{fit_peel}}
#'
#retro = function(model, n.peels = 7, ran = "log_NAA", do.sdrep = FALSE, n.newton = 0, MakeADFun.silent = FALSE, retro.silent = FALSE, save.input = FALSE)
retro = function(model, n.peels = 7, ran = NULL, use.mle = TRUE, do.sdrep = FALSE, n.newton = 0, MakeADFun.silent = FALSE, retro.silent = FALSE, save.input = FALSE, do.brps = FALSE, check.version = TRUE)
{
  data <- model$input$data
  par <- model$parList
  map <- model$input$map
  if(is.null(ran)) ran <- model$input$random

  if(check.version) verify_version(model)
  if(do.brps){
      if(any(data$can_move==1) & any(data$mig_type == 1)){
        warning("Cannot currently calculate standard errors of biological reference points internally when migration and movement are simultaneous for any stock.")
      } else {
        data$do_SPR_BRPs <- 1
        if(any(input$data$recruit_model %in% 3:4)) data$do_MSY_BRPs <- 1
      }
  } else {
    data$do_SPR_BRPs <- data$do_MSY_BRPs <- 0
  }
  temp <- model$input
  temp$random <- ran
  temp$data <- data
  if(use.mle) temp$par <- par
  peels <- list()
  if(n.peels>0) for(i in 1:n.peels) tryCatch(peels[[i]] <- 
    fit_peel(i, input = temp, do.sdrep = do.sdrep, n.newton = n.newton, MakeADFun.silent = MakeADFun.silent, retro.silent = retro.silent, 
      save.input = save.input), error = function(e) {peels[[i]]$err <<- conditionMessage(e)})
  # if(n.peels>0) for(i in 1:n.peels) peels[[i]] <- fit_peel(i, input = temp, do.sdrep = do.sdrep, n.newton = n.newton, 
  #   MakeADFun.silent = MakeADFun.silent, retro.silent = retro.silent, save.input = save.input)
  #if(n.peels>0) peels = list(fit_peel(1, input = temp, do.sdrep = do.sdrep, n.newton = n.newton, MakeADFun.silent = MakeADFun.silent, retro.silent = retro.silent, save.input = save.input))
  #if(n.peels>1) for(i in 2:n.peels) peels[[i]] = fit_peel(i, input = temp, do.sdrep = do.sdrep, n.newton = n.newton, MakeADFun.silent = MakeADFun.silent, retro.silent = retro.silent, save.input = save.input)
  return(peels)
}
