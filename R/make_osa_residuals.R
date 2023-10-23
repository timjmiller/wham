#' Calculate one-step-ahead residuals
#' 
#' Standard residuals are not appropriate for models with random effects. Instead, one-step-ahead (OSA) residuals
#' can be used for evaluating model goodness-of-fit (\href{https://doi.org/10.1007/s10651-017-0372-4}{Thygeson et al. (2017)},
#' implemented in \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}). OSA residual options
#' are passed to \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}} in a list \code{osa.opts}. Current options are method: 
#' oneStepGaussianOffMode (default), oneStepGaussian, or oneStepGeneric, and parallel: TRUE/FALSE. 
#' See \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}} for further details.
#' It is not recommended to run this function (or \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}) with any random effects and
#' mvtweedie age composition likelihoods due to extensive computational demand. An error will be thrown in such cases. 
#' See \href{https://doi.org/10.1016/j.fishres.2022.106487}{Trijoulet et al. (2023)} for OSA methods for age composition OSA residuals.
#'
#' @param model A fit WHAM model, output from \code{\link{fit_wham}}.
#'
#' @return the same fit TMB model with additional elements for osa residuals:
#'   \describe{
#'     \item{\code{$OSA.Ecov}}{data.frame returned by \code{\link[TMB]{TMB::oneStepPredict}} for environmental observations, if applicable.}
#'     \item{\code{$OSA.agregate}}{data.frame returned by \code{\link[TMB]{TMB::oneStepPredict}} for aggregate catch and index 
#'        observations conditional on any environmental observations, if applicable.}
#'     \item{\code{$OSA.agecomp}}{data.frame returned by \code{\link[TMB]{TMB::oneStepPredict}} for age composition observations 
#'        conditional on any aggregate catch or index, or environmental observations, if applicable.}
#'     \item{\code{$osa}}{One-step-ahead residuals (if \code{do.osa=TRUE})}
#'   }
#'
#' @export
#'
#' @seealso \code{\link{fit_wham}}
#'
#' @examples
#' \dontrun{
#' data("input4_SNEMAYT") # load SNEMA yellowtail flounder data and parameter settings
#' mod <- fit_wham(input4_SNEMAYT, do.osa =FALSE, do.retro =FALSE)
#' mod <- make_osa_residuals(mod) # calculate Mohn's rho
#' plot_wham_output(mod)
#' }
#' 
make_osa_residuals = function(model,osa.opts = list(method="oneStepGaussianOffMode", parallel=TRUE), sdrep_required = TRUE){
  verify_version(model)
  orig_vals <- c(model$env$data$do_SPR_BRPs,model$env$data$do_MSY_BRPs)
  model$env$data$do_SPR_BRPs <- model$env$data$do_MSY_BRPs <- 0

  # one-step-ahead residuals
  if(is.null(osa.opts$method)) osa.opts$method <- "oneStepGaussianOffMode"
  if(!osa.opts$method %in% c("oneStepGaussianOffMode","oneStepGaussian", "oneStepGeneric")){
    stop(paste0("Only osa methods allowed for TMB::oneStepPredict currently in WHAM are oneStepGaussianoffMode (default), oneStepGaussian, or oneStepGeneric"))
  }
  if(is.null(osa.opts$parallel)) osa.opts$parallel <- TRUE
  if(sdrep_required & !model$is_sdrep) stop(paste0("Only allowing OSA residuals for models with TMB::sdreport completed"))
  if(any(model$input$data$age_comp_model_fleets == 10)) stop("OSA residuals do not seem possible with mvtweedie age composition likelihoods.")
  if(any(model$input$data$age_comp_model_indices == 10)) stop("OSA residuals do not seem possible with mvtweedie age composition likelihoods.")
  
  cat("Doing OSA residuals...\n");
  input = model$input
  if(any(class(input$data$obs) == "list")) input$data$obs <- as.data.frame(input$data$obs) #simulated data$obs will be a list
  model$osa = input$data$obs
  model$osa$residual = NA
  #first do continuous obs, condition on obs without osa (probably none)
  subset.agecomp = which(model$osa$type %in% c("indexpaa", "catchpaa"))
  subset.ecov = which(model$osa$type %in% c("Ecov"))
  subset.aggregate = which(model$osa$type %in% c("logindex", "logcatch"))
  conditional. = NULL
  if(length(subset.ecov)){ #do Ecov first
    cat("Doing OSA for Ecov observations...\n")
    model$OSA.Ecov = suppressWarnings(TMB::oneStepPredict(
      obj = model,
      method = osa.opts$method,
      subset = subset.ecov,
      conditional = conditional.,
      observation.name = "obsvec",
      data.term.indicator = "keep",
      discrete = FALSE,
      parallel = osa.opts$parallel))
    model$osa$residual[subset.ecov] <- model$OSA.Ecov$residual
    conditional. = c(conditional., subset.ecov)
  }
  if(length(subset.aggregate)){
    cat("Doing OSA for aggregate catch and index observations...\n")
    model$OSA.aggregate = suppressWarnings(TMB::oneStepPredict(
      obj = model,
      method = osa.opts$method,
      subset = subset.aggregate,
      conditional = conditional.,
      observation.name = "obsvec",
      data.term.indicator = "keep",
      discrete = FALSE,
      parallel = osa.opts$parallel))
    model$osa$residual[subset.aggregate] <- model$OSA.aggregate$residual
    conditional. = c(conditional., subset.aggregate)
  }
  if(length(subset.agecomp)){
    cat("Doing OSA for catch and index age comp observations...\n")
    if(!is.null(input$data$condition_no_osa)) cat("OSA not available for logistic-normal-01-infl, logistic-normal-01-infl-2par, or mvtweedie age comp likelihoods...\n")
    conditional. = c(conditional., input$data$condition_no_osa)
    cat("Doing OSA for age comp observations...\n")
    model$OSA.agecomp = suppressWarnings(TMB::oneStepPredict(
      obj = model,
      method = osa.opts$method,
      subset = subset.agecomp,
      conditional = conditional.,
      observation.name = "obsvec",
      data.term.indicator = "keep",
      discrete = FALSE,
      parallel = osa.opts$parallel))
    #remove any 0 residuals associated with last ages of multinomial and D-M
    model$osa$residual[subset.agecomp] <- model$OSA.agecomp$residual
    conditional. = c(conditional., subset.agecomp)
    ind <- which(input$data$age_comp_model_fleets %in% c(1:2,11))
    if(length(ind)) {
      for(i in ind) {
        NAind <- which(model$osa$age == max(model$osa$age, na.rm =T) & model$osa$fleet == paste0("fleet_", i))
        model$osa$residual[NAind] <- NA
      }
    }
    ind <- which(input$data$age_comp_model_indices %in% c(1:2,11))
    if(length(ind)) {
      for(i in ind) {
        NAind <- which(model$osa$age == max(model$osa$age, na.rm =T) & model$osa$fleet == paste0("index_", i))
        model$osa$residual[NAind] <- NA
      }
    }
  }
  model$env$data$do_SPR_BRPs <- orig_vals[1]
  model$env$data$do_MSY_BRPs <- orig_vals[2]
  return(model)
}
