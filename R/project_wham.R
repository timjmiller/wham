#' Project a fit WHAM model
#'
#' Provides projections/forecasts for an existing (already fit) WHAM model.
#'
#' WHAM implements four options for handling fishing mortality in the projections.
#' Exactly one of these must be specified in \code{proj.opts}:
#'   \describe{
#'     \item{Use last year F (default)}{Set \code{proj.opts$use.lastF = TRUE}. WHAM will use F in the terminal model year for projections.}
#'     \item{Use F at X% SPR}{Set \code{proj.opts$use.FXSPR = TRUE}. WHAM will calculate F at X% SPR.}
#'     \item{Specify F}{Provide \code{proj.opts$proj.F}, an F vector with length = \code{n.yrs}.}
#'     \item{Specify catch}{Provide \code{proj.opts$proj.catch}, a vector of aggregate catch with length = \code{n.yrs}. WHAM will calculate F to get specified catch.}
#'   }
#'
#' @param model a previously fit wham model
#' @param proj.opts a named list with the following components:
#'   \describe{
#'     \item{\code{$n.yrs}}{integer, number of years to project/forecast. Default = \code{3}.}
#'     \item{\code{$use.lastF}}{T/F, use terminal year F for projections. Default = \code{TRUE}.}
#'     \item{\code{$use.FXSPR}}{T/F, calculate F at X% SPR for projections.}
#'     \item{\code{$proj.F}}{vector, user-specified fishing mortality for projections. Length must equal \code{n.yrs}.}
#'     \item{\code{$proj.catch}}{vector, user-specified aggregate catch for projections. Length must equal \code{n.yrs}.}
#'     \item{\code{$avg.yrs}}{vector, specify which years to average over for calculating reference points. Default = last 5 model years, \code{tail(model$years, 5)}.}
#'   }
#' @param n.newton integer, number of additional Newton steps after optimization. Passed to \code{\link{fit_tmb}}. Default = \code{0} for projections.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? See \code{\link[TMB]{sdreport}}. Default = \code{TRUE}.
#'
#' @return a projected WHAM model with additional output if specified:
#'   \describe{
#'     \item{\code{$rep}}{List of derived quantity estimates (see examples)}
#'     \item{\code{$sdrep}}{Parameter estimates (and standard errors if \code{do.sdrep=TRUE})}
#'     \item{\code{$peels}}{Retrospective analysis (if \code{do.retro=TRUE})}
#'     \item{\code{$osa}}{One-step-ahead residuals (if \code{do.osa=TRUE})}
#'   }
#'
#' @useDynLib wham
#' @export
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{fit_tmb}}
#'
#' @examples
#' \dontrun{
#' data("SNEMA_ytl") # load SNEMA yellowtail flounder data and parameter settings
#' mod <- fit_wham(input) # using default values (do.proj=T)
#'
#' mod2 <- fit_wham(input, do.retro=F, do.osa=F, do.proj=F) # fit model without projections, retro analysis, or OSA residuals
#' mod_proj <- project_wham(mod2) # add projections to previously fit model, using default values: use.lastF = TRUE, n.yrs = 3, avg.yrs = last 5 years
#'
#' names(mod_proj$rep) # list of derived quantities
#' tail(mod_proj$rep$SSB, 3) # get 3-year projected SSB estimates (weight, not numbers)

#' x = summary(mod_proj$sdrep)
#' unique(rownames(x))) # list of estimated parameters and derived quanitites with SE
#' x = x[rownames(x) == "log_SSB",] # SSB estimates with SE
#' ssb.mat = exp(cbind(x, x[,1] + qnorm(0.975)*cbind(-x[,2],x[,2])))/1000 # calculate 95% CI
#' colnames(ssb.mat) <- c("SSB","SSB_se","SSB_lower","SSB_upper")
#' tail(ssb.mat, 3) # 3-year projected SSB estimates with SE and 95% CI
#' }
project_wham = function(model, proj.opts=list(n.yrs=3, use.lastF=TRUE, use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL),
                        n.newton=0, do.sdrep=TRUE)
{
  # default: use average M, selectivity, etc. over last 5 model years to calculate ref points
  if(is.null(proj.opts$avg.yrs)) proj.opts$avg.yrs <- tail(model$years, 5)

  # check user options are valid 
  if(any(proj.opts$avg.yrs %in% model$years == FALSE)) stop(paste("","** Error setting up projections: **",
    "proj.opts$avg.yrs is not a subset of model years.","",sep='\n'))
  F.opt.ct <- sum(proj.opts$use.lastF, proj.opts$use.FXSPR, !is.null(proj.opts$proj.F), !is.null(proj.opts$proj.catch))
  if(F.opt.ct != 1) stop(paste("","** Error setting up projections: **",
    "Exactly one method of specifying F must be used (see ?project_wham).",
    "You have specified these in proj.opts:",
    capture.output(cat("  use.lastF = ",proj.opts$use.lastF)),
    capture.output(cat("  use.FXSPR = ",proj.opts$use.FXSPR)),
    capture.output(cat("  proj.F = ",proj.opts$proj.F)),
    capture.output(cat("  proj.catch = ",proj.opts$proj.catch)),"",sep='\n'))

  # fix parameters at previously estimated values
  # pad with NAs
  input2 <- prepare_projection(model)

  # 'refit' model, just estimates derived quantities in projection years (with uncertainty)
  mod <- fit_wham(input2, n.newton=n.newton, do.sdrep=do.sdrep, do.retro=F, do.osa=F, do.check=F)

  # pass along previously calculated retros, OSA residuals, and error messages
  if(!is.null(model$peels)) mod$peels <- model$peels # retrospective analysis
  if(!is.null(model$osa)) mod$osa <- model$osa # OSA residuals
  if(!is.null(model$err)) mod$err <- model$err # error messages 
  if(!is.null(model$err_retro)) mod$err_retro <- model$err_retro # error messages

  return(mod)
}
