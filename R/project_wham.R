#' Project a fit WHAM model
#'
#' Provides projections/forecasts for an existing (already fit) WHAM model.
#'
#' WHAM implements five options for handling fishing mortality in the projections.
#' Exactly one of these must be specified in \code{proj.opts}:
#'   \itemize{
#'     \item Use last year F (default). Set \code{proj.opts$use.last.F = TRUE}. WHAM will use F in the terminal model year for projections.
#'     \item Use average F. Set \code{proj.opts$use.avg.F = TRUE}. WHAM will use F averaged over \code{proj.opts$avg.yrs} for projections (as is done for M-, maturity-, and weight-at-age).
#'     \item Use F at X% SPR. Set \code{proj.opts$use.FXSPR = TRUE}. WHAM will calculate F at X% SPR.
#'     \item Specify F. Provide \code{proj.opts$proj.F}, an F vector with length = \code{n.yrs}.
#'     \item Specify catch. Provide \code{proj.opts$proj.catch}, a vector of aggregate catch with length = \code{n.yrs}. WHAM will calculate F to get specified catch.
#'   }
#'
#' \code{proj.opts$avg.yrs} controls which years will be averaged over in the projections.
#' The following quantities are averaged:
#'   \itemize{
#'     \item Maturity-at-age
#'     \item Weight-at-age
#'     \item Natural mortality-at-age
#'     \item Fishing mortality-at-age (if \code{proj.opts$use.avgF = TRUE})
#'   }
#'
#' WHAM implements four options for handling the environmental covariate(s) in the projections.
#' Exactly one of these must be specified in \code{proj.opts} if \code{Ecov} is in the model:
#'   \describe{
#'     \item{(Default) Continue Ecov process model (e.g. random walk, AR1)}{Set \code{$cont.Ecov = TRUE}. WHAM will estimate the Ecov process in projection years (i.e. continue the random walk / AR1 process).}
#'     \item{Use last year Ecov(s)}{Set \code{$use.last.Ecov = TRUE}. WHAM will use Ecov value from the terminal year (of population model) for projections.}
#'     \item{Use average Ecov(s)}{Provide \code{$avg.yrs.Ecov}, a vector specifying which years to average over the environmental covariate(s) for projections.}
#'     \item{Specify \code{Ecov}}{Provide \code{$proj.Ecov}, a matrix of user-specified environmental covariate(s) to use for projections. Dimensions must be # projection years (\code{proj.opts$n.yrs}) x # Ecovs (\code{ncols(Ecov$mean)}).}
#'   }
#' If the original model fit the Ecov in years beyond the population model, WHAM will use the already-fit
#' Ecov values for the projections. If the Ecov model extended at least \code{proj.opts$n.yrs} years
#' beyond the population model, then none of the above need be specified.
#'
#' @param model a previously fit wham model
#' @param proj.opts a named list with the following components:
#'   \itemize{
#'     \item \code{$n.yrs} (integer), number of years to project/forecast. Default = \code{3}.
#'     \item \code{$use.last.F} (T/F), use terminal year F for projections. Default = \code{TRUE}.
#'     \item \code{$use.FXSPR} (T/F), calculate F at X% SPR for projections.
#'     \item \code{$proj.F} (vector), user-specified fishing mortality for projections. Length must equal \code{n.yrs}.
#'     \item \code{$proj.catch} (vector), user-specified aggregate catch for projections. Length must equal \code{n.yrs}.
#'     \item \code{$avg.yrs} (vector), specify which years to average over for calculating reference points. Default = last 5 model years, \code{tail(model$years, 5)}.
#'     \item \code{$cont.Ecov} (T/F), continue Ecov process (e.g. random walk or AR1) for projections. Default = \code{TRUE}.
#'     \item \code{$use.last.Ecov} (T/F), use terminal year Ecov for projections.
#'     \item \code{$avg.Ecov.yrs} (vector), specify which years to average over the environmental covariate(s) for projections.
#'     \item \code{$proj.Ecov} (matrix), user-specified environmental covariate(s) for projections. \code{n.yrs} rows.
#'     \item \code{$cont.Mre} (T/F), continue M random effects (i.e. AR1_y or 2D AR1) for projections. Default = \code{TRUE}. If \code{FALSE}, M will be averaged over \code{$avg.yrs} (which defaults to last 5 model years).
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
project_wham = function(model, proj.opts=list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE, use.FXSPR=FALSE,
                                              proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
                                              cont.Ecov=TRUE, use.last.Ecov=FALSE, avg.Ecov.yrs=NULL, proj.Ecov=NULL, cont.Mre=NULL),
                        n.newton=3, do.sdrep=TRUE)
{
  # modify wham input (fix parameters at previously estimated values, pad with NAs)
  tryCatch(input2 <- prepare_projection(model, proj.opts)
    , error = function(e) {err <<- conditionMessage(e)})

  # refit model to estimate derived quantities in projection years
  if(!exists("err")) mod <- fit_wham(input2, n.newton=n.newton, do.sdrep=do.sdrep, do.retro=F, do.osa=F, do.check=F, do.proj=F)
  if(exists("err")){
    mod <- model # if error, still pass previous/full fit
    mod$err_proj <- err # store error message to print out in fit_wham
  }

  # pass along previously calculated retros, OSA residuals, error messages, and runtime
  if(!is.null(model$peels)) mod$peels <- model$peels # retrospective analysis
  if(!is.null(model$osa)) mod$osa <- model$osa # OSA residuals
  if(!is.null(model$err)) mod$err <- model$err # error messages
  if(!is.null(model$err_retro)) mod$err_retro <- model$err_retro # error messages
  mod$runtime <- model$runtime # runtime (otherwise would be just for projections)

  # print error message
  if(!is.null(mod$err_proj)) warning(paste("","** Error during projections. **",
    paste0("Check for issues with proj.opts, see ?project_wham."),"",mod$err_proj,"",sep='\n'))

  return(mod)
}
