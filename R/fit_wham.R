#' Fit WHAM model
#'
#' Fits the compiled WHAM model using \code{\link[TMB:MakeADFun]{TMB::MakeADFun}} and
#' \code{\link[stats:nlminb]{stats::nlminb}}. Runs retrospective analysis if specified.
#'
#' Standard residuals are not appropriate for models with random effects. Instead, one-step-ahead (OSA) residuals
#' can be used for evaluating model goodness-of-fit (\href{https://link.springer.com/article/10.1007/s10651-017-0372-4}{Thygeson et al. (2017)},
#' implemented in \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}). Additional OSA residual options
#' are passed to \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}} in a list \code{osa.opts}. For example,
#' to use the (much faster, ~1 sec instead of 2 min) full Gaussian approximation instead of the (default)
#' generic method, you can use \code{osa.opts=list(method="fullGaussian")}.
#'
#' @param input Named list with components:
#'   \describe{
#'     \item{\code{$data}}{Data to fit the assessment model to.}
#'     \item{\code{$par}}{Parameters, a list of all parameter objects required by the user template (both random and fixed effects). See \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{$map}}{Map, a mechanism for collecting and fixing parameters. See \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{$random}}{Character vector defining the parameters to treat as random effect. See \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{$years}}{Numeric vector of the years which the model spans. Not important for model fitting, but useful for plotting.}
#'     \item{\code{$model_name}}{Character, name of the model, e.g. \code{"Yellowtail flounder"}}
#'     \item{\code{$ages.lab}}{Character vector of the age labels, e.g. \code{c("1","2","3","4+").}}
#'   }
#' @param n.newton integer, number of additional Newton steps after optimization. Passed to \code{\link{fit_tmb}}. Default = \code{3}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? See \code{\link[TMB]{sdreport}}. Default = \code{TRUE}.
#' @param do.retro T/F, do retrospective analysis? Default = \code{TRUE}.
#' @param n.peels integer, number of peels to use in retrospective analysis. Default = \code{7}.
#' @param do.osa T/F, calculate one-step-ahead (OSA) residuals? Default = \code{TRUE}. See details. Returned
#'   as \code{mod$osa$residual}.
#' @param osa.opts list of options for calculating OSA residuals, passed to \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}.
#'   Default: \code{osa.opts = list(method="oneStepGeneric", parallel=TRUE)}.
#' @param model (optional), a previously fit wham model.
#' @param do.check T/F, check if model parameters are identifiable? Passed to \code{\link{fit_tmb}}. Runs \code{\link[TMBhelper::Check_Identifiable]{TMBhelper::Check_Identifiable}}. Default = \code{TRUE}.
#' @param do.proj T/F, do projections? Default = \code{FALSE}. If true, runs \code{\link{project_wham}}.
#' @param proj.opts list of options for projections
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
#'     \item \code{$proj.Ecov} (vector), user-specified environmental covariate(s) for projections. Length must equal \code{n.yrs}.
#'   }
#'
#' @return a fit TMB model with additional output if specified:
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
#' @seealso \code{\link{fit_tmb}}, \code{\link{retro}}, \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}, \code{\link{project_wham}}
#'
#' @examples
#' \dontrun{
#' data("SNEMA_ytl") # load SNEMA yellowtail flounder data and parameter settings
#' mod = fit_wham(input) # using default values
#' mod = fit_wham(input, do.retro=FALSE, osa.opts=list(method="fullGaussian")) # faster settings for initial model fitting
#'
#' names(mod$rep) # list of derived quantities
#' mod$rep$SSB # get SSB estimates (weight, not numbers)
#' m1$rep$NAA[,1] # get recruitment estimates (numbers, first column of numbers-at-age matrix)
#' m1$rep$F[,1] # get F estimates for fleet 1
#' }
fit_wham = function(input, n.newton = 3, do.sdrep = TRUE, do.retro = TRUE, n.peels = 7,
                    do.osa = TRUE, osa.opts = list(method="oneStepGeneric", parallel=TRUE), model=NULL, do.check = FALSE,
                    do.proj = FALSE, proj.opts=list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE, use.FXSPR=FALSE,
                                              proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
                                              cont.Ecov=TRUE, use.last.Ecov=FALSE, avg.Ecov.yrs=NULL, proj.Ecov=NULL))
{
  # fit model
  if(missing(model)){
    mod <- TMB::MakeADFun(input$data, input$par, DLL = "wham", random = input$random, map = input$map)
  } else {mod = model}

  mod <- fit_tmb(mod, n.newton = n.newton, do.sdrep = do.sdrep, do.check = do.check)
  mod$years <- input$years
  mod$years_full <- input$years_full
  mod$ages.lab <- input$ages.lab
  mod$model_name <- input$model_name

  # retrospective analysis
  if(do.retro) tryCatch(mod$peels <- retro(mod, ran = unique(names(mod$env$par[mod$env$random])), n.peels= n.peels)
    , error = function(e) {err <<- conditionMessage(e)})
  if(exists("err")) mod$err_retro <- err # store error message

  # one-step-ahead residuals
  if(do.osa){
    if(mod$is_sdrep){ # only do OSA residuals if sdrep ran
      cat("Doing OSA residuals...\n");
      OSA <- suppressWarnings(TMB::oneStepPredict(obj=mod, observation.name="obsvec",
                                  data.term.indicator="keep",
                                  method=osa.opts$method,
                                  discrete=FALSE, parallel=osa.opts$parallel))
      input$data$obs$residual <- OSA$residual;
      mod$osa <- input$data$obs
    } else warning(paste("","** Did not do OSA residual analyses. **",
    "Error during TMB::sdreport(). Check for unidentifiable parameters.","",sep='\n'))
  }
  mod$input <- input

  # projections
  if(do.proj) mod <- project_wham(mod, proj.opts=proj.opts) # calls prepare_projection + fit_wham(do.proj=F)

  # error message reporting
  if(!is.null(mod$err)) warning(paste("","** Error during Newton steps. **",
    "Check for unidentifiable parameters.","",mod$err,"",sep='\n'))
  if(!is.null(mod$err_retro)) warning(paste("","** Error during retrospective analysis. **",
    paste0("Check for issues with last ",n.peels," model years."),"",mod$err_retro,"",sep='\n'))

  return(mod)
}
