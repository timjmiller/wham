#' Project a fit WHAM model
#'
#' Provides projections/forecasts for an existing (already fit) WHAM model.
#'
#' WHAM implements four options for handling fishing mortality in the projections: 
#'   \describe{
#'     \item{Use last year F (default)}{Set \code{use.lastF = TRUE}. WHAM will use F in the terminal model year for projections.}
#'     \item{Use F at X% SPR}{Set \code{use.FXSPR = TRUE}. WHAM will calculate F at X% SPR.}
#'     \item{Specify F}{Provide \code{proj.F}, an F vector with length = \code{n.yrs}.}
#'     \item{Specify catch}{Provide \code{proj.catch}, a vector of aggregate catch with length = \code{n.yrs}. WHAM will calculate F to get specified catch.}
#'   }
#'
#' @param model a previously fit wham model
#' @param n.yrs integer, number of years to project/forecast. Default = \code{3}.
#' @param use.lastF T/F, use terminal year F for projections. Default = \code{TRUE}.
#' @param use.FXSPR T/F, calculate F at X% SPR for projections.
#' @param proj.F vector, user-specified fishing mortality for projections. Length must equal \code{n.yrs}.
#' @param proj.catch vector, user-specified aggregate catch for projections. Length must equal \code{n.yrs}.
#' @param avg.yrs vector, specify which years to average over for calculating reference points. Default = last 5 model years, \code{tail(mod$years, 5)}.
#' @param n.newton integer, number of additional Newton steps after optimization. Passed to \code{\link{fit_tmb}}. Default = \code{3}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? See \code{\link[TMB]{sdreport}}. Default = \code{TRUE}.
#' @param do.osa T/F, Calculate one-step-ahead (OSA) residuals? Default = \code{TRUE}. See details. Returned
#'   as \code{mod$osa$residual}.
#' @param osa.opts list of options for calculating OSA residuals, passed to \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}.
#'   Default: \code{osa.opts = list(method="oneStepGeneric", parallel=TRUE)}.
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
#' mod <- fit_wham(input) # using default values
#' mod <- fit_wham(input, do.retro=F, do.osa=F) # faster settings for initial model fitting
#'
#' mod_proj <- project_wham(mod) # using default values: use.lastF = TRUE, n.yrs = 3, avg.yrs = last 5 years
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
project_wham = function(model, n.yrs=3, use.lastF=TRUE, use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, 
  avg.yrs=NULL, n.newton=3, do.sdrep=TRUE, do.osa=TRUE, osa.opts=list(method="oneStepGeneric", parallel=TRUE))
{
  if(missing(model)){ 
    mod <- TMB::MakeADFun(input$data,input$par, DLL = "wham", random = input$random, map = input$map)
  } else {mod = model}
  
  mod <- fit_tmb(mod, n.newton = n.newton, do.sdrep = do.sdrep, do.check = do.check)
  mod$years <- input$years
  mod$ages.lab <- input$ages.lab
  mod$model_name <- input$model_name

  if(do.retro) tryCatch(mod$peels <- retro(mod, ran = unique(names(mod$env$par[mod$env$random])), n.peels= n.peels)
    , error = function(e) {err <<- conditionMessage(e)})
  if(exists("err")) mod$err_retro <- err # store error message

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

  if(!is.null(mod$err)) warning(paste("","** Error during Newton steps. **",
    "Check for unidentifiable parameters.","",mod$err,"",sep='\n'))

  if(!is.null(mod$err_retro)) warning(paste("","** Error during retrospective analysis. **",
    paste0("Check for issues with last ",n.peels," model years."),"",mod$err_retro,"",sep='\n'))  

  return(mod)
}
