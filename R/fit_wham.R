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
#' @param input Named list with several components:
#'   \describe{
#'     \item{\code{input$data}}{Data to fit the assessment model to.}
#'     \item{\code{input$par}}{Parameters, a list of all parameter objects required by the user template (both random and fixed effects). See \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{input$map}}{Map, a mechanism for collecting and fixing parameters. See \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{input$random}}{Character vector defining the parameters to treat as random effect. See \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{input$years}}{Numeric vector of the years which the model spans. Not important for model fitting, but useful for plotting.}
#'     \item{\code{input$map.n}}{I don't know what this is.}
#'     \item{\code{input$model_name}}{Character, name of the model, e.g. \code{"Yellowtail flounder"}}
#'     \item{\code{input$ages.lab}}{Character vector of the age labels, e.g. \code{c("1","2","3","4+").}}
#'   }
#' @param n.newton integer, number of additional Newton steps after optimization. Passed to \code{\link{fit_tmb}}. Default = \code{3}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? See \code{\link[TMB]{sdreport}}. Default = \code{TRUE}.
#' @param do.retro T/F, do retrospective analysis? Default = \code{TRUE}.
#' @param n.peels integer, number of peels to use in retrospective analysis. Default = \code{7}.
#' @param do.osa T/F, Calculate one-step-ahead (OSA) residuals? Default = \code{TRUE}. See details. Returned
#'   as \code{mod$osa$residual}.
#' @param osa.opts list of options for calculating OSA residuals, passed to \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}.
#'   Default: \code{osa.opts = list(method="oneStepGeneric", parallel=TRUE)}.
#' @param model (optional), a previously fit wham model.
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
#' @seealso \code{\link{fit_tmb}}, \code{\link{retro}}, \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}
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
fit_wham = function(input, n.newton = 3, do.sdrep = TRUE, do.retro = TRUE, n.peels = 7, do.osa = TRUE, osa.opts = list(method="oneStepGeneric", parallel=TRUE), model=NULL)
{
  # wham.dir <- find.package("wham")
  # dyn.load( paste0(wham.dir,"/libs/", TMB::dynlib(version)) )
  if(missing(model)){ 
    mod <- TMB::MakeADFun(input$data,input$par, DLL = "wham", random = input$random, map = input$map)
  } else {mod = model}
  
  mod <- fit_tmb(mod, n.newton = n.newton, do.sdrep = do.sdrep)
  mod$years <- input$years
  mod$ages.lab <- input$ages.lab
  mod$model_name <- input$model_name

  # only if no error
  if(is.null(mod$err)){
    if(do.retro) mod$peels = retro(mod, ran = unique(names(mod$env$par[mod$env$random])), n.peels= n.peels)
    if(do.osa){
      cat("Doing OSA residuals...\n");
      OSA <- TMB::oneStepPredict(obj=mod, observation.name="obsvec",
                                  data.term.indicator="keep",
                                  method=osa.opts$method,
                                  discrete=FALSE, parallel=osa.opts$parallel)
      input$data$obs$residual <- OSA$residual;
      mod$osa <- input$data$obs
    }
  } else warning(paste("","** Did not do sdrep, retro, or OSA residual analyses. **",
    "Error during Newton steps. Check for unidentifiable parameters.","",mod$err,"",sep='\n'))

  return(mod)
}
