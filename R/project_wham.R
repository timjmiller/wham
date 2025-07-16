#' Project a fit WHAM model
#'
#' Provides projections/forecasts for an existing (already fit) WHAM model.
#'
#' WHAM implements five options for handling fishing mortality in the projections.
#' Exactly one of these must be specified in \code{proj.opts}:
#'   \itemize{
#'     \item Use last year F (default). Set \code{proj.opts$use.last.F = TRUE}. WHAM will use F in the terminal model year for projections.
#'     \item Use average F. Set \code{proj.opts$use.avg.F = TRUE}. WHAM will use F averaged over \code{proj.opts$avg.yrs} for projections (as is done for M-, maturity-, and weight-at-age).
#'     \item Use F at X\% SPR. Set \code{proj.opts$use.FXSPR = TRUE}. WHAM will calculate F at X\% SPR.
#'     \item Specify F. Provide \code{proj.opts$proj.F}, an F vector with length = \code{n.yrs}.
#'     \item Specify catch. Provide \code{proj.opts$proj.catch}, a vector of aggregate catch with length = \code{n.yrs}. WHAM will calculate F to get specified catch.
#'   }
#'
#' \code{proj.opts$avg.yrs} controls which years the following will be averaged over in the projections:
#'   \itemize{
#'     \item Maturity-at-age
#'     \item Weight-at-age
#'     \item Natural mortality-at-age
#'     \item Fishing mortality-at-age (if \code{proj.opts$use.avgF = TRUE})
#'   }
#'
#' If fitting a model with recruitment estimated freely in each year, i.e. as fixed effects as in ASAP, WHAM handles recruitment 
#' in the projection years similarly to using the empirical cumulative distribution function. WHAM does this by calculating the mean
#' and standard deviation of log(R) over all model years (default) or a specified subset of years (\code{proj.opts$avg.rec.yrs}). WHAM then
#' treats recruitment in the projections as a random effect with this mean and SD, i.e. log(R) ~ N(meanlogR, sdlogR).
#' 
#' WHAM implements four options for handling the environmental covariate(s) in the projections.
#' Exactly one of these must be specified in \code{proj.opts} if \code{ecov} is in the model:
#'   \describe{
#'     \item{(Default) Continue ecov process model (e.g. random walk, AR1)}{Set \code{$cont.ecov = TRUE}. WHAM will estimate the ecov process in projection years (i.e. continue the random walk / AR1 process).}
#'     \item{Use last year ecov(s)}{Set \code{$use.last.ecov = TRUE}. WHAM will use ecov value from the terminal year (of population model) for projections.}
#'     \item{Use average ecov(s)}{Provide \code{$avg.yrs.ecov}, a vector specifying which years to average over the environmental covariate(s) for projections.}
#'     \item{Specify \code{ecov}}{Provide \code{$proj.ecov}, a matrix of user-specified environmental covariate(s) to use for projections. Dimensions must be # projection years (\code{proj.opts$n.yrs}) x # ecovs (\code{ncols(ecov$mean)}).}
#'   }
#' If the original model fit the ecov in years beyond the population model, WHAM will use the already-fit
#' ecov values for the projections. If the ecov model extended at least \code{proj.opts$n.yrs} years
#' beyond the population model, then none of the above need be specified.
#'
#' @param model a previously fit wham model
#' @param proj.opts a named list with the following components:
#'   \itemize{
#'     \item \code{$n.yrs} (integer), number of years to project/forecast. Default = \code{3}.
#'     \item \code{$use.last.F} (T/F), use terminal year F for projections. Default = \code{TRUE}.
#'     \item \code{$use.avg.F} (T/F), use average of F over certain years for projections. Default = \code{FALSE}. Years to average over determined by $avg.yrs defined below.
#'     \item \code{$use.FXSPR} (T/F), calculate and use F at X\% SPR for projections. Default = \code{FALSE}.
#'     \item \code{$use.FMSY} (T/F), calculate and use FMSY for projections. Default = \code{FALSE}.
#'     \item \code{$proj.F} (vector), user-specified fishing mortality for projections. Length must equal \code{n.yrs}.
#'     \item \code{$proj.catch} (vector), user-specified aggregate catch for projections. Length must equal \code{n.yrs}.
#'     \item \code{$avg.yrs} (vector), specify which years to use to average population attributes (MAA,FAA,WAA,maturity,movement) in projection years. Any BRPs calculated in projection years will also use these. Default = last 5 years, \code{tail(model$years, 5)}.
#'     \item \code{$cont.ecov} (T/F), continue ecov process (e.g. random walk or AR1) for projections. Default = \code{TRUE}.
#'     \item \code{$use.last.ecov} (T/F), use terminal year ecov for projections.
#'     \item \code{$avg.ecov.yrs} (vector), specify which years to average the environmental covariate(s) over for projections.
#'     \item \code{$proj.ecov} (matrix), user-specified environmental covariate(s) for projections. \code{n.yrs x n.ecov}.
#'     \item \code{$cont.M.re} (T/F), continue M random effects (i.e. AR1_y or 2D AR1) for projections. Default = \code{FALSE}. If \code{FALSE}, M will be averaged over \code{$avg.yrs.M} (which defaults to last 5 model years).
#'     \item \code{$cont.move.re} (T/F), continue any movement random effects for projections. Default = \code{FALSE}. If \code{FALSE}, movement parameters will be averaged over \code{$avg.yrs.move} (which defaults to last 5 model years).
#'     \item \code{$cont.L.re} (T/F), continue any L ("extra mortality rate") random effects for projections. Default = \code{FALSE}. If \code{FALSE}, L parameters will be averaged over \code{$avg.yrs.L} (which defaults to last 5 model years).
#'     \item \code{$avg.rec.yrs} (vector), specify which years to calculate the CDF of recruitment for use in projections. Default = all model years. Only used when recruitment is estimated as fixed effects (SCAA).
#'     \item \code{$percentFXSPR} (scalar), percent of F_XSPR to use for projections, only used if $use.FXSPR = TRUE. For example, to project with F = 75\% F_40\%SPR, \code{proj.opts$percentFXSPR = 75}. Default = 100.
#'     \item \code{$percentFMSY} (scalar), percent of F_MSY to use for projections, only used if $use.FMSY = TRUE and a stock-recruit relationship is assumed. Default = 100.
#'     \item \code{$proj_F_opt} (vector), integers specifying how to configure each year of the projection: 1: use terminal F, 2: use average F, 3: use F at X\% SPR, 4: use specified F, 5: use specified catch, 6: use Fmsy. Overrides any of the above specifications.
#'     \item \code{$proj_Fcatch} (vector or matrix), catch or F values to use each projection year: values are not used when using Fmsy, FXSPR, terminal F or average F. Overrides any of the above specifications of proj.F or proj.catch. if vector, total catch or F is supplied else matrix columns should be fleets for fleet-specific F to be found/used (\code{n.yrs} x 1 or n_fleets).
#'     \item \code{$proj_mature} (array), user-supplied maturity values for the projection years with dimensions (n_stocks x \code{n.yrs} x n_ages).
#'     \item \code{$proj_waa} (3-d array), user-supplied waa values for the projection years with first and third dimensions equal to that of \code{model$input$data$waa} (waa source x \code{n.yrs} x n_ages).
#'     \item \code{$proj_R_opt} (integer), 1: continue any RE processes for recruitment, 2: make projected recruitment consistent with average recruitment in SPR reference points and cancel any bias correction for NAA in projection years. 3: average recruitment deviations over $avg.yrs.R (if $sigma = "rec") 4: no recruitment deviations (if $sigma = "rec").
#'     \item \code{$proj_NAA_opt} (integer), 1: continue any RE processes for NAA, 2: average NAA deviations over $avg.yrs.NAA. 3: no NAA deviations.
#'     \item \code{$proj_NAA_init} (scalar), the default starting value for all NAA random effects in projection years is exp(10), which may not be large enough for some catch specification. Use this to change the default if a call to project_wham suggests it.
#'     \item \code{$proj_F_init} which F to initialize internal newton search for annual projected F for a given user-specifed catch. Default is 0.1
#'     \item \code{$avg.yrs.sel} list (length = n_fleets), years to average selectivity or FAA for each fleet for projection years. Any BRPs calculated in projection years will also use this. Default = last 5 years, \code{tail(model$years, 5)}.
#'     \item \code{$avg.yrs.waacatch} list (length = n_fleets), years to average weight at age for each fleet for projection years (if $proj_waa is NULL). Any BRPs calculated in projection years will also use this. Default = last 5 years, \code{tail(model$years, 5)}.
#'     \item \code{$avg.yrs.waassb} list (length = n_stocks), years to average weight at age for each stock SSB for projection years (if $proj_waa is NULL). Any BRPs calculated in projection years will also use this. Default = last 5 years, \code{tail(model$years, 5)}.
#'     \item \code{$avg.yrs.mature} list (length = n_stocks), years to average maturity at age for each stock for projection years (if $proj_mature is NULL). Any BRPs calculated in projection years will also use this. Default = last 5 years, \code{tail(model$years, 5)}.
#'     \item \code{$avg.yrs.L} list (length = n_regions), years to average extra mortality at age for each region for projection years (if $cond.L.re = FALSE). Any BRPs calculated in projection years will also use this. Default = last 5 years, \code{tail(model$years, 5)}.
#'     \item \code{$avg.yrs.M} list (length = n_stocks, each is a list with length = n_regions), years to average natural mortality at age for each stock and region for projection years (if $cont.M.re = FALSE). Any BRPs calculated in projection years will also use this. Default = last 5 years, \code{tail(model$years, 5)}.
#'     \item \code{$avg.yrs.move} list (length = n_stocks, each is a list with length = n_regions), years to average movement rates at age and season for each stock and region (at beginning of interval) for projection years (if $cont.move.re = FALSE). Any BRPs calculated in projection years will also use this. Default = last 5 years, \code{tail(model$years, 5)}.
#'     \item \code{$avg.yrs.R} list (length = n_stocks), years to average recruitment deviations for each stock and region for projection years (if $proj_R_opt = 3). Any BRPs calculated in projection years will also use this. Default = last 5 years, \code{tail(model$years, 5)}.
#'     \item \code{$avg.yrs.NAA} list (length = n_stocks, each is a list with length = n_regions), years to average NAA deviations for each stock and region for projection years (if $proj_NAA_opt = 2). Any BRPs calculated in projection years will also use this. Default = last 5 years, \code{tail(model$years, 5)}.
#'   }
#' @param n.newton integer, number of additional Newton steps after optimization. Passed to \code{\link{fit_tmb}}. Default = \code{0} for projections.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? See \code{\link[TMB]{sdreport}}. Default = \code{TRUE}.
#' @param MakeADFun.silent T/F, Passed to silent argument of \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}. Default = \code{FALSE}.
#' @param save.sdrep T/F, save the full \code{\link[TMB]{TMB::sdreport}} object? If \code{FALSE}, only save \code{\link[TMB:summary.sdreport]{summary.sdreport}} to reduce model object file size. Default = \code{TRUE}.
#' @param check.version T/F check whether version WHAM and TMB for fitted model match that of the version of WHAM using for projections. Default = \code{TRUE}.
#' @param TMB.bias.correct T/F whether to use the bias.correct feature of TMB::sdreport. Default = \code{FALSE}.
#' @param TMB.jointPrecision T/F whether to return the joint precision matrix for the fixed and random effects from TMB::sdreport. Default = \code{FALSE}.
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
#' data("input4_SNEMAYT") # load SNEMA yellowtail flounder input data and model settings
#' mod <- fit_wham(input4_SNEMAYT) # using default values (do.proj=T)
#'
#' mod2 <- fit_wham(input4_SNEMAYT, do.retro=F, do.osa=F, do.proj=F) # fit model without projections, retro analysis, or OSA residuals
#' mod_proj <- project_wham(mod2) # add projections to previously fit model, using default values: use.lastF = TRUE, n.yrs = 3, avg.yrs = last 5 years
#'
#' names(mod_proj$rep) # list of derived quantities
#' tail(mod_proj$rep$SSB, 3) # get 3-year projected SSB estimates (weight, not numbers)
#'
#' x <- summary(mod_proj$sdrep)
#' unique(rownames(x)) # list of estimated parameters and derived quanitites with SE
#' x <- x[rownames(x) == "log_SSB",] # SSB estimates with SE
#' ssb.mat <- exp(cbind(x, x[,1] + qnorm(0.975)*cbind(-x[,2],x[,2])))/1000 # calculate 95% CI
#' colnames(ssb.mat) <- c("SSB","SSB_se","SSB_lower","SSB_upper")
#' tail(ssb.mat, 3) # 3-year projected SSB estimates with SE and 95% CI
#' }
project_wham <- function(model, 
  proj.opts=list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE, use.FXSPR=FALSE, use.FMSY=FALSE,
    cont.ecov=TRUE, use.last.ecov=FALSE, percentFXSPR=100, percentFMSY=100),
  n.newton=3, do.sdrep=TRUE, MakeADFun.silent=FALSE, save.sdrep=TRUE, check.version = TRUE, TMB.bias.correct=FALSE, TMB.jointPrecision = FALSE)
{
  if(check.version) verify_version(model)
  # modify wham input (fix parameters at previously estimated values, pad with NAs)
  tryCatch(input2 <- prepare_projection(model, proj.opts, check.version = check.version)
    , error = function(e) {model$err_proj <<- conditionMessage(e)})

  if("err_proj" %in% names(model)) stop(model$err_proj)
  else{# refit model to estimate derived quantities in projection years
  #if(!exists("err")) 
    tryCatch(proj_mod <- TMB::MakeADFun(input2$data, input2$par, DLL = "wham", random = input2$random, map = input2$map, silent = MakeADFun.silent),
      error = function(e) {model$err_MakeADFun <<- conditionMessage(e)})
    if(!is.null(proj_mod$err_MakeADFun)) stop(model$err_proj)
    proj_mod$years <- input2$years
    proj_mod$years_full <- input2$years_full
    proj_mod$ages.lab <- input2$ages.lab
    proj_mod$model_name <- input2$model_name
    proj_mod$call <- match.call()
    proj_mod$input <- input2
    wham_commit <- packageDescription("wham")$GithubSHA1
    proj_mod$wham_commit <- ifelse(is.null(wham_commit), "local install", paste0("Github (timjmiller/wham@", wham_commit, ")")) 
    wham_version <- packageDescription("wham")$Version
    proj_mod$wham_version <- paste0(wham_version, " / ", proj_mod$wham_commit)
    TMB_commit <- packageDescription("TMB")$GithubSHA1
    proj_mod$TMB_commit <- ifelse(is.null(TMB_commit), "local install", paste0("Github (kaskr/adcomp@", TMB_commit, ")")) 
    TMB_version <- packageDescription("TMB")$Version
    proj_mod$TMB_version <- paste0(TMB_version, " / ", proj_mod$TMB_commit, ")")
    
    #If model has not been fitted (i.e., for setting up an operating model/mse), then we do not want to find the Emp. Bayes Posteriors for the random effects.
    is.fit <- !is.null(model$opt)
    if(!is.fit) mle <- model$par
    else {
      mle <- model$opt$par
      proj_mod$fn(mle)
    }
    proj_mod$marg_nll <- proj_mod$fn(mle) #to make sure it is the same as the base model
    if(is.fit & !is.na(proj_mod$marg_nll)) if(abs(proj_mod$marg_nll - model$opt$obj)>1e-5) warning(paste0("Difference between projection model nll and base model nll is ", proj_mod$marg_nll - model$opt$obj))
    proj_mod$rep <- proj_mod$report(proj_mod$env$last.par.best)
    proj_mod$parList <- proj_mod$env$parList(x=mle)
    proj_mod <- check_projF(proj_mod) #projections added.
    if(is.fit & do.sdrep) # only do sdrep if no error and the model has been previously fitted.
    {
      if(model$is_sdrep){ #much faster than old way below
        proj_mod$sdrep <- try(TMB::sdreport(proj_mod, par.fixed = mle, hessian.fixed = solve(model$sdrep$cov.fixed), bias.correct = TMB.bias.correct, getJointPrecision = TMB.jointPrecision))
      } else {
        proj_mod$sdrep <- try(TMB::sdreport(proj_mod, bias.correct = TMB.bias.correct, getJointPrecision = TMB.jointPrecision))
      }
      proj_mod$is_sdrep <- !is.character(proj_mod$sdrep)
      proj_mod$na_sdrep <- ifelse(proj_mod$is_sdrep, any(is.na(summary(proj_mod$sdrep,"fixed")[,2])), NA)
      # if(proj_mod$is_sdrep) proj_mod$na_sdrep <- any(is.na(summary(proj_mod$sdrep,"fixed")[,2])) else proj_mod$na_sdrep = NA
      if(!save.sdrep) proj_mod$sdrep <- summary(proj_mod$sdrep) # only save summary to reduce model object size
    } else {
      proj_mod$is_sdrep <- FALSE
      proj_mod$na_sdrep <- NA
    }
  }

  # pass along previously calculated retros, OSA residuals, error messages, and runtime
  noproj_elements <- names(model)[!names(model) %in% names(proj_mod)] #if anything, should just be OSA.aggregate and OSA.agecomp
  
  proj_mod[noproj_elements] <- model[noproj_elements]
  proj_mod[c("years","years_full","ages.lab")] <- proj_mod$input[c("years","years_full","ages.lab")]
  proj_mod$date <- Sys.time()
  if(!do.sdrep) proj_mod$sdrep <- NULL #remove sdrep of unprojected model

  # print error message
  if(!is.null(model$err_proj))
  {
    proj_mod$err_proj <- model$err_proj
    warning(paste("","** Error during projections. **",
      paste0("Check for issues with proj.opts, see ?project_wham."),"",mod$err_proj,"",sep='\n'))
  }
  return(proj_mod)
}