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
#'   Default: \code{osa.opts = list(method="cdf", parallel=TRUE)}. Discrete versions used for multinomial and Dirichlet-multinomial age composition observations.
#' @param do.post.samp T/F, obtain sample from posterior of random effects? Default = \code{TRUE}. NOT YET IMPLEMENTED.
#' @param model (optional), a previously fit wham model.
#' @param do.check T/F, check if model parameters are identifiable? Passed to \code{\link{fit_tmb}}. Runs internal function \code{check_estimability}, originally provided by https://github.com/kaskr/TMB_contrib_R/TMBhelper. Default = \code{TRUE}.
#' @param MakeADFun.silent T/F, Passed to silent argument of \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}. Default = \code{FALSE}.
#' @param retro.silent T/F, Passed to argument of internal retro function. Determines whether peel number is printed to screen. Default = \code{FALSE}.
#' @param do.proj T/F, do projections? Default = \code{FALSE}. If true, runs \code{\link{project_wham}}.
#' @param proj.opts a named list with the following components:
#'   \itemize{
#'     \item \code{$n.yrs} (integer), number of years to project/forecast. Default = \code{3}.
#'     \item \code{$use.last.F} (T/F), use terminal year F for projections. Default = \code{TRUE}.
#'     \item \code{$use.FXSPR} (T/F), calculate and use F at X% SPR for projections.
#'     \item \code{$use.FMSY} (T/F), calculate and use FMSY for projections.
#'     \item \code{$proj.F} (vector), user-specified fishing mortality for projections. Length must equal \code{n.yrs}.
#'     \item \code{$proj.catch} (vector), user-specified aggregate catch for projections. Length must equal \code{n.yrs}.
#'     \item \code{$avg.yrs} (vector), specify which years to average over for calculating reference points. Default = last 5 model years, \code{tail(model$years, 5)}.
#'     \item \code{$cont.ecov} (T/F), continue ecov process (e.g. random walk or AR1) for projections. Default = \code{TRUE}.
#'     \item \code{$use.last.ecov} (T/F), use terminal year ecov for projections.
#'     \item \code{$avg.ecov.yrs} (vector), specify which years to average over the environmental covariate(s) for projections.
#'     \item \code{$proj.ecov} (vector), user-specified environmental covariate(s) for projections. Length must equal \code{n.yrs}.
#'     \item \code{$cont.Mre} (T/F), continue M random effects (i.e. AR1_y or 2D AR1) for projections. Default = \code{TRUE}. If \code{FALSE}, M will be averaged over \code{$avg.yrs} (which defaults to last 5 model years).
#'     \item \code{$avg.rec.yrs} (vector), specify which years to calculate the CDF of recruitment for use in projections. Default = all model years.
#'     \item \code{$percentFXSPR} (scalar), percent of F_XSPR to use for calculating catch in projections, only used if proj.opts$use.FXSPR = TRUE. For example, GOM cod uses F = 75% F_40%SPR, so \code{proj.opts$percentFXSPR = 75}. Default = 100.
#'     \item \code{$percentFMSY} (scalar), percent of F_MSY to use for calculating catch in projections, only used if $use.FMSY = TRUE.
#'   }
#' @param do.fit T/F, fit the model using \code{fit_tmb}. Default = \code{TRUE}.
#' @param save.sdrep T/F, save the full \code{\link[TMB]{TMB::sdreport}} object? If \code{FALSE}, only save \code{\link[TMB:summary.sdreport]{summary.sdreport}} to reduce model object file size. Default = \code{TRUE}.
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
#' data("input4_SNEMAYT") # load SNEMA yellowtail flounder data and parameter settings
#' mod = fit_wham(input4_SNEMAYT) # using default values
#' mod = fit_wham(input4_SNEMAYT, do.retro=FALSE, osa.opts=list(method="oneStepGeneric")) # slower OSA method. 
#'
#' names(mod$rep) # list of derived quantities
#' mod$rep$SSB # get SSB estimates (weight, not numbers)
#' m1$rep$NAA[,1] # get recruitment estimates (numbers, first column of numbers-at-age matrix)
#' m1$rep$F[,1] # get F estimates for fleet 1
#' }
fit_wham = function(input, n.newton = 3, do.sdrep = TRUE, do.retro = TRUE, n.peels = 7,
                    do.osa = TRUE, osa.opts = list(method="cdf", parallel=TRUE), do.post.samp = TRUE,
                    model=NULL, do.check = FALSE, MakeADFun.silent=FALSE, retro.silent = FALSE, do.proj = FALSE,
                    proj.opts=list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE, use.FXSPR=FALSE, proj.F=NULL, 
                      proj.catch=NULL, avg.yrs=NULL, cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, 
                      proj.ecov=NULL, cont.Mre=NULL, avg.rec.yrs=NULL, percentFXSPR=100),
                    do.fit = TRUE, save.sdrep=TRUE)
{

  # fit model
  if(missing(model)){
    mod <- TMB::MakeADFun(input$data, input$par, DLL = "wham", random = input$random, map = input$map, silent = MakeADFun.silent)
  } else {mod = model}

  mod$years <- input$years
  mod$years_full <- input$years_full
  mod$ages.lab <- input$ages.lab
  mod$model_name <- input$model_name
  mod$input <- input
  ver <- sessioninfo::package_info() %>% as.data.frame %>% dplyr::filter(package=="wham") %>% dplyr::select(loadedversion, source) %>% unname
  mod$wham_version <- paste0(ver, collapse=" / ")
  if(do.fit){
    btime <- Sys.time()
    #mod <- fit_tmb(mod, n.newton = n.newton, do.sdrep = do.sdrep, do.check = do.check, save.sdrep = save.sdrep)
    mod <- fit_tmb(mod, n.newton = n.newton, do.sdrep = FALSE, do.check = do.check, save.sdrep = save.sdrep)
    mod$runtime <- round(difftime(Sys.time(), btime, units = "mins"),2) # don't count retro or proj in runtime
    mod <- check_which_F_age(mod)
    mod <- check_FXSPR(mod)
    if(mod$env$data$do_proj==1) mod <- check_projF(mod) #projections added.
    if(do.sdrep) mod <- do_sdrep(mod, save.sdrep = save.sdrep)

    # retrospective analysis
    if(do.retro){
      tryCatch(mod$peels <- retro(mod, ran = unique(names(mod$env$par[mod$env$random])), n.peels= n.peels,
        MakeADFun.silent = MakeADFun.silent, retro.silent = retro.silent), error = function(e) {mod$err_retro <<- conditionMessage(e)})
      #assigning mod$err_retro accomplishes the below if statement
      #if(exists("err")){
      #  mod$err_retro <- err # store error message to print out in fit_wham
      #  rm("err")
      #}
    }

    # one-step-ahead residuals
    if(do.osa){
      if(mod$is_sdrep){ # only do OSA residuals if sdrep ran
        cat("Doing OSA residuals...\n");
        mod$env$data$do_osa = 1
        full_set = 1:length(input$data$obsvec)
        input$data$obs$residual = NA
        if(!is.null(input$data$condition_no_osa)) cat("OSA not available for some age comp likelihoods...\n")
        #first do continuous obs, condition on obs without osa (probably none)
        subset. = setdiff(full_set, c(input$data$subset_discrete_osa, input$data$conditional_no_osa))
        OSA.continuous <- suppressWarnings(TMB::oneStepPredict(obj=mod, observation.name="obsvec",
                                    data.term.indicator="keep",
                                    method=osa.opts$method,
                                    discrete=FALSE, parallel=osa.opts$parallel,
                                    subset = subset., conditional = input$data$conditional_no_osa))
        input$data$obs$residual[subset.] <- OSA.continuous$residual;
        mod$OSA.continuous = OSA.continuous
        if(!is.null(input$data$subset_discrete_osa)) {
          cat("Doing OSA for discrete age comp likelihoods...\n")
          conditional = union(input$data$condition_no_osa, subset.) #all with continuous and without osa 
          subset. = input$data$subset_discrete_osa
          #first do continuous
          OSA.discrete <- suppressWarnings(TMB::oneStepPredict(obj=mod, observation.name="obsvec",
                                      data.term.indicator="keep",
                                      method= osa.opts$method,
                                      discrete=TRUE, parallel=osa.opts$parallel,
                                      conditional = conditional))
          input$data$obs$residual[subset.] <- OSA.discrete$residual;
          mod$OSA.discrete = OSA.discrete
        }
        mod$osa <- input$data$obs
        mod$env$data$do_osa = 0 #set this back to not using OSA likelihoods
      } else warning(paste("","** Did not do OSA residual analyses. **",
      "Error during TMB::sdreport(). Check for unidentifiable parameters.","",sep='\n'))
    }

    # projections, calls prepare_projection + fit_wham(do.proj=F)
    if(do.proj) mod <- project_wham(mod, proj.opts=proj.opts, MakeADFun.silent = MakeADFun.silent, do.sdrep = do.sdrep, save.sdrep = save.sdrep)

    # error message reporting
    if(!is.null(mod$err)) warning(paste("","** Error during model fit. **",
      "Check for unidentifiable parameters.","",mod$err,"",sep='\n'))
    if(!is.null(mod$err_retro)) warning(paste("","** Error during retrospective analysis. **",
      paste0("Check for issues with last ",n.peels," model years."),"",mod$err_retro,"",sep='\n'))
  }
  else { #model not fit, but generate report and parList so project_wham can be used without fitted model.
    mod$rep = mod$report() #par values don't matter because function has not been evaluated
    mod$parList = mod$env$parList()
    mod <- check_which_F_age(mod)
    mod <- check_FXSPR(mod)
  }

  return(mod)
}


check_which_F_age = function(mod)
{
  mod$env$data$which_F_age[] = apply(mod$rep$FAA_tot,1, function(x) which(x == max(x))[1])
  mod$rep = mod$report()
  return(mod)
}
do_sdrep = function(model, save.sdrep = TRUE)
{
  model$sdrep <- try(TMB::sdreport(model))
  model$is_sdrep <- !is.character(model$sdrep)
  if(model$is_sdrep) model$na_sdrep <- any(is.na(summary(model$sdrep,"fixed")[,2])) else model$na_sdrep = NA
  if(!save.sdrep) model$sdrep <- summary(model$sdrep) # only save summary to reduce model object size
  return(model)
}

check_FXSPR = function(mod)
{
  #If model has not been fitted (i.e., for setting up an operating model/mse), then we do not want to find the Emp. Bayes Posteriors for the random effects.
  if(is.null(mod$opt)) mle = mod$par
  else {
    mle = mod$opt$par
  }
  percentSPR_out = exp(mod$rep$log_SPR_FXSPR - mod$rep$log_SPR0)
  ind = which(round(percentSPR_out,4) != round(mod$env$data$percentSPR/100,4))
  if(length(ind))
  {
    for(i in 1:2) #two tries to fix initial FXSPR value
    {
      if(mod$env$data$n_years_proj) years = mod$years_full
      else years = mod$years
      redo_SPR_years = years[ind]
      warning(paste0("Changing initial values for estimating FXSPR for years ", paste(redo_SPR_years, collapse = ","), "."))
      mod$env$data$FXSPR_init[ind] = mod$env$data$FXSPR_init[ind]*0.5
      mod$fn(mle)
      mod$rep = mod$report()
      percentSPR_out = exp(mod$rep$log_SPR_FXSPR - mod$rep$log_SPR0)
      ind = which(round(percentSPR_out,4) != round(mod$env$data$percentSPR/100,4))
      if(!length(ind)) break
    }
  }
  if(length(ind)) warning(paste0("Still bad initial values and estimates of FXSPR for years ", paste(years[ind], collapse = ","), "."))

  percentSPR_out_static = exp(mod$rep$log_SPR_FXSPR_static - mod$rep$log_SPR0_static)
  ind = which(round(percentSPR_out_static,4) != round(mod$env$data$percentSPR/100,4))
  if(length(ind))
  {
    for(i in 1:2) #two tries to fix initial FXSPR value
    {
      warning(paste0("Changing initial values for estimating static FXSPR."))
      mod$env$data$static_FXSPR_init = mod$env$data$static_FXSPR_init*0.5
      mod$fn(mle)
      mod$rep = mod$report()
      percentSPR_out_static = exp(mod$rep$log_SPR_FXSPR_static - mod$rep$log_SPR0_static)
      ind = which(round(percentSPR_out_static,4) != round(mod$env$data$percentSPR/100,4))
      if(!length(ind)) break
    }
  }
  if(length(ind)) warning(paste0("Still bad initial values and estimates of static FXSPR."))

  return(mod)     
}

check_projF = function(mod)
{
  if(is.null(mod$opt)) mle = mod$par
  else {
    mle = mod$opt$par
  }
  proj_F_opt = mod$env$data$proj_F_opt
  ind = which(proj_F_opt == 3) #FXSPR
  if(length(ind))
  {
    y = mod$env$data$n_years_model + ind
    correct_F = round(mod$env$data$percentFXSPR * exp(mod$rep$log_FXSPR[y])/100, 4)
    used_F = round(mod$rep$FAA_tot[cbind(y,mod$env$data$which_F_age[y])],4)
    bad = which(correct_F != used_F)
    if(length(bad))
    {
      redo_SPR_years = mod$years_full[y[bad]]
      warning(paste0("Changing initial values for estimating FXSPR used to define F in projection years ", paste(redo_SPR_years, collapse = ","), "."))
      mod$env$data$F_init_proj[ind[bad]] = mod$env$data$FXSPR_init_proj[y[bad]]
      mod$fn(mle)
      mod$rep = mod$report()
      correct_F = round(mod$env$data$percentFXSPR * exp(mod$rep$log_FXSPR[y])/100, 4)
      used_F = round(mod$rep$FAA_tot[cbind(y,mod$env$data$which_F_age[y])],4)
      bad = which(correct_F != used_F)
    }
    y_bad_FXSPR = mod$years_full[y[bad]]
    if(length(bad)) warning(paste0("Still bad initial values and estimates of FXSPR used to define F in projection years ", paste(y_bad_FXSPR, collapse = ","), "."))
  }
  ind = which(proj_F_opt == 5) #Find F from catch
  if(length(ind))
  {
    y = mod$env$data$n_years_model + ind
    #print(y)
    #print(dim(mod$rep$pred_catch))
    bad = which(round(mod$env$data$proj_Fcatch[ind],4) != round(rowSums(mod$rep$pred_catch[y,,drop=F]),4))
    if(length(bad))
    {
      for(i in 1:2)
      {
        redo_Catch_years = mod$years_full[y[bad]]
        #print(redo_Catch_years)
        warning(paste0("Changing initial values for finding F from Catch in projection years ", paste(redo_Catch_years, collapse = ","), , "."))
        mod$env$data$F_init_proj[ind[bad]] = mod$env$data$F_init_proj[ind[bad]]*0.5
        mod$fn(mle)
        mod$rep = mod$report()
        bad = which(round(mod$env$data$proj_Fcatch[ind],4) != round(sum(mod$rep$pred_catch[y,,drop=F]),4))
        if(!length(bad)) break
      }
    }
    y_bad_Fcatch = mod$years_full[y[bad]]
    if(length(bad)) warning(paste0("Still bad initial values for finding F from Catch in projection years ", paste(y_bad_Fcatch, collapse = ","), , "."))
  }
  return(mod)
}
