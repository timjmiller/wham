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
#' @param osa.opts list of 2 options (method, parallel) for calculating OSA residuals, passed to \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}.
#'   Default: \code{osa.opts = list(method="oneStepGaussianOffMode", parallel=TRUE)}. See \code{\link{make_osa_residuals}}.
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
#'     \item \code{$proj.ecov} (matrix), user-specified environmental covariate(s) for projections. \code{n.yrs x n_Ecov}.
#'     \item \code{$cont.Mre} (T/F), continue M random effects (i.e. AR1_y or 2D AR1) for projections. Default = \code{TRUE}. If \code{FALSE}, M will be averaged over \code{$avg.yrs} (which defaults to last 5 model years).
#'     \item \code{$avg.rec.yrs} (vector), specify which years to calculate the CDF of recruitment for use in projections. Default = all model years.
#'     \item \code{$percentFXSPR} (scalar), percent of F_XSPR to use for calculating catch in projections, only used if proj.opts$use.FXSPR = TRUE. For example, GOM cod uses F = 75% F_40%SPR, so \code{proj.opts$percentFXSPR = 75}. Default = 100.
#'     \item \code{$percentFMSY} (scalar), percent of F_MSY to use for calculating catch in projections, only used if $use.FMSY = TRUE.
#'   }
#' @param do.fit T/F, fit the model using \code{fit_tmb}. Default = \code{TRUE}.
#' @param save.sdrep T/F, save the full \code{\link[TMB]{TMB::sdreport}} object? If \code{FALSE}, only save \code{\link[TMB:summary.sdreport]{summary.sdreport}} to reduce model object file size. Default = \code{TRUE}.
#' @param do.brps T/F, calculate and report biological reference points. Default = \code{TRUE}.
#' @param fit.tmb.control list of optimizer controlling attributes passed to \code{\link[wham]{fit_tmb}}. Default is \code{list(use.optim = FALSE, opt.control = list(iter.max = 1000, eval.max = 1000))}, so stats::nlminb is used to opitmize.
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
fit_wham <- function(input, n.newton = 3, do.sdrep = TRUE, do.retro = TRUE, n.peels = 7,
                    do.osa = TRUE, osa.opts = list(method="oneStepGaussianOffMode", parallel=TRUE), do.post.samp = TRUE,
                    model=NULL, do.check = FALSE, MakeADFun.silent=FALSE, retro.silent = FALSE, do.proj = FALSE,
                    proj.opts=list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE, use.FXSPR=FALSE, proj.F=NULL, 
                      proj.catch=NULL, avg.yrs=NULL, cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, 
                      proj.ecov=NULL, cont.Mre=NULL, avg.rec.yrs=NULL, percentFXSPR=100),
                    do.fit = TRUE, save.sdrep=TRUE, do.brps = TRUE, fit.tmb.control = NULL)
{

  # fit model
  if(missing(model)){
    mod <- TMB::MakeADFun(input$data, input$par, DLL = "wham", random = input$random, map = input$map, silent = MakeADFun.silent)
  } else {
    verify_version(model)
    mod <- model
  }

  mod$years <- input$years
  mod$years_full <- input$years_full
  mod$ages.lab <- input$ages.lab
  mod$model_name <- input$model_name
  mod$input <- input
  mod$call <- match.call()
  mod <- check_which_F_age(mod) #can be an issue if estimated full F is at age with 0 selectivity
  #mod$rep <- mod$report() #not needed because check_which_F_age calls mod$report()
  #ver <- sessioninfo::package_info() %>% as.data.frame %>% dplyr::filter(package=="wham") %>% dplyr::select(loadedversion, source) %>% unname
  #mod$wham_version <- paste0(ver, collapse=" / ")
  wham_commit <- packageDescription("wham")$GithubSHA1
  mod$wham_commit <- ifelse(is.null(wham_commit), "local install", paste0("Github (timjmiller/wham@", wham_commit, ")")) 
  wham_version <- packageDescription("wham")$Version
  mod$wham_version <- paste0(wham_version, " / ", mod$wham_commit)

  if(do.fit){
    btime <- Sys.time()
    #mod <- fit_tmb(mod, n.newton = n.newton, do.sdrep = do.sdrep, do.check = do.check, save.sdrep = save.sdrep)
    use.optim <- FALSE
    opt.control <- NULL
    if(!is.null(fit.tmb.control$use.optim)) use.optim <- fit.tmb.control$use.optim
    if(!is.null(fit.tmb.control$opt.control)) opt.control <- fit.tmb.control$opt.control
    mod <- fit_tmb(mod, n.newton = n.newton, do.sdrep = FALSE, do.check = do.check, save.sdrep = save.sdrep, use.optim=use.optim, opt.control = opt.control)
    mod$runtime <- round(difftime(Sys.time(), btime, units = "mins"),2) # don't count retro or proj in runtime
    if(do.brps){
      mod <- do_reference_points(mod)
      # if(any(mod$input$data$can_move==1) & any(mod$input$data$mig_type == 1)){
      #   warning("Cannot currently calculate standard errors of biological reference points internally when survival and movement are simultaneous for any stock.")
      # } else {
      #   mod <- check_which_F_age(mod)
      #   mod$input$data$do_SPR_BRPs <- mod$env$data$do_SPR_BRPs <- 1
      #   if(any(input$data$recruit_model %in% 3:4)) mod$input$data$do_MSY_BRPs <- mod$env$data$do_MSY_BRPs <- 1
      #   mod$rep = mod$report() #par values don't matter because function has not been evaluated
      #   mod <- check_FXSPR(mod)
      # }
    }
    
    if(mod$env$data$do_proj==1) mod <- check_projF(mod) #projections added.
    if(do.sdrep) mod <- do_sdreport(mod, save.sdrep = save.sdrep)

    # retrospective analysis
    if(do.retro){
      tryCatch(mod$peels <- retro(mod, n.peels= n.peels,
        MakeADFun.silent = MakeADFun.silent, retro.silent = retro.silent), error = function(e) {mod$err_retro <<- conditionMessage(e)})
    }
    
    # one-step-ahead residuals
    if(do.osa & mod$is_sdrep){
      mod <- make_osa_residuals(mod, osa.opts = osa.opts)
    } else if(do.osa) warning(paste("","** Did not do OSA residual analyses. **",
        "If do.sdrep = TRUE, then there was an error during TMB::sdreport(), and so should check for unidentifiable parameters.","",
        sep='\n'))

    # projections, calls prepare_projection + fit_wham(do.proj=F)
    if(do.proj) mod <- project_wham(mod, proj.opts=proj.opts, MakeADFun.silent = MakeADFun.silent, do.sdrep = do.sdrep, save.sdrep = save.sdrep)

    # error message reporting
    if(!is.null(mod$err)) warning(paste("","** Error during model fit. **",
      "Check for unidentifiable parameters.","",mod$err,"",sep='\n'))
    if(!is.null(mod$err_retro)) warning(paste("","** Error during retrospective analysis. **",
      paste0("Check for issues with last ",n.peels," model years."),"",mod$err_retro,"",sep='\n'))
  } else { #model not fit, but generate report and parList so project_wham can be used without fitted model.
    if(do.brps){
      mod <- do_reference_points(mod)
      # mod$input$data$do_SPR_BRPs <- mod$env$data$do_SPR_BRPs <- 1
      # if(any(input$data$recruit_model %in% 3:4)) input$data$do_MSY_BRPs <- mod$env$data$do_MSY_BRPs <- 1
    }
    mod$parList <- mod$env$parList()
    mod <- check_which_F_age(mod)
    #mod$rep <- mod$report() #not needed because check_which_F_age calls mod$report()
    if(is.null(mod$TMB_commit)){
      TMB_commit <- packageDescription("TMB")$GithubSHA1
      mod$TMB_commit <- ifelse(is.null(TMB_commit), "local install", paste0("Github (kaskr/adcomp@", TMB_commit, ")")) 
    }
    if(is.null(mod$TMB_version)){
      TMB_version <- packageDescription("TMB")$Version
      mod$TMB_version <- paste0(TMB_version, " / ", mod$TMB_commit, ")")
    }
    if(is.null(mod$is_sdrep)) mod$is_sdrep = FALSE

    #mod <- check_FXSPR(mod)
  }

  return(mod)
}


check_which_F_age <- function(mod)
{
  if(is.null(mod$opt)) mle <- mod$par
  else {
    mle <- mod$opt$par
  }
  mod$fn(mle)
  mod$rep <- mod$report()
  for(y in 1:dim(mod$rep$FAA)[2]){
    temp <- apply(rbind(mod$rep$FAA[,y,]),2,sum)
    mod$env$data$which_F_age[y] <- mod$input$data$which_F_age[y] <- which(temp == max(temp))[1]
  }
  if(mod$env$data$do_SPR_BRPs == 1){
    temp <- apply(mod$rep$FAA_static,2,sum)
    mod$env$data$which_F_age_static <- mod$input$data$which_F_age_static <- as.numeric(which(temp == max(temp))[1])
  }
  mod$retape()
  mod$fn(mle)
  mod$rep <- mod$report()
  return(mod)
}

# do_sdreport <- function(model, save.sdrep = TRUE) {
#   model$sdrep <- try(TMB::sdreport(model))
#   model$is_sdrep <- !is.character(model$sdrep)
#   if(model$is_sdrep) model$na_sdrep <- any(is.na(summary(model$sdrep,"fixed")[,2])) else model$na_sdrep = NA
#   if(!save.sdrep) model$sdrep <- summary(model$sdrep) # only save summary to reduce model object size
#   return(model)
# }

check_FXSPR <- function(mod)
{
  #If model has not been fitted (i.e., for setting up an operating model/mse), then we do not want to find the Emp. Bayes Posteriors for the random effects.
  if(is.null(mod$opt)) mle <- mod$par
  else {
    mle <- mod$opt$par
  }
  #allow calculation of X%SPR BRPs
  if(mod$env$data$do_SPR_BRPs != 1) stop("in check_FXSPR: model$env$data$do_SPR_BRPs is not equal to 1")
  #mod$env$data$do_SPR_BRPs <- mod$input$data$do_SPR_BRPs <- 1
  mod$retape()
  mod$fn(mle)
  mod$rep <- mod$report()


  percentSPR_out <- exp(cbind(mod$rep$log_SPR_FXSPR - mod$rep$log_SPR0)[,mod$input$data$n_stocks+1])
  ind = which(round(percentSPR_out,4) != round(mod$env$data$percentSPR/100,4))
  years <- mod$years
  if(length(ind))
  {
    for(i in 1:2) #two tries to fix initial FXSPR value
    {
      redo_SPR_years <- years[ind]
      warning(paste0("Changing initial values for estimating FXSPR for years ", paste(redo_SPR_years, collapse = ","), "."))
      mod$env$data$FXSPR_init[ind] <- mod$input$data$FXSPR_init[ind] <- mod$env$data$FXSPR_init[ind]*2
      mod$retape()
      mod$fn(mle)
      mod$rep <- mod$report()
      percentSPR_out <- exp(cbind(mod$rep$log_SPR_FXSPR - mod$rep$log_SPR0)[,mod$input$data$n_stocks+1])
      ind <- which(round(percentSPR_out,4) != round(mod$env$data$percentSPR/100,4))
      if(!length(ind)) break
    }
  }
  if(length(ind)) warning(paste0("Still bad initial values and estimates of FXSPR for years ", paste(years[ind], collapse = ","), "."))

  percentSPR_out_static <- exp(mod$rep$log_SPR_FXSPR_static - mod$rep$log_SPR0_static)[mod$input$data$n_stocks+1]
  ind <- which(round(percentSPR_out_static,4) != round(mod$env$data$percentSPR/100,4))
  if(length(ind))
  {
    for(i in 1:2) #two tries to fix initial FXSPR value
    {
      warning(paste0("Changing initial values for estimating static FXSPR."))
      mod$env$data$FXSPR_static_init <- mod$input$data$FXSPR_static_init <- mod$env$data$FXSPR_static_init*0.5
      mod$retape()
      mod$fn(mle)
      mod$rep <- mod$report()
      percentSPR_out_static <- exp(mod$rep$log_SPR_FXSPR_static - mod$rep$log_SPR0_static)[mod$input$data$n_stocks+1]
      ind = which(round(percentSPR_out_static,4) != round(mod$env$data$percentSPR/100,4))
      if(!length(ind)) break
    }
  }
  if(length(ind)) warning(paste0("Still bad initial values and estimates of static FXSPR."))

  return(mod)
}

check_projF <- function(mod)
{
  if(is.null(mod$opt)) mle <- mod$par
  else {
    mle <- mod$opt$par
  }
  proj_F_opt <- mod$env$data$proj_F_opt
  ind = which(proj_F_opt == 3) #FXSPR
  if(length(ind))
  {
    y <- mod$env$data$n_years_model + ind
    correct_F <- round(mod$env$data$percentFXSPR * exp(mod$rep$log_FXSPR[y])/100, 2)
    FAA_tot <- apply(mod$rep$FAA,2:3, sum)
    used_F <- round(FAA_tot[cbind(y,mod$env$data$which_F_age[y])],2)
    # print(used_F)
    bad <- which(correct_F != used_F)
    if(length(bad))
    {
      redo_SPR_years <- mod$years_full[y[bad]]
      warning(paste0("Changing initial values for estimating FXSPR used to define F in projection years ", paste(redo_SPR_years, collapse = ","), "."))
      mod$env$data$F_proj_init[ind[bad]] <- mod$input$data$F_proj_init[ind[bad]] <- mod$env$data$FXSPR_init[y[bad]]
      mod$retape()
      mod$fn(mle)
      mod$rep <- mod$report()
      correct_F <- round(mod$env$data$percentFXSPR * exp(mod$rep$log_FXSPR[y])/100, 2)
      used_F <- round(FAA_tot[cbind(y,mod$env$data$which_F_age[y])],2)
      bad <- which(correct_F != used_F)
    }
    y_bad_FXSPR <- mod$years_full[y[bad]]
    if(length(bad)) warning(paste0("Still bad initial values and estimates of FXSPR used to define F in projection years ", paste(y_bad_FXSPR, collapse = ","), "."))
  }
  ind <- which(proj_F_opt == 5) #Find F from catch
  if(length(ind))
  {
    y <- mod$env$data$n_years_model + ind
    by_fleet <- NCOL(mod$env$data$proj_Fcatch) == mod$env$data$n_fleets
    if(!by_fleet) bad <- which(round(mod$env$data$proj_Fcatch[ind,1],4) != round(rowSums(mod$rep$pred_catch[y,,drop=F]),4))
    else bad <- which(any(round(mod$env$data$proj_Fcatch[ind,],4) != round(mod$rep$pred_catch[y,],4)))
    if(length(bad))
    {
      # print(mod$marg_nll)
      if(is.na(mod$marg_nll)){
        # print(mod$rep$pred_NAA[1,1,which(!mod$years_full %in% mod$years),])
        FAA_tot <- exp(mod$rep$log_FAA_tot[which(!mod$years_full %in% mod$years),, drop = FALSE])
        y_na <- is.na(apply(FAA_tot,1,sum))
        # print(y_na)
        #if(any(proj_F_opt[y_na]==5)) 
        # print(mod$rep$log_FAA_tot[which(!mod$years_full %in% mod$years),])
        if(any(y_na)) stop("Need to change initial log_NAA parameter values in projection years. See ?project_wham and proj.opts$proj_NAA_init")
      }
      for(i in 1:2)
      {
        redo_Catch_years <- mod$years_full[y[bad]]

        warning(paste0("Changing initial values for finding F from Catch in projection years ", paste(redo_Catch_years, collapse = ","), "."))
        mod$env$data$F_proj_init[ind[bad]] <- mod$input$data$F_proj_init[ind[bad]] <- mod$env$data$F_proj_init[ind[bad]]*0.5
        mod$retape()
        mod$fn(mle)
        mod$rep <- mod$report()
        if(!by_fleet) bad <- which(round(mod$env$data$proj_Fcatch[ind,1],4) != round(rowSums(mod$rep$pred_catch[y,,drop=F]),4))
        else bad <- which(any(round(mod$env$data$proj_Fcatch[ind,],4) != round(mod$rep$pred_catch[y,],4)))
        # bad <- which(round(mod$env$data$proj_Fcatch[ind],4) != round(sum(mod$rep$pred_catch[y,,drop=F]),4))
        if(!length(bad)) break
      }
    }
    y_bad_Fcatch <- mod$years_full[y[bad]]
    if(length(bad)) warning(paste0("Still bad initial values for finding F from Catch in projection years ", paste(y_bad_Fcatch, collapse = ","), "."))
  }
  return(mod)
}