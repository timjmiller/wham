#' Prepare input data and parameters to project an already fit WHAM model
#'
#' \code{prepare_projection} is an internal function called by \code{\link{project_wham}},
#' which in turn is called by \code{\link{fit_wham}} if \code{do.proj = TRUE}.
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
#'     \item \code{$avg.yrs} (vector), specify which years to average over for calculating reference points. Default = last 5 model years, \code{tail(model$years, 5)}.
#'     \item \code{$cont.ecov} (T/F), continue ecov process (e.g. random walk or AR1) for projections. Default = \code{TRUE}.
#'     \item \code{$use.last.ecov} (T/F), use terminal year ecov for projections.
#'     \item \code{$avg.ecov.yrs} (vector), specify which years to average over the environmental covariate(s) for projections.
#'     \item \code{$proj.ecov} (matrix), user-specified environmental covariate(s) for projections. \code{n.yrs x n.ecov}.
#'     \item \code{$cont.M.re} (T/F), continue M random effects (i.e. AR1_y or 2D AR1) for projections. Default = \code{FALSE}. If \code{FALSE}, M will be averaged over \code{$avg.yrs} (which defaults to last 5 model years).
#'     \item \code{$cont.move.re} (T/F), continue any movement random effects for projections. Default = \code{FALSE}. If \code{FALSE}, movement parameters will be averaged over \code{$avg.yrs} (which defaults to last 5 model years).
#'     \item \code{$cont.L.re} (T/F), continue any movement random effects for projections. Default = \code{FALSE}. If \code{FALSE}, movement parameters will be averaged over \code{$avg.yrs} (which defaults to last 5 model years).
#'     \item \code{$avg.rec.yrs} (vector), specify which years to calculate the CDF of recruitment for use in projections. Default = all model years. Only used when recruitment is estimated as fixed effects (SCAA).
#'     \item \code{$percentFXSPR} (scalar), percent of F_XSPR to use for calculating catch in projections, only used if $use.FXSPR = TRUE. For example, GOM cod uses F = 75\% F_40\%SPR, so \code{proj.opts$percentFXSPR = 75}. Default = 100.
#'     \item \code{$percentFMSY} (scalar), percent of F_MSY to use for calculating catch in projections, only used if $use.FMSY = TRUE.
#'     \item \code{$proj_F_opt} (vector), integers specifying how to configure each year of the projection: 1: use terminal F, 2: use average F, 3: use F at X\% SPR, 4: use specified F, 5: use specified catch, 6: use Fmsy. Overrides any of the above specifications.
#'     \item \code{$proj_Fcatch} (vector or matrix), catch or F values to use each projection year: values are not used when using Fmsy, FXSPR, terminal F or average F. Overrides any of the above specifications of proj.F or proj.catch. if vector, total catch or F is supplied else matrix columns should be fleets for fleet-specific F to be found/used (\code{n.yrs} x 1 or n_fleets).
#'     \item \code{$proj_mature} (array), user-supplied maturity values for the projection years with dimensions (n_stocks x \code{n.yrs} x n_ages).
#'     \item \code{$proj_waa} (3-d array), user-supplied waa values for the projection years with first and third dimensions equal to that of \code{model$input$data$waa} (waa source x \code{n.yrs} x n_ages).
#'     \item \code{$proj_R_opt} (integer), 1: continue any RE processes for recruitment, 2: make projected recruitment consistent with average recruitment in SPR reference points and cancel any bias correction for NAA in projection years.
#'     \item \code{$proj_NAA_init} (scalar), the default starting value for all NAA random effects in projection years is exp(10), which may not be large enough for some catch specification. Use this to change the default if a call to project_wham suggests it.
#'   }
#' @param check.version T/F check whether version WHAM and TMB for fitted model match that of the version of WHAM using for projections. Default = \code{TRUE}.
#'
#' @return same as \code{\link{prepare_wham_input}}, a list ready for \code{\link{fit_wham}}:
#'   \describe{
#'     \item{\code{data}}{Named list of data, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{par}}{Named list of parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{map}}{Named list of factors that determine which parameters are estimated, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{random}}{Character vector of parameters to treat as random effects, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{years}}{Numeric vector of representing (non-projection) model years of WHAM model}
#'     \item{\code{years_full}}{Numeric vector of representing all model and projection years of WHAM model}
#'     \item{\code{ages.lab}}{Character vector of age labels, ending with plus-group}
#'     \item{\code{model_name}}{Character, name of stock/model (specified in call to \code{prepare_wham_input})}
#'   }
#'
#' @seealso \code{\link{prepare_wham_input}}, \code{\link{project_wham}}
#'
#'
#' @export
prepare_projection = function(model, proj.opts, check.version=FALSE) {
# if(is.null(proj.opts)) proj.opts=list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE, use.FXSPR=FALSE, use.FMSY=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
#                                       cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL, cont.M.re=NULL, avg.rec.yrs=NULL, percentFXSPR=100,
#                                       percentFMSY=100, proj_F_opt = NULL, proj_Fcatch = NULL)
  if(is.null(proj.opts)) proj.opts <- list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE, use.FXSPR=FALSE, use.FMSY=FALSE,
    cont.ecov=TRUE, use.last.ecov=FALSE, percentFXSPR=100, percentFMSY=100)
  if(check.version) verify_version(model)
  # default: 3 projection years
  # peel <- 0
  # if(!is.null(model$peel)) peel <- model$peel # projecting off of a peel

  input <- model$input
  if(is.null(proj.opts$n.yrs)) proj.opts$n.yrs <- 3
  # default: use average M, selectivity, etc. over last 5 model years to calculate ref points
  if(is.null(proj.opts$avg.yrs)) proj.opts$avg.yrs <- input$years[model$env$data$avg_years_ind+1] #tail(model$years, 5)  
  if(any(proj.opts$avg.yrs<0)) stop("negative years are specified in proj.opts$avg.yrs to average projection inputs")
  if(all(is.null(proj.opts$use.last.F), is.null(proj.opts$use.avg.F), is.null(proj.opts$use.FXSPR), is.null(proj.opts$use.FMSY), is.null(proj.opts$proj.F), is.null(proj.opts$proj.catch))){
    proj.opts$use.last.F=TRUE; proj.opts$use.avg.F=FALSE; proj.opts$use.FXSPR=FALSE; proj.opts$use.FMSY=FALSE; proj.opts$proj.F=NULL; proj.opts$proj.catch=NULL
  }
  if(is.null(proj.opts$percentFXSPR)) proj.opts$percentFXSPR = 100
  if(is.null(proj.opts$percentFMSY)) proj.opts$percentFMSY = 100

  data <- input$data
  # check options for F/catch are valid
  if(any(!(proj.opts$avg.yrs %in% input$years))) stop(paste("","** Error setting up projections: **",
    "proj.opts$avg.yrs is not a subset of model years.","",sep='\n'))
  F.opt.ct <- sum(proj.opts$use.avg.F, proj.opts$use.last.F, proj.opts$use.FXSPR, proj.opts$use.FMSY, !is.null(proj.opts$proj.F), !is.null(proj.opts$proj.catch))
  if(F.opt.ct != 1) stop(paste("","** Error setting up projections: **",
    "Exactly one method of specifying F must be used (see ?project_wham).",
    "You have specified these in 'proj.opts':",
    capture.output(cat("  $use.last.F = ",proj.opts$use.last.F)),
    capture.output(cat("  $use.avg.F = ",proj.opts$use.avg.F)),
    capture.output(cat("  $use.FXSPR = ",proj.opts$use.FXSPR)),
    capture.output(cat("  $use.FMSY = ",proj.opts$use.FMSY)),
    capture.output(cat("  $proj.F = ",proj.opts$proj.F)),
    capture.output(cat("  $proj.catch = ",proj.opts$proj.catch)),"",sep='\n'))

  # add new data objects for projections
  data$do_proj = 1
  data$n_years_model <- data$n_years_model
  data$n_years_proj = proj.opts$n.yrs
  years_full <- head(input$years, length(input$years))
  years_full <- c(years_full, tail(years_full,1) + 1:proj.opts$n.yrs)
  input$years_full <- years_full
  #input$years_full = c(input$years, tail(input$years,1) +1:proj.opts$n.yrs)
  proj_yrs_ind <- data$n_years_model + 1:data$n_years_proj
  avg.yrs.ind <- match(proj.opts$avg.yrs, input$years)
  data$avg_years_ind = avg.yrs.ind - 1 # c++ indices start at 0
  
  #new option for making long term projections consistent with prevailing spr-based reference points.
  if(!is.null(proj.opts$proj_R_opt)){
    if(length(proj.opts$proj_R_opt)!= 1 | !(proj.opts$proj_R_opt %in% 1:2)) stop("proj.opts$proj_R_opt must be either 1 or 2.\n")
    data$proj_R_opt = proj.opts$proj_R_opt
  } else data$proj_R_opt = 1 #continue using mean of RE process for predicted R by default

  data$proj_F_opt = rep(0,data$n_years_proj) 
  if(!is.null(proj.opts$use.last.F)) if(proj.opts$use.last.F) data$proj_F_opt[] = 1
  if(!is.null(proj.opts$use.avg.F)) if(proj.opts$use.avg.F) data$proj_F_opt[] = 2
  if(!is.null(proj.opts$use.FXSPR)) if(proj.opts$use.FXSPR) data$proj_F_opt[] = 3
  if(!is.null(proj.opts$use.FMSY)) if(proj.opts$use.FMSY){
    if(data$recruit_model > 2) data$proj_F_opt[] = 6
    else{
      warning("Trying to use FMSY in projections but there is no stock-recruit model assumed. Will just use FXSPR in this case")
      data$proj_F_opt[] = 3
    }
  }
  data$proj_Fcatch = cbind(rep(0,data$n_years_proj))
  if(!is.null(proj.opts$proj.F)){
    data$proj_F_opt[] = 4
    if(length(proj.opts$proj.F) != data$n_years_proj) stop("length of proj.F is not = number of projection years")
    data$proj_Fcatch[] = proj.opts$proj.F
  }
  if(!is.null(proj.opts$proj.catch)){
    data$proj_F_opt[] = 5
    if(length(proj.opts$proj.catch) != data$n_years_proj) stop("length of proj.catch is not = number of projection years")
    data$proj_Fcatch[] = proj.opts$proj.catch
  }
  if(!is.null(proj.opts$proj_F_opt)){
    if(length(proj.opts$proj_F_opt) != data$n_years_proj) stop("length of proj_F_opt is not = number of projection years")
    data$proj_F_opt = proj.opts$proj_F_opt
  }
  if(!is.null(proj.opts$proj_Fcatch)){
    if(!is.matrix(proj.opts$proj_Fcatch)){
      if(length(proj.opts$proj_Fcatch) != data$n_years_proj) stop("length of proj_Fcatch is not = number of projection years")
      data$proj_Fcatch[] = proj.opts$proj_Fcatch
    } else{
      if(dim(proj.opts$proj_Fcatch)[1] != data$n_years_proj) stop("number of rows for proj_Fcatch is not = number of projection years")
      if(dim(proj.opts$proj_Fcatch)[2] != data$n_fleets) stop("number of cols for proj_Fcatch is not = number of fleets")
      data$proj_Fcatch = proj.opts$proj_Fcatch
    }
  }
  data$proj_Fcatch[which(!data$proj_F_opt %in% 4:5),] = 0

  if(any(data$proj_F_opt == 3)) data$percentFXSPR = proj.opts$percentFXSPR
  if(any(data$proj_F_opt == 6)) data$percentFMSY = proj.opts$percentFMSY
  
  data$FXSPR_init = c(data$FXSPR_init,rep(data$FXSPR_init[data$n_years_model], data$n_years_proj))
  data$FMSY_init = c(data$FMSY_init,rep(data$FMSY_init[data$n_years_model], data$n_years_proj))
  data$F_proj_init = rep(0.1, data$n_years_proj)
  data$F_proj_init[which(data$proj_F_opt == 3)] = data$FXSPR_init[proj_yrs_ind][which(data$proj_F_opt == 3)]
  data$F_proj_init[which(data$proj_F_opt == 6)] = data$FMSY_init[proj_yrs_ind][which(data$proj_F_opt == 6)]
  #define age for full F in projections
  FAA_proj <- colMeans(apply(model$rep$FAA[,avg.yrs.ind,,drop = FALSE],2:3, sum))
  #FAA_proj = colMeans(rbind(model$rep$FAA_tot[avg.yrs.ind,]))
  data$which_F_age = c(data$which_F_age, rep(which.max(FAA_proj), data$n_years_proj))

  # modify data objects for projections (pad with average over avg.yrs): mature, fracyr_SSB, waa
  avg_cols = function(x) apply(x, 2, mean, na.rm=TRUE)
  if(!is.null(proj.opts$proj_mature)){
    dims.check <- c(data$n_stocks, proj.opts$n.yrs, dim(data$mature)[2])
    cat("\nUsing user-suplied maturity values for projected maturity.\n")
    if(length(dim(proj.opts$proj_mature)) != length(dims.check)) {
      stop(paste0("\n** Error setting up projections: **\n",
                 "proj.opts$proj_mature must be an array with dimensions: ", paste(dims.check,collapse = ','), ".\n"))      
    }
    if(any(dim(proj.opts$proj_mature)!= dims.check)){
      stop(paste0("\n** Error setting up projections: **\n",
                 "proj.opts$proj_mature must be an array with dimensions: ", paste(dims.check,collapse = ','), ".\n"))      
    }
    data$mature_proj <- proj.opts$proj_mature
  } else {
    data$mature_proj <- array(0, dim = c(1,1,1)) #tests length on c++ side for whether to use it.
  }

  # proj_waa dims are (dim(input$data$waa)[1] x n.yrs x n_age)
  if(!is.null(proj.opts$proj_waa)){
    dims.check <- c(dim(data$waa)[1],proj.opts$n.yrs, dim(data$waa)[3])
    cat("\nUsing user-suplied WAA values for projected WAA.\n")
    if(length(dim(proj.opts$proj_waa)) != length(dims.check)) {
      stop(paste0("\n** Error setting up projections: **\n",
                 "proj.opts$proj_waa must be a 3d-array with dimensions: ", paste(dims.check,collapse = ','), ".\n"))      
    }
    if(any(dim(proj.opts$proj_waa)!= dims.check)){
      stop(paste0("\n** Error setting up projections: **\n",
                 "proj.opts$proj_waa must be a 3d-array with dimensions: ", paste(dims.check,collapse = ','), ".\n"))      
    }
    data$waa_proj <- proj.opts$proj_waa
  } else {
    data$waa_proj <- array(0, dim = c(1,1,1))  #tests length on c++ side for whether to use it.
  }

  
  # initialize pars at previously estimated values
  par <- model$parList
  # fill_vals <- function(x){as.factor(rep(NA, length(x)))}
  map <- input$map
  random <- input$random
  # map <- lapply(par, fill_vals)

  # SCAA (fixed effect Rec devs): set up logR_proj to treat recruitment as random effects in projections
  # will need to add this to likelihood for all stocks whether SCAA is used or not, but it will only be used in projections for SCAA stocks?
  if(any(data$NAA_re_model == 0)){ # SCAA
    if(any(proj.opts$avg.rec.yrs %in% model$years == FALSE)) stop(paste("","** Error setting up projections: **",
      "proj.opts$avg.rec.yrs is not a subset of model years.","",sep='\n'))
    if(is.null(proj.opts$avg.rec.yrs)) proj.opts$avg.rec.yrs <- model$years#[-1] # option for which model years to use for mean and SD calc
    avg.rec.yrs <- match(proj.opts$avg.rec.yrs, model$years)#[-1])
    data$logR_mean <- sapply(1:data$n_stocks, function(x) mean(log(model$rep$NAA[x,data$spawn_regions[x],avg.rec.yrs,1])))
    data$logR_sd <- sapply(1:data$n_stocks, function(x) sd(log(model$rep$NAA[x,data$spawn_regions[x],avg.rec.yrs,1])))
    par$logR_proj <- matrix(rep(data$logR_mean, each = data$n_years_proj), data$n_years_proj, data$n_stocks)
    map$logR_proj <- factor(1:length(par$logR_proj)) # turn on estimation of logR_proj
    random = c(random, "logR_proj")
  }

  if(all(data$NAA_re_model>0) & !is.null(proj.opts$avg.rec.yrs)) stop(paste("","** Error setting up projections: **",
    "proj.opts$avg.rec.yrs should only be used for SCAA model projections.
    This model already treats recruitment deviations as random effects.","",sep='\n'))

  # # in case projecting off of a peel
  # data$use_agg_catch <- cbind(data$use_agg_catch[1:data$n_years_model,])
  # data$use_catch_paa <- cbind(data$use_catch_paa[1:data$n_years_model,])
  # data$use_indices <- cbind(data$use_indices[1:data$n_years_model,])
  # data$use_index_paa <- cbind(data$use_index_paa[1:data$n_years_model,])

  # expand NAA_re
  input_NAA <- input
  input_NAA$asap3 <- NULL
  input_NAA$data$n_years_model <- data$n_years_model + data$n_years_proj
  input_NAA <- set_NAA(input_NAA, input$options$NAA_re) #use same machinery to map NAA and now we have options saved
  map$log_NAA <- input_NAA$map$log_NAA
  dims <- dim(input_NAA$par$log_NAA)
  log_NAA_init <- 10
  if(!is.null(proj.opts$proj_NAA_init)){
    if(length(proj.opts$proj_NAA_init)>1) stop("proj_NAA_init should be a single value")
    if(proj.opts$proj_NAA_init<0) stop("proj_NAA_init should be >0")
    log_NAA_init <- log(proj.opts$proj_NAA_init)
  }
  tmp <- array(log_NAA_init, dim = dims) #a large number for projections is particularly important when user is specifying catch.
  #for(s in 1:data$n_stocks) for(r in 1:data$n_regions) tmp[s,r,,] <- par$log_NAA[s,r,data$n_years_model-1,1] #fill in with last recruitment
  tmp[,,1:(data$n_years_model-1),] <- par$log_NAA[,,1:(data$n_years_model-1),]
  par$log_NAA <- tmp

  # check ecov options are valid, all Ecov will be projected if there observations do not occur in projection years
  # projection options are for how they are used in effects on population and are now done on c++ side
  # default is to continue Ecov in projection years, but will not matter if effects on M/mu and M/mu is averaged.
  # Ecov_out_R/M/mu/q will have projection years as determined. Ecov_x will still project any Ecov

  #Ecov, if everything is null, continue ecov processes
  if(all(is.null(proj.opts$cont.ecov), is.null(proj.opts$use.last.ecov), is.null(proj.opts$avg.ecov.yrs), is.null(proj.opts$proj.ecov))){
    data$proj_Ecov_opt = rep(1, data$n_Ecov)
    proj.opts$cont.ecov=TRUE 
    proj.opts$use.last.ecov=FALSE
  }
  if(is.null(proj.opts$cont.ecov)) proj.opts$cont.ecov=FALSE #one of the other ecov options is not null

  if(any(input$data$Ecov_model > 0)){
    end_model <- tail(input$years_full,1) #now need to go to the end of projection years
    end_Ecov <- tail(input$years_Ecov, 1)
    if(end_Ecov < end_model){
      print("prepare_projection: Ecov last year is before model last projection year. Padding Ecov...")
      map$Ecov_obs_logsigma <- matrix(as.integer(map$Ecov_obs_logsigma), ncol = data$n_Ecov)
      data$Ecov_obs <- rbind(data$Ecov_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      par$Ecov_obs_logsigma <- rbind(par$Ecov_obs_logsigma, matrix(par$Ecov_obs_logsigma[NROW(par$Ecov_obs_logsigma),], nrow = end_model-end_Ecov, ncol = data$n_Ecov, byrow=T))
      map$Ecov_obs_logsigma <- rbind(map$Ecov_obs_logsigma, matrix(NA, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_use_obs <- rbind(data$Ecov_use_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      
      n_years_proj_Ecov <- end_model - end_Ecov
      data$Ecov_use_proj <- matrix(0, nrow = n_years_proj_Ecov, ncol=data$n_Ecov)
      #data$years_use_Ecov <- 1:(data$n_years_Ecov + n_years_proj_Ecov) - 1

      # pad Ecov_re for projections
      map$Ecov_re = matrix(as.integer(map$Ecov_re), data$n_years_Ecov, data$n_Ecov)
      par$Ecov_re = rbind(par$Ecov_re, matrix(0,n_years_proj_Ecov, data$n_Ecov))

      tmp.re <- matrix(NA, n_years_proj_Ecov, data$n_Ecov)
      if(!is.null(proj.opts$proj.ecov)){ #projection ecov values are supplied, so use fixed values instead of RE
        if(NCOL(proj.opts$proj.ecov) != data$n_Ecov) stop("number of columns of proj.opts$proj.ecov is not equal to n_Ecov")
        if(NROW(proj.opts$proj.ecov) != n_years_proj_Ecov) stop(paste0("number of rows of proj.opts$proj.ecov should be ", n_years_proj_Ecov))
        for(i in 1:data$n_Ecov) if(data$Ecov_model[i]>0) {
          data$Ecov_use_proj[,i] <- proj.opts$proj.ecov[,i]
          # doesn't matter. Ecov_use_proj is just used to fill out Ecov_out_R, Ecov_out_M, etc. in projection years
          # if(data$Ecov_model[i] == 1) data$Ecov_use_proj[,i] <- proj.opts$proj.ecov[,i]  # random walk
          # if(data$Ecov_model[i] == 2) data$Ecov_use_proj[,i] <- proj.opts$proj.ecov[,i] - par$Ecov_process_pars[1,i] # AR(1)
        }
      }
      for(i in 1:data$n_Ecov) if(data$Ecov_model[i]>0) {
        tmp.re[,i] = 1
      }
      if(sum(!is.na(tmp.re))) tmp.re[which(!is.na(tmp.re))] <- max(map$Ecov_re, na.rm = TRUE) + 1:sum(!is.na(tmp.re))
      map$Ecov_re <- factor(rbind(map$Ecov_re, tmp.re))
      
      input$years_Ecov <- c(input$years_Ecov, seq(end_Ecov+1, end_model))
      map$Ecov_obs_logsigma = factor(map$Ecov_obs_logsigma)
      end_Ecov <- end_model
      data$n_years_Ecov <- length(input$years_Ecov)

      Ecov.opt.ct <- sum(proj.opts$cont.ecov, proj.opts$use.last.ecov, !is.null(proj.opts$avg.ecov.yrs), !is.null(proj.opts$proj.ecov))
      if(Ecov.opt.ct == 0) proj.opts$cont.ecov = TRUE; Ecov.opt.ct = 1;
      if(Ecov.opt.ct != 1) stop(paste("","** Error setting up projections: **",
        "Exactly one method of specifying ecov must be used (see ?project_wham).",
        "You have specified these in 'proj.opts':",
        capture.output(cat("  $cont.ecov = ",proj.opts$cont.ecov)),
        capture.output(cat("  $use.last.ecov = ",proj.opts$use.last.ecov)),
        capture.output(cat("  $avg.ecov.yrs = ",proj.opts$avg.ecov.yrs)),
        capture.output(cat("  $proj.ecov = ",proj.opts$proj.ecov)),"",sep='\n'))
      if(!is.null(proj.opts$avg.ecov.yrs) & any(proj.opts$avg.ecov.yrs %in% model$years == FALSE)) stop(paste("","** Error setting up projections: **",
        "proj.opts$avg.ecov.yrs is not a subset of model years.","",sep='\n'))

      if(!is.null(proj.opts$avg.ecov.yrs)) {
        data$avg_years_Ecov <- match(proj.opts$avg.ecov.yrs, model$years) - 1
        data$proj_Ecov_opt <- rep(2, data$n_Ecov)
      }
      if(!is.null(proj.opts$proj.ecov)) data$proj_Ecov_opt <- rep(4, data$n_Ecov)
      if(proj.opts$use.last.ecov) data$proj_Ecov_opt <- rep(3, data$n_Ecov)
      if(proj.opts$cont.ecov) data$proj_Ecov_opt <- rep(1, data$n_Ecov)

    } else{
      print(paste0("Ecovs are already fit through projection years."))
    }
  }
    # n.beyond <- end.beyond <- integer()
    # for(i in 1:data$n_Ecov) {
    #   n.beyond[i] = data$n_years_Ecov-1-max(data$ind_Ecov_out_end_R[i,],data$ind_Ecov_out_end_mu[i,,,,,],
    #     data$ind_Ecov_out_end_M[i,,,],data$ind_Ecov_out_end_q[i,])
    #   end.beyond[i] <- min(n.beyond[i], data$n_years_proj)
    #   stop()
    #   if(end.beyond[i] == data$n_years_proj) print(paste0("ecov ",i," already fit through projection years. Using fit ecov ",i," for projections..."))
    # }

  #   if(all(end.beyond < data$n_years_proj)){ # if Ecov proj options ARE necessary, check they are valid
  #     Ecov.opt.ct <- sum(proj.opts$cont.ecov, proj.opts$use.last.ecov, !is.null(proj.opts$avg.ecov.yrs), !is.null(proj.opts$proj.ecov))
  #     if(Ecov.opt.ct == 0) proj.opts$cont.ecov = TRUE; Ecov.opt.ct = 1;
  #     if(Ecov.opt.ct != 1) stop(paste("","** Error setting up projections: **",
  #       "Exactly one method of specifying ecov must be used (see ?project_wham).",
  #       "You have specified these in 'proj.opts':",
  #       capture.output(cat("  $cont.ecov = ",proj.opts$cont.ecov)),
  #       capture.output(cat("  $use.last.ecov = ",proj.opts$use.last.ecov)),
  #       capture.output(cat("  $avg.ecov.yrs = ",proj.opts$avg.ecov.yrs)),
  #       capture.output(cat("  $proj.ecov = ",proj.opts$proj.ecov)),"",sep='\n'))
  #     if(!is.null(proj.opts$avg.ecov.yrs) & any(proj.opts$avg.ecov.yrs %in% model$years == FALSE)) stop(paste("","** Error setting up projections: **",
  #       "proj.opts$avg.ecov.yrs is not a subset of model years.","",sep='\n'))
  #   }
  # } else { # need to create objects if no Ecov
  #   end.beyond = rep(data$n_years_proj, data$n_Ecov) # effectively say that Ecov already extends # of proj years
  # }
  #data$n_years_proj_Ecov = max(data$n_years_proj-end.beyond)

  # pad Ecov_re for projections
  # Ecov.proj <- matrix(0, nrow = data$n_years_proj_Ecov, ncol=data$n_Ecov)
  # if(any(data$Ecov_model > 0)){
  #   if(data$n_years_proj_Ecov>0){ # need to pad Ecov_re
  #     par$Ecov_re <- rbind(par$Ecov_re[1:data$n_years_Ecov,,drop=F], Ecov.proj) # pad Ecov_re if necessary
  #     data$Ecov_use_re <- rbind(data$Ecov_use_re, matrix(0, nrow=data$n_years_proj_Ecov, ncol=data$n_Ecov))
  #     map$Ecov_re = matrix(as.integer(map$Ecov_re), data$n_years_Ecov, data$n_Ecov)
  #     tmp.re <- matrix(NA, data$n_years_proj_Ecov, data$n_Ecov)
  #     for(i in 1:data$n_Ecov) if(data$Ecov_model[i]>0) {
  #       tmp.re[,i] = 1
  #       data$Ecov_use_re[,i] = c(data$Ecov_use_re[1:data$n_years_Ecov,i], tmp.re[,i])
  #     }
  #     if(sum(!is.na(tmp.re))) tmp.re[which(!is.na(tmp.re))] <- max(map$Ecov_re, na.rm = TRUE) + 1:sum(!is.na(tmp.re)) 
  #     map$Ecov_re <- factor(rbind(map$Ecov_re, tmp.re))
  #   }
  # }
  
  #M
  # options for M in projections, data$proj_M_opt:
  #   1 = continue random effects (if they exist) - need to pad M_re
  #   2 = use average
  if(!is.null(proj.opts$cont.M.re)){
    if(proj.opts$cont.M.re & !"M_re" %in% input$random){
      stop(paste("","** Error setting up projections **",
        "proj.opts$cont.M.re = TRUE but no random effects on M"))
    }
    data$proj_M_opt <- ifelse(proj.opts$cont.M.re, 1, 2) # 1 = continue M_re, 2 = average
  } else { # if NULL, default is to continue M random effects (if they exist!)
    data$proj_M_opt <- 2 #default is to use average M
    #data$proj_M_opt <- ifelse(model$env$data$M_re_model %in% c(2,4,5), 1, 2) # 2 = IID, 4 = AR1_y, 5 = 2D AR1
  }
  # expand M_re
  input_M <- input
  input_M$asap3 <- NULL
  input_M$data$n_years_model <- data$n_years_model + data$n_years_proj
  input_M <- set_M(input_M, input$options$M) #use same machinery to map M and now we have options saved
  map$M_re <- input_M$map$M_re
  dims <- dim(par$M_re)
  dims[3] <- data$n_years_model + data$n_years_proj
  tmp <- array(0, dim = dims)
  tmp[,,1:data$n_years_model,] <- par$M_re[,,1:data$n_years_model,]
  #for(a in 1:dims[1]) for(b in 1:dims[2]) for(c in 1:dims[4]) tmp[a,b,,c] <- c(par$M_re[a,b,,c], rep(0,data$n_years_proj))
  par$M_re <- tmp
      # print("here6")

  # options for L in projections, data$proj_L_opt:
  #   1 = continue random effects (if they exist) - need to pad L_re
  #   2 = use average
  if(!is.null(proj.opts$cont.L.re)){
    if(proj.opts$cont.L.re & !"L_re" %in% input$random){
      stop(paste("","** Error setting up projections **",
        "proj.opts$cont.L.re = TRUE but no random effects on L"))
    }
    data$proj_L_opt <- ifelse(proj.opts$cont.L.re, 1, 2) # 1 = continue L_re, 2 = average
  } else { # if NULL, default is to continue L random effects (if they exist!)
    data$proj_L_opt <- 2 #default is to use average L
  }
  #expand L_re
  input_L <- input
  input_L$asap3 <- NULL
  input_L$data$n_years_model <- data$n_years_model + data$n_years_proj
  input_L <- set_L(input_L, input$options$L) #use same machinery to map M and now we have options saved
  map$L_re <- input_L$map$L_re
  dims <- dim(par$L_re)
  dims[1] <- data$n_years_model + data$n_years_proj
  tmp <- matrix(0, dims[1], dims[2])
  tmp[1:data$n_years_model,] <- par$L_re[1:data$n_years_model,]
  par$L_re <- tmp

  # options for mu in projections, data$proj_mu_opt:
  #   1 = continue random effects (if they exist) - need to pad mu_re
  #   2 = use average
  if(!is.null(proj.opts$cont.move.re)){
    if(data$n_regions == 1) stop("proj.opts$cont.move.re is specified but there is only one region.")
    if(proj.opts$cont.move..re & !"mu_re" %in% input$random){
      stop(paste("","** Error setting up projections **",
        "proj.opts$cont.move.re = TRUE but no random effects on movement"))
    }
    data$proj_mu_opt <- ifelse(proj.opts$cont.move.re, 1, 2) # 1 = continue mu_re, 2 = average
  } else { # if NULL, default is to continue mu random effects (if they exist!)
    data$proj_mu_opt <- 2 #default is to use average mu
  }
  # expand mu_re
  #if(data$n_regions>1){
    input_mu <- input
    input_mu$asap3 <- NULL
    input_mu$data$n_years_model <- data$n_years_model + data$n_years_proj
    input_mu <- set_move(input_mu, input$options$move) #use same machinery to map M and now we have options saved
    map$mu_re <- input_mu$map$mu_re
    dims <- dim(par$mu_re)
    dims[4] <- data$n_years_model + data$n_years_proj
    tmp <- array(0, dim = dims)
    tmp[,,,1:data$n_years_model,,] <- par$mu_re[,,,1:data$n_years_model,,]
    par$mu_re <- tmp
  #}

  # expand q_re
  input_q <- input
  input_q$asap3 <- NULL
  input_q$data$n_years_model <- data$n_years_model + data$n_years_proj
  input_q <- set_q(input_q, input$options$q) #use same machinery to map q and now we have options saved
  map$q_re <- input_q$map$q_re
  dims <- dim(par$q_re)
  dims[2] <- data$n_years_model + data$n_years_proj
  tmp <- matrix(0, nrow = data$n_years_model + data$n_years_proj, ncol = NCOL(par$q_re))
  tmp[1:data$n_years_model,] <- par$q_re[1:data$n_years_model,]
  par$q_re <- tmp

  data$years_use <- 1:(data$n_years_model + data$n_years_proj) -1
  data$years_use_Ecov <-  1:data$n_years_Ecov - 1
  input$data <- data
  input$par <- par
  input$map <- map
  input$random <- random
  input$options$proj <- proj.opts
  attr(input$par, 'check.passed') = NULL
  attr(input$data, 'check.passed') = NULL
  return(input)
}
