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
#'     \item \code{$use.FXSPR} (T/F), calculate and use F at X% SPR for projections. Default = \code{FALSE}.
#'     \item \code{$use.FMSY} (T/F), calculate and use FMSY for projections. Default = \code{FALSE}.
#'     \item \code{$proj.F} (vector), user-specified fishing mortality for projections. Length must equal \code{n.yrs}.
#'     \item \code{$proj.catch} (vector), user-specified aggregate catch for projections. Length must equal \code{n.yrs}.
#'     \item \code{$avg.yrs} (vector), specify which years to average over for calculating reference points. Default = last 5 model years, \code{tail(model$years, 5)}.
#'     \item \code{$cont.ecov} (T/F), continue ecov process (e.g. random walk or AR1) for projections. Default = \code{TRUE}.
#'     \item \code{$use.last.ecov} (T/F), use terminal year ecov for projections.
#'     \item \code{$avg.ecov.yrs} (vector), specify which years to average over the environmental covariate(s) for projections.
#'     \item \code{$proj.ecov} (matrix), user-specified environmental covariate(s) for projections. \code{n.yrs x n.ecov}.
#'     \item \code{$cont.Mre} (T/F), continue M random effects (i.e. AR1_y or 2D AR1) for projections. Default = \code{TRUE}. If \code{FALSE}, M will be averaged over \code{$avg.yrs} (which defaults to last 5 model years).
#'     \item \code{$avg.rec.yrs} (vector), specify which years to calculate the CDF of recruitment for use in projections. Default = all model years. Only used when recruitment is estimated as fixed effects (SCAA).
#'     \item \code{$percentFXSPR} (scalar), percent of F_XSPR to use for calculating catch in projections, only used if $use.FXSPR = TRUE. For example, GOM cod uses F = 75% F_40%SPR, so \code{proj.opts$percentFXSPR = 75}. Default = 100.
#'     \item \code{$percentFMSY} (scalar), percent of F_MSY to use for calculating catch in projections, only used if $use.FMSY = TRUE.
#'     \item \code{$proj_F_opt} (vector), integers specifying how to configure each year of the projection: 1: use terminal F, 2: use average F, 3: use F at X% SPR, 4: use specified F, 5: use specified catch, 6: use Fmsy. Overrides any of the above specifications.
#'     \item \code{$proj_Fcatch} (vector), catch or F values to use each projection year: values are not used when using Fmsy, FXSPR, terminal F or average F. Overrides any of the above specifications of proj.F or proj.catch.
#'   }
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
prepare_projection = function(model, proj.opts)
{
  if(is.null(proj.opts)) proj.opts=list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE, use.FXSPR=FALSE, use.FMSY=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
                                       cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL, cont.Mre=NULL, avg.rec.yrs=NULL, percentFXSPR=100,
                                       percentFMSY=100, proj_F_opt = NULL, proj_Fcatch = NULL)
  # default: 3 projection years
  if(is.null(proj.opts$n.yrs)) proj.opts$n.yrs <- 3
  # default: use average M, selectivity, etc. over last 5 model years to calculate ref points
  if(is.null(proj.opts$avg.yrs)) proj.opts$avg.yrs <- tail(model$years, 5)  
  if(all(is.null(proj.opts$use.last.F), is.null(proj.opts$use.avg.F), is.null(proj.opts$use.FXSPR), is.null(proj.opts$use.FMSY), is.null(proj.opts$proj.F), is.null(proj.opts$proj.catch))){
    proj.opts$use.last.F=TRUE; proj.opts$use.avg.F=FALSE; proj.opts$use.FXSPR=FALSE; proj.opts$use.FMSY=FALSE; proj.opts$proj.F=NULL; proj.opts$proj.catch=NULL
  }
  if(all(is.null(proj.opts$cont.ecov), is.null(proj.opts$use.last.ecov), is.null(proj.opts$avg.ecov.yrs), is.null(proj.opts$proj.ecov))){
    proj.opts$cont.ecov=TRUE; proj.opts$use.last.ecov=FALSE; proj.opts$avg.ecov.yrs=NULL; proj.opts$proj.ecov=NULL
  }
  if(is.null(proj.opts$cont.ecov)) proj.opts$cont.ecov=FALSE
  if(is.null(proj.opts$percentFXSPR)) proj.opts$percentFXSPR = 100
  if(is.null(proj.opts$percentFMSY)) proj.opts$percentFMSY = 100

  input1 <- model$input
  data <- input1$data
  # options for M in projections, data$proj_M_opt:
  #   1 = continue random effects (if they exist) - need to pad M_re
  #   2 = use average
  if(!is.null(proj.opts$cont.Mre)){
    if(proj.opts$cont.Mre & model$env$data$M_re_model %in% c(1,3)){ # 1 = none, 3 = AR1_a (not time-varying)
      stop(paste("","** Error setting up projections **",
        "proj.opts$cont.Mre = TRUE but no time-varying random effects on M"))
    }
    data$proj_M_opt <- ifelse(proj.opts$cont.Mre, 1, 2) # 1 = continue M_re, 2 = average
  } else { # if NULL, default is to continue M random effects (if they exist!)
    data$proj_M_opt <- ifelse(model$env$data$M_re_model %in% c(2,4,5), 1, 2) # 2 = IID, 4 = AR1_y, 5 = 2D AR1
  }

  # check options for F/catch are valid
  if(any(proj.opts$avg.yrs %in% model$years == FALSE)) stop(paste("","** Error setting up projections: **",
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

  # check ecov options are valid
  if(any(input1$data$Ecov_model > 0)){
    n_effects = dim(input1$par$Ecov_beta)[1]
    # if Ecovs already extend beyond population model, Ecov proj options are unnecessary
    n.beyond <- end.beyond <- integer()
    for(i in 1:data$n_Ecov) {
      n.beyond[i] = data$n_years_Ecov-1-max(data$ind_Ecov_out_end[i,])
      end.beyond[i] <- min(n.beyond[i], proj.opts$n.yrs)     
      if(end.beyond[i] == proj.opts$n.yrs) print(paste0("ecov ",i," already fit through projection years. Using fit ecov ",i," for projections..."))
    }
    Ecov.proj <- matrix(rep(NA, max(proj.opts$n.yrs-end.beyond)*data$n_Ecov), ncol=data$n_Ecov)
    #this isn't used anywhere
    #Ecov.map <- matrix(rep(NA, max(proj.opts$n.yrs-end.beyond)*data$n_Ecov), ncol=data$n_Ecov)
    if(all(end.beyond < proj.opts$n.yrs)){ # if Ecov proj options ARE necessary, check they are valid
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
    }
  } else { # need to create objects if no Ecov
    end.beyond = rep(proj.opts$n.yrs, data$n_Ecov) # effectively say that Ecov already extends # of proj years
  }

  # add new data objects for projections
  data$do_proj = 1
  data$n_years_proj = proj.opts$n.yrs
  data$n_years_proj_Ecov = max(proj.opts$n.yrs-end.beyond)
  avg.yrs.ind <- match(proj.opts$avg.yrs, input1$years)
  data$avg_years_ind = avg.yrs.ind - 1 # c++ indices start at 0
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
  data$proj_Fcatch = rep(0,data$n_years_proj)
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
    if(length(proj.opts$proj_Fcatch) != data$n_years_proj) stop("length of proj_Fcatch is not = number of projection years")
    data$proj_Fcatch = proj.opts$proj_Fcatch
  }
  data$proj_Fcatch[which(!data$proj_F_opt %in% 4:5)] = 0

  if(any(data$proj_F_opt == 3)) data$percentFXSPR = proj.opts$percentFXSPR
  if(any(data$proj_F_opt == 6)) data$percentFMSY = proj.opts$percentFMSY
#     else {
# stop(paste("","** Error setting up projections: **",
#         "percentFXSPR is not used if FXSPR is not used to calculate catch in projections.",
#         "You have chosen ",
#         capture.output(cat(c("use.last.F","use.avg.F","use.FXSPR","proj.F","proj.catch")[data$proj_F_opt])),"",sep='\n'))      
#     }

  data$FXSPR_init = c(data$FXSPR_init,rep(data$FXSPR_init[data$n_years_model], data$n_years_proj))
  data$FMSY_init = c(data$FMSY_init,rep(data$FMSY_init[data$n_years_model], data$n_years_proj))
  data$F_proj_init = rep(0.1, data$n_years_proj)
  data$F_proj_init[which(data$proj_F_opt == 3)] = data$FXSPR_init[which(data$proj_F_opt == 3)]
  data$F_proj_init[which(data$proj_F_opt == 6)] = data$FMSY_init[which(data$proj_F_opt == 6)]
  #define age for full F in projections
  FAA_proj = colMeans(rbind(model$rep$FAA_tot[avg.yrs.ind,]))
  data$which_F_age = c(data$which_F_age, rep(which(FAA_proj == max(FAA_proj))[1], data$n_years_proj))

  # modify data objects for projections (pad with average over avg.yrs): mature, fracyr_SSB, waa
  avg_cols = function(x) apply(x, 2, mean, na.rm=TRUE)
  data$mature <- rbind(data$mature[1:data$n_years_model,], matrix(rep(avg_cols(data$mature[avg.yrs.ind,]), proj.opts$n.yrs), nrow=proj.opts$n.yrs, byrow=TRUE))
  data$fracyr_SSB <- c(data$fracyr_SSB[1:data$n_years_model], rep(mean(data$fracyr_SSB[avg.yrs.ind]), proj.opts$n.yrs))
  toadd <- apply(data$waa[,avg.yrs.ind,], c(1,3), mean)
  if(dim(data$waa)[2] > data$n_years_model) data$waa <- data$waa[,1:data$n_years_model,]
  tmp <- array(NA, dim = dim(data$waa) + c(0,proj.opts$n.yrs,0))
  tmp[,1:data$n_years_model,] <- data$waa
  for(i in seq(data$n_years_model+1,data$n_years_model+proj.opts$n.yrs)) tmp[,i,] <- toadd
  data$waa <- tmp

  # initialize pars at previously estimated values
  par <- model$parList
  # fill_vals <- function(x){as.factor(rep(NA, length(x)))}
  map <- input1$map
  random <- input1$random
  # map <- lapply(par, fill_vals)

  # SCAA (fixed effect Rec devs): set up logR_proj to treat recruitment as random effects in projections
  if(data$n_NAA_sigma == 0){ # SCAA
    if(any(proj.opts$avg.rec.yrs %in% model$years == FALSE)) stop(paste("","** Error setting up projections: **",
      "proj.opts$avg.rec.yrs is not a subset of model years.","",sep='\n'))
    if(is.null(proj.opts$avg.rec.yrs)) proj.opts$avg.rec.yrs <- model$years[-1] # option for which model years to use for mean and SD calc
    avg.rec.yrs <- match(proj.opts$avg.rec.yrs, model$years[-1])
    data$logR_mean <- mean(par$log_NAA[avg.rec.yrs,1])
    data$logR_sd <- sd(par$log_NAA[avg.rec.yrs,1])
    par$logR_proj <- rep(data$logR_mean, proj.opts$n.yrs)
    map$logR_proj <- factor(1:proj.opts$n.yrs) # turn on estimation of logR_proj
    random = c(random, "logR_proj")

    # pad log_NAA (even though projection years not used)
    # par$log_NAA <- rbind(par$log_NAA, matrix(10, nrow=proj.opts$n.yrs, ncol=data$n_ages, byrow=TRUE))
    # map$log_NAA <- factor(c(map$log_NAA, rep(NA, proj.opts$n.yrs*data$n_ages)))
    par$log_NAA <- rbind(par$log_NAA[1:(data$n_years_model-1),], matrix(10, nrow=proj.opts$n.yrs, ncol=data$n_ages, byrow=TRUE))
    map$log_NAA <- factor(rbind(matrix(as.numeric(map$log_NAA), data$n_years_model-1, data$n_ages)[1:(data$n_years_model-1),], matrix(NA, proj.opts$n.yrs, data$n_ages)))
  } else {
    # recruit RE: pad log_NAA, map ages > 1 to NA
    # all NAA RE: pad log_NAA, estimate all
    if(!is.null(proj.opts$avg.rec.yrs)) stop(paste("","** Error setting up projections: **",
      "proj.opts$avg.rec.yrs should only be used for SCAA model projections.
      This model already treats recruitment deviations as random effects.","",sep='\n'))

    # pad log_NAA
    log_NAA_mean <- apply(par$log_NAA[1:(data$n_years_model-1),], 2, mean) # start at mean for each age
    par$log_NAA <- rbind(par$log_NAA[1:(data$n_years_model-1),], matrix(log_NAA_mean, nrow=proj.opts$n.yrs, ncol=data$n_ages, byrow=TRUE))
    tmp <- par$log_NAA
    if(data$n_NAA_sigma == 1) tmp[,-1] <- NA # don't estimate NAA_devs for ages > 1 if RE only on recruitment
    ind.notNA <- which(!is.na(tmp))
    tmp[ind.notNA] <- 1:length(ind.notNA)
    map$log_NAA = factor(tmp)
  }

    # pad q_re
    par$q_re <- rbind(par$q_re[1:data$n_years_model,,drop=F], matrix(0, data$n_years_proj, data$n_indices))
    map$q_re <- rbind(matrix(as.integer(map$q_re), data$n_years_model, data$n_indices), matrix(NA, data$n_years_proj, data$n_indices))
    if(any(data$use_q_re == 1)){
      ind = which(data$use_q_re)
      map$q_re[,ind] <- 1:(length(ind) * (data$n_years_model + data$n_years_proj)) #turn on appropriate columns of q_re
    }
    map$q_re <- factor(map$q_re)


  # pad Ecov_re for projections
  if(any(data$Ecov_model > 0)){
    if(any(end.beyond < proj.opts$n.yrs)){ # need to pad Ecov_re
      for(j in 1:data$n_Ecov){
        # for(i in 1:(proj.opts$n.yrs-end.beyond[j])){
        for(i in 1:max(proj.opts$n.yrs-end.beyond)){
          if(!is.null(proj.opts$use.last.ecov)) if(proj.opts$use.last.ecov){ # use last Ecov (pad Ecov_re but map to NA)
            Ecov.proj[i,j] <- model$rep$Ecov_re[max(data$ind_Ecov_out_end[j,])+1+end.beyond[j],j]
          }
          if(!is.null(proj.opts$avg.ecov.yrs)){ # use average Ecov (pad Ecov_re but map to NA)
            ecov.yrs <- data$year1_Ecov+0:(data$n_years_Ecov-1+data$n_years_proj_Ecov)
            avg.yrs.ind.ecov <- match(proj.opts$avg.ecov.yrs, ecov.yrs)
            Ecov.proj[i,j] <- avg_cols(as.matrix(model$rep$Ecov_re[avg.yrs.ind.ecov,j]))
          }
          if(!is.null(proj.opts$cont.ecov)) if(proj.opts$cont.ecov){ # continue Ecov process (pad Ecov_re and estimate)
            Ecov.proj[i,j] <- 0
          }
          if(!is.null(proj.opts$proj.ecov)){ # use specified Ecov, may have to back-calculate Ecov_re from Ecov_x
              if(data$Ecov_model[j] == 1) Ecov.proj[i,j] <- proj.opts$proj.ecov[i,j] # random walk
              if(data$Ecov_model[j] == 2) Ecov.proj[i,j] <- proj.opts$proj.ecov[i,j] - par$Ecov_process_pars[1,j] # AR(1)
          }
        }
      }
      par$Ecov_re <- rbind(par$Ecov_re[1:data$n_years_Ecov,,drop=F], Ecov.proj) # pad Ecov_re if necessary

      # pad map$Ecov_re
      tmp.re <- matrix(1:length(par$Ecov_re), dim(par$Ecov_re)[1], data$n_Ecov, byrow=FALSE)
      data$Ecov_use_re <- matrix(0, nrow=data$n_years_Ecov + data$n_years_proj_Ecov, ncol=data$n_Ecov)
      for(i in 1:data$n_Ecov){
        tmp.re[,i] <- if(data$Ecov_model[i]==0) rep(NA,dim(par$Ecov_re)[1]) else tmp.re[,i]
        if(data$Ecov_model[i]==1) tmp.re[1,i] <- NA # if Ecov is a rw, first year of Ecov_re is not used bc Ecov_x[1] uses Ecov1 (fixed effect)
      }
      if(!proj.opts$cont.ecov){
        for(i in 1:data$n_Ecov){
          tmp.re[1:(proj.opts$n.yrs-end.beyond[i])+data$n_years_Ecov,] <- NA
        }
      } 
      ind.notNA <- which(!is.na(tmp.re))
      tmp.re[ind.notNA] <- 1:length(ind.notNA)
      data$Ecov_use_re[ind.notNA] <- 1 # don't want to add to NLL in projection years (= 0)
      map$Ecov_re = factor(tmp.re)
    }
  }

  # only need to pad M_re if continuing M_re in projections (otherwise will average)
  # only M_re models 2, 4, 5
  if(data$proj_M_opt == 1){ 
    # par$M_re <- rbind(par$M_re, matrix(NA, nrow=proj.opts$n.yrs, ncol=data$n_ages))
    # tmp <- par$M_re
    # ind.proj <- which(is.na(tmp))
    # par$M_re[ind.proj] <- 0

    par$M_re <- rbind(par$M_re, matrix(0, nrow=proj.opts$n.yrs, ncol=data$n_ages))
    tmp <- par$M_re
    if(data$M_re_model %in% c(2,5)){ # iid / 2d ar1
      tmp[] = 1:(dim(tmp)[1]*dim(tmp)[2]) # all y,a estimated
    }
    if(data$M_re_model == 4){ # ar1_y (devs by year, constant by age)
      for(i in 1:dim(tmp)[1]) tmp[i,] = i
    }
    map$M_re <- factor(tmp)

    # tmp[-ind.proj] <- NA
    # if(data$M_re_model %in% c(2,5)) tmp[ind.proj] <- 1:length(ind.proj)
    # if(data$M_re_model == 4){ # ar1_y (devs by year, constant by age)
    #   for(i in 1:proj.opts$n.yrs) tmp[data$n_years_model+i,] = i
    # }
    # map$M_re <- factor(tmp)
  }

  input2 <- list(data=data, par = par, map = map, random = random,
    years = input1$years, years_full = c(input1$years, tail(input1$years,proj.opts$n.yrs) + proj.opts$n.yrs),
    ages.lab = input1$ages.lab, model_name = input1$model_name)
  attr(input2$par, 'check.passed') = NULL
  attr(input2$data, 'check.passed') = NULL
  return(input2)
}
