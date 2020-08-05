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
#'     \item \code{$use.FXSPR} (T/F), calculate F at X% SPR for projections.
#'     \item \code{$proj.F} (vector), user-specified fishing mortality for projections. Length must equal \code{n.yrs}.
#'     \item \code{$proj.catch} (vector), user-specified aggregate catch for projections. Length must equal \code{n.yrs}.
#'     \item \code{$avg.yrs} (vector), specify which years to average over for calculating reference points. Default = last 5 model years, \code{tail(model$years, 5)}.
#'     \item \code{$cont.ecov} (T/F), continue ecov process (e.g. random walk or AR1) for projections. Default = \code{TRUE}.
#'     \item \code{$use.last.ecov} (T/F), use terminal year ecov for projections.
#'     \item \code{$avg.ecov.yrs} (vector), specify which years to average over the environmental covariate(s) for projections.
#'     \item \code{$proj.ecov} (vector), user-specified environmental covariate(s) for projections. Length must equal \code{n.yrs}.
#'     \item \code{$cont.Mre} (T/F), continue M random effects (i.e. AR1_y or 2D AR1) for projections. Default = \code{TRUE}. If \code{FALSE}, M will be averaged over \code{$avg.yrs} (which defaults to last 5 model years).
#'   }
#'
#' @return same as \code{\link{prepare_wham_input}}, a list ready for \code{\link{fit_wham}}:
#'   \describe{
#'     \item{\code{data}}{Named list of data, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{par}}{Named list of parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{map}}{not sure what this does, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{random}}{Character vector of parameters to treat as random effects, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{years}}{Numeric vector of years to fit WHAM model (specified in ASAP3 .dat file)}
#'     \item{\code{ages.lab}}{Character vector of age labels, ending with plus-group (specified in ASAP3 .dat file)}
#'     \item{\code{model_name}}{Character, name of stock/model (specified in call to \code{prepare_wham_input})}
#'   }
#'
#' @seealso \code{\link{prepare_wham_input}}, \code{\link{project_wham}}
#'
prepare_projection = function(model, proj.opts)
{
# library(wham)
# write.dir <- "/home/bstock/Documents/wham/sandbox/ex3_projections"
# setwd(write.dir)
# model = readRDS("m6.rds")
# proj.opts=list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE, use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL,
#                 avg.yrs=NULL, cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL)
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

  # default: use average M, selectivity, etc. over last 5 model years to calculate ref points
  if(is.null(proj.opts$avg.yrs)) proj.opts$avg.yrs <- tail(model$years, 5)

  # default: 3 projection years
  if(is.null(proj.opts$n.yrs)) proj.opts$n.yrs <- 3

  # check options for F/catch are valid
  if(any(proj.opts$avg.yrs %in% model$years == FALSE)) stop(paste("","** Error setting up projections: **",
    "proj.opts$avg.yrs is not a subset of model years.","",sep='\n'))
  F.opt.ct <- sum(proj.opts$use.avg.F, proj.opts$use.last.F, proj.opts$use.FXSPR, !is.null(proj.opts$proj.F), !is.null(proj.opts$proj.catch))
  if(F.opt.ct != 1) stop(paste("","** Error setting up projections: **",
    "Exactly one method of specifying F must be used (see ?project_wham).",
    "You have specified these in 'proj.opts':",
    capture.output(cat("  $use.last.F = ",proj.opts$use.last.F)),
    capture.output(cat("  $use.avg.F = ",proj.opts$use.avg.F)),
    capture.output(cat("  $use.FXSPR = ",proj.opts$use.FXSPR)),
    capture.output(cat("  $proj.F = ",proj.opts$proj.F)),
    capture.output(cat("  $proj.catch = ",proj.opts$proj.catch)),"",sep='\n'))

  # check ecov options are valid
  if(any(input1$data$Ecov_model > 0)){
    # if Ecov already extends beyond population model, Ecov proj options are unnecessary
    n.beyond <- data$n_years_Ecov-1-data$ind_Ecov_out_end
    end.beyond <- min(n.beyond, proj.opts$n.yrs)
    Ecov.proj <- matrix(rep(NA, (proj.opts$n.yrs-end.beyond)*dim(model$rep$Ecov_re)[2]), ncol=dim(model$rep$Ecov_re)[2], byrow=TRUE)
    Ecov.map <- matrix(rep(NA, (proj.opts$n.yrs-end.beyond)*dim(model$rep$Ecov_re)[2]), ncol=dim(model$rep$Ecov_re)[2], byrow=TRUE)
    if(end.beyond == proj.opts$n.yrs){ print("ecov already fit through projection years. Using fit ecov for projections...")
    } else { # if Ecov proj options ARE necessary, check they are valid
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
    end.beyond = proj.opts$n.yrs # effectively say that Ecov already extends # of proj years
  }

  # add new data objects for projections
  data$do_proj = 1
  data$n_years_proj = proj.opts$n.yrs
  data$n_years_proj_Ecov = proj.opts$n.yrs-end.beyond
  avg.yrs.ind <- match(proj.opts$avg.yrs, input1$years)
  data$avg_years_ind = avg.yrs.ind - 1 # c++ indices start at 0
  data$proj_F_opt = rep(0,data$n_years_proj) 
  if(proj.opts$use.last.F) data$proj_F_opt[] = 1
  if(proj.opts$use.avg.F) data$proj_F_opt[] = 2
  if(proj.opts$use.FXSPR) data$proj_F_opt[] = 3
  if(!is.null(proj.opts$proj.F)){
    data$proj_F_opt[] = 4
    data$proj_Fcatch = proj.opts$proj.F
  }
  if(!is.null(proj.opts$proj.catch)){
    data$proj_F_opt[] = 5
    data$proj_Fcatch = proj.opts$proj.catch
  }
  if(data$proj_F_opt[1] %in% 1:3) data$proj_Fcatch = rep(0, proj.opts$n.yrs)

  #define age for full F in projections
  FAA_proj = colMeans(rbind(model$rep$FAA_tot[avg.yrs.ind,]))
  data$which_F_age = which(FAA_proj == max(FAA_proj))[1]

  # modify data objects for projections (pad with average over avg.yrs): mature, fracyr_SSB, waa
  avg_cols = function(x) apply(x, 2, mean, na.rm=TRUE)
  data$mature <- rbind(data$mature, matrix(rep(avg_cols(data$mature[avg.yrs.ind,]), proj.opts$n.yrs), nrow=proj.opts$n.yrs, byrow=TRUE))
  data$fracyr_SSB <- c(data$fracyr_SSB, rep(mean(data$fracyr_SSB[avg.yrs.ind]), proj.opts$n.yrs))
  toadd <- apply(data$waa[,avg.yrs.ind,], c(1,3), mean)
  tmp <- array(NA, dim = dim(data$waa) + c(0,proj.opts$n.yrs,0))
  tmp[,1:data$n_years_model,] <- data$waa
  for(i in seq(data$n_years_model+1,data$n_years_model+proj.opts$n.yrs)) tmp[,i,] <- toadd
  data$waa <- tmp

  # initialize pars at previously estimated values
  par <- model$parList
  # fill_vals <- function(x){as.factor(rep(NA, length(x)))}
  map <- input1$map
  # map <- lapply(par, fill_vals)

  # SCAA: pad NAA_re with 0, map all ages to NA
  # Rec: pad NAA_re with 0, map ages > 1 to NA (estimate Rec devs as RE)
  # all: pad NAA_re with 0, estimate all ages as RE
  # if(data$n_NAA_sigma > 0){
  par$log_NAA <- rbind(par$log_NAA, matrix(NA, nrow=proj.opts$n.yrs, ncol=data$n_ages))
  tmp <- par$log_NAA
  ind.NA <- which(is.na(tmp))
  tmp[-ind.NA] <- NA
  tmp[ind.NA] <- 1:length(ind.NA)
  par$log_NAA[ind.NA] <- 10 

  tmp <- par$log_NAA
  if(data$n_NAA_sigma < 2) tmp[,-1] <- NA # don't estimate NAA_devs for ages > 1 if SCAA or RE on recruitment
  if(data$n_NAA_sigma == 0) tmp[which(tmp[,1] == 10),1] <- NA # for SCAA, don't estimate fixed effect rec devs in projection years
  ind.notNA <- which(!is.na(tmp))
  tmp[ind.notNA] <- 1:length(ind.notNA)
  map$log_NAA = factor(tmp)
  # }

  # pad Ecov_re for projections
  if(any(data$Ecov_model > 0)){
    if(end.beyond < proj.opts$n.yrs){ # need to pad Ecov_re
      for(i in 1:(proj.opts$n.yrs-end.beyond)){
        if(proj.opts$use.last.ecov){ # use last Ecov (pad Ecov_re but map to NA)
          Ecov.proj[i,] <- model$rep$Ecov_re[data$ind_Ecov_out_end+1+end.beyond,]
        }
        if(!is.null(proj.opts$avg.ecov.yrs)){ # use average Ecov (pad Ecov_re but map to NA)
          avg.yrs.ind.ecov <- match(proj.opts$avg.ecov.yrs, input1$years)
          Ecov.proj[i,] <- avg_cols(as.matrix(model$rep$Ecov_re[avg.yrs.ind.ecov,]))
        }
        if(proj.opts$cont.ecov){ # continue Ecov process (pad Ecov_re and estimate)
          Ecov.proj[i,] <- rep(0, data$n_Ecov)
        }
        if(!is.null(proj.opts$proj.ecov)){ # use specified Ecov, have to back-calculate Ecov_re from Ecov_x
          for(j in 1:data$n_Ecov){ 
            #random walk
            if(data$Ecov_model[j] == 1) Ecov.proj[i,j] <- proj.opts$proj.ecov[i,j]
            #AR(1)
            if(data$Ecov_model[j] == 2) Ecov.proj[i,j] <- proj.opts$proj.ecov[i,j] - par$Ecov_process_pars[1,j] 
          }
        }
      }
      par$Ecov_re <- rbind(par$Ecov_re, Ecov.proj) # pad Ecov_re if necessary

      # pad map$Ecov_re
      tmp.re <- matrix(1:length(par$Ecov_re), dim(par$Ecov_re)[1], data$n_Ecov, byrow=FALSE)
      data$Ecov_use_re <- matrix(0, nrow=data$n_years_Ecov + data$n_years_proj_Ecov, ncol=data$n_Ecov)
      for(i in 1:data$n_Ecov){
        tmp.re[,i] <- if(data$Ecov_model[i]==0) rep(NA,dim(par$Ecov_re)[1]) else tmp.re[,i]
        if(data$Ecov_model[i]==1) tmp.re[1,i] <- NA # if Ecov is a rw, first year of Ecov_re is not used bc Ecov_x[1] uses Ecov1 (fixed effect)
      }
      if(!proj.opts$cont.ecov) tmp.re[1:(proj.opts$n.yrs-end.beyond)+data$n_years_Ecov,] <- NA
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

  return(list(data=data, par = par, map = map, random = input1$random,
    years = input1$years, years_full = c(input1$years, tail(input1$years,proj.opts$n.yrs) + proj.opts$n.yrs),
    ages.lab = input1$ages.lab, model_name = input1$model_name))
}
