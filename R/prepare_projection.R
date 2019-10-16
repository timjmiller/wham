#' Prepare input data and parameters to project an already fit WHAM model
#'
#' \code{prepare_projection} is an internal function called by \code{\link{project_wham}},
#' which in turn is called by \code{\link{fit_wham}} if \code{do.proj = TRUE}.
#'
#' @param model a previously fit wham model
#' @param proj.opts a named list with the following components:
#'   \describe{
#'     \item{\code{$n.yrs}}{integer, number of years to project/forecast. Default = \code{3}.}
#'     \item{\code{$use.lastF}}{T/F, use terminal year F for projections. Default = \code{TRUE}.}
#'     \item{\code{$use.avgF}}{T/F, use average F for projections.}
#'     \item{\code{$use.FXSPR}}{T/F, calculate F at X% SPR for projections.}
#'     \item{\code{$proj.F}}{vector, user-specified fishing mortality for projections. Length must equal \code{n.yrs}.}
#'     \item{\code{$proj.catch}}{vector, user-specified aggregate catch for projections. Length must equal \code{n.yrs}.}
#'     \item{\code{$avg.yrs}}{vector, specify which years to average over for calculating reference points. Default = last 5 model years, \code{tail(model$years, 5)}.}
#'     \item{\code{$cont.Ecov}}{T/F, continue Ecov process (e.g. random walk or AR1) for projections. Default = \code{TRUE}.}
#'     \item{\code{$use.last.Ecov}}{T/F, use terminal year Ecov for projections.}
#'     \item{\code{$avg.Ecov.yrs}}{vector, specify which years to average over the environmental covariate(s) for projections.}
#'     \item{\code{$proj.Ecov}}{vector, user-specified environmental covariate(s) for projections. Length must equal \code{n.yrs}.}
#'   }
#'
#' @return same as \code{prepare_wham_input}, a list ready for \code{fit_wham}:
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
#                 avg.yrs=NULL, cont.Ecov=TRUE, use.last.Ecov=FALSE, avg.Ecov.yrs=NULL, proj.Ecov=NULL)

  # default: use average M, selectivity, etc. over last 5 model years to calculate ref points
  if(is.null(proj.opts$avg.yrs)) proj.opts$avg.yrs <- tail(model$years, 5)

  # check options for F/catch are valid 
  if(any(proj.opts$avg.yrs %in% model$years == FALSE)) stop(paste("","** Error setting up projections: **",
    "proj.opts$avg.yrs is not a subset of model years.","",sep='\n'))
  F.opt.ct <- sum(proj.opts$use.last.F, proj.opts$use.FXSPR, !is.null(proj.opts$proj.F), !is.null(proj.opts$proj.catch))
  if(F.opt.ct != 1) stop(paste("","** Error setting up projections: **",
    "Exactly one method of specifying F must be used (see ?project_wham).",
    "You have specified these in 'proj.opts':",
    capture.output(cat("  $use.last.F = ",proj.opts$use.last.F)),
    capture.output(cat("  $use.avg.F = ",proj.opts$use.avg.F)),
    capture.output(cat("  $use.FXSPR = ",proj.opts$use.FXSPR)),
    capture.output(cat("  $proj.F = ",proj.opts$proj.F)),
    capture.output(cat("  $proj.catch = ",proj.opts$proj.catch)),"",sep='\n'))

  # check Ecov options are valid 
  input1 <- model$input
  data <- input1$data
  if(any(input1$data$Ecov_model > 0)){
    # if Ecov already extends beyond population model, Ecov proj options are unnecessary
    n.beyond <- data$n_years_Ecov-1-data$ind_Ecov_out_end
    end.beyond <- min(n.beyond, proj.opts$n.yrs)
    Ecov.proj <- matrix(rep(NA, (proj.opts$n.yrs-end.beyond)*dim(model$rep$Ecov_re)[2]), ncol=dim(model$rep$Ecov_re)[2], byrow=TRUE)
    if(end.beyond == proj.opts$n.yrs){ print("Ecov already fit through projection years. Using fit Ecov for projections...")
    } else { # if Ecov proj options ARE necessary, check they are valid
      Ecov.opt.ct <- sum(proj.opts$cont.Ecov, proj.opts$use.last.Ecov, !is.null(proj.opts$avg.Ecov.yrs), !is.null(proj.opts$proj.Ecov))
      if(Ecov.opt.ct != 1) stop(paste("","** Error setting up projections: **",
        "Exactly one method of specifying Ecov must be used (see ?project_wham).",
        "You have specified these in 'proj.opts':",
        capture.output(cat("  $cont.Ecov = ",proj.opts$cont.Ecov)),
        capture.output(cat("  $use.last.Ecov = ",proj.opts$use.last.Ecov)),
        capture.output(cat("  $avg.Ecov.yrs = ",proj.opts$avg.Ecov.yrs)),
        capture.output(cat("  $proj.Ecov = ",proj.opts$proj.Ecov)),"",sep='\n'))
      if(!is.null(proj.opts$avg.Ecov.yrs) & any(proj.opts$avg.Ecov.yrs %in% model$years == FALSE)) stop(paste("","** Error setting up projections: **",
        "proj.opts$avg.Ecov.yrs is not a subset of model years.","",sep='\n'))
    } 
  }

  # add new data objects for projections
  data$do_proj = 1
  data$n_years_proj = proj.opts$n.yrs
  avg.yrs.ind <- match(proj.opts$avg.yrs, input1$years)
  data$avg_years_ind = avg.yrs.ind - 1 # c++ indices start at 0
  if(proj.opts$use.last.F) data$proj_F_opt = 1
  if(proj.opts$use.avg.F) data$proj_F_opt = 2
  if(proj.opts$use.FXSPR) data$proj_F_opt = 3
  if(!is.null(proj.opts$proj.F)){
    data$proj_F_opt = 4
    data$proj_Fcatch = proj.opts$proj.F
  } 
  if(!is.null(proj.opts$proj.catch)){
    data$proj_F_opt = 5
    data$proj_Fcatch = proj.opts$proj.catch
  }
  if(data$proj_F_opt %in% 1:3) data$proj_Fcatch = rep(0, proj.opts$n.yrs)

  # modify data objects for projections (pad with average over avg.yrs): mature, fracyr_SSB, waa
  avg_cols = function(x) apply(x, 2, mean, na.rm=TRUE)
  data$mature <- rbind(data$mature, matrix(rep(avg_cols(data$mature[avg.yrs.ind,]), proj.opts$n.yrs), nrow=proj.opts$n.yrs, byrow=TRUE))
  data$fracyr_SSB <- c(data$fracyr_SSB, rep(mean(data$fracyr_SSB[avg.yrs.ind]), proj.opts$n.yrs))
  toadd <- apply(data$waa[,avg.yrs.ind,], c(1,3), mean)
  tmp <- array(NA, dim = dim(data$waa) + c(0,proj.opts$n.yrs,0))
  tmp[,1:data$n_years_model,] <- data$waa
  for(i in seq(data$n_years_model+1,data$n_years_model+proj.opts$n.yrs)) tmp[,i,] <- toadd
  data$waa <- tmp

  # initialize and fix pars at previously estimated values
  par <- model$parList
  fill_vals <- function(x){as.factor(rep(NA, length(x)))}
  map <- lapply(par, fill_vals)
  
  # pad parameters for projections: log_NAA and Ecov_re
  par$log_NAA <- rbind(par$log_NAA, matrix(NA, nrow=proj.opts$n.yrs, ncol=data$n_ages))
  map$log_NAA <- as.factor(c(map$log_NAA, seq(1:(proj.opts$n.yrs*data$n_ages))))
  if(any(data$Ecov_model > 0)){
    if(end.beyond < proj.opts$n.yrs){ # need to pad Ecov_re
      for(i in 1:(proj.opts$n.yrs-end.beyond)){
        if(proj.opts$use.last.Ecov){ # use last Ecov
          Ecov.proj[i,] <- model$rep$Ecov_re[data$ind_Ecov_out_end+1+end.beyond,]
        }
        if(!is.null(proj.opts$avg.Ecov.yrs)){ # use average Ecov
          avg.yrs.ind.Ecov <- match(proj.opts$avg.Ecov.yrs, input1$years)
          Ecov.proj[i,] <- avg_cols(model$rep$Ecov_re[avg.yrs.ind.Ecov,])
        }
        if(proj.opts$cont.Ecov){ # continue Ecov process
          Ecov.proj[i,] <- rep(0, data$n_Ecov) 
        }
        # if(!is.null(proj.opts$proj.Ecov)){ # use specified Ecov, have to back-calculate Ecov_re from Ecov_x
        # }
      }
    }
    par$Ecov_re <- rbind(model$rep$Ecov_re, Ecov.proj) # pad Ecov_re if necessary
    map$Ecov_re <- as.factor(c(map$Ecov_re, seq(1:length(Ecov.proj))))
  }

  # remove random effects if they are not estimated (all mapped to NA)
  check_allNA <- function(x){ifelse(length(levels(map[[x]])) > 0, FALSE, TRUE)}
  random <- input1$random[!sapply(input1$random, check_allNA)]

  return(list(data=data, par = par, map = map, random = random, 
    years = c(input1$years, tail(input1$years,proj.opts$n.yrs) + proj.opts$n.yrs),
    ages.lab = input1$ages.lab, model_name = input1$model_name))
}