#' Fit model peeling off \emph{i} years of data
#'
#' Internal function called by \code{\link{retro}} for \emph{i} in 1--\code{n.peels}.
#' Fits the model peeling off \emph{i} years of data (calls \code{\link{fit_tmb}}).
#'
#' @param peel Integer, number of years of data to remove before model fitting.
#' @param input input with same structure as that provided by \code{\link{prepare_wham_input}}. May want to use input$par = model$parList to start at MLEs.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? Default = \code{FALSE}.
#' @param n.newton integer, number of additional Newton steps after optimizafit_tmbtion for each peel. Default = \code{3}.
#' @param MakeADFun.silent T/F, Passed to silent argument of \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}. Default = \code{FALSE}.
#' @param retro.silent T/F, Passed to argument of internal fit_peel function. Determines whether peel number is printed to screen. Default = \code{FALSE}.
#' @param save.input T/F, should modified input list be saved? Necessary to project from a peel but increases model object size. Default = \code{FALSE}.
#'
#' @return \code{out}, output of \code{\link{fit_tmb}} for peel \emph{i}
#'
#' @export
#' 
#' @seealso \code{\link{fit_wham}}, \code{\link{retro}}, \code{\link{fit_tmb}}
#'
fit_peel = function(peel, input, do.sdrep = FALSE, n.newton = 3, MakeADFun.silent = FALSE, retro.silent = FALSE, save.input = FALSE)
{
  out = list()
  if(!retro.silent) print(paste0("Retro Peel: ", peel))
  temp = input
  data <- temp$data
  par <- temp$par
  map <- temp$map
  #n_years = data$n_years_model = temp$data$n_years_model - peel  
  n_years <- temp$data$n_years_model - peel
  #data$which_F_age = data$which_F_age[1:n_years,]
  data$use_agg_catch[n_years + 1:peel,] <- 0
  data$use_catch_paa[n_years + 1:peel,] <- 0
  data$use_indices[n_years + 1:peel,] <- 0
  data$use_index_paa[n_years + 1:peel,] <- 0
  #new year indicator for which random effects to include in likelihoods used primarily for peels (not selectivity though?)
  data$years_use <- 1:n_years - 1
  data$do_annual_SPR_BRPs <- 0 #don't do these calculations

#print("in fit peel")
  # peeling ecov is tricky bc ecov_years can be different than model_years - make sure to peel to same year
  if(any(data$Ecov_model != 0)){
    n.beyond <- tail(data$Ecov_year,1) - tail(temp$years,1) # n years ecov extends beyond model (if > 1 ecov, obs dim = longest)
    print(n.beyond)
    peel.ecov <- peel + max(n.beyond, 0) # cannot be less than peel bc prepare_wham_input pads Ecov_year if ecov ends before model
    print(peel.ecov)
    n_years_Ecov <- data$n_years_Ecov - peel.ecov
    #new year indicator for which random effects to include in likelihoods used primarily for peels
    data$years_use_Ecov <- 1:n_years_Ecov - 1
    if(!is.null(map$Ecov_re)) map$Ecov_re <- matrix(as.integer(map$Ecov_re), NROW(par$Ecov_re), NCOL(par$Ecov_re))
    else map$Ecov_re <- matrix(NA, NROW(par$Ecov_re), NCOL(par$Ecov_re))
    #Ecov_re_na_ind = matrix(NA, data$n_years_Ecov, data$n_Ecov)
    k <- 0
    for(i in 1:length(data$Ecov_model)){
      #Ecov_re for RW models apparently were always fixed before? 
      if(data$Ecov_model[i] == 2) {
        map$Ecov_re[1:n_years_Ecov,i] = k  + 1:n_years_Ecov
        k <- k + n_years_Ecov
      }
      if(data$Ecov_model[i] == 1){ #first Ecov_re not used
        map$Ecov_re[2:n_years_Ecov,i] = k  + 2:n_years_Ecov
        k <- k + n_years_Ecov - 1
      }
      map$Ecov_re[n_years_Ecov + 1:peel.ecov,i] <- NA
    }
    map$Ecov_re = factor(map$Ecov_re)
    #Ecov_re_na_ind = rbind(Ecov_re_na_ind, matrix(NA, peel.ecov, data$n_Ecov))
    #if(sum(!is.na(Ecov_re_na_ind))) Ecov_re_na_ind[!is.na(Ecov_re_na_ind)] = 1:sum(!is.na(Ecov_re_na_ind))
    #temp$map$Ecov_re = factor(Ecov_re_na_ind)
    #data$ind_Ecov_out_end = data$ind_Ecov_out_end - peel # reduce by model dim, not ecov dim
    print(dim(data$Ecov_use_obs))
    print(data$n_years_Ecov)
    data$Ecov_use_obs[n_years_Ecov + 1:peel.ecov, ] <- 0
    #data$Ecov_use_obs[(data$n_years_Ecov+1):(data$n_years_Ecov+peel.ecov), ] <- 0
  }

  #if(any(data$Ecov_obs_sigma_opt %in% c(3,4))){
  if(any(data$Ecov_obs_sigma_opt == 3)){
    temp$map$Ecov_obs_logsigma = factor(rbind(head(matrix(as.numeric(as.character(temp$map$Ecov_obs_logsigma)), ncol=data$n_Ecov), -peel.ecov), matrix(NA, ncol=data$n_Ecov, nrow=peel.ecov)))
    #temp$map$Ecov_obs_logsigma_re = factor(rbind(head(matrix(as.numeric(as.character(temp$map$Ecov_obs_logsigma_re)), ncol=data$n_Ecov), -peel.ecov), matrix(NA, ncol=data$n_Ecov, nrow=peel.ecov)))
  }
#print("in fit peel")

  # #peel any q random effects
  map$q_re <- matrix(as.integer(map$q_re), NROW(par$q_re), NCOL(par$q_re))
  map$q_re[n_years + 1:peel] <- NA
  map$q_re <- factor(map$q_re)
#print("in fit peel")
  
  #peel any NAA fixed effects
  #if(any(data$NAA_re_model == 0)){
    map$log_NAA <- array(as.integer(map$log_NAA), dim = dim(par$log_NAA))
    map$log_NAA[,,n_years-1 + 1:peel,] <- NA
    #print(dim(map$log_NAA))
    #print(dim(par$log_NAA))
    map$log_NAA <- factor(map$log_NAA)
  #}
#print("in fit peel")

  #peel any M random effects
  map$M_re <- array(as.integer(map$M_re), dim = dim(par$M_re))
  map$M_re[,,n_years + 1:peel,] <- NA
  print(dim(map$M_re))
  print(dim(par$M_re))
  map$M_re <- factor(map$M_re)

  # #peel any mu random effects
  map$mu_re <- array(as.integer(map$mu_re), dim = dim(par$mu_re))
  map$mu_re[,,,n_years + 1:peel,,] <- NA
  print(dim(map$mu_re))
  print(dim(par$mu_re))
  map$mu_re <- factor(map$mu_re)

  map$F_pars <- rbind(matrix(1:(data$n_fleets * n_years), n_years), matrix(rep(NA, peel * data$n_fleets), peel))
  # print(dim(map$F_pars))
  # print(dim(par$F_pars))
  map$F_pars <- factor(map$F_pars)
# print("in fit peel")

  temp$data <- data
  temp$map <- map
  temp$par <- par
  # print(lapply(par, length))
  # print("before TMB call")
  temp.mod <- TMB::MakeADFun(temp$data, temp$par, DLL="wham", random = temp$random, map = temp$map, silent = MakeADFun.silent)
  print(" in 1")
  print(temp.mod$fn())
  # trep <- temp.mod$report()
  # print(sapply(trep[grep("nll_", names(trep))], sum))
  # print(trep$nll_agg_catch)
  # print(trep$pred_catch)
  # print(trep$pred_log_catch)
  out = fit_tmb(temp.mod, do.sdrep = do.sdrep, n.newton = n.newton, do.check=FALSE)
#  print(" in 2")
  if(save.input){
    out$input <- temp
    # out$years <- head(input$years, length(input$years) - peel)
    # out$years_full <- head(input$years_full, length(input$years_full) - peel)
    # out$input$years <- out$years
    # out$input$years_full <- out$years_full
    out$input$model_name <- paste0(input$model_name, " peel ",peel)
    out$ages.lab <- input$ages.lab
    out$model_name <- paste0(input$model_name, " peel ",peel)
  }
  return(out)
}
