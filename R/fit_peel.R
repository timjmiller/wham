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
  if(!retro.silent) print(peel)
  temp = input
  n_years = temp$data$n_years_catch = temp$data$n_years_indices = temp$data$n_years_model = temp$data$n_years_model - peel
  temp$data$which_F_age = temp$data$which_F_age[1:n_years]
  temp$data$avg_years_ind = temp$data$avg_years_ind - peel
  temp$data$avg_years_ind <- temp$data$avg_years_ind[which(temp$data$avg_years_ind>0)] #in case peels going all the way back near initial year.
  peel_R_avg_yrs <- which((temp$data$XSPR_R_avg_yrs+1) %in% (1:n_years))
  if(length(peel_R_avg_yrs)) {
    temp$data$XSPR_R_avg_yrs <- temp$data$XSPR_R_avg_yrs[peel_R_avg_yrs]
  } else{
    warning("Years specified to use for average recruitment in SPR-based reference points are not included in peel ", peel, " so all years in peel will be used.")
    temp$data$XSPR_R_avg_yrs <- 1:n_years - 1
  }
  #temp$data$XSPR_R_avg_yrs = temp$data$XSPR_R_avg_yrs[1:n_years]
  
  # peeling ecov is tricky bc ecov_years can be different than model_years - make sure to peel to same year
  if(any(temp$data$Ecov_model != 0)){
    # temp$data$n_years_Ecov = model$dat$n_years_Ecov - peel
    n.beyond <- tail(temp$data$Ecov_year,1) - tail(temp$years,1) # n years ecov extends beyond model (if > 1 ecov, obs dim = longest)
    peel.ecov <- peel + max(n.beyond, 0) # cannot be less than peel bc prepare_wham_input pads Ecov_year if ecov ends before model
    temp$data$n_years_Ecov <- temp$data$n_years_Ecov - peel.ecov
    Ecov_re_na_ind = matrix(NA, temp$data$n_years_Ecov, temp$data$n_Ecov)
    for(i in 1:length(temp$data$Ecov_model)){
      if(temp$data$Ecov_model[i] != 1) Ecov_re_na_ind[,i] = 1
      #if(temp$data$Ecov_model[i] == 1) Ecov_re_na_ind[,i] = rbind(matrix(rep(NA, temp$data$n_Ecov)), matrix(1:((temp$data$n_years_Ecov-peel-1)*temp$data$n_Ecov), nrow=temp$data$n_years_Ecov-peel-1, ncol=temp$data$n_Ecov), matrix(rep(NA, peel*temp$data$n_Ecov), nrow=peel))
    }
    Ecov_re_na_ind = rbind(Ecov_re_na_ind, matrix(NA, peel.ecov, temp$data$n_Ecov))
    if(sum(!is.na(Ecov_re_na_ind))) Ecov_re_na_ind[!is.na(Ecov_re_na_ind)] = 1:sum(!is.na(Ecov_re_na_ind))
    temp$map$Ecov_re = factor(Ecov_re_na_ind)
    temp$data$ind_Ecov_out_end = temp$data$ind_Ecov_out_end - peel # reduce by model dim, not ecov dim
    temp$data$Ecov_use_obs[(temp$data$n_years_Ecov+1):(temp$data$n_years_Ecov+peel.ecov), ] <- 0
  }

  #peel any q random effects
  n_q_re <- sum(temp$data$use_q_re >0)
  if(n_q_re >0) {
    ind = which(temp$data$use_q_re >0)
    tmp = matrix(NA, n_years + peel, temp$data$n_indices)
    tmp[1:n_years,ind] = 1:(n_years*n_q_re)
    temp$map$q_re = factor(tmp)
  }
  
  # if("log_NAA" %in% temp$random){
    tmp <- rbind(matrix(1:(temp$data$n_ages*(n_years-1)), n_years-1), matrix(rep(NA, peel*temp$data$n_ages), peel))
    if(temp$data$n_NAA_sigma < 2) tmp[,-1] <- NA # always estimate Rec devs (col 1), whether random effect or not
    ind.notNA <- which(!is.na(tmp))
    tmp[ind.notNA] <- 1:length(ind.notNA)
    temp$map$log_NAA = factor(tmp)
  # }

  F_devs_na_ind = rbind(matrix(1:(temp$data$n_fleets * (n_years-1)), n_years-1), matrix(rep(NA, peel * temp$data$n_fleets), peel))
  temp$map$F_devs = factor(F_devs_na_ind)

  if(any(temp$data$Ecov_obs_sigma_opt %in% c(3,4))){
    temp$map$Ecov_obs_logsigma = factor(rbind(head(matrix(as.numeric(as.character(temp$map$Ecov_obs_logsigma)), ncol=temp$data$n_Ecov), -peel.ecov), matrix(NA, ncol=temp$data$n_Ecov, nrow=peel.ecov)))
    #temp$map$Ecov_obs_logsigma_re = factor(rbind(head(matrix(as.numeric(as.character(temp$map$Ecov_obs_logsigma_re)), ncol=temp$data$n_Ecov), -peel.ecov), matrix(NA, ncol=temp$data$n_Ecov, nrow=peel.ecov)))
  }

  temp.mod <- TMB::MakeADFun(temp$data, temp$par, DLL="wham", random = temp$random, map = temp$map, silent = MakeADFun.silent)
  out = fit_tmb(temp.mod, do.sdrep = do.sdrep, n.newton = n.newton, do.check=FALSE)
  if(save.input){
    out$input <- temp
    out$years <- head(input$years, length(input$years) - peel)
    out$years_full <- head(input$years_full, length(input$years_full) - peel)
    out$input$years <- out$years
    out$input$years_full <- out$years_full
    out$input$model_name <- paste0(input$model_name, " peel ",peel)
    out$ages.lab <- input$ages.lab
    out$model_name <- paste0(input$model_name, " peel ",peel)
  }
  return(out)
}
