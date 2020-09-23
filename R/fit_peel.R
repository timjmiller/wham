#' Fit model peeling off \emph{i} years of data
#'
#' Internal function called by \code{\link{retro}} for \emph{i} in 1--\code{n.peels}.
#' Fits the model peeling off \emph{i} years of data (calls \code{\link{fit_tmb}}).
#'
#' @param peel Integer, number of years of data to remove before model fitting.
#' @param model Output from \code{\link{fit_tmb}}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? Default = \code{FALSE}.
#' @param n.newton integer, number of additional Newton steps after optimization for each peel. Default = \code{3}.
#' @param MakeADFun.silent T/F, Passed to silent argument of \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}. Default = \code{FALSE}.
#' @param retro.silent T/F, Passed to argument of internal fit_peel function. Determines whether peel number is printed to screen. Default = \code{FALSE}.
#'
#' @return \code{out}, output of \code{\link{fit_tmb}} for peel \emph{i}
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{retro}}, \code{\link{fit_tmb}}
#'
fit_peel = function(peel, model, do.sdrep = FALSE, n.newton = 3, MakeADFun.silent = FALSE, retro.silent = FALSE)
{
  out = list()
  if(!retro.silent) print(peel)
  temp = model
  n_years = temp$dat$n_years_catch = temp$dat$n_years_indices = temp$dat$n_years_model = temp$dat$n_years_model - peel
  if(temp$dat$Ecov_model != 0){
    # temp$dat$n_years_Ecov = model$dat$n_years_Ecov - peel
    if(temp$dat$Ecov_model != 1) Ecov_re_na_ind = rbind(matrix(1:((temp$dat$n_years_Ecov-peel)*temp$dat$n_Ecov), nrow=temp$dat$n_years_Ecov-peel, ncol=temp$dat$n_Ecov), matrix(rep(NA, peel*temp$dat$n_Ecov), nrow=peel))
    if(temp$dat$Ecov_model == 1) Ecov_re_na_ind = rbind(matrix(rep(NA, temp$dat$n_Ecov)), matrix(1:((temp$dat$n_years_Ecov-peel-1)*temp$dat$n_Ecov), nrow=temp$dat$n_years_Ecov-peel-1, ncol=temp$dat$n_Ecov), matrix(rep(NA, peel*temp$dat$n_Ecov), nrow=peel))
    temp$map$Ecov_re = factor(Ecov_re_na_ind)
    temp$dat$ind_Ecov_out_end = model$dat$ind_Ecov_out_end - peel
    temp$dat$Ecov_use_obs[(temp$dat$n_years_Ecov-peel+1):temp$dat$n_years_Ecov, ] <- 0
  }
  
  if("log_NAA" %in% temp$random){
    tmp <- rbind(matrix(1:(temp$dat$n_ages*(n_years-1)), n_years-1), matrix(rep(NA, peel*temp$dat$n_ages), peel))
    if(temp$dat$n_NAA_sigma < 2) tmp[,-1] <- NA # always estimate Rec devs (col 1), whether random effect or not
    ind.notNA <- which(!is.na(tmp))
    tmp[ind.notNA] <- 1:length(ind.notNA)
    temp$map$log_NAA = factor(tmp)
  }

  F_devs_na_ind = rbind(matrix(1:(temp$dat$n_fleets * (n_years-1)), n_years-1), matrix(rep(NA, peel * temp$dat$n_fleets), peel))
  temp$map$F_devs = factor(F_devs_na_ind)

  if(temp$dat$Ecov_obs_sigma_opt %in% c(3,4)){
    temp$map$Ecov_obs_logsigma = factor(rbind(head(matrix(as.numeric(as.character(temp$map$Ecov_obs_logsigma)), ncol=temp$dat$n_Ecov), -peel), matrix(NA, ncol=temp$dat$n_Ecov, nrow=peel)))
  }

  temp.mod <- TMB::MakeADFun(temp$dat, temp$par, DLL="wham", random = temp$random, map = temp$map, silent = MakeADFun.silent)
  out = fit_tmb(temp.mod, do.sdrep = do.sdrep, n.newton = n.newton, do.check=FALSE)
  return(out)
}
