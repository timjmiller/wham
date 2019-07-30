#' Fit model peeling off \emph{i} years of data
#'
#' Internal function called by \code{\link{retro}} for \emph{i} in 1--\code{n.peels}.
#' Fits the model peeling off \emph{i} years of data (calls \code{\link{fit_tmb}}).
#'
#' @param peel Integer, number of years of data to remove before model fitting.
#' @param model Output from \code{\link{fit_tmb}}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? Default = \code{FALSE}.
#' @param n.newton integer, number of additional Newton steps after optimization for each peel. Default = \code{3}.
#'
#' @return \code{out}, output of \code{\link{fit_tmb}} for peel \emph{i}
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{retro}}, \code{\link{fit_tmb}}
#'
fit_peel = function(peel, model, do.sdrep = FALSE, n.newton = 3)
{
  out = list()
  print(peel)
  temp = model
  n_years = temp$dat$n_years_catch = temp$dat$n_years_indices = temp$dat$n_years_model = temp$dat$n_years_model - peel
  if(temp$dat$Ecov_model != 0){
    temp$dat$n_years_Ecov = model$dat$n_years_Ecov - peel
    Ecov_re_na_ind = rbind(matrix(1:(temp$dat$n_years_Ecov*temp$dat$n_Ecov), nrow=temp$dat$n_years_Ecov, ncol=temp$dat$n_Ecov), matrix(rep(NA, peel*temp$dat$n_Ecov), nrow=peel))
    temp$map$Ecov_re = factor(Ecov_re_na_ind)
    temp$dat$ind_Ecov_out_end = model$dat$ind_Ecov_out_end - peel
  }
  log_NAA_na_ind = rbind(matrix(1:(temp$dat$n_ages*(n_years-1)), n_years-1), matrix(rep(NA, peel*temp$dat$n_ages), peel))
  F_devs_na_ind = rbind(matrix(1:(temp$dat$n_fleets * (n_years-1)), n_years-1), matrix(rep(NA, peel * temp$dat$n_fleets), peel))
  log_R_na_ind = c(1:(n_years-1), rep(NA, peel))
  if("log_R" %in% model$random) temp$map$log_R = factor(log_R_na_ind)
  if("log_NAA" %in% model$random) temp$map$log_NAA = factor(log_NAA_na_ind)
  temp$map$F_devs = factor(F_devs_na_ind)
  temp.mod <- TMB::MakeADFun(temp$dat, temp$par, DLL="wham", random = temp$random, map = temp$map)
  out = fit_tmb(temp.mod, do.sdrep = do.sdrep, n.newton = n.newton)
  return(out)
}
