#' Check convergence of WHAM model
#'
#' Access quick convergence checks from `TMB` and `nlminb`.
#'
#' @param mod output from \code{\link{fit_wham}}
#' @param ret T/F, return list? Default = FALSE, just prints to console
#'
#' @return a list with at least the first three of these components:
#'   \describe{
#'     \item{\code{$convergence}}{From \code{\link[stats:nlminb]{stats::nlminb}}, "0 indicates successful convergence for nlminb"}
#'     \item{\code{$maxgr}}{Max absolute gradient value, from `max(abs(mod$gr(mod$opt$par)))`}
#'     \item{\code{$maxgr_par}}{Name of parameter with max gradient}
#'     \item{\code{$is_sdrep}}{If \code{\link[TMB:sdreport]{TMB::sdreport}} was performed
#'     for this model, this indicates whether it performed without error}
#'     \item{\code{$na_sdrep}}{If \code{\link[TMB:sdreport]{TMB::sdreport}} was performed
#'     without error for this model, this indicates which (if any) components of the
#'     diagonal of the inverted hessian were returned as NA}
#'   }
#'
#' @export
#' @seealso \code{\link{fit_wham}}, \code{\link{fit_tmb}}, \code{\link[stats:nlminb]{stats::nlminb}}
#' @examples
#' \dontrun{
#' data("input4_SNEMAYT") # load SNEMA yellowtail flounder data and parameter settings
#' mod = fit_wham(input4_SNEMAYT) # using default values
#' check_convergence(mod)
#' }
check_convergence <- function(mod, ret=FALSE, f = ""){
  app = f != ""
  res <- list()
  res$convergence <- mod$opt$convergence
  res$maxgr <- max(abs(mod$final_gradient))
  res$maxgr_par <- names(mod$par)[which.max(mod$final_gradient)]
  # print to screen
  if(res$convergence == 0){
    cat("stats:nlminb thinks the model has converged: mod$opt$convergence == 0\n", file = f)
  } else {
    cat("stats:nlminb thinks the model has NOT converged: mod$opt$convergence != 0\n", file = f, append = app)
  }
  cat("Maximum gradient component:",formatC(res$maxgr, format = "e", digits = 2),"\n", file = f, append = app)
  cat("Max gradient parameter:",res$maxgr_par,"\n", file = f, append = app)
  if("sdrep" %in% names(mod)){
    res$is_sdrep = mod$is_sdrep
    if(res$is_sdrep){
      res$na_sdrep = mod$na_sdrep
      if(!(res$na_sdrep)) cat("TMB:sdreport() was performed successfully for this model\n", file = f, append = app)
      else cat("TMB:sdreport() was performed for this model, but it appears hessian was not invertible\n", file = f, append = app)
    }
    else cat("TMB:sdreport() was performed for this model, but it appears hessian was not invertible\n", file = f, append = app)
  }
  else cat("TMB:sdreport() was not performed for this model\n", file = f, append = app)

  if(ret) return(res)
}
