#' Check convergence of WHAM model
#'
#' Access quick convergence checks from `TMB` and `nlminb`.
#'
#' @param mod output from \code{\link{fit_wham}}
#' @param ret T/F, return list? Default = FALSE, just prints to console
#'
#' @return a list with three components:
#'   \describe{
#'     \item{\code{$convergence}}{From \code{\link[stats:nlminb]{stats::nlminb}}, "0 indicates successful convergence"}
#'     \item{\code{$maxgr}}{Max gradient, from `max(mod$gr(mod$opt$par))`}
#'     \item{\code{$maxgr_par}}{Name of parameter with max gradient}
#'   }
#'
#' @export
#' @seealso \code{\link{fit_wham}}, \code{\link{fit_tmb}}, \code{\link[stats:nlminb]{stats::nlminb}}
#' @examples
#' \dontrun{
#' data("SNEMA_ytl") # load SNEMA yellowtail flounder data and parameter settings
#' mod = fit_wham(input) # using default values
#' check_convergence(mod)
#' }
check_convergence <- function(mod, ret=FALSE){
  res <- list()
  res$convergence <- mod$opt$convergence
  res$maxgr <- max(mod$final_gradient)
  res$maxgr_par <- names(mod$par)[which.max(mod$final_gradient)]

  # print to screen
  if(res$convergence == 0){
    cat("The model appears converged: mod$opt$convergence == 0\n")
  } else {
    cat("The model does NOT appear converged: mod$opt$convergence != 0\n")
  }
  cat("Maximum gradient component:",formatC(res$maxgr, format = "e", digits = 2),"\n")
  cat("Max gradient parameter:",res$maxgr_par,"\n")

  if(ret) return(res)
}
