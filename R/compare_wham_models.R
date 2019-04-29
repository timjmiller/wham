#' Compare WHAM models
#'
#' After fitting multiple WHAM models, \code{compare_wham_models} produces a table
#' of AIC and Mohn's rho to aid model comparison.
#'
#' @param mods (named) list of fit WHAM models. If no names are given, "m1", "m2", ...
#' will be used.
#' @param fname character, filename to save CSV results table (.csv will be appended). Default = "model_comparison".
#' @param sort T/F, sort by AIC? Default = TRUE.
#'
#' @return a list with the following components:
#'   \describe{
#'     \item{\code{aic}}{Vector of AIC by model}
#'     \item{\code{rho}}{Matrix of Mohn's rho by model}
#'     \item{\code{best}}{Name of best model (lowest AIC)}
#'     \item{\code{tab}}{Results table of AIC and Mohn's rho}
#'   }
#'
#' @seealso \code{\link{fit_wham}}
#'
#' @examples
#' \dontrun{
#' m1 <- fit_wham(input1)
#' m2 <- fit_wham(input2)
#' mods <- list(m1=m1, m2=m2)
#' res <- compare_wham_models(mods)
#' res$best
#' res$tab
#' }
#'
#' @export
compare_wham_models <- function(mods, fname = "model_comparison", sort = TRUE){
  if(is.null(names(mods))) names(mods) <- paste0("m",1:length(mods))
  aic <- sapply(mods, function(x){
    2*(x$opt$obj + length(x$opt$par))
  })
  aic <- round(aic, 1)
  rho <- t(sapply(mods, function(x){
    mohns_rho(x)
  }))[ ,c("R","SSB","Fbar")]
  rho <- round(rho, 4)
  # apply(rho, 1, function(y) mean(abs(y)))
  tab <- cbind(AIC = aic, rho = rho)

  best <- names(mods)[which(aic == min(aic))]
  ord <- order(aic)
  aic <- aic[ord]
  rho <- rho[ord,]
  tab <- tab[ord,]
  write.csv(tab, file = paste0(file.path(getwd(),fname),".csv"))

  return(list(aic=aic, rho=rho, best=best, tab=tab))
}
