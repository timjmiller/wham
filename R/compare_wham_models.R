#' Compare WHAM models
#'
#' After fitting multiple WHAM models, \code{compare_wham_models} produces a table
#' of AIC and Mohn's rho to aid model comparison.
#'
#' @param mods (named) list of fit WHAM models. If no names are given, "m1", "m2", ...
#' will be used.
#' @param fname character, filename to save CSV results table (.csv will be appended). Default = "model_comparison".
#' @param sort T/F, sort by AIC? Default = TRUE.
#' @param mohns.rho T/F, calculate Mohn's rho? Retrospective analysis must have been run for all modes. Default = TRUE.
#'
#' @return a list with the following components:
#'   \describe{
#'     \item{\code{daic}}{Vector of delta-AIC by model}
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
compare_wham_models <- function(mods, fname = "model_comparison", sort = TRUE, calc.rho = TRUE, calc.aic = TRUE){
  if(is.null(names(mods))) names(mods) <- paste0("m",1:length(mods))
  aic.tab <- aic <- daic <- NULL
  if(calc.aic){
    if(sum(mapply(function(x) x$env$data$Ecov_model==0, mods)) %in% c(0,length(mods))){
      ecov.obs <- lapply(mods, function(x) x$env$data$Ecov_use_obs)
      all.identical <- function(l) all(mapply(identical, head(l, 1), tail(l, -1)))
      if(!all.identical(ecov.obs)){
        stop("Different env covariate data used, cannot compare models' AIC.
             Set 'calc.aic = FALSE' to compare only Mohn's rho, or only select
             models fit to the same data.")
      }
      aic <- sapply(mods, function(x){
        k = length(x$opt$par)
        2*(x$opt$obj + k) # AIC
        # 2*(x$opt$obj + k + k*(k+1)/(n-k-1)) # AICc
      })
      aic <- round(aic, 1)
      daic <- round(aic - min(aic), 1)
    } else {
      stop("Env covariate in some model(s) but not all. Cannot compare AIC
           for models with different data (here, some have environmental data
           and some do not).")
    }
    aic.tab <- cbind(daic, aic)
    colnames(aic.tab) <- c("dAIC","AIC")
  }
  rho <- NULL
  if(calc.rho){
    if(any(mapply(function(x) is.null(x$peels), mods))){
      stop("Not all models have peels --> Cannot compare Mohn's rho.
           Set 'calc.rho = FALSE' to compare only AICc, or re-run models
           with 'fit_wham(do.retro = TRUE)'.")
    }
    rho <- t(sapply(mods, function(x){
      mohns_rho(x)
    }))[ ,c("R","SSB","Fbar")]
    rho <- round(rho, 4)
    colnames(rho) <- paste0("rho_",c("R","SSB","Fbar"))
    # apply(rho, 1, function(y) mean(abs(y)))
  }
  tab <- cbind(aic.tab, rho)

  best <- NULL
  if(calc.aic) best <- names(mods)[which(aic == min(aic))]
  if(sort){
    ord <- order(aic)
    daic <- daic[ord]
    aic <- aic[ord]
    rho <- rho[ord,]
    tab <- tab[ord,]
  }
  write.csv(tab, file = paste0(file.path(getwd(),fname),".csv"))

  print(tab) # print to console
  return(list(daic=daic, aic=aic, rho=rho, best=best, tab=tab))
}
