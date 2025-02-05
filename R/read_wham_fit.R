#' Read WHAM fit
#'
#' Gets output from a fit WHAM model for plotting with other models.
#' Internal function, called within \code{\link{compare_wham_models}}.
#'
#' @param mod output from \code{\link{fit_wham}}
#' @param alphaCI (1-alpha)\% confidence intervals will be calculated. Default = 0.05 for 95\% CI.  
#'
#' @return a named list with the following elements:
#'   \describe{
#'     \item{\code{$years}}{numeric vector, model years only, e.g. \code{1972:2020}}
#'     \item{\code{$years_full}}{numeric vector, model + proj years, e.g. \code{1972:2022}}
#'     \item{\code{$selAA}}{list of length(n_selblocks), first the fleet blocks then indices, i.e. if 4 fleet blocks and 3 indices, \code{selAA[[5]]} is for index 1. Each element is a matrix, years (rows) x ages (cols), selectivity at age}
#'     \item{\code{$selblock_pointer_fleets}}{matrix, years x fleets, indices of selAA used by each fleet in each year}
#'     \item{\code{$selblock_pointer_indices}}{matrix, years x indices, indices of selAA used by each index in each year}
#'     \item{\code{$MAA}}{array, stocks x regions x years x ages, natural mortality}
#'     \item{\code{$log_SSB}}{matrix, years x 2, log-scale spawning stock biomass. 1st col = MLE, 2nd col = SE.}
#'     \item{\code{$log_F}}{matrix, years x 2, log-scale fully-selected F. 1st col = MLE, 2nd col = SE.}
#'     \item{\code{$log_NAA_rep}}{array, stocks x regions x years x ages, numbers at age}
#'     \item{\code{$NAA_CV}}{array, stocks x regions x years x ages, CV of numbers at age}
#'     \item{\code{$percentSPR}}{scalar, X\% SPR used to calculate reference points, default = 40}
#'     \item{\code{$log_Y_FXSPR}}{matrix, years x 2, log-scale yield at FXSPR. 1st col = MLE, 2nd col = SE.}
#'     \item{\code{$log_FXSPR}}{matrix, years x 2, log-scale FXSPR. 1st col = MLE, 2nd col = SE.}
#'     \item{\code{$log_SSB_FXSPR}}{matrix, years x 2, log-scale SSB at FXSPR. 1st col = MLE, 2nd col = SE.}
#'     \item{\code{$log_rel_ssb_F_cov}}{list, length n_years, each element is a 2x2 covariance matrix with SSB/SSB_FXSPR first and F/F_FXSPR second}
#'   }
#'
#' @export
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{read_asap3_fit}}, \code{\link{compare_wham_models}}
#'
read_wham_fit <- function(mod, alphaCI=0.05){
  # if sdreport succeeded but didn't save full sdreport object in mod, recalculate it here
  if(mod$is_sdrep & class(mod$sdrep)[1] != "sdreport"){ #only summary saved
    mod$sdrep <- TMB::sdreport(mod)
  }
  n_ages <- mod$env$data$n_ages
  n_years <- length(mod$years_full)

  x <- list()

  if(!is.null(mod$sdrep)){
    std <- list(est = TMB:::as.list.sdreport(mod$sdrep, report=T, what = "Est"), se = TMB:::as.list.sdreport(mod$sdrep, report=T, what = "Std"))
    inds <- list()
    std.summ <- summary(mod$sdrep, "report")
  	# inds$Y.t <- matrix(which(rownames(std) == "log_Y_FXSPR"), ncol = all_catch)
  	inds$F.t <- which(rownames(std.summ) == "log_FXSPR")
  	inds$SSB.t <- matrix(which(rownames(std.summ) == "log_SSB_FXSPR"), ncol = NCOL(std$est$log_SSB_FXSPR))[,NCOL(std$est$log_SSB_FXSPR)]
  	inds$ssb <- which(rownames(std.summ) == "log_SSB_all")
  	inds$full.f <- which(rownames(std.summ) == "log_F_tot")
    rep_names <- names(std$est)
    # rep_names <- c("log_FMSY", "log_SSB_MSY", "log_MSY", "log_SSB", "log_SSB_all", "log_F_tot", "log_NAA_rep", "log_FXSPR", "log_Y_FXSPR", "log_SSB_FXSPR")
    cov <- mod$sdrep$cov
    x$log_rel_ssb_F_cov <- lapply(1:n_years, function(x){
      K <- cbind(c(1,-1,0,0),c(0,0,1,-1))
      ind <- c(inds$ssb[x],inds$SSB.t[x],inds$full.f[x],inds$F.t[x])
      tcov <- cov[ind,ind]
      ests <- std.summ[ind,1]
      diff <- t(K) %*% ests
      return(list(diff, t(K) %*% tcov %*% K))
    })
    for(i in 1:length(rep_names)){
      ci.type <- "I"
      if(length(grep("log_", rep_names[i]))) ci.type <- "exp"
      x[[rep_names[i]]] <- list(est = std$est[[rep_names[i]]], se = std$se[[rep_names[i]]], 
        ci = get.ci(std$est[[rep_names[i]]], std$se[[rep_names[i]]], type = ci.type, alpha.ci = alphaCI, asap = FALSE))
    }
    if("log_FMSY" %in% rownames(std.summ) & "log_SSB_MSY" %in% rownames(std.summ)){
      inds$SSBmsy <- matrix(which(rownames(std.summ) == "log_SSB_MSY"), ncol = NCOL(std$est$log_SSB_MSY))[,NCOL(std$est$log_SSB_MSY)]
      inds$Fmsy <- which(rownames(std) == "log_FMSY")
      x$log_rel_ssb_F_cov_msy <- lapply(1:n_years, function(x){
        K <- cbind(c(1,-1,0,0),c(0,0,1,-1))
        ind <- c(inds$ssb[x],inds$SSBmsy[x],inds$full.f[x],inds$Fmsy[x])
        tcov <- cov[ind,ind]
        ests <- std.summ[ind,1]
        diff <- t(K) %*% ests
        return(list(diff, t(K) %*% tcov %*% K))
      })    
    }
  } else {
    x[names(mod$rep)] <- lapply(mod$rep, function(x) {
      out <- list(est = x)
      out$se <- out$est
      out$se[] <- 0
      out$ci <- list(lo = out$se, hi = out$se)
      return(out)
    })
  }
  x$years <- mod$years
  x$years_full <- mod$years_full
  x$selAA <- mod$rep$selAA
  x$selblock_pointer_fleets <- mod$env$data$selblock_pointer_fleets
  x$selblock_pointer_indices <- mod$env$data$selblock_pointer_indices
  x$MAA <- mod$rep$MAA
  x$percentSPR <- mod$env$data$percentSPR
  return(x)
}

