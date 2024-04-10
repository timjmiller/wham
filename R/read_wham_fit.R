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
  if(mod$is_sdrep & class(mod$sdrep)[1] != "sdreport"){
    mod$sdrep <- TMB::sdreport(mod)
  }
  n_ages <- mod$env$data$n_ages
  n_years <- length(mod$years_full)

  x <- list()
  x$years <- mod$years
  x$years_full <- mod$years_full
  x$selAA <- mod$rep$selAA
  x$selblock_pointer_fleets <- mod$env$data$selblock_pointer_fleets
  x$selblock_pointer_indices <- mod$env$data$selblock_pointer_indices
  x$MAA <- mod$rep$MAA

  # std = summary(mod$sdrep)
  # all_catch = mod$input$data$n_fleets+mod$input$data$n_regions+1
  # all_stocks = mod$input$data$n_stocks+1 
  # std <- summary(mod$sdrep, "report")
  std <- list(est = TMB:::as.list.sdreport(mod$sdrep, report=T, what = "Est"), se = TMB:::as.list.sdreport(mod$sdrep, report=T, what = "Std"))
  inds <- list()
  std.summ <- summary(mod$sdrep, "report")
	# inds$Y.t <- matrix(which(rownames(std) == "log_Y_FXSPR"), ncol = all_catch)
	inds$F.t <- which(rownames(std.summ) == "log_FXSPR")
	inds$SSB.t <- matrix(which(rownames(std.summ) == "log_SSB_FXSPR"), ncol = NCOL(std$est$log_SSB_FXSPR))[,NCOL(std$est$log_SSB_FXSPR)]
	inds$ssb <- which(rownames(std.summ) == "log_SSB_all")
	inds$full.f <- which(rownames(std.summ) == "log_F_tot")
  
  # inds <- list(Y.t = which(rownames(std) == "log_Y_FXSPR"))
  # inds$F.t <- which(rownames(std) == "log_FXSPR")
  # inds$SSB.t <- which(rownames(std) == "log_SSB_FXSPR")
  # inds$ssb <- which(rownames(std) == "log_SSB")
  # inds$faa <- which(rownames(std) == "log_FAA_tot")
  # log.faa <- matrix(std[inds$faa,1], n_years, n_ages)
  # age.full.f <- apply(log.faa,1, function(x) max(which(x == max(x))))
  # inds$full.f <- (age.full.f-1)*n_years + 1:n_years  + min(inds$faa) - 1 #cbind(1:n_years, age.full.f)
  #inds$naa <- array(which(rownames(std) == "log_NAA_rep"), dim = dim(mod$rep$NAA))
  rep_names <- names(std$est)
  # rep_names <- c("log_FMSY", "log_SSB_MSY", "log_MSY", "log_SSB", "log_SSB_all", "log_F_tot", "log_NAA_rep", "log_FXSPR", "log_Y_FXSPR", "log_SSB_FXSPR")
  for(i in 1:length(rep_names)){
      ci.type <- "I"
      if(length(grep("log_", rep_names[i]))) ci.type <- "exp"
      x[[rep_names[i]]] <- list(est = std$est[[rep_names[i]]], se = std$se[[rep_names[i]]], 
        ci = get.ci(std$est[[rep_names[i]]], std$se[[rep_names[i]]], type = ci.type, alpha.ci = alphaCI, asap = FALSE))
  }
  # if("log_FMSY" %in% names(std$est)){
  # # if("log_FMSY" %in% rownames(std)){
  #   # inds$Fmsy <- which(rownames(std) == "log_FMSY")
  #   x$log_FMSY <- list(est = std$est$log_FMSY, se = std$se$log_FMSY, ci = get.ci(std$est$log_FMSY, std$se$log_FMSY, type = "exp"))
  #   # x$log_FMSY <- cbind(std[inds$Fmsy,1:2], get.ci(std[inds$Fmsy,1:2], alpha=alphaCI))
  #   # colnames(x$log_FMSY) <- c("log_est","log_se","est","lo","hi")
  # }
  #if("log_SSB_MSY" %in% names(std$est)){
  # # if("log_SSB_MSY" %in% rownames(std)){
  #  inds$SSBmsy <- matrix(which(rownames(std.summ) == "log_SSB_MSY"), ncol = NCOL(std))
  #   x$log_SSB_MSY <- list(est = std$est$log_SSB_MSY, se = std$se$log_SSB_MSY, ci = get.ci(std$est$log_FMSY, std$se$log_FMSY, type = "exp"))
  #   # x$log_SSB_MSY <- cbind(std[inds$SSBmsy[,all_stocks],1:2], get.ci(std[inds$SSBmsy[,all_stocks],1:2], alpha=alphaCI))
  #   # colnames(x$log_SSB_MSY) <- c("log_est","log_se","est","lo","hi")
  # }
  # if("log_MSY" %in% names(std$est)){
  # # if("log_MSY" %in% rownames(std)){
  #   inds$msy <- array(which(rownames(std) == "log_MSY"), dim = dim(mod$rep$log_MSY))
  #   x$log_MSY <- cbind(std[inds$msy[all_catch,all_stocks,],1:2], get.ci(std[inds$msy[all_catch,all_stocks,],1:2], alpha=alphaCI))
  #   colnames(x$log_MSY) <- c("log_est","log_se","est","lo","hi")
  # }
  # inds$catch <- which(rownames(std) == "log_pred_catch")
  # x$log_pred_catch <- cbind(std[inds$catch,1:2], get.ci(std[inds$catch,1:2], alpha=alphaCI))
  # colnames(x$log_pred_catch) <- c("log_est","log_se","est","lo","hi")

  # x$log_SSB <- cbind(std[inds$ssb,1:2], get.ci(std[inds$ssb,1:2], alpha=alphaCI))
  # colnames(x$log_SSB) <- c("log_est","log_se","est","lo","hi")
  # x$SSB_CV <- std$se$log_SSB_all
  # x$log_F <- cbind(std[inds$full.f,1:2], get.ci(std[inds$full.f,1:2], alpha=alphaCI))
  # colnames(x$log_F) <- c("log_est","log_se","est","lo","hi")
  # # x$F_CV <- std[inds$full.f,2]
  # x$log_NAA <- array(std[inds$naa,1], dim = dim(mod$rep$NAA))
  #x$log_NAA <- std$est$log_NAA_rep
  #x$NAA_CV <- std$se$log_NAA_rep #array(std[inds$naa,2], dim = dim(mod$rep$NAA))
  # x$log_NAA_lo <- exp(x$log_NAA - qnorm(1-alphaCI/2)*x$NAA_CV)
  # x$log_NAA_hi <- exp(x$log_NAA + qnorm(1-alphaCI/2)*x$NAA_CV)

  x$percentSPR <- mod$env$data$percentSPR
  # x$log_Y_FXSPR <- cbind(std[inds$Y.t[,all_catch],1:2], get.ci(std[inds$Y.t[,all_catch],1:2], alpha=alphaCI))
  # colnames(x$log_Y_FXSPR) <- c("log_est","log_se","est","lo","hi")
  # x$log_FXSPR <- cbind(std[inds$F.t,1:2], get.ci(std[inds$F.t,1:2], alpha=alphaCI))
  # colnames(x$log_FXSPR) <- c("log_est","log_se","est","lo","hi")
  # x$log_SSB_FXSPR <- cbind(std[inds$SSB.t[,all_stocks],1:2], get.ci(std[inds$SSB.t[,all_stocks],1:2], alpha=alphaCI))
  # colnames(x$log_SSB_FXSPR) <- c("log_est","log_se","est","lo","hi")
  # x$Y_FXSPR_CV <- std[inds$Y.t,2]
  # x$FXSPR_CV <- std[inds$F.t,2]
  # x$SSB_FXSPR_CV <- std[inds$SSB.t,2]

  cov <- mod$sdrep$cov
  x$log_rel_ssb_F_cov <- lapply(1:n_years, function(x){
    K <- cbind(c(1,-1,0,0),c(0,0,1,-1))
    ind <- c(inds$ssb[x],inds$SSB.t[x],inds$full.f[x],inds$F.t[x])
    tcov <- cov[ind,ind]
    ests <- std.summ[ind,1]
    diff <- t(K) %*% ests
    return(list(diff, t(K) %*% tcov %*% K))
  })
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
  return(x)
}

