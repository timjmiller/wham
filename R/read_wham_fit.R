#' Read WHAM fit
#'
#' Gets output from a fit WHAM model for plotting.
#'
#' @param mod output from \code{\link{fit_wham}}
#'
#' @return a named list with the following elements:
#'   \describe{
#'     \item{\code{$years}}{numeric vector, model years only, e.g. 1972:2020}
#'     \item{\code{$years_full}}{numeric vector, model + proj years, e.g. 1972:2022}
#'     \item{\code{$selAA}}{list of length(n_selblocks), first the fleet blocks then indices, i.e. if 4 fleet blocks and 3 indices, selAA[[5]] is for index 1. Each element is a matrix, years (rows) x ages (cols), selectivity at age}
#'     \item{\code{$selblock_pointer_fleets}}{matrix, years x fleets, indices of selAA used by each fleet in each year}
#'     \item{\code{$selblock_pointer_indices}}{matrix, years x indices, indices of selAA used by each index in each year}
#'     \item{\code{$MAA}}{matrix, years x ages, natural mortality}
#'     \item{\code{$log_SSB}}{vector, years, log-scale spawning stock biomass}
#'     \item{\code{$SSB_CV}}{vector, years, CV of spawning stock biomass (i.e. 95% CI calculated as exp(log_SSB +/- 1.96*SSB_CV)}
#'     \item{\code{$log_F}}{vector, years, log-scale fully-selected F}
#'     \item{\code{$F_CV}}{vector, years, CV of full F}
#'     \item{\code{$log_NAA}}{matrix, years x ages, numbers at age}
#'     \item{\code{$NAA_CV}}{matrix, years x ages, CV of numbers at age}
#'     \item{\code{$percentSPR}}{scalar, X% SPR used to calculate reference points, default = 40}
#'     \item{\code{$log_Y_FXSPR}}{vector, years, log-scale yield at FXSPR}
#'     \item{\code{$log_FXSPR}}{vector, years, log-scale FXSPR}
#'     \item{\code{$log_SSB_FXSPR}}{vector, years, log-scale SSB at FXSPR}
#'     \item{\code{$Y_FXSPR_CV}}{vector, years, CV of yield at FXSPR}
#'     \item{\code{$FXSPR_CV}}{vector, years, CV of FXSPR}
#'     \item{\code{$SSB_FXSPR_CV}}{vector, years, CV of SSB at FXSPR}
#'     \item{\code{$log_rel_ssb_F_cov}}{list, length n_years, each element is a 2x2 covariance matrix with SSB/SSB_FXSPR first and F/F_FXSPR second}
#'   }
#'
#' @export
#'
#' @seealso \code{\link{fit_wham}}, \code{wham_plots_tables}
#'
read_wham_fit <- function(mod){
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
  std <- summary(mod$sdrep, "report")
  inds <- list(Y.t = which(rownames(std) == "log_Y_FXSPR"))
  inds$F.t <- which(rownames(std) == "log_FXSPR")
  inds$SSB.t <- which(rownames(std) == "log_SSB_FXSPR")
  inds$ssb <- which(rownames(std) == "log_SSB")
  inds$faa <- which(rownames(std) == "log_FAA_tot")
  log.faa <- matrix(std[inds$faa,1], n_years, n_ages)
  age.full.f <- apply(log.faa,1, function(x) max(which(x == max(x))))
  inds$full.f <- (age.full.f-1)*n_years + 1:n_years  + min(inds$faa) - 1 #cbind(1:n_years, age.full.f)
  inds$naa <- which(rownames(std) == "log_NAA_rep")

  x$log_SSB <- std[inds$ssb,1]
  x$SSB_CV <- std[inds$ssb,2]
  x$log_F <- std[inds$full.f,1]
  x$F_CV <- std[inds$full.f,2]
  x$log_NAA <- matrix(std[inds$naa,1], n_years, n_ages)
  x$NAA_CV <- matrix(std[inds$naa,2], n_years, n_ages)

  x$percentSPR <- mod$env$data$percentSPR
  x$log_Y_FXSPR <- std[inds$Y.t,1]
  x$log_FXSPR <- std[inds$F.t,1]
  x$log_SSB_FXSPR <- std[inds$SSB.t,1]
  x$Y_FXSPR_CV <- std[inds$Y.t,2] 
  x$FXSPR_CV <- std[inds$F.t,2]
  x$SSB_FXSPR_CV <- std[inds$SSB.t,2]

  cov <- mod$sdrep$cov
  x$log_rel_ssb_F_cov <- lapply(1:n_years, function(x){
    K <- cbind(c(1,-1,0,0),c(0,0,1,-1))
    ind <- c(inds$ssb[x],inds$SSB.t[x],inds$full.f[x],inds$F.t[x])
    tcov <- cov[ind,ind]
    return(t(K) %*% tcov %*% K)
  })

  return(x)
}
  
