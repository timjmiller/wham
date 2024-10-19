#' Jitter starting values of a fitted WHAM model
#'
#' Refits a previously fitted WHAM model a number of times with different initial parameter values
#'
#' Because the likelihood surface of a specific model may not be monotone, optimization may result in a local maximum which can be dependent on the parameter
#' values used to start the optimization. It is considered good practice to check whether the optimization is consistent for a range of starting values.
#' There does not appear to be an approach that is both best and computationally efficient to do this check so, this function allows the user to provide the 
#' alternative starting values for all fixed effects parameters to use in the jittering. When \code{initial_vals} is not provided, a default method is used
#' where n_jitter random draws from a MVN distribution with mean equal to the MLE from the fitted model and covariance matrix equal to that estimated for the
#' the MLEs.
#'
#' @param fit_RDS (required) location of RDS file with fitted WHAM model.
#' @param n_jitter the number of different starting locations on the likelihood surface to refit the model from. Default = 10.
#' @param initial_vals (optional) matrix (n_jitter x n parameters) of starting values to use for jittering.
#' @param which_rows (optional) which rows of intial_vals to use for jittering. Useful if doing jitter fits in stages.
#' @param do_parallel T/F whether to do jitter fits in parallel. Requires snowfall and parallel packages to be installed. Default = TRUE.
#' @param n_cores (optional) the number of cores to use for parallel fitting.
#' @param res_dir directory where to save individual files for each jitter refit. Useful if jittering for some reason crashes R. If not provided, no files will be saved.
#' @param wham_location (optional) location of WHAM package. Useful if not using the WHAM installation in the standard library location.
#' @param test_dir (optional) directory for package repository. To be used when the function is being called during package testing rather than an installed version of WHAM.
#'
#' @return a list of two elements. First is the results which is a list (length = n_jitter) of lists with 3 elements: minimized negative log-likelihood, jitter MLEs, and gradient.
#'  Second element is the matrix of initial values for the jitters.
#'
#' @seealso \code{\link{fit_wham}}
#' @export
#'
jitter_wham <- function(fit_RDS = NULL, n_jitter = 10, initial_vals = NULL, which_rows = NULL, do_parallel = TRUE, n_cores  = NULL, res_dir = NULL, wham_location = NULL, test_dir = NULL){
  
  if(is.null(fit_RDS)) stop("Provide fit_RDS, an RDS file name for a fitted WHAM model.")
  if(!is.null(res_dir)) {
    cat("res_dir is provided, so jitter files will be saved to ", res_dir, ". \n")
  }
  #if(is.null(wham_location)) wham_location <- system.file(package="wham")
  if(is.null(which_rows)) which_rows <- 1:n_jitter
  is_snowfall <- nchar(system.file(package="snowfall"))>0
  is_parallel <- nchar(system.file(package="parallel"))>0
  
  mod <- readRDS(fit_RDS)

  if(is.null(initial_vals)) {
    if(is.null(mod$sdrep)) stop("If initial_vals are not provided, the model provided by fit_RDS must have do_sdreport() completed.")
    cov <- mod$sdrep$cov.fixed
    chol.L <- t(chol(cov))
    set.seed(8675309)
    initial_vals <- t(sapply(1:n_jitter, function(x) mod$opt$par + chol.L %*% cbind(rnorm(n= NCOL(cov)))))
  }
  if(!all(which_rows %in% 1:NROW(initial_vals))) stop("some of which_rows are outside the rows of initial_vals.")
  if(!is.null(res_dir)) saveRDS(initial_vals, file.path(res_dir, "initial_values.RDS"))
  if(do_parallel){
    if(is_snowfall & is_parallel){
      if(is.null(n_cores)) n_cores = parallel::detectCores()/2
      if(!is.null(res_dir)) snowfall::sfInit(parallel=TRUE, cpus=n_cores, slaveOutfile=file.path(res_dir,"jitter_log.txt"))
      snowfall::sfExportAll()
      jit_res <- snowfall::sfLapply(which_rows, function(row_i){
        if(is.null(test_dir)) library(wham, lib.loc = wham_location)
        else pkgload::load_all(test_dir)
        jit_mod <- readRDS(fit_RDS)
        jit_mod$env$data$do_SPR_BRPs[] <-  0
        jit_mod$env$data$do_MSY_BRPs[] <- 0
        jit_mod$env$inner.control$trace <- FALSE #to silence inner optimization in log file
        jit_mod$par[] <- initial_vals[row_i,]
        # x <- try(nlminb(initial_vals[row_i,], jit_mod$fn, jit_mod$gr))
        x <- try(fit_tmb(jit_mod, n.newton = 0, do.sdrep = FALSE))
        # snowfall::sfCat(names(x))
        # print(x)
    		out <- list(obj = NA, 
          par = rep(NA,length(jit_mod$par)), 
          grad = rep(NA,length(jit_mod$par)))
        out$initial_vals_row <- row_i
        out$initial_vals <- initial_vals[row_i,]
        out$rep <- NULL
        out$parList <- NULL
        out$opt <- NULL
        out$last.par.best <- rep(NA, length(jit_mod$env$last.par.best))
    		# if(!(is.character(x) | is.null(x))){
        if(!(is.character(x$opt) | is.null(x$opt))){
          out$opt <- x$opt
          out$rep <- x$rep
          out$parList <- x$parList
          out$last.par.best <- x$env$last.par.best
          out$obj <- x$opt$obj
          out$par <- x$opt$par
          out$grad <- x$gr(x$opt$par)
    		}
        if(!is.null(res_dir)){
      		saveRDS(out, file.path(res_dir, paste0("jitter_sim_", row_i, ".RDS")))
        }
    		return(out)
      })
      snowfall::sfStop()
    } else stop("To do jitter fits in parallel, install the snowfall and parallel packages. Otherwise, set do_parallel = FALSE.")
  } else{
    if(!(is_snowfall & is_parallel)) cat("If snowfall and parallel packages are installed, jitters can be fit in parallel. \n")
    jit_res <- list()
    for(i in 1:length(which_rows)){
      jit_mod <- readRDS(fit_RDS)
      jit_mod$env$data$do_SPR_BRPs[] <- 0
      jit_mod$par[] <- initial_vals[which_rows[i],]
      # jit_fit <- try(nlminb(initial_vals[row_i,], jit_mod$fn, jit_mod$gr))
      x <- try(fit_tmb(jit_mod, n.newton = 0, do.sdrep = FALSE))
      jit_res[[i]] <- list(obj = NA, par = rep(NA,length(jit_mod$par)), grad = rep(NA,length(jit_mod$par)))
      jit_res[[i]]$initial_vals_row <- row_i
      jit_res[[i]]$initial_vals <- initial_vals[row_i,]
      jit_res[[i]]$rep <- NULL
      jit_res[[i]]$parList <- NULL
      jit_res[[i]]$opt <- NULL
      jit_res[[i]]$last.par.best <- rep(NA, length(jit_mod$env$last.par.best))
      # if(!(is.character(x) | is.null(x))){
      if(!(is.character(x$opt) | is.null(x$opt))){
        jit_res[[i]]$opt <- x$opt
        jit_res[[i]]$rep <- x$rep
        jit_res[[i]]$parList <- x$parList
        jit_res[[i]]$last.par.best <- x$env$last.par.best
        jit_res[[i]]$obj <- x$opt$obj
        jit_res[[i]]$par <- x$opt$par
        jit_res[[i]]$grad <- x$gr(x$opt$par)
        # jit_res[[i]]$grad <- x$gr(x$par)
      }
      if(!is.null(res_dir)){
        saveRDS(jit_res[[i]], file.path(res_dir, paste0("jitter_sim_", which_rows[i], ".RDS")))
      }
    }
  }
  return(list(jitter_results = jit_res, initial_values = initial_vals))
}
