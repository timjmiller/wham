#' Perform self-test simulation and estimation of a fitted WHAM model
#'
#' Generate simulated data from fitted model and refit the model to the simulated data
#'
#' @param fit_RDS (required) location of RDS file with fitted WHAM model.
#' @param n the number of simulated data sets to create and fit. Default = 10.
#' @param seeds (optional) vector of seeds to use to generate each simulated data set.
#' @param which_seeds (optional) which \code{seeds} to use for simualted data set. Useful if doing self test in stages.
#' @param conditional T/F whether to fix rather than simulate estimated random effects. Deafult = TRUE.
#' @param map_change (optional) list of input$map elements for altering mapping assumptions of fitted model.
#' @param res_dir directory where to save individual files for each self test fit. Useful if doing self test in stages. If not provided, no files will be saved.
#' @param do_parallel T/F whether to do self-test fits in parallel. Requires snowfall and parallel packages to be installed. Default = TRUE.
#' @param n_cores (optional) the number of cores to use for parallel fitting.
#' @param wham_location (optional) location of WHAM package. Useful if not using the WHAM installation in the standard library location.
#' @param test_dir (optional) directory for package repository. To be used when the function is being called during package testing rather than an installed version of WHAM.
#' @param save_inputs T/F whether to save the simulated inputs in res_dir. Default = FALSE.
#'
#' @return a list of two elements. First is the results which is a list (length = n) of lists with 6 elements: minimized negative log-likelihood, MLEs, gradient, SSB, F, abundance at age.
#'  Second element is the vector of seeds used for self-test simulations.
#'
#' @seealso \code{\link{fit_wham}}
#' @export
#'
self_test <- function(fit_RDS = NULL, n = 10, seeds = NULL, which_seeds = NULL, conditional = TRUE, map_change = NULL, do_parallel = TRUE, n_cores  = NULL, 
  res_dir = NULL, wham_location = NULL, test_dir = NULL, save_inputs = FALSE){
  
  if(is.null(fit_RDS)) stop("Provide fit_RDS, an RDS file name for a fitted WHAM model.")
  if(!is.null(res_dir)) {
    cat("res_dir is provided, so jitter files will be saved to ", res_dir, ". \n")
  }
  #if(is.null(wham_location)) wham_location <- system.file(package="wham")
  if(is.null(which_seeds)) which_seeds <- 1:n
  is_snowfall <- nchar(system.file(package="snowfall"))>0
  is_parallel <- nchar(system.file(package="parallel"))>0
  
  mod <- readRDS(fit_RDS)
  sim_type <- "sim"
  if(conditional) sim_type <- "cond_sim"

  if(is.null(seeds)) {
    set.seed(8675309)
    seeds <- sample(0:1e9, n)
  }
  if(!all(which_seeds %in% 1:length(seeds))) stop("some of which_seeds are outside 1 to n.")
  
  if(do_parallel){
    if(is_snowfall & is_parallel){
      if(is.null(n_cores)) n_cores <- parallel::detectCores()/2
      if(!is.null(res_dir)) snowfall::sfInit(parallel=TRUE, cpus=n_cores, slaveOutfile=file.path(res_dir,"self_test_log.txt"))
      snowfall::sfExportAll()
      sim_res <- snowfall::sfLapply(which_seeds, function(i){
        if(is.null(test_dir)) library(wham, lib.loc = wham_location)
        else pkgload::load_all(test_dir)
        fit <- readRDS(fit_RDS)
        sim_input <- fit$input
        sim_input$par <- fit$parList
        sim_input$data$do_SPR_BRPs[] <- 0
        sim_input$random <- NULL # so fit_wham doesn't try to do inner optimization
        if(conditional){
          temp <- grep(glob2rx("do_simulate_*_re"), names(sim_input$data), value = TRUE)
          # temp <- c("do_simulate_Ecov_re", "do_simulate_L_re", "do_simulate_M_re", "do_simulate_mu_prior_re", "do_simulate_mu_re", "do_simulate_N_re",
          #           "do_simulate_q_prior_re", "do_simulate_q_re", "do_simulate_sel_re")
          sim_input$data[temp] <- rep(list(0), length(temp))
          #sim_input$data[temp] <- lapply(temp, function(x) sim_input$data[[x]][] <- 0)
        }
        if(!is.null(map_change)) sim_input$map[names(map_change)] <- map_change
        sim_mod <- fit_wham(sim_input, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)
        set.seed(seeds[i])
        sim_input$data <- sim_mod$simulate(complete=TRUE)
        sim_input$random <- fit$input$random #set random correctly for estimation
        x <- try(fit_wham(sim_input, do.sdrep = FALSE, do.retro = FALSE, do.osa = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE))
        out <- list(obj = NA, 
          par = rep(NA,length(sim_mod$par)), 
          grad = rep(NA, length(sim_mod$par)), 
          SSB = matrix(NA,NROW(sim_mod$rep$SSB),NCOL(sim_mod$rep$SSB)), 
          F = rep(NA,length(sim_mod$rep$log_F_tot)), 
          NAA = array(NA, dim = dim(sim_mod$rep$NAA)))
        if(!is.null(x$opt)){
          out$obj <- x$opt$obj
          out$par <- x$opt$par
          out$grad <- x$final_gradient
          out$SSB <- x$rep$SSB
          out$F <- exp(x$rep$log_F_tot)
          out$NAA <- x$rep$NAA
        }
        out$seed <- seeds[i]
        if(!is.null(res_dir)){
          if(save_inputs) saveRDS(sim_input, file.path(res_dir, paste0(sim_type, "_input_", i, ".RDS")))
          saveRDS(out, file.path(res_dir, paste0(sim_type, "_", i, ".RDS")))
        }
        return(out)
      })
      snowfall::sfStop()
    } else stop("To do self test fits in parallel, install the snowfall and parallel packages. Otherwise, set do_parallel = FALSE.")
  } else{
    if(!(is_snowfall & is_parallel)) cat("If snowfall and parallel packages are installed, self test can be fit in parallel. \n")
    sim_res <- list()
    for(i in 1:length(which_seeds)){
      fit <- readRDS(fit_RDS)
      sim_input <- fit$input
      sim_input$par <- fit$parList
      sim_input$data$do_SPR_BRPs[] <- 0
      sim_input$random <- NULL # so fit_wham doesn't try to do inner optimization
      if(conditional){
          temp <- grep(glob2rx("do_simulate_*_re"), names(sim_input$data), value = TRUE)
        # temp <- c("do_simulate_Ecov_re", "do_simulate_L_re", "do_simulate_M_re", "do_simulate_mu_prior_re", "do_simulate_mu_re", "do_simulate_N_re",
        #           "do_simulate_q_prior_re", "do_simulate_q_re", "do_simulate_sel_re")
        sim_input$data[temp] <- rep(list(0), length(temp))
        # sim_input$data[temp] <- lapply(temp, function(x) sim_input$data[[x]][] <- 0)
      }
      if(!is.null(map_change)) sim_input$map[names(map_change)] <- map_change
      sim_mod <- fit_wham(sim_input, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)
      set.seed(seeds[i])
      sim_input$data <- sim_mod$simulate(complete=TRUE)
      sim_input$random <- fit$input$random #set random correctly for estimation
      x <- try(fit_wham(sim_input, do.sdrep = FALSE, do.retro = FALSE, do.osa = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE))
      out <- list(obj = NA, 
        par = rep(NA,length(sim_mod$par)), 
        grad = rep(NA, length(sim_mod$par)), 
        SSB = matrix(NA,NROW(sim_mod$rep$SSB),NCOL(sim_mod$rep$SSB)), 
        F = rep(NA,length(sim_mod$rep$log_F_tot)), 
        NAA = array(NA, dim = dim(sim_mod$rep$NAA)))
      if(!(is.character(x) | is.null(x$opt))){
        out$obj <- x$opt$obj
        out$par <- x$opt$par
        out$grad <- x$final_gradient
        out$SSB <- x$rep$SSB
        out$F <- exp(x$rep$log_F_tot)
        out$NAA <- x$rep$NAA
      }
      out$seed <- seeds[i]
      sim_res[[i]] <- out
      if(!is.null(res_dir)){
        if(save_inputs) saveRDS(sim_input, file.path(res_dir, paste0(sim_type, "_input_", i, ".RDS")))
        saveRDS(out, file.path(res_dir, paste0(sim_type, "_", i, ".RDS")))
      }
    }
  }
  return(list(self_test_results = sim_res, seeds = seeds))
}
