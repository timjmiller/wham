#' Specify configuration for fully-selected fishing mortality
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param F_opts (optional) named list of initial values for annual fully-selected fishing mortality and configuration method for estimation.
#' 
#' \code{F_opts} specifies a few as well as the effect on the population. Environmental covariate data need not span
#' the same years as the fisheries data. It can be \code{NULL} if no environmental data are to be fit.
#' Otherwise, it must be a named list with the following components:
#'   \describe{
#'     \item{$F}{matrix (n_years x n_fleets) of (initial) values for fully-selected fishing morality.}
#'     \item{$F_config}{integer 1: (default) configure F parameters (on log scale) as an F in the initial year and then deviations from one year to the next,
#'        or 2: configure F parameters as (log) annual values.}
#'  }

set_F = function(input, F_opts = NULL)
{
  asap3 = input$asap3
  input$data$F_config <- 1 #log_F1, F_devs
  input$par$F_pars <- matrix(0,input$data$n_years_model, input$data$n_fleets)
  #F_devs = matrix(0, input$data$n_years_model-1, input$data$n_fleets)
	if(!is.null(asap3)) {
    k = 1
    for(i in 1:length(asap3))for(j in 1:length(asap3[[i]]$F1_ini)){
      input$par$F_pars[1,k] = log(asap3[[i]]$F1_ini[j]) # use F1_ini values from asap3 file  
      k = k + 1
    }
	} else {
  	input$par$F_pars[1,] = log(0.2) # old
  }
  if(!is.null(F_opts$F_config)) input$data$F_config <- F_opts$F_config
  if(!is.null(F_opts$F)) {
    if(input$data$F_config == 2) {
      input$par$F_pars[] = log(F_opts$F[])
    } else { # F_config = 1
      input$par$F_pars[1,] <- log(F_opts$F[1,])
      for(f in 1:data$n_fleets) input$par$F_pars[-1,f] <- diff(log(F_opts$F[,f]))
    }
  }
  return(input)
}