set_growth = function(input, growth)
{
  data = input$data
  par = input$par
  map = input$map
  
  n_par = 5 # 5 parameters: K, Linf, t0, CV1, CV2
  n_re_par = 2  
	
  data$n_growth_a = data$n_ages # all ages but only used for K so far
  data$growth_re_model = rep(1, times = n_par) # constant all
  data$n_growth = n_par # 5 parameters to estimate: K, Linf, t0, CV1, CV2, in that order
  data$growth_model = 1 # 1: vB-classic, 2: vB-K_age
  data$growth_est <- rep(0, times = n_par) # default = don't estimate growth parameters
  
  growth_re_ini = matrix(NA, nrow = data$n_years_model, ncol = data$n_growth)
  
  if(is.null(input$asap3)) asap3 = NULL
  else {
    stop('Growth module does not work with ASAP3 input')
  }
  # growth default options:
  if(is.null(growth)) {
    data$growth_model = 1
    data$growth_re_model = rep(1, times = n_par) # default = no RE / 'none'
    data$growth_est <- rep(0, times = n_par) # default = don't estimate M
  }

  # prepare growth options:
  if(!is.null(growth)){

    if(!is.null(growth$model)){ # growth model to be used
      if(!(growth$model %in% c("vB_classic"))) stop("growth$model must be 'vB_classic'")
      data$growth_model <- match(growth$model, c("vB_classic")) # 1
      if(growth$model %in% c("vB_classic")){
        data$n_growth = n_par # 5 parameters to estimate
        data$growth_est = rep(0, times = n_par) # estimate?
        data$growth_ini = c(log(0.2), log(100), 0.01, log(0.1), log(0.1)) # taken from growth argument?
      }
    }

    if(!is.null(growth$re)){
      if(length(growth$re) != n_par) stop("Length(growth$re) must be equal to the number of parameters of the chosen growth model.")
      if(!(growth$re[1] %in% c("none","iid","ar1_y"))) stop("growth$re[1] (k) must be one of the following: 'none','iid','ar1_y'")
      if(!(growth$re[2] %in% c("none","iid","ar1_y"))) stop("growth$re[2] (Linf) must be one of the following: 'none','iid','ar1_y'")
      if(!(growth$re[3] %in% c("none","iid","ar1_y"))) stop("growth$re[3] (t0) must be one of the following: 'none','iid','ar1_y'")
      if(!(growth$re[4] %in% c("none","iid","ar1_y"))) stop("growth$re[4] (CV1) must be one of the following: 'none','iid','ar1_y'")
      if(!(growth$re[5] %in% c("none","iid","ar1_y"))) stop("growth$re[5] (CV2) must be one of the following: 'none','iid','ar1_y'")
      data$growth_re_model[1] <- match(growth$re[1], c("none","iid","ar1_y")) # Respect this order to create array later
      data$growth_re_model[2] <- match(growth$re[2], c("none","iid","ar1_y")) # Respect this order to create array later
      data$growth_re_model[3] <- match(growth$re[3], c("none","iid","ar1_y")) # Respect this order to create array later
    }

    if(!is.null(growth$init_vals)){
      if(length(growth$init_vals) != data$n_growth) stop("Length(growth$init_vals) must be equal to the number of parameters of the chosen growth model.")
      growth_ini <- c(log(growth$init_vals[1:2]), growth$init_vals[3], log(growth$init_vals[4:5]))
	  growth_re_ini[] <- 0#  initialize yearly deviations at 0
    }
	
	if(!is.null(growth$est_pars)){
      if(length(growth$est_pars) > data$n_growth) stop("Should be equal or less than the number of parameters of the chosen growth model.")
      data$growth_est[growth$est_pars] = 1
	  growth_re_ini[] <- 0#  initialize yearly deviations at 0
	}

  }
  data$n_growth_est <- sum(data$growth_est)

  # growth pars --------------------------
  
  par$growth_a = growth_ini
  par$growth_re = growth_re_ini
  par$growth_repars = matrix(0, ncol = n_re_par, nrow = data$n_growth)
  par$growth_repars[,1] = log(0.1) # start sigma at 0.1, rho at 0


  # --------------------------------
  # Prepare data for growth:
  tmp1 <- par$growth_a

  tmp1[data$growth_est==0] = NA
  ind.notNA <- which(!is.na(tmp1))
  tmp1[ind.notNA] <- 1:length(ind.notNA)
  map$growth_a <- factor(tmp1)

  # RE info:
  map$growth_re = NULL
  map$growth_repars = NULL
  for(i in 1:data$n_growth) {

	  # growth_re: "none","iid","ar1_y"
	  tmp1 <- rep(NA, times = data$n_years_model)
	  if(data$growth_re_model[i] == 1) tmp1[] = NA # no RE (either estimate RE for all ages or none at all)
	  if(data$growth_re_model[i] %in% c(2,3)){ # iid - only y
		  tmp1[] = (1:length(tmp1)) + data$n_years_model*(i-1) # all y estimated
	  }
	  map$growth_re <- c(map$growth_re, tmp1)

	  # K_repars: sigma_M, rho_M_y
	  if(data$growth_re_model[i] == 1) tmp1 <- rep(NA,n_re_par) # no RE pars to estimate
	  if(data$growth_re_model[i] == 2) tmp1 <- c(1,NA) # estimate sigma
	  if(data$growth_re_model[i] == 3) tmp1 <- c(1,2) # ar1_y: estimate sigma, rho_y

	  map$growth_repars = c(map$growth_repars, tmp1)

  }
  
  map$growth_re = as.factor(map$growth_re)
  nNoNA = which(!is.na(map$growth_repars))
  map$growth_repars[nNoNA] = 1:length(nNoNA)
  map$growth_repars = as.factor(map$growth_repars)
  
  # End section

  input$data = data
  input$par = par
  input$map = map
  return(input)

}