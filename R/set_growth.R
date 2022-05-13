set_growth = function(input, growth)
{
  data = input$data
  par = input$par
  map = input$map
    
  n_re_par = 2 # number of parameters for RE

  if(is.null(input$asap3)) asap3 = NULL
  else {
    stop('Growth module does not work with ASAP3 input')
  }
  # growth default options:
  n_par_def = 5 # 5 parameters: K, Linf, t0, CV1, CV2 
  data$growth_model = 1 # 1: vB-classic, 2: use LAA as input
  data$growth_re_model = rep(1, times = n_par_def) # default = no RE / 'none'
  data$growth_est <- rep(0, times = n_par_def) # default = don't estimate M
  data$n_growth_par = n_par_def # 5 parameters to estimate: K, Linf, t0, CV1, CV2, in that order
  growth_re_ini = array(0, dim = c(data$n_years_model, data$n_ages, n_par_def))
  growth_ini = c(log(0.2), log(100), log(0.01*-1 + 1), log(0.1), log(0.1)) 
  
  # prepare growth options:
  if(!is.null(growth)){

    n_par = c(5, 1) # for model 1 and 2, respectively

    if(!is.null(growth$model)){ # growth model to be used
      if(!(growth$model %in% c("vB_classic", "input_LAA"))) stop("growth$model must be 'vB_classic' or 'input_LAA'")
      data$growth_model <- match(growth$model, c("vB_classic", "input_LAA")) # 1
      if(growth$model %in% c("vB_classic", "input_LAA")){
        data$n_growth_par = n_par[data$growth_model] # 5 parameters to estimate
        data$growth_est = rep(0, times = n_par[data$growth_model]) # estimate?
        data$growth_re_model = rep(1, times = n_par[data$growth_model]) # default = no RE / 'none'
        growth_re_ini = array(0, dim = c(data$n_years_model, data$n_ages, n_par[data$growth_model]))
      }
    }

    if(is.null(growth$input_LAA) & data$growth_model == 2) stop("'input_LAA' must be provided when growth model is 2.")

    if(!is.null(growth$re)){
      if(length(growth$re) != n_par[data$growth_model]) stop("Length(growth$re) must be equal to the number of parameters of the chosen growth model (equal to 1 for 'input_LAA').")
      
      if(data$growth_model == 1) {
        if(!(growth$re[1] %in% c("none","iid_y","iid_c","ar1_y","ar1_c"))) stop("growth$re[1] (k) must be one of the following: 'none','iid_y','iid_c','ar1_y','ar1_c'")
        if(!(growth$re[2] %in% c("none","iid_y","iid_c","ar1_y","ar1_c"))) stop("growth$re[2] (Linf) must be one of the following: 'none','iid_y','iid_c','ar1_y','ar1_c'")
        if(!(growth$re[3] %in% c("none","iid_y","iid_c","ar1_y","ar1_c"))) stop("growth$re[3] (t0) must be one of the following: 'none','iid_y','iid_c','ar1_y','ar1_c'")
        if(!(growth$re[4] %in% c("none","iid_y","iid_c","ar1_y","ar1_c"))) stop("growth$re[4] (CV1) must be one of the following: 'none','iid_y','iid_c','ar1_y','ar1_c'")
        if(!(growth$re[5] %in% c("none","iid_y","iid_c","ar1_y","ar1_c"))) stop("growth$re[5] (CV2) must be one of the following: 'none','iid_y','iid_c','ar1_y','ar1_c'")
        data$growth_re_model[1] <- match(growth$re[1], c("none","iid_y","iid_c","ar1_y","ar1_c")) # Respect this order to create array later
        data$growth_re_model[2] <- match(growth$re[2], c("none","iid_y","iid_c","ar1_y","ar1_c")) # Respect this order to create array later
        data$growth_re_model[3] <- match(growth$re[3], c("none","iid_y","iid_c","ar1_y","ar1_c")) # Respect this order to create array later
        data$growth_re_model[4] <- match(growth$re[4], c("none","iid_y","iid_c","ar1_y","ar1_c")) # Respect this order to create array later
        data$growth_re_model[5] <- match(growth$re[5], c("none","iid_y","iid_c","ar1_y","ar1_c")) # Respect this order to create array later
      }
      if(data$growth_model == 2) {
        if(!identical(dim(growth$input_LAA), c(data$n_years_model, data$n_ages))) stop("Dimensions of 'input_LAA' must be the number of years and ages, respectively.")
        if(!(growth$re[1] %in% c("none","iid_y","iid_c","ar1_y","ar1_c"))) stop("growth$re[1] must be one of the following: 'none','iid_y','iid_c','ar1_y','ar1_c'")
        data$growth_re_model[1] <- match(growth$re[1], c("none","iid_y","iid_c","ar1_y","ar1_c")) # Respect this order to create array later
        growth_ini = growth$input_LAA # LAA
      }
    }

    if(!is.null(growth$init_vals) & data$growth_model == 2) stop("'init_vals' must be omitted when growth model is 'input_LAA'.")
    if(!is.null(growth$est_pars) & data$growth_model == 2) stop("'est_pars' must be omitted when growth model is 'input_LAA'.")

    if(!is.null(growth$init_vals)){
      if(length(growth$init_vals) != data$n_growth_par) stop("Length(growth$init_vals) must be equal to the number of parameters of the chosen growth model.")
      growth_ini <- c(log(growth$init_vals[1:2]), log(growth$init_vals[3]*-1 + 1), log(growth$init_vals[4:5]))
    }
	
  	if(!is.null(growth$est_pars)){
        if(length(growth$est_pars) > data$n_growth_par) stop("Should be equal or less than the number of parameters of the chosen growth model.")
        data$growth_est[growth$est_pars] = 1
  	}

  }
  data$n_growth_est <- sum(data$growth_est)

  # growth pars --------------------------
  
  par$growth_a = as.matrix(growth_ini)
  par$growth_re = growth_re_ini
  par$growth_repars = matrix(0, ncol = n_re_par, nrow = data$n_growth_par)
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
  tmp.g.repars <- matrix(NA, ncol = n_re_par, nrow = data$n_growth_par)
  max_val_par = 0
  active_sum = -1
  this_max = 0
  for(i in 1:data$n_growth_par) {

    # growth_re: "none","iid","ar1_y"
    tmp1 <- matrix(NA, nrow = data$n_years_model, ncol = data$n_ages)
    if(data$growth_re_model[i] %in% c(2,4)){ # iid ar1- only y
      tmp1[] = rep(1:nrow(tmp1), times = ncol(tmp1))  # all y estimated
      active_sum = active_sum + 1
      max_val_par = max_val_par + this_max * min(1, active_sum)
      this_max = max(tmp1, na.rm = TRUE)
    }
    if(data$growth_re_model[i] %in% c(3,5)){ # iid ar1 - only c
      loop_row = rep(0, times = ncol(tmp1) + nrow(tmp1) - 1)
      loop_row[1:ncol(tmp1)] = (ncol(tmp1) - 1):0
      loop_col = rep(ncol(tmp1) - 1, times = ncol(tmp1) + nrow(tmp1) - 1)
      loop_col[(length(loop_col) - ncol(tmp1) + 1):length(loop_col)] = (ncol(tmp1) - 1):0

      for(j in seq_along(loop_col)) {
        tmp1[(loop_row[j]:loop_col[j])*(nrow(tmp1) + 1) + (j - ncol(tmp1) + 1)] <- j
      }

      active_sum = active_sum + 1
      max_val_par = max_val_par + this_max * min(1, active_sum)
      this_max = max(tmp1, na.rm = TRUE)
    }

    map$growth_re <- c(map$growth_re, as.vector(tmp1) + max_val_par)

	  # K_repars: sigma_M, rho_M_y
	  if(data$growth_re_model[i] == 1) tmp.g.repars[i,] <- rep(NA,n_re_par) # no RE pars to estimate
	  if(data$growth_re_model[i] == 2) tmp.g.repars[i,] <- c(1,NA) # estimate sigma y
	  if(data$growth_re_model[i] == 3) tmp.g.repars[i,] <- c(1,NA) # estimate sigma c
    if(data$growth_re_model[i] == 4) tmp.g.repars[i,] <- c(1,2) # ar1_y: estimate sigma y, rho_y
    if(data$growth_re_model[i] == 5) tmp.g.repars[i,] <- c(1,2) # ar1_c: estimate sigma c, rho_c

  }

  map$growth_re = as.factor(map$growth_re)
  ind.notNA <- which(!is.na(tmp.g.repars))
  tmp.g.repars[ind.notNA] <- 1:length(ind.notNA)
  map$growth_repars = factor(tmp.g.repars)

  
  # End section

  input$data = data
  input$par = par
  input$map = map
  return(input)

}