set_LAA = function(input, LAA)
{
  data = input$data
  par = input$par
  map = input$map
    
  n_re_par = 2 # number of parameters for RE

  if(is.null(input$asap3)) asap3 = NULL
  else {
    stop('LAA does not work with ASAP3 input')
  }

  # LAA default options:
  LAA_re_ini = matrix(0, ncol = data$n_ages, nrow = data$n_years_model)
  data$LAA_re_model = 1
  data$LAA_est = rep(0, data$n_ages)
  LAA_ini = log( 100 + (3 - 100)*exp(-0.2*(1:data$n_ages - 1)) )
  if(!is.null(LAA)) {
    LAA_ini = LAA$LAA_vals
    if(!is.null(LAA$LAA_est)) data$LAA_est[LAA$LAA_est] = 1
    if(!is.null(LAA$re))  {
      if(!(LAA$re %in% c("none","iid_y","iid_c","ar1_y","ar1_c"))) stop("LAA$re must be one of the following: 'none','iid_y','iid_c','ar1_y','ar1_c'")
      data$LAA_re_model <- match(LAA$re, c("none","iid_y","iid_c","ar1_y","ar1_c")) # Respect this order to create array later
    }
  }

  # growth pars --------------------------
  
  par$LAA_a = log(LAA_ini)
  par$LAA_re = LAA_re_ini
  par$LAA_repars = matrix(0, ncol = n_re_par, nrow = 1)
  par$LAA_repars[,1] = log(0.1) # start sigma at 0.1, rho at 0


  # --------------------------------
  # Prepare data for growth:
  tmp1 <- par$LAA_a
  tmp1[data$LAA_est==0] = NA
  ind.notNA <- which(!is.na(tmp1))
  tmp1[ind.notNA] <- 1:length(ind.notNA)
  map$LAA_a <- factor(tmp1)

  # RE info:
  map$LAA_re = NULL
  map$LAA_repars = NULL
  tmp.g.repars <- matrix(NA, ncol = n_re_par, nrow = 1)
  max_val_par = 0
  active_sum = -1
  this_max = 0

    i = 1
    # growth_re: "none","iid","ar1_y"
    tmp1 <- matrix(NA, nrow = data$n_years_model, ncol = data$n_ages)
    if(data$LAA_re_model[i] %in% c(2,4)){ # iid ar1- only y
      tmp1[] = rep(1:nrow(tmp1), times = ncol(tmp1))  # all y estimated
      active_sum = active_sum + 1
      max_val_par = max_val_par + this_max * min(1, active_sum)
      this_max = max(tmp1, na.rm = TRUE)
    }
    if(data$LAA_re_model[i] %in% c(3,5)){ # iid ar1 - only c
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

    map$LAA_re <- c(map$LAA_re, as.vector(tmp1) + max_val_par)

	  # K_repars: sigma_M, rho_M_y
	  if(data$LAA_re_model[i] == 1) tmp.g.repars[i,] <- rep(NA,n_re_par) # no RE pars to estimate
	  if(data$LAA_re_model[i] == 2) tmp.g.repars[i,] <- c(1,NA) # estimate sigma y
	  if(data$LAA_re_model[i] == 3) tmp.g.repars[i,] <- c(1,NA) # estimate sigma c
    if(data$LAA_re_model[i] == 4) tmp.g.repars[i,] <- c(1,2) # ar1_y: estimate sigma y, rho_y
    if(data$LAA_re_model[i] == 5) tmp.g.repars[i,] <- c(1,2) # ar1_c: estimate sigma c, rho_c


  map$LAA_re = as.factor(map$LAA_re)
  ind.notNA <- which(!is.na(tmp.g.repars))
  tmp.g.repars[ind.notNA] <- 1:length(ind.notNA)
  map$LAA_repars = factor(tmp.g.repars)

  
  # End section

  input$data = data
  input$par = par
  input$map = map
  return(input)

}