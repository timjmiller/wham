set_LW = function(input, LW)
{
  data = input$data
  par = input$par
  map = input$map

  n_re_par = 2 # number of parameters for RE

  # LW default options:
  n_par_def = 2 # 2 parameters: a and b: W = a*L^b
  data$LW_re_model = rep(1, times = n_par_def) # default = no RE / 'none'
  data$n_LW_par = n_par_def # 
  data$LW_est <- rep(0, times = n_par_def) # default = don't estimate M
  LW_re_ini = array(0, dim = c(data$n_years_model, data$n_ages, n_par_def))
  LW_ini = c(log(5e-06), log(3)) 

  # prepare LW options:
  if(!is.null(LW)){

    data$weight_model = 2 # use LW

    if(!is.null(LW$re)){
      
      if(length(LW$re) != data$n_LW_par) stop("Number of 're' must be equal to the number of L-W parameters.")
      for(k in 1:data$n_LW_par) {
          if(!(LW$re[k] %in% c("none","iid_y","iid_c","ar1_y","ar1_c"))) stop(paste0("LW$re[", k, "] must be one of the following: 'none','iid_y','iid_c','ar1_y','ar1_c'"))
          data$LW_re_model[k] <- match(LW$re[k], c("none","iid_y","iid_c","ar1_y","ar1_c")) # Respect this order to create array later
      }
      
    }

    if(!is.null(LW$init_vals)){
      if(length(LW$init_vals) != data$n_LW_par) stop("Length(LW$init_vals) must be 2.")
      LW_ini <- log(LW$init_vals)
    }
  
    if(!is.null(LW$est_pars)){
        if(length(LW$est_pars) > data$n_LW_par) stop("LW$est_pars should contain values equal or less than 2.")
        data$LW_est[LW$est_pars] = 1
    }

  }
  data$n_LW_est <- sum(data$LW_est)

  # LW pars --------------------------
  
  par$LW_a = as.matrix(LW_ini)
  par$LW_re = LW_re_ini
  par$LW_repars = matrix(0, ncol = n_re_par, nrow = data$n_LW_par)
  par$LW_repars[,1] = log(0.1) # start sigma at 0.1, rho at 0


  # --------------------------------
  # Prepare data for LW:
  tmp1 <- par$LW_a

  tmp1[data$LW_est==0] = NA
  ind.notNA <- which(!is.na(tmp1))
  tmp1[ind.notNA] <- 1:length(ind.notNA)
  map$LW_a <- factor(tmp1)

  # RE info:
  map$LW_re = NULL
  map$LW_repars = NULL
  tmp.g.repars <- matrix(NA, ncol = n_re_par, nrow = data$n_LW_par)
  max_val_par = 0
  active_sum = -1
  this_max = 0
  for(i in 1:data$n_LW_par) {

    # LW_re: "none","iid","ar1_y"
    tmp1 <- matrix(NA, nrow = data$n_years_model, ncol = data$n_ages)
    if(data$LW_re_model[i] %in% c(2,4)){ # iid ar1- only y
      tmp1[] = rep(1:nrow(tmp1), times = ncol(tmp1))  # all y estimated
      active_sum = active_sum + 1
      max_val_par = max_val_par + this_max * min(1, active_sum)
      this_max = max(tmp1, na.rm = TRUE)
    }
    if(data$LW_re_model[i] %in% c(3,5)){ # iid ar1 - only c
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

    map$LW_re <- c(map$LW_re, as.vector(tmp1) + max_val_par)

    # K_repars: sigma_M, rho_M_y
    if(data$LW_re_model[i] == 1) tmp.g.repars[i,] <- rep(NA,n_re_par) # no RE pars to estimate
    if(data$LW_re_model[i] == 2) tmp.g.repars[i,] <- c(1,NA) # estimate sigma y
    if(data$LW_re_model[i] == 3) tmp.g.repars[i,] <- c(1,NA) # estimate sigma c
    if(data$LW_re_model[i] == 4) tmp.g.repars[i,] <- c(1,2) # ar1_y: estimate sigma y, rho_y
    if(data$LW_re_model[i] == 5) tmp.g.repars[i,] <- c(1,2) # ar1_c: estimate sigma c, rho_c

  }

  map$LW_re = as.factor(map$LW_re)
  ind.notNA <- which(!is.na(tmp.g.repars))
  tmp.g.repars[ind.notNA] <- 1:length(ind.notNA)
  map$LW_repars = factor(tmp.g.repars)

  # End section

  input$data = data
  input$par = par
  input$map = map
  return(input)

}