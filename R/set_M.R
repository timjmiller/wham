set_M = function(input, M)
{
  data = input$data
  par = input$par
  map = input$map
  data$n_M_a = data$n_ages
  data$M_model = 2
  data$use_b_prior = 0
  data$M_re_model = 1 # default = no RE / 'none'
  data$M_est <- rep(0, data$n_M_a) # default = don't estimate M
  M_first_est = NA  
  M_re_ini = matrix(NA, data$n_years_model, data$n_M_a)
  
  if(is.null(input$asap3)) asap3 = NULL
  else {
    asap3 = input$asap3
    M_a_ini <- log(asap3$M[1,])
    M_re_ini[] <- matrix(log(asap3$M), data$n_years_model, data$n_M_a) - matrix(M_a_ini, data$n_years_model, data$n_M_a, byrow=T)
  }

  # natural mortality options, default = use values from ASAP file, no estimation
  if(is.null(M))
  {
    data$M_model = 2
    data$use_b_prior = 0
    data$M_re_model = 1 # default = no RE / 'none'
    data$M_est <- rep(0, data$n_M_a) # default = don't estimate M
    M_first_est = NA
    # Use ASAP file M values
    if(!is.null(asap3)){
    }
    else {
      M_a_ini <- log(rep(0.2, data$n_ages))
      M_re_ini[] <- 0
    }
  }
  if(!is.null(M)){
    if(!is.null(M$model)){ # M model options
      if(!(M$model %in% c("constant","age-specific","weight-at-age"))) stop("M$model must be either 'constant', 'age-specific', or 'weight-at-age'")
      if(!is.null(M$re)) if(M$model == "age-specific" & M$re == "ar1_a") stop("Cannot estimate age-specific mean M and AR1 deviations M_a.
If you want an AR1 process on M-at-age, set M$model = 'constant' and M$re = 'ar1_a'.")
      data$M_model <- match(M$model, c("constant","age-specific","weight-at-age"))
      if(M$model %in% c("constant","weight-at-age")){
        data$n_M_a = 1
        data$M_est = 0
        if(!is.null(asap3)) M_a_ini = log(asap3$M[1,1])
        else M_a_ini = log(0.2)
        if(!is.null(asap3)) {
          if(is.null(M$initial_means) & length(unique(asap3$M[1,])) > 1) warning("Constant or weight-at-age M specified (so only 1 mean M parameter),
but first row of MAA matrix has > 1 unique value.
Initializing M at age-1 MAA values. To avoid this warning
without changing ASAP file, specify M$initial_means.")
        }
        if(!is.null(M$logb_prior)){
          if(!is.logical(M$logb_prior)) stop("M$logb_prior must be either TRUE or FALSE")
          if(M$logb_prior) data$use_b_prior = 1
        }
      }
    }
    if(!is.null(M$re)){
      if(length(M$re) != 1) stop("Length(M$re) must be 1")
      if(!(M$re %in% c("none","iid","ar1_a","ar1_y","2dar1"))) stop("M$re must be one of the following: 'none','iid','ar1_a','ar1_y','2dar1'")
      data$M_re_model <- match(M$re, c("none","iid","ar1_a","ar1_y","2dar1"))
    }
    if(!is.null(M$initial_means)){
      if(length(M$initial_means) != data$n_M_a) stop("Length(M$initial_means) must be # ages (if age-specific M) or 1 (if constant or weight-at-age M)")
      M_a_ini <- log(M$initial_means)
      # overwrite ASAP file M values
      M_re_ini[] <- 0# if estimating mean M for any ages, initialize yearly deviations at 0
    }
    if(!is.null(M$est_ages)){
      if(!all(M$est_ages %in% 1:data$n_M_a)) stop("All M$est_ages must be in 1:n.ages (if age-specific M) or 1 (if constant or weight-at-age M)")
      data$M_est[M$est_ages] = 1 # turn on estimation for specified M-at-age
      M_first_est <- M$est_ages[1]
      M_re_ini[] <- 0# if estimating mean M for any ages, initialize yearly deviations at 0
    }
  }
  data$n_M_est <- sum(data$M_est)

  # natural mortality pars
  par$M_a <- M_a_ini # deviations by age
  par$M_re <- M_re_ini # deviations from mean M_a on log-scale, PARAMETER_ARRAY
  par$M_repars <- rep(0, 3)
  par$M_repars[1] <- log(0.1) # start sigma at 0.1, rho at 0
  if(data$M_re_model == 3) par$M_repars[3] <- 0 # if ar1 over ages only, fix rho_y = 0
  if(data$M_re_model == 4) par$M_repars[2] <- 0 # if ar1 over years only, fix rho_a = 0
  # check if only 1 estimated mean M (e.g. because weight-at-age M or if all but 1 age is fixed), can't estimate rho_a
  # if(data$n_M_est < 2) par$M_repars[2] <- 0
  par$log_b = log(0.305)

  tmp <- par$M_a
  tmp[data$M_est==0] = NA
  ind.notNA <- which(!is.na(tmp))
  tmp[ind.notNA] <- 1:length(ind.notNA)
  map$M_a <- factor(tmp)
  if(data$M_model != 3) map$log_b = factor(rep(NA,length(par$log_b)))

  # M_re: "none","iid","ar1_a","ar1_y","2dar1"
  tmp <- par$M_re
  if(data$M_re_model == 1) tmp[] = NA # no RE (either estimate RE for all ages or none at all)
  if(data$M_re_model %in% c(2,5)){ # 2d ar1
    tmp[] = 1:(dim(tmp)[1]*dim(tmp)[2]) # all y,a estimated
  }
  if(data$M_re_model == 3){ # ar1_a (devs by age, constant by year)
    for(i in 1:dim(tmp)[2]) tmp[,i] = i
  }
  if(data$M_re_model == 4){ # ar1_y (devs by year, constant by age)
    for(i in 1:dim(tmp)[1]) tmp[i,] = i
  }
  map$M_re <- factor(tmp)

  # M_repars: sigma_M, rho_M_a, rho_M_y
  if(data$M_re_model == 1) tmp <- rep(NA,3) # no RE pars to estimate
  if(data$M_re_model == 2) tmp <- c(1,NA,NA) # estimate sigma
  if(data$M_re_model == 3) tmp <- c(1,2,NA) # ar1_a: estimate sigma, rho_a
  if(data$M_re_model == 4) tmp <- c(1,NA,2) # ar1_y: estimate sigma, rho_y
  if(data$M_re_model == 5) tmp <- 1:3 # 2dar1: estimate all
  map$M_repars = factor(tmp)

  input$data = data
  input$par = par
  input$map = map
  return(input)

}