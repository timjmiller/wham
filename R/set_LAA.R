set_LAA = function(input, LAA)
{
  data = input$data
  par = input$par
  map = input$map
    
  n_re_par = 3 # number of parameters for RE

  # if(is.null(input$asap3)) asap3 = NULL
  # else {
  #   stop('LAA does not work with ASAP3 input')
  # }

  # LAA default options:
  LAA_re_ini = matrix(0, ncol = data$n_ages, nrow = data$n_years_model)
  data$LAA_re_model = 1
  data$LAA_est = rep(0, data$n_ages)
  LAA_ini = log( 100 + (3 - 100)*exp(-0.2*(1:data$n_ages - 1)) )
  if(!is.null(LAA)) {
    LAA_ini = LAA$LAA_vals
    if(!is.null(LAA$LAA_est)) data$LAA_est[LAA$LAA_est] = 1
    if(!is.null(LAA$re))  {
      if(!(LAA$re %in% c("none","iid","iid_a","ar1_a","2dar1"))) stop("LAA$re must be one of the following: 'none','iid','iid_a','ar1_a','2dar1'")
      data$LAA_re_model <- match(LAA$re, c("none","iid","iid_a","ar1_a","2dar1")) # Respect this order to create array later
    }
  }

  # growth pars --------------------------
  
  par$LAA_a = log(LAA_ini)
  par$LAA_re = LAA_re_ini
  par$LAA_repars = rep(0, times = n_re_par)
  par$LAA_repars[1] = log(0.1) # start sigma at 0.1, rho at 0


  # --------------------------------
  # Prepare data for growth:
  tmp1 <- par$LAA_a
  tmp1[data$LAA_est==0] = NA
  ind.notNA <- which(!is.na(tmp1))
  tmp1[ind.notNA] <- 1:length(ind.notNA)
  map$LAA_a <- factor(tmp1)

  # RE info:
  # LAA_re: "none","iid","iid_a","ar1_a","2dar1"
  tmp <- par$LAA_re
  if(data$LAA_re_model == 1) tmp[] = NA # no RE (either estimate RE for all ages or none at all)
  if(data$LAA_re_model %in% c(2,5)){ # 2d ar1
    tmp[] = 1:(dim(tmp)[1]*dim(tmp)[2]) # all y,a estimated
  }
  if(data$LAA_re_model %in% c(3,4)){ # ar1_a (devs by age, constant by year)
    for(i in 1:dim(tmp)[2]) tmp[,i] = i
  }
  map$LAA_re <- factor(tmp)

  # LAA_repars: 
  if(data$LAA_re_model == 1) tmp <- rep(NA,3) # no RE pars to estimate
  if(data$LAA_re_model == 2) tmp <- c(1,NA,NA) # estimate sigma
  if(data$LAA_re_model == 3) tmp <- c(1,NA,NA) # estimate sigma over ages
  if(data$LAA_re_model == 4) tmp <- c(1,2,NA) # ar1_a: estimate sigma, rho_a
  if(data$LAA_re_model == 5) tmp <- 1:3 # 2dar1: estimate all
  map$LAA_repars = factor(tmp)

  # End section

  input$data = data
  input$par = par
  input$map = map
  return(input)

}