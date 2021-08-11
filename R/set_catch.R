set_catch = function(input, catch_opts= NULL)
{
  data = input$data
  if(is.null(input$asap3)){
    asap3 = NULL
    print("here")
    print(is.null(catch_opts$n_fleets))
    if(is.null(catch_opts$n_fleets)) data$n_fleets = 1
    else data$n_fleets = catch_opts$n_fleets
  } else {
    asap3 = input$asap3
    data$n_fleets = asap3$n_fleets
  }
  print(data)

	data$agg_catch = matrix(NA, data$n_years_model, data$n_fleets)
  data$catch_paa = array(NA, dim = c(data$n_fleets, data$n_years_model, data$n_ages))
  data$use_agg_catch = matrix(1, data$n_years_model, data$n_fleets)
  data$use_catch_paa = matrix(1, data$n_years_model, data$n_fleets)
	if(!is.null(asap3))
	{
    data$n_fleets = asap3$n_fleets
    asap3$use_catch_acomp <- rep(1,asap3$n_fleets) #default is to use age comp for catch
  	for(i in 1:data$n_fleets) data$agg_catch[,i] = asap3$CAA_mats[[i]][,data$n_ages + 1]
  	data$agg_catch_sigma = asap3$catch_cv
    data$agg_catch_sigma[which(data$agg_catch_sigma < 1e-15)] = 100
    data$agg_catch_sigma = sqrt(log(data$agg_catch_sigma^2 + 1))
	  for(i in 1:data$n_fleets)
	  {
	    temp = asap3$CAA_mats[[i]][,1:data$n_ages]
	    temp[which(is.na(temp))] = 0
	    temp[which(temp<0)] = 0
	    data$catch_paa[i,,] = temp/apply(temp,1,sum)
	    for(y in 1:data$n_years_model) if(asap3$CAA_mats[[i]][y,data$n_ages+1] < 1e-15) data$use_agg_catch[y,i] = 0
	    if(asap3$use_catch_acomp[i] != 1) data$use_catch_paa[,i] = 0
	    else for(y in 1:data$n_years_model) if(asap3$catch_Neff[y,i] < 1e-15 | sum(data$catch_paa[i,y,] > 1e-15)<2) data$use_catch_paa[y,i] = 0
	  }
	  data$catch_Neff = asap3$catch_Neff	  
  }
  else
  {
    data$n_fleets = 1
    if(is.null(catch_opts$agg_catch)) data$agg_catch[] = 100
    else data$agg_catch = catch_opts$agg_catch
  	
  	if(is.null(catch_opts$catch_paa)) data$catch_paa[] = 1/data$n_ages
    else data$catch_paa[] = catch_opts$catch_paa

		if(is.null(catch_opts$agg_catch_sigma)) data$agg_catch_sigma = matrix(sqrt(log(0.1^2 + 1)), data$n_years_model, data$n_fleets)
    else data$agg_catch_sigma = catch_opts$agg_catch_sigma
	  
    if(is.null(catch_opts$catch_Neff)) data$catch_Neff = matrix(200, data$n_years_model, data$n_fleets)	  
    else data$catch_Neff = catch_opts$catch_Neff

    for(i in 1:data$n_fleets) for(y in 1:data$n_years_model){ 
      if(data$catch_Neff[y,i] < 1e-15 | sum(data$catch_paa[i,y,] > 1e-15)<2 | any(is.na(data$catch_paa[i,y,]))) data$use_catch_paa[y,i] = 0
    }
  }

  data$catch_paa[is.na(data$catch_paa)] = 0

  data$catch_aref = matrix(NA, data$n_years_model, data$n_fleets)
  for(i in 1:data$n_fleets) data$catch_aref[,i] = get_aref_fn(data$catch_paa[i,,])

  input$par$log_catch_sig_scale = rep(0, data$n_fleets)
  input$map$log_catch_sig_scale = factor(rep(NA, data$n_fleets))

  input$data = data
  return(input)
}