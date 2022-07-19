set_catch = function(input, catch_opts= NULL)
{
  data = input$data
  if(is.null(input$asap3)){
    asap3 = NULL
    # print("here")
    # print(is.null(catch_opts$n_fleets))
    if(is.null(catch_opts$n_fleets)) data$n_fleets = 1
    else data$n_fleets = catch_opts$n_fleets
  } else {
    asap3 = input$asap3
    data$n_fleets = asap3$n_fleets
  }
  # print(data)

	data$agg_catch = matrix(NA, data$n_years_model, data$n_fleets)
  data$catch_paa = array(NA, dim = c(data$n_fleets, data$n_years_model, data$n_ages))
  data$use_agg_catch = matrix(1, data$n_years_model, data$n_fleets)
  data$use_catch_paa = matrix(0, data$n_years_model, data$n_fleets)
  data$catch_pal = array(NA, dim = c(data$n_fleets, data$n_years_model, data$n_lengths))
  data$use_catch_pal = matrix(0, data$n_years_model, data$n_fleets)  
  data$catch_caal = array(NA, dim = c(data$n_fleets, data$n_years_model, data$n_lengths, data$n_ages))
  data$use_catch_caal = array(0, dim = c(data$n_years_model, data$n_fleets))
  
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
	    if(asap3$use_catch_acomp[i] != 1){
        data$use_catch_paa[,i] = 0
      } else { # use catch paa in at least some years - not necessarily all, have to go through year by year
        for(y in 1:data$n_years_model){
          if(is.na(sum(data$catch_paa[i,y,] > 1e-15))){ # handle negative or NA paa
            data$use_catch_paa[y,i] = 0
          } else {
            if(asap3$catch_Neff[y,i] < 1e-15 | sum(data$catch_paa[i,y,] > 1e-15)<2) data$use_catch_paa[y,i] = 0
          }
        } 
      }
	  }
	  data$catch_Neff = asap3$catch_Neff	  
    data$selblock_pointer_fleets = cbind(sapply(asap3$sel_block_assign, function(x) return(x)))
  }
  else
  {
    data$n_fleets = 1
    if(is.null(catch_opts$agg_catch)) data$agg_catch[] = 1
    else data$agg_catch = catch_opts$agg_catch
  	
  	if(is.null(catch_opts$catch_paa)) data$catch_paa[] = 1/data$n_ages
    else data$catch_paa[] = catch_opts$catch_paa

  	if(is.null(catch_opts$catch_pal)) data$catch_pal[] = 1/data$n_lengths
    else data$catch_pal[] = catch_opts$catch_pal
	
    if(is.null(catch_opts$catch_caal)) data$catch_caal[] = 1/data$n_ages
    else data$catch_caal[] = catch_opts$catch_caal

	if(is.null(catch_opts$catch_cv)) data$agg_catch_sigma = matrix(sqrt(log(0.1^2 + 1)), data$n_years_model, data$n_fleets)
    else data$agg_catch_sigma = matrix(sqrt((log(catch_opts$catch_cv^2 + 1))), data$n_years_model, data$n_fleets)
	  
    if(is.null(catch_opts$catch_Neff)) data$catch_Neff = matrix(0, data$n_years_model, data$n_fleets)	  
    else data$catch_Neff = catch_opts$catch_Neff

    if(is.null(catch_opts$catch_NeffL)) data$catch_NeffL = matrix(0, data$n_years_model, data$n_fleets)	  
    else data$catch_NeffL = catch_opts$catch_NeffL

    if(is.null(catch_opts$catch_caal_Neff)) data$catch_caal_Neff = array(0, dim = c(data$n_years_model, data$n_fleets, data$n_lengths))
    else data$catch_caal_Neff = catch_opts$catch_caal_Neff

    if(!is.null(catch_opts$use_catch_paa)) data$use_catch_paa[] = catch_opts$use_catch_paa
    if(!is.null(catch_opts$use_catch_pal)) data$use_catch_pal[] = catch_opts$use_catch_pal
    if(!is.null(catch_opts$use_catch_caal)) data$use_catch_caal[] = catch_opts$use_catch_caal

    for(i in 1:data$n_fleets) for(y in 1:data$n_years_model){ 
      if(data$catch_Neff[y,i] < 1e-15 | sum(data$catch_paa[i,y,] > 1e-15)<2 | any(is.na(data$catch_paa[i,y,]))) data$use_catch_paa[y,i] = 0
    }
    for(i in 1:data$n_fleets) for(y in 1:data$n_years_model){ 
      if(data$catch_NeffL[y,i] < 1e-15 | sum(data$catch_pal[i,y,] > 1e-15)<2 | any(is.na(data$catch_pal[i,y,]))) data$use_catch_pal[y,i] = 0
    }
    if(is.null(catch_opts$selblock_pointer_fleets)) data$selblock_pointer_fleets = matrix(rep(1:data$n_fleets, each = data$n_years_model), data$n_years_model, data$n_fleets)
    else data$selblock_pointer_fleets = catch_opts$selblock_pointer_fleets
  }

  data$catch_paa[is.na(data$catch_paa)] = 0
  data$catch_pal[is.na(data$catch_pal)] = 0

  input$par$log_catch_sig_scale = rep(0, data$n_fleets)
  input$map$log_catch_sig_scale = factor(rep(NA, data$n_fleets))

  input$data = data
  return(input)
}