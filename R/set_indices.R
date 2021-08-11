set_indices = function(input, index_opts=NULL)
{
	data = input$data
	if(is.null(input$asap3)) {
		asap3 = NULL
    if(is.null(index_opts$n_indices)) data$n_indices = 1
    else data$n_indices = index_opts$n_indices
	}
  else {
  	asap3 = input$asap3
 	  which_indices <- which(asap3$use_index ==1)
	  asap3$n_indices = length(which_indices)
	  asap3$survey_index_units <- asap3$index_units[which_indices]
	  asap3$survey_acomp_units <- asap3$index_acomp_units[which_indices]
	  asap3$survey_WAA_pointers <- asap3$index_WAA_pointers[which_indices]
	  asap3$survey_month <- matrix(asap3$index_month[which_indices], asap3$n_years, asap3$n_indices, byrow = TRUE)
	  asap3$use_survey_acomp <- asap3$use_index_acomp[which_indices]
	  asap3$index_sel_option <- asap3$index_sel_option[which_indices]
	  asap3$index_sel_ini = asap3$index_sel_ini[which_indices]
	  asap3$index_WAA_pointers = asap3$index_WAA_pointers[which_indices]
	  asap3$IAA_mats <- asap3$IAA_mats[which_indices]
	  asap3$use_survey <- asap3$use_index[which_indices]
	  data$n_indices <- asap3$n_indices
	} 
  data$agg_indices = matrix(NA, data$n_years_model, data$n_indices)
  data$use_indices = matrix(1, data$n_years_model, data$n_indices)
  data$agg_index_sigma = matrix(NA, data$n_years_model, data$n_indices)
  data$index_paa = array(NA, dim = c(data$n_indices, data$n_years_model, data$n_ages))
  data$use_index_paa = matrix(1, data$n_years_model, data$n_indices)
  data$index_Neff = matrix(NA, data$n_years_model, data$n_indices)
  data$q_lower <- rep(0,data$n_indices)
  data$q_upper <- rep(1000,data$n_indices)
	if(!is.null(asap3))
	{
	  data$units_indices <- asap3$survey_index_units
	  data$fracyr_indices = (asap3$survey_month-1)/12 #make sure that this is right
  	for(i in 1:data$n_indices) data$agg_indices[,i] = asap3$IAA_mats[[i]][,2]
	  for(i in 1:data$n_indices)
	  {
	    for(y in 1:data$n_years_model) if(asap3$IAA_mats[[i]][y,2] < 1e-15) data$use_indices[y,i] = 0
	  }
	  for(i in 1:data$n_indices) data$agg_index_sigma[,i] = asap3$IAA_mats[[i]][,3]
	  for(i in 1:data$n_indices)
	  {
	    temp = asap3$IAA_mats[[i]][,3 + 1:data$n_ages]
	    temp[which(is.na(temp))] = 0
	    temp[which(temp<0)] = 0
	    data$index_paa[i,,] = temp/apply(temp,1,sum)
	    #data$index_paa[i,,][sign(asap3$IAA_mats[[i]][,3 + 1:data$n_ages]) == -1] = 0
	  }
	  for(i in 1:data$n_indices)
	  {
	    if(asap3$use_survey_acomp[i] != 1) data$use_index_paa[,i] = 0
	    else for(y in 1:data$n_years_model) if(asap3$IAA_mats[[i]][y,4 + data$n_ages] < 1e-15 | sum(data$index_paa[i,y,] > 1e-15)<2) data$use_index_paa[y,i] = 0
	  }
	  data$units_index_paa <- asap3$survey_acomp_units
	  for(i in 1:data$n_indices) data$index_Neff[,i] = asap3$IAA_mats[[i]][,4 + data$n_ages]

  	input$par$logit_q = gen.logit(asap3$q_ini, data$q_lower, data$q_upper) # use q_ini values from asap3 file
	}
	else
	{
		if(is.null(index_opts$units_indices)) data$units_indices = rep(1,data$n_indices) #biomass
		else data$units_indices = index_opts$units_indices
		
		if(is.null(index_opts$fracyr_indices)) data$fracyr_indices =matrix(0.5, data$n_years_model, data$n_indices)
		else data$fracyr_indices = index_opts$fracyr_indices
		
		if(is.null(index_opts$agg_indices)) data$agg_indices[] = 10
		else 

		if(is.null(index_opts$agg_index_sigma)) data$agg_index_sigma[] = sqrt(log(0.3^2 + 1))
		else 

		if(is.null(index_opts$index_paa)) data$index_paa[] = 1/data$n_ages
		else 
		if(is.null(index_opts$units_index_paa)) data$units_index_paa = rep(2,data$n_indices) #numbers
		
		if(is.null(index_opts$index_Neff)) data$index_Neff[] = 100

		#These are specified at the top, so only change if needed
		if(!is.null(index_opts$q_lower)) data$q_lower = index_opts$q_lower
		if(!is.null(index_opts$q_upper)) data$q_upper = index_opts$q_upper
		
		if(is.null(index_opts$q)) input$par$logit_q[] = gen.logit(0.3, data$q_lower, data$q_upper)
		else input$par$logit_q[] = gen.logit(index_opts$q, data$q_lower, data$q_upper)
	}

  data$agg_index_sigma[which(data$agg_index_sigma < 1e-15)] = 100
  data$agg_index_sigma = sqrt(log(data$agg_index_sigma^2 + 1))
  data$index_paa[is.na(data$index_paa)] = 0
  
  data$index_aref = matrix(NA, data$n_years_model, data$n_indices)
  for(i in 1:data$n_indices) data$index_aref[,i] = get_aref_fn(data$index_paa[i,,])  

  input$par$log_index_sig_scale = rep(0, data$n_indices)
  input$map$log_index_sig_scale = factor(rep(NA, data$n_indices))

	
  input$data = data
  return(input)
}