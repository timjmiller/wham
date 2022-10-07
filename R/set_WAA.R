set_WAA = function(input, waa_opts = NULL)
{
	# Weight-at-age
  asap3 = if(is.null(input$asap3)) asap3 = NULL
  else asap3 = input$asap3
	data = input$data
	if(!is.null(asap3))
	{
	  i <- c(seq(1,(asap3$n_fleets+1)*2-1,2),(asap3$n_fleets+1)*2 + 1:2)
	  WAA_pointers <- asap3$WAA_pointers[i] #wham has no discard data, so remove those WAA matrices
	  data$waa = array(NA, dim = c(length(asap3$WAA_mats), data$n_years_model, data$n_ages))
	  for(i in 1:length(asap3$WAA_mats)) data$waa[i,,] = asap3$WAA_mats[[i]]
	  data$waa_pointer_indices = asap3$index_WAA_pointers
	  data$waa_type = 1 # 1 = waa info present and used
	  data$waa_cv = array(0, dim = dim(data$waa))
	  data$use_catch_waa = matrix(0, nrow = data$n_years_model, ncol = data$n_fleets)
	  data$use_index_waa = matrix(0, nrow = data$n_years_model, ncol = data$n_indices)
	} else {
    	if(is.null(waa_opts$waa)){
	    	L = 100*(1-exp(-0.3*(1:data$n_ages - 0)))
		    W = rep(exp(-11)*L^3, each = data$n_years_model)
			dim_WAA = c(data$n_fleets + data$n_indices + 2, data$n_years_model, data$n_ages)
			data$waa = array(2, dim = dim_WAA) # 2 kg for all ages, this will be replaced in WHAM	
			data$waa_type = 2 # 2 = waa info not provided. use LW
			data$waa_cv = array(0, dim = dim(data$waa))		
	  		data$use_catch_waa = matrix(0, nrow = data$n_years_model, ncol = data$n_fleets)
	  		data$use_index_waa = matrix(0, nrow = data$n_years_model, ncol = data$n_indices)
		} else {
			data$waa = waa_opts$waa
			data$waa_type = 1
			data$waa_cv = array(0, dim = dim(data$waa))
			data$use_catch_waa = matrix(0, nrow = data$n_years_model, ncol = data$n_fleets)
			data$use_index_waa = matrix(0, nrow = data$n_years_model, ncol = data$n_indices)
		}
		
		if(is.null(waa_opts$waa_pointer_fleets)) data$waa_pointer_fleets = 1:data$n_fleets
		else data$waa_pointer_fleets = waa_opts$waa_pointer_fleets

		if(is.null(waa_opts$waa_pointer_indices)) data$waa_pointer_indices = 1:data$n_indices + data$n_fleets
		else data$waa_pointer_indices = waa_opts$waa_pointer_indices		

		if(is.null(waa_opts$waa_pointer_totcatch)) data$waa_pointer_totcatch = data$n_fleets+data$n_indices + 1
		else data$waa_pointer_totcatch = waa_opts$waa_pointer_totcatch

		if(is.null(waa_opts$waa_pointer_ssb)) data$waa_pointer_ssb = data$n_fleets+data$n_indices + 2
		else data$waa_pointer_ssb = waa_opts$waa_pointer_ssb

		if(is.null(waa_opts$waa_pointer_jan1)) data$waa_pointer_jan1 = data$n_fleets+data$n_indices + 2
		else data$waa_pointer_jan1 = waa_opts$waa_pointer_jan1

		if(!is.null(waa_opts$waa_type))  data$waa_type = waa_opts$waa_type

		if(!is.null(waa_opts$waa_cv)) data$waa_cv = waa_opts$waa_cv

		if(!is.null(waa_opts$use_catch_waa)) data$use_catch_waa = waa_opts$use_catch_waa

		if(!is.null(waa_opts$use_index_waa)) data$use_index_waa = waa_opts$use_index_waa
	}

  input$data = data
  return(input)
}