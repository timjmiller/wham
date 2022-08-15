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
	  data$waa_info = 1 # 1 = waa info present and used
	}
	else
	{
    if(is.null(waa_opts$waa)){
    	L = 100*(1-exp(-0.3*(1:data$n_ages - 0)))
	    W = rep(exp(-11)*L^3, each = data$n_years_model)
			dim_WAA = c(data$n_fleets + data$n_indices + 2, data$n_years_model, data$n_ages)
			WAA_pointers = c(1:data$n_fleets,data$n_fleets+data$n_indices + c(1,2,2))
			data$waa = array(2, dim = dim_WAA) # 2 kg for all ages, this will be replaced in WHAM
			data$waa_pointer_indices = rep(1,data$n_indices) + data$n_fleets
			data$waa_info = 0 # 0 = waa info not provided. use LW
		}
		else {
			data$waa = waa_opts$waa
			dim_waa = dim(data$waa)
			WAA_pointers = c(1:data$n_fleets,data$n_fleets+data$n_indices + c(1,2,2)) #Jan1 = SSB
			data$waa_pointer_indices = rep(1,data$n_indices) + data$n_fleets
			data$waa_info = 1 # 1 = waa info present and used
		}
	}

  data$waa_pointer_fleets = WAA_pointers[1:data$n_fleets]
  data$waa_pointer_totcatch = WAA_pointers[data$n_fleets + 1]
  data$waa_pointer_ssb = WAA_pointers[data$n_fleets + 2]
  data$waa_pointer_jan1 = WAA_pointers[data$n_fleets + 3]

  input$data = data
  return(input)
}