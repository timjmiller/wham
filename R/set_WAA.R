set_WAA = function(input)
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
	}
	else
	{
    L = 100*(1-exp(-0.3*(1:data$n_ages - 0)))
    W = rep(exp(-11)*L^3, each = data$n_years_model)
		dim_WAA = c(1, data$n_years_model, data$n_ages)
		WAA_pointers = rep(1,data$n_fleets + 3)
		data$waa = array(W, dim = dim_WAA)
		data$waa_pointer_indices = rep(1,data$n_indices)
	}

	data$waa_pointer_fleets = WAA_pointers[1:data$n_fleets]
  data$waa_pointer_totcatch = WAA_pointers[data$n_fleets + 1]
  data$waa_pointer_ssb = WAA_pointers[data$n_fleets + 2]
  data$waa_pointer_jan1 = WAA_pointers[data$n_fleets + 3]

  input$data = data
  return(input)
}