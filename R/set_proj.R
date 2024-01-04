set_proj = function(input, proj.opts = NULL)
{
	data = input$data
	if(is.null(proj.opts))
	{
	  data$do_proj <- 0
	  data$n_years_proj <- 0
	  data$n_years_proj_Ecov <- 0
  	#data$avg_years_ind <- data$n_years_model - (5:1) # set in prepare_wham_input, needed even without projections
	  #data$avg_years_ind <- 0
	  data$proj_F_opt <- 0
	  data$proj_Fcatch <- 0
	  data$proj_M_opt <- 0
	  data$proj_R_opt <- 0
	  data$logR_mean <- 0 # only used for SCAA projections
	  data$logR_sd <- 0 # only used for SCAA projections
	  data$FXSPR_init = rep(0.1, data$n_years_model + data$n_years_proj)
	  data$FMSY_init = rep(0.1, data$n_years_model + data$n_years_proj)
	  data$F_proj_init = 0
	}
	else {

		#do nothing?
		#prepare_projection requires a model object returned by fit_wham

	}
	input$data = data
	return(input)
}