set_proj = function(input, proj.opts = NULL)
{
	data = input$data
	if(is.null(proj.opts))
	{
  	input$par$logR_proj <- matrix(0, data$n_stocks, 1) # will be set by prepare_projection if SCAA
  	input$map$logR_proj <- factor(rep(NA, data$n_stocks))
  
	  data$do_proj <- 0
	  data$n_years_proj <- 0
	  #data$n_years_proj_Ecov <- 0
  	data$avg_years_Ecov <- data$n_years_model - (5:1) # c++ indices start at 0
	  data$proj_F_opt <- 1
	  data$proj_Fcatch <- matrix(0)
	  data$proj_M_opt <- 1
	  data$proj_R_opt <- 1
		data$proj_L_opt <- 1
		data$proj_mu_opt <- 1
		data$mature_proj <- array(0,dim =c(1,1,1))
		data$waa_proj <- array(0,dim =c(1,1,1))
		if(!is.null(data$n_Ecov)){ #sometimes set_proj is called within other set_X functions before set_ecov has been called.
			data$proj_Ecov_opt <- rep(1, data$n_Ecov)
			data$Ecov_use_proj <- matrix(0, 1, data$n_Ecov)
		}
	  data$logR_mean <- 0 # only used for SCAA projections
	  data$logR_sd <- 0 # only used for SCAA projections
	  data$FXSPR_init = rep(0.5, data$n_years_model + data$n_years_proj)
	  data$FMSY_init = rep(0.5, data$n_years_model + data$n_years_proj)
	  data$F_proj_init = 0.1
		data$percentFMSY = 100
		data$percentFXSPR = 100
	}
	else {

		#do nothing?
		#prepare_projection requires a model object returned by fit_wham

	}
	input$data = data
	input$options$proj <- proj.opts
	return(input)
}