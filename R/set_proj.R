set_proj <- function(input, proj.opts = NULL)
{
	data <- input$data
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
		avg_years_ind <- tail(1:data$n_years_model,5) - 1
		# data$avg_years_ind <- avg_years_ind #default values to average FAA, M, movement, WAA, and maturity for projections (need to separate these)
		
		data$avg_years_ind_L <- matrix(0, data$n_years_model+1, data$n_regions) 
		data$avg_years_ind_L[1:(length(avg_years_ind) +1),] <- c(length(avg_years_ind),avg_years_ind) #default values to average L for projections
		
		tmat <- matrix(0, data$n_years_model+1, data$n_fleets) 
		tmat[1:(length(avg_years_ind) +1),] <- c(length(avg_years_ind),avg_years_ind) 
		data$avg_years_ind_sel <- data$avg_years_ind_waacatch <- tmat

		ind <- as.matrix(expand.grid(1:(length(avg_years_ind) +1), 1:data$n_stocks,1:data$n_regions))
		tarray <- array(0, dim = c(data$n_stocks,data$n_regions,data$n_years_model+1)) 
		tarray[ind[,c(2,3,1)]] <- c(length(avg_years_ind),avg_years_ind)
		data$avg_years_ind_M <- data$avg_years_ind_mat <- data$avg_years_ind_waassb <- data$avg_years_ind_move <- tarray
		
		ind <- as.matrix(expand.grid(1:(length(avg_years_ind) +1), 1:data$n_stocks,1:data$n_regions, 1:data$n_ages))
		tarray2 <- array(0, dim = c(data$n_stocks,data$n_regions,data$n_ages,data$n_years_model+1)) 
		tarray2[ind[,c(2,3,4,1)]] <- c(length(avg_years_ind),avg_years_ind)
		data$avg_years_ind_NAA <- tarray2
	  
	  data$proj_NAA_opt <- 1
	  
	  data$logR_mean <- 0 # only used for SCAA projections
	  data$logR_sd <- 0 # only used for SCAA projections
	  data$FXSPR_init <- c(data$FXSPR_init, rep(data$FXSPR_init[data$n_years_model], data$n_years_proj))
	  data$FMSY_init <- c(data$FMSY_init, rep(data$FMSY_init[data$n_years_model], data$n_years_proj))
	  data$F_proj_init <- 0.1
		data$percentFMSY <- 100
		data$percentFXSPR <- 100
	}
	else {

		#do nothing?
		#prepare_projection requires a model object returned by fit_wham

	}
	input$data <- data
	input$options$proj <- proj.opts
	return(input)
}