#' Reduce the years of the model
#'
#' Internal function called by \code{\link{fit_peel}} for \emph{i} in 1--\code{n.peels}.
#' Creates the input for the model peeling off \emph{i} years of data (calls \code{\link{fit_tmb}}).
#' Reduces the year dimension of all data, parameters, maps, so that one can project correctly from the peeled model.
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{prepare_wham_input}}). NOT from prepare_projection or project_wham
#' @param years_peeled which of input$years to peel from the model input
#' @param retro (T/F) whether this is for a retro peel (Default = TRUE)
#'
reduce_input <- function(input,years_peeled, retro = TRUE){
	
	if(input$data$n_years_proj>0) stop("cannot run wham:::reduce_input on an input for a projection")
	new <- input
	data <- new$data
	par <- new$par
	map <- new$map
	ind <- which(!(input$years %in% years_peeled))
	years <- input$years[ind]
	n_years <- length(years)
	new$years <- new$years_full <- years


	data$n_years_model <- n_years
	data$years_use <- data$years_use[ind]
	data$XSPR_R_avg_yrs <- data$XSPR_R_avg_yrs[which(data$XSPR_R_avg_yrs %in% (ind-1))]
	if(length(data$XSPR_R_avg_yrs)<1) {
		warning("Years specified to use for average recruitment in SPR-based reference points are not included in years of reduced input so all years will be used.")
		data$XSPR_R_avg_yrs <- 1:n_years - 1
	}
	data$avg_years_ind <- data$n_years_model - input$data$n_years_model + data$avg_years_ind#shift which years to average back given the number of years to peel.
	if(all(data$avg_years_ind<0)) data$avg_years_ind <- 1:data$n_years_model - 1 #shouldn't really ever happen
	data$avg_years_ind <- data$avg_years_ind[which(data$avg_years_ind>=0)] #reduce if number of years used is more than the number available in the peel.

	data$which_F_age <- data$which_F_age[ind]
	data$FXSPR_init <- data$FXSPR_init[ind]
	data$FMSY_init <- data$FMSY_init[ind]

	data$fracyr_SSB <- data$fracyr_SSB[ind,,drop=F]
	data$mature <- data$mature[,ind,,drop=F]
	data$catch_Neff <- data$catch_Neff[ind,,drop=F]
	data$agg_catch <- data$agg_catch[ind,,drop=F]
	data$agg_catch_sigma <- data$agg_catch_sigma[ind,,drop=F]
	data$catch_paa <- data$catch_paa[,ind,,drop=F]
	data$use_agg_catch <- data$use_agg_catch[ind,,drop=F]
	data$use_catch_paa <- data$use_catch_paa[ind,,drop=F]
	data$selblock_pointer_fleets <- data$selblock_pointer_fleets[ind,,drop=F]
	data$index_Neff <- data$index_Neff[ind,,drop=F]
	data$agg_index_sigma <- data$agg_index_sigma[ind,,drop=F]
	data$agg_indices <- data$agg_indices[ind,,drop=F]
	data$index_paa <- data$index_paa[,ind,,drop=F]
	data$use_indices <- data$use_indices[ind,,drop=F]
	data$use_index_paa <- data$use_index_paa[ind,,drop=F]
	data$selblock_pointer_indices <- data$selblock_pointer_indices[ind,,drop=F]
	data$fracyr_indices <- data$fracyr_indices[ind,,drop=F]
	data$waa <- data$waa[,ind,,drop=F]
	data$selblock_years <- data$selblock_years[ind,,drop=F]
	data$n_years_selblocks <- apply(data$selblock_years,2,sum)

	data$keep_C <- data$keep_C[ind,,drop=F]
	data$keep_I <- data$keep_I[ind,,drop=F]
	data$keep_Cpaa <- data$keep_Cpaa[,ind,,drop=F]
	data$keep_Ipaa <- data$keep_Ipaa[,ind,,drop=F]

	if(!is.null(map$q_re)) {
		tmp <- par$q_re
		tmp[] <- as.integer(map$q_re)
		map$q_re <- factor(tmp[ind,])	
	}
	par$q_re <- par$q_re[ind,,drop=F]

	if(!is.null(map$mu_re)) {
		tmp <- par$mu_re
		tmp[] <- as.integer(map$mu_re)
		map$mu_re <- factor(tmp[,,,ind,,])	
	}
	par$mu_re <- par$mu_re[,,,ind,,,drop=F]

	if(!is.null(map$M_re)) {
		tmp <- par$M_re
		tmp[] <- as.integer(map$M_re)
		map$M_re <- factor(tmp[,,ind,])	
	}
	par$M_re <- par$M_re[,,ind,,drop=F]
	
	if(!is.null(map$log_NAA)) {
		tmp <- par$log_NAA
		tmp[] <- as.integer(map$log_NAA)
		map$log_NAA <- factor(tmp[,,head(ind,n_years-1),])	
	}
	par$log_NAA<- par$log_NAA[,,head(ind,n_years-1),,drop = F] 
	
	if(!is.null(map$L_re)) {
		tmp <- par$L_re
		tmp[] <- as.integer(map$L_re)
		map$L_re <- factor(tmp[ind,])	
	}
	par$L_re <- par$L_re[ind,,drop=F]
	
	if(!is.null(map$selpars_re)){
		tmp <- par$selpars_re
		tmp[] <- as.integer(map$selpars_re)
		map$selpars_re <- factor(tmp[,ind,]) #n_selblocks x n_years x n_ages
  }
  par$selpars_re <- par$selpars_re[,ind,,drop=F]

	map$F_pars <- par$F_pars <- par$F_pars[ind,,drop=F]
	map$F_pars[] <- 1:(length(ind)*data$n_fleets)
	map$F_pars <- as.factor(map$F_pars)

  dims <- dim(par$q_re)
  dims[2] <- data$n_years_model + data$n_years_proj
  tmp <- matrix(0, nrow = data$n_years_model + data$n_years_proj, ncol = NCOL(par$q_re))
  tmp[1:data$n_years_model,] <- par$q_re[1:data$n_years_model,]
  par$q_re <- tmp

	new$par <- par
	new$map <- map
	new$data <- data

	#adjust Ecov dimensions
	if(!is.null(new$options$ecov)){
		ecov <- new$options$ecov
		if(retro) ind_E <- which(ecov$year < min(years_peeled))
		else ind_E <- which(!(ecov$year %in% years_peeled))
	  ecov$year <- ecov$year[ind_E]
		ecov$mean <- ecov$mean[ind_E,,drop=F]
		ecov$use_obs <- ecov$use_obs[ind_E,,drop=F]
		if(is.list(ecov$logsigma)) if(is.matrix(ecov$logsigma[[1]])){
			ecov$logsigma[[1]] <- ecov$logsigma[[1]][ind_E,,drop = FALSE]
		}
		if(is.matrix(ecov$logsigma)) ecov$logsigma <- ecov$logsigma[ind_E,,drop = FALSE]
		new <- set_ecov(new,ecov)
		#else set_ecov should define the pars correctly
	}

	#reset the obs_vec
	new <- set_osa_obs(new)
	return(new)

}