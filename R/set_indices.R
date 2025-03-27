#' Specify index selectivity blocks and aggregate and age composition observations for indices
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{prepare_wham_input}})
#' @param index_info (optional) list specifying various aspects about catch by indices (see details)
#' 
#' \code{index_info} specifies observations, and various configuration options for index-specific catch observations and will overwrite attributes specified in the ASAP data file.
#' If \code{NULL}, all settings from the ASAP data file or basic_info are used.
#' \code{index_info} is a list with any of the following entries:
#'   \describe{
#'     \item{$n_indices}{number of indices}
#'     \item{$index_regions}{vector (n_indices) of regions where each fleet operates.}
#'     \item{$index_seasons}{vector (n_indices) of 0/1 values flagging which seasons each index occurs.}
#'     \item{$index_names}{character vector (n_indices) of names for indices Used for naming results in plots and tables.}
#'     \item{$agg_indices}{matrix (n_years_model x n_indices) of annual aggregate index catches.}
#'     \item{$agg_index_cv}{matrix (n_years_model x n_indices) of CVs for annual aggregate index catches.}
#'     \item{$fracyr_indices}{matrix (n_years_model x n_indices) of fractions of year at which index occurs within the season (difference between time of survey and time at start of season).}
#'     \item{$use_indices}{matrix (n_years_model x n_indices) of 0/1 values flagging whether to use aggregate observations.}
#'     \item{$units_indices}{matrix (n_years_model x n_indices) of 1/2 values flagging whether aggregate observations are biomass (1) or numbers (2).}
#'     \item{$index_paa}{array (n_indices x n_years_model x n_ages) of annual catch proportions at age by index.}
#'     \item{$use_index_paa}{matrix (n_years_model x n_indices) of 0/1 values flagging whether to use proportions at age observations.}
#'     \item{$units_index_paa}{matrix (n_years_model x n_indices) of 1/2 values flagging whether composition observations are biomass (1) or numbers (2).}
#'     \item{$index_Neff}{matrix (n_years_model x n_indices) of effective sample sizes for proportions at age observations.}
#'     \item{$waa_pointer_indices}{vector (n_indices) of itegers indicated waa to use for each index.}
#'     \item{$selblock_pointer_indices}{matrix (n_years_model x n_indices) of itegers indicated selblocks to use.}
#'     \item{$initial_index_sd_scale}{vector (n_indices) of scalar multipliers of annual log-observation standard deviation. Default = 1.}
#'     \item{$map_index_sd_scale}{integer vector (n_indices) specifying which sd scalar parameters to fix. Use \code{NA} to fix a parameter and integers to estimate. 
#'         Use the same integer for multiple indices to estimate a shared scalar parameter.}
#'   }
#'
#' @export
set_indices <- function(input, index_info=NULL) {
  if(is.null(input$use_asap3)){
  	input$use_asap3 <- TRUE
  	if(is.null(input$asap3)) input$use_asap3 <- FALSE
  }
	data <- input$data
	if(is.null(data)) data <- list()
  input$log$indices <- list()
	if(input$use_asap3) {
  	asap3 <- input$asap3
		for(i in 1:length(asap3)){
			which_indices <- which(asap3[[i]]$use_index ==1)
			asap3[[i]]$n_indices <- length(which_indices)
		}
    n_indices_per_region <- sapply(asap3, function(x) x$n_indices)
    data$n_indices <- sum(n_indices_per_region)
	} else{ #no asap and no number of indices provided
		if(is.null(data$n_indices)) data$n_indices <- 1
	}
	if(is.null(index_info)) index_info <- list()
	if(!is.null(index_info$n_indices)) data$n_indices <- index_info$n_indices
  # index_data_names <- c("n_indices", "agg_indices","index_paa","use_indices","use_index_paa","units_indices", "units_index_paa","index_Neff", "agg_index_sigma",
  # 	"index_regions", "index_seasons","fracyr_indices", "selblock_pointer_indices", "waa_pointer_indices")
	# if(is.null(asap3)) {
	#   # data$n_indices <- 1
	# } #else {
	if(input$use_asap3){
    for(i in 1:length(asap3)) {
			which_indices <- which(asap3[[i]]$use_index ==1)
			asap3[[i]]$n_indices <- length(which_indices)
			asap3[[i]]$survey_index_units <- asap3[[i]]$index_units[which_indices]
			asap3[[i]]$survey_acomp_units <- asap3[[i]]$index_acomp_units[which_indices]
			asap3[[i]]$survey_WAA_pointers <- asap3[[i]]$index_WAA_pointers[which_indices]
			asap3[[i]]$survey_month <- asap3[[i]]$index_month[which_indices]
			#asap3[[i]]$survey_month <- matrix(asap3[[i]]$index_month[which_indices], asap3[[i]]$n_years, asap3[[i]]$n_indices, byrow = TRUE)
			asap3[[i]]$use_survey_acomp <- asap3[[i]]$use_index_acomp[which_indices]
			asap3[[i]]$index_WAA_pointers = asap3[[i]]$index_WAA_pointers[which_indices]
			asap3[[i]]$IAA_mats <- asap3[[i]]$IAA_mats[which_indices]
			asap3[[i]]$use_survey <- asap3[[i]]$use_index[which_indices]
			if(!is.null(asap3[[i]]$index.names)) input$index_names <- c(input$index_names, asap3[[i]]$index.names[which_indices])
    }
	  input$asap3 <- asap3
	  data$index_regions <- 0
	  data$index_seasons <- 0
    data$fracyr_indices <- matrix(0, data$n_years_model, data$n_indices)
    data$agg_indices <- matrix(0, data$n_years_model, data$n_indices)
    data$agg_index_sigma <- matrix(0, data$n_years_model, data$n_indices)
    data$index_paa <- array(0, c(data$n_indices,data$n_years_model,data$n_ages))
    data$units_indices <- 0
    data$units_index_paa <- 0
    data$use_indices <- matrix(1, data$n_years_model, data$n_indices)
    data$use_index_paa <- matrix(1, data$n_years_model, data$n_indices)
    data$index_Neff <- matrix(1, data$n_years_model, data$n_indices)
    data$index_Neff <- matrix(1, data$n_years_model, data$n_indices)
    data$selblock_pointer_indices <- matrix(1, data$n_years_model, data$n_indices)
    k <- 1
    for(i in 1:length(asap3)) {
      for(j in 1:asap3[[i]]$n_indices) {
        data$index_regions[k] <- i #each asap file is a separate region
		  	data$units_indices[k] <- asap3[[i]]$survey_index_units[j]
		  	tmp <- (asap3[[i]]$survey_month[j]-1)/12 #make sure that this is right
				int_starts <- cumsum(c(0,data$fracyr_seasons))
		  	ind <- max(which(int_starts <= tmp))
		  	data$index_seasons[k] <- ind
		  	data$fracyr_indices[,k] <- tmp - int_starts[ind]
				data$agg_indices[,k] <- asap3[[i]]$IAA_mats[[j]][,2]
	    	for(y in 1:data$n_years_model) if(asap3[[i]]$IAA_mats[[j]][y,2] < 1e-15) data$use_indices[y,k] <- 0
	  		data$agg_index_sigma[,k] <- asap3[[i]]$IAA_mats[[j]][,3]
		    temp <- asap3[[i]]$IAA_mats[[j]][,3 + 1:data$n_ages]
	    	temp[which(is.na(temp))] <- 0
	    	temp[which(temp<0)] <- 0
	    	data$index_paa[k,,] <- temp/apply(temp,1,sum) #all 0s will make NaN
				if(asap3[[i]]$use_survey_acomp[j] != 1){
					data$use_index_paa[,k] <- 0
				} else {
					for(y in 1:data$n_years_model) {
						flag <- asap3[[i]]$IAA_mats[[j]][y,4 + data$n_ages] < 1e-15 | sum(data$index_paa[k,y,] > 1e-15) < 2 | any(is.na(data$index_paa[k,y,]))
						if(flag) data$use_index_paa[y,k] <- 0
					}
				}
				data$units_index_paa[k] <- asap3[[i]]$survey_acomp_units[j]
				data$index_Neff[,k] <- asap3[[i]]$IAA_mats[[j]][,4 + data$n_ages]
        data$selblock_pointer_indices[,k] <- max(data$selblock_pointer_fleets) + k #set_catch already called
        k <- k + 1
			}
		}
		data$waa_pointer_indices <- integer()
		#fill with index waa pointers
		i <- 1
		for(k in 1:length(asap3)) {
			x <- asap3[[k]]
			for(f in 1:x$n_indices){
				data$waa_pointer_indices[i] <- data$n_fleets + i
				i <- i + 1
			}
		}
		if(is.null(input$index_names)) input$index_names <- paste0("index_", 1:data$n_indices)
	
	}
  if(!is.null(index_info$index_names)){
    if(is.character(index_info$index_names) & length(index_info$index_names)== data$n_indices) input$index_names <- index_info$index_names
    else stop("index_info$index_names is either not a character vector or its length is not equal to the number of indices.")
  } else {
  	if(is.null(input$index_names) | (length(input$index_names) != data$n_indices)) {
  		input$index_names <- paste0("index_", 1:data$n_indices)
  	}
  }
  for(name in c("index_regions", "index_seasons","waa_pointer_indices")){
	  if(!is.null(index_info[[name]])){
	    if(length(index_info[[name]]) != data$n_indices){
	      stop(paste0("length of index_info$", name, " is not equal to the number of indices.\n"))
	    }
	  	data[[name]] <- index_info[[name]]
	  } else {
	  	if(is.null(data[[name]]) | (length(data[[name]]) != data$n_indices)) {
	  		data[[name]] <- rep(1, data$n_indices)
	  	}
	  }
	}
  for(name in c("units_indices", "units_index_paa")){
	  if(!is.null(index_info[[name]])){
	    if(length(index_info[[name]]) != data$n_indices){
	      stop(paste0("length of index_info$", name, " is not equal to the number of indices.\n"))
	    }
	  	data[[name]] <- index_info[[name]]
	  } else {
	  	if(is.null(data[[name]]) | (length(data[[name]]) != data$n_indices)) {
	  		data[[name]] <- rep(2, data$n_indices)
	  	}
	  }
	}
	vals <- c(10,0.3,100,1,1,data$fracyr_seasons[1]*0.5)
	index_info$agg_index_sigma <- index_info$agg_index_cv

	obs <- c("agg_indices", "agg_index_sigma","index_Neff", "use_indices", "use_index_paa","fracyr_indices")
  for(i in 1:length(obs) ){
	  if(!is.null(index_info[[obs[i]]])){
	  	data[[obs[i]]] <- index_info[[obs[i]]]
	  } else {
	  	if(is.null(data[[obs[i]]]) | !identical(dim(data[[obs[i]]]), as.integer(c(data$n_years_model, data$n_indices)))) {
	  		data[[obs[i]]] <- matrix(vals[i], data$n_years_model, data$n_indices)
	  	}
	  }
	}
  if(!is.null(index_info$selblock_pointer_indices)){
  	data$selblock_pointer_indices <- index_info$selblock_pointer_indices
  } else {
  	if(is.null(data$selblock_pointer_indices) | !identical(dim(data$selblock_pointer_indices), as.integer(c(data$n_years_model, data$n_indices)))) {
  		data$selblock_pointer_indices <- matrix(1:data$n_indices, data$n_years_model, data$n_indices, byrow = TRUE)
	  	if(!is.null(data$selblock_pointer_fleets)) data$selblock_pointer_indices <- data$selblock_pointer_indices + max(data$selblock_pointer_fleets)
  	}
  }

  if(!is.null(index_info$index_paa)) {
  	data$index_paa <- index_info$index_paa
  } else {
		if(is.null(data$index_paa) | !identical(dim(data$index_paa), as.integer(c(data$n_indices, data$n_years_model, data$n_ages)))){
			data$index_paa <- array(1/data$n_ages, c(data$n_indices, data$n_years_model, data$n_ages))
		} 
  }

  if(is.null(input$by_pwi)){
    if(is.null(data$waa)) stop("basic_info argument does not include an array of weight at age. Add that with appropriate dimensions before calling set_index with index_info$waa_pointer_indices.")
    if(any(!(data$waa_pointer_indices %in% 1:dim(data$waa)[1]))){
      stop("some index_info$waa_pointer_indices are outside the number of waa matrices.\n")
    }
  }
  # if(!is.null(index_info$waa_pointer_indices)){
  #   if(is.null(input$by_pwi)){
  #     if(is.null(data$waa)) stop("basic_info argument does not include an array of weight at age. Add that with appropriate dimensions before calling set_index with index_info$waa_pointer_indices.")
  #     if(any(!(index_info$waa_pointer_indices %in% 1:dim(data$waa)[1]))){
  #       stop("some index_info$waa_pointer_indices are outside the number of waa matrices.\n")
  #     }
  #   }
    # if(length(index_info$waa_pointer_indices) != data$n_indices){
    #   stop("length of index_info$waa_pointer_indices is not equal to the number of fleets.\n")
    # }
    # data$waa_pointer_indices <- index_info$waa_pointer_indices
  # }
  if(is.null(index_info$waa_pointer_indices) & input$use_asap3) input$log$indices <- c(input$log$indices, "waa_pointer_indices determined from ASAP file(s). \n")

  ################################################################################
  # for plotting, in years where index is not used set sigma = avg of used years
  tmp <- data$agg_index_sigma
  tmp[data$use_indices == 0] <- NA
  mean_agg_ind_sigma <- apply(tmp, 2, mean, na.rm=T)
  for(i in 1:data$n_indices) data$agg_index_sigma[data$use_indices[,i] == 0,i] <- mean_agg_ind_sigma[i]
  ################################################################################

	data$index_paa[is.na(data$index_paa)] <- 0
  data$agg_index_sigma[which(data$agg_index_sigma < 1e-15)] <- 100  
  data$agg_index_sigma <- sqrt(log(data$agg_index_sigma^2 + 1))

  input$par$log_index_sig_scale <- rep(0, data$n_indices)
  if(!is.null(index_info$initial_index_sd_scale)) input$par$log_index_sig_scale[] <- log(index_info$initial_index_sd_scale)
  input$map$log_index_sig_scale <-factor(rep(NA, data$n_indices))
  if(!is.null(index_info$map_index_sd_scale)) input$map$log_index_sig_scale <- factor(index_info$map_index_sd_scale)
  
  input$data <- data
  if(is.null(input$by_pwi)) { #check whether called by prepare_wham_input
  	input <- set_q(input, input$options$catchability)
  	input <- set_selectivity(input, input$options$selectivity)
		if(any(!input$data$selblock_pointer_indices %in% 1:input$data$n_selblocks)){
			input$log$indices <- c(input$log$indices, 
				paste0("NOTE: Some index selblock_pointers are outside of the number of selectivity blocks currently specified.\n",
				"\tMake sure to run set_selectivity with appropriate options before using the input.\n"))
		}
		input <- set_ecov(input, input$options$ecov) #may need to resize dimensions for catchability if number of indices changes
  	input <- set_age_comp(input, input$options$age_comp)
  	input <- set_osa_obs(input)
  }
  input$log$indices <- c(input$log$indices, paste0("Selectivity blocks for each index:\n",
    paste0(input$index_names, ": ", apply(data$selblock_pointer_indices, 2, function(x) paste0(sort(unique(x)), collapse = ", ")), collapse ="\n"), "\n\n")
  )

 	input$options$index <- index_info
  if(length(input$log$indices))  input$log$indices <- c(
    "--Indices----------------------------------------------------------------------------------------------------------------------------",
  	"\n", 
  	input$log$indices,
    "-------------------------------------------------------------------------------------------------------------------------------------",
  "\n\n")

  if(is.null(input$by_pwi)) message(unlist(input$log$indices, recursive=T))
 return(input)
}