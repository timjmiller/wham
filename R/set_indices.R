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
set_indices = function(input, index_info=NULL) {
	data = input$data
  asap3 = input$asap3
  input$log$indices <- list()
	if(is.null(asap3)) {
	  data$n_indices = 1
	} else {
		input$index_names <- NULL
    for(i in 1:length(asap3)) {
			which_indices <- which(asap3[[i]]$use_index ==1)
			asap3[[i]]$n_indices = length(which_indices)
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
    n_indices_per_region = sapply(asap3, function(x) x$n_indices)
    data$n_indices = sum(n_indices_per_region)
	} 
	if(!is.null(index_info$n_indices)) data$n_indices = index_info$n_indices
	if(is.null(input$index_names)) input$index_names <- paste0("index_", 1:data$n_indices)

  data$index_regions = rep(1, data$n_indices)
  data$index_seasons = rep(1, data$n_indices)
	
	data$agg_indices = data$agg_index_sigma = data$index_Neff = matrix(NA, data$n_years_model, data$n_indices)
	data$index_paa = array(NA, dim = c(data$n_indices, data$n_years_model, data$n_ages))
	data$use_indices = matrix(1, data$n_years_model, data$n_indices)
	data$use_index_paa = matrix(1, data$n_years_model, data$n_indices)
  data$selblock_pointer_indices = matrix(0, data$n_years_model, data$n_indices)
  data$units_indices = rep(2,data$n_indices)
  data$units_index_paa = rep(2,data$n_indices)
  data$fracyr_indices = matrix(data$fracyr_seasons[1]*0.5, data$n_years_model, data$n_indices)

	if(!is.null(asap3)) {
    k <- 1
    for(i in 1:length(asap3)) {
      for(j in 1:asap3[[i]]$n_indices) {
        data$index_regions[k] = i #each asap file is a separate region
		  	data$units_indices[k] <- asap3[[i]]$survey_index_units[j]
		  	tmp = (asap3[[i]]$survey_month[j]-1)/12 #make sure that this is right
				int_starts <- cumsum(c(0,data$fracyr_seasons))
		  	ind = max(which(int_starts <= tmp))
		  	data$index_seasons[k] = ind
		  	data$fracyr_indices[,k] = tmp - int_starts[ind]
		  	
				data$agg_indices[,k] = asap3[[i]]$IAA_mats[[j]][,2]
	    	for(y in 1:data$n_years_model) if(asap3[[i]]$IAA_mats[[j]][y,2] < 1e-15) data$use_indices[y,k] = 0
	  		data$agg_index_sigma[,k] = asap3[[i]]$IAA_mats[[j]][,3]
		    
		    temp = asap3[[i]]$IAA_mats[[j]][,3 + 1:data$n_ages]
	    	temp[which(is.na(temp))] = 0
	    	temp[which(temp<0)] = 0
	    	data$index_paa[k,,] = temp/apply(temp,1,sum) #all 0s will make NaN
				if(asap3[[i]]$use_survey_acomp[j] != 1){
					data$use_index_paa[,k] = 0
				} else {
					for(y in 1:data$n_years_model) {
						flag <- asap3[[i]]$IAA_mats[[j]][y,4 + data$n_ages] < 1e-15 | sum(data$index_paa[k,y,] > 1e-15) < 2 | any(is.na(data$index_paa[k,y,]))
						if(flag) data$use_index_paa[y,k] = 0
					}
				}
				data$units_index_paa[k] <- asap3[[i]]$survey_acomp_units[j]
				data$index_Neff[,k] = asap3[[i]]$IAA_mats[[j]][,4 + data$n_ages]
        data$selblock_pointer_indices[,k] = max(data$selblock_pointer_fleets) + k #set_catch already called
        k <- k + 1
			}
		}
	}
	else {
		data$agg_indices[] = 10
		data$agg_index_sigma[] = 0.3
		data$index_paa[] = 1/data$n_ages
		data$index_Neff[] = 100
		data$selblock_pointer_indices[] = rep(1:data$n_indices, each = data$n_years_model) + max(data$selblock_pointer_fleets)
    input$index_names <- paste0("Index ", 1:data$n_indices)
	}
	if(!is.null(index_info$use_indices)) data$use_indices[] = index_info$use_indices
	if(!is.null(index_info$use_index_paa)) data$use_index_paa[] = index_info$use_index_paa
	if(!is.null(index_info$units_indices)) data$units_indices[] = index_info$units_indices
	if(!is.null(index_info$fracyr_indices)) data$fracyr_indices[] = index_info$fracyr_indices
	if(!is.null(index_info$agg_indices)) data$agg_indices[] = index_info$agg_indices
	if(!is.null(index_info$agg_index_cv)) data$agg_index_sigma[] = index_info$agg_index_cv
	if(!is.null(index_info$index_paa)) data$index_paa[] = index_info$index_paa
	if(!is.null(index_info$units_index_paa)) data$units_index_paa[] = index_info$units_index_paa
	if(!is.null(index_info$index_Neff)) data$index_Neff[] = index_info$index_Neff
  
  if(!is.null(index_info$waa_pointer_indices)){
    if(!is_internal_call()){
      if(is.null(data$waa)) stop("basic_info argument does not include an array of weight at age. Add that with appropriate dimensions before calling set_index with index_info$waa_pointer_indices.")
      if(any(!(index_info$waa_pointer_indices %in% 1:dim(data$waa)[1]))){
        stop("some index_info$waa_pointer_indices are outside the number of waa matrices.\n")
      }
    }
    if(length(index_info$waa_pointer_indices) != data$n_indices){
      stop("length of index_info$waa_pointer_indices is not equal to the number of fleets.\n")
    }
    data$waa_pointer_indices <- index_info$waa_pointer_indices
  } else{
		if(!is.null(asap3)) {
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
      input$log$indices <- c(input$log$indices, "waa_pointer_indices determined from ASAP file(s). \n")
    } else{ #no asap and no waa_pointer provided
      input$log$indices <- c(input$log$indices, "index_info$waa_pointer_indices was not provided, so the first waa matrix will be used for all indices. \n")
      data$waa_pointer_indices <- rep(1,data$n_indices)
    }
	}

	if(!is.null(index_info$selblock_pointer_indices)) data$selblock_pointer_indices[] = index_info$selblock_pointer_indices
	if(!is.null(index_info$index_seasons)) data$index_seasons[] = index_info$index_seasons
	if(!is.null(index_info$index_regions)) data$index_regions[] = index_info$index_regions

  ################################################################################
  # for plotting, in years where index is not used set sigma = avg of used years
  tmp <- data$agg_index_sigma
  tmp[data$use_indices == 0] = NA
  mean_agg_ind_sigma <- apply(tmp, 2, mean, na.rm=T)
  for(i in 1:data$n_indices) data$agg_index_sigma[data$use_indices[,i] == 0,i] = mean_agg_ind_sigma[i]
  ################################################################################

	data$index_paa[is.na(data$index_paa)] = 0
  data$agg_index_sigma[which(data$agg_index_sigma < 1e-15)] = 100  
  data$agg_index_sigma = sqrt(log(data$agg_index_sigma^2 + 1))

  input$par$log_index_sig_scale = rep(0, data$n_indices)
  if(!is.null(index_info$initial_index_sd_scale)) input$par$log_index_sig_scale[] <- log(index_info$initial_index_sd_scale)
  input$map$log_index_sig_scale <-rep(NA, data$n_indices)
  if(!is.null(index_info$map_index_sd_scale)) input$map$log_index_sig_scale <- index_info$map_index_sd_scale
  input$map$log_index_sig_scale = factor(input$map$log_index_sig_scale)
  
  input$asap3 <- asap3

  input$data = data
  if(length(input$log$indices))  input$log$indices <- c("Indices: \n", input$log$indices)
 	input$options$index <- index_info
  if(!is_internal_call()) { #check whether called by prepare_wham_input
  	input <- set_selectivity(input, input$options$selectivity)
  	input <- set_age_comp(input, input$options$age_comp)
  	input <- set_osa_obs(input)
    cat(unlist(input$log$indices, recursive=T))
  }
  return(input)
}