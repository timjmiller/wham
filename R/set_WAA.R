set_WAA <- function(input, waa_info = NULL) {
	data <- input$data
  asap3 <- input$asap3
	# Weight-at-age
	data$waa <- array(NA,dim = c(data$n_fleets + data$n_regions + data$n_indices + data$n_stocks, data$n_years_model, data$n_ages))
	if(!is.null(asap3)) {
		waa_temp = list()
		n_waa_total = sum(sapply(asap3, function(x) x$n_fleets)) + 2 * length(asap3) + sum(sapply(asap3, function(x) length(x$index_WAA_pointers)))
		data$waa <- array(NA,dim = c(n_waa_total, data$n_years_model, data$n_ages))
		data$waa_pointer_fleets <- integer()
		data$waa_pointer_ssb <- integer()
		data$waa_pointer_indices <- integer()
		data$waa_pointer_totcatch <- integer()
		j = 1
    for(a in 1:length(asap3)) {
			waa_pointer_fleets <- 2*(1:asap3[[a]]$n_fleets) - 1 #wham has no discard data, so remove those WAA matrices
			for(k in 1:length(waa_pointer_fleets)) {
				data$waa[j,,] = asap3[[a]]$WAA_mats[[waa_pointer_fleets[k]]]
				data$waa_pointer_fleets = c(data$waa_pointer_fleets,j)
				j <- j + 1
			}
			waa_pointer_totcatch <- asap3[[a]]$n_fleets*2+1
			data$waa[j,,] = asap3[[a]]$WAA_mats[[waa_pointer_totcatch]]
			data$waa_pointer_totcatch = c(data$waa_pointer_totcatch,j)
			j <- j + 1
			waa_pointer_ssb <- asap3[[a]]$n_fleets*2+3
			data$waa[j,,] = asap3[[a]]$WAA_mats[[waa_pointer_ssb]]
			data$waa_pointer_ssb = c(data$waa_pointer_ssb,j)
			j <- j + 1
			waa_pointer_indices <- asap3[[a]]$index_WAA_pointers
			for(k in 1:length(waa_pointer_indices)) {
				data$waa[j,,] = asap3[[a]]$WAA_mats[[waa_pointer_indices[k]]]
				data$waa_pointer_indices = c(data$waa_pointer_indices,j)
				j <- j + 1
			}
		}
	} else {
  	L = 100*(1-exp(-0.3*(1:data$n_ages - 0)))
    W = rep(exp(-11)*L^3, each = data$n_years_model)
		data$waa = array(W, dim = c(1, data$n_years_model, data$n_ages))
		data$waa_pointer_fleets = rep(1,data$n_fleets)
		data$waa_pointer_totcatch = 1
		data$waa_pointer_ssb = 1
		data$waa_pointer_indices = rep(1,data$n_indices)
	}

  if(!is.null(waa_opts$waa)){
		data$waa = waa_opts$waa
		dim_waa = dim(data$waa)
		if(is.null(waa_opts$waa_pointer_fleets)){
			cat("waa_opts$waa is provided without waa_opts$waa_pointer_fleets so the first waa matrix will be used for all fleets. \n")
			data$waa_pointer_fleets = rep(1,data$n_fleets)	
		}
	}
		
	if(!is.null(waa_opts$waa_pointer_fleets)){
		if(max(waa_opts$waa_pointer_fleets) > dim(data$waa)[1]){
			stop("some waa_opts$waa_pointer_fleets are outside the number of waa matrices.\n")
		}
		if(length(waa_opts$waa_pointer_fleets) != data$n_fleets){
			stop("length of waa_opts$waa_pointer_fleets is not equal to the number of fleets.\n")
		}
		data$waa_pointer_fleets = waa_opts$waa_pointer_fleets
	}
	if(!is.null(waa_opts$waa_pointer_indices)){
		if(max(waa_opts$waa_pointer_indices) > dim(data$waa)[1]){
			stop("some waa_opts$waa_pointer_indices are outside the number of waa matrices.\n")
		}
		if(length(waa_opts$waa_pointer_indices) != data$n_indices){
			stop("length of waa_opts$waa_pointer_indices is not equal to the number of indices.\n")
		}
		data$waa_pointer_indices = waa_opts$waa_pointer_indices
	}
	if(!is.null(waa_opts$waa_pointer_ssb)){
		if(max(waa_opts$waa_pointer_ssb) > dim(data$waa)[1]){
			stop("some waa_opts$waa_pointer_ssb are outside the number of waa matrices.\n")
		}
		if(length(waa_opts$waa_pointer_ssb) != data$n_stocks){
			stop("length of waa_opts$waa_pointer_ssb is not equal to the number of stocks.\n")
		}
		data$waa_pointer_ssb = waa_opts$waa_pointer_ssb
	}
	if(!is.null(waa_opts$waa_pointer_totcatch)){
		if(max(waa_opts$waa_pointer_totcatch) > dim(data$waa)[1]){
			stop("some waa_opts$waa_pointer_totcatch are outside the number of waa matrices.\n")
		}
		if(length(waa_opts$waa_pointer_totcatch) != data$n_regions){
			stop("length of waa_opts$waa_pointer_totcatch is not equal to the number of regions.\n")
		}
		data$waa_pointer_totcatch = waa_opts$waa_pointer_totcatch
	}

  input$data = data
  return(input)
}