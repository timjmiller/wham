set_WAA <- function(input, waa_info = NULL) {
	data <- input$data
  asap3 <- input$asap3
	# Weight-at-age

  if(!is.null(waa_info$waa)){
		data$waa = waa_info$waa
		dim_waa = dim(data$waa)
		if(length(dim_waa) != 3) stop("basic_info$waa must be a 3d array. second index is number of years, third is number of ages.")
	} else {
		if(!is.null(asap3)) {
			# data$waa <- array(NA,dim = c(data$n_fleets + data$n_regions + data$n_indices + data$n_stocks, data$n_years_model, data$n_ages))
			data$waa <- array(NA,dim = c(data$n_fleets + data$n_indices + data$n_stocks, data$n_years_model, data$n_ages))
			data$waa_pointer_ssb <- integer()
			# data$waa_pointer_fleets <- integer()
			# data$waa_pointer_indices <- integer()
			# waa_pointer_totcatch <- integer()
			
			#fill with fleet catch waa
			i <- 1
			for(k in 1:length(asap3)) {
				x <- asap3[[k]]
				for(f in 1:x$n_fleets){
					data$waa[i,,] <- x$WAA_mats[[x$WAA_pointers[2*f-1]]]
			# 		data$waa_pointer_fleets[i] <- i
					i <- i + 1
				}
			}
			# fill with total catch waa (region)
			# for(k in 1:length(asap3)) {
			# 	x <- asap3[[k]]
			# 	data$waa[data$n_fleets + k,,] <- x$WAA_mats[[x$WAA_pointers[2*x$n_fleets+1]]]
			# # 	waa_pointer_totcatch[k] <- data$n_fleets + k
			# }
			#fill with index waa
			i <- 1
			for(k in 1:length(asap3)) {
				x <- asap3[[k]]
				for(f in 1:x$n_indices){
					# data$waa[data$n_fleets + data$n_regions + i,,] <- x$WAA_mats[[x$index_WAA_pointers[f]]]
					data$waa[data$n_fleets + i,,] <- x$WAA_mats[[x$index_WAA_pointers[f]]]
			# 		data$waa_pointer_indices[i] <- data$n_fleets + data$n_regions + i
					i <- i + 1
				}
			}
			#fill with ssb waa (stocks)
			for(k in 1:length(asap3)) {
				x <- asap3[[k]]
				# data$waa[data$n_fleets + data$n_regions + data$n_indices + k,,] <- x$WAA_mats[[x$WAA_pointers[2*x$n_fleets+3]]]
				data$waa[data$n_fleets + data$n_indices + k,,] <- x$WAA_mats[[x$WAA_pointers[2*x$n_fleets+3]]]
				data$waa_pointer_ssb[k] <- data$n_fleets + data$n_indices + k
			}
		} else { #no asap and no waa provided
			L = 100*(1-exp(-0.3*(1:data$n_ages - 0)))
	    W = rep(exp(-11)*L^3, each = data$n_years_model)
			data$waa = array(W, dim = c(1, data$n_years_model, data$n_ages))
			input$log$waa <- c(input$log$waa, "input$data$waa filled with fake data. \n")
		}
	}
	

	input$log$waa <- list()
	if(!is.null(waa_info$waa_pointer_ssb)){
		if(any(!(waa_info$waa_pointer_ssb %in% 1:dim(data$waa)[1]))){
			stop("some basic_info$waa_pointer_ssb are outside the number of waa matrices.\n")
		}
		if(length(waa_info$waa_pointer_ssb) != data$n_stocks){
			stop("length of basic_info$waa_pointer_ssb is not equal to the number of stocks.\n")
		}
		data$waa_pointer_ssb = waa_info$waa_pointer_ssb
	} else{
		if(is.null(asap3)) {
			data$waa_pointer_ssb = rep(1,data$n_stocks)
			input$log$waa <- c(input$log$waa, "basic_info$waa_pointer_ssb was not provided and no asap3 file(s), so the first waa matrix will be used. \n")
		}
	}
	if(is.null(waa_info$waa_pointer_M)){
		input$log$waa <- c(input$log$waa, "basic_info$waa_pointer_M was not provided, so waa_pointer_ssb will be used. \n")
		data$waa_pointer_M = rep(data$waa_pointer_ssb,data$n_stocks)
	} else{
		if(any(!(waa_info$waa_pointer_M %in% 1:dim(data$waa)[1]))){
			stop("some basic_info$waa_pointer_M are outside the number of waa matrices.\n")
		}
		if(length(waa_info$waa_pointer_M) != data$n_stocks){
			stop("length of basic_info$waa_pointer_M is not equal to the number of regions.\n")
		}
		data$waa_pointer_M[] = waa_info$waa_pointer_M
	}
	if(!is.null(data$waa_pointer_fleets)){
    if(any(!(data$waa_pointer_fleets %in% 1:dim(data$waa)[1]))){
      stop("some data$waa_pointer_fleets are outside the number of waa matrices.\n")
    }
	}
	if(!is.null(data$waa_pointer_indices)){
    if(any(!(data$waa_pointer_indices %in% 1:dim(data$waa)[1]))){
    	print(data$waa_pointer_indices)
    	print(dim(data$waa))
      stop("some data$waa_pointer_indices are outside the number of waa matrices.\n")
    }
	}
	if(length(input$log$waa))	input$log$waa <- c("WAA: \n", input$log$waa)

  input$data = data
  input$options$waa <- waa_info
  if(!is_internal_call()) cat(unlist(input$log$waa, recursive=T))
 	return(input)
}