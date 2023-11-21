set_WAA <- function(input, waa_info = NULL) {
	data <- input$data
  asap3 <- input$asap3
	# Weight-at-age
	data$waa <- array(NA,dim = c(data$n_fleets + data$n_regions + data$n_indices + data$n_stocks, data$n_years_model, data$n_ages))
	if(!is.null(asap3)) {
		# waa_temp = list()
		data$waa_pointer_fleets <- integer()
		data$waa_pointer_ssb <- integer()
		data$waa_pointer_indices <- integer()
		data$waa_pointer_totcatch <- integer()
		#fill with fleet catch waa
		i <- 1
		for(k in 1:length(asap3)) {
			x <- asap3[[k]]
			for(f in 1:x$n_fleets){
				data$waa[i,,] <- x$WAA_mats[[x$WAA_pointers[2*f-1]]]
				data$waa_pointer_fleets[i] <- i
				i <- i + 1
			}
		}
		#fill with total catch waa (region)
		for(k in 1:length(asap3)) {
			x <- asap3[[k]]
			data$waa[data$n_fleets + k,,] <- x$WAA_mats[[x$WAA_pointers[2*x$n_fleets+1]]]
			data$waa_pointer_totcatch[k] <- data$n_fleets + k
		}
		#fill with index waa
		i <- 1
		for(k in 1:length(asap3)) {
			x <- asap3[[k]]
			for(f in 1:x$n_indices){
				data$waa[data$n_fleets + data$n_regions + i,,] <- x$WAA_mats[[x$index_WAA_pointers[f]]]
				data$waa_pointer_indices[i] <- data$n_fleets + data$n_regions + i
				i <- i + 1
			}
		}
		#fill with ssb waa (stocks)
		for(k in 1:length(asap3)) {
			x <- asap3[[k]]
			data$waa[data$n_fleets + data$n_regions + data$n_indices + k,,] <- x$WAA_mats[[x$WAA_pointers[2*x$n_fleets+3]]]
			data$waa_pointer_ssb[k] <- data$n_fleets + data$n_regions + data$n_indices + k
		}
		# n_waa_total = sapply(asap3, function(x){
		# 	i <- c(seq(1,(x$n_fleets+1)*2-1,2),(x$n_fleets+1)*2 + 1:2)
	 #  	WAA_pointers <- x$WAA_pointers[i] #wham has no discard data, so remove those WAA matrices
		# 	WAA_pointers <- c(WAA_pointers,x$index_WAA_pointers) #need any for indices too
		# 	length(unique(WAA_pointers))
		# })
		#print(sum(n_waa_total))

		# data$waa_pointer_fleets <- integer()
		# data$waa_pointer_ssb <- integer()
		# data$waa_pointer_indices <- integer()
		# data$waa_pointer_totcatch <- integer()
		# data$waa <- array(NA,dim = c(sum(n_waa_total), data$n_years_model, data$n_ages))
		# j = 1
		# n_waa <- 0
  #   for(a in 1:length(asap3)) {
		# 	#fleet 1 catch, fleet 2 catch, ..., fleet n catch, totcatch, ssb
		# 	x <- asap3[[a]]
		# 	i <- c(seq(1,(x$n_fleets+1)*2-1,2),(x$n_fleets+1)*2 + 2)
	 #  	WAA_pointers <- c(x$WAA_pointers[i],x$index_WAA_pointers) #wham has no discard data, so remove those WAA matrices
		# 	for(k in unique(WAA_pointers)){ #just retain the needed WAA matrices
		# 		data$waa[j,,] <- x$WAA_mats[[k]]  #note order is changing
		# 		j <- j + 1
		# 	}
		# 	#print(WAA_pointers)
		# 	new_pointer <- n_waa + (1:length(unique(WAA_pointers)))
		# 	#print(new_pointer)
		# 	new_pointer <- new_pointer[match(WAA_pointers,unique(WAA_pointers))]
		# 	#print(new_pointer)
		# 	data$waa_pointer_fleets <- c(data$waa_pointer_fleets,new_pointer[1:x$n_fleets])
		# 	#print(data$waa_pointer_fleets)
		# 	data$waa_pointer_totcatch = c(data$waa_pointer_totcatch,new_pointer[x$n_fleets+1])
		# 	#print(data$waa_pointer_totcatch)
		# 	data$waa_pointer_ssb = c(data$waa_pointer_ssb,new_pointer[x$n_fleets+2])
		# 	data$waa_pointer_indices = c(data$waa_pointer_indices,new_pointer[x$n_fleets+2 + 1:x$n_indices])
		# 	#print(data$waa_pointer_indices)
		# 	n_waa <- n_waa + n_waa_total[a]
		# }
	} else {
  	L = 100*(1-exp(-0.3*(1:data$n_ages - 0)))
    W = rep(exp(-11)*L^3, each = data$n_years_model)
		data$waa = array(W, dim = c(1, data$n_years_model, data$n_ages))
		data$waa_pointer_fleets = rep(1,data$n_fleets)
		data$waa_pointer_totcatch = 1
		data$waa_pointer_ssb = 1
		data$waa_pointer_indices = rep(1,data$n_indices)
	}
	input$log$waa <- list()
  if(!is.null(waa_info$waa)){
		data$waa = waa_info$waa
		dim_waa = dim(data$waa)
		if(length(dim_waa) != 3) stop("waa_info$waa must be a 3d array. second index is number of years, third is number of ages.")
		if(is.null(waa_info$waa_pointer_fleets)){
			input$log$waa <- c(input$log$waa, "waa_info$waa is provided without waa_info$waa_pointer_fleets so the first waa matrix will be used for all fleets. \n")
			data$waa_pointer_fleets = rep(1,data$n_fleets)
		}
		if(is.null(waa_info$waa_pointer_indices)){
			input$log$waa <- c(input$log$waa, "waa_info$waa is provided without waa_info$waa_pointer_totcatch so the first waa matrix will be used. \n")
			data$waa_pointer_indices = rep(1,data$n_indices)
		}
		if(is.null(waa_info$waa_pointer_ssb)){
			input$log$waa <- c(input$log$waa, "waa_info$waa is provided without waa_info$waa_pointer_ssb so the first waa matrix will be used. \n")
			data$waa_pointer_ssb = rep(1,data$n_stocks)
		}
		if(is.null(waa_info$waa_pointer_M)){
			input$log$waa <- c(input$log$waa, "waa_info$waa is provided without waa_info$waa_pointer_M so the first waa matrix will be used. \n")
			data$waa_pointer_M = rep(1,data$n_stocks)
		}
	}
		
	if(!is.null(waa_info$waa_pointer_fleets)){
		if(any(!(waa_info$waa_pointer_fleets %in% 1:dim(data$waa)[1]))){
			stop("some waa_info$waa_pointer_fleets are outside the number of waa matrices.\n")
		}
		if(length(waa_info$waa_pointer_fleets) != data$n_fleets){
			stop("length of waa_info$waa_pointer_fleets is not equal to the number of fleets.\n")
		}
		data$waa_pointer_fleets = waa_info$waa_pointer_fleets
	}
	if(!is.null(waa_info$waa_pointer_indices)){
		if(any(!(waa_info$waa_pointer_indices %in% 1:dim(data$waa)[1]))){
			stop("some waa_info$waa_pointer_indices are outside the number of waa matrices.\n")
		}
		if(length(waa_info$waa_pointer_indices) != data$n_indices){
			stop("length of waa_info$waa_pointer_indices is not equal to the number of indices.\n")
		}
		data$waa_pointer_indices = waa_info$waa_pointer_indices
	}
	if(!is.null(waa_info$waa_pointer_ssb)){
		if(any(!(waa_info$waa_pointer_ssb %in% 1:dim(data$waa)[1]))){
			stop("some waa_info$waa_pointer_ssb are outside the number of waa matrices.\n")
		}
		if(length(waa_info$waa_pointer_ssb) != data$n_stocks){
			stop("length of waa_info$waa_pointer_ssb is not equal to the number of stocks.\n")
		}
		data$waa_pointer_ssb = waa_info$waa_pointer_ssb
	}
	if(!is.null(waa_info$waa_pointer_totcatch)){
		if(any(!(waa_info$waa_pointer_totcatch %in% 1:dim(data$waa)[1]))){
			stop("some waa_info$waa_pointer_totcatch are outside the number of waa matrices.\n")
		}
		if(length(waa_info$waa_pointer_totcatch) != data$n_regions){
			stop("length of waa_info$waa_pointer_totcatch is not equal to the number of regions.\n")
		}
		data$waa_pointer_totcatch = waa_info$waa_pointer_totcatch
	}
	data$waa_pointer_M <- data$waa_pointer_ssb
	if(!is.null(waa_info$waa_pointer_M)){
		if(any(!(waa_info$waa_pointer_M %in% 1:dim(data$waa)[1]))){
			stop("some waa_info$waa_pointer_M are outside the number of waa matrices.\n")
		}
		if(length(waa_info$waa_pointer_M) != data$n_stocks){
			stop("length of waa_info$waa_pointer_totcatch is not equal to the number of regions.\n")
		}
		data$waa_pointer_M[] = waa_info$waa_pointer_M
	}
	if(length(input$log$waa))	input$log$waa <- c("WAA: \n", input$log$waa)

  input$data = data
  input$options$waa <- waa_info
 	return(input)
}