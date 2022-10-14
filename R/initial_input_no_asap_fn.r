
initial_input_no_asap_fn <- function(input, basic_info){
	
	input$years = 1975:2014	
	if(!is.null(basic_info$years)) {
		if(!is.integer(basic_info$years)) stop("basic_info$years has been specified, but it is not an integer vector")
		input$years = basic_info$years
	}
	input$data$n_years_model = length(input$years)

	input$data$n_stocks = 1
	if(!is.null(basic_info$n_stocks)) {
		if(!is.integer(basic_info$n_stocks)) stop("basic_info$n_stocks has been specified, but it is not an integer")
		input$data$n_stocks = basic_info$n_stocks
	}
	
	input$data$n_regions = 1
	if(!is.null(basic_info$n_regions)) {
		if(!is.integer(basic_info$n_regions)) stop("basic_info$n_regions has been specified, but it is not an integer.")
		input$data$n_regions = basic_info$n_regions
	}
	
	input$data$n_seasons = 1
	if(!is.null(basic_info$n_seasons)) {
		if(!is.integer(basic_info$n_seasons)) stop("basic_info$n_seasons has been specified, but it is not an integer.")
		input$data$n_seasons = basic_info$n_seasons
	}
	
	input$data$fracyr_seasons = 1
	if(!is.null(basic_info$fracyr_seasons)) {
		if(!is.numeric(basic_info$fracyr_seasons)) stop("basic_info$fracyr_seasons has been specified, but it is not numeric.")
		if(length(basic_info$fracyr_seasons) != input$data$n_seasons) stop("length of basic_info$fracyr_seasons not equal to n_seasons.")
		input$data$fracyr_seasons = basic_info$fracyr_seasons
		if(length(basic_info$fracyr_seasons) != 1) cat("basic_info$fracyr_seasons implies more than one season. Ensure that 
			input$data$fracyr_SSB, input$data$spawn_seasons, input$data$fracyr_indices, input$data$index_seasons, are specified apropriately. \n ")
	}

	input$data$n_ages = 10
	if(!is.null(basic_info$ages)) {
		if(!is.integer(basic_info$ages)) stop("basic_info$ages has been specified, but it is not an integer vector")
		input$data$n_ages = length(basic_info$ages)
	}

	input$data$fracyr_SSB = matrix(rep(0.25, length(input$years)), input$data$n_years_model, input$data$n_stocks)
	if(!is.null(basic_info$fracyr_SSB)){
		if(!is.matrix(basic_info$fracyr_SSB) | NROW(basic_info$fracyr_SSB) != input$data$n_years_model | NCOL(basic_info$fracyr_SSB) != input$data$n_stocks){
			stop("basic_info$fracyr_SSB has been specified, but it's not a matrix or its dimensions are incorrect")	
		} 
		input$data$fracyr_SSB[] = basic_info$fracyr_SSB
	}
	input$data$mature = array(NA, dim = c(1,input$data$n_years_model, input$data$n_ages))
	input$data$mature[1,,] = t(matrix(1/(1 + exp(-1*(1:input$data$n_ages - input$data$n_ages/2))), input$data$n_ages, length(input$years)))
	if(!is.null(basic_info$maturity)){
		if(!(length(basic_info$maturity) %in% c(1,input$data$n_stocks * input$data$n_ages*length(input$years)))) stop("basic_info$mature has been specified, but it's length is not 1 or length(ages)*length(years)")
		else input$data$mature[] = basic_info$maturity
	}		
	input$data$Fbar_ages = 1:input$data$n_ages
	if(!is.null(basic_info$Fbar_ages)) {
		if(!is.integer(basic_info$Fbar_ages)) stop("basic_info$Fbar_ages has been specified, but it is not an integer vector")
		else input$data$Fbar_ages = basic_info$Fbar_ages
	}

  #indicate which regions each stock and age can be in on Jan 1 and which regions spawning occurs and recruitment originates on Jan 1.
	if(data$n_regions!=data$n_stocks & is.null(basic_info$NAA_where)) {
		stop("basic_info$NAA_where must be provided when data$n_regions is not equal to data$n_stocks. \n")
	}
  data$NAA_where = array(0,dim = c(data$n_stocks,data$n_regions,data$n_ages))
	for(s in 1:data$n_stocks) data$NAA_where[s,s,] = 1
  
	if(!is.null(basic_info$NAA_where)) {
    if(!all(dim(NAA_re$NAA_where) == dim(data$NAA_where))) stop("Dimensions of basic_info$NAA_where are not correct. \n")
    data$NAA_where[] = NAA_re$NAA_where
  }
	for(s in 1:n_stocks) {
		sr = which(data$NAA_where[s,,1] == 1)
		if(length(sr)>1) stop("NAA_where needs to have age 1 occuring in only 1 region on Jan 1.")
	} 
	data$spawn_regions = rep(1:n_stocks, data$n_stocks)
  if(is.null(basic_info$spawn_regions)) {
		cat("input$data$spawn_regions will be defined from input$data$NAA_where for age 1.\n")
		for(s in 1:n_stocks) {
			sr = which(data$NAA_where[s,,1] == 1)
			data$spawn_regions[s] = sr
		}
  } else {
    if(!all(length(NAA_re$spawn_regions) == data$n_stocks)) stop("length of NAA_re$spawn_regions is not equal to data$n_stocks. \n")
    if(max(NAA_re$spawn_regions) > data$n_regions)) stop("maximum value of NAA_re$spawn_regions is greater than data$n_regions. \n")
		for(s in 1:n_stocks) {
			sr = which(data$NAA_where[s,,1] == 1)
			if(sr != basic_info$spawn_regions[s]) stop("basic_info$spawn_regions not consistent with input$dataNAA_where.")
			data$spawn_regions[] = basic_info$spawn_regions
		}
  }
	
	#print(input)
	return(input)
}
