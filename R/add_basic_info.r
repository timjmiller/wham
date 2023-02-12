add_basic_info <- function(input, basic_info){
	input$years = 1975:2014
	if(!is.null(basic_info$years)) {
		if(!is.integer(basic_info$years)) stop("basic_info$years has been specified, but it is not an integer vector")
		else input$years = basic_info$years
	}	
	input$data$n_ages = 10
	if(!is.null(basic_info$ages)) {
		if(!is.integer(basic_info$ages)) stop("basic_info$ages has been specified, but it is not an integer vector")
		else {
			input$data$n_ages = length(basic_info$ages)
		}
	}
	input$data$lengths = seq(from = 2, to = 82, by = 4)
	input$data$n_lengths = length(input$data$lengths)
	if(!is.null(basic_info$lengths)) {
		input$data$lengths = basic_info$lengths
		input$data$n_lengths = length(basic_info$lengths)
	}
	input$data$fracyr_SSB = rep(0.25, length(input$years))
	if(!is.null(basic_info$fracyr_SSB)){
		if(!(length(basic_info$fracyr_SSB) %in% c(1,length(input$years)))) stop("basic_info$fracyr_SSB has been specified, but it's length is not 1 or length(years)")
		else input$data$fracyr_SSB[] = basic_info$fracyr_SSB
	}
	#age-based maturity:
	input$data$mature = t(matrix(1, input$data$n_ages, length(input$years))) # Imporant for growth WHAM
	if(!is.null(basic_info[["maturity"]])){
		if(!(length(basic_info$maturity) %in% c(1,input$data$n_ages*length(input$years)))) stop("basic_info$mature has been specified, but it's length is not 1 or length(ages)*length(years)")
		else input$data$mature[] = basic_info$maturity
	}
	#len-based maturity:
	# input$data$mature_len = t(matrix(1, input$data$n_lengths, length(input$years)))
	# if(!is.null(basic_info[["maturity_len"]])){
	# 	if(!(length(basic_info$maturity_len) %in% c(1,input$data$n_lengths*length(input$years)))) stop("basic_info$mature_len has been specified, but it's length is not 1 or n_lengths*length(years)")
	# 	else input$data$mature_len[] = basic_info$maturity_len
	# }
	# Continue...
	input$data$Fbar_ages = 1:input$data$n_ages
	if(!is.null(basic_info$Fbar_ages)) {
		if(!is.integer(basic_info$Fbar_ages)) stop("basic_info$Fbar_ages has been specified, but it is not an integer vector")
		else input$data$Fbar_ages = basic_info$Fbar_ages
	}
	# For growth:
	input$data$age_L1 = 1
	if(!is.null(basic_info$age_L1)) input$data$age_L1 = basic_info$age_L1 # assuming age = 1 at L1
  	if(input$data$age_L1 < 1) stop("'age_L1' cannot be younger than 1")
  	input$data$age_L1_ceil = as.integer(ceiling(input$data$age_L1)) # very important for parametric growth.

	return(input)
}
