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
	input$data$fracyr_SSB = rep(0.25, length(input$years))
	if(!is.null(basic_info$fracyr_SSB)){
		if(!(length(basic_info$fracyr_SSB) %in% c(1,length(input$years)))) stop("basic_info$fracyr_SSB has been specified, but it's length is not 1 or length(years)")
		else input$data$fracyr_SSB[] = basic_info$fracyr_SSB
	}
	input$data$mature = t(matrix(1/(1 + exp(-1*(1:input$data$n_ages - input$data$n_ages/2))), input$data$n_ages, length(input$years)))
	if(!is.null(basic_info$maturity)){
		if(!(length(basic_info$maturity) %in% c(1,input$data$n_ages*length(input$years)))) stop("basic_info$mature has been specified, but it's length is not 1 or length(ages)*length(years)")
		else input$data$mature[] = basic_info$maturity
	}		
	input$data$Fbar_ages = 1:input$data$n_ages
	if(!is.null(basic_info$Fbar_ages)) {
		if(!is.integer(basic_info$Fbar_ages)) stop("basic_info$Fbar_ages has been specified, but it is not an integer vector")
		else input$data$Fbar_ages = basic_info$Fbar_ages
	}
	#print(input)
	return(input)
}
