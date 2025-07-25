#' Specify catch selectivity blocks and aggregate and age composition observations for catch
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{prepare_wham_input}})
#' @param catch_info (optional) list specifying various aspects about catch by fleet (see details)
#' 
#' \code{catch_info} specifies observations, and various configuration options for fleet-specific catch observations and will overwrite attributes specified in the ASAP data file.
#' If \code{NULL}, all settings from the ASAP data file or basic_info are used.
#' \code{catch_info} is a list with any of the following entries:
#'   \describe{
#'     \item{$n_fleets}{number of fleets}
#'     \item{$fleet_regions}{vector (n_fleets) of regions where each fleet operates.}
#'     \item{$fleet_seasons}{matrix (n_fleets x n_seasons) of 0/1 values flagging which seasons each fleet operates.}
#'     \item{$fleet_names}{character vector (n_fleets) of names for fleets. Used for naming results in plots and tables.}
#'     \item{$agg_catch}{matrix (n_years_model x n_fleets) of annual aggregate catches by fleet.}
#'     \item{$agg_catch_cv}{matrix (n_years_model x n_fleets) of CVs for annual aggregate catches by fleet.}
#'     \item{$catch_paa}{array (n_fleets x n_years_model x n_ages) of annual catch proportions at age by fleet.}
#'     \item{$use_catch_paa}{matrix (n_years_model x n_fleets) of 0/1 values flagging whether to use proportions at age observations.}
#'     \item{$catch_Neff}{matrix (n_years_model x n_fleets) of effective sample sizes for proportions at age observations.}
#'     \item{$waa_pointer_fleets}{vector (n_fleets) of itegers indicated waa to use for each fleet.}
#'     \item{$selblock_pointer_fleets}{matrix (n_years_model x n_fleets) of itegers indicated selblocks to use.}
#'     \item{$initial_catch_sd_scale}{vector (n_fleets) of scalar multipliers of annual log-observation standard deviation. Default = 1.}
#'     \item{$map_catch_sd_scale}{integer vector (n_fleets) specifying which sd scalar parameters to fix. Use \code{NA} to fix a parameter and integers to estimate. 
#'         Use the same integer for multiple fleets to estimate a shared scalar parameter.}
#'   }
#'
#' @return a named list with same elements as the input provided with catch observations and fleet options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3)
#` newcatch <- matrix(500, input$data$n_years_model, input$data$n_fleets)
#' input <- set_catch(input, catch_info = list(agg_catch = newcatch)) #constant catch of 500 mt
#' }
#'
#' @export
set_catch <- function(input, catch_info= NULL) {
  if(is.null(input[["use_asap3"]])){
    input[["use_asap3"]] <- TRUE
    if(is.null(input[["asap3"]])) input[["use_asap3"]] <- FALSE
  }
  data <- input[["data"]]
  if(is.null(data)) data <- list()
  asap3 <- input$asap3
  input$log$catch <- list()
  if(!input[["use_asap3"]]) {
  # if(is.null(asap3)){
    data$n_fleets <- 1
  } else {
    for(i in 1:length(asap3)) asap3[[i]]$use_catch_acomp <- rep(1,asap3[[i]]$n_fleets) #default is to use age comp for catch
    n_fleets_per_region <- sapply(asap3, function(x) x$n_fleets)
    data$n_fleets <- sum(n_fleets_per_region)
		for(i in 1:length(asap3)) if(!is.null(asap3[[i]]$fleet.names)) input$fleet_names <- c(input$fleet_names, asap3[[i]]$fleet.names)
  }
  if(!is.null(catch_info$n_fleets)) data$n_fleets <- catch_info$n_fleets 
  if(!is.null(catch_info$fleet_names)) {
    if(is.character(catch_info$fleet_names) & length(catch_info$fleet_names)== data$n_fleets) input$fleet_names <- catch_info$fleet_names
    else stop("catch_info$fleet_names is either not a character vector or its length is not equal to the number of fleets.")
    # if(length(unique(input$fleet_names)) != data$n_fleets) stop("Some catch_info$fleet_names are repeated. Provide unique names.")
  }
  if(is.null(input$fleet_names)) input$fleet_names <- paste0("fleet_", 1:data$n_fleets)

  data$fleet_regions <- rep(1, data$n_fleets)
  data$fleet_seasons <- matrix(1, data$n_fleets, data$n_seasons)

	data$agg_catch <- data$agg_catch_sigma <- data$catch_Neff <- matrix(NA, data$n_years_model, data$n_fleets)
  data$catch_paa <- array(NA, dim = c(data$n_fleets, data$n_years_model, data$n_ages))
  data$use_agg_catch <- matrix(1, data$n_years_model, data$n_fleets)
  data$use_catch_paa <- matrix(1, data$n_years_model, data$n_fleets)
  data$selblock_pointer_fleets <- matrix(0, data$n_years_model, data$n_fleets)
	
  if(input[["use_asap3"]]) {
  # if(!is.null(asap3)) {
    #data$n_fleets <- asap3$n_fleets
    k <- 1
    cum_n_selblocks <- 0
    for(i in 1:length(asap3)) {
      asap3[[i]]$use_catch_acomp <- rep(1,asap3[[i]]$n_fleets) #default is to use age comp for catch
      allselblocks_i <- sort(unique(unlist(asap3[[i]]$sel_block_assign)))
      for(j in 1:asap3[[i]]$n_fleets) {
        data$fleet_regions[k] <- i #each asap file is a separate region
        data$agg_catch[,k] <- asap3[[i]]$CAA_mats[[j]][,data$n_ages + 1]
        data$agg_catch_sigma[,k] <- asap3[[i]]$catch_cv[,j]
        temp <- asap3[[i]]$CAA_mats[[j]][,1:data$n_ages]
        temp[which(is.na(temp))] <- 0
        temp[which(temp<0)] <- 0
        data$catch_paa[k,,] <- temp/apply(temp,1,sum)
        for(y in 1:data$n_years_model) if(asap3[[i]]$CAA_mats[[j]][y,data$n_ages+1] < 1e-15) data$use_agg_catch[y,k] <- 0
        if(asap3[[i]]$use_catch_acomp[j] != 1){
          data$use_catch_paa[,k] <- 0
        } else { # use catch paa in at least some years - not necessarily all, have to go through year by year
          for(y in 1:data$n_years_model){
            if(is.na(sum(data$catch_paa[k,y,] > 1e-15))){ # handle negative or NA paa
              data$use_catch_paa[y,k] <- 0
            } else {
              if(asap3[[i]]$catch_Neff[y,j] < 1e-15 | sum(data$catch_paa[k,y,] > 1e-15)<2) data$use_catch_paa[y,k] <- 0
            }
          } 
        }
        data$catch_Neff[,k] <- asap3[[i]]$catch_Neff[,j]
        temp <- asap3[[i]]$sel_block_assign[[j]]
        temp <- match(temp,allselblocks_i) #index unique values
        data$selblock_pointer_fleets[,k] <- cum_n_selblocks + temp #max grows each time
        k <- k + 1
      }
      cum_n_selblocks <- cum_n_selblocks + length(allselblocks_i)
    }
    data$agg_catch_sigma[which(data$agg_catch_sigma < 1e-15)] <- 100
    data$agg_catch_sigma <- sqrt(log(data$agg_catch_sigma^2 + 1))
  }
  else
  {
    #data$n_fleets <- 1
    data$agg_catch[] <- 1000  	
  	data$catch_paa[] <- 1/data$n_ages
		data$agg_catch_sigma[] <- sqrt(log(0.1^2 + 1))
    data$catch_Neff[] <- 200	  
    for(i in 1:data$n_fleets) for(y in 1:data$n_years_model){ 
      if(data$catch_Neff[y,i] < 1e-15 | sum(data$catch_paa[i,y,] > 1e-15)<2 | any(is.na(data$catch_paa[i,y,]))) data$use_catch_paa[y,i] <- 0
    }
    data$selblock_pointer_fleets[] <- rep(1:data$n_fleets, each = data$n_years_model)
    input$fleet_names <- paste0("Fleet ", 1:data$n_fleets)
  }

  if(!is.null(catch_info$agg_catch)) data$agg_catch[] <- catch_info$agg_catch
  if(!is.null(catch_info$catch_paa)) data$catch_paa[] <- catch_info$catch_paa
  if(!is.null(catch_info$agg_catch_cv)) data$agg_catch_sigma[] <- sqrt((log(catch_info$agg_catch_cv^2 + 1)))
  if(!is.null(catch_info$catch_Neff)) data$catch_Neff[] <- catch_info$catch_Neff
  if(!is.null(catch_info$use_catch_paa)) data$use_catch_paa[] <- catch_info$use_catch_paa
  if(!is.null(catch_info$waa_pointer_fleets)){
    if(is.null(input$by_pwi)){
      if(is.null(data$waa)) stop("basic_info argument does not include an array of weight at age. Add that with appropriate dimensions before calling set_catch with catch_info$waa_pointer_fleets.")
      if(any(!(catch_info$waa_pointer_fleets %in% 1:dim(data$waa)[1]))){
        stop("some catch_info$waa_pointer_fleets are outside the number of waa matrices.\n")
      }
    }
    if(length(catch_info$waa_pointer_fleets) != data$n_fleets){
      stop("length of catch_info$waa_pointer_fleets is not equal to the number of fleets.\n")
    }
    data$waa_pointer_fleets <- catch_info$waa_pointer_fleets
  } else{
    # if(!is.null(asap3)) {
    if(input[["use_asap3"]]) {
      data$waa_pointer_fleets <- integer()
      #fill with fleet catch waa pointers
      i <- 1
      for(k in 1:length(asap3)) {
        x <- asap3[[k]]
        for(f in 1:x$n_fleets){
          data$waa_pointer_fleets[i] <- i
          i <- i + 1
        }
      }
      input$log$catch <- c(input$log$catch, "waa_pointer_fleets determined from ASAP file(s). \n")
    } else{ #no asap and no waa_pointer provided
      input$log$catch <- c(input$log$catch, "catch_info$waa_pointer_fleets was not provided, so the first waa matrix will be used for all fleets. \n")
      data$waa_pointer_fleets <- rep(1,data$n_fleets)
    }
  }
  if(!is.null(catch_info$selblock_pointer_fleets)) data$selblock_pointer_fleets[] <- catch_info$selblock_pointer_fleets
  if(!is.null(catch_info$fleet_seasons)) data$fleet_seasons[] <- catch_info$fleet_seasons
  if(!is.null(catch_info$fleet_regions)) data$fleet_regions[] <- catch_info$fleet_regions

  data$catch_paa[is.na(data$catch_paa)] <- 0

  input$par$log_catch_sig_scale <- rep(0, data$n_fleets)
  if(!is.null(catch_info$initial_catch_sd_scale)) input$par$log_catch_sig_scale[] <- log(catch_info$initial_catch_sd_scale)
  input$map$log_catch_sig_scale <-rep(NA, data$n_fleets)
  if(!is.null(catch_info$map_catch_sd_scale)) input$map$log_catch_sig_scale <- catch_info$map_catch_sd_scale
  input$map$log_catch_sig_scale <- factor(input$map$log_catch_sig_scale)

  input$data <- data
  input$asap3 <- asap3

  if(is.null(input$by_pwi)) { #check whether called by prepare_wham_input
     input <- set_F(input, input[["options"]][["F"]])
     input <- set_selectivity(input, input$options$selectivity)
    if(any(!input[["data"]][["selblock_pointer_fleets"]] %in% 1:input[["data"]][["n_selblocks"]])){
      input[["log"]][["catch"]] <- c(input[["log"]][["catch"]], 
        paste0("NOTE: Some fleet selblock_pointers are outside of the number of selectivity blocks currently specified.\n",
        "\tMake sure to run set_selectivity with appropriate options before using the input.\n"))
    }
    input <- set_age_comp(input, input$options$age_comp)
    input <- set_osa_obs(input)
  }
  input$log$catch <- c(input$log$catch, paste0("Selectivity blocks for each fleet:\n",
    paste0(input$fleet_names, ": ", apply(data$selblock_pointer_fleets, 2, function(x) paste0(sort(unique(x)), collapse = ", ")), collapse ="\n"), "\n\n")
  )
  input$options$catch <- catch_info
  if(length(input$log$catch))  input$log$catch <- c(
    "--Catch------------------------------------------------------------------------------------------------------------------------------",
    " \n", 
    input$log$catch,
    "-------------------------------------------------------------------------------------------------------------------------------------",
  "\n\n")
 
  if(is.null(input$by_pwi)) message(unlist(input$log$catch, recursive=T))
  return(input)
}