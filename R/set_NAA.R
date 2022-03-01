set_NAA = function(input, NAA_re=NULL)
{

  data = input$data
  par = input$par
  map = input$map
  asap3 = if(is.null(input$asap3)) asap3 = NULL
  else asap3 = input$asap3

  #set up initial NAA
  if(is.null(NAA_re$N1_model)) {
    data$N1_model = 0 #0: just age-specific numbers at age
  } else {
    data$N1_model = NAA_re$N1_model 
  }
  if(is.null(NAA_re$N1_pars)){
    if(data$N1_model == 0){
      if(!is.null(asap3)) par$log_N1_pars = log(asap3$N1_ini) # use N1_ini values from asap3 file
      else par$log_N1_pars = log(exp(10)*exp(-(0:(data$n_ages-1))*0.2))
    }
    if(data$N1_model == 1) par$log_N1_pars = c(10,log(0.1)) # allowed in wham.cpp but no option to set here (must be hard-coded after calling prepare_wham_input)
  }
  else par$log_N1_pars = log(NAA_re$N1_pars)

  if(data$N1_model>1) stop("NAA_re$N1_model can only be 0 or 1 currently")

  # NAA_re options
  if(is.null(NAA_re$sigma)){ # default = SCAA
    data$n_NAA_sigma <- 0
    data$NAA_sigma_pointers <- rep(1,data$n_ages)
  } 
  if(!is.null(NAA_re$sigma)) {
    if(NAA_re$sigma == "rec"){
      data$n_NAA_sigma <- 1
      data$NAA_sigma_pointers <- rep(1,data$n_ages)
    } else {
      if(NAA_re$sigma == "rec+1"){ # default state-space model with two NAA_sigma (one for recruits, one for ages > 1)
        data$n_NAA_sigma <- 2
        data$NAA_sigma_pointers <- c(1,rep(2,data$n_ages-1))
      } else {
        if(length(NAA_re$sigma) != data$n_ages) stop("NAA_re$sigma must either be 'rec' (random effects on recruitment only), 
'rec+1' (random effects on all NAA with ages > 1 sharing sigma_a,
or a vector with length == n.ages specifying which sigma_a to use for each age.")
        if(length(NAA_re$sigma) == data$n_ages){
          if(is.na(unique(NAA_re$sigma))){
            data$n_NAA_sigma <- 0
            data$NAA_sigma_pointers <- rep(1,data$n_ages)            
          } else {
            data$n_NAA_sigma <- max(unique(NAA_re$sigma), na.rm=T)
            data$NAA_sigma_pointers <- NAA_re$sigma
          }
        }
      }
    }
  }
  if(data$recruit_model > 2 & data$n_NAA_sigma == 0) warning("SCAA model specified, yearly recruitment deviations estimated as fixed effects. Stock-recruit function also specified. WHAM will fit the SCAA model but without estimating a stock-recruit function.
This message will not appear if you set recruit_model = 2 (random about mean).")

  # NAA_re pars
  if(data$n_NAA_sigma == 0){
    par$log_NAA_sigma <- 0
    map$log_NAA_sigma <- factor(NA)
  } else {
    par$log_NAA_sigma <- rep(0,data$n_NAA_sigma)
  }
  par$trans_NAA_rho <- c(0,0)
  par$log_NAA = matrix(10, data$n_years_model-1, data$n_ages)
  par$logR_proj <- 0 # will be set by prepare_projection if SCAA
  map$logR_proj <- factor(NA)
  
  # NAA_re and NAA_rho map
  if(!is.null(NAA_re$cor)){
    if(!NAA_re$cor %in% c("iid","ar1_a","ar1_y","2dar1")) stop("NAA_re$cor must be one of 'iid','ar1_a','ar1_y','2dar1'")
  } else {
    NAA_re$cor <- 'iid'
  }
  tmp <- par$trans_NAA_rho
  if(NAA_re$cor %in% c("iid","ar1_y")) tmp[1] = NA 
  if(NAA_re$cor %in% c("iid","ar1_a")) tmp[2] = NA
  ind.notNA <- which(!is.na(tmp))
  tmp[ind.notNA] <- 1:length(ind.notNA)
  map$trans_NAA_rho = factor(tmp)

  tmp <- par$log_NAA
  if(data$n_NAA_sigma < 2) tmp[,-1] <- NA # always estimate Rec devs (col 1), whether random effect or not
  ind.notNA <- which(!is.na(tmp))
  tmp[ind.notNA] <- 1:length(ind.notNA)
  map$log_NAA = factor(tmp)

  #set up recruitment
  if(!is.null(NAA_re$recruit_model)) {
    data$recruit_model = NAA_re$recruit_model #overrides recruit_model argument to wham::prepare_wham_input
    if(data$recruit_model > 1 & length(NAA_re$sigma) == 0) stop("NAA_re$recruit_model > 1 has been specified, but NAA_re$sigma must either be 'rec' (random effects on recruitment only), 
      'rec+1' (random effects on all NAA with ages > 1 sharing sigma_a,
      or a vector with length == n.ages specifying which sigma_a to use for each age.")
  }
  par$mean_rec_pars = numeric(c(0,1,2,2)[data$recruit_model])

  if(!is.null(NAA_re$use_steepness)) data$use_steepness = sum(NAA_re$use_steepness)
  else data$use_steepness = 0 #use regular SR parameterization by default, steepness still can be estimated as derived par.
  
  if(!is.null(NAA_re$recruit_pars)){
    if(data$recruit_model == 2) par$mean_rec_pars[] = log(NAA_re$recruit_pars[1])
    if(data$recruit_model %in% 3:4){
      if(data$use_steepness == 1){
        if(data$recruit_model == 3) par$mean_rec_pars[1] = log(NAA_re$recruit_pars[1] - 0.2) - log(1-NAA_re$recruit_pars[1])
        if(data$recruit_model == 4) par$mean_rec_pars[1] = log(NAA_re$recruit_pars[1] - 0.2)
        par$mean_rec_pars[2] = log(NAA_re$recruit_pars[2])
      }
      else par$mean_rec_pars[] = log(NAA_re$recruit_pars)
    }
  }
  else{
    par$mean_rec_pars[] = 0	
    if(data$recruit_model==2) {
  		if(!is.null(asap3)) par$mean_rec_pars = log(asap3$N1_ini[1]) # initialize R0 at initial age-1
  		else par$mean_rec_pars = 10
  	}  
    if(data$recruit_model==4) par$mean_rec_pars[2] = -10
    if(data$recruit_model == 1) map$mean_rec_pars = factor(rep(NA, length(par$mean_rec_pars)))
    if(data$n_NAA_sigma == 0) map$mean_rec_pars = factor(rep(NA, length(par$mean_rec_pars))) #SCAA with recruitments as fixed effects
  }

  input$data = data
  input$par = par
  input$map = map
  return(input)
}