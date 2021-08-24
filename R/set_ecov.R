set_ecov = function(input, ecov)
{
  data = input$data
  par = input$par
  map = input$map
  #random = input$random

  # --------------------------------------------------------------------------------
  # Environmental covariate data
  if(is.null(ecov)){
    data$Ecov_obs <- matrix(1, nrow=1, ncol=1)
    par.Ecov.obs.logsigma <- matrix(-2.3, nrow=1, ncol=1)
    map.Ecov.obs.logsigma <- factor(NA)
    par.Ecov.obs.sigma.par <- matrix(0, nrow=1, ncol=1)
    map.Ecov.obs.sigma.par <- factor(NA)
    data$Ecov_obs_sigma_opt <- 1
    data$n_Ecov <- 1
    data$Ecov_use_obs <- matrix(0, nrow=1, ncol=1)
    data$Ecov_year <- matrix(0, nrow=1, ncol=1)
    data$year1_Ecov <- 0
    data$year1_model <- input$years[1]
    data$Ecov_lag <- 0
    data$Ecov_model <- 0
    data$Ecov_where <- 0
    data$Ecov_how <- 0
    data$Ecov_poly <- 1
    data$n_years_Ecov <- 1
    data$ind_Ecov_out_start <- data$ind_Ecov_out_end <- 0
    data$Ecov_label <- "none"
    data$Ecov_use_re <- matrix(0, nrow=1, ncol=1)
  } else {
    if(class(ecov$mean)[1] == "matrix") {data$Ecov_obs <- ecov$mean} else{
      warning("ecov$mean is not a matrix. Coercing to a matrix...")
      data$Ecov_obs <- as.matrix(ecov$mean)
    }
    if(class(ecov$use_obs)[1] == "matrix"){
      data$Ecov_use_obs <- ecov$use_obs
    } else{
      warning("ecov$use_obs is not a matrix with same dimensions as ecov$mean. Coercing to a matrix...")
      data$Ecov_use_obs <- as.matrix(as.integer(ecov$use_obs))
    }
    if(!identical(dim(data$Ecov_use_obs), dim(data$Ecov_obs))) stop("Dimensions of ecov$use_obs != dimensions of ecov$mean")

    # Handle Ecov sigma options
    data$n_Ecov <- dim(data$Ecov_obs)[2] # num Ecovs
    n_Ecov_obs <- dim(data$Ecov_obs)[1] # num Ecov obs
    if(class(ecov$logsigma)[1] == "matrix"){
      data$Ecov_obs_sigma_opt = 1
      par.Ecov.obs.logsigma <- ecov$logsigma
      if(!identical(dim(par.Ecov.obs.logsigma), dim(data$Ecov_obs))) stop("Dimensions of ecov$mean != dimensions of ecov$logsigma")
    }
    if(class(ecov$logsigma)[1] == 'numeric'){
      data$Ecov_obs_sigma_opt = 1
      print("ecov$logsigma is numeric. Coercing to a matrix...")
      if(length(ecov$logsigma) == data$n_Ecov) par.Ecov.obs.logsigma <- matrix(rep(ecov$logsigma, each=n_Ecov_obs), ncol=data$n_Ecov)
      if(length(ecov$logsigma) == n_Ecov_obs && data$n_Ecov == 1) par.Ecov.obs.logsigma <- matrix(ecov$logsigma, ncol=1)
      if(length(ecov$logsigma) != data$n_Ecov && length(ecov$logsigma) != n_Ecov_obs) stop("ecov$logsigma is numeric but length is not equal to # of ecovs or ecov observations")
    }
    if(class(ecov$logsigma)[1] == 'character'){
      if(ecov$logsigma == 'est_1'){ # estimate 1 Ecov obs sigma for each Ecov
        data$Ecov_obs_sigma_opt = 2
        par.Ecov.obs.logsigma <- matrix(-1.3, nrow=n_Ecov_obs, ncol=data$n_Ecov)
        map.Ecov.obs.logsigma <- matrix(rep(1:data$n_Ecov, each=n_Ecov_obs), ncol=data$n_Ecov)
        par.Ecov.obs.sigma.par <- matrix(-1.3, nrow=2, ncol=data$n_Ecov)
        map.Ecov.obs.sigma.par <- matrix(NA, nrow=2, ncol=data$n_Ecov) # turn off RE pars
      }
      if(ecov$logsigma == 'est_fe'){ # estimate Ecov obs sigma for each Ecov obs
        data$Ecov_obs_sigma_opt = 3
        par.Ecov.obs.logsigma <- matrix(-1.3, nrow=n_Ecov_obs, ncol=data$n_Ecov) # fixed effect inits
        map.Ecov.obs.logsigma <- matrix(1:(n_Ecov_obs*data$n_Ecov), nrow=n_Ecov_obs, ncol=data$n_Ecov) # est all
        par.Ecov.obs.sigma.par <- matrix(-1.3, nrow=2, ncol=data$n_Ecov)
        map.Ecov.obs.sigma.par <- matrix(NA, nrow=2, ncol=data$n_Ecov) # turn off RE pars
      }
      if(ecov$logsigma == 'est_re'){
        data$Ecov_obs_sigma_opt = 4
        par.Ecov.obs.logsigma <- matrix(-1.3, nrow=n_Ecov_obs, ncol=data$n_Ecov) # random effect inits
        map.Ecov.obs.logsigma <- matrix(1:(n_Ecov_obs*data$n_Ecov), nrow=n_Ecov_obs, ncol=data$n_Ecov) # est all
        par.Ecov.obs.sigma.par <- matrix(c(rep(-1.3, data$n_Ecov), rep(-2.3, data$n_Ecov)), ncol=data$n_Ecov, byrow=TRUE) # random effect pars
        map.Ecov.obs.sigma.par <- matrix(1:(2*data$n_Ecov), nrow=2, ncol=data$n_Ecov)
      }
    }
    if(data$Ecov_obs_sigma_opt == 1){ # Ecov sigma given, initialized at given values
      map.Ecov.obs.logsigma <- matrix(NA, nrow=n_Ecov_obs, ncol=data$n_Ecov) # turn off estimation
      par.Ecov.obs.sigma.par <- matrix(-1.3, nrow=2, ncol=data$n_Ecov)
      map.Ecov.obs.sigma.par <- matrix(NA, nrow=2, ncol=data$n_Ecov) # turn off RE pars
    }

    if(length(ecov$year) != n_Ecov_obs) stop("ecov$year is not the same length as # rows in ecov$mean")
    data$Ecov_year <- as.numeric(ecov$year)
    data$year1_Ecov <- ecov$year[1]
    data$year1_model <- input$years[1]
    end_model <- tail(input$years,1)
    end_Ecov <- tail(ecov$year,1)
    if(length(ecov$label) == data$n_Ecov){
      data$Ecov_label <- ecov$label
    } else {
      warning("Number of Ecov labels not equal to number of Ecovs
              Setting Ecov labels = 'Ecov 1', 'Ecov 2', ...")
      data$Ecov_label = paste0("Ecov ",1:data$n_Ecov)
    }

    # # check that Ecov year vector doesn't have missing gaps
    # pad Ecov if it starts after model year1 - max(lag)
    if(data$year1_Ecov > data$year1_model - max(ecov$lag)){
      print("ecov does not start by model year 1 - max(lag). Padding ecov...")
      data$Ecov_obs <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(ecov$lag)), ncol = data$n_Ecov), data$Ecov_obs)
      par.Ecov.obs.logsigma <- rbind(matrix(-1.3, nrow = data$year1_Ecov-(data$year1_model-max(ecov$lag)), ncol = data$n_Ecov), par.Ecov.obs.logsigma)
      map.Ecov.obs.logsigma <- rbind(matrix(NA, nrow = data$year1_Ecov-(data$year1_model-max(ecov$lag)), ncol = data$n_Ecov), map.Ecov.obs.logsigma)
      data$Ecov_use_obs <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(ecov$lag)), ncol = data$n_Ecov), data$Ecov_use_obs)
      data$Ecov_year <- c(seq(data$year1_model - max(ecov$lag), data$year1_Ecov-1), data$Ecov_year)
      data$year1_Ecov <- data$year1_model - max(ecov$lag)
    }

    # pad Ecov if it ends before last model year
    if(end_Ecov < end_model){
      print("Ecov last year is before model last year. Padding Ecov...")
      data$Ecov_obs <- rbind(data$Ecov_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      par.Ecov.obs.logsigma <- rbind(par.Ecov.obs.logsigma, matrix(-1.3, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      map.Ecov.obs.logsigma <- rbind(map.Ecov.obs.logsigma, matrix(NA, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_use_obs <- rbind(data$Ecov_use_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_year <- c(data$Ecov_year, seq(end_Ecov+1, end_model))
      end_Ecov <- end_model
    }
    data$n_years_Ecov <- dim(data$Ecov_obs)[1] # num years Ecov to model (padded)
    data$Ecov_use_re <- matrix(1, nrow=data$n_years_Ecov, ncol=data$n_Ecov)

    # get index of Ecov_x to use for Ecov_out (Ecovs can have diff lag)
    data$ind_Ecov_out_start <- data$ind_Ecov_out_end <- rep(NA, data$n_Ecov)
    for(i in 1:data$n_Ecov){
      data$ind_Ecov_out_start[i] <- which(data$Ecov_year==data$year1_model)-ecov$lag[i]-1 # -1 is for cpp indexing
      data$ind_Ecov_out_end[i] <- which(data$Ecov_year==end_model)-ecov$lag[i]-1 # -1 is for cpp indexing
    }

    if(!identical(length(ecov$lag), length(ecov$label), data$n_Ecov)) stop("Length of Ecov_lag and Ecov_label vectors not equal to # Ecov")
    data$Ecov_lag <- ecov$lag
    if(!all(ecov$process_model %in% c(NA,"rw", "ar1"))){
      stop("ecov$process_model must be 'rw' (random walk), 'ar1', or NA (do not fit)")
    }
    if(is.na(ecov$process_model) && ecov$how !=0){
      stop("ecov$process_model not chosen (NA) but ecov$how specified.
       Either 1) choose an ecov process model ('rw' or 'ar1'),
              2) turn off ecov (set ecov$how = 0 and ecov$process_model = NA),
           or 3) fit ecov but with no effect on population (ecov$how = 0, ecov$process_model 'rw' or 'ar1').")
    }
    data$Ecov_model <- sapply(ecov$process_model, match, c("rw", "ar1"))

    if(any(ecov$where == 'recruit') & data$n_NAA_sigma == 0){
      stop("Cannot estimate ecov effect on recruitment when
      recruitment in each year is estimated freely as fixed effect parameters.
      Either remove ecov-recruit effect or estimate recruitment
      (or all numbers-at-age) as random effects.")
    }
    if(is.null(ecov$where)) stop("ecov$where must be specified, 'recruit' or 'M'")
    if(!any(ecov$where %in% c('recruit','M','none'))){
      stop("Only ecov effects on recruitment and M currently implemented.
      Set ecov$where = 'recruit', 'M', or 'none'.")
    }
    data$Ecov_where <- sapply(ecov$where, match, c('none','recruit','M')) - 1

    if(is.null(ecov$how)) stop("ecov$how must be specified")
    if(length(ecov$how) != data$n_Ecov) stop("ecov$how must be a vector of length(n.ecov)")
    for(i in 1:data$n_Ecov){
      if(data$Ecov_where[i] == 2) if(!ecov$how[i] %in% c(0,1)){
        stop("Sorry, only the following ecov effects on M are currently implemented.
        Set ecov$how = 0 (no effect) or 1 (effect on mean M, shared across ages).")
      }
      if(data$Ecov_where[i] == 1) if(!ecov$how[i] %in% c(0,1,2,4)){
        stop("Sorry, only the following ecov effects on recruitment are currently implemented.
        Set ecov$how = 0 (no effect), 1 (controlling), 2 (limiting, Bev-Holt only), or 4 (masking).")
      }
      if(data$Ecov_where[i] == 1 & data$recruit_model == 1){
        stop("Random walk recruitment cannot have an ecov effect on recruitment.
        Either choose a different recruit_model (2, 3, or 4), or remove the Ecov effect.")
      }
      if(data$Ecov_where[i] == 1) if(data$recruit_model == 4 & ecov$how[i] == 2){
        stop("'Limiting' ecov effect on Ricker recruitment not implemented.
        Either set ecov$how = 0 (no effect), 1 (controlling), or 4 (masking)...
        Or set recruit_model = 3 (Bev-Holt).")
      }      
    }
    data$Ecov_how <- ecov$how
    data$Ecov_poly <- rep(1,data$n_Ecov)
    ecov_str <- as.list(rep('linear',data$n_Ecov))
    if(!is.null(ecov$link_model)){
      if(!is.na(ecov$link_model)){
        for(i in 1:data$n_Ecov){
          ecov_str[[i]] <- strsplit(ecov$link_model[i],"-")[[1]]
          if(!ecov_str[[i]][1] %in% c('linear','poly')) stop("Only 'linear' or 'poly-x' (x = 1, 2, ...) ecov link models implemented.")
          if(ecov_str[[i]][1]=='linear') data$Ecov_poly[i] <- 1
          if(ecov_str[[i]][1]=='poly') data$Ecov_poly[i] <- as.numeric(ecov_str[[i]][2])
          ecov_str[[i]] = ecov$link_model[i]
        }
      }
    }
    if(!is.null(ecov$ages)){
      if(length(ecov$ages) != data$n_Ecov) stop("ecov$ages must be a list of length n.ecov")
      for(i in 1:data$n_Ecov){
        if(!all(ecov$ages[[i]] %in% 1:data$n_ages)) stop("All ecov$ages must be in 1:n.ages")
      }
    } else {
      ecov$ages <- vector("list", data$n_Ecov)
      for(i in 1:data$n_Ecov) ecov$ages[[i]] <- 1:data$n_ages # default: ecov affects all ages
    }

    cat(paste0("Please check that the environmental covariates have been loaded
and interpreted correctly.

Model years: ", data$year1_model, " to ", end_model,"
Ecov years: ", data$year1_Ecov, " to ", end_Ecov,"

"))
    for(i in 1:data$n_Ecov){
      years <- data$Ecov_year[as.logical(data$Ecov_use_obs[,i])]

      if(data$Ecov_where[i] == 1){ # recruitment
        cat(paste0("Ecov ",i,": ",ecov$label[i],"
",c('*NO*','Controlling','Limiting','Lethal','Masking','Directive')[data$Ecov_how[i]+1]," (",ecov_str[[i]],") effect on: ", c('recruitment','M')[data$Ecov_where[i]],"

In model years:
"))
      }
      if(data$Ecov_where[i] == 2){ # M
        cat(paste0("Ecov ",i,": ",ecov$label[i],"
",c('*NO*',ecov_str[[i]])[data$Ecov_how[i]+1]," effect on: ", c('recruitment','M')[data$Ecov_where[i]],"

In model years:
"))
      }

cat(years, fill=TRUE)
lastyr <- tail(years,1)
cat(paste0("Lag: ",data$Ecov_lag[i],"
Ex: ",ecov$label[i]," in ",years[1]," affects ", c('recruitment','M')[data$Ecov_where[i]]," in ",years[1+data$Ecov_lag[i]],"
    ",ecov$label[i]," in ",lastyr," affects ", c('recruitment','M')[data$Ecov_where[i]]," in ",lastyr+data$Ecov_lag[i],"

"))
    }
    data$Ecov_label <- list(data$Ecov_label)
  } # end load Ecov

  # Ecov pars
  par$Ecov_re = matrix(rnorm(data$n_years_Ecov*data$n_Ecov), data$n_years_Ecov, data$n_Ecov)
  max.poly <- max(data$Ecov_poly)
  par$Ecov_beta = array(0, dim=c(max.poly, data$n_Ecov, data$n_ages)) # beta_R in eqns 4-5, Miller et al. (2016)
  par$Ecov_process_pars = matrix(0, 3, data$n_Ecov) # nrows = RW: 2 par (Ecov1, log_sig), AR1: 3 par (mu, log_sig, phi); ncol = N_ecov
  par$Ecov_process_pars[1,] = -1.3 # start sig_ecov at 0.27
  par$Ecov_obs_logsigma <- par.Ecov.obs.logsigma
  par$Ecov_obs_sigma_par <- par.Ecov.obs.sigma.par

  # turn off 3rd Ecov par if it's a RW
  tmp.pars <- par$Ecov_process_pars
  for(i in 1:data$n_Ecov) tmp.pars[3,i] <- ifelse(data$Ecov_model[i]==1, NA, 0)
  ind.notNA <- which(!is.na(tmp.pars))
  tmp.pars[ind.notNA] <- 1:length(ind.notNA)

  # turn off Ecov_beta to fit Ecov process model without effect on population
  tmp <- array(NA, dim=dim(par$Ecov_beta))
  ct = 1
  for(j in 1:data$n_Ecov){
    for(i in 1:max.poly){
      for(a in 1:data$n_ages){
        if(data$Ecov_how[j] > 0 & i <= data$Ecov_poly[j] & a %in% ecov$ages[[j]]) tmp[i,j,a] = ct # default share ecov_beta across ages
      }
      ct = ct+1
    }
  }
  map$Ecov_beta = factor(tmp)

  # turn off Ecov pars if no Ecov (re, process)
  # for any Ecov_model = NA, ecov$how must be 0 and beta is already turned off
  data$Ecov_model[is.na(data$Ecov_model)] = 0 # turn any NA into 0
  tmp.re <- matrix(1:length(par$Ecov_re), data$n_years_Ecov, data$n_Ecov, byrow=FALSE)
  for(i in 1:data$n_Ecov){
    tmp.pars[,i] <- if(data$Ecov_model[i]==0) rep(NA,3) else tmp.pars[,i]
    tmp.re[,i] <- if(data$Ecov_model[i]==0) rep(NA,data$n_years_Ecov) else tmp.re[,i]
    if(data$Ecov_model[i]==1) tmp.re[1,i] <- NA # if Ecov is a rw, first year of Ecov_re is not used bc Ecov_x[1] uses Ecov1 (fixed effect)
  }
  ind.notNA <- which(!is.na(tmp.re))
  tmp.re[ind.notNA] <- 1:length(ind.notNA)
  map$Ecov_re = factor(tmp.re)
  ind.notNA <- which(!is.na(tmp.pars))
  tmp.pars[ind.notNA] <- 1:length(ind.notNA)
  map$Ecov_process_pars = factor(tmp.pars)
  map$Ecov_obs_logsigma <- factor(map.Ecov.obs.logsigma)
  map$Ecov_obs_sigma_par <- factor(map.Ecov.obs.sigma.par)

  input$data = data
  input$par = par
  input$map = map
  #input$random = random
  return(input)
}
