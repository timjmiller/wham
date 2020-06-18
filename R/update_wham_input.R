#' Update input data and parameters from simulated data with projection
#'
#' When data are simulated from a wham model with projections, this function can be used to update the input model years.
#' This is helpful for using WHAM as an assessment model in a management strategy evaluation. See \code{prepare_wham_input}, 
#' \code{prepare_wham_om_input}, \code{prepare_wham_om_proj}, and \code{project_wham} for further details on data and parameter structure.
#'
#' @param simres list containing data and parameters (output from model$simulate(complete=TRUE), where model is output from \code{fit_wham})
#' @param model object used to produce simres.
#' @param n_years_add number of projection years to convert to model_years
#'
#' @return a named list with the following components:
#'   \describe{
#'     \item{\code{data}}{Named list of data, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{par}}{Named list of parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{map}}{Named list defining how to optionally collect and fix parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{random}}{Character vector of parameters to treat as random effects, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{years}}{Numeric vector of years to fit WHAM model}
#'     \item{\code{years_full}}{Numeric vector of years to fit and project WHAM model}
#'     \item{\code{ages.lab}}{Character vector of age labels, ending with plus-group}
#'     \item{\code{model_name}}{Character, name of stock/model (specified in call to \code{prepare_wham_input})}
#'   }
#'
#' @seealso \code{\link{read_asap3_dat}}, \code{\link{fit_wham}}, \href{https://www.nefsc.noaa.gov/nft/ASAP.html}{ASAP}, \href{https://www.sciencedirect.com/science/article/pii/S1385110197000221}{Iles & Beverton (1998)}
#'
#' @examples
#' \dontrun{
#' input = prepare_wham_om_input(asap3)
#' mod = fit_wham(input, do.fit = FALSE, do.proj = TRUE)
#' sim = mod$simulate(complete=TRUE)
#' sim_input = update_wham_input(sim, mod, 3)
#' simfit = fit_wham(sim_input, MakeADFun.silent = TRUE)
#' plot_wham_output(simfit)
#' }
#'
#' @export
update_wham_input <- function(simres, model, n_years_add){
  if(n_years_add> simres$n_years_proj) stop("the number of years to update is greater than the number of projection years in the simulated data.")
  map = model$input$map
  if(!is.null(model$parList)) par = model$parList 
  else par = model$input$par
  random = model$input$random
  
  nyo = simres$n_years_model
  nya = n_years_add
  nym = ny0 + nya
  data = list(n_years_model = nym)
  data$n_years_catch = nym
  data$n_years_indices = nym
  data$n_ages = simres$n_ages
  data$n_fleets = simres$n_fleets
  data$n_indices <- simres$n_indices

  # Selectivity
  data$n_selblocks <- simres$n_selblocks
  data$selblock_models <- simres$selblock_models
  data$selblock_models_re <- simres$selblock_models_re
  data$n_selpars <- simres$n_selpars
  data$selblock_pointer_fleets = rbind(cbind(simres$selblock_pointer_fleets), t(matrix(simres$selblock_pointer_fleets[nyo,],data$n_fleets,nya)))
  data$selblock_pointer_indices = rbind(cbind(simres$selblock_pointer_indices), t(matrix(simres$selblock_pointer_indices[nyo,],data$n_indices,nya)))
  selblock_pointers <- cbind(data$selblock_pointer_fleets, data$selblock_pointer_indices)
  data$selblock_years <- matrix(0, nrow=data$n_years_model, ncol=data$n_selblocks)
  for(b in 1:data$n_selblocks) data$selblock_years[,b] <- apply(selblock_pointers, 1, function(x) b %in% x)
  data$n_years_selblocks <- apply(data$selblock_years, 2, sum)

  # Age composition model
  data$age_comp_model_fleets = simres$age_comp_model_fleets
  data$n_age_comp_pars_fleets = simres$n_age_comp_pars_fleets
  data$age_comp_model_indices = simres$age_comp_model_indices
  data$n_age_comp_pars_indices = simres$n_age_comp_pars_indices

  data$fracyr_SSB = simres$fracyr_SSB[1:nym]
  data$mature = simres$mature[1:nym,]

  # Weight-at-age
  data$waa_pointer_fleets = inputt$waa_pointer_fleets
  data$waa_pointer_totcatch = inputt$waa_pointer_totcatch
  data$waa_pointer_indices = inputt$waa_pointer_indices
  data$waa_pointer_ssb = inputt$waa_pointer_ssb
  data$waa_pointer_jan1 = inputt$waa_pointer_jan1
  data$waa_pointer_totcatch = asap3$WAA_pointers[data$n_fleets + 1]
  data$waa_pointer_indices = asap3$index_WAA_pointers
  data$waa_pointer_ssb = asap3$WAA_pointers[data$n_fleets + 2]
  data$waa_pointer_jan1 = asap3$WAA_pointers[data$n_fleets + 3]
  data$waa = simres$waa[,1:nym,,drop=FALSE]
  
  # Catch
  data$agg_catch = rbind(cbind(simres$agg_catch),cbind(simres$agg_catch_proj[1:nya,]))
  data$agg_catch_sigma = rbind(cbind(simres$catch_cv), t(matrix(simres$catch_cv[nyo,], data$n_fleets, nya)))
  data$catch_paa = array(NA, dim = c(data$n_fleets, data$n_years_model, data$n_ages))
  data$catch_paa[,1:nyo,] = simres$catch_paa
  data$catch_paa[,nyo+1:nya,] = simres$catch_paa_proj
  data$use_agg_catch = rbind(cbind(simres$use_agg_catch), t(matrix(simres$use_agg_catch[nyo,],data$n_fleets, nya)))
  data$use_catch_paa = rbind(cbind(simres$use_catch_paa), t(matrix(simres$use_catch_paa[nyo,],data$n_fleets, nya)))
  data$catch_Neff = rbind(cbind(simres$catch_Neff), t(matrix(simres$catch_Neff[nyo,],data$n_fleets, nya)))
  data$catch_aref = matrix(NA, data$n_years_model, data$n_fleets)
  for(i in 1:data$n_fleets) data$catch_aref[,i] = get_aref_fn(data$catch_paa[i,,])

  # Indices/surveys
  data$units_indices <- simres$units_indices
  data$fracyr_indices = rbind(cbind(simres$fracyr_indices), t(matrix(simres$fracyr_indices[nyo,],data$n_indices, nya)))
  data$agg_indices = rbind(cbind(simres$agg_indices),cbind(simres$agg_indices_proj[1:nya,]))
  data$use_indices = rbind(cbind(simres$use_indices), t(matrix(simres$use_indices[nyo,],data$n_indices, nya)))
  data$agg_index_sigma = rbind(cbind(simres$agg_index_sigma), t(matrix(simres$agg_index_sigma[nyo,],data$n_indices, nya)))
  data$units_index_paa <- simres$units_index_paa
  data$index_paa = array(NA, dim = c(data$n_indices, data$n_years_model, data$n_ages))
  data$index_paa[,1:nyo,] = simres$index_paa
  data$index_paa[,nyo+1:nya,] = simres$index_paa_proj
  data$use_index_paa = rbind(cbind(simres$use_index_paa), t(matrix(simres$use_index_paa[nyo,],data$n_indices, nya)))
  data$index_Neff = rbind(cbind(simres$index_Neff), t(matrix(simres$index_Neff[nyo,],data$n_indices, nya)))
  data$index_aref = matrix(NA, data$n_years_model, data$n_indices)
  for(i in 1:data$n_indices) data$index_aref[,i] = get_aref_fn(data$index_paa[i,,])
  data$q_lower <- simres$q_lower
  data$q_upper <- simres$q_upper

  data$selpars_est = simres$selpars_est
  data$n_selpars_est <- apply(data$selpars_est > 0, 1, sum)
  data$selpars_lower = simres$selpars_lower
  data$selpars_upper = simres$selpars_upper
    
  data$n_NAA_sigma = simres$n_NAA_sigma
  data$NAA_sigma_pointers = simres$NAA_sigma_pointers
  data$recruit_model = simres$recruit_model
  
  data$n_M_a = simres$n_M_a
  data$M_model = simres$M_model
  data$M_re_model = simres$M_re_model
  data$M_est = simres$M_est
  data$n_M_est <- sum(data$M_est)

  data$use_b_prior = simres$use_b_prior
  
  data$N1_model = simres$N1_model
  data$which_F_age = simres$which_F_age
  data$use_steepness = simres$use_steepness
  data$bias_correct_pe = simres$bias_correct_pe
  data$bias_correct_oe = simres$bias_correct_oe
  data$Fbar_ages = simres$Fbar_ages
  data$simulate_state = simres$simulate_state
  data$simulate_data = simres$simulate_data
  data$simulate_period = simres$simulate_period
  data$percentSPR = simres$percentSPR
  data$XSPR_R_opt = simres$XSPR_R_opt
  data$XSPR_R_avg_yrs = simres$XSPR_R_avg_yrs
  

  model_years <- model$years_full[1:nym]

  data$n_Ecov <- simres$n_Ecov
  data$Ecov_obs_sigma_opt <- simres$Ecov_obs_sigma_opt
  data$year1_Ecov <- simres$year1_Ecov
  data$Ecov_lag <- simres$Ecov_lag
  data$year1_model <- simres$year1_model
  data$Ecov_model <- simres$Ecov_model
  data$Ecov_where <- simres$Ecov_where
  data$Ecov_how <- simres$Ecov_how
  data$Ecov_poly <- simres$Ecov_poly
  data$Ecov_recruit <- simres$Ecov_recruit
  data$Ecov_growth <- simres$Ecov_growth
  data$Ecov_mortality <- simres$Ecov_mortality
  data$Ecov_label <- simres$Ecov_label
  data$Ecov_obs <- simres$Ecov_obs
  data$Ecov_use_obs <- simres$Ecov_use_obs
  data$Ecov_year <- simres$Ecov_year
  data$n_years_Ecov <- nyeo <-simres$n_years_Ecov
  #nyeo is the number of years of Ecov random effects in unupdated model
  data$ind_Ecov_out_start <- simres$ind_Ecov_out_start
  data$ind_Ecov_out_end <- simres$ind_Ecov_out_end
  data$Ecov_use_re <- simres$Ecov_use_re

  # --------------------------------------------------------------------------------
  # Environmental covariate data
  if(NROW(data$Ecov_obs)>1){
    # Handle Ecov sigma options
    n_Ecov_obs <- dim(data$Ecov_obs)[1] # num Ecov obs

    end_model <- tail(model_years,1) #last updated model year
    end_Ecov <- model_years[nyo] #end of model before update

    # pad Ecov if it ends before last model year
    if(end_Ecov < end_model){
      print("Ecov last year is before model last year. Padding Ecov...")
      # warning("Ecov last year is before model last year. Padding Ecov...")
      #IFF simulate() for wham_v0.cpp is changed to include simulating Ecov_obs in the projection years, then change the next few lines
      data$Ecov_obs <- rbind(data$Ecov_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      #par.Ecov.obs.logsigma <- rbind(par.Ecov.obs.logsigma, matrix(-1.3, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      #map.Ecov.obs.logsigma <- rbind(map.Ecov.obs.logsigma, matrix(NA, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_use_obs <- rbind(data$Ecov_use_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_year <- c(data$Ecov_year, seq(end_Ecov+1, end_model))
      end_Ecov <- end_model
    }
    data$n_years_Ecov <- NROW(data$Ecov_obs) # num years Ecov to model (padded)
    data$Ecov_use_re <- matrix(1, nrow=data$n_years_Ecov, ncol=data$n_Ecov)

    # get index of Ecov_x to use for Ecov_out (Ecovs can have diff lag)
    data$ind_Ecov_out_start <- data$ind_Ecov_out_end <- rep(NA, data$n_Ecov)
    for(i in 1:data$n_Ecov){
      data$ind_Ecov_out_start[i] <- which(data$Ecov_year==data$year1_model)-data$Ecov_lag[i]-1 # -1 is for cpp indexing
      data$ind_Ecov_out_end[i] <- which(data$Ecov_year==end_model)-data$Ecov_lag[i]-1 # -1 is for cpp indexing
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
",c('*NO*','Controlling','Limiting','Lethal','Masking','Directive')[data$Ecov_how+1]," (",ecov_str[[i]],") effect on: ", c('recruitment','M')[data$Ecov_where[i]],"

In model years:
"))
      }
      if(data$Ecov_where[i] == 2){ # M
        cat(paste0("Ecov ",i,": ",ecov$label[i],"
",ecov_str[[i]]," effect on: ", c('recruitment','M')[data$Ecov_where[i]],"

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
  } # end load Ecov

  # add vector of all observations for one step ahead residuals ==========================
  # 5 components: fleet catch (log), index catch (log), Ecov, paa catch, paa index
  obs.colnames <- c("year","fleet","age","type","val")
  obs <- data.frame(matrix(ncol = length(obs.colnames), nrow = 0))
  colnames(obs) <- obs.colnames

  # 1. log fleet catch
  x <- as.data.frame(data$agg_catch)
  x[data$use_agg_catch==0] <- NA # can't fit to fleets/years with 0 catch
  colnames(x) <- paste0("fleet_", 1:data$n_fleets)
  x$year <- 1:data$n_years_catch
  tmp <- tidyr::gather(x, fleet, val, -year)
  tmp <- tmp[complete.cases(tmp),]  
  tmp$val <- log(tmp$val) # all obs of 0 catch should have use_agg_catch==0, turned to NA, and removed
  tmp$age <- NA
  tmp$type <- "logcatch"
  obs <- rbind(obs, tmp[, obs.colnames])

  # 2. log index catch
  x <- as.data.frame(data$agg_indices)
  x[data$use_indices==0] <- NA # only include index data to fit in obsvec
  colnames(x) <- paste0("index_", 1:data$n_indices)
  x$year <- 1:data$n_years_indices # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s
  tmp <- tidyr::gather(x, fleet, val, -year)
  tmp <- tmp[complete.cases(tmp),]
  tmp$val <- log(tmp$val) # all obs of 0 catch should have use_indices==0, turned to NA, and already removed
  tmp$age <- NA
  tmp$type <- "logindex"
  obs <- rbind(obs, tmp[, obs.colnames])

  # 3. Ecov
  if(!all(data$Ecov_use_obs==0)){
    x <- as.data.frame(data$Ecov_obs)
    x[data$Ecov_use_obs==0] <- NA # only include index data to fit in obsvec
    colnames(x) <- paste0("Ecov_", 1:data$n_Ecov)
    # x$year <- 1:data$n_years_Ecov # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s
    x$year <- seq(from=data$year1_Ecov-data$year1_model+1, length.out=data$n_years_Ecov) # don't assume Ecov and model years are the same
    tmp <- tidyr::gather(x, fleet, val, -year)
    tmp <- tmp[complete.cases(tmp),]
    tmp$age <- NA
    tmp$type <- "Ecov"
    obs <- rbind(obs, tmp[, obs.colnames])
}

  # # 4. paa catch
  # dimnames(data$catch_paa) <- list(fleet=paste0("fleet_", 1:data$n_fleets),
  #                                  year=1:data$n_years_catch,
  #                                  age=1:data$n_ages)
  # x <- as.data.frame(dplyr::as.tbl_cube(data$catch_paa, met_name = "val"))
  # x$type <- "paacatch"
  # obs <- rbind(obs, x[, obs.colnames])

  # # 5. paa index
  # dimnames(data$index_paa) <- list(fleet=paste0("index_", 1:data$n_indices),
  #                                  year=1:data$n_years_indices,
  #                                  age=1:data$n_ages)
  # x <- as.data.frame(dplyr::as.tbl_cube(data$index_paa, met_name = "val"))
  # x$type <- "paaindex"
  # obs <- rbind(obs, x[, obs.colnames])

  # order by year, fleet, age, type
  o <- order(as.numeric(obs$year), obs$fleet, as.numeric(obs$age), obs$type)
  obs <- obs[o,]

  # calculate obsvec indices in keep arrays
  obs$ind <- 1:dim(obs)[1]
  # data$keep_C <- matrix(subset(obs, type=='logcatch')$ind, nrow=data$n_years_catch, ncol=data$n_fleets, byrow=TRUE)
  data$keep_C <- matrix(NA, nrow=data$n_years_catch, ncol=data$n_fleets)
  xl <- lapply(seq_len(nrow(data$use_agg_catch)), function(r) which(data$use_agg_catch[r,]==1))
  Col <- unlist(xl)
  Row <- rep(1:data$n_years_catch, times=sapply(xl, length))
  data$keep_C[cbind(Row,Col)] <- subset(obs, type=='logcatch')$ind

  data$keep_I <- matrix(NA, nrow=data$n_years_indices, ncol=data$n_indices)
  # data$keep_I[data$use_indices==1] <- subset(obs, type=='logindex')$ind
  # xl <- apply(data$use_indices,1,function(r) which(r==1))
  xl <- lapply(seq_len(nrow(data$use_indices)), function(r) which(data$use_indices[r,]==1))
  Col <- unlist(xl)
  Row <- rep(1:data$n_years_indices, times=sapply(xl, length))
  data$keep_I[cbind(Row,Col)] <- subset(obs, type=='logindex')$ind

  data$keep_E <- matrix(NA, nrow=data$n_years_Ecov, ncol=data$n_Ecov)
  # data$keep_E[data$Ecov_use_obs==1] <- subset(obs, type=='Ecov')$ind
  xl <- lapply(seq_len(nrow(data$Ecov_use_obs)), function(r) which(data$Ecov_use_obs[r,]==1))
  # xl <- apply(data$Ecov_use_obs,1,function(r) which(r==1))
  Col <- unlist(xl)
  Row <- rep(1:data$n_years_Ecov, times=sapply(xl, length))
  data$keep_E[cbind(Row,Col)] <- subset(obs, type=='Ecov')$ind

  data$keep_Cpaa <- array(NA, dim=c(data$n_fleets, data$n_years_catch, data$n_ages))
  for(i in 1:data$n_fleets) data$keep_Cpaa[i,,] <- matrix(subset(obs, type=='paacatch' & fleet==paste0("fleet_",i))$ind, nrow=data$n_years_catch, ncol=data$n_ages, byrow=TRUE)
  data$keep_Ipaa <- array(NA, dim=c(data$n_indices, data$n_years_indices, data$n_ages))
  for(i in 1:data$n_indices) data$keep_Ipaa[i,,] <- matrix(subset(obs, type=='paaindex' & fleet==paste0("index_",i))$ind, nrow=data$n_years_indices, ncol=data$n_ages, byrow=TRUE)
  # subtract 1 bc TMB indexes from 0
  data$keep_C <- data$keep_C - 1
  data$keep_I <- data$keep_I - 1
  data$keep_E <- data$keep_E - 1
  data$keep_Cpaa <- data$keep_Cpaa - 1
  data$keep_Ipaa <- data$keep_Ipaa - 1

  data$obs <- obs
  data$obsvec <- obs$val

  # projection data will always be modified by 'prepare_projection'
  data$do_proj <- 0
  data$n_years_proj <- 0
  data$n_years_proj_Ecov <- 0
  data$avg_years_ind <- 0
  data$proj_F_opt <- 0
  data$proj_Fcatch <- 0
  data$proj_M_opt <- 0

  # data$obsvec[data$keep_I[data$use_indices==1]+1] - log(data$agg_indices[data$use_indices==1])
  # data$obsvec[data$keep_E[data$Ecov_use_obs==1]+1] - data$Ecov_obs[data$Ecov_use_obs==1]

  # -------------------------------------------------------------------
  # Parameters
  #par$mean_rec_pars = numeric(c(0,1,2,2)[recruit_model])
  #if(recruit_model==2) par$mean_rec_pars = 10
  #if(recruit_model==4) par$mean_rec_pars[2] = -10
  #par$logit_q = rep(-8, data$n_indices)
  #par$log_F1 = rep(-2, data$n_fleets)
  F = matrix(NA, nym,data$n_fleets)
  par$F_devs = matrix(0, nym-1, data$n_fleets)
  for(i in 1:data$n_fleets) par$F_devs[,i] = diff(log(simres$FAA[i,1:nym,data$which_F_age]))
   
  #par$F_devs = matrix(0, data$n_years_model-1, data$n_fleets)
  #if(data$N1_model == 1) par$log_N1_pars = c(10,log(0.1))
  #if(data$N1_model == 0) par$log_N1_pars = rep(10,data$n_ages)
  
  # NAA_re pars
  par$log_NAA = cbind(simres$log_NAA[1:(nym-1),])
  

  tmp <- par$log_NAA
  if(data$n_NAA_sigma < 2) tmp[,-1] <- NA # always estimate Rec devs (col 1), whether random effect or not
  ind.notNA <- which(!is.na(tmp))
  tmp[ind.notNA] <- 1:length(ind.notNA)
  map$log_NAA = factor(tmp)

  # selectivity pars
  #par$logit_selpars = log(selpars_ini-selpars_lo) - log(selpars_hi - selpars_ini)
  #par$logit_selpars[!is.na(map$logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars<0] = -10
  #par$logit_selpars[!is.na(map$logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars>0] = 10
  # number of estimated selpars per block * number of years per block (only if that block has re)
  if(any(data$selblock_models_re > 1)){
    tmp_vec = c()
    ind = 0
    for(b in 1:n_selblocks) 
    {
      tmp = matrix(NA, n_years_selblocks[b], n_selpars_est[b])
      for(j in 1: n_selpars_est(b)) {
        tmp[,j] = c(tmp, par$selpars_re[ind + 1:(n_years_selblocks(b)-nya)], rep(0, nya))
        ind = ind + n_years_selblocks(b) - nya
      }
      tmp_vec = c(tmp_vec, tmp)
    }
    par$selpars_re = tmp_vec
    tmp_vec <- c()
    ct <- 0
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] > 1){
        tmp <- matrix(0, nrow=data$n_years_selblocks[b], ncol=data$n_selpars_est[b])
        if(data$selblock_models_re[b] %in% c(2,5)){ # 2d ar1
          tmp[] = 1:(dim(tmp)[1]*dim(tmp)[2]) + ct # all y,a estimated
        }
        if(data$selblock_models_re[b] == 3){ # ar1_a (devs by age, constant by year)
          for(i in 1:dim(tmp)[2]) tmp[,i] = (i + ct)
        }
        if(data$selblock_models_re[b] == 4){ # ar1_y (devs by year, constant by age)
          for(i in 1:dim(tmp)[1]) tmp[i,] = (i + ct)
        }
        ct = max(tmp)
        tmp_vec = c(tmp_vec, as.vector(tmp))
      }
    }
    map$selpars_re <- factor(tmp_vec)
  } else {
    #par$selpars_re <- matrix(0)
    #map$selpars_re <- factor(NA)
  }
  #par$sel_repars <- matrix(0, nrow=data$n_selblocks, ncol=3)
  #par$sel_repars[,1] <- log(0.1) # start sigma at 0.1, rho at 0
  for(b in 1:data$n_selblocks){
    #if(data$selblock_models_re[b] == 3) par$sel_repars[b,3] <- 0 # if ar1 over ages only, fix rho_y = 0
    #if(data$selblock_models_re[b] == 4) par$sel_repars[b,2] <- 0 # if ar1 over years only, fix rho = 0
    # check if only 1 estimated sel par (e.g. because all but 1 age is fixed), can't estimate rho
    #if(data$n_selpars_est[b] < 2) par$sel_repars[b,2] <- 0
  }

  # age comp pars
  n_catch_acomp_pars = c(0,1,1,3,1,2)[data$age_comp_model_fleets[which(apply(data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2)[data$age_comp_model_indices[which(apply(data$use_index_paa,2,sum)>0)]]
  #par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  #par$index_paa_pars = rep(0, sum(n_index_acomp_pars))

  # natural mortality pars
  #par$M0 <- M0_ini # mean M
  #par$M_a <- M_a_ini # deviations by age
  par$M_re <- rbind(par$M_re, matrix(0,nya,NCOL(par$M_re))) # deviations from mean M_a on log-scale, PARAMETER_ARRAY
  #par$M_repars <- rep(0, 3)
  #par$M_repars[1] <- log(0.1) # start sigma at 0.1, rho at 0
  #if(data$M_re_model == 3) par$M_repars[3] <- 0 # if ar1 over ages only, fix rho_y = 0
  #if(data$M_re_model == 4) par$M_repars[2] <- 0 # if ar1 over years only, fix rho_a = 0
  # check if only 1 estimated mean M (e.g. because weight-at-age M or if all but 1 age is fixed), can't estimate rho_a
  # if(data$n_M_est < 2) par$M_repars[2] <- 0
  #par$log_b = log(0.305)
  #par$log_catch_sig_scale = rep(0, data$n_fleets)
  #par$log_index_sig_scale = rep(0, data$n_indices)

  # Ecov pars
  #par$Ecov_re = matrix(0, data$n_years_Ecov, data$n_Ecov)
  if(data$n_years_Ecov>nyeo) par$Ecov_re = rbind(par$Ecov_re, matrix(rnorm((data$n_years_Ecov-nyeo)*data$n_Ecov), data$n_years_Ecov-nyeo, data$n_Ecov))
  #par$Ecov_re = matrix(rnorm(data$n_years_Ecov*data$n_Ecov), data$n_years_Ecov, data$n_Ecov)
  max.poly <- max(data$Ecov_poly)
  #par$Ecov_beta = matrix(0, nrow=max.poly, ncol=data$n_Ecov) # beta_R in eqns 4-5, Miller et al. (2016)
  #par$Ecov_process_pars = matrix(0, 3, data$n_Ecov) # nrows = RW: 2 par (log_sig, Ecov1), AR1: 3 par (mu, phi, log_sig); ncol = N_ecov
  #par.Ecov.obs.logsigma <- rbind(par.Ecov.obs.logsigma, matrix(-1.3, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
  #map.Ecov.obs.logsigma <- rbind(map.Ecov.obs.logsigma, matrix(NA, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
  if(data$n_years_Ecov>nyeo) {
    if(data$Ecov_obs_sigma_opt %in% c(1,3:4)) {
      par$Ecov_obs_logsigma <- rbind(cbind(par$Ecov_obs_logsigma), matrix(-1.3, data$n_years_Ecov-nyeo, data$n_Ecov))
      if(data$Ecov_obs_sigma_opt == 1) map$Ecov_obs_logsigma = factor(rep(NA, data$n_years_Ecov*data$n_Ecov))
      if(data$Ecov_obs_sigma_opt %in% 3:4) map$Ecov_obs_logsigma = factor(rbind(matrix(as.integer(map$Ecov_obs_logsigma), nyeo, data$n_Ecov), matrix(NA,data$n_years_Ecov-neo, data$n_Ecov)))
    }
  }

  # turn off Ecov pars if no Ecov (re, process)
  # for any Ecov_model = NA, ecov$how must be 0 and beta is already turned off
  data$Ecov_model[is.na(data$Ecov_model)] = 0 # turn any NA into 0
  tmp.re <- matrix(as.integer(map$Ecov_re), nyeo, data$n_Ecov)
  new.re <- matrix(NA, data$n_years_Ecov-nyeo, data$n_Ecov)
  if(any(!is.na(tmp.re))) {
    if(data$n_years_Ecov>nyeo) new.re[] <- max(tmp.re, na.rm = TRUE) + 1:((data$n_years_Ecov-nyeo)*data$n_Ecov)
    #tmp.re <- matrix(1:length(par$Ecov_re), data$n_years_Ecov, data$n_Ecov, byrow=FALSE)
    for(i in 1:data$n_Ecov){
      new.re[,i] <- if(data$Ecov_model[i]==0) rep(NA,data$n_years_Ecov)
    }
    ind.notNA <- which(!is.na(new.re))
    new.re[ind.notNA] <- max(tmp.re, na.rm = TRUE) + 1:length(ind.notNA)
  }
  map$Ecov_re = factor(rbind(tmp.re,new.re))

  # M_re: "none","iid","ar1_a","ar1_y","2dar1"
  #tmp <- par$M_re
  oldmap = matrix(as.integer(map$M_re), nyo, NCOL(par$M_re))
  tmp <- matrix(NA, nya, NCOL(par$M_re))
  if(any(!is.na(oldmap))) maxmf = max(oldmap, na.rm = TRUE)
  else maxmf = integer()
  #if(data$M_re_model == 1) tmp[] = NA # no RE (either estimate RE for all ages or none at all)
  if(length(maxmf))
  {
    if(data$M_re_model %in% c(2,5)){ # 2d ar1
      tmp[] = maxmf + 1:(dim(tmp)[1]*dim(tmp)[2]) # all y,a estimated
    }
    if(data$M_re_model == 4){ # ar1_y (devs by year, constant by age)
      tmp[] = maxmf + 1:NROW(tmp) #repeats and fills by columns
      #for(i in 1:dim(tmp)[1]) tmp[i,] = maxmf + i
    }
    map$M_re <- factor(rbind(oldmap, tmp))
  }

  return(list(data=data, par = par, map = map, random = random, years = model_years, years_full = model_years,
    ages.lab = model$ages.lab, model_name = model$model_name))
}
