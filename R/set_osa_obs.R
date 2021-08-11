set_osa_obs = function(input)
{
  data = input$data

  # 5 components: fleet catch (log), index catch (log), Ecov, paa catch, paa index
  obs.colnames <- c("year","fleet","age","type","val")
  obs <- data.frame(matrix(ncol = length(obs.colnames), nrow = 0))
  colnames(obs) <- obs.colnames

  # 1. log fleet catch
  x <- as.data.frame(data$agg_catch)
  x[data$use_agg_catch==0] <- NA # can't fit to fleets/years with 0 catch
  colnames(x) <- paste0("fleet_", 1:data$n_fleets)
  obs_levels <- paste0("fleet_", 1:data$n_fleets)
  x$year <- 1:data$n_years_catch
  tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="fleet")
  tmp <- tmp[complete.cases(tmp),]  
  tmp$val <- log(tmp$val) # all obs of 0 catch should have use_agg_catch==0, turned to NA, and removed
  tmp$age <- NA
  tmp$type <- "logcatch"
  obs <- rbind(obs, tmp[, obs.colnames])

  # 2. log index catch
  x <- as.data.frame(data$agg_indices)
  x[data$use_indices==0] <- NA # only include index data to fit in obsvec
  colnames(x) <- paste0("index_", 1:data$n_indices)
  obs_levels <- c(obs_levels, paste0("index_", 1:data$n_indices))
  x$year <- 1:data$n_years_indices # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s
  tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="fleet")
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
    obs_levels <- c(obs_levels, paste0("Ecov_", 1:data$n_Ecov))
    x$year <- seq(from=data$year1_Ecov-data$year1_model+1, length.out=data$n_years_Ecov) # don't assume Ecov and model years are the same
    tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="fleet")
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
  obs$fleet <- factor(obs$fleet, levels=obs_levels)
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

  input$data = data
  return(input)
}