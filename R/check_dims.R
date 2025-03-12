check_dims <- function(input){
	
	input_data_dims <- lapply(input$data, dim)
	vecs <- which(sapply(input_data_dims, is.null))
	input_data_dims[vecs] <- lapply(input$data[vecs], length)

	#################################################
	#basics are single integers
	basics <- c("n_ages", "n_regions", "n_stocks", "n_seasons", "n_fleets", "n_indices", "n_years_model", "n_ages", "n_selblocks", "n_years_Ecov",
		"n_Ecov", "n_years_proj")
	should_dims <- lapply(1:length(input_data_dims), function(x) integer(0))
	names(should_dims) <- names(input_data_dims)
	for(i in basics) should_dims[[i]] <- 1L
	good_dims <- rep(FALSE, length(should_dims))
	names(good_dims) <- names(should_dims)

	good_dims[basics] <- sapply(basics, function(x) identical(should_dims[[x]], input_data_dims[[x]]))
 	if(any(!good_dims[basics])) stop(paste0("These should be a single integer: ", paste0(basics[!good_dims[basics]], collapse = ", ")))
	for(i in basics) assign(i, as.integer(input$data[[i]]))


	#################################################
	#flags are single integers
	do_flags <-	c("do_MSY_BRPs", "do_osa", "do_post_samp_Ecov", "do_post_samp_M", "do_post_samp_mu", "do_post_samp_N", "do_post_samp_q", "do_post_samp_sel", "do_proj", 
	"do_simulate_data", "do_simulate_Ecov_re", "do_simulate_L_re", "do_simulate_M_re", "do_simulate_mu_prior_re",  "do_simulate_mu_re", "do_simulate_N_re",
  "do_simulate_q_prior_re", "do_simulate_q_re", "do_simulate_sel_re", "do_SPR_BRPs")

	for(i in do_flags) should_dims[[i]] <- 1L
 	
 	# if(any(!good_dims[do_flags])) stop(paste0("These should be a single integer: ", paste0(do_flags[!good_dims[do_flags]], collapse = ", ")))
	
	other_1s <- c("M_model", "F_config", "decouple_recruitment", "use_b_prior", "log_b_model", 	"bias_correct_brps", "bias_correct_oe", "bias_correct_pe", "n_regions_is_small",
		"SPR_weight_type", "XSPR_R_opt", "which_F_age_static", "use_alt_AR1", "proj_M_opt", "proj_R_opt", "proj_mu_opt", "proj_L_opt", "percentSPR", "FXSPR_static_init", "FMSY_static_init",
		"percentFMSY", "percentFXSPR")

	for(i in other_1s) should_dims[[i]] <- 1L

	#################################################
	#vectors
	
 	# if(input_data_dims[["do_simulate_period"]] !=should_dims["do_simulate_period"])) stop("input$data$do_simulate_period should be length 2")
 	two <- "do_simulate_period"
 	three <- "do_simulate_data"
	n_stocks_vecs <- c("spawn_regions", "spawn_seasons", "mig_type", "waa_pointer_ssb",  "logR_mean",  "logR_sd", "recruit_model", "waa_pointer_M", "N1_model", "NAA_re_model", 
		"SPR_weights")
 	n_fleets_vecs <- c("waa_pointer_fleets", "fleet_regions", "age_comp_model_fleets")
	n_indices_vecs <- c("index_seasons", "waa_pointer_indices", "index_regions", "age_comp_model_indices", "units_indices", "units_index_paa", "q_lower", "q_upper", "use_q_prior", "logit_q_prior_sigma", "use_q_re")
	n_seasons_vecs <- "fracyr_seasons"
	n_regions_vecs <- "L_model"
	n_Ecov_vecs <- c("Ecov_model", "Ecov_obs_sigma_opt", "Ecov_use_re", "proj_Ecov_opt")
	n_blocks_vecs <- c("selblock_models", "selblock_models_re", "n_selpars", "n_selpars_est", "n_years_selblocks")
	n_y_m_vecs <- c("years_use")
	n_y_p_vecs <- c("F_proj_init", "proj_F_opt")
	n_y_mp_vecs <- c("which_F_age", "FXSPR_init", "FMSY_init") #n_years_model + n_years_proj
	n_years_E_vecs <- c("years_use_Ecov")


	for(i in two) should_dims[[i]] <- 2L
	for(i in three) should_dims[[i]] <- 3L
	for(i in n_stocks_vecs) should_dims[[i]] <- n_stocks
	for(i in n_fleets_vecs) should_dims[[i]] <- n_fleets
	for(i in n_indices_vecs) should_dims[[i]] <- n_indices
	for(i in n_seasons_vecs) should_dims[[i]] <- n_seasons
	for(i in n_regions_vecs) should_dims[[i]] <- n_regions
	for(i in n_Ecov_vecs) should_dims[[i]] <- n_Ecov
	for(i in n_blocks_vecs) should_dims[[i]] <- n_selblocks
	for(i in n_y_m_vecs) should_dims[[i]] <- n_years_model
	for(i in n_y_p_vecs) should_dims[[i]] <- n_years_proj
	for(i in n_y_mp_vecs) should_dims[[i]] <- n_years_model + n_years_proj
	for(i in n_years_E_vecs) should_dims[[i]] <- n_years_Ecov

	n_Fbar_ages_vecs <- "Fbar_ages"
	n_obs_vecs <- c("obsvec", "agesvec")
	n_XSPR_R_avg_yrs_vecs <- "XSPR_R_avg_yrs"
	n_avg_years_ind_static_vecs <- "avg_years_ind_static"
	n_avg_years_ind_vecs <- "avg_years_ind"
	n_avg_years_Ecov_vecs <- "avg_years_Ecov"
	data_not_checked <- c(n_obs_vecs, n_Fbar_ages_vecs, n_XSPR_R_avg_yrs_vecs, n_avg_years_ind_static_vecs, n_avg_years_ind_vecs,n_avg_years_Ecov_vecs,n_y_p_vecs)

	#################################################
 	#matrices

 	n_years_n_stocks <- c("fracyr_SSB")
 	n_fleets_n_seasons <- c("fleet_seasons")
 	n_years_n_fleets <- c("agg_catch", "agg_catch_sigma", "use_agg_catch", "use_catch_paa", "catch_Neff", "keep_C", "selblock_pointer_fleets")
 	n_years_n_indices <- c("agg_indices", "agg_index_sigma", "use_indices", "use_index_paa", "index_Neff", "fracyr_indices", "keep_I", "selblock_pointer_indices")
 	n_years_n_Ecov <- c("Ecov_use_obs", "Ecov_obs", "keep_E")
 	n_years_n_blocks <- c("selblock_years")
 	n_Ecov_n_stocks <- c("Ecov_how_R","ind_Ecov_out_start_R","ind_Ecov_out_end_R", "n_poly_Ecov_R")
 	n_Ecov_n_indices <- c("Ecov_how_q","ind_Ecov_out_start_q","ind_Ecov_out_end_q","n_poly_Ecov_q")
 	n_blocks_n_pars <- c("selpars_est", "selpars_lower", "selpars_upper")
 	n_stocks_n_regions <- c("n_M_re", "M_re_model")
 	n_regions_n_reg_m1 <- c("mu_model")
 	
	for(i in n_years_n_stocks) should_dims[[i]] <- c(n_years_model, n_stocks)
	for(i in n_fleets_n_seasons) should_dims[[i]] <- c(n_fleets, n_seasons)
	for(i in n_years_n_fleets) should_dims[[i]] <- c(n_years_model, n_fleets)
	for(i in n_years_n_indices) should_dims[[i]] <- c(n_years_model, n_indices)
	for(i in n_years_n_Ecov) should_dims[[i]] <- c(n_years_Ecov, n_Ecov)
	for(i in n_years_n_blocks) should_dims[[i]] <- c(n_years_model, n_selblocks)
	for(i in n_Ecov_n_stocks) should_dims[[i]] <- c(n_Ecov, n_stocks)
	for(i in n_Ecov_n_indices) should_dims[[i]] <- c(n_Ecov, n_indices)
	for(i in n_blocks_n_pars) should_dims[[i]] <- c(n_selblocks, n_ages + 6L)
	for(i in n_stocks_n_regions) should_dims[[i]] <- c(n_stocks, n_regions)
	for(i in n_regions_n_reg_m1) should_dims[[i]] <- c(n_regions, n_regions - 1L)
	# for(i in n_proj_n_Ecov) should_dims[[i]] <- c(n_years_proj, n_Ecov)
	# for(i in n_proj_1) should_dims[[i]] <- c(n_years_proj, 1)

 	n_obs_7 <- c("obs")
 	n_proj_n_Ecov <- c("Ecov_use_proj")
 	n_proj_1 <- "proj_Fcatch" #or n_years_proj x n_fleets
 	n_proj_n_Ecov <- c("Ecov_use_proj")
 	n_proj_1 <- "proj_Fcatch" #or n_years_proj x n_fleets
	data_not_checked <- c(data_not_checked, n_obs_7, n_proj_n_Ecov, n_proj_1)


	#################################################

	#arrays
	n_fleets_n_y_2 <- "keep_Cpaa"
	for(i in n_fleets_n_y_2) should_dims[[i]] <- c(n_fleets, n_years_model,2L)
	n_indices_n_y_2 <- "keep_Ipaa"
	for(i in n_indices_n_y_2) should_dims[[i]] <- c(n_indices, n_years_model,2L)
	n_stocks_n_years_n_ages <- c("mature")
	for(i in n_stocks_n_years_n_ages) should_dims[[i]] <- c(n_stocks, n_years_model,n_ages)
	n_unk_n_years_n_ages <- c("waa")
	for(i in n_unk_n_years_n_ages) should_dims[[i]] <- c(input_data_dims[["waa"]][1], n_years_model,n_ages)
	n_stocks_n_regions_n_ages <- c("NAA_where", "M_re_index")
	for(i in n_stocks_n_regions_n_ages) should_dims[[i]] <- c(n_stocks, n_regions,n_ages)

	n_fleets_n_years_n_ages <- c("catch_paa")
	for(i in n_fleets_n_years_n_ages) should_dims[[i]] <- c(n_fleets, n_years_model,n_ages)
	n_indices_n_years_n_ages <- c("index_paa")
	for(i in n_indices_n_years_n_ages) should_dims[[i]] <- c(n_indices, n_years_model,n_ages)

	n_stocks_n_seasons_n_regions_n_regions <- c("can_move")
	for(i in n_stocks_n_seasons_n_regions_n_regions) should_dims[[i]] <- c(n_stocks, n_seasons,n_regions,n_regions)
	n_stocks_n_seasons_n_regions <- c("must_move")
	for(i in n_stocks_n_seasons_n_regions) should_dims[[i]] <- c(n_stocks, n_seasons,n_regions)
	n_stocks_n_seasons_n_regions_n_reg_m1 <- c("use_mu_prior","trans_mu_prior_sigma")
	for(i in n_stocks_n_seasons_n_regions_n_reg_m1) should_dims[[i]] <- c(n_stocks, n_seasons, n_regions, n_regions-1L)

	n_Ecov_n_stocks_n_ages_n_regions <- c("Ecov_how_M", "ind_Ecov_out_start_M", "ind_Ecov_out_end_M", "n_poly_Ecov_M")
	for(i in n_Ecov_n_stocks_n_ages_n_regions) should_dims[[i]] <- c(n_Ecov,n_stocks, n_ages,n_regions)
	n_Ecov_n_stocks_n_ages_n_seasons_n_regions_n_reg_m1 <- c("Ecov_how_mu", "ind_Ecov_out_start_mu", "ind_Ecov_out_end_mu", "n_poly_Ecov_mu")
	for(i in n_Ecov_n_stocks_n_ages_n_seasons_n_regions_n_reg_m1) should_dims[[i]] <- c(n_Ecov, n_stocks, n_ages, n_seasons, n_regions, n_regions-1L)

	n_stocks_n_proj_n_ages <- c("mature_proj", "waa_proj")
	for(i in n_stocks_n_proj_n_ages) should_dims[[i]] <- c(n_stocks, n_years_proj,n_ages)

	data_not_checked <- c(data_not_checked, n_stocks_n_proj_n_ages)

	should_data_dims <- should_dims
	# name <- "keep_Cpaa"
	# print(input_data_dims[[name]])
	# print(should_dims[[name]])
	# print(class(input_data_dims[[name]]))
	# print(class(should_dims[[name]]))
	# print(identical(should_dims[[name]], input_data_dims[[name]]))
	# stop()
	bad_data_dims <- sapply(names(should_data_dims), function(x) !identical(should_data_dims[[x]],input_data_dims[[x]]))
	bad_data_dims <- bad_data_dims[which(bad_data_dims)]
	# print(bad_data_dims)
	# stop()
	bad_data_dims <- bad_data_dims[which(!names(bad_data_dims) %in% data_not_checked)]
	print(bad_data_dims)
	#################################################
 
 	#pars
	input_par_dims <- lapply(input$par, dim)
	vecs <- which(sapply(input_par_dims, is.null))
	input_par_dims[vecs] <- lapply(input$par[vecs], length)
	# print(input_par_dims)
	# stop()

	#################################################

	should_dims <- lapply(1:length(input_par_dims), function(x) integer(0))
	names(should_dims) <- names(input_par_dims)

 	n_fleets_vec <- "log_catch_sig_scale"
	for(i in n_fleets_vec) should_dims[[i]] <- n_fleets
 	n_indices_vec <- c("logit_q", "q_prior_re", "log_index_sig_scale")
	for(i in n_indices_vec) should_dims[[i]] <- n_indices


 	#matrices

 	n_stocks_2 <- "mean_rec_pars"
	for(i in n_stocks_2) should_dims[[i]] <- c(n_stocks, 2L)
 	n_years_n_indices <- "q_re"
	for(i in n_years_n_indices) should_dims[[i]] <- c(n_years_model, n_indices)
 	n_indices_2 <- "q_repars"
	for(i in n_indices_2) should_dims[[i]] <- c(n_indices, 2L)
 	n_years_n_fleets <- "F_pars"
	for(i in n_years_n_fleets) should_dims[[i]] <- c(n_years_model, n_fleets)
 	n_blocks_n_pars <- c("logit_selpars")
	for(i in n_blocks_n_pars) should_dims[[i]] <- c(n_selblocks, n_ages+6L)
 	n_blocks_3 <- c("sel_repars")
	for(i in n_blocks_3) should_dims[[i]] <- c(n_selblocks, 3L)
	n_fleets_3 <- "catch_paa_pars"
	for(i in n_fleets_3) should_dims[[i]] <- c(n_fleets, 3L)
  n_indices_3 <- "index_paa_pars"
	for(i in n_indices_3) should_dims[[i]] <- c(n_indices, 3L)
	n_regions_3 <- "L_repars"
	for(i in n_regions_3) should_dims[[i]] <- c(n_regions, 3L)
	n_stocks_n_regions <- "log_b"
	for(i in n_stocks_n_regions) should_dims[[i]] <- c(n_stocks, n_regions)
	n_years_n_regions <- "L_re"
	for(i in n_years_n_regions) should_dims[[i]] <- c(n_years_model, n_regions)
	n_y_E_n_E <- c("Ecov_re", "Ecov_obs_logsigma", "Ecov_obs_logsigma_re")
	for(i in n_y_E_n_E) should_dims[[i]] <- c(n_years_Ecov, n_Ecov)
	three_n_E <- "Ecov_process_pars"
	for(i in three_n_E) should_dims[[i]] <- c(3L, n_Ecov)
	two_n_E <- "Ecov_obs_sigma_par"
	for(i in two_n_E) should_dims[[i]] <- c(2L, n_Ecov)
	n_py_n_stocks <- "logR_proj"
	for(i in n_py_n_stocks) should_dims[[i]] <- c(n_years_proj, n_stocks)

	par_not_checked <- n_py_n_stocks
	#arrays

	n_stocks_n_seasons_n_regions_n_reg_m1 <- c("mu_prior_re", "trans_mu")
	for(i in n_stocks_n_seasons_n_regions_n_reg_m1) should_dims[[i]] <- c(n_stocks, n_seasons, n_regions, n_regions-1L)
	n_stocks_n_seasons_n_regions_n_reg_m1_3 <- "mu_repars"
	for(i in n_stocks_n_seasons_n_regions_n_reg_m1_3) should_dims[[i]] <- c(n_stocks, n_seasons, n_regions, n_regions-1L, 3L)
	n_stocks_n_ages_n_seasons_n_y_n_regions_n_reg_m1 <- "mu_re"
	for(i in n_stocks_n_ages_n_seasons_n_y_n_regions_n_reg_m1) should_dims[[i]] <- c(n_stocks, n_ages, n_seasons, n_years_model, n_regions, n_regions-1L)
	n_stocks_n_regions_3 <- c("N1_repars", "trans_NAA_rho", "M_repars")
	for(i in n_stocks_n_regions_3) should_dims[[i]] <- c(n_stocks, n_regions, 3L)
	n_stocks_n_regions_n_ages <- c("log_N1", "log_NAA_sigma", "Mpars")
	for(i in n_stocks_n_regions_n_ages) should_dims[[i]] <- c(n_stocks, n_regions, n_ages)
	n_stocks_n_regions_n_y_m1_n_ages <- "log_NAA"
	for(i in n_stocks_n_regions_n_y_m1_n_ages) should_dims[[i]] <- c(n_stocks, n_regions, n_years_model-1L, n_ages)
	n_stocks_n_regions_n_y_n_ages <- "M_re"
	for(i in n_stocks_n_regions_n_y_n_ages) should_dims[[i]] <- c(n_stocks, n_regions, n_years_model, n_ages)
	n_blocks_n_y_n_ages <- "selpars_re"
	for(i in n_blocks_n_y_n_ages) should_dims[[i]] <- c(n_selblocks, n_years_model, n_ages)
	n_stocks_n_E_np <- "Ecov_beta_R"
	for(i in n_stocks_n_E_np) should_dims[[i]] <- c(n_stocks, n_Ecov, as.integer(max(input$data[["n_poly_Ecov_R"]])))
	n_stocks_n_ages_n_regions_n_E_np <- "Ecov_beta_M"
	for(i in n_stocks_n_ages_n_regions_n_E_np) should_dims[[i]] <- c(n_stocks, n_ages, n_regions, n_Ecov, as.integer(max(input$data[["n_poly_Ecov_M"]])))
	n_stocks_n_ages_n_seasons_n_regions_n_reg_m1_n_E_np <- "Ecov_beta_mu"
	for(i in n_stocks_n_ages_n_seasons_n_regions_n_reg_m1_n_E_np) should_dims[[i]] <- c(n_stocks, n_ages, n_seasons, n_regions, n_regions-1L, n_Ecov, as.integer(max(c(0,input$data[["n_poly_Ecov_mu"]]))))
	n_indices_n_E_np <- "Ecov_beta_q"
	for(i in n_indices_n_E_np) should_dims[[i]] <- c(n_indices, n_Ecov, as.integer(max(input$data[["n_poly_Ecov_q"]])))

	# name <- "Ecov_beta_mu"
	# print(input_par_dims[[name]])
	# print(should_dims[[name]])
	# print(class(input_par_dims[[name]]))
	# print(class(should_dims[[name]]))
	# print(identical(should_dims[[name]], input_par_dims[[name]]))

	should_par_dims <- should_dims
	bad_par_dims <- sapply(names(should_par_dims), function(x) !identical(should_par_dims[[x]],input_par_dims[[x]]))
	bad_par_dims <- bad_par_dims[which(bad_par_dims)]
	bad_par_dims <- bad_par_dims[which(!names(bad_par_dims) %in% par_not_checked)]
	return(list(bad_data_dims = bad_data_dims, bad_par_dims = bad_par_dims))

}