#define TMB_LIB_INIT R_init_wham
#include <TMB.hpp>
#include "all.hpp"


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  
  DATA_INTEGER(n_years_model);
  DATA_IVECTOR(years_use); //years to use for evaluating likelihoods (and simulating values). normally = 0,....,n_years_model-1. used for peels.
  DATA_INTEGER(n_seasons);
  DATA_VECTOR(fracyr_seasons); //length of intervals for seasons
  //int n_seasons = fracyr_seasons.size();
  DATA_INTEGER(n_regions);
  DATA_INTEGER(n_stocks);
  DATA_INTEGER(n_fleets);
  DATA_INTEGER(n_indices);
  DATA_INTEGER(n_ages);
  //DATA_IVECTOR(n_ages_fleet);
  //DATA_IVECTOR(n_ages_indices);
  DATA_IVECTOR(mig_type); //n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
  DATA_MATRIX(fracyr_SSB); //n_years x n_stocks:  size of interval from beginning of season to time of spawning within that season
  //int n_stocks = fracyr_SSB.cols();
  DATA_IVECTOR(spawn_regions); //n_stocks
  DATA_IVECTOR(spawn_seasons); //n_stocks
  //DATA_IARRAY(n_mu) //n_stocks x n_years x x n_seasons x n_ages
  //DATA_IVECTOR(mu_row) //length the total number of migration parameters
  //DATA_IVECTOR(mu_col) //length the total number of migration parameters
  //DATA_IVECTOR(mu_pointer) //n_stocks * n_years * n_ages
  //DATA_INTEGER(n_seasons_recruited); //easiest if this is = n_seasons 
  DATA_ARRAY(mature); //n_stocks x n_years x n_ages
  //int n_ages = mature.dim(2);
  DATA_IVECTOR(waa_pointer_fleets); //n_fleets indicator which waa to use for each fleet
  //int n_fleets = waa_pointer_fleets.size();
  //DATA_IVECTOR(waa_pointer_totcatch); //n_regions
  DATA_IVECTOR(waa_pointer_indices);
  //int n_indices = waa_pointer_indices.size();
  DATA_IVECTOR(waa_pointer_ssb); //n_stocks
  DATA_IVECTOR(waa_pointer_M); //n_stocks, possibly used with M_model = 2, M = a W^b
  DATA_ARRAY(waa); //(n_fleets + n_indices + n_stocks + 1(totcatch)) x n_years x n_ages_model
  //M is a parameter
  //DATA_ARRAY(MAA); // n_years x n_ages x n_regions
  
  DATA_IVECTOR(fleet_regions); //length = n_fleets
  DATA_IMATRIX(fleet_seasons); //n_fleets x n_seasons; indicator with length n_seasons, e.g., 1,1,1.. would mean operating all year.
  DATA_MATRIX(agg_catch); //n_years x n_fleets 
  DATA_IMATRIX(use_agg_catch); //n_years x n_fleets 
  DATA_MATRIX(agg_catch_sigma); //n_years x n_fleets 
  DATA_ARRAY(catch_paa); //n_fleets x n_years x n_ages_model 
  DATA_IMATRIX(use_catch_paa); //n_years x n_fleets 
  DATA_MATRIX(catch_Neff); //n_years x n_fleets 
  DATA_IVECTOR(age_comp_model_fleets); //length = n_fleets
  DATA_INTEGER(F_config); //1: F_pars is log_F1, F_devs, 2: F_pars is log_F.
  DATA_IVECTOR(index_regions); //length = n_indices
  DATA_IVECTOR(index_seasons); //n_indices x n_seasons. indicator must only be for one season
  DATA_IVECTOR(units_indices); //length = n_indices
  DATA_MATRIX(fracyr_indices); //n_years x n_indices: size of interval from beginning of season to time of surve within that season
  DATA_MATRIX(agg_indices); //n_years x n_indices
  DATA_IMATRIX(use_indices);  //n_years x n_indices
  DATA_MATRIX(agg_index_sigma); //n_years x n_indices
  DATA_IVECTOR(units_index_paa); //length = n_indices
  DATA_ARRAY(index_paa); //n_indices x n_years x n_ages_model
  DATA_IMATRIX(use_index_paa); //n_years x n_indices
  DATA_MATRIX(index_Neff); //n_years x n_indices
  DATA_IVECTOR(age_comp_model_indices); //length = n_indices

  // data for environmental covariate(s), Ecov
  DATA_INTEGER(n_Ecov); // also = 1 if no Ecov
  DATA_INTEGER(n_years_Ecov); // num years in Ecov  process model
  DATA_IVECTOR(years_use_Ecov); //years to use for evaluating likelihoods (and simulating values). normally = 0,....,n_years_Ecov-1. used for peels.
  DATA_IMATRIX(Ecov_use_obs); // all 0 if no Ecov
  DATA_MATRIX(Ecov_obs);
  DATA_IMATRIX(Ecov_how_R); // n_Ecov x n_stocks: specific to recruitment effects. 0 = no effect, 1 = controlling, 2 = limiting, 3 = lethal, 4 = masking, 5 = directive
  DATA_IARRAY(Ecov_how_M); // n_Ecov x n_stocks x n_ages x n_regions: 0/1 values indicating to use effects on natural mortality at age.
  DATA_IMATRIX(Ecov_how_q); // n_Ecov x n_indices: 0/1 values indicating to use effects on catchability for each index.
  DATA_IARRAY(Ecov_how_mu); // n_Ecov x n_stocks x n_ages x n_seasons x n_regions x n_regions-1: 0/1 values indicating to use effects on migration for each stock for each region (less 1).
  //DATA_IMATRIX(Ecov_where); // n_Ecov x 3+n_indices. 0/1 values with columns corresponding to recruit, mortality, migration, indices in that order
  DATA_IVECTOR(Ecov_model); // 0 = no Ecov, 1 = RW, 2 = AR1
  //DATA_IMATRIX(ind_Ecov_out_start); // n_Ecov x (2 + n_indices) index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects specific the multiple types of effects each Ecov can have)
  DATA_IMATRIX(ind_Ecov_out_start_R); // n_Ecov x n_stocks: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IARRAY(ind_Ecov_out_start_M); // n_Ecov x n_stocks x n_ages x n_regions: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IMATRIX(ind_Ecov_out_start_q); // n_Ecov x n_indices: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IARRAY(ind_Ecov_out_start_mu); // n_Ecov x n_stocks x n_ages x n_seasons x n_regions x n_regions-1: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  //DATA_IMATRIX(ind_Ecov_out_end); // n_Ecov x (2 + n_indices) index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects specific the multiple types of effects each Ecov can have)
  DATA_IMATRIX(ind_Ecov_out_end_R); // n_Ecov x n_stocks: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IARRAY(ind_Ecov_out_end_M); // n_Ecov x n_stocks x n_ages x n_regions: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IMATRIX(ind_Ecov_out_end_q); // n_Ecov x n_indices: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IARRAY(ind_Ecov_out_end_mu); // n_Ecov x n_stocks x n_ages x n_seasons x n_regions x n_regions-1: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IVECTOR(Ecov_obs_sigma_opt); // n_Ecov, 1 = given, 2 = estimate 1 value, shared among obs, 3 = estimate for each obs, 4 = estimate for each obs as random effects
  DATA_IMATRIX(n_poly_Ecov_R); // dim = n_ecov x n_stocks, order of orthogonal polynomial to use for effect of each covariate on each stock
  DATA_IARRAY(n_poly_Ecov_M); // dim = n_ecov x n_stocks x n_ages x n_regions, order of orthogonal polynomial to use for effect of each covariate on each stock
  DATA_IARRAY(n_poly_Ecov_mu); // dim = n_ecov x n_stocks x n_ages x n_seasons x n_regions x (n_regions-1), order of orthogonal polynomial to use for effect of each covariate on each stock
  DATA_IMATRIX(n_poly_Ecov_q); // dim = n_ecov x n_indices, order of orthogonal polynomial to use for effect of each covariate on each index
  
  DATA_IVECTOR(Ecov_use_re); // n_Ecov: 0/1: use Ecov_re? If yes, add to nll.

  DATA_VECTOR(q_lower); //length = n_indices
  DATA_VECTOR(q_upper); //length = n_indices
  DATA_IVECTOR(use_q_prior); //length = n_indices
  DATA_VECTOR(logit_q_prior_sigma); //length = n_indices
  DATA_IVECTOR(use_q_re);  //length = n_indices, 0= no re, >0 = use re
  
  DATA_IVECTOR(selblock_models);
  int n_selblocks = selblock_models.size();
  DATA_IVECTOR(selblock_models_re); // for each block: 1 = none, 2 = IID, 3 = ar1, 4 = ar1_y, 5 = 2dar1
  DATA_IVECTOR(n_selpars); //length = n_selblocks
  DATA_IMATRIX(selpars_est); // n_blocks x (n_pars + n_ages), is the selpar estimated in this block?
  DATA_IVECTOR(n_selpars_est); // of the selpars, how many are actually estimated (not fixed at 0 or 1)
  DATA_IVECTOR(n_years_selblocks); // for each block, number of years the block covers
  DATA_IMATRIX(selblock_years); // n_years_model x n_selblocks, = 1 if block covers year, = 0 if not
  DATA_IMATRIX(selblock_pointer_fleets);
  DATA_IMATRIX(selblock_pointer_indices);
  DATA_MATRIX(selpars_lower);
  DATA_MATRIX(selpars_upper);
  
  DATA_IVECTOR(recruit_model); //length = n_stocks
  DATA_IVECTOR(N1_model); //n_stocks, 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
  DATA_IVECTOR(NAA_re_model); //n_stocks, 0 SCAA, 1 "rec", 2 "rec+1"
  DATA_IARRAY(NAA_where); //n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
  //int n_regions = NAA_where.dim(1);
  //DATA_IVECTOR(n_M_a); //length = n_regions
  DATA_IMATRIX(n_M_re); // n_stocks x n_regions how many time-varying RE each year? n_ages? 1? n_est_M? max(n_M_re) <= n_ages)
  DATA_IARRAY(M_re_index); // n_stocks x n_regions x n_ages, indicators of which M_re to use for which age. length(unique(M_re_index[s,r,])) == n_M_re[s,r]
  DATA_INTEGER(M_model); // 1: M = f(age), 2: M = f(WAA)
  DATA_IMATRIX(M_re_model); //n_stocks x n_regions; 1 = none, 2 = ar1_a, 3 = ar1_y, 4 = 2dar1
  DATA_INTEGER(use_b_prior); //for M_model = 2: M = a*W^b model, use prior (and re) for log_b
  DATA_INTEGER(log_b_model); //1: constant, 2: differ by stock, 3: differ by region, 4: differ by both
  DATA_IVECTOR(L_model); //length = n_regions; L_model = 0: don't use L (extra mortality); L_model = 1: use constant L; L_model = 2: use iid re; L_model = 3: use ar1 re 
  DATA_IARRAY(can_move); //n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
  DATA_IARRAY(must_move); //n_stocks x n_seasons x n_regions: 0/1 determining if it must leave the region
  DATA_ARRAY(trans_mu_prior_sigma); //n_stocks x n_region(from) x n_regions-1 (to); sd for mu parameters on transformed (-inf,inf) scale for 
  DATA_IARRAY(use_mu_prior); //n_stocks x n_seasons x n_regions x n_regions-1: 0/1 whether to apply prior for each movement parameter
  DATA_IMATRIX(mu_model); //n_regions x n_regions - 1 
  // 1 = constant across stocks, ages, time (1 fixed effect for r,rr). 
  // 2 = differ by age (1 fixed effect, n_ages random effects for r,rr). 
  // 3 = differ by year, (1 fixed effect, n_years random effects for r,rr)
  // 4 = differ by age,year (1 fixed effect, n_years x n_ages random effects for r,rr)
  // 5 = differ by stock (n_stocks fixed effects for r,rr). 
  // 6 = differ by stock, age (n_stocks fixed effects, n_ages random effects for r,rr). 
  // 7 = differ by stock, year (n_stocks fixed effects, n_years random effects for r,rr)
  // 8 = differ by stock, age,year (n_stocks fixed effects, n_years x n_ages random effects for r,rr)
  // 9 = differ by season (n_seasons fixed effects for r,rr). 
  // 10 = differ by season,age (n_seasons fixed effects, n_ages random effects for r,rr). 
  // 11 = differ by season,year (n_seasons fixed effects, n_years random effects for r,rr)
  // 12 = differ by season,age,year (n_seasons fixed effects, n_years x n_ages random effects for r,rr)
  // 13 = differ by stock, season (n_stocks x n_seasons fixed effects for r,rr). 
  // 14 = differ by stock, season, age (n_stocks x n_seasons fixed effects, n_ages random effects for r,rr). 
  // 15 = differ by stock, season, year (n_stocks x n_seasons fixed effects, n_years random effects for r,rr)
  // 16 = differ by stock, season, age,year (n_stocks x n_seasons fixed effects, n_years x n_ages random effects for r,rr)
  //DATA_IARRAY(use_mu_re); //n_stocks x ages x n_seasons x n_years_model x n_regions x n_regions-1: 0/1 whether to use temporal RE for each movement parameter

  //DATA_IMATRIX(which_F_age); // (n_years_model + n_years_proj x 2); age, fleet for which F to use for max F for msy/ypr calculations and projections
  DATA_IVECTOR(which_F_age); // (n_years_model + n_years_proj); age for which F to use for max Fmsy/Fxspr calculations and projections
  //DATA_IVECTOR(which_F_season); // (n_years_model + n_years_proj); which season of F to use for max F for msy/ypr calculations and projections
  //DATA_IVECTOR(which_F_fleet); // (n_years_model + n_years_proj); which fleet of F to use for max F for msy/ypr calculations and projections
  //DATA_INTEGER(use_steepness); // which parameterization to use for BH/Ricker S-R, if needed.
  DATA_INTEGER(bias_correct_pe); //bias correct lognormal process error?
  DATA_INTEGER(bias_correct_oe); //bias correct lognormal observation error?
  DATA_IVECTOR(Fbar_ages);
  
  //DATA_IMATRIX(R1_pointer); //n_stocks x n_regions. Tells where log_R1 parameters are used.
  //DATA_IMATRIX(N1_sigma_pointer); //n_stocks x n_regions. Tells which N1_sigma_par to use where.
  //DATA_IARRAY(NAA_sigma_pointers); //n_stocks x n_ages x n_regions
  //DATA_IARRAY(NAA_re_indicator); //n_stocks x (n_years-1) x n_ages x n_regions will estimate a random effect where indicator is not 0. sum = n_NAA_re.

  //DATA_IVECTOR(simulate_state); //vector (0/1) if 1 then state parameters (NAA, MAA, sel, Ecov, q) in that order) will be simulated.
  //DATA_INTEGER(do_simulate_NAA_re); //(0/1) if 1 then simulate numbers at age random effects.
  DATA_INTEGER(do_simulate_Ecov_re); //(0/1) if 1 then simulate selectivity random effects.
  DATA_INTEGER(do_simulate_sel_re); //(0/1) if 1 then simulate selectivity random effects.
  DATA_INTEGER(do_simulate_M_re); //(0/1) if 1 then simulate natural mortality random effects.
  DATA_INTEGER(do_simulate_q_re); //(0/1) if 1 then simulate catchability random effects.
  DATA_INTEGER(do_simulate_q_prior_re); //(0/1) if 1 then simulate q prior random effects.
  DATA_INTEGER(do_simulate_mu_re); //(0/1) if 1 then simulate mu (migration) random effects.
  DATA_INTEGER(do_simulate_mu_prior_re); //(0/1) if 1 then simulate mu prior random effects.
  DATA_INTEGER(do_simulate_L_re); //(0/1) if 1 then simulate L (extra mortality) random effects.
  DATA_INTEGER(do_simulate_N_re); //(0/1) if 1 then simulate N1 and NAA random effects.
  DATA_IVECTOR(do_simulate_data); //vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.

  // data for one-step-ahead (OSA) residuals
  DATA_INTEGER(do_osa); //whether to do osa residuals. For efficiency reasons with age comp likelihoods.
  DATA_VECTOR(obsvec); // vector of all observations for OSA residuals
  DATA_IVECTOR(agesvec); // vector of ages associated with paa observations for OSA residuals. same length as obsvec
  DATA_VECTOR_INDICATOR(keep, obsvec); // for OSA residuals
  DATA_IMATRIX(keep_C); // indices for catch obs, can loop years/fleets with keep(keep_C(y,f))
  DATA_IMATRIX(keep_I);
  DATA_IMATRIX(keep_E); // Ecov
  DATA_IARRAY(keep_Cpaa);
  DATA_IARRAY(keep_Ipaa);
  DATA_INTEGER(do_post_samp_N); //whether to ADREPORT posterior residuals for NAA re. 
  DATA_INTEGER(do_post_samp_M); //whether to ADREPORT posterior residuals for M re. 
  DATA_INTEGER(do_post_samp_mu); //whether to ADREPORT posterior residuals for movement re. 
  DATA_INTEGER(do_post_samp_sel); //whether to ADREPORT posterior residuals for selectivity re. 
  DATA_INTEGER(do_post_samp_Ecov); //whether to ADREPORT posterior residuals for Ecov re. 
  DATA_INTEGER(do_post_samp_q); //whether to ADREPORT posterior residuals for q re. 
  int sum_do_post_samp = do_post_samp_N + do_post_samp_M + do_post_samp_mu + do_post_samp_sel + do_post_samp_Ecov + do_post_samp_q;
  //reference points
  //DATA_IVECTOR(do_simulate_period); //vector (0/1) if 1 then period (model years, projection years) will be simulated.
  DATA_INTEGER(do_SPR_BRPs); //whether to calculate and adreport reference points. 
  DATA_INTEGER(do_MSY_BRPs); //whether to calculate and adreport reference points. 
  DATA_INTEGER(SPR_weight_type); //0 = use average recruitment for each stock for weighting, 1= use SPR_weights 
  DATA_VECTOR(SPR_weights); //n_stocks; weights to use for stock-specific SPR for weight some. should sum to 1.
  DATA_SCALAR(percentSPR); // percentage to use for SPR-based reference points. Default = 40.
  DATA_IVECTOR(XSPR_R_avg_yrs); // model year indices (TMB, starts @ 0) to use for averaging recruitment when defining SSB_XSPR (if XSPR_R_opt = 2,4)
  DATA_VECTOR(FXSPR_init); // annual initial values to use for newton steps to find FXSPR (n_years_model+n_proj_years)
  DATA_VECTOR(FMSY_init); // annual initial values to use for newton steps to find FMSY (n_years_model+n_proj_years)
  DATA_INTEGER(n_regions_is_small) //is the number of regions "small"? determines different matrix inversion methods in TMB
  
  //static brp info
  DATA_SCALAR(FXSPR_static_init); // initial value to use for newton steps to find FXSPR_static
  DATA_SCALAR(FMSY_static_init); // initial value to use for newton steps to find FXSPR_static
  DATA_INTEGER(which_F_age_static); // which age,fleet of F to use for full total F for static brps (max of average FAA_tot over avg_years_ind)
  DATA_INTEGER(XSPR_R_opt); //1(3): use annual R estimates(predictions) for annual SSB_XSPR, 2(4): use average R estimates(predictions). See XSPR_R_avg_yrs for years to average over.
  
  DATA_INTEGER(use_alt_AR1) //0: use density namespace, 1: use nll calculated by "hand".
  
  // data for projections
  DATA_INTEGER(n_years_proj); // number of years to project  
  DATA_IVECTOR(avg_years_ind); // model year indices (TMB, starts @ 0) to use for averaging MAA, waa, maturity, and F (if use.avgF = TRUE)
  DATA_IVECTOR(proj_Ecov_opt); // if any, how to use each ecov in pop projections: 1 = continue Ecov_re, 2 = average Ecov (over avg_years_ind), 3 = terminal year Ecov, 4 = user-specified
  DATA_MATRIX(Ecov_use_proj); // n_years_proj x n_Ecov matrix of fixed user-supplied values to use in projections if proj_Ecov_opt = 4
  DATA_IVECTOR(avg_years_Ecov); // model year indices (TMB, starts @ 0) to use for averaging ecov for projections if proj_Ecov_opt = 2
  DATA_IVECTOR(proj_F_opt); // for each projection year, 1 = last year F (default), 2 = average F, 3 = F at X% SPR, 4 = user-specified F, 5 = calculate F from user-specified catch
  DATA_VECTOR(proj_Fcatch); // user-specified F or catch in projection years, only used if proj_F_opt = 4 or 5
  DATA_INTEGER(proj_M_opt); // 1 = continue M_re (check for time-varying M_re on R side), 2 = average M (over avg_years_ind)
  DATA_INTEGER(proj_mu_opt); // 1 = continue mu_re (check for time-varying M_re on R side), 2 = average mu (over avg_years_ind)
  DATA_INTEGER(proj_L_opt); // 1 = continue mu_re (check for time-varying M_re on R side), 2 = average mu (over avg_years_ind)
  DATA_VECTOR(logR_mean); // (n_stocks) empirical mean recruitment in model years, used for SCAA recruit projections
  DATA_VECTOR(logR_sd); //  (n_stocks) empirical sd recruitment in model years, used for SCAA recruit projections
  DATA_VECTOR(F_proj_init); // annual initial values  to use for newton steps to find F for use in projections  (n_years_proj)
  DATA_SCALAR(percentFMSY); // percent of FMSY to use for calculating catch in projections.
  DATA_SCALAR(percentFXSPR); // percent of F_XSPR to use for calculating catch in projections. For example, GOM cod uses F = 75% F_40%SPR, so percentFXSPR = 75 and percentSPR = 40. Default = 100.
  

  // parameters - general
  PARAMETER_MATRIX(mean_rec_pars); //n_stocks x 2
  PARAMETER_VECTOR(logit_q);
  PARAMETER_VECTOR(q_prior_re); //n_indices (if a prior is used for q, this is used instead of logit_q)
  PARAMETER_MATRIX(q_re); //n_years x n_indices (time series of)
  PARAMETER_MATRIX(q_repars) //n_indices x 2 (sigma, rho)
  PARAMETER_MATRIX(F_pars); //n_years_model x n_fleets
  //PARAMETER_MATRIX(F_devs);
  //PARAMETER_VECTOR(mu); //migration parameters
  PARAMETER_ARRAY(mu_prior_re); //n_stocks x n_seasons x n_regions x n_regions-1
  PARAMETER_ARRAY(trans_mu); //n_stocks x n_seasons x n_regions x n_regions-1 (mean) migration parameters
  PARAMETER_ARRAY(mu_re); //n_stocks x n_ages x n_seasons x n_y x n_regions x n_regions-1 RE for migration
  PARAMETER_ARRAY(mu_repars); //n_stocks x n_seasons x n_regions x n_regions-1 x 3 (sig, rho_a, rho_y)
  //N1 might need some tweaking. for example, if all fish are forced to be in spawning region at the beginning of the year, then there should be no N1 in other regions.
  //PARAMETER_VECTOR(log_R1); //length must be consistent with R1_pointer above.
  PARAMETER_ARRAY(N1_repars); // (n_stocks x n_regions x 3) mean, sig, rho
  PARAMETER_ARRAY(log_N1); // (n_stocks x n_regions x n_ages)
  PARAMETER_ARRAY(log_NAA_sigma); // (n_stocks x n_regions x n_ages) sigmas for NAA RE
  PARAMETER_ARRAY(trans_NAA_rho); // (n_stocks x n_regions x 3) rho_a, rho_y, recruits rho_y
  //Just have annual NAA currently
  PARAMETER_ARRAY(log_NAA); //(n_stocks x n_regions x nyears-1 x n_ages) 
  
  PARAMETER_MATRIX(logR_proj); // (n_proj_years x n_stocks) recruitment (random effects) in proj years, only if SCAA
  PARAMETER_MATRIX(logit_selpars); // mean selectivity, dim = n_selblocks x n_ages + 6 (n_ages for by age, 2 for logistic, 4 for double-logistic)
  PARAMETER_VECTOR(selpars_re);    // deviations in selectivity parameters (random effects), length = sum(n_selpars)*n_years per block
  PARAMETER_MATRIX(sel_repars);    // fixed effect parameters controlling selpars_re, dim = n_blocks, 3 (sigma, rho, rho_y)
  PARAMETER_MATRIX(catch_paa_pars); //n_fleets x 3
  PARAMETER_MATRIX(index_paa_pars); //n_indices x 3
  PARAMETER_VECTOR(log_catch_sig_scale); //n_fleets
  PARAMETER_VECTOR(log_index_sig_scale); //n_indices
  
  PARAMETER_ARRAY(Mpars); // (n_stocks x n_regions x n_ages) mean log M-at-age, fixed effects
  PARAMETER_ARRAY(M_re); // random effects for year- and age-varying M deviations from mean Mpars, dim = n_stocks x n_regions x n_years x n_ages
  PARAMETER_ARRAY(M_repars); // parameters controlling M_re, (n_stocks x n_regions x length = 3 (sigma_M, rho_M_a, rho_M_y))
  PARAMETER_MATRIX(log_b); //n_stocks x n_regions (for M = a * W^b model)
  PARAMETER_MATRIX(L_re); // random effects for year- and region-varying extra/unknown mortality source, dim = n_years x n_regions
  PARAMETER_MATRIX(L_repars); // parameters controlling L_re, dim = n_regions x 3 (mu, sigma_L, rho_L_y)

  // parameters - environmental covariate ("Ecov")
  PARAMETER_MATRIX(Ecov_re); // nrows = n_years_Ecov, ncol = N_Ecov
  //PARAMETER_ARRAY(Ecov_beta); // dim = (2 + n_stocks*n_regions + n_indices) x n_poly x n_ecov x n_ages
  PARAMETER_ARRAY(Ecov_beta_R); // dim = n_stocks x n_ecov x max(n_poly_R), effects on recruitment, beta_R in eqns 4-5, Miller et al. (2016)
  PARAMETER_ARRAY(Ecov_beta_M); // dim = n_stocks x n_ages x n_regions x n_ecov x max(n_poly_M), effects on natural mortality
  PARAMETER_ARRAY(Ecov_beta_mu); // dim = n_stocks x n_ages x n_seasons x n_regions x (n_regions-1) x n_ecov x max(n_poly_mu), effects on movement
  PARAMETER_ARRAY(Ecov_beta_q); // dim = n_indices x n_ecov x max(n_poly_q), effects on catchability
  PARAMETER_MATRIX(Ecov_process_pars); // nrows = RW: 2 par (Ecov1, sig), AR1: 3 par (mu, sig, phi); ncol = N_ecov
  PARAMETER_MATRIX(Ecov_obs_logsigma); // N_Ecov_years x n_Ecov. options: just given (data), or fixed effect(s)
  PARAMETER_MATRIX(Ecov_obs_logsigma_re); // N_Ecov_years x n_Ecov. columns of random effects used if Ecov_obs_sigma_opt = 4 
  PARAMETER_MATRIX(Ecov_obs_sigma_par); // ncol = N_Ecov, nrows = 2 (mean, sigma of random effects)
  Type nll = 0.0; //negative log-likelihood
  //int trace = 1;
  int no_trace = 0;
  int n_years_pop = n_years_model + n_years_proj;

  //make expanded (or not) fracyr_SSB, waa, maturity
  matrix<Type> fracyr_SSB_all(n_years_pop, n_stocks);
  fracyr_SSB_all.setZero();
  for(int y = 0; y < n_years_model; y ++) fracyr_SSB_all.row(y) = vector<Type> (fracyr_SSB.row(y));
  array<Type> mature_all(n_stocks,n_years_pop, n_ages);
  mature_all.setZero();
  array<Type> waa_ssb(n_stocks,n_years_pop,n_ages);
  waa_ssb.setZero();
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years_model; y ++) for(int a = 0; a < n_ages; a++) {
    mature_all(s,y,a) = mature(s,y,a); 
    waa_ssb(s,y,a) = waa(waa_pointer_ssb(s)-1,y,a);
  }
  array<Type> waa_catch(n_fleets,n_years_pop,n_ages);
  waa_catch.setZero();
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_years_model; y ++) for(int a = 0; a < n_ages; a++) {
    waa_catch(f,y,a) = waa(waa_pointer_fleets(f)-1,y,a);
  }

  /////////////////////////////////////////
  // Environmental covariate process model --------------------------------------

  // 'true' estimated Ecov (x_t in Miller et al. 2016 CJFAS)
  matrix<Type> Ecov_x = get_Ecov(Ecov_model, Ecov_process_pars, Ecov_re, Ecov_use_re);
  if(Ecov_model.sum()>0) {
    matrix<Type> nll_Ecov = get_nll_Ecov(Ecov_model, Ecov_process_pars, Ecov_re, Ecov_use_re, years_use_Ecov);
    nll += nll_Ecov.sum();
    REPORT(nll_Ecov);
    SIMULATE if(do_simulate_Ecov_re == 1){
      Ecov_re = simulate_Ecov_re(Ecov_model, Ecov_process_pars, Ecov_re, Ecov_use_re, years_use_Ecov);
      Ecov_x = get_Ecov(Ecov_model, Ecov_process_pars, Ecov_re, Ecov_use_re);
      REPORT(Ecov_re);
    }
    if(do_post_samp_Ecov == 1) ADREPORT(Ecov_re);
  }
  REPORT(Ecov_x);
  if(do_post_samp_Ecov == 1) ADREPORT(Ecov_re);
  /////////////////////////////////////////
  //see(Ecov_x);

  /////////////////////////////////////////
  // Environmental covariate observation model -------------------------------------
  //TODO: Ecov obs are not yet simulated in projection years!!!!!!!!
  matrix<Type> nll_Ecov_obs(n_years_Ecov, n_Ecov);
  nll_Ecov_obs.setZero();
  //Type nll_Ecov_obs = Type(0);
  Type nll_Ecov_obs_sig = Type(0); // Ecov obs sigma random effects (opt = 4)
  matrix<Type> Ecov_obs_sigma(n_years_Ecov, n_Ecov);
  for(int i = 0; i < n_Ecov; i++){
    for(int y = 0; y < n_years_Ecov; y++){
      if(Ecov_obs_sigma_opt(i) == 4){
        Type mu_logsigma = Ecov_obs_sigma_par(0,i);
        Type sd_logsigma = exp(Ecov_obs_sigma_par(1,i));
        nll_Ecov_obs_sig -= dnorm(Ecov_obs_logsigma_re(y,i), mu_logsigma, sd_logsigma, 1);
        SIMULATE if(do_simulate_data(2)) {
          Ecov_obs_logsigma_re(y,i) = rnorm(mu_logsigma, sd_logsigma);
        }
        Ecov_obs_sigma(y,i) = exp(Ecov_obs_logsigma_re(y,i));
      } else{
        Ecov_obs_sigma(y,i) = exp(Ecov_obs_logsigma(y,i));
      }
      if(Ecov_use_obs(y,i) == 1){
        nll_Ecov_obs(y,i) -= keep(keep_E(y,i)) * dnorm(obsvec(keep_E(y,i)), Ecov_x(y,i), Ecov_obs_sigma(y,i), 1);
        nll_Ecov_obs(y,i) -= keep.cdf_lower(keep_E(y,i)) * log(squeeze(pnorm(obsvec(keep_E(y,i)), Ecov_x(y,i), Ecov_obs_sigma(y,i))));
        nll_Ecov_obs(y,i) -= keep.cdf_upper(keep_E(y,i)) * log(1.0 - squeeze(pnorm(obsvec(keep_E(y,i)), Ecov_x(y,i), Ecov_obs_sigma(y,i))));
        SIMULATE if(do_simulate_data(2)) {
          Ecov_obs(y,i) = rnorm(Ecov_x(y,i), Ecov_obs_sigma(y,i));
          obsvec(keep_E(y,i)) = Ecov_obs(y,i);
        }
      }
    }
  }
  nll += nll_Ecov_obs_sig;
  nll += nll_Ecov_obs.sum();
  REPORT(nll_Ecov_obs);
  REPORT(nll_Ecov_obs_sig);
  REPORT(Ecov_obs_sigma);
  SIMULATE if(do_simulate_data(2)) {
    REPORT(Ecov_obs);
    REPORT(Ecov_obs_logsigma);
  }
  if(Ecov_model.sum() > 0){
    matrix<Type> Ecov_resid = Ecov_obs.array() - Ecov_x.block(0,0,n_years_Ecov,n_Ecov).array();
    if(sum_do_post_samp == 0){
      ADREPORT(Ecov_x);
      ADREPORT(Ecov_resid);
    }
  }

  /////////////////////////////////////////
  
    
  /////////////////////////////////////////////////////////
  //////next is setting up Ecov_out, Ecov_lm for R, M, mu, q
  /////////////////////////////////////////////////////////
  // Lag environmental covariates -------------------------------------
  // Then use Ecov_out_*(t) for processes in year t, instead of Ecov_x
  //Recruit
  array<Type> Ecov_out_R(n_stocks, n_years_pop, n_Ecov);
  Ecov_out_R.setZero();
  for(int s = 0; s < n_stocks; s++){
    vector<int> t_ind_s = ind_Ecov_out_start_R.col(s);
    vector<int> t_ind_e = ind_Ecov_out_end_R.col(s);
    //see(t_ind_s);
    //see(t_ind_e);
    matrix<Type> tmp = get_Ecov_out(Ecov_x, n_years_model, n_years_proj, t_ind_s, t_ind_e, proj_Ecov_opt, avg_years_Ecov, Ecov_use_proj);
    //see(tmp);
    for(int j = 0; j < tmp.rows(); j++) for(int k = 0; k < n_Ecov; k++) Ecov_out_R(s,j,k) = tmp(j,k);
  }
  REPORT(Ecov_out_R); 
  //see(Ecov_out_R);
  array<Type> Ecov_lm_R(n_stocks, n_years_pop, n_Ecov); 
  int max_n_poly_R = Ecov_beta_R.dim(2); // now a 3D array dim: (n_stocks,n_Ecov,max(n_poly_Ecov_R))
  //see(max_n_poly_R);
  for(int s = 0; s < n_stocks; s++) {
    matrix<Type> Ecov_beta_R_s(n_Ecov,max_n_poly_R);
    matrix<Type> Ecov_out_R_s(n_years_pop,n_Ecov);
    vector<int> n_poly_Ecov_R_s(n_Ecov);
    Ecov_beta_R_s.setZero();
    for(int i = 0; i <n_Ecov; i++) {
      n_poly_Ecov_R_s(i) = n_poly_Ecov_R(i,s);
      //see(n_poly_Ecov_R_s(i));
      //see(n_years_pop);
      //see(Ecov_out_R.dim);
      //see(Ecov_beta_R_s);
      for(int y = 0; y < n_years_pop; y++) Ecov_out_R_s(y,i) = Ecov_out_R(s,y,i);
      //see("here");
      for(int j = 0; j < max_n_poly_R; j++) Ecov_beta_R_s(i,j) = Ecov_beta_R(s,i,j);
    }
    //see(n_poly_Ecov_R_s);
    //see(Ecov_beta_R_s);
    //see(Ecov_out_R_s);
    //see("after");
    matrix<Type> Ecov_lm_R_s = get_Ecov_lm(Ecov_beta_R_s,Ecov_out_R_s, n_years_model, n_years_proj, n_poly_Ecov_R_s);
    for(int y = 0; y < n_years_pop; y++) for(int i = 0; i <n_Ecov; i++) Ecov_lm_R(s,y,i) = Ecov_lm_R_s(y,i);
  }
  //see(Ecov_lm_R)
  REPORT(Ecov_lm_R);
  ///////////////////////

  //M
  array<Type> Ecov_out_M(n_stocks, n_ages, n_regions, n_years_pop, n_Ecov);
  Ecov_out_M.setZero();
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++){
    vector<int> t_ind_s(n_Ecov), t_ind_e(n_Ecov);
    for(int i = 0; i < n_Ecov; i++){
      t_ind_s(i) = ind_Ecov_out_start_M(i,s,a,r);
      t_ind_e(i) = ind_Ecov_out_end_M(i,s,a,r);
    }
    matrix<Type> tmp = get_Ecov_out(Ecov_x, n_years_model, n_years_proj, t_ind_s, t_ind_e, proj_Ecov_opt, avg_years_Ecov, Ecov_use_proj);
    for(int j = 0; j < tmp.rows(); j++) for(int i = 0; i < n_Ecov; i++) Ecov_out_M(s,a,r,j,i) = tmp(j,i);
  }
  REPORT(Ecov_out_M);
  array<Type> Ecov_lm_M(n_stocks, n_regions, n_ages, n_years_pop, n_Ecov); 
  int max_n_poly_M = Ecov_beta_M.dim(4); // now a 5D array dim: (n_stocks, n_ages, n_regions, n_Ecov, max(n_poly_Ecov_M))
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++){
    matrix<Type> Ecov_beta_M_s_r_a(n_Ecov,max_n_poly_M);
    matrix<Type> Ecov_out_M_s_r_a(n_years_pop, n_Ecov);
    vector<int> n_poly_Ecov_M_s_r_a(n_Ecov);
    for(int i = 0; i <n_Ecov; i++) {
      n_poly_Ecov_M_s_r_a(i) = n_poly_Ecov_M(i,s,a,r);
      for(int y = 0; y < n_years_pop; y++) Ecov_out_M_s_r_a(y,i) = Ecov_out_M(s,a,r,y,i);
      for(int j = 0; j < max_n_poly_M; j++) Ecov_beta_M_s_r_a(i,j) = Ecov_beta_M(s,a,r,i,j);
    }
    matrix<Type> Ecov_lm_M_s_r_a = get_Ecov_lm(Ecov_beta_M_s_r_a,Ecov_out_M_s_r_a, n_years_model, n_years_proj, n_poly_Ecov_M_s_r_a);
    for(int y = 0; y < n_years_pop; y++) for(int i = 0; i <n_Ecov; i++) Ecov_lm_M(s,r,a,y,i) = Ecov_lm_M_s_r_a(y,i);
  }
  REPORT(Ecov_lm_M);
  
  //q
  array<Type> Ecov_out_q(n_indices, n_years_pop, n_Ecov);
  for(int i = 0; i < n_indices; i++){
    vector<int> t_ind_s = ind_Ecov_out_start_q.col(i);
    vector<int> t_ind_e = ind_Ecov_out_end_q.col(i);
    matrix<Type> tmp = get_Ecov_out(Ecov_x,n_years_model, n_years_proj, t_ind_s, t_ind_e, proj_Ecov_opt, avg_years_Ecov, Ecov_use_proj);
    for(int j = 0; j < tmp.rows(); j++) for(int k = 0; k < n_Ecov; k++) Ecov_out_q(i,j,k) = tmp(j,k);
  }
  //see("q");
  REPORT(Ecov_out_q);
  array<Type> Ecov_lm_q(n_indices,n_years_pop, n_Ecov); 
  int max_n_poly_q = Ecov_beta_q.dim(2); // now a 3D array dim: (n_indices, n_Ecov, n_poly)
  //see(Ecov_lm_q.dim);
  //see(Ecov_out_q.dim);
  //see(n_years_pop);
  for(int i = 0; i < n_indices; i++){
  //see("q1");
    matrix<Type> Ecov_beta_q_i(n_Ecov,max_n_poly_q);
    matrix<Type> Ecov_out_q_i(n_years_pop, n_Ecov);
    vector<int> n_poly_Ecov_q_i(n_Ecov);
    for(int k = 0; k <n_Ecov; k++) {
      n_poly_Ecov_q_i(k) = n_poly_Ecov_q(k,i); 
      for(int y = 0; y < n_years_pop; y++) Ecov_out_q_i(y,k) = Ecov_out_q(i,y,k);
      for(int j = 0; j < max_n_poly_q; j++) Ecov_beta_q_i(k,j) = Ecov_beta_q(i,k,j);
    }
    //see("q2");
    //see(n_poly_Ecov_q_i);

    matrix<Type> Ecov_lm_q_i = get_Ecov_lm(Ecov_beta_q_i,Ecov_out_q_i, n_years_model, n_years_proj, n_poly_Ecov_q_i);
    //see(Ecov_lm_q_i.rows());
    //see(Ecov_lm_q_i.cols());
    //see("q3");
    //see(Ecov_lm_q_i);
    for(int y = 0; y < n_years_pop; y++) for(int k = 0; k < n_Ecov; k++) Ecov_lm_q(i,y,k) = Ecov_lm_q_i(y,k);
    //see("q4");
  }

  REPORT(Ecov_lm_q);
  
  //mu
  array<Type> Ecov_out_mu(n_stocks, n_ages, n_seasons, n_regions, n_regions-1, n_years_pop, n_Ecov);
  for(int s = 0; s < n_stocks; s++)  for(int a = 0; a < n_ages; a++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions-1; rr++){
    vector<int> t_ind_s(n_Ecov), t_ind_e(n_Ecov);
    for(int i = 0; i < n_Ecov; i++) {
      t_ind_s(i) = ind_Ecov_out_start_mu(i,s,a,t,r,rr);
      t_ind_e(i) = ind_Ecov_out_end_mu(i,s,a,t,r,rr);
    }
    matrix<Type> tmp = get_Ecov_out(Ecov_x, n_years_model, n_years_proj, t_ind_s, t_ind_e, proj_Ecov_opt, avg_years_Ecov, Ecov_use_proj);
    for(int j = 0; j < tmp.rows(); j++) for(int i = 0; i < n_Ecov; i++) Ecov_out_mu(s,a,t,r,rr,j,i) = tmp(j,i);
  }
  REPORT(Ecov_out_mu);
  array<Type> Ecov_lm_mu(n_stocks, n_ages, n_seasons, n_regions, n_regions-1, n_years_pop, n_Ecov);
  int max_n_poly_mu = Ecov_beta_mu.dim(6); // a 7D array dim: (n_stocks, n_ages, n_seasons, n_regions, n_regions-1, n_Ecov, max(n_poly))
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions-1; rr++){
    matrix<Type> Ecov_beta_mu_s_a_t_r_rr(n_Ecov,max_n_poly_mu);
    matrix<Type> Ecov_out_mu_s_a_t_r_rr(n_years_pop, n_Ecov);
    vector<int> n_poly_Ecov_mu_s_a_t_r_rr(n_Ecov);
    for(int i = 0; i <n_Ecov; i++) {
      n_poly_Ecov_mu_s_a_t_r_rr(i) = n_poly_Ecov_mu(i,s,a,t,r,rr);
      for(int y = 0; y < n_years_pop; y++) Ecov_out_mu_s_a_t_r_rr(y,i) = Ecov_out_mu(s,a,t,r,rr,y,i);
      for(int j = 0; j <max_n_poly_mu; j++) Ecov_beta_mu_s_a_t_r_rr(i,j) = Ecov_beta_mu(s,a,t,r,rr,i,j);
    }
    matrix<Type> Ecov_lm_s_a_t_r_rr = get_Ecov_lm(Ecov_beta_mu_s_a_t_r_rr,Ecov_out_mu_s_a_t_r_rr, n_years_model, n_years_proj, n_poly_Ecov_mu_s_a_t_r_rr);
    for(int y = 0; y < n_years_pop; y++) for(int i = 0; i <n_Ecov; i++) Ecov_lm_mu(s,a,t,r,rr,y,i) = Ecov_lm_s_a_t_r_rr(y,i);
  }
  REPORT(Ecov_lm_mu);
  /////////////////////////////////////////
  
  /////////////////////////////////////////
  // Selectivity --------------------------------------------------------------
  if(selblock_models_re.sum()>0) {
    vector<Type> nll_sel = get_nll_sel(selblock_models_re, n_years_selblocks, n_selpars_est, selpars_re, sel_repars);
    nll += nll_sel.sum();
    REPORT(nll_sel);
    SIMULATE if(do_simulate_sel_re){
      selpars_re = simulate_selpars_re(selblock_models_re, n_years_selblocks, n_selpars_est, selpars_re, sel_repars);
      REPORT(selpars_re);
    }
    if(do_post_samp_sel) ADREPORT(selpars_re);
  }
  vector<matrix<Type> > selpars_re_mats = get_selpars_re_mats(n_selpars, selblock_years, selpars_est, 
    n_years_model, selpars_re, selblock_models, selblock_models_re);
  REPORT(selpars_re_mats); //can't report a vector<array<Type>> ?

  vector<matrix<Type> > selpars = get_selpars(selblock_models, n_selpars, logit_selpars, 
    selpars_re_mats, selpars_lower, selpars_upper, n_years_model);
  REPORT(selpars);

  vector<matrix<Type> > selAA = get_selAA(n_years_model, n_ages, n_selblocks, selpars, selblock_models);
  REPORT(selAA);
  /////////////////////////////////////////
 
  /////////////////////////////////////////
  //extra mortality parameter: missing catch?
  if(L_model.sum()>0) {
    vector<Type> nll_L = get_nll_L(L_model, L_repars, L_re);
    nll += nll_L.sum();
    REPORT(nll_L);
    SIMULATE if(do_simulate_L_re==1){
      L_re = simulate_L_re(L_model, L_repars, L_re);
      REPORT(L_re);
    }
  }
  //n_years_model x n_regions
  matrix<Type> L = get_L(L_model, L_repars, L_re, n_years_model, n_years_proj, proj_L_opt, avg_years_ind);
  REPORT(L);
  /////////////////////////////////////////

  /////////////////////////////////////////
  //catchability
  if(use_q_prior.sum()>0) {
    vector<Type> nll_q_prior = get_nll_q_prior(q_prior_re, logit_q, logit_q_prior_sigma, use_q_prior);
    nll += nll_q_prior.sum();
    REPORT(nll_q_prior);
    SIMULATE if(do_simulate_q_prior_re==1){
      q_prior_re = simulate_q_prior_re(q_prior_re, logit_q, logit_q_prior_sigma, use_q_prior);
      REPORT(q_prior_re);
    }
  }
  if(use_q_re.sum()>0) {
    matrix<Type> nll_q_re = get_nll_q_re(q_repars, q_re, use_q_re, years_use);
    nll += nll_q_re.sum();
    REPORT(nll_q_re);
    SIMULATE if(do_simulate_q_re ==1){
      q_re = simulate_q_re(q_repars, q_re, use_q_re, years_use);
    }
  }
  if(do_post_samp_q) ADREPORT(q_re);
  
  matrix<Type> logit_q_mat = get_logit_q_mat(logit_q, q_re, q_prior_re, use_q_prior, use_q_re, Ecov_lm_q, Ecov_how_q);
  REPORT(logit_q_mat);
  
  if((use_q_re.sum()>0) | (Ecov_how_q.sum() > 0)) if(sum_do_post_samp == 0) ADREPORT(logit_q_mat);
  
  //n_years_pop x n_indices;
  matrix<Type> q = get_q(logit_q_mat, q_lower, q_upper);
  REPORT(q);
  //n_indices x n_years_pop x n_ages;
  array<Type> QAA = get_QAA(q,selAA, selblock_pointer_indices, n_years_model, n_ages);
  REPORT(QAA);
  /////////////////////////////////////////
  
  /////////////////////////////////////////
  //natural mortality 
  //RE and log_M (possibly updated in time steps)
  matrix<Type> nll_M = get_nll_M(M_repars, M_re_model, M_model, M_re, n_M_re, years_use);
  nll += nll_M.sum();
  REPORT(nll_M);
  SIMULATE if(do_simulate_M_re){
    M_re = simulate_M_re(M_repars, M_re_model, M_model, M_re, n_M_re, years_use);
    REPORT(M_re);
  }
  if(do_post_samp_M) ADREPORT(M_re);

  //log_b prior for M_model = 2
  if((M_model == 2) & use_b_prior){
    matrix<Type> nll_log_b = get_nll_log_b(log_b_model, log_b, bias_correct_pe);
    nll += nll_log_b.sum();
    REPORT(nll_log_b);
    SIMULATE if(do_simulate_M_re){
      log_b = simulate_log_b(log_b_model, log_b, bias_correct_pe);
      REPORT(log_b);
    }
    if(do_post_samp_M) ADREPORT(log_b);
  }
  
  //n_stocks x n_regions x n_years x n_ages
  array<Type> log_M = get_log_M(M_re, M_re_index, M_model, n_years_model, Mpars, log_b, waa, waa_pointer_M, Ecov_lm_M, Ecov_how_M, n_years_proj, 
    proj_M_opt, avg_years_ind);
  array<Type> MAA = get_MAA(log_M);
  REPORT(log_M);
  REPORT(MAA);
  /////////////////////////////////////////
  
  /////////////////////////////////////////
  //movement
  // mu_model: 
  // 1 = constant across stocks, ages, time (1 fixed effect for r,rr). 
  // 2 = differ by age (1 fixed effect, n_ages random effects for r,rr). 
  // 3 = differ by year, (1 fixed effect, n_years random effects for r,rr)
  // 4 = differ by age,year (1 fixed effect, n_years x n_ages random effects for r,rr)
  // 5 = differ by stock (n_stocks fixed effects for r,rr). 
  // 6 = differ by stock, age (n_stocks fixed effects, n_ages random effects for r,rr). 
  // 7 = differ by stock, year (n_stocks fixed effects, n_years random effects for r,rr)
  // 8 = differ by stock, age,year (n_stocks fixed effects, n_years x n_ages random effects for r,rr)
  // 9 = differ by season (n_seasons fixed effects for r,rr). 
  // 10 = differ by season,age (n_seasons fixed effects, n_ages random effects for r,rr). 
  // 11 = differ by season,year (n_seasons fixed effects, n_years random effects for r,rr)
  // 12 = differ by season,age,year (n_seasons fixed effects, n_years x n_ages random effects for r,rr)
  // 13 = differ by stock, season (n_stocks x n_seasons fixed effects for r,rr). 
  // 14 = differ by stock, season, age (n_stocks x n_seasons fixed effects, n_ages random effects for r,rr). 
  // 15 = differ by stock, season, year (n_stocks x n_seasons fixed effects, n_years random effects for r,rr)
  // 16 = differ by stock, season, age,year (n_stocks x n_seasons fixed effects, n_years x n_ages random effects for r,rr)
  //priors, RE, get full lm link for migration parameters
  //n_stocks x n_ages x n_seasons x n_years_pop X n_regions x n_regions-1
  //continues random processes in projection years!
  if(n_regions>1){//only do this mess if the number of regions is greater than 1.
    if(use_mu_prior.sum()>0) {
      array<Type> nll_mu_prior = get_nll_mu_prior(mu_prior_re, trans_mu, trans_mu_prior_sigma, use_mu_prior, mu_model);
      nll += nll_mu_prior.sum();
      REPORT(nll_mu_prior);
      SIMULATE if(do_simulate_mu_prior_re){
        mu_prior_re = simulate_mu_prior_re(mu_prior_re, trans_mu, trans_mu_prior_sigma, use_mu_prior, mu_model);
        REPORT(mu_prior_re);
      }
    }
    //if((mu_model != 1) & (mu_model != 5) & (mu_model != 9) & (mu_model != 13)){ //some type of random effects
    //see(11.2);
    array<Type> nll_mu_re = get_nll_mu(mu_repars, mu_re, mu_model, can_move, years_use);
    nll += nll_mu_re.sum();
    REPORT(nll_mu_re);
    SIMULATE if(do_simulate_mu_re){
      mu_re = simulate_mu_re(mu_repars, mu_re, mu_model, can_move, years_use);
    }
    if(do_post_samp_mu) ADREPORT(mu_re);
  }
    
  array<Type> trans_mu_base = get_trans_mu_base(trans_mu, mu_re, mu_prior_re, use_mu_prior, mu_model, Ecov_lm_mu, Ecov_how_mu);
  REPORT(trans_mu_base);
  //n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions - 1
  //rows sum to 1 for mig_type = 0 (prob move), rows sum to 0 for mig_type 1 (instantaneous)
  array<Type> mu = get_mu(trans_mu_base, can_move, must_move, mig_type, n_years_proj, n_years_model, proj_mu_opt, avg_years_ind);
  REPORT(mu);
  /////////////////////////////////////////

  /////////////////////////////////////////
  // Construct fishing mortality-at-age (FAA)
  matrix<Type> log_F = get_log_F(F_pars, F_config, n_years_pop);
  //n_fleets x n_years_pop x n_ages  (projection years not yet populated)
  array<Type> FAA = get_FAA(log_F, selAA, selblock_pointer_fleets, n_ages, n_years_model);
  /////////////////////////////////////////

  /////////////////////////////////////////
  //Population model and likelihoods
  //First: get everything needed to generate expected numbers at age
  
  //int P_dim = n_regions + n_fleets + 1; // probablity transition matrix is P_dim x P_dim
  
  //get probability transition matrices for yearly survival, movement, capture...
  array<Type> annual_Ps = get_annual_Ps(n_years_model, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons, FAA, log_M, mu, L);
  //seasonal PTMs for last year, just for inspection
  array<Type> seasonal_Ps_terminal_year = get_seasonal_Ps_y(n_years_model-1,fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons, 
    FAA, log_M, mu, L);
  REPORT(seasonal_Ps_terminal_year);
  //just survival categories for spawning
  array<Type> annual_SAA_spawn = get_annual_SAA_spawn(n_years_model, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons, fracyr_SSB, 
    spawn_seasons, FAA, log_M, mu, L); 

  //get annual stock-recruit pars if needed
  matrix<Type> log_SR_a = get_SR_log_a(recruit_model, mean_rec_pars, Ecov_lm_R, Ecov_how_R);
  matrix<Type> log_SR_b = get_SR_log_b(recruit_model, mean_rec_pars, Ecov_lm_R, Ecov_how_R);

  bool any_N1_re = false;
  for(int s = 0; s < n_stocks; s++) if(N1_model(s) ==2) any_N1_re = true;
  if(any_N1_re) { //Initial numbers at age are random effects
    matrix<Type> nll_N1 = get_nll_N1(N1_model, log_N1, N1_repars, NAA_where);
    nll += nll_N1.sum();
    //see(nll);
    REPORT(nll_N1);
    SIMULATE if(do_simulate_N_re){
      log_N1 = simulate_log_N1(N1_model, log_N1, N1_repars, NAA_where);
      REPORT(log_N1);
    }
    if(do_post_samp_N) ADREPORT(log_N1);
  }

  //initial realized numbers at age
  //n_stocks x n_regions x n_ages
  array<Type> N1 = get_NAA_1(N1_model,log_N1, NAA_where, log_M, FAA, which_F_age, 
   spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, mu, L, fracyr_seasons, 
   avg_years_ind, n_regions_is_small);
  REPORT(N1);
  // vector<array<Type>> N1_components = get_eq_NAA_components(N1_model,log_N1, NAA_where, log_M, FAA, which_F_age, spawn_regions, fleet_regions, fleet_seasons, 
  //     can_move, mig_type, mu, L, fracyr_seasons, avg_years_ind, n_regions_is_small);
  // // see(N1_components.size());
  // array<Type> FAA11 = N1_components(0);
  // // see(FAA11);
  // REPORT(FAA11);
  // array<Type> SAA11 = N1_components(1);
  // REPORT(SAA11);
  // if(n_stocks == 2){
  //   array<Type> FAA12 = N1_components(2);
  //   // see(FAA12);
  //   REPORT(FAA12);
  //   array<Type> SAA12 = N1_components(3);
  //   REPORT(SAA12);
  // }
  // array<Type> selex = N1_components(n_stocks*2);
  // REPORT(selex);

  //initial predicted numbers at age
  //n_stocks x n_regions x n_ages
  array<Type> pred_N1 = get_pred_N1(N1_model, N1, NAA_where, N1_repars);
  REPORT(pred_N1);

  //should work for SCAA and RE models
  array<Type> all_NAA = get_all_NAA(NAA_re_model, N1_model, N1, N1_repars, log_NAA, NAA_where, 
   mature_all, waa_ssb, recruit_model, mean_rec_pars, log_SR_a, log_SR_b, 
   Ecov_how_R, Ecov_lm_R, spawn_regions,  annual_Ps, annual_SAA_spawn, n_years_model,0); //log_NAA should be mapped accordingly to exclude NAA=0 e.g., recruitment by region.
  array<Type> NAA = extract_NAA(all_NAA);
  //This will use get_all_NAA, get_SSB, and get_pred_NAA to form devs and calculate likelihoods
  matrix<Type> R_XSPR = get_RXSPR(all_NAA, spawn_regions, n_years_model, n_years_proj, XSPR_R_opt, XSPR_R_avg_yrs);
  

  //need to do projections before evaluating nll component for NAA
  vector<Type> fracyr_ssb_y = get_avg_ssbfrac(fracyr_SSB,avg_years_ind); 
  matrix<Type> mat_y = get_avg_mat(mature,avg_years_ind);
  matrix<Type> waa_ssb_y = get_avg_waa(waa,avg_years_ind,waa_pointer_ssb);
  matrix<Type> waa_catch_y = get_avg_waa(waa,avg_years_ind,waa_pointer_fleets);
  if(n_years_proj > 0){

    for(int y = n_years_model; y < n_years_pop; y++){
      fracyr_SSB_all.row(y) = fracyr_ssb_y;
      for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) {
        mature_all(s,y,a) = mat_y(s,a);
        waa_ssb(s,y,a) = waa_ssb_y(s,a);
      }
      for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) waa_catch(f,y,a) = waa_catch_y(f,a);

      // see("yproj");
      // see(y);
      // see(annual_Ps.dim);
      all_NAA = update_all_NAA(y, all_NAA, NAA_re_model, N1_model, N1, N1_repars, log_NAA, NAA_where, 
        mature_all, waa_ssb, recruit_model, mean_rec_pars, log_SR_a, log_SR_b, 
        Ecov_how_R, Ecov_lm_R, spawn_regions,  annual_Ps, annual_SAA_spawn, n_years_model, logR_proj, no_trace);
      NAA = extract_NAA(all_NAA);
      R_XSPR = get_RXSPR(all_NAA, spawn_regions, n_years_model, n_years_proj, XSPR_R_opt, XSPR_R_avg_yrs);
      //There are many options for defining F in projection years so a lot of inputs
      FAA = update_FAA_proj(y, proj_F_opt, FAA, NAA, log_M, mu, L, mat_y, waa_ssb_y, waa_catch_y, fleet_regions, fleet_seasons, 
        fracyr_ssb_y, spawn_regions, can_move, must_move, mig_type, avg_years_ind, n_years_model, which_F_age, fracyr_seasons, 
            n_regions_is_small, percentSPR, proj_Fcatch, percentFXSPR, percentFMSY, R_XSPR,
        FXSPR_init, FMSY_init, F_proj_init, log_SR_a, log_SR_b, spawn_seasons, recruit_model, SPR_weights, SPR_weight_type, no_trace);
        // if(trace) see(y);
        // if(trace) for(int a = 0; a < n_ages; a++) see(FAA(0,y,a));
      annual_Ps = update_annual_Ps(y, annual_Ps, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons, FAA, log_M, mu, L);
      annual_SAA_spawn = update_annual_SAA_spawn(y, annual_SAA_spawn, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons, 
        fracyr_SSB_all, spawn_seasons, FAA, log_M, mu, L);
    }
    //if(trace) std::exit(EXIT_FAILURE);
  }

  NAA = extract_NAA(all_NAA);
  array<Type> pred_NAA = extract_pred_NAA(all_NAA);
  array<Type> NAA_devs = get_NAA_devs(all_NAA, NAA_where, NAA_re_model);

  if(NAA_re_model.sum() > 0){ //at least one stock is not SCAA
    matrix<Type> nll_NAA = get_NAA_nll(NAA_re_model, all_NAA, log_NAA_sigma, trans_NAA_rho, NAA_where, spawn_regions, years_use, bias_correct_pe,
      use_alt_AR1);
    // matrix<Type> nll_NAA = get_NAA_nll(N1, N1_model, N1_repars, NAA_re_model, log_NAA, log_NAA_sigma, trans_NAA_rho, NAA_where, recruit_model, mean_rec_pars, log_SR_a, log_SR_b, 
    //   Ecov_how_R, Ecov_lm_R, spawn_regions, annual_Ps, annual_SAA_spawn, waa_ssb, mature, years_use, bias_correct_pe, n_years_model);
    nll += nll_NAA.sum();
    //see(nll);
    REPORT(nll_NAA);
    SIMULATE if(do_simulate_N_re){
      array<Type> NAA_devs_sim = simulate_NAA_devs(NAA, NAA_re_model, log_NAA_sigma, trans_NAA_rho, NAA_where, spawn_regions, years_use, 
        bias_correct_pe);
      //repopulate log_NAA, NAA, pred_NAA, SSB,etc.
      log_NAA = get_simulated_log_NAA(N1_model, N1, N1_repars, NAA_re_model, NAA_devs_sim, log_NAA, NAA_where, recruit_model, mean_rec_pars,
        log_SR_a, log_SR_b, Ecov_how_R, Ecov_lm_R, spawn_regions, annual_Ps, annual_SAA_spawn, waa_ssb, mature_all, n_years_model, logR_proj);
      all_NAA = get_all_NAA(NAA_re_model, N1_model, N1, N1_repars, log_NAA, NAA_where, 
        mature_all, waa_ssb, recruit_model, mean_rec_pars, log_SR_a, log_SR_b, 
        Ecov_how_R, Ecov_lm_R, spawn_regions,  annual_Ps, annual_SAA_spawn, n_years_model,no_trace);
      R_XSPR = get_RXSPR(all_NAA, spawn_regions, n_years_model, n_years_proj, XSPR_R_opt, XSPR_R_avg_yrs);
      if(n_years_proj > 0){

        for(int y = n_years_model; y < n_years_pop; y++){
          log_NAA = get_simulated_log_NAA(N1_model, N1, N1_repars, NAA_re_model, NAA_devs_sim, log_NAA, NAA_where, recruit_model, mean_rec_pars,
            log_SR_a, log_SR_b, Ecov_how_R, Ecov_lm_R, spawn_regions, annual_Ps, annual_SAA_spawn, waa_ssb, mature_all, n_years_model, logR_proj);
          all_NAA = update_all_NAA(y, all_NAA, NAA_re_model, N1_model, N1, N1_repars, log_NAA, NAA_where, 
            mature_all, waa_ssb, recruit_model, mean_rec_pars, log_SR_a, log_SR_b, 
            Ecov_how_R, Ecov_lm_R, spawn_regions,  annual_Ps, annual_SAA_spawn, n_years_model, logR_proj, no_trace);
            
          R_XSPR = get_RXSPR(all_NAA, spawn_regions, n_years_model, n_years_proj, XSPR_R_opt, XSPR_R_avg_yrs);
          NAA = extract_NAA(all_NAA);
          //There are many options for defining F in projection years so a lot of inputs
          FAA = update_FAA_proj(y, proj_F_opt, FAA, NAA, log_M, mu, L, mat_y, waa_ssb_y, waa_catch_y, fleet_regions, fleet_seasons, 
            fracyr_ssb_y, spawn_regions, can_move, must_move, mig_type, avg_years_ind, n_years_model, which_F_age, fracyr_seasons, 
            n_regions_is_small, percentSPR, proj_Fcatch, percentFXSPR, percentFMSY, R_XSPR, FXSPR_init, FMSY_init, F_proj_init, 
            log_SR_a, log_SR_b, spawn_seasons, recruit_model, SPR_weights, SPR_weight_type, no_trace);
          annual_Ps = update_annual_Ps(y, annual_Ps, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons, FAA, log_M, mu, L);
          annual_SAA_spawn = update_annual_SAA_spawn(y, annual_SAA_spawn, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons, 
            fracyr_SSB_all, spawn_seasons, FAA, log_M, mu, L);
        }

      }

      NAA = extract_NAA(all_NAA);
      pred_NAA = extract_pred_NAA(all_NAA);
      NAA_devs = get_NAA_devs(all_NAA, NAA_where, NAA_re_model);
      REPORT(log_NAA);
      REPORT(NAA_devs_sim);
    }
    if(do_post_samp_N) ADREPORT(log_NAA);
  } else { // SCAA models
    // likelihood of NAA deviations
    if(n_years_proj > 0){ // SCAA treats recruitment in proj years as random effects with fixed mean, SD
      matrix<Type> nll_Rproj(n_years_proj, n_stocks);
      nll_Rproj.setZero();
      for(int y = 0; y < n_years_proj; y++) for(int s = 0; s < n_stocks; s++) if(NAA_re_model(s) == 0){
        nll_Rproj(y,s) -= dnorm(logR_proj(y,s), logR_mean(s), logR_sd(s), 1);
        SIMULATE if(do_simulate_N_re) logR_proj(y,s) = rnorm(logR_mean(s), logR_sd(s));
      }
      REPORT(nll_Rproj);
      nll += nll_Rproj.sum();
      SIMULATE if(do_simulate_N_re) REPORT(logR_proj);
    }
  }
  //need to do this
  //log_F = update_log_F(log_F, FAA, which_F_age);
  //matrix<Type> F(n_years_pop,n_fleets); //n_years_pop x n_fleets (projection years not yet populated)
  //for(int f = 0; f < log_F.cols(); f++) F.col(f) = exp(vector<Type> (log_F.col(f)));
  // see(F);
  REPORT(annual_Ps);
  REPORT(annual_SAA_spawn);
  //REPORT(F);
  REPORT(all_NAA);
  // REPORT(log_F);
  REPORT(NAA);
  REPORT(pred_NAA);
  REPORT(NAA_devs);
  REPORT(FAA);
  REPORT(R_XSPR);
  REPORT(mature_all);
  REPORT(waa_ssb);
  REPORT(waa_catch);
  REPORT(fracyr_SSB_all);

  //Now get annual NAA at spawning and SSB.
  array<Type> NAA_spawn = get_NAA_spawn(NAA, annual_SAA_spawn, spawn_regions);
  REPORT(NAA_spawn);
  matrix<Type> SSB = get_SSB(NAA_spawn,waa_ssb,mature_all);
  REPORT(SSB);
  matrix<Type> log_SSB(n_years_pop,n_stocks);
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years_pop; y++) log_SSB(y,s) = log(SSB(y,s));
  //need to be careful here about log_NAA random effects in regions where there will be zero predicted due to migration parameterization
  
  /////////////////////////////////////////
  //catch observations
  array<Type> pred_stock_CAA = get_pred_stock_CAA(NAA, annual_Ps);
  REPORT(pred_stock_CAA);
  array<Type> pred_CAA = get_pred_CAA(pred_stock_CAA);
  REPORT(pred_CAA);
  array<Type> pred_catch_paa = get_pred_catch_paa(pred_CAA, n_years_model);
  REPORT(pred_catch_paa);
  array<Type> pred_stock_catch = get_pred_stock_catch(pred_stock_CAA,waa_catch);
  REPORT(pred_stock_catch);
  matrix<Type> pred_catch = get_pred_catch(pred_stock_catch);
  REPORT(pred_catch);
  matrix<Type> pred_log_catch = get_pred_log_catch(pred_catch, agg_catch_sigma, log_catch_sig_scale, bias_correct_oe);
  REPORT(pred_log_catch);
  
  matrix<Type> nll_agg_catch = get_nll_agg_catch(pred_log_catch, agg_catch_sigma, log_catch_sig_scale, obsvec,
    use_agg_catch, keep_C, keep);
  nll += nll_agg_catch.sum();
  //see(nll);
  REPORT(nll_agg_catch);
  SIMULATE if(do_simulate_data(0)){
    agg_catch = simulate_agg_catch(pred_log_catch, agg_catch, agg_catch_sigma, log_catch_sig_scale, use_agg_catch);
    REPORT(agg_catch);
    obsvec = sim_agg_catch_in_obsvec(obsvec,keep_C,agg_catch, use_agg_catch);
  }

  matrix<Type> nll_catch_acomp = get_nll_catch_acomp(pred_catch_paa, use_catch_paa, catch_paa,
    catch_Neff, age_comp_model_fleets, catch_paa_pars, keep_Cpaa, keep, obsvec, agesvec, do_osa);
  nll += nll_catch_acomp.sum();
  REPORT(nll_catch_acomp);
  matrix<Type> catch_Neff_out = get_Neff_out(catch_Neff, age_comp_model_fleets, catch_paa_pars);
  REPORT(catch_Neff_out);
  //see(nll);
  SIMULATE if(do_simulate_data(0)){
    obsvec = simulate_catch_paa_in_obsvec(obsvec, agesvec, pred_catch_paa, use_catch_paa,  keep_Cpaa, catch_Neff, 
      age_comp_model_fleets, catch_paa_pars, no_trace);
    catch_paa = sim_obsvec_in_catch_paa(obsvec, agesvec, catch_paa, use_catch_paa, keep_Cpaa, age_comp_model_fleets, no_trace);
    REPORT(catch_paa);
  }
  /////////////////////////////////////////

  
  /////////////////////////////////////////
  //index observations
  array<Type> NAA_index = get_NAA_index(NAA, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons, fracyr_indices, index_seasons,
    index_regions, FAA, log_M, mu, L, n_years_model);
  REPORT(NAA_index);
  array<Type> pred_IAA = get_pred_IAA(QAA, NAA_index);
  REPORT(pred_IAA);
  array<Type> pred_index_paa = get_pred_index_paa(pred_IAA, units_index_paa, waa, waa_pointer_indices);
  REPORT(pred_index_paa);

  matrix<Type> pred_indices = get_pred_indices(pred_IAA, units_indices, waa, waa_pointer_indices);
  REPORT(pred_indices);
  matrix<Type> pred_log_indices = get_pred_log_indices(pred_indices, agg_index_sigma, log_index_sig_scale, bias_correct_oe);
  REPORT(pred_log_indices);

  matrix<Type> nll_agg_indices = get_nll_agg_indices(pred_log_indices, agg_index_sigma, log_index_sig_scale, obsvec,
    use_indices, keep_I, keep);
  nll += nll_agg_indices.sum();
  //see(nll);
  REPORT(nll_agg_indices);
  SIMULATE if(do_simulate_data(1)){
    agg_indices = simulate_agg_indices(pred_log_indices, agg_indices, agg_index_sigma, log_index_sig_scale, use_indices);
    REPORT(agg_indices);
    obsvec = sim_agg_indices_in_obsvec(obsvec,keep_I,agg_indices, use_indices);
  }

  matrix<Type> nll_index_acomp = get_nll_index_acomp(pred_index_paa, use_index_paa, index_paa,
    index_Neff, age_comp_model_indices, index_paa_pars, keep_Ipaa, keep, obsvec, agesvec, do_osa);
  nll += nll_index_acomp.sum();
  REPORT(nll_index_acomp);
  matrix<Type> index_Neff_out = get_Neff_out(index_Neff, age_comp_model_indices, index_paa_pars);
  REPORT(index_Neff_out);
  //see(nll);
  SIMULATE if(do_simulate_data(1)){
    obsvec = simulate_index_paa_in_obsvec(obsvec, agesvec, pred_index_paa, use_index_paa,  keep_Ipaa, index_Neff, 
      age_comp_model_indices, index_paa_pars);
    index_paa = sim_obsvec_in_index_paa(obsvec, agesvec, index_paa, use_index_paa, keep_Ipaa, age_comp_model_indices);
    REPORT(index_paa);
    // vector<Type> tf_paa_check = obsvec.segment(keep_Ipaa(0,5,0),keep_Ipaa(0,5,1));
    // REPORT(tf_paa_check);
    // vector<Type> t_pred_paa(n_ages);
    // for(int a = 0; a < n_ages; a++) t_pred_paa(a) = pred_index_paa(0,5,a);    
    // vector<int> age_check = agesvec.segment(keep_Ipaa(0,5,0), keep_Ipaa(0,5,1));
    // REPORT(age_check);
    // vector<Type> tf_paa_check = sim_acomp(t_pred_paa, index_Neff(5,0), age_check, age_comp_model_indices(0), 
    //   vector<Type>(index_paa_pars.row(0)));
    // REPORT(tf_paa_check);
    // vector<Type> paa_check = make_paa(tf_paa_check, age_comp_model_indices(0), age_check, n_ages);
    // //vector<Type> paa_check = obsvec_to_paa(0, 5, obsvec, agesvec, use_index_paa, keep_Ipaa, age_comp_model_indices, n_ages);
    // REPORT(paa_check);
  }
  /////////////////////////////////////////
  SIMULATE if(sum(do_simulate_data) > 0) REPORT(obsvec);
      //see(log_M);
  REPORT(nll);


  if(do_SPR_BRPs){
    vector< array<Type>> annual_SPR_res = get_annual_SPR_res(SPR_weights, log_M, FAA, spawn_seasons,  
      spawn_regions, fleet_regions, fleet_seasons, fracyr_seasons, can_move, must_move, mig_type, trans_mu_base, 
      L, which_F_age, waa_ssb, waa_catch, mature_all, percentSPR, NAA, fracyr_SSB_all, FXSPR_init, 
      R_XSPR, n_regions_is_small, SPR_weight_type, 0, 10);
    
    array<Type> log_FAA_XSPR = annual_SPR_res(0);
    REPORT(log_FAA_XSPR);
    array<Type> log_SSB_FXSPR = annual_SPR_res(1);
    REPORT(log_SSB_FXSPR);
    array<Type> log_Y_FXSPR = annual_SPR_res(2);
    REPORT(log_Y_FXSPR);
    array<Type> log_SPR_FXSPR = annual_SPR_res(3);
    REPORT(log_SPR_FXSPR);
    array<Type> log_SPR0 = annual_SPR_res(4);
    REPORT(log_SPR0);
    array<Type> log_YPR_FXSPR = annual_SPR_res(5);
    REPORT(log_YPR_FXSPR);
    array<Type> log_FXSPR_iter = annual_SPR_res(6); 
    REPORT(log_FXSPR_iter);
    vector<Type> log_FXSPR = log_FXSPR_iter.matrix().col(9);
    REPORT(log_FXSPR);

    vector< array<Type>> static_SPR_res =  get_SPR_res(SPR_weights, log_M, FAA, spawn_seasons,  
      spawn_regions, fleet_regions, fleet_seasons, fracyr_seasons, can_move, must_move, mig_type, trans_mu_base, 
      L, which_F_age_static, waa_ssb, waa_catch, mature_all, percentSPR, NAA, fracyr_SSB_all, FXSPR_static_init, 
      avg_years_ind, avg_years_ind, avg_years_ind, avg_years_ind, avg_years_ind, avg_years_ind, avg_years_ind, 
      vector<Type> (R_XSPR.row(n_years_model-1)), //This will be constant across years if XSPR_R_opt = 2 or 4
      n_regions_is_small, SPR_weight_type, 0, 10);
    
    //matrix<Type> log_FAA_XSPR_static = static_SPR_res(0).matrix();
    array<Type> log_FAA_XSPR_static = static_SPR_res(0); //(n_fleets + n_regions + 1) x n_ages
    // matrix<Type> log_FAA_XSPR_static(n_fleets+1,n_ages);
    // for(int f = 0; f <= n_fleets; f++) for(int a = 0; a < n_ages; a++) log_FAA_XSPR_static(f,a) = static_SPR_res(0)(n_ages*f + a);
    vector<Type> log_SSB_FXSPR_static = static_SPR_res(1).matrix().row(0);
    //array<Type> log_SSB_FXSPR_static = static_SPR_res(1);
    vector<Type> log_Y_FXSPR_static = static_SPR_res(2).matrix().row(0);
    //array<Type> log_Y_FXSPR_static = static_SPR_res(2);
    vector<Type> log_SPR_FXSPR_static = static_SPR_res(3).matrix().row(0);
    //array<Type> log_SPR_FXSPR_static = static_SPR_res(3);
    vector<Type> log_SPR0_static = static_SPR_res(4).matrix().row(0);
    //array<Type> log_SPR0_static = static_SPR_res(4);
    //matrix<Type> log_YPR_FXSPR_static = static_SPR_res(5).matrix();
    array<Type> log_YPR_FXSPR_static = static_SPR_res(5);
    // matrix<Type> log_YPR_FXSPR_static(n_stocks,n_fleets+1);
    // for(int s = 0; s < n_stocks; s++) {
    //   for(int f = 0; f < n_fleets; f++) log_YPR_FXSPR_static(s,f) = static_SPR_res(5)(s*n_fleets + f);
    //   log_YPR_FXSPR_static(s,n_fleets) = static_SPR_res(5)(n_stocks*n_fleets + s);
    // }
    //matrix<Type> log_FXSPR_iter_static = static_SPR_res(6).matrix();
    //array<Type> log_FXSPR_iter_static = static_SPR_res(6);
    vector<Type> log_FXSPR_iter_static = static_SPR_res(6).matrix().row(0);
    Type log_FXSPR_static = log_FXSPR_iter_static(log_FXSPR_iter_static.size()-1);
    REPORT(log_FAA_XSPR_static);
    REPORT(log_SSB_FXSPR_static);
    REPORT(log_Y_FXSPR_static);
    REPORT(log_SPR_FXSPR_static);
    REPORT(log_SPR0_static);
    REPORT(log_YPR_FXSPR_static);
    REPORT(log_FXSPR_static);
    REPORT(log_FXSPR_iter_static);

    array<Type> annual_SPR0AA = get_annual_SPR0_at_age(log_M, spawn_seasons, fracyr_seasons, can_move, must_move,
      mig_type, trans_mu_base, L, waa_ssb,  mature_all, fracyr_SSB_all, n_regions_is_small);
    REPORT(annual_SPR0AA);

    if((sum_do_post_samp == 0) & (mig_type.sum() == 0)) {
      ADREPORT(log_FXSPR);
      ADREPORT(log_SSB_FXSPR);
      ADREPORT(log_Y_FXSPR);
      ADREPORT(log_SPR0);
      ADREPORT(log_FAA_XSPR_static);
      ADREPORT(log_FXSPR_static);
      ADREPORT(log_SSB_FXSPR_static);
      ADREPORT(log_SPR0_static);
      ADREPORT(log_Y_FXSPR_static);
    }
  }
  int is_SR = 0;
  for(int s = 0; s < n_stocks; s++) if((recruit_model(s) == 3) | (recruit_model(s) == 4)) is_SR++;
  if((is_SR> 0) & do_MSY_BRPs) {

    vector< array <Type> > annual_MSY_res = get_annual_MSY_res(recruit_model,
      log_SR_a, log_SR_b, log_M, FAA, spawn_seasons, spawn_regions, fleet_regions,
      fleet_seasons, fracyr_seasons, can_move, must_move, mig_type, trans_mu_base, 
      L, which_F_age, waa, waa_pointer_ssb, waa_pointer_fleets, mature, fracyr_SSB, FMSY_init, 
      n_regions_is_small, 1, 10);
    
    array<Type> log_SSB_MSY = annual_MSY_res(0);
    REPORT(log_SSB_MSY);
    array<Type> log_R_MSY = annual_MSY_res(1);
    REPORT(log_R_MSY);
    array<Type> log_SPR_MSY = annual_MSY_res(2);
    REPORT(log_SPR_MSY);
    array<Type> log_FAA_MSY = annual_MSY_res(3);
    REPORT(log_FAA_MSY);
    array<Type> log_MSY = annual_MSY_res(4);
    REPORT(log_MSY);
    array<Type> log_YPR_MSY = annual_MSY_res(5);
    REPORT(log_YPR_MSY);
    array<Type> log_FMSY_iter = annual_MSY_res(6);
    REPORT(log_FMSY_iter);
    vector<Type> log_FMSY = log_FMSY_iter.matrix().col(9);
    REPORT(log_FMSY);

    vector<Type> log_FMSY_alt = get_log_FMSY(FAA, fleet_regions, fleet_seasons, spawn_seasons, spawn_regions, can_move, mig_type, 
      fracyr_seasons, which_F_age, recruit_model, log_SR_a, log_SR_b, fracyr_SSB_all, log_M, mu, L, waa_ssb, waa_catch, mature_all, n_regions_is_small,
      FMSY_init,0);

    vector< matrix<Type>> static_MSY_res =  get_MSY_res(recruit_model,
      log_SR_a, log_SR_b, log_M, FAA, spawn_seasons, spawn_regions, fleet_regions,
      fleet_seasons, fracyr_seasons, can_move, must_move, mig_type, trans_mu_base, 
      L, which_F_age_static, waa, waa_pointer_ssb, waa_pointer_fleets, mature, fracyr_SSB, FMSY_static_init, 
      avg_years_ind, avg_years_ind, avg_years_ind, avg_years_ind, avg_years_ind, avg_years_ind, avg_years_ind, avg_years_ind,
      n_regions_is_small, 0, 10);
      vector<Type> log_SSB_MSY_static = static_MSY_res(0).col(0);
      REPORT(log_SSB_MSY_static);
      vector<Type> log_R_MSY_static = static_MSY_res(1).col(0);
      REPORT(log_R_MSY_static);
      vector<Type> log_SPR_MSY_static = static_MSY_res(2).col(0);
      REPORT(log_SPR_MSY_static);
      matrix<Type> log_FAA_MSY_static = static_MSY_res(3);
      REPORT(log_FAA_MSY_static);
      matrix<Type> log_MSY_static = static_MSY_res(4);
      REPORT(log_MSY_static);
      matrix<Type> log_YPR_MSY_static = static_MSY_res(5);
      REPORT(log_YPR_MSY_static);
      matrix<Type> log_FMSY_iter_static = static_MSY_res(6);
      REPORT(log_FMSY_iter_static);
      vector<Type> log_FMSY_static = log_FMSY_iter_static.col(9);
      REPORT(log_FMSY_static);

    if(sum_do_post_samp == 0) if((n_regions == 1) | (mig_type.sum() == 0)) {
      ADREPORT(log_FMSY);
      ADREPORT(log_SSB_MSY);
      ADREPORT(log_R_MSY);
      ADREPORT(log_MSY);
      ADREPORT(log_FMSY_static);
      ADREPORT(log_SSB_MSY_static);
      ADREPORT(log_R_MSY_static);
      ADREPORT(log_MSY_static);
    }

    REPORT(log_SR_a);
    REPORT(log_SR_b);
    if(sum_do_post_samp == 0){
      ADREPORT(log_SR_a);
      ADREPORT(log_SR_b);
    }  
  }
  matrix<Type> log_index_resid(n_years_model, n_indices), log_catch_resid(n_years_model, n_fleets);
  log_index_resid.setZero();
  log_catch_resid.setZero();
  for(int y = 0; y < n_years_model; y++){
    for(int i = 0; i < n_indices; i++){
      if(use_indices(y,i) == 1) log_index_resid(y,i) = log(agg_indices(y,i)) - pred_log_indices(y,i);
    }
    for(int f = 0; f < n_fleets; f++) log_catch_resid(y,f) = log(agg_catch(y,f)) - pred_log_catch(y,f);
  }
  REPORT(log_catch_resid);
  REPORT(log_index_resid);
  
  //if(reportMode==0){
  array<Type> log_FAA = get_log_FAA(FAA);
  array<Type> FAA_by_region = get_FAA_by_region(FAA, fleet_regions, n_regions);
  REPORT(FAA_by_region);
  array<Type> log_FAA_by_region = get_log_FAA(FAA_by_region); //FAA_by_region is also a 3-d array (region,year,age)
  array<Type> log_NAA_rep = get_log_NAA_rep(NAA, NAA_where);
  matrix<Type> Fbar(FAA_by_region.dim(1),FAA_by_region.dim(0));
  Fbar.setZero();
  int n_Fbar_ages = Fbar_ages.size();
  matrix<Type> log_FAA_tot(FAA_by_region.dim(1),FAA_by_region.dim(2));
  vector<Type> log_F_tot(FAA_by_region.dim(1));
  vector<Type> log_SSB_all(SSB.rows());
  log_SSB_all.setZero();
  log_FAA_tot.setZero(); log_F_tot.setZero();
  for(int y = 0; y < log_SSB_all.size(); y++) log_SSB_all(y) = log(SSB.row(y).sum());
  for(int y = 0; y < FAA_by_region.dim(1); y++) {
    for(int r = 0; r < FAA_by_region.dim(0); r++) for(int a = 0; a < n_Fbar_ages; a++) {
      Fbar(y,r) += FAA_by_region(r,y,Fbar_ages(a)-1)/Type(n_Fbar_ages);
    }
    for(int a = 0; a < FAA_by_region.dim(2); a++){
      for(int r = 0; r < FAA_by_region.dim(0); r++) log_FAA_tot(y,a) += FAA_by_region(r,y,a);
      log_FAA_tot(y,a) = log(log_FAA_tot(y,a));
    }
    log_F_tot(y) = log_FAA_tot(y,which_F_age(y)-1);
  }
  matrix<Type> log_Fbar = log(Fbar.array());
  REPORT(Fbar);

  REPORT(q);
  REPORT(QAA);
  REPORT(log_F_tot);
  REPORT(log_FAA_tot);
  REPORT(log_FAA_by_region);
  //}
  if(sum_do_post_samp == 0){
    ADREPORT(log_NAA_rep);
    ADREPORT(log_SSB);
    ADREPORT(log_SSB_all);
    //ADREPORT(log_F);
    ADREPORT(log_FAA);
    ADREPORT(log_FAA_tot);
    ADREPORT(log_FAA_by_region);
    ADREPORT(log_F_tot);
    ADREPORT(log_Fbar);
    //ADREPORT(log_catch_resid);
    //ADREPORT(log_index_resid);
    if((n_regions>1) & (sum(can_move)>0)){ //only adreport the necessary parameters
      array<int> mu_sdrep_index = get_mu_sdrep_indices(mu_model, trans_mu_base);
      REPORT(mu_sdrep_index);
      vector<Type> trans_mu_base_sdrep = get_trans_mu_base_sdrep(trans_mu_base, mu_model, mu_sdrep_index);
      REPORT(trans_mu_base_sdrep);
      ADREPORT(trans_mu_base_sdrep);
    }
  }
  //see(nll);
  return nll;
}

