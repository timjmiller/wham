0#define TMB_LIB_INIT R_init_wham
#include <TMB.hpp>
#include "all.hpp"


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  
  DATA_INTEGER(n_years_model); 
  DATA_INTEGER(n_seasons);
  DATA_VECTOR(fracyr_seasons); //length of intervals for seasons
  DATA_INTEGER(n_regions);
  DATA_INTEGER(n_stocks);
  DATA_INTEGER(n_fleets);
  DATA_INTEGER(n_indices);
  DATA_INTEGER(n_ages);
  DATA_IVECTOR(n_ages_fleet);
  DATA_IVECTOR(n_ages_indices);
  DATA_IVECTOR(mig_type); //n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
  
  DATA_MATRIX(fracyr_SSB); //n_years x n_stocks:  size of interval from beginning of season to time of spawning within that season
  DATA_IVECTOR(spawn_regions); //n_stocks
  DATA_IVECTOR(spawn_seasons); //n_stocks
  //DATA_IARRAY(n_mu) //n_stocks x n_years x x n_seasons x n_ages
  //DATA_IVECTOR(mu_row) //length the total number of migration parameters
  //DATA_IVECTOR(mu_col) //length the total number of migration parameters
  //DATA_IVECTOR(mu_pointer) //n_stocks * n_years * n_ages
  //DATA_INTEGER(n_seasons_recruited); //easiest if this is = n_seasons 
  DATA_ARRAY(mature); //n_stocks x n_years x n_ages
  DATA_IVECTOR(waa_pointer_fleets);
  DATA_IVECTOR(waa_pointer_totcatch); //n_regions
  DATA_IVECTOR(waa_pointer_indices);
  DATA_IVECTOR(waa_pointer_ssb); //n_stocks
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
  DATA_INTEGER(N1_model); //0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
  DATA_IMATRIX(NAA_where); //n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
  DATA_IVECTOR(n_M_a); //length = n_regions
  DATA_IMATRIX(n_M_re); // n_stocks x n_regions how many time-varying RE each year? n_ages? 1? n_est_M?
  DATA_IVECTOR(M_model); //length = n_regions; 1: "all ages all stocks", 2: "by age all stocks", 3: "weight-at-age all stocks", 4: "all ages by stock", 5: "by age by stock", 6: "weight-at-age by stock"
  DATA_IVECTOR(M_re_model); //length = n_regions; 1 = none, 2 = IID, 3 = ar1_a, 4 = ar1_y, 5 = 2dar1
  DATA_INTEGER(use_b_prior); //for M_model = 3: M = a*W^b model, use prior (and re) for log_b
  DATA_INTEGER(log_b_model); //1: constant, 2: differ by stock, 3: differ by region, 4: differ by both
  DATA_IVECTOR(L_model); //length = n_regions; L_model = 0: don't use L (extra mortality); L_model = 1: use constant L; L_model = 2: use iid re; L_model = 3: use ar1 re 
  DATA_IARRRAY(can_move); //n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
  DATA_IARRRAY(must_move); //n_stocks x ages x n_seasons x n_regions: 0/1 determining if it must leave the region
  DATA_IARRAY(use_mu_prior); //n_stocks x ages x n_seasons x n_regions x n_regions-1: 0/1 whether to apply prior for each movement parameter
  DATA_INTEGER(mu_model); 
  // 1 = constant across stocks, ages, time (n_seasons x n_regions x (n_regions -1) pars). 
  // 2 = differ by age (n_seasons x n_regions x (n_regions -1) fixed effects, n_ages AR1 random effects for each). 
  // 3 = differ by year (n_seasons x n_regions x (n_regions -1) fixed effects, n_years AR1 random effects for each)
  // 4 = differ by stock (n_seasons x n_stocks x n_regions x (n_regions -1) pars). 
  // 5 = differ by stock, age (n_seasons x n_stocks x n_regions x (n_regions -1) fixed effects, n_ages AR1 random effects for each). 
  // 6 = differ by stock, year (n_seasons x n_stocks x n_regions x (n_regions -1) fixed effects, n_years AR1 random effects for each)
  //DATA_IARRAY(use_mu_re); //n_stocks x ages x n_seasons x n_years_model x n_regions x n_regions-1: 0/1 whether to use temporal RE for each movement parameter
  
  DATA_IVECTOR(which_F_age); // (n_years_model + n_years_proj) x n_stocks; which age of F to use for full total F for msy/ypr calculations and projections
  DATA_INTEGER(use_steepness); // which parameterization to use for BH/Ricker S-R, if needed.
  DATA_INTEGER(bias_correct_pe); //bias correct lognormal process error?
  DATA_INTEGER(bias_correct_oe); //bias correct lognormal observation error?
  DATA_IVECTOR(Fbar_ages);
  
  //DATA_IMATRIX(R1_pointer); //n_stocks x n_regions. Tells where log_R1 parameters are used.
  //DATA_IMATRIX(N1_sigma_pointer); //n_stocks x n_regions. Tells which N1_sigma_par to use where.
  //DATA_IARRAY(NAA_sigma_pointers); //n_stocks x n_ages x n_regions
  //DATA_IARRAY(NAA_re_indicator); //n_stocks x (n_years-1) x n_ages x n_regions will estimate a random effect where indicator is not 0. sum = n_NAA_re.

  //DATA_IVECTOR(simulate_state); //vector (0/1) if 1 then state parameters (NAA, MAA, sel, Ecov, q) in that order) will be simulated.
  DATA_INTEGER(do_simulate_NAA_re); //(0/1) if 1 then simulate numbers at age random effects.
  DATA_INTEGER(do_simulate_Ecov_re); //(0/1) if 1 then simulate selectivity random effects.
  DATA_INTEGER(do_simulate_sel_re); //(0/1) if 1 then simulate selectivity random effects.
  DATA_INTEGER(do_simulate_M_re); //(0/1) if 1 then simulate natural mortality random effects.
  DATA_INTEGER(do_simulate_q_re); //(0/1) if 1 then simulate catchability random effects.
  DATA_INTEGER(do_simulate_q_prior_re); //(0/1) if 1 then simulate q prior random effects.
  DATA_INTEGER(do_simulate_mu_re); //(0/1) if 1 then simulate mu (migration) random effects.
  DATA_INTEGER(do_simulate_L_re); //(0/1) if 1 then simulate L (extra mortality) random effects.
  DATA_INTEGER(do_simulate_data); //vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
  //DATA_IVECTOR(do_simulate_period); //vector (0/1) if 1 then period (model years, projection years) will be simulated.
  DATA_SCALAR(percentSPR); // percentage to use for SPR-based reference points. Default = 40.
  DATA_SCALAR(percentFXSPR); // percent of F_XSPR to use for calculating catch in projections. For example, GOM cod uses F = 75% F_40%SPR, so percentFXSPR = 75 and percentSPR = 40. Default = 100.
  DATA_SCALAR(percentFMSY); // percent of FMSY to use for calculating catch in projections.
  DATA_INTEGER(XSPR_R_opt); //1(3): use annual R estimates(predictions) for annual SSB_XSPR, 2(4): use average R estimates(predictions). See next line for years to average over.
  DATA_IVECTOR(XSPR_R_avg_yrs); // model year indices (TMB, starts @ 0) to use for averaging recruitment when defining SSB_XSPR (if XSPR_R_opt = 2,4)
  DATA_VECTOR(FXSPR_init); // annual initial values to use for newton steps to find FXSPR (n_years_model+n_proj_years)
  DATA_VECTOR(FMSY_init); // annual initial values to use for newton steps to find FMSY (n_years_model+n_proj_years)
  
  // data for one-step-ahead (OSA) residuals
  DATA_INTEGER(do_osa); //whether to do osa residuals. For efficiency reasons with age comp likelihoods.
  DATA_VECTOR(obsvec); // vector of all observations for OSA residuals
  DATA_VECTOR_INDICATOR(keep, obsvec); // for OSA residuals
  DATA_IMATRIX(keep_C); // indices for catch obs, can loop years/fleets with keep(keep_C(y,f))
  DATA_IMATRIX(keep_I);
  DATA_IMATRIX(keep_E); // Ecov
  DATA_IARRAY(keep_Cpaa);
  DATA_IARRAY(keep_Ipaa);
  DATA_IVECTOR(do_post_samp_N); //length = 5, whether to ADREPORT posterior residuals for NAA, M, selectivity, Ecov, q. 
  DATA_IVECTOR(do_post_samp_M); //length = 5, whether to ADREPORT posterior residuals for NAA, M, selectivity, Ecov, q. 
  DATA_IVECTOR(do_post_samp_sel); //length = 5, whether to ADREPORT posterior residuals for NAA, M, selectivity, Ecov, q. 
  DATA_IVECTOR(do_post_samp_Ecov); //length = 5, whether to ADREPORT posterior residuals for NAA, M, selectivity, Ecov, q. 
  DATA_IVECTOR(do_post_samp_q); //length = 5, whether to ADREPORT posterior residuals for NAA, M, selectivity, Ecov, q. 

  // data for environmental covariate(s), Ecov
  DATA_INTEGER(n_Ecov); // also = 1 if no Ecov
  DATA_INTEGER(n_years_Ecov); // num years in Ecov  process model
  DATA_IMATRIX(Ecov_use_obs); // all 0 if no Ecov
  DATA_MATRIX(Ecov_obs);
  DATA_IMATRIX(Ecov_how_R); // n_Ecov x n_stocks: specific to recruitment effects. 0 = no effect, 1 = controlling, 2 = limiting, 3 = lethal, 4 = masking, 5 = directive
  DATA_IMATRIX(use_Ecov_R); // n_Ecov x n_stocks: 0/1 values indicating to use effects on recruitment
  DATA_IARRAY(use_Ecov_M); // n_Ecov x n_stocks x n_ages x n_regions: 0/1 values indicating to use effects on natural mortality at age.
  DATA_IMATRIX(use_Ecov_q); // n_Ecov x n_indices: 0/1 values indicating to use effects on catchability for each index.
  DATA_IARRAY(use_Ecov_mu); // n_Ecov x n_stocks x n_ages x n_seasons x n_regions x n_regions-1: 0/1 values indicating to use effects on migration for each stock for each region (less 1).
  //DATA_IMATRIX(Ecov_where); // n_Ecov x 3+n_indices. 0/1 values with columns corresponding to recruit, mortality, migration, indices in that order
  DATA_IVECTOR(Ecov_model); // 0 = no Ecov, 1 = RW, 2 = AR1
  //DATA_IMATRIX(ind_Ecov_out_start); // n_Ecov x (2 + n_indices) index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects specific the multiple types of effects each Ecov can have)
  DATA_IMATRIX(ind_Ecov_out_start_R); // n_Ecov x n_stocks: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IARRAY(ind_Ecov_out_start_M); // n_Ecov x n_stocks x n_ages x n_regions: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IMATRIX(ind_Ecov_out_start_q); // n_Ecov x n_indices: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IARRAY(ind_Ecov_out_start_mu); // n_Ecov x n_stocks x n_regions-1: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  //DATA_IMATRIX(ind_Ecov_out_end); // n_Ecov x (2 + n_indices) index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects specific the multiple types of effects each Ecov can have)
  DATA_IVECTOR(ind_Ecov_out_end_R); // n_Ecov: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IARRAY(ind_Ecov_out_end_M); // n_Ecov x n_ages: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IMATRIX(ind_Ecov_out_end_q); // n_Ecov x n_indices: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IARRAY(ind_Ecov_out_end_mu); // n_Ecov x n_stocks x n_regions-1: index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects each Ecov can have)
  DATA_IVECTOR(Ecov_obs_sigma_opt); // n_Ecov, 1 = given, 2 = estimate 1 value, shared among obs, 3 = estimate for each obs, 4 = estimate for each obs as random effects
  DATA_IVECTOR(Ecov_use_re); // n_Ecov: 0/1: use Ecov_re? If yes, add to nll.

  // data for projections
  DATA_INTEGER(do_proj); // 1 = yes, 0 = no
  DATA_INTEGER(n_years_proj); // number of years to project
  DATA_INTEGER(n_years_proj_Ecov); // number of years to project Ecov
  DATA_IVECTOR(avg_years_ind); // model year indices (TMB, starts @ 0) to use for averaging MAA, waa, maturity, and F (if use.avgF = TRUE)
  DATA_IVECTOR(proj_F_opt); // for each projection year, 1 = last year F (default), 2 = average F, 3 = F at X% SPR, 4 = user-specified F, 5 = calculate F from user-specified catch
  DATA_VECTOR(proj_Fcatch); // user-specified F or catch in projection years, only used if proj_F_opt = 4 or 5
  DATA_INTEGER(proj_M_opt); // 1 = continue M_re (check for time-varying M_re on R side), 2 = average M (over avg_years_ind)
  DATA_VECTOR(logR_mean); // (n_stocks) empirical mean recruitment in model years, used for SCAA recruit projections
  DATA_VECTOR(logR_sd); //  (n_stocks) empirical sd recruitment in model years, used for SCAA recruit projections
  DATA_VECTOR(F_proj_init); // annual initial values  to use for newton steps to find F for use in projections  (n_years_proj)
  
  //static brp info
  DATA_INTEGER(which_F_age_static); // which age of F to use for full total F for static brps (max of average FAA_tot over avg_years_ind)
  DATA_SCALAR(static_FXSPR_init); // initial value to use for newton steps to find FXSPR_static
 
  // parameters - general
  PARAMETER_MATRIX(mean_rec_pars); //n_stocks x n_rec_pars (determined by recruit_model)
  PARAMETER_VECTOR(logit_q);
  PARAMETER_VECTOR(q_prior_re); //n_indices (if a prior is used for q, this is used instead of logit_q)
  PARAMETER_ARRAY(q_re); //n_years x n_indices (time series of)
  PARAMETER_MATRIX(q_repars) //n_indices x 2 (sigma, rho)
  PARAMETER_MATRIX(log_F); //n_years_model x n_fleets
  //PARAMETER_MATRIX(F_devs);
  //PARAMETER_VECTOR(mu); //migration parameters
  PARAMETER_ARRAY(mu_prior_re); //n_stocks x n_ages x n_seasons x n_regions x n_regions-1
  PARAMETER_ARRAY(trans_mu); //n_stocks x n_ages x n_seasons x n_regions x n_regions-1 (mean) migration parameters
  PARAMETER_ARRAY(mu_re); //n_stocks x n_ages x n_seasons x n_y x n_regions x n_regions-1 RE for migration
  PARAMETER_ARRAY(mu_repars); //n_stocks x n_ages x n_seasons x n_regions x 4 (sig, rho_a, rho_y, rho_r)
  //N1 might need some tweaking. for example, if all fish are forced to be in spawning region at the beginning of the year, then there should be no N1 in other regions.
  //PARAMETER_VECTOR(log_R1); //length must be consistent with R1_pointer above.
  PARAMETER_MATRIX(N1_repars); // (n_stocks x 3) mean, sig, rho
  PARAMETER_MATRIX(log_N1); // (n_stocks x n_ages)
  PARAMETER_MATRIX(log_NAA_sigma); // (n_stocks x n_ages) vector sigmas used with NAA_sigma_pointers
  PARAMETER_MATRIX(trans_NAA_rho); // (n_stocks x 2) rho_a, rho_y (length = 2)
  //Just have annual NAA right now
  PARAMETER_ARRAY(log_NAA); //(n_stocks x n_regions x nyears-1 x n_ages) 
  
  PARAMETER_MATRIX(logR_proj); // (n_stocks x n_proj_years) recruitment (random effects) in proj years, only if SCAA
  PARAMETER_MATRIX(logit_selpars); // mean selectivity, dim = n_selblocks x n_ages + 6 (n_ages for by age, 2 for logistic, 4 for double-logistic)
  PARAMETER_VECTOR(selpars_re);    // deviations in selectivity parameters (random effects), length = sum(n_selpars)*n_years per block
  PARAMETER_MATRIX(sel_repars);    // fixed effect parameters controlling selpars_re, dim = n_blocks, 3 (sigma, rho, rho_y)
  PARAMETER_MATRIX(catch_paa_pars); //n_fleets x 3
  PARAMETER_MATRIX(index_paa_pars); //n_indices x 3
  PARAMETER_VECTOR(log_catch_sig_scale) //n_fleets
  PARAMETER_VECTOR(log_index_sig_scale) //n_indices
  
  PARAMETER_ARRAY(M_a); // (n_stocks x n_regions x n_ages) mean M-at-age, fixed effects, length = n_ages if M_model = 2 (age-specific), length = 1 if M_model = 1 (constant) or 3 (weight-at-age M)
  PARAMETER_ARRAY(M_re); // random effects for year- and age-varying M deviations from mean M_a, dim = n_stocks x n_regions x n_years x n_ages
  PARAMETER_ARRAY(M_repars); // parameters controlling M_re, (n_stocks x n_regions x length = 3 (sigma_M, rho_M_a, rho_M_y))
  PARAMETER_MATRIX(log_b); //n_stocks x n_regions (for M = a * W^b model)
  PARAMETER_ARRAY(L_re); // random effects for year- and region-varying extra/unknown mortality source, dim = n_years x n_regions
  PARAMETER_ARRAY(L_repars); // random effects for year- and region-varying extra/unknown mortality source, dim = n_years x n_regions

  // parameters - environmental covariate ("Ecov")
  PARAMETER_MATRIX(Ecov_re); // nrows = n_years_Ecov, ncol = N_Ecov
  //PARAMETER_ARRAY(Ecov_beta); // dim = (2 + n_stocks*n_regions + n_indices) x n_poly x n_ecov x n_ages
  PARAMETER_ARRAY(Ecov_beta_R); // dim = n_stocks x n_regions x n_ecov x n_poly, effects on recruitment, beta_R in eqns 4-5, Miller et al. (2016)
  PARAMETER_ARRAY(Ecov_beta_M); // dim = n_stocks x n_ages x n_regions x n_ecov x n_poly, effects on natural mortality
  PARAMETER_ARRAY(Ecov_beta_mu); // dim = n_stocks x n_ages x n_seasons x n_regions x (n_regions-1) x n_ecov x n_poly, effects on movement
  PARAMETER_ARRAY(Ecov_beta_q); // dim = n_indices x n_ecov x n_poly, effects on catchability
  PARAMETER_MATRIX(Ecov_process_pars); // nrows = RW: 2 par (Ecov1, sig), AR1: 3 par (mu, sig, phi); ncol = N_ecov
  PARAMETER_MATRIX(Ecov_obs_logsigma); // N_Ecov_years x n_Ecov. options: just given (data), or fixed effect(s)
  PARAMETER_MATRIX(Ecov_obs_logsigma_re); // N_Ecov_years x n_Ecov. columns of random effects used if Ecov_obs_sigma_opt = 4 
  PARAMETER_MATRIX(Ecov_obs_sigma_par); // ncol = N_Ecov, nrows = 2 (mean, sigma of random effects)
  Type nll = 0.0; //negative log-likelihood
  int n_years_pop = n_years_model + n_years_proj;
  
  
  /////////////////////////////////////////
  // Environmental covariate process model --------------------------------------
  matrix<Type> Ecov_x(n_years_Ecov + n_years_proj_Ecov, n_Ecov); // 'true' estimated Ecov (x_t in Miller et al. 2016 CJFAS)

  Ecov_x.setZero();
  matrix<Type> Ecov_x = get_Ecov(Ecov_model, Ecov_pars, Ecov_re, Ecov_use_re);
  if(Ecov_model.sum()>0) {
    vector<Type> nll_Ecov = get_nll_Ecov(Ecov_model, Ecov_pars, Ecov_re, Ecov_use_re);
    nll += nll_Ecov.sum();
    REPORT(nll_Ecov);
    SIMULATE if(do_simulate_Ecov_re == 1){
      Ecov_re = simulate_Ecov_re(Ecov_model, Ecov_pars, Ecov_re, Ecov_use_re);
      Ecov_x = get_Ecov(Ecov_model, Ecov_pars, Ecov_re, Ecov_use_re);
      REPORT(Ecov_re);
    }
    if(do_post_samp_Ecov == 1) ADREPORT(Ecov_re);
  }
  REPORT(Ecov_x);
  if(do_post_samp_Ecov == 1) ADREPORT(Ecov_re);
  /////////////////////////////////////////

  /////////////////////////////////////////
  // Environmental covariate observation model -------------------------------------
  //TODO: Ecov obs are not yet simulated in projection years!!!!!!!!
  Type nll_Ecov_obs = Type(0);
  Type nll_Ecov_obs_sig = Type(0); // Ecov obs sigma random effects (opt = 4)
  matrix<Type> Ecov_obs_sigma(n_years_Ecov, n_Ecov);
  for(int i = 0; i < n_Ecov; i++){
    for(int y = 0; y < n_years_Ecov; y++){
      if(Ecov_obs_sigma_opt(i) == 4){
        Type mu_logsigma = Ecov_obs_sigma_par(0,i);
        Type sd_logsigma = exp(Ecov_obs_sigma_par(1,i));
        nll_Ecov_obs_sig -= dnorm(Ecov_obs_logsigma_re(y,i), mu_logsigma, sd_logsigma, 1);
        SIMULATE if(do_simulate_Ecov_obs == 1) {
          Ecov_obs_logsigma_re(y,i) = rnorm(mu_logsigma, sd_logsigma);
        }
        Ecov_obs_sigma(y,i) = exp(Ecov_obs_logsigma_re(y,i));
      } else{
        Ecov_obs_sigma(y,i) = exp(Ecov_obs_logsigma(y,i));
      }
      if(Ecov_use_obs(y,i) == 1){
        nll_Ecov_obs -= keep(keep_E(y,i)) * dnorm(obsvec(keep_E(y,i)), Ecov_x(y,i), Ecov_obs_sigma(y,i), 1);
        nll_Ecov_obs -= keep.cdf_lower(keep_E(y,i)) * log(squeeze(pnorm(obsvec(keep_E(y,i)), Ecov_x(y,i), Ecov_obs_sigma(y,i))));
        nll_Ecov_obs -= keep.cdf_upper(keep_E(y,i)) * log(1.0 - squeeze(pnorm(obsvec(keep_E(y,i)), Ecov_x(y,i), Ecov_obs_sigma(y,i))));
        SIMULATE if(do_simulate_Ecov_obs == 1) {
          Ecov_obs(y,i) = rnorm(Ecov_x(y,i), Ecov_obs_sigma(y,i));
          obsvec(keep_E(y,i)) = Ecov_obs(y,i);
        }
      }
    }
  }
  nll += nll_Ecov_obs_sig;
  nll += nll_Ecov_obs;
  REPORT(nll_Ecov_obs_sig);
  SIMULATE if(do_simulate_Ecov_obs ==1) {
    REPORT(Ecov_obs);
    REPORT(Ecov_obs_logsigma);
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
    matrix<Type> tmp = get_Ecov_out(Ecov_x, n_years_model, n_years_proj, t_ind_s, t_ind_e);
    for(int j = 0; j < tmp.rows(); j++) for(int k = 0; k < n_Ecov; k++) Ecov_out_R(s,j,k) = tmp(j,k);
  }
  REPORT(Ecov_out_R);
  
  array<Type> Ecov_lm_R(n_stocks, n_regions, n_years_pop, n_Ecov); 
  int n_poly_R = Ecov_beta_R.dim(3); // now a 4D array dim: (n_stocks,n_regions, n_Ecov, n_poly)
  matrix<Type> Ecov_beta_R_s_r(n_Ecov,n_poly_R);
  for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++){
    for(int i = 0; i <n_Ecov; i++) for(int j = 0; i <n_poly_R; j++) Ecov_beta_R_s_r(i,j) = Ecov_beta_R(s,r,i,j);
    matrix<Type> tmp = get_Ecov_lm(Ecov_beta_R_s_r,Ecov_out_R);
    for(int y = 0; y < tmp.dim(0); y++) for(int i = 0; i <n_Ecov; i++) Ecov_lm_R(s,r,y,i) = tmp(y,i);
  }
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
    matrix<Type> tmp = get_Ecov_out(Ecov_x, n_years_model, n_years_proj, t_ind_s, t_ind_e);
    for(int j = 0; j < tmp.rows(); j++) for(int i = 0; i < n_Ecov; k++) Ecov_out_M(s,a,r,j,i) = tmp(j,i);
  }
  REPORT(Ecov_out_M);
  //M
  array<Type> Ecov_lm_M(n_stocks, n_regions, n_ages, n_years_pop, n_Ecov); 
  int n_poly_M = Ecov_beta_M.dim(4); // now a 5D array dim: (n_stocks, n_ages, n_regions, n_Ecov, n_poly)
  for(int a = 0; a < n_ages_model; a++){
    matrix<Type> Ecov_beta_M_s_r(n_Ecov,n_poly_M);
    for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++){
      for(int i = 0; i <n_Ecov; i++) for(int j = 0; i <n_poly_M; j++) Ecov_beta_M_s_r(i,j) = Ecov_beta_M(s,a,r,i,j);
      matrix<Type> tmp = get_Ecov_lm(Ecov_beta_M_s_r,Ecov_out_M);
      for(int y = 0; y < tmp.dim(0); y++) for(int i = 0; i <n_Ecov; i++) Ecov_lm_M(s,r,a,y,i) = tmp(y,i);
    }
  }
  REPORT(Ecov_lm_M);
  
  //q
  array<Type> Ecov_out_q(n_indices, n_years_pop, n_Ecov);
  for(int i = 0; i < n_indices; i++){
    vector<int> t_ind_s = ind_Ecov_out_start_q.col(i);
    vector<int> t_ind_e = ind_Ecov_out_end_q.col(i);
    matrix<Type> tmp = get_Ecov_out(Ecov_x,n_years_model, n_years_proj, t_ind_s, t_ind_e);
    for(int j = 0; j < tmp.rows(); j++) for(int k = 0; k < n_Ecov; k++) Ecov_out_q(i,j,k) = tmp(j,k);
  }
  REPORT(Ecov_out_q);
  array<Type> Ecov_lm_q(n_indices,n_years_pop, n_Ecov); 
  int n_poly_q = Ecov_beta_q.dim(2); // now a 3D array dim: (n_indices, n_Ecov, n_poly)
  for(int a = 0; a < n_indices; a++){
    matrix<Type> Ecov_beta_q_a(n_Ecov,n_poly_q);
    //for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++){
      for(int i = 0; i <n_Ecov; i++) for(int j = 0; i <n_poly_q; j++) Ecov_beta_q_a(i,j) = Ecov_beta_q(a,i,j);
      matrix<Type> tmp = get_Ecov_lm(Ecov_beta_q_a,Ecov_out_q);
      for(int y = 0; y < tmp.dim(0); y++) for(int i = 0; i <n_Ecov; i++) Ecov_lm_q(a,y,i) = tmp(y,i);
    //}
  }
  REPORT(Ecov_lm_q);
  
  //mu
  array<Type> Ecov_out_mu(n_stocks, n_ages, n_seasons, n_regions, n_regions-1, n_years_pop, n_Ecov);
  for(int h = 0; h < n_stocks; h++) for(int i = 0; i < n_regions-1; i++){
    vector<int> t_ind_s(n_Ecov), t_int_e(n_Ecov);
    for(int l = 0; l < n_Ecov; l++) t_ind_s(l) = ind_Ecov_out_start_mu(i,j,l);
    for(int l = 0; l < n_Ecov; l++) t_ind_e(l) = ind_Ecov_out_end_mu(i,j,l);
    matrix<Type> tmp = get_Ecov_out(Ecov_x, n_years_model, n_years_proj, t_ind_s, t_ind_e);
    for(int j = 0; j < tmp.rows(); j++) for(int i = 0; i < n_Ecov; i++) Ecov_out_mu(h,i,j,k) = tmp(j,k);
  }
  REPORT(Ecov_out_mu);
  array<Type> Ecov_lm_mu(n_stocks, n_ages, n_seasons, n_regions, n_regions-1, n_years_pop, n_Ecov);
  int n_poly_mu = Ecov_beta_mu.dim(5); // now a 6D array dim: (n_stocks, n_ages, n_seasons, n_regions, n_regions-1, n_Ecov, n_poly)
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages_model; a++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions-1; rr++){
    matrix<Type> Ecov_beta_mu_s_r(n_Ecov,n_poly_mu);
    for(int i = 0; i <n_Ecov; i++) for(int j = 0; i <n_poly_mu; j++) Ecov_beta_mu_s_r(i,j) = Ecov_beta_mu(s,a,t,r,rr,i,j);
    matrix<Type> tmp = get_Ecov_lm(Ecov_beta_mu_s_r,Ecov_out_mu);
    for(int y = 0; y < tmp.dim(0); y++) for(int i = 0; i <n_Ecov; i++) Ecov_lm_mu(s,a,t,r,rr,y,i) = tmp(y,i);
  }
  REPORT(Ecov_lm_mu);
  
  // Calculate ecov link model (b1*ecov + b2*ecov^2 + ...) --------------------
  // ecov_beta is now 4D array, dim = (2 + n_indices) x n_poly x n_ecov x n_ages
  //int n_poly = Ecov_beta.dim(1); // now a 4D array dim: (n_effects,n_poly,n_Ecov,n_ages) is second dimension
  //vector<matrix<Type>> Ecov_lm(n_Ecov)(n_effects); // ecov linear model for each Ecov, dim = n_years_pop, n_ages
  // Ecov_lm.setZero();
  // Ecov_lm stores the linear models for each Ecov and where it is used. dim = n_Ecov, n_effects, n_years_pop, n_ages
  // n_effects dimension is: 0: recruitment, 1: M, 2-1+n_indices: which catchability it affects
  //array<Type> Ecov_lm(n_Ecov, n_effects,n_years_pop, n_ages); 
  /////////////////////////////////////////
  
  vector<int> any_index_age_comp(n_indices);
  vector<int> any_fleet_age_comp(n_fleets);
  for(int i = 0; i < P_dim; i++) I_mat(i,i) = 1.0;

  for(int i = 0; i < n_indices; i++)
  {
    any_index_age_comp(i) = 0;
    for(int y = 0; y < n_years; y++) if(use_index_paa(y,i) == 1) any_index_age_comp(i) = 1;
  }
  for(int i = 0; i < n_fleets; i++)
  {
    any_fleet_age_comp(i) = 0;
    for(int y = 0; y < n_years; y++) if(use_catch_paa(y,i) == 1) any_fleet_age_comp(i) = 1;
  }


  /////////////////////////////////////////
  // Selectivity --------------------------------------------------------------
  // need to change to use selectivity.hpp
  vector<array<Type> > selpars_re_mats(n_selblocks); // gets selectivity deviations (RE vector, selpars_re) as vector of matrices (nyears x npars), one for each block
  vector<matrix<Type> > selpars(n_selblocks); // selectivity parameter matrices for each block, nyears x npars
  Type nll_sel = 0.0;
  int istart = 0;
  int ct = 0;
  for(int b = 0; b < n_selblocks; b++){
    array<Type> tmp2(n_years_model,n_selpars(b));
    tmp2.setZero();
    selpars_re_mats(b) = tmp2;

    int jstart = 0; // offset for indexing selectivity pars, depends on selectivity model for block b: n_ages (age-specific) + 2 (logistic) + 4 (double-logistic)
    if(selblock_models(b) == 2) jstart = n_ages;
    if(selblock_models(b) == 3) jstart = n_ages + 2;

    if(selblock_models_re(b) > 1){
      // fill in sel devs from RE vector, selpars_re (fixed at 0 if RE off)
      array<Type> tmp(n_years_selblocks(b), n_selpars_est(b));
      for(int j=0; j<n_selpars_est(b); j++){
        tmp.col(j) = selpars_re.segment(istart,n_years_selblocks(b));
        istart += n_years_selblocks(b);
      }

      //question: is it faster here to just work on the selectivity parameters as re rather than the deviations?
      // likelihood of RE sel devs (if turned on)
      Type sigma; // sd selectivity deviations (fixed effect)
      Type rho; // among-par correlation selectivity deviations (fixed effect)
      Type rho_y; // among-year correlation selectivity deviations (fixed effect)
      Type Sigma_sig_sel;
      sigma = exp(sel_repars(b,0));
      rho = rho_trans(sel_repars(b,1)); // rho_trans ensures correlation parameter is between -1 and 1, see helper_functions.hpp
      rho_y = rho_trans(sel_repars(b,2)); // rho_trans ensures correlation parameter is between -1 and 1, see helper_functions.hpp
      
      if((selblock_models_re(b) == 2) | (selblock_models_re(b) == 5)){
        // 2D AR1 process on selectivity parameter deviations
        Sigma_sig_sel = pow(pow(sigma,2) / ((1-pow(rho_y,2))*(1-pow(rho,2))),0.5);
        nll_sel += SCALE(SEPARABLE(AR1(rho),AR1(rho_y)), Sigma_sig_sel)(tmp);
        SIMULATE if(do_simulate_sel_re == 1) SEPARABLE(AR1(rho),AR1(rho_y)).simulate(tmp);
      } else {
        // 1D AR1 process on selectivity parameter deviations
        if(selblock_models_re(b) == 3){ // ar1 across parameters in selblock, useful for age-specific pars.
          vector<Type> tmp0 = tmp.matrix().row(0); //random effects are constant across years 
          Sigma_sig_sel = pow(pow(sigma,2) / (1-pow(rho,2)),0.5);
          nll_sel += SCALE(AR1(rho), Sigma_sig_sel)(tmp0);
          SIMULATE if(do_simulate_sel_re == 1)
          {
            AR1(rho).simulate(tmp0);
            for(int i = 0; i < tmp0.size(); i++) tmp(0,i) = tmp0(i);
          }
        } else { // selblock_models_re(b) = 4, ar1_y
          vector<Type> tmp0 = tmp.matrix().col(0); //random effects are constant across years 
          Sigma_sig_sel = pow(pow(sigma,2) / (1-pow(rho_y,2)),0.5);
          //Sigma_sig_sel = sigma;
          nll_sel += SCALE(AR1(rho_y), Sigma_sig_sel)(tmp0);
          SIMULATE if(do_simulate_sel_re == 1)
          {
            AR1(rho_y).simulate(tmp0);
            tmp.col(0) = tmp0;
          }
        }
      }
      SIMULATE if(do_simulate_sel_re == 1) {
        tmp = tmp * Sigma_sig_sel;
        istart -= n_selpars_est(b) * n_years_selblocks(b); //bring it back to the beginning for this selblock
        for(int j=0; j<n_selpars_est(b); j++){
          for(int y = 0; y < n_years_selblocks(b); y++){
            selpars_re(istart) = tmp(y,j);
            istart++;
          }
        }
      }

      // construct deviations array with full dimensions (n_years_model instead of n_years_selblocks, n_selpars instead of n_selpars_est)
      for(int j=0; j<n_selpars(b); j++){
        for(int y=0; y<n_years_model; y++){
          if((selblock_years(y,b) == 1) & (selpars_est(b,j+jstart) > 0)){
            selpars_re_mats(b)(y,j) = selpars_re(ct);
            ct++;
          }
        }
      }
    }

    // get selpars = mean + deviations
    matrix<Type> tmp1(n_years_model, n_selpars(b));
    for(int j=jstart; j<(jstart+n_selpars(b)); j++){ // transform from logit-scale
      for(int i=0; i<n_years_model; i++){
        tmp1(i,j-jstart) = selpars_lower(b,j) + (selpars_upper(b,j) - selpars_lower(b,j)) / (1.0 + exp(-(logit_selpars(b,j) + selpars_re_mats(b).matrix()(i,j-jstart))));
      }
    }
    selpars(b) = tmp1;
  }
  REPORT(selpars);
  REPORT(sel_repars);
  REPORT(selpars_re); //even if not simulated
  if(do_post_samp_sel == 1) ADREPORT(selpars_re);
  REPORT(logit_selpars);
  REPORT(nll_sel);
  vector<matrix<Type>> selAA(n_selblocks); // selAA(b)(y,a) gives selectivity by block, year, age; selAA(b) is matrix with dim = n_years x n_ages;
  selAA = get_selectivity(n_years_model, n_ages, n_selblocks, selpars, selblock_models); // Get selectivity by block, age, year
  nll += nll_sel;
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
  matrix<Type> L = get_L(L_model, L_pars, L_re);
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
    matrix<Type> nll_q_re = get_nll_q_re(q_repars, q_re, use_q_re);
    nll += nll_q_re.sum();
    REPORT(nll_q_re);
    SIMULATE if(do_simulate_q_re ==1){
      q_re = simulate_q_re(q_repars, q_re, use_q_re);
    }
  }
  if(do_post_samp_q == 1) ADREPORT(q_re);
  matrix<Type> logit_q_mat = get_logit_q_mat(logit_q, q_re, q_prior_re, Ecov_lm_q, use_q_prior, use_q_re);
  REPORT(logit_q_mat);
  if(use_q_re.sum()>0) if(do_post_samp_q==0) ADREPORT(logit_q_mat);
  
  //n_years_pop x n_indices;
  matrix<Type> q = get_q(logit_q_mat, q_lower, q_upper);
  REPORT(q);
  //n_indices x n_years_pop x n_ages;
  array<Type> QAA = get_QAA(q,selAA, selblock_pointer_indices, q_lower, q_upper, n_years_model, n_ages);
  REPORT(QAA);
  /////////////////////////////////////////
  
  /////////////////////////////////////////
  //natural mortality 
  //RE and log_M_base (possibly updated in time steps)
  matrix<Type> nll_M = get_nll_M(M_repars, M_re_model, M_model, M_re, n_M_re);
  nll += nll_M.sum();
  REPORT(nll_M);
  SIMULATE if(do_simulate_M_re == 1){
    M_re = simulate_M_re(M_repars, M_re_model, M_model, M_re, n_M_re);
    REPORT(M_re);
  }
  if(do_post_samp_M == 1) ADREPORT(M_re);

  //log_b prior for M_model = 3
  if(use_log_b_prior==1){
    matrix<Type> nll_log_b = get_nll_log_b(log_b_model, log_b, bias_correct_pe);
    nll += nll_log_b.sum();
    REPORT(nll_log_b);
    SIMULATE if(do_simulate_M_re == 1){
      log_b = simulate_log_b(log_b_model, log_b, bias_correct);
      REPORT(log_b);
    }
    if(do_post_samp_M == 1) ADREPORT(log_b);
  }
  
  //n_stocks x n_regions x n_years x n_ages
  array<Type> log_M_base = get_log_M_base(M_re, M_model, n_years_model, M_a, log_b, waa, waa_pointer_M, Ecov_lm, use_Ecov, do_proj, proj_M_opt);
  REPORT(log_M_base);
  /////////////////////////////////////////
  
  /////////////////////////////////////////
  //movement
  // mu_model: 
  // 1 = constant across stocks, ages, time (n_seasons x n_regions x (n_regions -1) pars). 
  // 2 = differ by age (n_seasons x n_regions x (n_regions -1) fixed effects, n_ages AR1 random effects for each). 
  // 3 = differ by year (n_seasons x n_regions x (n_regions -1) fixed effects, n_years AR1 random effects for each)
  // 4 = differ by stock (n_seasons x n_stocks x n_regions x (n_regions -1) pars). 
  // 5 = differ by stock, age (n_seasons x n_stocks x n_regions x (n_regions -1) fixed effects, n_ages AR1 random effects for each). 
  // 6 = differ by stock, year (n_seasons x n_stocks x n_regions x (n_regions -1) fixed effects, n_years AR1 random effects for each)
  //priors, RE, get full lm link for migration parameters
  if(n_regions>1){//only do this mess if the number of regions is greater than 1.
    if(use_mu_prior.sum()>0) {
      array<Type> nll_mu_prior = get_nll_mu_prior(mu_prior_re, trans_mu, trans_mu_prior_sigma, use_mu_prior, mu_model);
      nll += nll_mu_prior.sum();
      REPORT(nll_mu_prior);
      SIMULATE if(do_simulate_mu_prior_re==1){
        mu_prior_re = simulate_mu_prior_re(mu_prior_re, trans_mu, trans_mu_prior_sigma, use_mu_prior, mu_model);
        REPORT(mu_prior_re);
      }
    }
    if(use_mu_re.sum()>0) {
      array<Type> nll_mu_re = get_nll_mu_re(mu_repars, mu_re, mu_model);
      nll += nll_mu_re.sum();
      REPORT(nll_mu_re);
      SIMULATE if(do_simulate_mu_re ==1){
        mu_re = simulate_mu_re(mu_repars, mu_re, mu_model);
      }
    }
    if(do_post_samp_mu == 1) ADREPORT(mu_re);
    //n_stocks x n_ages x n_seasons x n_years_pop X n_regions x n_regions-1
    //continues random processes in projection years!
    array<Type> trans_mu_base = get_trans_mu_base(trans_mu, mu_re, mu_prior_re, use_mu_prior, mu_model, Ecov_lm_mu, use_Ecov_mu);
    REPORT(trans_mu_base);
    //n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions
    //rows sum to 1 for mig_type = 0 (prob move), rows sum to 0 for mig_type 1 (instantaneous)
    array<Type> mu = get_mu(trans_mu_base, can_move, must_move, mig_type);
    REPORT(mu);
    //if(use_mu_re.sum()>0) if(do_post_samp_mu==0) ADREPORT(trans_mu_full);
  }
  /////////////////////////////////////////


  /////////////////////////////////////////
  // Construct fishing mortality-at-age (FAA)
  matrix<Type> F = exp(log_F); //n_years_model x n_fleets
  REPORT(F);
  //n_fleets x n_years_model x n_seasons x n_ages
  array<Type> FAA = get_FAA(F,fleet_seasons, selAA, selblock_pointer_fleets, n_ages, n_seasons);
  REPORT(FAA);
  //n_regions x n_years_model x n_seasons x n_ages
  array<Type> FAA_tot = get_FAA_tot(FAA, fleet_regions, n_regions);
  REPORT(FAA_tot);
  /////////////////////////////////////////
  
  /////////////////////////////////////////
  //Population model and likelihoods
  //First: initial numbers at age
  if(N1_model ==2) { //Initial numbers at age are random effects
    vector<Type> nll_N1 = get_nll_N1(log_N1, N1_repars);
    nll += nll_N1.sum();
    REPORT(nll_N1);
    SIMULATE if(do_simulate_N==1){
      log_N1 = simulate_N1(N1_model, log_N1, N1_repars);
      REPORT(log_N1);
    }
  }
  //n_stocks x n_regions x n_ages
  array<Type> pred_N1 = get_pred_N1(N1_model, log_N1, N1_repars);
  REPORT(pred_N1);

  //fill out NAA (they are all parameters)
  //n_stocks x n_regions x n_years_pop x n_ages
  array<Type> NAA = get_NAA(N1_model, log_N1, log_NAA, NAA_where) //log_NAA should be mapped accordingly to exclude NAA=0 e.g., recruitment by region.
  REPORT(NAA);
  
  //fill out initial year predicted numbers at age
  array<Type> pred_NAA = NAA;
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++){
    pred_NAA(s,r,0,a) = pred_N1(s,r,a);
  }
  
  //get probability transition matrices for yearly survival, movement, capture...
  //also get annual NAA at spawning and NAA for each index along the way.
  array<Type> annual_Ps(n_stocks,n_years_pop,n_ages,P_dim,P_dim);
  array<Type> annual_Ps_SSB(n_stocks,n_years_pop,n_ages,n_regions,n_regions); //just survival categories
  array<Type> P_season_terminal(n_seasons,n_stocks,n_ages,P_dim,P_dim);
  matrix<Type> I_mat(P_dim,P_dim);
  for(int i = 0; i < P_dim; i++) I_mat(i,i) = 1.0;
  array<Type> NAA_SSB(n_stocks,n_years_pop,n_ages);
  NAA_SSB.setZero();
  array<Type> NAA_index(n_stocks,n_indices,n_years_pop,n_ages);
  NAA_index.setZero();
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years_model; y++) for(int a = 0; a < n_ages; a++) 
  {
    matrix<Type> P_y = I_mat; //reset for each year, age, stock
    for(int t = 0; t < n_seasons; t++) 
    {
      //get numbers at age a for stock s in each region at time of spawning
      if(t == spawn_seasons(s)-1)
      {
        //P(0,t) x P(t_s-t): PTM over interval from beginning of year to time of spawning within the season
        matrix<Type> P_SSB = P_y * get_P(a, y, s, t, fleet_regions, can_move, mig_type, fracyr_SSB(y,s), FAA, log_M_base, mu, L);
        for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) annual_Ps_SSB(s,y,a,i,j) = P_SSB(i,j);
        //spawners only spawning in stock region...
        for(int r = 0; r < n_regions; r++) NAA_SSB(s,y,a) += NAA(s,r,y,a) * P_SSB(r,spawn_regions(s)-1);
      }
      for(int i = 0; i < n_indices; i++) 
      {
        if(t == index_seasons(i)-1)
        { 
          //P(0,t) x P(t_i-t): PTM over interval from beginning of year to time of index within the season
          matrix<Type> P_index = P_y * get_P(a, y, s, t, fleet_regions, can_move, mig_type, fracyr_indices(y,i), FAA, log_M_base, mu, L);
          for(int r = 0; r < n_regions; r++) NAA_index(s,i,y,a) += P_index(r,index_regions(i)-1) * NAA(s,y,a,r);
        }
      }
      //P(t,u): PTM over entire season interval
      matrix<Type> P_t = get_P(a, y, s, t, fleet_regions, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
      if(y == n_years_model-1)//collect each seasonal matrix for the last model year?
      {
        for(int d = 0; d < P_dim; d++) for(int dd = 0; dd < P_dim; dd++) P_season_terminal(t,s,a,d,dd) = P_t(d,dd);
      }
      P_y = P_y * P_t;
    }
    for(int i = 0; i < P_dim; i++) for(int j = 0; j < P_dim; j++) annual_Ps(s,y,a,i,j) = P_y(i,j);
  }
  //std::cout << "here 5"  << "\n";

  //get SSB from NAA_SSB, waa and mature
  //n_years_pop x n_stocks
  //make sure mature has projection years if necessary
  matrix<Type> SSB = get_SSB(NAA_ssb,waa,waa_pointer_ssb,mature);
  
  array<Type> pred_NAA = get_pred_NAA(NAA, SSB, spawn_regions, NAA_where, mean_rec_pars, annual_Ps);

  //need to be careful here about log_NAA random effects in regions where there will be zero predicted due to migration parameterization
  //FIX ME:
  array<Type> nll_NAA = get_NAA_nll(NAA, pred_NAA, NAA_repars, NAA_re_model, NAA_where);
  REPORT(nll_NAA);
  nll += nll_NAA.sum();

  vector<Type> nll_agg_catch(n_fleets);
  vector<Type> nll_catch_acomp(n_fleets);
  nll_agg_catch.setZero();
  nll_catch_acomp.setZero();
  //vector<Type> sigma2_log_NAA = exp(log_NAA_sigma*two);
  array<Type> pred_CAA(n_years_pop,n_fleets,n_ages);
  pred_CAA.setZero();
  array<Type> pred_stock_CAA(n_years_pop,n_stocks,n_ages,n_regions);
  pred_stock_CAA.setZero(); //(n_years,n_stocks,n_ages,n_regions);
  array<Type> pred_stock_catch(n_years_pop,n_stocks,n_regions);
  pred_stock_catch.setZero(); //(n_years,n_stocks,n_regions);
  matrix<Type> pred_catch_region(n_years_pop,n_regions);
  pred_catch_region.setZero(); //(n_years,n_regions);
  array<Type> pred_stock_prop_catch(n_years_pop,n_stocks,n_regions);
  pred_stock_prop_catch.setZero(); //(n_years,n_stocks,n_regions);
  array<Type> pred_catch_paa(n_years_mode + n_years_projl,n_fleets,n_ages);
  matrix<Type> pred_catch(n_years_pop,n_fleets);
  pred_catch.setZero();
  matrix<Type> log_pred_catch(n_years_pop,n_fleets);
  vector<Type> t_paa(n_ages)
  vector<Type> t_pred_paa(n_ages);
  //int P_dim = n_regions + n_fleets + 1;

  for(int y = 0; y < n_years; y++)
  {
    int acomp_par_count = 0;
    for(int f = 0; f < n_fleets; f++)
    {
      //pred_catch(y,f) = zero;
      Type tsum = zero;
      for(int a = 0; a < n_ages; a++) 
      {
        for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) 
        {
          for(int rr = 0; rr < n_regions; rr++) if(fleet_regions(f)-1 == rr)
          {
            pred_stock_CAA(y,s,a,rr) += NAA(s,y,a,r) * annual_Ps(s,y,a,r,n_regions + f);
            pred_stock_catch(y,s,rr) += waa(waa_pointer_fleets(f)-1,y,a) * NAA(s,y,a,r) * annual_Ps(s,y,a,r,n_regions + f);
            pred_catch_region(y,rr) += NAA(s,y,a,r) * annual_Ps(s,y,a,r,n_regions + f) * waa(waa_pointer_fleets(f)-1,y,a);
          }
          pred_CAA(y,f,a) +=  NAA(s,y,a,r) * all_P(s,y,a,r,n_regions + f);
        }
        pred_catch(y,f) += waa(waa_pointer_fleets(f)-1,y,a) * pred_CAA(y,f,a);
        tsum += pred_CAA(y,f,a);
      }
      nll_agg_catch(f) -= dnorm(log(agg_catch(y,f)),log(pred_catch(y,f)),agg_catch_sigma(y,f),1);
      if(any_fleet_age_comp(f) == 1)
      {
        vector<Type> acomp_pars(n_age_comp_pars_fleets(f));
        for(int j = 0; j < n_age_comp_pars_fleets(f); j++) 
        {
          acomp_pars(j) = catch_paa_pars(acomp_par_count);
          acomp_par_count++;
        }
        if(use_catch_paa(y,f) == 1) 
        {
          for(int a = 0; a < n_ages; a++)
          {
            pred_catch_paa(y,f,a) = pred_CAA(y,f,a)/tsum;
            t_pred_paa(a) = pred_catch_paa(y,f,a);
            t_paa(a) = catch_paa(f * n_years + y,a);
          }
          //std::cout << "for y = " << y << ", get_acomp_ll: " << get_acomp_ll(y, n_ages, catch_Neff(y,f), age_comp_model_fleets(f), t_paa, t_pred_paa, acomp_pars, catch_aref(y,f)) << "\n";
          nll_catch_acomp(f) -= get_acomp_ll(y, n_ages, catch_Neff(y,f), age_comp_model_fleets(f), t_paa, t_pred_paa, acomp_pars, catch_aref(y,f));
        }
      }
    }
  }
  for(int y = 0; y < n_years; y++) for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++)
  {    
    pred_stock_prop_catch(y,s,r) = pred_stock_catch(y,s,r)/pred_catch_region(y,r);
  }
  
  std::cout << "nll_agg_catch: " << "\n" << nll_agg_catch << "\n";
  nll += nll_agg_catch.sum();
  std::cout << "nll_catch_acomp: " << "\n" << nll_catch_acomp << "\n";
  nll += nll_catch_acomp.sum();

  vector<Type> nll_agg_indices(n_indices), nll_index_acomp(n_indices);
  nll_agg_indices.setZero();
  nll_index_acomp.setZero();
  array<Type> pred_IAA(n_years_pop,n_indices,n_ages);
  pred_IAA.setZero();
  array<Type> pred_index_paa(n_years_pop,n_indices,n_ages);
  matrix<Type> pred_indices(n_years_pop,n_indices);
  pred_indices.setZero();
  for(int y = 0; y < n_years; y++) 
  {
    int acomp_par_count = 0;
    for(int i = 0; i < n_indices; i++) 
    {
      Type tsum = zero;
      for(int a = 0; a < n_ages; a++) 
      {
        //get numbers at age a for stock s in each region at time of spawning
        for(int s = 0; s < n_stocks; s++) pred_IAA(y,i,a) += QAA(y,i,a) * NAA_index(s,i,y,a);
        if(units_indices(i) == 1) pred_indices(y,i) += waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        else pred_indices(y,i) += pred_IAA(y,i,a);    
      }
      for(int a = 0; a < n_ages; a++) 
      {
        if(units_index_paa(i) == 1) pred_IAA(y,i,a) = waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        tsum += pred_IAA(y,i,a);
      }
      if(use_indices(y,i) == 1)
      {
        nll_agg_indices(i) -= dnorm(log(agg_indices(y,i)),log(pred_indices(y,i)),agg_index_sigma(y,i),1);
      }
      if(any_index_age_comp(i) == 1)
      {
        vector<Type> acomp_pars(n_age_comp_pars_indices(i));
        for(int j = 0; j < n_age_comp_pars_indices(i); j++) 
        {
          acomp_pars(j) = index_paa_pars(acomp_par_count);
          acomp_par_count++;
        }
        if(use_index_paa(y,i) > 0)
        {
          for(int a = 0; a < n_ages; a++)
          {
            pred_index_paa(y,i,a) = pred_IAA(y,i,a)/tsum;
            t_pred_paa(a) = pred_index_paa(y,i,a);
            t_paa(a) = index_paa(i * n_years + y,a);
          }
          /*if(i == 1) 
          {
            Type acomp_ll = get_acomp_ll(y, n_ages, index_Neff(y,i), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(y,i));
            std::cout << "y: " << y << "\n";
            std::cout << index_Neff(y,i) << "\n";
            std::cout << age_comp_model_indices(i) << "\n";
            std::cout << "acomp_pars: " << "\n" << acomp_pars << "\n";
            std::cout << "index_aref(y,i): " << index_aref(y,i) << "\n";
            std::cout << "t_paa: " << "\n" << t_paa << "\n";
            std::cout << "t_pred_paa: " << "\n" << t_pred_paa << "\n";
            std::cout << "acomp_ll: " << acomp_ll << "\n";
          }*/
          nll_index_acomp(i) -= get_acomp_ll(y, n_ages, index_Neff(y,i), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(y,i));
        }
      }
    }
  }
  std::cout << "nll_agg_indices: " << "\n" << nll_agg_indices << "\n";
  nll += nll_agg_indices.sum();
  std::cout << "nll_index_acomp: " << "\n" << nll_index_acomp << "\n";
  nll += nll_index_acomp.sum();
  
  matrix<Type> log_SSB(n_years_pop,n_stocks);
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years_pop; y++) log_SSB(y,s) = log(SSB(y,s));
  
  //if(reportMode==0){
    REPORT(NAA);
    REPORT(pred_NAA);
    REPORT(SSB);
    REPORT(selblocks);
    REPORT(q);
    REPORT(F);
    REPORT(pred_catch);
    REPORT(pred_indices);
    REPORT(NAA_index);
    REPORT(QAA);
    ADREPORT(log_SSB);
    REPORT(pred_CAA);
    REPORT(pred_stock_CAA);
    REPORT(pred_stock_catch);
    REPORT(pred_stock_prop_catch);
    REPORT(pred_catch_region);
    REPORT(all_P);
    REPORT(P_season_terminal);
  //}
  
  return nll;
}

