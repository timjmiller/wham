template <class Type>
matrix<Type> get_nll_catch_acomp(array<Type> pred_catch_paa, matrix<int> use_catch_paa, array<Type> catch_paa,
  matrix<Type> catch_Neff, vector<int> age_comp_model_fleets, matrix<Type> catch_paa_pars, 
  array<int> keep_Cpaa, data_indicator<vector<Type>, Type> keep, vector<Type> obsvec, vector<int> agesvec, int do_osa){
  int n_fleets = pred_catch_paa.dim(0);
  int n_y = catch_paa.dim(1);
  int n_ages = pred_catch_paa.dim(2);
  matrix<Type> nll_catch_acomp(n_y,n_fleets);
  nll_catch_acomp.setZero();

  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_y; y++)if(use_catch_paa(y,f)) {
    vector<Type> paa_obs_y(n_ages);
    vector<Type> t_pred_paa(n_ages);
    for(int a = 0; a < n_ages; a++){
      t_pred_paa(a) = pred_catch_paa(f,y,a);
      paa_obs_y(a) = catch_paa(f,y,a);
    }
    //NB: indexing in obsvec MUST be: keep_Cpaa(i,y,0),...,keep_Cpaa(i,y,0) + keep_Cpaa(i,y,1) - 1
    //keep_Cpaa(i,y,0) is first val, keep_Cpaa(i,y,1) is the length of the vector
    vector<Type> tf_paa_obs = obsvec.segment(keep_Cpaa(f,y,0), keep_Cpaa(f,y,1));
    vector<int> ages_obs_y = agesvec.segment(keep_Cpaa(f,y,0), keep_Cpaa(f,y,1));
    nll_catch_acomp(y,f) -= get_acomp_ll(tf_paa_obs, t_pred_paa, catch_Neff(y,f), ages_obs_y, age_comp_model_fleets(f), 
      vector<Type>(catch_paa_pars.row(f)), keep.segment(keep_Cpaa(f,y,0),keep_Cpaa(f,y,1)), do_osa, paa_obs_y);
  }
  return nll_catch_acomp;
}

template <class Type>
vector<Type> simulate_catch_paa_in_obsvec(vector<Type> obsvec, vector<int> agesvec, array<Type> pred_catch_paa, matrix<int> use_catch_paa,
  array<int> keep_Cpaa, matrix<Type> catch_Neff, vector<int> age_comp_model_fleets, matrix<Type> catch_paa_pars, int trace = 0){
  if(trace) see("in simulate_catch_paa_in_obsvec");
  int n_fleets = pred_catch_paa.dim(0);
  if(trace) see(n_fleets);
  int n_y = keep_Cpaa.dim(1);
  if(trace) see(n_y);
  int n_ages = pred_catch_paa.dim(2);
  if(trace) see(n_ages);
  vector<Type> obsvec_out = obsvec;
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_y; y++) if(use_catch_paa(y,f)) {
    vector<Type> t_pred_paa(n_ages);
    for(int a = 0; a < n_ages; a++) t_pred_paa(a) = pred_catch_paa(f,y,a);
    vector<int> ages_obs_y = agesvec.segment(keep_Cpaa(f,y,0), keep_Cpaa(f,y,1));
    vector<Type> tf_paa_obs = sim_acomp(t_pred_paa, catch_Neff(y,f), ages_obs_y, age_comp_model_fleets(f), 
      vector<Type>(catch_paa_pars.row(f)));
    obsvec_out.segment(keep_Cpaa(f,y,0),keep_Cpaa(f,y,1)) = tf_paa_obs;
    if((f == 0) & (y == 0)){
      if(trace) see(f);
      if(trace) see(y);
      if(trace) see(t_pred_paa);
      if(trace) see(ages_obs_y);
      if(trace) see(tf_paa_obs);
      if(trace) see(keep_Cpaa(f,y,0));
      if(trace) see(keep_Cpaa(f,y,1));
      if(trace) see(obsvec.segment(keep_Cpaa(f,y,0),keep_Cpaa(f,y,1)));
      if(trace) see(obsvec_out.segment(keep_Cpaa(f,y,0),keep_Cpaa(f,y,1)));
    }
  }
  if(trace) see("end simulate_catch_paa_in_obsvec");
  return obsvec_out;
}

template <class Type>
array<Type> sim_obsvec_in_catch_paa(vector<Type> obsvec, vector<int> agesvec, array<Type> catch_paa, matrix<int> use_catch_paa, array<int> keep_Cpaa, 
  vector<int> age_comp_model_fleets, int trace = 0){
  if(trace) see("in sim_obsvec_in_catch_paa");
  int n_fleets = catch_paa.dim(0);
  if(trace) see(n_fleets);
  int n_y = catch_paa.dim(1);
  if(trace) see(n_y);
  int n_ages = catch_paa.dim(2);
  if(trace) see(n_ages);
  array<Type> catch_paa_out = catch_paa;
  if(trace) see(catch_paa_out.dim);
  vector<Type> paa_obs_y(n_ages);
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_y; y++) if(use_catch_paa(y,f)) {
    if(trace) see(f);
    if(trace) see(y);
    vector<Type> tf_paa_obs = obsvec.segment(keep_Cpaa(f,y,0),keep_Cpaa(f,y,1));
    if(trace) see(tf_paa_obs);
    vector<int> ages_obs_y = agesvec.segment(keep_Cpaa(f,y,0), keep_Cpaa(f,y,1));
    if(trace) see(ages_obs_y);
    paa_obs_y = make_paa(tf_paa_obs, age_comp_model_fleets(f), ages_obs_y, n_ages);
    if(trace) see(paa_obs_y);
    for(int a = 0; a < n_ages; a++) catch_paa_out(f,y,a) = paa_obs_y(a);
  }
  if(trace) see("end sim_obsvec_in_catch_paa");
  return catch_paa_out;
}

template <class Type>
vector<Type> sim_agg_catch_in_obsvec(vector<Type> obsvec, matrix<int> keep_C, matrix<Type> agg_catch, matrix<int> use_agg_catch){
  int n_y = agg_catch.rows();
  int n_fleets = agg_catch.cols();
  vector<Type> obsvec_out = obsvec;
  for(int y = 0; y < n_y; y++) for(int f = 0; f < n_fleets; f++){
    if(use_agg_catch(y,f)) obsvec_out(keep_C(y,f)) = log(agg_catch(y,f));
  }
  return obsvec_out;
}


template <class Type>
matrix<Type> simulate_agg_catch(matrix<Type> pred_log_catch, matrix<Type> agg_catch, matrix<Type> agg_catch_sigma, 
  vector<Type> log_catch_sig_scale, matrix<int> use_agg_catch){
  int n_y = agg_catch.rows();
  int n_fleets = agg_catch.cols();
  matrix<Type> agg_catch_out = agg_catch;
  for(int y = 0; y < n_y; y++) for(int f = 0; f < n_fleets; f++) if(use_agg_catch(y,f)) {
    Type sig = agg_catch_sigma(y,f)*exp(log_catch_sig_scale(f));
    agg_catch(y,f) = exp(rnorm(pred_log_catch(y,f), sig));
  }
  return agg_catch;
}


template <class Type>
matrix<Type> get_nll_agg_catch(matrix<Type> pred_log_catch, matrix<Type> agg_catch_sigma, vector<Type> log_catch_sig_scale,
  vector<Type> obsvec, matrix<int> use_agg_catch, matrix<int> keep_C, data_indicator<vector<Type>, Type> keep){
  int n_y = agg_catch_sigma.rows();
  int n_fleets = agg_catch_sigma.cols();

  matrix<Type> nll_agg_catch(n_y,n_fleets);
  nll_agg_catch.setZero();

  for(int y = 0; y < n_y; y++) for(int f = 0; f < n_fleets; f++) if(use_agg_catch(y,f)){
    Type sig = agg_catch_sigma(y,f)*exp(log_catch_sig_scale(f));
    nll_agg_catch(y,f) -= keep(keep_C(y,f)) * dnorm(obsvec(keep_C(y,f)), pred_log_catch(y,f), sig,1);
    nll_agg_catch(y,f) -= keep.cdf_lower(keep_C(y,f)) * log(squeeze(pnorm(obsvec(keep_C(y,f)), pred_log_catch(y,f), sig)));
    nll_agg_catch(y,f) -= keep.cdf_upper(keep_C(y,f)) * log(1.0 - squeeze(pnorm(obsvec(keep_C(y,f)), pred_log_catch(y,f), sig)));
  }
  return nll_agg_catch;
}

template <class Type>
matrix<Type> get_pred_log_catch(matrix<Type> pred_catch, matrix<Type> agg_catch_sigma, vector<Type> log_catch_sig_scale, 
  int bias_correct_oe){
  int n_y = agg_catch_sigma.rows();
  int n_fleets = agg_catch_sigma.cols();
  matrix<Type> pred_log_catch(n_y,n_fleets);
  pred_log_catch.setZero();

  for(int y = 0; y < n_y; y++) for(int f = 0; f < n_fleets; f++) {
    pred_log_catch(y,f) = log(pred_catch(y,f));
    Type sig = agg_catch_sigma(y,f)*exp(log_catch_sig_scale(f));
    if(bias_correct_oe) pred_log_catch(y,f) -= 0.5*exp(2*log(sig));
  }
  return pred_log_catch;
}


//returns annual values of 
template <class Type>
vector<Type> get_log_FXSPR(Type percentSPR, array<Type> FAA, vector<int> fleet_regions, matrix<int> fleet_seasons, 
  vector<int> spawn_seasons, vector<int> spawn_regions, array<int> can_move, vector<int> mig_type, vector<Type> fracyr_seasons, 
  vector<int> which_F_age, 
  matrix<Type> fracyr_SSB, array<Type> log_M, array<Type> mu, matrix<Type> L, matrix<Type> log_SPR0, array<Type> waa_ssb, array<Type> mature, 
  vector<Type> SPR_weights, int SPR_weight_type, int small_dim, matrix<Type> R_XSPR, vector<Type> FXSPR_init, int trace = 0){

  int n_years_pop = waa_ssb.dim(1);
  vector<int> yvec(1);
  vector<Type> log_FXSPR(n_years_pop);

  for(int y = 0; y < n_years_pop; y++){
    yvec(0) = y;
    if(trace) see(y);
    array<Type> waa_ssb_y = get_matrix_y(waa_ssb, y);
    //array<Type> waa_catch_y = get_matrix_y(waa_catch, y);
    array<Type> mature_y = get_matrix_y(mature, y);
    array<Type> sel_y = get_avg_fleet_sel(FAA, yvec, which_F_age(y));
    if(trace) see(sel_y);
    vector<Type> L_y = L.row(y);
    if(trace) see(L_y);
    array<Type> log_M_y = get_log_M_y(y,log_M);
    if(trace) see(log_M_y);
    array<Type> mu_y = get_mu_y(y,mu);
    if(trace) see(mu_y);
    if(trace) see(mu_y.dim);

    vector<Type> FXSPR = get_FXSPR(spawn_seasons, spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, vector<Type> (fracyr_SSB.row(y)),
      sel_y, log_M_y, mu_y, L_y, mature_y,  waa_ssb_y, fracyr_seasons, vector<Type> (R_XSPR.row(y)), vector<Type> (log_SPR0.row(y)),
      percentSPR, SPR_weights, SPR_weight_type, small_dim, FXSPR_init(y), 10, trace);
    log_FXSPR(y) = log(FXSPR(0));
    if(trace) see(log_FXSPR(y));
  }

  return log_FXSPR;
}

template <class Type>
Type get_log_FXSPR_static(Type percentSPR, vector<int> avg_years_ind, array<Type> FAA, vector<int> fleet_regions, matrix<int> fleet_seasons, 
  vector<int> spawn_seasons, vector<int> spawn_regions, array<int> can_move, array<int> must_move, vector<int> mig_type, vector<Type> fracyr_seasons, 
  int which_F_age, matrix<Type> fracyr_SSB, array<Type> log_M, array<Type> trans_mu_base, matrix<Type> L, vector<Type> log_SPR0_static, array<Type> waa, 
  vector<int> waa_pointer_ssb, array<Type> mature, vector<Type> SPR_weights, int SPR_weight_type, int small_dim, array<Type> NAA, 
  vector<int> XSPR_R_avg_yrs, Type FXSPR_init, int trace = 0){
  //returns static FXSPR with inputs averaged over avg_years_ind years. Recruitment averaged over  XSPR_R_avg_yrs years

  array<Type> waa_ssb = get_avg_waa(waa, avg_years_ind, waa_pointer_ssb);
  if(trace) see(waa_ssb);
  array<Type> avg_mat = get_avg_mat(mature, avg_years_ind);
  if(trace) see(avg_mat);
  array<Type> avg_sel = get_avg_fleet_sel(FAA, avg_years_ind, which_F_age);
  if(trace) see(avg_sel);
  vector<Type> avg_L = get_avg_L(L, avg_years_ind, 0);
  if(trace) see(avg_L);
  array<Type> log_avg_M = get_avg_M(log_M, avg_years_ind, 1);
  if(trace) see(log_avg_M);
  array<Type> avg_mu = get_avg_mu(trans_mu_base,avg_years_ind, mig_type, can_move, must_move);
  if(trace) see(avg_mu);
  vector<Type> avg_fracyr_SSB = get_avg_ssbfrac(fracyr_SSB, avg_years_ind);
  if(trace) see(avg_fracyr_SSB);
  vector<Type> R_XSPR(NAA.dim(0));
  R_XSPR.setZero();
  for(int s = 0; s < NAA.dim(0); s++) for(int y = 0; y < XSPR_R_avg_yrs.size(); y++){
    R_XSPR(s) += NAA(s,spawn_regions(s)-1,XSPR_R_avg_yrs(y),0)/XSPR_R_avg_yrs.size();
  }

  vector<Type> FXSPR = get_FXSPR(spawn_seasons, spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, avg_fracyr_SSB,
    avg_sel, log_avg_M, avg_mu, avg_L, avg_mat, waa_ssb, fracyr_seasons, R_XSPR, log_SPR0_static,
    percentSPR, SPR_weights, SPR_weight_type, small_dim, FXSPR_init, 10, trace);
  Type log_FXSPR = log(FXSPR(0));
  if(trace) see(log_FXSPR);
  return log_FXSPR;
}

template <class Type>
array<Type> get_FAA_from_log_F(vector<Type> log_F, vector<int> which_F_age, array<Type> FAA){

  vector<int> yvec(1);
  int n_years_pop = log_F.size();
  int n_fleets = FAA.dim(0);
  int n_ages = FAA.dim(2);
  array<Type> FAA_from_F(n_fleets, n_years_pop, n_ages);
  for(int y = 0; y < n_years_pop; y++){
    yvec(0) = y;
    matrix<Type> sel_y = get_avg_fleet_sel(FAA, yvec, which_F_age(y));
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) FAA_from_F(f,y,a) = exp(log_F(y)) * sel_y(f,a);
  }
  return FAA_from_F;
}


template <class Type>
matrix<Type> get_log_YPR(matrix<int> fleet_seasons, vector<int> fleet_regions, 
  array<int> can_move, vector<int> mig_type, array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L, 
  array<Type> waa_catch, vector<Type> fracyr_seasons, int small_dim, int trace = 0){
  /* 
    returns a matrix n_years_pop x n_stocks of annual equilibrium spawning biomass per recruit (at age) by stock and region. 
    If movement is set up approriately all fish can be made to return to a single spawning region for each stock.
       spawn_seasons: vector of indicators telling which season spawning occurs for each stock
       spawn_regions: vector of indicators telling which season spawning occurs for each stock
       fleat_seasons: n_fleets x n_seasons matrix of indicators telling which season fleet is operating
       fleet_regions: n_fleets vector of indicators telling which region fleet is operating in
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
                FAA: n_fleets x n_years_pop x n_ages
              log_M: log M (density-independent components): n_stocks x n_regions x n_years_pop x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions array of movement matrices
                  L: n_years_pop x n_regions. "extra" unobserved mortality
            waa_catch: n_fleets x n_years_pop x n_ages_model. weight at age
     fracyr_seasons: n_seasons: length of intervals for each season
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_fleets = FAA.dim(0);
  //int n_seasons = can_move.dim(1);
  int n_years_pop = FAA.dim(1);
  int n_ages = FAA.dim(2);
  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  if(trace) see(fleet_regions);
  matrix<Type> FAA_y(n_fleets,n_ages);
  FAA_y.setZero();
  matrix<Type> log_YPR(n_years_pop,n_fleets);
  log_YPR.setZero();
  for(int y = 0; y < n_years_pop; y++){
    if(trace) see(y);
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) FAA_y(f,a) = FAA(f,y,a);
    if(trace) see(FAA_y);
    matrix<Type> waa_catch_y = get_matrix_y(waa_catch, y);
    if(trace) see(waa_catch_y);
    if(trace) see(waa_catch.dim);
    array<Type> mu_y = get_mu_y(y,mu);
    if(trace) see(mu_y);
    if(trace) see(mu.dim);
    if(trace) see(L);
    if(trace) see(L.rows());
    vector<Type> L_y = L.row(y);
    if(trace) see(L_y);
    array<Type> log_M_y = get_log_M_y(y,log_M);
    if(trace) see(log_M_y);
    if(trace) see(log_M.dim);

    array<Type> YPR_srf = get_YPR_srf(fleet_regions, fleet_seasons, can_move, mig_type, 
      FAA_y, log_M_y, mu_y, L_y, waa_catch_y, fracyr_seasons, 0, small_dim);
    if(trace) see(YPR_srf);
    for(int f = 0; f < n_fleets; f++) {
      for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) log_YPR(y,f) += YPR_srf(s,r,f);
      log_YPR(y,f) = log(log_YPR(y,f));
    } 
    if(trace) see(log_YPR.row(y));
  }
  return log_YPR;
}

template <class Type>
matrix<Type> get_log_SPR(vector<int> spawn_seasons, vector<int> spawn_regions, matrix<int> fleet_seasons, vector<int> fleet_regions, 
  array<int> can_move, vector<int> mig_type, matrix<Type> fracyr_SSB, array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L, 
  array<Type> mature, array<Type> waa_ssb, vector<Type> fracyr_seasons, int small_dim, int trace = 0){
  /* 
    returns a matrix n_years_pop x n_stocks of annual equilibrium spawning biomass per recruit (at age) by stock and region. 
    If movement is set up approriately all fish can be made to return to a single spawning region for each stock.
       spawn_seasons: vector of indicators telling which season spawning occurs for each stock
       spawn_regions: vector of indicators telling which season spawning occurs for each stock
       fleat_seasons: n_fleets x n_seasons matrix of indicators telling which season fleet is operating
       fleet_regions: n_fleets vector of indicators telling which region fleet is operating in
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_years_pop x n_stocks:  size of interval from beginning of season to time of spawning within that season
                FAA: n_fleets x n_years_pop x n_ages
              log_M: log M (density-independent components): n_stocks x n_regions x n_years_pop x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions array of movement matrices
                  L: n_years_pop x n_regions. "extra" unobserved mortality
             mature: n_stocks x n_years_pop x n_ages proportion mature at age
            waa_ssb: (n_stocks) x n_years_pop x n_ages_model. weight at age
     fracyr_seasons: n_seasons: length of intervals for each season
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_fleets = FAA.dim(0);
  //int n_seasons = can_move.dim(1);
  int n_ages = mature.dim(2);
  int n_stocks = mature.dim(0);
  int n_years_pop = mature.dim(1);
  if(trace) see(fleet_regions);
  matrix<Type> FAA_y(n_fleets,n_ages);
  FAA_y.setZero();
  matrix<Type> log_SPR(n_years_pop,n_stocks);
  log_SPR.setZero();
  for(int y = 0; y < n_years_pop; y++){
    if(trace) see(y);
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) FAA_y(f,a) = FAA(f,y,a);
    if(trace) see(FAA_y);
    matrix<Type> mature_y = get_matrix_y(mature, y);
    if(trace) see(mature_y);
    if(trace) see(mature.dim);
    matrix<Type> waa_ssb_y = get_matrix_y(waa_ssb, y);
    if(trace) see(waa_ssb_y);
    if(trace) see(waa_ssb.dim);
    array<Type> mu_y = get_mu_y(y,mu);
    if(trace) see(mu_y);
    if(trace) see(mu.dim);
    if(trace) see(L);
    if(trace) see(L.rows());
    vector<Type> L_y = L.row(y);
    if(trace) see(L_y);
    array<Type> log_M_y = get_log_M_y(y,log_M);
    if(trace) see(log_M_y);
    if(trace) see(log_M.dim);
    array<Type> SPR = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, 
      vector<Type> (fracyr_SSB.row(y)), FAA_y, log_M_y, mu_y, L_y, mature_y, waa_ssb_y, fracyr_seasons, 0, small_dim, 0);
    for(int s = 0; s < n_stocks; s++) log_SPR(y,s) = log(SPR(s,spawn_regions(s)-1,spawn_regions(s)-1)); 
    if(trace) see(log_SPR.row(y));
  }
  return log_SPR;
}

template <class Type>
matrix<Type> get_log_SPR0(vector<int> spawn_seasons, vector<int> spawn_regions, array<int> can_move, 
  vector<int> mig_type, matrix<Type> fracyr_SSB, array<Type> log_M, array<Type> mu, matrix<Type> L, 
  array<Type> mature, array<Type> waa_ssb, vector<Type> fracyr_seasons, int small_dim, int trace = 0){
  /* 
    returns a matrix n_years_pop x n_stocks of annual unfished equilibrium spawning biomass per recruit (at age) by stock and region. 
    If movement is set up approriately all fish can be made to return to a single spawning region for each stock.
       spawn_season: vector of indicators telling which season spawning occurs for each stock
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_years_pop x n_stocks:  size of interval from beginning of season to time of spawning within that season
         log_M: log M (density-independent components): n_stocks x n_regions x n_years_pop x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions array of movement matrices
                  L: n_years_pop x n_regions. "extra" unobserved mortality
             mature: n_stocks x n_years_pop x n_ages proportion mature at age
            waa_ssb: (n_stocks) x n_years_pop x n_ages_model. weight at age
     fracyr_seasons: n_seasons: length of intervals for each season
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_seasons = can_move.dim(1);
  int n_ages = mature.dim(2);
  int n_stocks = mature.dim(0);
  int n_years_pop = mature.dim(1);
  vector<int> fleet_regions(1);
  fleet_regions(0) = 1;
  if(trace) see(fleet_regions);
  matrix<int> fleet_seasons(1,n_seasons);
  fleet_seasons.setZero();
  matrix<Type> FAA0(1,n_ages);
  FAA0.setZero();
  if(trace) see(FAA0);
  matrix<Type> log_SPR0(n_years_pop,n_stocks);
  log_SPR0.setZero();
  for(int y = 0; y < n_years_pop; y++){
    if(trace) see(y);
    matrix<Type> mature_y = get_matrix_y(mature, y);
    if(trace) see(mature_y);
    if(trace) see(mature.dim);
    matrix<Type> waa_ssb_y = get_matrix_y(waa_ssb, y);
    if(trace) see(waa_ssb_y);
    if(trace) see(waa_ssb.dim);
    array<Type> mu_y = get_mu_y(y,mu);
    if(trace) see(mu_y);
    if(trace) see(mu.dim);
    if(trace) see(L);
    if(trace) see(L.rows());
    vector<Type> L_y = L.row(y);
    if(trace) see(L_y);
    array<Type> log_M_y = get_log_M_y(y,log_M);
    if(trace) see(log_M_y);
    if(trace) see(log_M.dim);
    array<Type> SPR_0 = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, 
      vector<Type> (fracyr_SSB.row(y)), FAA0, log_M_y, mu_y, L_y, mature_y, waa_ssb_y, fracyr_seasons, 0, small_dim, 0);
    for(int s = 0; s < n_stocks; s++) log_SPR0(y,s) = log(SPR_0(s,spawn_regions(s)-1,spawn_regions(s)-1)); 
    if(trace) see(log_SPR0.row(y));
  }
  return log_SPR0;
}

template <class Type>
vector<Type> get_log_SPR0_static(vector<int> avg_years_ind, vector<int> spawn_seasons, vector<int> spawn_regions, 
  array<int> can_move, array<int> must_move, vector<int> mig_type, matrix<Type> fracyr_SSB, array<Type> log_M, 
  array<Type> trans_mu_base, matrix<Type> L, array<Type> mature, array<Type> waa, vector<int> waa_pointer_ssb, 
  vector<Type> fracyr_seasons, int small_dim, int trace = 0){
  /* 
    returns a vector (n_stocks) of unfished equilibrium spawning biomass per recruit (at age) by stock. 
    If movement is set up approriately all fish can be made to return to a single spawning region for each stock.
       spawn_season: vector of indicators telling which season spawning occurs for each stock
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_years_pop x n_stocks:  size of interval from beginning of season to time of spawning within that season
         log_M: log M (density-independent components): n_stocks x n_regions x n_years_pop x n_ages
                 trans_mu_base: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions-1 array of movement matrices
                  L: n_years_pop x n_regions. "extra" unobserved mortality
             mature: n_stocks x n_years_pop x n_ages proportion mature at age
            waa_ssb: (n_stocks) x n_years_pop x n_ages_model. weight at age
     fracyr_seasons: n_seasons: length of intervals for each season
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  vector<int> fleet_regions(1);
  fleet_regions(0) = 1;
  if(trace) see(fleet_regions);
  matrix<int> fleet_seasons(1,fracyr_seasons.size());
  fleet_seasons.setZero();
  matrix<Type> FAA0(1,mature.dim(2));
  FAA0.setZero();
  if(trace) see(FAA0);
  matrix<Type> waa_ssb = get_avg_waa(waa, avg_years_ind, waa_pointer_ssb);
  if(trace) see(waa_ssb);
  matrix<Type> avg_mat = get_avg_mat(mature, avg_years_ind);
  if(trace) see(avg_mat);
  vector<Type> avg_L = get_avg_L(L, avg_years_ind, 0);
  if(trace) see(avg_L);
  array<Type> log_avg_M = get_avg_M(log_M, avg_years_ind, 1);
  if(trace) see(log_avg_M);
  array<Type> avg_mu = get_avg_mu(trans_mu_base,avg_years_ind, mig_type, can_move, must_move);
  if(trace) see(avg_mu);
  vector<Type> avg_fracyr_SSB = get_avg_ssbfrac(fracyr_SSB, avg_years_ind);
  if(trace) see(avg_fracyr_SSB);
  array<Type> SPR_0 = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, 
    avg_fracyr_SSB, FAA0, log_avg_M, avg_mu, avg_L, avg_mat, waa_ssb, fracyr_seasons, 0, small_dim, 0);
  vector<Type> log_SPR0(spawn_regions.size());
  for(int s = 0; s < spawn_regions.size(); s++) log_SPR0(s) = log(SPR_0(s,spawn_regions(s)-1,spawn_regions(s)-1)); 
  if(trace) see(log_SPR0);
  return log_SPR0;
}

template <class Type>
matrix<Type> get_P(int age, int year, int stock, int season, vector<int> fleet_regions, matrix<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, Type time, array<Type> FAA, array<Type>log_M, 
  array<Type> mu, matrix<Type> L) {
  /*
    produce the probability transition matrix for a given stock,, age, season, year
                age: which age
               year: which year
              stock: which stock
             season: which season
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
               time: time interval between beginning and end of the interval (units = years)
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                  L: n_years_model x n_regions; "extra" mortality rate
  */
  int n_regions = log_M.dim(1);
  int n_fleets = fleet_regions.size();
  int dim = n_regions+n_fleets+1;
  matrix<Type> P(dim,dim);
  P.setZero(); //zero it out.
  vector<Type> F(n_fleets), M(n_regions), Z(n_regions);
  F.setZero();
  M.setZero();
  Z.setZero();
  for(int r = 0; r < n_regions; r++) {
    M(r) = exp(log_M(stock,r,year,age)) + L(year,r); //add in "extra" mortality
    Z(r) = M(r);
  }
  for(int f = 0; f < n_fleets; f++) {
    if(fleet_seasons(f,season)) F(f) = FAA(year,f,age);
    Z(fleet_regions(f)-1) += F(f);
  }
  if(n_regions == 1){ //usual Baranov
    if(time < 1e-15) {
      P(0,0) = 1.0;
    } else {
      P(0,0) = exp(-Z(0)* time);
      for(int f = 0; f < n_fleets; f++) {
        P(0,1+f) = F(f) * (1 - exp(-Z(0) * time))/Z(0);
      }
      P(0,1+n_fleets) = M(0) * (1 - exp(-Z(0) * time))/Z(0);
    }
    for(int i = 1; i < dim; i++) P(i,i) = 1.0;
  } else{ //n_regions > 1
    matrix<int> can_move_now(n_regions,n_regions);
    can_move_now.setZero();
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) can_move_now(i,j) = can_move(stock,season,i,j);
    //if(n_mu(stock,year,season,age)>0) //migration is happening
    if(can_move_now.sum()>0) { //migration is happening
      if(mig_type(stock) == 0) { //migration is instantaneous after survival and mortality, so P is easy.
        if(time < 1e-15) {
        //prob of survival is 1 when interval is zero
        for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) P(i,j) =  mu(stock,age,season,year,i,j);
        } else {
          for(int f = 0; f < n_fleets; f++) P(fleet_regions(f)-1,n_regions + f) = F(f) * (1.0 - exp(-Z(fleet_regions(f)-1) * time))/Z(fleet_regions(f)-1);
          for(int i = 0; i < n_regions; i++) {
            //probs of survival and mortality and then moving between regions
            for(int j = 0; j < n_regions; j++) 
            {
              //survival only a function of region starting in, then modify survival in in each region by prob of moving.
              P(i,j) = exp(-Z(i) * time) * mu(stock,age,season,year,i,j); 
            }
            //caught only in region starting in since migration happens after survival 
            
            //all other dead
            P(i,dim-1) = M(i) * (1.0 - exp(-Z(i) * time))/Z(i);
          }
        }
        //fish and nat mort state: if dead, must stay dead
        for(int i = n_regions; i < dim; i++) P(i,i) = 1.0; 
      }
      if(mig_type(stock) == 1) {//migration occurs continuously during interval, so P is not so easy.
        //prob of survival is 1 when interval is zero
        if(time < 1e-15) {
          for(int i = 0; i < n_regions; i++) P(i,i) = 1.0;
        } else {
          matrix<Type> A(dim,dim);
          A.setZero(); //zero it out.
          for(int i = 0; i < n_regions; i++) A(i,dim-1) += M(i); //other dead
          for(int f = 0; f < n_fleets; f++) A(fleet_regions(f)-1,n_regions + f) += F(f);
          for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
            //need to not put in mu diagonal because it has -sum of mig rates in there already.
            if(i != j) A(i,j) += mu(stock,age,season,year,i,j); // transition intensities
          }
          for(int i = 0; i < n_regions; i++) A(i,i) = -(A.row(i)).sum(); //hazard
          A = A * time;
          P = expm(A);
        }
      }
    }
    else {//no migration during this interval, so P is easy.
      //prob of survival is 1 when interval is zero
      if(time < 1e-15) {
        for(int i = 0; i < dim; i++) P(i,i) = 1.0; 
      } else {
        for(int i = 0; i < n_regions; i++) 
        {
          P(i,i) = exp(-Z(i) * time);
          //other dead
          P(i,dim-1) = M(i) * (1.0 - exp(-Z(i) * time))/Z(i);
        }
        for(int f = 0; f < n_fleets; f++) 
        {
          P(fleet_regions(f)-1,n_regions+f) = F(f) * (1.0 - exp(-Z(fleet_regions(f)-1) * time))/Z(fleet_regions(f)-1);
        } 
      }
    }
  } //end n_regions > 1
  return(P);
}

template <class Type>
matrix<Type> get_P(int age, int stock, int season, vector<int> fleet_regions, matrix<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, Type time, matrix<Type> FAA, array<Type>log_M, 
  array<Type> mu, vector<Type> L) {
  /*
    produce the probability transition matrix for a given stock, age, season
                age: which age
              stock: which stock
             season: which season
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
               time: time interval between beginning and end of the interval (units = years)
                FAA: fishing mortality: n_fleets x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_regions x n_regions; movement rates
                  L: n_regions; "extra" mortality rate
  */
  array<Type> FAA_d(FAA.rows(),1, FAA.cols());
  FAA_d.setZero();
  for(int f = 0; f < FAA.rows(); f++) for(int a = 0; a < FAA.cols(); a++) FAA_d(f,0,a) = FAA(f,a);
  array<Type> log_M_d(log_M.dim(0),log_M.dim(1),1,log_M.dim(2));
  log_M_d.setZero();
  for(int s = 0; s < log_M.dim(0); s++) for(int a = 0; a < log_M.dim(2); a++) for(int r = 0; r < log_M.dim(1); r++) {
    log_M_d(s,r,0,a) = log_M(s,r,a);
  }
  array<Type> mu_d(mu.dim(0), mu.dim(1), mu.dim(2), 1, mu.dim(3),mu.dim(4));
  mu_d.setZero();
  for(int s = 0; s < mu.dim(0); s++) for(int a = 0; a < mu.dim(1); a++) for(int t = 0; t < mu.dim(2); t++) {
    for(int r = 0; r < mu.dim(3); r++) for(int rr = 0; rr < mu.dim(4); rr++){
      mu_d(s,a,t,0,r,rr) = mu_d(s,a,t,r,rr);
    }
  }
  matrix<Type> L_d(1,L.size());
  L_d.setZero();
  for(int r = 0; r < mu.dim(3); r++) L_d(0,r) = L(r);
  matrix<Type> P = get_P(age, 0, stock, season, fleet_regions, fleet_seasons, can_move, mig_type, time, FAA_d, log_M_d, mu_d, L_d);
  return(P);
}

template <class Type>
matrix<Type> get_S(matrix<Type> P, int n_regions){
  /*
    extract the submatrix from a PTM that contains the proportions surviving in each region
              P: the probablity transition matrix
      n_regions: the number of regions
  */
  matrix<Type> S(n_regions,n_regions);
  S.setZero();
  for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) S(i,j) = P(i,j);
  return(S);
}


//NOTE get_P here is defined as class T instead of Type and is not distiguishable from Type defined function.
//This is why the get_P function(s) in PTM are currently commented out. Not sure if this affects expected model performance.
template <class T>
matrix<T> get_P(int age, int year, int stock, int season, vector<int> fleet_regions, matrix<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, T time, array<T> FAA, array<T>log_M, 
  array<T> mu, matrix<T> L, int trace = 0) {
  /*
    produce the probability transition matrix for a given stock,, age, season, year
                age: which age
               year: which year
              stock: which stock
             season: which season
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
               time: time interval between beginning and end of the interval (units = years)
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                  L: n_years_model x n_regions; "extra" mortality rate
  */
  int n_regions = log_M.dim(1);
  int n_fleets = fleet_regions.size();
  int dim = n_regions+n_fleets+1;
  if(trace) see("inside get_P");
  matrix<T> P(dim,dim);
  P.setZero(); //zero it out.
  vector<T> F(n_fleets), M(n_regions), Z(n_regions);
  F.setZero();
  M.setZero();
  Z.setZero();
  for(int r = 0; r < n_regions; r++) {
    M(r) = exp(log_M(stock,r,year,age)) + L(year,r); //add in "extra" mortality
    Z(r) = M(r);
  }
  if(trace) see(M);
  for(int f = 0; f < n_fleets; f++) {
    if(fleet_seasons(f,season)) F(f) = FAA(f,year,age);
    Z(fleet_regions(f)-1) += F(f);
  }
  if(trace) see(F);
  if(trace) see(Z);
  if(trace) see(time);
  if(n_regions == 1){ //usual Baranov
    if(time < 1e-15) {
      P(0,0) = 1.0;
    } else {
      P(0,0) = exp(-Z(0)* time);
      for(int f = 0; f < n_fleets; f++) {
        P(0,1+f) = F(f) * (1 - exp(-Z(0) * time))/Z(0);
      }
      P(0,1+n_fleets) = M(0) * (1 - exp(-Z(0) * time))/Z(0);
    }
    for(int i = 1; i < dim; i++) P(i,i) = 1.0;
    if(trace) see(P);
  } else{ //n_regions > 1
    matrix<int> can_move_now(n_regions,n_regions);
    can_move_now.setZero();
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) can_move_now(i,j) = can_move(stock,season,i,j);
    //if(n_mu(stock,year,season,age)>0) //migration is happening
    if(can_move_now.sum()>0) { //migration is happening
      if(mig_type(stock) == 0) { //migration is instantaneous after survival and mortality, so P is easy.
        if(time < 1e-15) {
        //prob of survival is 1 when interval is zero
        for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) P(i,j) =  mu(stock,age,season,year,i,j);
        } else {
          for(int f = 0; f < n_fleets; f++) P(fleet_regions(f)-1,n_regions + f) = F(f) * (1.0 - exp(-Z(fleet_regions(f)-1) * time))/Z(fleet_regions(f)-1);
          for(int i = 0; i < n_regions; i++) {
            //probs of survival and mortality and then moving between regions
            for(int j = 0; j < n_regions; j++) 
            {
              //survival only a function of region starting in, then modify survival in in each region by prob of moving.
              P(i,j) = exp(-Z(i) * time) * mu(stock,age,season,year,i,j); 
            }
            //caught only in region starting in since migration happens after survival 
            
            //all other dead
            P(i,dim-1) = M(i) * (1.0 - exp(-Z(i) * time))/Z(i);
          }
        }
        //fish and nat mort state: if dead, must stay dead
        for(int i = n_regions; i < dim; i++) P(i,i) = 1.0; 
      }
      if(mig_type(stock) == 1) {//migration occurs continuously during interval, so P is not so easy.
        //prob of survival is 1 when interval is zero
        if(time < 1e-15) {
          for(int i = 0; i < n_regions; i++) P(i,i) = 1.0;
        } else {
          matrix<T> A(dim,dim);
          A.setZero(); //zero it out.
          for(int i = 0; i < n_regions; i++) A(i,dim-1) += M(i); //other dead
          for(int f = 0; f < n_fleets; f++) A(fleet_regions(f)-1,n_regions + f) += F(f);
          for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
            //need to not put in mu diagonal because it has -sum of mig rates in there already.
            if(i != j) A(i,j) += mu(stock,age,season,year,i,j); // transition intensities
          }
          for(int i = 0; i < n_regions; i++) A(i,i) = -(A.row(i)).sum(); //hazard
          A = A * time;
          P = expm(A);
        }
      }
    }
    else {//no migration during this interval, so P is easy.
      //prob of survival is 1 when interval is zero
      if(time < 1e-15) {
        for(int i = 0; i < dim; i++) P(i,i) = 1.0; 
      } else {
        for(int i = 0; i < n_regions; i++) 
        {
          P(i,i) = exp(-Z(i) * time);
          //other dead
          P(i,dim-1) = M(i) * (1.0 - exp(-Z(i) * time))/Z(i);
        }
        for(int f = 0; f < n_fleets; f++) 
        {
          P(fleet_regions(f)-1,n_regions+f) = F(f) * (1.0 - exp(-Z(fleet_regions(f)-1) * time))/Z(fleet_regions(f)-1);
        } 
      }
    }
  } //end n_regions > 1
  if(trace) see("end inside get_P");

  return(P);
}

//age, 0, stock, season, fleet_regions, fleet_seasons, can_move, mig_type, time
template <class T>
matrix<T> get_P_sta(int age, int stock, int season, vector<int> fleet_regions, matrix<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, T time, matrix<T> FAA, array<T>log_M, 
  array<T> mu, vector<T> L, int trace = 0) {
  /*
    this expands the dimensions of log_M, FAA, mu, L so that the get_P function defined above can be used for reference points
    produce the probability transition matrix for a given stock, age, season
                age: which age
              stock: which stock
             season: which season
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
               time: time interval between beginning and end of the interval (units = years)
                FAA: fishing mortality: n_fleets x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_regions x n_regions; movement rates
                  L: n_regions; "extra" mortality rate
  */
  array<T> FAA_d(FAA.rows(),1, FAA.cols());
  FAA_d.setZero();
  for(int f = 0; f < FAA.rows(); f++) for(int a = 0; a < FAA.cols(); a++) FAA_d(f,0,a) = FAA(f,a);
  if(trace) see(FAA_d);
  array<T> log_M_d(log_M.dim(0),log_M.dim(1),1,log_M.dim(2));
  log_M_d.setZero();
  for(int s = 0; s < log_M.dim(0); s++) for(int a = 0; a < log_M.dim(2); a++) for(int r = 0; r < log_M.dim(1); r++) {
    log_M_d(s,r,0,a) = log_M(s,r,a);
  }
  if(trace) see(log_M_d);
  array<T> mu_d(mu.dim(0), mu.dim(1), mu.dim(2), 1, mu.dim(3),mu.dim(4));
  mu_d.setZero();
  if(trace) see(mu.dim);
  if(trace) see(mu_d.dim);
  for(int s = 0; s < mu.dim(0); s++) for(int a = 0; a < mu.dim(1); a++) for(int t = 0; t < mu.dim(2); t++) {
    for(int r = 0; r < mu.dim(3); r++) for(int rr = 0; rr < mu.dim(4); rr++){
      mu_d(s,a,t,0,r,rr) = mu(s,a,t,r,rr);
    }
  }
  if(trace) see(mu_d);
  matrix<T> L_d(1,L.size());
  L_d.setZero();
  for(int r = 0; r < mu.dim(3); r++) L_d(0,r) = L(r);
  if(trace) see(L_d);
  matrix<T> P = get_P(age, 0, stock, season, fleet_regions, fleet_seasons, can_move, mig_type, time, FAA_d, log_M_d, mu_d, L_d, trace);
  if(trace) see(P);
  return P;
}

template <class Type>
array<Type> get_YPR(vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, matrix<Type> FAA, array<Type> log_M, array<Type> mu, vector<Type> L, matrix<Type> waacatch, 
  vector<Type> fracyr_seasons, int age_specific, int small_dim){
  /* 
    calculate equilibrium yield per recruit (at age) by stock and region.
        fleet_regions: vector of indicators telling which region each fleet is operating
        fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
             can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
             mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
                  FAA: fishing mortality: n_fleets x n_ages
                log_M: n_stocks x n_regions x n_ages
                   mu: n_stocks x n_ages x n_seasons x n_regions x n_regions array of movement matrices
                    L: n_regions. "extra" unobserved mortality
             waacatch: n_fleets x n_ages. weight at age
       fracyr_seasons: n_seasons: length of intervals for each season
         age_specific: 0/1 telling whether to return SSB/R by age or not. dimensions of returned arrays are different. 
            small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */
  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  int n_ages = log_M.dim(2);
  int n_fleets = FAA.rows();
  int n_seasons = can_move.dim(1);
  int P_dim = n_regions + n_fleets + 1;
  array<Type> YPRAA(n_stocks,n_ages,n_regions,n_fleets); //Yield/R at age in each fleet column, given recruited in region row
  YPRAA.setZero();
  matrix<Type> W(n_fleets,n_fleets);
  W.setZero();

  matrix<Type> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < P_dim; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<Type> S_ya = get_S(I, n_regions);
    for(int a = 0; a < n_ages-1; a++) {
      for(int i = 0; i < n_fleets; i++) W(i,i) = waacatch(i,a);
      matrix<Type> P_ya = I; //PTM for year and age
      for(int t = 0; t < n_seasons; t++) {
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      }
      // YPR at year and age = prob alive to up to age a-1 x prob caught at age a x waa
      matrix<Type> YPR_ya = S_ya * get_D(P_ya, n_regions,n_fleets) * W; //should be n_regions x n_fleets
      S_ya = S_ya * get_S(P_ya, n_regions);  //accumulate for next age
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) YPRAA(s,a,i,j) = YPR_ya(i,j);
    }
    //now plus group
    for(int i = 0; i < n_fleets; i++) W(i,i) = waacatch(i,n_ages-1);
    matrix<Type> P_ya = I; //PTM for year and age
    for(int t = 0; t < n_seasons; t++) {
      //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
      P_ya = P_ya * get_P(n_ages-1, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
    }
    matrix<Type> fundm(n_regions,n_regions);
    fundm.setZero();
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      fundm(i,j) = -P_ya(i,j);
      if(i==j) fundm(i,j) += 1; //I - S_y
    }
    if(small_dim) fundm = fundm.inverse(); else fundm = atomic::matinv(fundm);
    //for plus group S_ya = S_y,a-1 x (I - S_y,+)^-1
    S_ya = S_ya * fundm;
    // YPR at year and age = (prob alive to up to age a-1 + prob alive at older ages) x prob caught at these ages x waa
    matrix<Type> YPR_ya = S_ya * get_D(P_ya,n_regions,n_fleets) * W; //should be n_regions x n_fleets
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) YPRAA(s,n_ages-1,i,j) = YPR_ya(i,j);
  }
  if(age_specific) {
    return YPRAA; 
  } else { //sum across age
    array<Type> YPR(n_stocks,n_regions,n_fleets);
    YPR.setZero();
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) {
      YPR(s,i,j) += YPRAA(s,a,i,j);
    }
    return YPR;
  }
}

template <class Type>
array<Type> get_SPR(vector<int> spawn_seasons, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, vector<Type> fracyr_SSB, matrix<Type> FAA, array<Type> log_M, array<Type> mu, vector<Type> L, 
  matrix<Type> mature, matrix<Type> waa_ssb, vector<Type> fracyr_seasons, int age_specific, int small_dim){
  /* 
    calculate equilibrium spawning biomass per recruit (at age) by stock and region. If movement is set up approriately 
    all fish can be made to return to a single spawning region for each stock.
    NB: This version of get_SPR does not allow year indices. It is assumed that the yearly (or averaged) values are provided.
       spawn_season: vector of indicators telling which season spawning occurs for each stock
      fleet_regions: vector of indicators telling which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_stocks:  size of interval from beginning of season to time of spawning within that season
                FAA: fishing mortality: n_fleets x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_regions x n_regions array of movement matrices
                  L: n_regions. "extra" unobserved mortality
             mature: n_stocks x n_ages proportion mature at age
            waa_ssb: n_stocks x n_ages. weight at age
     fracyr_seasons: n_seasons: length of intervals for each season
       age_specific: 0/1 telling whether to return SSB/R by age or not. dimensions of returned arrays are different. 
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  int n_ages = log_M.dim(3);
  int n_fleets = FAA.rows();
  int n_seasons = can_move.dim(1);
  int P_dim = n_regions + n_fleets + 1;
  array<Type> SPRAA(n_stocks,n_ages,n_regions,n_regions); //SSB/R at age in each region column, given recruited in region row
  SPRAA.setZero();

  matrix<Type> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < P_dim; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<Type> S_ya = get_S(I, n_regions);
    for(int a = 0; a < n_ages-1; a++) {
      matrix<Type> P_ya = I, P_spawn = I; //PTM for age and up to time of spawning
      for(int t = 0; t < n_seasons; t++) {
        if(t == spawn_seasons(s)-1) {
          //P(0,t_spawn): PTM over entire year up to time of spawning
          P_spawn = P_ya * get_P_t(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(s), FAA, log_M, mu, L);
        }
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P_t(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      } 
      // SSB/R at year and age = prob alive to up to age a-1 x prob spawn at age a x waa x mature
      //should be n_regions x n_regions
      matrix<Type> SPR_ya = S_ya * get_S(P_spawn, n_regions) * mature(s,a) * waa_ssb(s,a);
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
        //SSB per Recruit in each region (should only have positive values in spawn_regions(s)-1?)
        SPRAA(s,a,i,j) += SPR_ya(i,j); 
      }
      S_ya = S_ya * get_S(P_ya, n_regions); //accumulate for next age
    }
    //now plus group
    matrix<Type> P_ya = I, P_spawn = I; //PTM for year and age and up to time of spawning
    for(int t = 0; t < n_seasons; t++) {
      if(t == spawn_seasons(s)-1) {
        //P(0,t_spawn): PTM over entire year up to time of spawning
        P_spawn = P_ya * get_P(n_ages-1, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(s), FAA, log_M, mu, L);
      }
      //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
      P_ya = P_ya * get_P(n_ages-1, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
    }
    matrix<Type> fundm(n_regions,n_regions);
    fundm.setZero();
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      fundm(i,j) = -P_ya(i,j);
      if(i==j) fundm(i,j) += 1;
    }
    if(small_dim) fundm = fundm.inverse(); else fundm = atomic::matinv(fundm);
    //for plus group S_ya = S_y,a-1 x (I - S_y,+)^-1
    S_ya = S_ya * fundm;
    // SSB/R at year and age = prob alive to up to age a-1 x prob spawn at age a x waa x mature
    //should be n_regions x n_regions
    matrix<Type> SPR_ya = S_ya * get_S(P_spawn, n_regions) * mature(s,n_ages-1) * waa_ssb(s,n_ages-1);
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      //SSB per Recruit in each region (should only have positive values in spawn_regions(s)-1?)
      SPRAA(s,n_ages-1,i,j) += SPR_ya(i,j); 
    }
  }
  if(age_specific) {
    return SPRAA; 
  } else {//sum over ages
    array<Type> SPR(n_stocks,n_regions,n_regions);
    SPR.setZero();
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      SPR(s,i,j) += SPRAA(s,a,i,j);
    }
    return SPR;
  }
}

template <class Type>
Type get_SPR(Type log_F, vector<Type> M, vector<Type> sel, vector<Type> mat, vector<Type> waassb, Type fracyearSSB)
{
  int n_ages = M.size();
  Type SPR = 0.0, ntemp = 1.0;
  vector<Type> F(n_ages), Z(n_ages);

  F = exp(log_F) * sel;
  Z = F + M;
  for(int age=0; age<n_ages-1; age++)
  {
    SPR += ntemp * mat(age) * waassb(age) * exp(-fracyearSSB * Z(age));
    ntemp *= exp(-Z(age));
  }
  ntemp /= 1.0-exp(-Z(n_ages-1));
  SPR += ntemp * mat(n_ages-1) * waassb(n_ages-1) * exp(-fracyearSSB*Z(n_ages-1));

  return SPR;
}

template <class Type>
Type get_SPR(Type log_F, vector<Type> M, matrix<Type> sel, vector<Type> mat, vector<Type> waassb, Type fracyearSSB)
{
  int n_ages = M.size();
  int n_fleets = sel.rows();
  Type SPR = 0.0, ntemp = 1.0;
  vector<Type> Z = M;
  matrix<Type> F = exp(log_F) * sel;
  //Z = F + M;
  for(int age=0; age<n_ages-1; age++)
  {
    for(int f = 0; f < n_fleets;f++) Z(age) += F(f,age);
    SPR += ntemp * mat(age) * waassb(age) * exp(-fracyearSSB * Z(age));
    ntemp *= exp(-Z(age));
  }
  for(int f = 0; f < n_fleets;f++) Z(n_ages-1) += F(f,n_ages-1);
  ntemp /= 1.0-exp(-Z(n_ages-1));
  SPR += ntemp * mat(n_ages-1) * waassb(n_ages-1) * exp(-fracyearSSB*Z(n_ages-1));

  return SPR;
}


template <class Type>
vector<Type> get_SPRAA(Type log_F, vector<Type> M, vector<Type> sel, vector<Type> mat, vector<Type> waassb, Type fracyearSSB)
{
  int n_ages = M.size();
  Type ntemp = 1.0;
  vector<Type> F(n_ages), Z(n_ages), SPR(n_ages);

  F = exp(log_F) * sel;
  Z = F + M;
  for(int age=0; age<n_ages-1; age++)
  {
    SPR(age) = ntemp * mat(age) * waassb(age) * exp(-fracyearSSB * Z(age));
    ntemp *= exp(-Z(age));
  }
  ntemp /= 1.0-exp(-Z(n_ages-1));
  SPR(n_ages-1) = ntemp * mat(n_ages-1) * waassb(n_ages-1) * exp(-fracyearSSB*Z(n_ages-1));

  return SPR;
}

template <class Type>
Type get_YPR(Type log_F, vector<Type> M, vector<Type> sel, vector<Type> waacatch)
{
  int n_ages = M.size();
  Type YPR = 0.0, ntemp = 1.0;
  vector<Type> F(n_ages), Z(n_ages);

  F = exp(log_F) * sel;
  Z = F + M;
  for(int age=0; age<n_ages-1; age++)
  {
    YPR += ntemp * F(age) * waacatch(age) * (1.0 - exp(-Z(age)))/Z(age);
    ntemp *= exp(-Z(age));
  }
  ntemp /= 1.0 - exp(-Z(n_ages-1));
  YPR += ntemp * F(n_ages-1) * waacatch(n_ages-1) * (1.0 - exp(-Z(n_ages-1)))/Z(n_ages-1);

  return YPR;
}

template <class Type>
Type get_YPR(Type log_F, vector<Type> M, matrix<Type> sel, matrix<Type> waacatch)
{
  int n_ages = M.size();
  int n_fleets = sel.rows();
  Type YPR = 0.0, ntemp = 1.0;
  vector<Type> Z = M;

  matrix<Type> F = exp(log_F) * sel;
  for(int age=0; age<n_ages; age++) for(int f = 0; f < n_fleets; f ++) Z(age) += F(f,age);
  for(int age=0; age<n_ages-1; age++) {
    for(int f = 0; f < n_fleets; f ++) YPR += ntemp * F(f,age) * waacatch(f,age) * (1.0 - exp(-Z(age)))/Z(age);
    ntemp *= exp(-Z(age));
  }
  ntemp /= 1.0 - exp(-Z(n_ages-1));
  for(int f = 0; f < n_fleets; f ++) YPR += ntemp * F(f,n_ages-1) * waacatch(f,n_ages-1) * (1.0 - exp(-Z(n_ages-1)))/Z(n_ages-1);
  return YPR;
}

template <class Type>
Type get_SPR_0(vector<Type> M, vector<Type> mat, vector<Type> waassb, Type fracyearSSB)
{
  int n_ages = M.size();
  Type SPR_0 = Type(0.0);
  Type ntemp0 = Type(1.0);
  for (int a = 0; a < n_ages - 1; a++)
  {
    SPR_0 += ntemp0 * mat(a) * waassb(a) * exp(-fracyearSSB * M(a));
    ntemp0 *= exp(-M(a));
  }
  ntemp0 /= Type(1.0)-exp(-M(n_ages-1));
  SPR_0 += ntemp0 * mat(n_ages-1) * waassb(n_ages-1) * exp(-fracyearSSB * M(n_ages-1));
  return SPR_0;
}

/* calculate beverton-holt or ricker equilibrium yield */
template<class Type>
struct sr_yield {
  /* Data and parameter objects for yield calculation: */
  Type SR_a;
  Type SR_b;
  vector<Type> M;
  vector<Type> sel;
  vector<Type> mat;
  vector<Type> waassb;
  vector<Type> waacatch;
  Type fracyearSSB;
  int sr_type;

  /* Constructor */
  sr_yield(Type SR_a_, Type SR_b_,
  vector<Type> M_,
  vector<Type> sel_,
  vector<Type> mat_,
  vector<Type> waassb_,
  vector<Type> waacatch_,
  Type fracyearSSB_, int sr_type_) :
    SR_a(SR_a_), SR_b(SR_b_), M(M_), sel(sel_), mat(mat_),
  waassb(waassb_), waacatch(waacatch_), fracyearSSB(fracyearSSB_), sr_type(sr_type_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_ages = M.size();
    T YPR = 0, SPR = 0, ntemp = 1, R;
    vector<T> F(n_ages), Z(n_ages);

    F = exp(log_F(0)) * sel.template cast<T>();
    Z = F + M.template cast<T>();
    for(int age=0; age<n_ages-1; age++)
    {
      YPR += ntemp * F(age) * T(waacatch(age)) * (1- exp(-Z(age)))/Z(age);
      SPR += ntemp * T(mat(age) * waassb(age)) * exp(-T(fracyearSSB) * Z(age));
      ntemp *= exp(-Z(age));
    }
    ntemp /= 1 - exp(-Z(n_ages-1));
    YPR += ntemp * F(n_ages-1) * T(waacatch(n_ages-1)) * (1 - exp(-Z(n_ages-1)))/Z(n_ages-1);
    SPR += ntemp * T(mat(n_ages-1) * waassb(n_ages-1)) * exp(-T(fracyearSSB)*Z(n_ages-1));

    //Type SPR = get_SPR(x, M, sel, mat, waassb, fracyearSSB);
    if(sr_type == 0) R = (T(SR_a) - 1/SPR) / T(SR_b); //beverton-holt
    if(sr_type == 1) R = log(T(SR_a) * SPR)/(T(SR_b) * SPR); //ricker
    T Y = YPR * R;
    return Y;
  }
};

/* calculate beverton-holt or ricker equilibrium yield */
/* uses selectivity and waa by fleet */
template<class Type>
struct sr_yield_fleet {
  /* Data and parameter objects for yield calculation: */
  Type SR_a;
  Type SR_b;
  vector<Type> M;
  matrix<Type> sel;
  vector<Type> mat;
  vector<Type> waassb;
  matrix<Type> waacatch;
  Type fracyearSSB;
  int sr_type;

  /* Constructor */
  sr_yield_fleet(Type SR_a_, Type SR_b_,
  vector<Type> M_,
  matrix<Type> sel_,
  vector<Type> mat_,
  vector<Type> waassb_,
  matrix<Type> waacatch_,
  Type fracyearSSB_, int sr_type_) :
    SR_a(SR_a_), SR_b(SR_b_), M(M_), sel(sel_), mat(mat_),
  waassb(waassb_), waacatch(waacatch_), fracyearSSB(fracyearSSB_), sr_type(sr_type_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_ages = M.size();
    int n_fleets = sel.rows();
    T YPR = 0, SPR = 0, ntemp = 1, R;
    vector<T> Z = M.template cast<T>();

    matrix<T> F = exp(log_F(0)) * sel.template cast<T>();
    //Z = F + M.template cast<T>();
    for(int age=0; age<n_ages; age++) for(int f = 0; f < n_fleets; f++) Z(age) += F(f,age);
    for(int age=0; age<n_ages-1; age++) {
      for(int f = 0; f < n_fleets; f++) YPR += ntemp * F(f,age) * T(waacatch(f,age)) * (1- exp(-Z(age)))/Z(age);
      SPR += ntemp * T(mat(age) * waassb(age)) * exp(-T(fracyearSSB) * Z(age));
      ntemp *= exp(-Z(age));
    }
    ntemp /= 1 - exp(-Z(n_ages-1));
    for(int f = 0; f < n_fleets; f++) YPR += ntemp * F(f,n_ages-1) * T(waacatch(f,n_ages-1)) * (1 - exp(-Z(n_ages-1)))/Z(n_ages-1);
    SPR += ntemp * T(mat(n_ages-1) * waassb(n_ages-1)) * exp(-T(fracyearSSB)*Z(n_ages-1));

    //Type SPR = get_SPR(x, M, sel, mat, waassb, fracyearSSB);
    if(sr_type == 0) R = (T(SR_a) - 1/SPR) / T(SR_b); //beverton-holt
    if(sr_type == 1) R = log(T(SR_a) * SPR)/(T(SR_b) * SPR); //ricker
    T Y = YPR * R;
    return Y;
  }
};

/* calculate SPR at F */
template<class Type>
struct spr_F {
  /* Data and parameter objects for calculation: */
  vector<Type> M;
  vector<Type> sel;
  vector<Type> mat;
  vector<Type> waassb;
  Type fracyearSSB;

  /* Constructor */
  spr_F(
  vector<Type> M_,
  vector<Type> sel_,
  vector<Type> mat_,
  vector<Type> waassb_,
  Type fracyearSSB_) :
    M(M_), sel(sel_), mat(mat_),
  waassb(waassb_), fracyearSSB(fracyearSSB_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_ages = M.size();
    T SPR = 0, ntemp = 1;
    vector<T> F(n_ages), Z(n_ages);

    F = exp(log_F(0)) * sel.template cast<T>();
    Z = F + M.template cast<T>();
    for(int age=0; age<n_ages-1; age++)
    {
      SPR += ntemp * T(mat(age) * waassb(age)) * exp(-T(fracyearSSB) * Z(age));
      ntemp *= exp(-Z(age));
    }
    ntemp /= 1 - exp(-Z(n_ages-1));
    SPR += ntemp * T(mat(n_ages-1) * waassb(n_ages-1)) * exp(-T(fracyearSSB)*Z(n_ages-1));

    return SPR;
  }
};

/* calculate SPR at F */
/* uses selectivity and by fleet */
template<class Type>
struct spr_F_fleet {
  /* Data and parameter objects for calculation: */
  vector<Type> M;
  matrix<Type> sel;
  vector<Type> mat;
  vector<Type> waassb;
  Type fracyearSSB;

  /* Constructor */
  spr_F_fleet(
  vector<Type> M_,
  matrix<Type> sel_,
  vector<Type> mat_,
  vector<Type> waassb_,
  Type fracyearSSB_) :
    M(M_), sel(sel_), mat(mat_),
  waassb(waassb_), fracyearSSB(fracyearSSB_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_ages = M.size();
    int n_fleets = sel.rows();
    T SPR = 0, ntemp = 1;
    vector<T> Z = M.template cast<T>();

    matrix<T> F = exp(log_F(0)) * sel.template cast<T>();
    for(int age=0; age<n_ages; age++) for(int f = 0; f < n_fleets; f++) Z(age) += F(f,age);
    //Z = F + M.template cast<T>();
    for(int age=0; age<n_ages-1; age++) {
      SPR += ntemp * T(mat(age) * waassb(age)) * exp(-T(fracyearSSB) * Z(age));
      ntemp *= exp(-Z(age));
    }
    ntemp /= 1 - exp(-Z(n_ages-1));
    SPR += ntemp * T(mat(n_ages-1) * waassb(n_ages-1)) * exp(-T(fracyearSSB)*Z(n_ages-1));

    return SPR;
  }
};

template <class Type>
matrix<Type> get_SPR_res(matrix<Type> MAA, array<Type> FAA, vector<int> which_F_age, array<Type> waa, int waa_pointer_ssb, vector<int> waa_pointer_fleets,
  matrix<Type> mature, Type percentSPR, vector<Type> predR, vector<Type> fracyr_SSB, vector<Type> log_SPR0, vector<Type> F_init)
{
  int n = 10;
  int ny = MAA.rows();
  int na = MAA.cols();
  int nf = waa_pointer_fleets.size();
  matrix<Type> res(ny, 5+n);
  vector<Type> log_FXSPR(ny), log_FXSPR_i(1), log_SPR_FXSPR(ny);
  vector<Type> log_SSB_FXSPR(ny), log_Y_FXSPR(ny);
  matrix<Type> log_FXSPR_iter(ny,n);
  //matrix<Type> catchWAA = extract_matrix_array3(waa, waa_pointer_tot_catch-1);
  matrix<Type> ssbWAA = extract_matrix_array3(waa, waa_pointer_ssb-1);
  vector<Type> log_YPR_FXSPR(ny), waassb(na), mat(na), M(na);
  matrix<Type> sel(nf,na), waacatch(nf, na);
  for(int y = 0; y < ny; y++)
  {
    log_FXSPR_iter(y,0) = res(y,5) = log(F_init(y));
    M = MAA.row(y);
    waassb = ssbWAA.row(y);
    waacatch = get_waacatch_y(waa, y, na, waa_pointer_fleets);
    //waacatch = catchWAA.row(y);
    mat = mature.row(y);
    Type F_full = 0;
    for(int f = 0; f < nf; f++) F_full += FAA(y,f,which_F_age(y)-1);
    for(int f = 0; f < nf; f++) for(int a = 0; a < na; a++) sel(f,a) = FAA(y,f,a)/F_full;
    spr_F_fleet<Type> sprF(M, sel, mat, waassb, fracyr_SSB(y));
    for (int i=0; i<n-1; i++)
    {
      log_FXSPR_i(0) = log_FXSPR_iter(y,i);
      vector<Type> grad_spr_F = autodiff::gradient(sprF,log_FXSPR_i);
      log_FXSPR_iter(y,i+1) = log_FXSPR_iter(y,i) - (sprF(log_FXSPR_i) - 0.01*percentSPR*exp(log_SPR0(y)))/grad_spr_F(0);// /hess_sr_yield(0,0);
      res(y,5+i+1) = log_FXSPR_iter(y,i+1);
    }
    log_FXSPR(y) = log_FXSPR_iter(y,n-1);
    log_SPR_FXSPR(y) = log(get_SPR(log_FXSPR(y), M, sel, mat, waassb, fracyr_SSB(y)));
    log_YPR_FXSPR(y) = log(get_YPR(log_FXSPR(y), M, sel, waacatch));
    log_SSB_FXSPR(y) = log(predR(y)) + log_SPR_FXSPR(y);
    log_Y_FXSPR(y) = log(predR(y)) + log_YPR_FXSPR(y);
    res(y,0) = log_FXSPR(y);
    res(y,1) = log_SSB_FXSPR(y);
    res(y,2) = log_Y_FXSPR(y);
    res(y,3) = log_SPR_FXSPR(y);
    res(y,4) = log_YPR_FXSPR(y);
  }
  return res;
}

template <class Type>
vector<Type> get_static_SPR_res(matrix<Type> MAA, array<Type> FAA, int which_F_age, array<Type> waa, int waa_pointer_ssb, vector<int> waa_pointer_fleets,
  matrix<Type> mature, Type percentSPR, matrix<Type> NAA, vector<Type> fracyr_SSB, Type F_init, 
  vector<int> years_M, vector<int> years_mat, vector<int> years_sel, vector<int> years_waa_ssb, vector<int> years_waa_catch, vector<int> years_R)
{
  int n = 10;
  int na = MAA.cols();
  int nf = waa_pointer_fleets.size();
  //matrix<Type> waacatch = extract_matrix_array3(waa, waa_pointer_tot_catch-1);
  matrix<Type> waassb = extract_matrix_array3(waa, waa_pointer_ssb-1);
  matrix<Type> sel(nf,na), waa_c(nf,na);
  vector<Type> M(na), mat(na), waa_s(na);
  sel.setZero(); M.setZero(); mat.setZero(); waa_s.setZero(), waa_c.setZero();
  Type ssbfrac = 0, R = 0;
  //get average inputs over specified years
  for(int y = 0; y < years_R.size(); y++) R += NAA(years_R(y),0)/years_R.size();
  for(int y = 0; y < years_mat.size(); y++) for(int a = 0; a < na; a++) mat(a) += mature(years_mat(y),a)/years_mat.size();
  for(int y = 0; y < years_M.size(); y++) for(int a = 0; a < na; a++) M(a) += MAA(years_M(y),a)/years_M.size();
  for(int y = 0; y < years_waa_catch.size(); y++) for(int f = 0; f < nf; f++) for(int a = 0; a < na; a++) {
    waa_c(f,a) += waa(waa_pointer_fleets(f)-1, years_waa_catch(y),a)/years_waa_catch.size();
  }
  for(int y = 0; y < years_waa_ssb.size(); y++) {
   ssbfrac += fracyr_SSB(years_waa_ssb(y))/years_waa_ssb.size();
   for(int a = 0; a < na; a++) waa_s(a) += waassb(years_waa_ssb(y),a)/years_waa_ssb.size();
  }
  for(int y = 0; y < years_sel.size(); y++) for(int a = 0; a < na; a++) for(int f = 0; f < nf; f++) {
    sel(f,a) += FAA(years_sel(y),f,a)/years_sel.size();
  }
  Type F_mult= 0;
  for(int f = 0; f < nf; f++) F_mult += sel(f,which_F_age-1);
  sel = sel/F_mult;
  Type spr0 = get_SPR_0(M, mat, waa_s, ssbfrac); 

  vector<Type> res(6+n), log_FXSPR_i(1), log_FXSPR_iter(n);
  log_FXSPR_iter(0) = res(6) = log(F_init);
  spr_F_fleet<Type> sprF(M, sel, mat, waa_s, ssbfrac);
  for (int i=0; i<n-1; i++)
  {
    log_FXSPR_i(0) = log_FXSPR_iter(i);
    vector<Type> grad_spr_F = autodiff::gradient(sprF,log_FXSPR_i);
    log_FXSPR_iter(i+1) = log_FXSPR_iter(i) - (sprF(log_FXSPR_i) - 0.01*percentSPR * spr0)/grad_spr_F(0);// /hess_sr_yield(0,0);
    res(6+i+1) = log_FXSPR_iter(i+1);
  }
  res(0) = log_FXSPR_iter(n-1); //log_FXSPR
  res(3) = log(get_SPR(res(0), M, sel, mat, waa_s, ssbfrac)); //log_SPR_FXSPR, should be log(X*SPR0/100)
  res(1) = log(R) + res(3); //log_SSB_FXSPR
  res(4) = log(get_YPR(res(0), M, sel, waa_c)); //log_YPR_FXSPR
  res(2) = log(R) + res(4); //log_Y_FXSPR
  res(5) = log(spr0); //log_SPR0
  return res;
}

/* calculate Catch at F */
template<class Type>
struct catch_F {
  /* Data and parameter objects for calculation: */
  vector<Type> M;
  vector<Type> sel;
  vector<Type> waacatch;
  vector<Type> NAA;

  /* Constructor */
  catch_F(
  vector<Type> M_,
  vector<Type> sel_,
  vector<Type> waacatch_,
  vector<Type> NAA_) :
    M(M_), sel(sel_), waacatch(waacatch_), NAA(NAA_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it achieves required catch
    int n_ages = M.size();
    T Catch = 0;
    vector<T> F(n_ages), Z(n_ages);

    F = exp(log_F(0)) * sel.template cast<T>();
    Z = F + M.template cast<T>();
    for(int age=0; age<n_ages; age++)
    {
      Catch += T(waacatch(age)) * T(NAA(age)) * F(age) *(1 - exp(-Z(age)))/Z(age);
    }
    return Catch;
  }
};

/* calculate Catch at F */
template<class Type>
struct log_catch_F {
  /* Data and parameter objects for calculation: */
  vector<Type> M;
  vector<Type> sel;
  vector<Type> waacatch;
  vector<Type> NAA;

  /* Constructor */
  log_catch_F(
  vector<Type> M_,
  vector<Type> sel_,
  vector<Type> waacatch_,
  vector<Type> NAA_) :
    M(M_), sel(sel_), waacatch(waacatch_), NAA(NAA_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it achieves required catch
    int n_ages = M.size();
    T Catch = 0;
    vector<T> F(n_ages), Z(n_ages);

    F = exp(log_F(0)) * sel.template cast<T>();
    Z = F + M.template cast<T>();
    for(int age=0; age<n_ages; age++)
    {
      Catch += T(waacatch(age)) * T(NAA(age)) * F(age) *(1 - exp(-Z(age)))/Z(age);
    }
    return log(Catch);
  }
};

/* calculate total catch at F across multiple fleets*/
template<class Type>
struct log_catch_fleets_F {
  /* Data and parameter objects for calculation: */
  vector<Type> M;
  matrix<Type> sel; //n_fleets x n_ages
  matrix<Type> waacatch; //n_fleets x n_ages
  vector<Type> NAA;

  /* Constructor */
  log_catch_fleets_F(
  vector<Type> M_,
  matrix<Type> sel_,
  matrix<Type> waacatch_,
  vector<Type> NAA_) :
    M(M_), sel(sel_), waacatch(waacatch_), NAA(NAA_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it achieves required catch
    int n_ages = M.size();
    int n_fleets = waacatch.rows();
    T Catch = 0;
    vector<T> Z = M.template cast<T>();

    matrix<T> F = exp(log_F(0)) * sel.template cast<T>();
    //Z = F + M.template cast<T>();
    for(int age=0; age<n_ages; age++) {
      for(int f = 0; f < n_fleets; f++) Z(age) += F(f,age);
      for(int f = 0; f < n_fleets; f++){
        Catch += T(waacatch(f,age)) * T(NAA(age)) * F(f,age) *(1 - exp(-Z(age)))/Z(age);
      }
    }
    return log(Catch);
  }
};

template <class Type>
Type get_F_from_Catch(Type Catch, vector<Type> NAA, vector<Type> M, vector<Type> sel, vector<Type> waacatch, Type F_init)
{
  int n = 10;
  vector<Type> log_F_i(1);
  vector<Type> log_F_iter(n);
  log_F_iter.fill(log(F_init)); //starting value
  catch_F<Type> catchF(M, sel, waacatch, NAA);
  for (int i=0; i<n-1; i++)
  {
    log_F_i(0) = log_F_iter(i);
    vector<Type> grad_catch_F = autodiff::gradient(catchF,log_F_i);
    log_F_iter(i+1) = log_F_iter(i) - (catchF(log_F_i) - Catch)/grad_catch_F(0);
  }
  Type res = exp(log_F_iter(n-1));
  return res;
}

template <class Type>
Type get_F_from_log_Catch(Type Catch, vector<Type> NAA, vector<Type> M, vector<Type> sel, vector<Type> waacatch, Type F_init)
{
  int n = 10;
  vector<Type> log_F_i(1);
  vector<Type> log_F_iter(n);
  log_F_iter.fill(log(F_init)); //starting value
  log_catch_F<Type> logcatchF(M, sel, waacatch, NAA);
  for (int i=0; i<n-1; i++)
  {
    log_F_i(0) = log_F_iter(i);
    vector<Type> grad_log_catch_F = autodiff::gradient(logcatchF,log_F_i);
    log_F_iter(i+1) = log_F_iter(i) - (logcatchF(log_F_i) - log(Catch))/grad_log_catch_F(0);
  }
  Type res = exp(log_F_iter(n-1));
  return res;
}

//multiple fleets
template <class Type>
Type get_F_from_log_Catch(Type Catch, vector<Type> NAA, vector<Type> M, matrix<Type> sel, matrix<Type> waacatch, Type F_init)
{
  int n = 10;
  vector<Type> log_F_i(1);
  vector<Type> log_F_iter(n);
  log_F_iter.fill(log(F_init)); //starting value
  log_catch_fleets_F<Type> logcatchF(M, sel, waacatch, NAA);
  for (int i=0; i<n-1; i++)
  {
    log_F_i(0) = log_F_iter(i);
    vector<Type> grad_log_catch_F = autodiff::gradient(logcatchF,log_F_i);
    log_F_iter(i+1) = log_F_iter(i) - (logcatchF(log_F_i) - log(Catch))/grad_log_catch_F(0);
  }
  Type res = exp(log_F_iter(n-1));
  return res;
}

template <class Type>
Type get_FXSPR(vector<Type> M, vector<Type> sel, vector<Type> waassb,
  vector<Type> mat, Type percentSPR, Type fracyr_SSB, Type log_SPR0, Type F_init)
{
  int n = 10;
  vector<Type> log_FXSPR_i(1);
  vector<Type> log_FXSPR_iter(n);
  log_FXSPR_iter.fill(log(F_init)); //starting value
  spr_F<Type> sprF(M, sel, mat, waassb, fracyr_SSB);
  for (int i=0; i<n-1; i++)
  {
    log_FXSPR_i(0) = log_FXSPR_iter(i);
    vector<Type> grad_spr_F = autodiff::gradient(sprF,log_FXSPR_i);
    log_FXSPR_iter(i+1) = log_FXSPR_iter(i) - (sprF(log_FXSPR_i) - 0.01*percentSPR*exp(log_SPR0))/grad_spr_F(0);// /hess_sr_yield(0,0);
  }
  Type res = exp(log_FXSPR_iter(n-1));
  return res;
}

template <class Type>
Type get_FMSY(Type log_a, Type log_b, vector<Type> M, vector<Type> sel, vector<Type> waacatch, vector<Type> waassb,
  vector<Type> mat, Type fracyr_SSB, Type log_SPR0, int recruit_model, Type F_init)
{    
  int n = 10;
  vector<Type> log_FMSY_i(1);
  vector<Type> log_FMSY_iter(n);
  log_FMSY_iter.fill(log(F_init)); //starting value
  Type a = exp(log_a);
  Type b = exp(log_b);
  int sr_type = 0; //recruit_model = 3, B-H
  if(recruit_model == 4) sr_type = 1; //recruit_model = 4, Ricker
  sr_yield<Type> sryield(a, b, M, sel, mat, waassb, waacatch, fracyr_SSB, sr_type);
  for (int i=0; i<n-1; i++)
  {
    log_FMSY_i(0) = log_FMSY_iter(i);
    vector<Type> grad_sr_yield = autodiff::gradient(sryield,log_FMSY_i);
    matrix<Type> hess_sr_yield = autodiff::hessian(sryield,log_FMSY_i);
    log_FMSY_iter(i+1) = log_FMSY_iter(i) - grad_sr_yield(0)/hess_sr_yield(0,0);
  }
  Type FMSY = exp(log_FMSY_iter(n-1));
  return FMSY;
}

//fleet-specific FAA
template <class Type>
Type get_FMSY(Type log_a, Type log_b, vector<Type> M, matrix<Type> sel, matrix<Type> waacatch, vector<Type> waassb,
  vector<Type> mat, Type fracyr_SSB, Type log_SPR0, int recruit_model, Type F_init)
{    
  int n = 10;
  vector<Type> log_FMSY_i(1);
  vector<Type> log_FMSY_iter(n);
  log_FMSY_iter.fill(log(F_init)); //starting value
  Type a = exp(log_a);
  Type b = exp(log_b);
  int sr_type = 0; //recruit_model = 3, B-H
  if(recruit_model == 4) sr_type = 1; //recruit_model = 4, Ricker
  sr_yield<Type> sryield(a, b, M, sel, mat, waassb, waacatch, fracyr_SSB, sr_type);
  for (int i=0; i<n-1; i++)
  {
    log_FMSY_i(0) = log_FMSY_iter(i);
    vector<Type> grad_sr_yield = autodiff::gradient(sryield,log_FMSY_i);
    matrix<Type> hess_sr_yield = autodiff::hessian(sryield,log_FMSY_i);
    log_FMSY_iter(i+1) = log_FMSY_iter(i) - grad_sr_yield(0)/hess_sr_yield(0,0);
  }
  Type FMSY = exp(log_FMSY_iter(n-1));
  return FMSY;
}

template <class Type>
Type get_pred_recruit_y(int y, int recruit_model, vector<Type> mean_rec_pars, vector<Type> SSB, matrix<Type> NAA, vector<Type> log_SR_a, 
  vector<Type> log_SR_b, matrix<int> Ecov_where, vector<int> Ecov_how, array<Type> Ecov_lm){

  /*
   * y: year (between 1 and n_years_model+n_years_proj)
   * recruit_model: which recruitment model (1-4)
   * mean_rec_pars: vector of any recruitment parameters (defined in main code)
   * SSB: vector of yearly SSB (uses y-1 for any S-R relationship)
   * NAA: matrix of numbers at age
   * log_SR_a: yearly "a" parameters for SR function
   * log_SR_b: yearly "b" parameters for SR function
   * Ecov_where: matrix of 0/1 with first column determining if Ecov is affecting recruitment
   * Ecov_how: integer vector with an element that tells how the Ecov is affecting recruitment
   * Ecov_lm: array that holds linear predictor for Ecov
   */
  //recruit_model == 1, random walk
  Type pred_recruit = NAA(y-1,0);
  if(recruit_model == 1) // random walk
  {
    //pred_NAA(y,0) = NAA(y-1,0);
  }
  else
  {
    if(recruit_model == 2) // random about mean
    {
      pred_recruit = exp(mean_rec_pars(0));
      int nE = Ecov_where.rows();
      for(int i=0; i < nE; i++){
        if(Ecov_where(i,0) == 1) if(Ecov_how(i) == 1) pred_recruit *= exp(Ecov_lm(i,0,y,0));
      }
      //pred_NAA(y,0) = exp(mean_rec_pars(0));
      //if(Ecov_recruit > 0) if(Ecov_how(Ecov_recruit-1) == 1) pred_NAA(y,0) *= exp(Ecov_lm(y,Ecov_recruit-1));
    }
    else
    {
      if(recruit_model == 3) // BH stock recruit (if ecov effect, already modified SR_a and/or SR_b)
      {
        pred_recruit = exp(log_SR_a(y)) * SSB(y-1)/(1 + exp(log_SR_b(y))*SSB(y-1));
        //pred_NAA(y,0) = exp(log_SR_a(y)) * SSB(y-1)/(1 + exp(log_SR_b(y))*SSB(y-1));
      }
      else // recruit_model = 4, Ricker stock recruit (if ecov effect, already modified SR_a and/or SR_b)
      {
        pred_recruit = exp(log_SR_a(y)) * SSB(y-1) * exp(-exp(log_SR_b(y)) * SSB(y-1));
        //pred_NAA(y,0) = exp(log_SR_a(y)) * SSB(y-1) * exp(-exp(log_SR_b(y)) * SSB(y-1));
      }
    }
  }
  return(pred_recruit);
}

template <class Type>
vector<Type> get_pred_NAA_y(int y, int recruit_model, vector<Type> mean_rec_pars, vector<Type> SSB, matrix<Type> NAA, vector<Type> log_SR_a, 
  vector<Type> log_SR_b, matrix<int> Ecov_where, vector<int> Ecov_how, array<Type> Ecov_lm, matrix<Type> ZAA){

  /*
   * y: year (between 1 and n_years_model+n_years_proj)
   * recruit_model: which recruitment model (1-4)
   * mean_rec_pars: vector of any recruitment parameters (defined in main code)
   * SSB: vector of yearly SSB (uses y-1 for any S-R relationship)
   * NAA: matrix of numbers at age
   * log_SR_a: yearly "a" parameters for SR function
   * log_SR_b: yearly "b" parameters for SR function
   * Ecov_where: matrix of 0/1 integer with first column determining if Ecov i affects recruitment (= 1)
   * Ecov_how: integer vector with an element that tells how the Ecov is affecting recruitment
   * Ecov_lm: array that holds linear predictor for Ecov
   * ZAA: matrix of total mortality rate by year and age
   */
  int n_ages = NAA.cols();
  vector<Type> pred_NAA(n_ages);
  
  // Expected recruitment
  pred_NAA(0) = get_pred_recruit_y(y, recruit_model, mean_rec_pars, SSB, NAA, log_SR_a, 
    log_SR_b, Ecov_where, Ecov_how, Ecov_lm);
    
  // calculate pred_NAA for ages after recruitment
  for(int a = 1; a < n_ages-1; a++) pred_NAA(a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
  pred_NAA(n_ages-1) = NAA(y-1,n_ages-2) * exp(-ZAA(y-1,n_ages-2)) + NAA(y-1,n_ages-1) * exp(-ZAA(y-1,n_ages-1));
  
  return(pred_NAA);
}

template <class Type>
Type get_SSB(matrix<Type> NAA, matrix<Type> ZAA, array<Type> waa, matrix<Type> mature, int y, int ssbpointer, vector<Type> fracyr_SSB){
  int na = mature.cols();
  Type SSB = 0;
  vector<Type> waay = get_waa_y(waa,y,na,ssbpointer);
  for(int a = 0; a < na; a++) SSB += NAA(y,a) * waay(a) * mature(y,a) * exp(-ZAA(y,a)*fracyr_SSB(y));
  return(SSB);
}
 
template <class Type>
matrix<Type> get_F_proj(int y, int n_fleets, vector<int> proj_F_opt, array<Type> FAA, matrix<Type> NAA, matrix<Type> MAA, matrix<Type> mature, 
  matrix<Type> waacatch, vector<Type> waassb, vector<Type> fracyr_SSB, vector<Type> log_SPR0, vector<int> avg_years_ind, 
  int n_years_model, vector<int> which_F_age, Type percentSPR, vector<Type> proj_Fcatch, Type percentFXSPR, Type F_init, 
  vector<Type> log_a, vector<Type> log_b, int recruit_model, Type percentFMSY){
    /* 
     get F to project for next time step
              y:  year of projection (>n_years_model)
     proj_F_opt:  for each year, how to specify F for projection (1 to 6)
            FAA:  FAA array from main code.
            NAA:  NAA matrix from main code
            MAA:  MAA matrix from main code.
         mature:  maturity matrix from main code.
       waacatch:  matrix (n_fleets x n_ages) weight at age in catch to use for year y (use function get_waacatch_y defined above)
         waassb:  weight at age for SSB to use in year y (use function get_waa_y defined above)
     fracyr_SSB:  vector of yearly fractions of the year when spawning occurs
       log_SPR0:  vector of yearly log(unfished SSB/R) 
  avg_years_ind:  integer vector of years to average F for projection
  n_years_model:  number of years before projection begins
    which_F_age:  define which age has max F
     percentSPR:  percentage (0-100) of unfished spawning potential to determine F_percentSPR
    proj_Fcatch:  vector (n_years_proj) of user specified Fishing mortality rates to project
   percentFXSPR:  percentage (0-100) of F_percentSPR to use in catch, e.g. GOM cod uses F = 75% F_40%SPR
         F_init:  initial value to use for FXSPR or FMSY newton method
          log_a:  annual log(a) for stock-recruit relationship
          log_b:  annual log(b) for stock-recruit relationship
  recruit_model:  integer for which type of recruit model is assumed (= 3 or 4 for using Fmsy)
    percentFMSY:  percentage (0-100) of FMSY to use in catch.
    */
  int n_toavg = avg_years_ind.size();
  int n_ages = waacatch.cols();
  int proj_F_opt_y = proj_F_opt(y-n_years_model);

  //proj_F_opt == 1, last year F (default)
  matrix<Type> FAA_proj(n_fleets, n_ages);
  vector<Type> FAA_tot_proj(n_ages);
  FAA_tot_proj.setZero();
  FAA_proj.setZero();
  if(proj_F_opt_y == 1){ // last year F (default)
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) FAA_proj(f,a) = FAA(n_years_model-1,f,a);
  }
  else { //proj_F_opt_y>1
    //option 1: average F is by fleet and Ftot is sum of fleet averages
    //when there is more than 1 fleet, the sum of predicted catch across fleets will not generally equal the total catch using FAA_tot and waa_totcatch.
    for(int f = 0; f < n_fleets; f++) 
    {
      for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_toavg; i++){
        FAA_proj(f,a) += FAA(avg_years_ind(i),f,a);
      }
      FAA_proj.row(f) /= Type(n_toavg);
    }
      
    //proj_F_opt_y == 2, average F
    if(proj_F_opt_y == 2){
      //already defined
    }
    //get selectivity using average over avg.yrs
    matrix<Type> sel_proj(n_fleets,n_ages);
    vector<Type> sel_tot_proj(n_ages);
    vector<Type> FAA_tot_proj = FAA_proj.colwise().sum();
    for(int a = 0; a < n_ages; a++) sel_tot_proj(a) = FAA_tot_proj(a)/FAA_tot_proj(which_F_age(y)-1);
    //selectivity at age and fleet as a proportion of full total F (sum across fleets and ages = 1).
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) sel_proj(f,a) = FAA_proj(f,a)/FAA_tot_proj(which_F_age(y)-1);
    vector<Type> M = MAA.row(y);
    vector<Type> NAA_y = NAA.row(y);
    vector<Type> mat = mature.row(y);
    Type fracyr_SSB_y = fracyr_SSB(y);
    Type log_SPR0_y = log_SPR0(y);
    Type log_a_y = log_a(y);
    Type log_b_y = log_b(y);
    if(proj_F_opt_y == 3){ // F at X% SPR
      FAA_tot_proj = get_FXSPR(M, sel_tot_proj, waassb, mat, percentSPR, fracyr_SSB_y, log_SPR0_y, F_init) * sel_tot_proj * 0.01*percentFXSPR;
      FAA_proj = FAA_tot_proj(which_F_age(y)-1) * sel_proj;
    }
    if(proj_F_opt_y == 4){ // user-specified F
      if(proj_Fcatch(y-n_years_model) < 1e-10){ // if F = 0, sel_proj is NaN
        FAA_proj.setZero();
      } else {
        FAA_proj = Type(proj_Fcatch(y-n_years_model)) * sel_proj;
        FAA_tot_proj = FAA_proj.colwise().sum();
        //FAA_tot_proj = Type(proj_Fcatch(y-n_years_model)) * sel_tot_proj;
        //FAA_proj = FAA_tot_proj(which_F_age(y)-1) * sel_proj;
      }
    }
    if(proj_F_opt_y == 5){ // calculate F from user-specified catch
      Type thecatch = proj_Fcatch(y-n_years_model);
      if(thecatch < 1e-10){ // if catch = 0, F = 0 and sel_proj is NaN
        FAA_proj.setZero();
      } else {
        //FAA_tot_proj = get_F_from_log_Catch(thecatch, NAA_y, M, sel_tot_proj, waacatch, F_init) * sel_tot_proj;
        FAA_proj = get_F_from_log_Catch(thecatch, NAA_y, M, sel_proj, waacatch, F_init) * sel_proj;
        FAA_tot_proj = FAA_proj.colwise().sum();
      }
    }
    if(proj_F_opt_y == 6){ //Fmsy
      FAA_proj = get_FMSY(log_a_y, log_b_y, M, sel_proj, waacatch, waassb, mat, fracyr_SSB_y, log_SPR0_y, recruit_model, F_init) * sel_proj * 0.01* percentFMSY;
      FAA_tot_proj = FAA_proj.colwise().sum();
      //FAA_proj = FAA_tot_proj(which_F_age(y)-1) * sel_proj;
    }
  }
  return(FAA_proj);
}


template <class Type>
matrix<Type> sim_pop(array<Type> NAA_devs, int recruit_model, vector<Type> mean_rec_pars, vector<Type> SSBin, matrix<Type> NAAin, vector<Type> log_SR_a, 
  vector<Type> log_SR_b, matrix<int> Ecov_where, vector<int> Ecov_how, array<Type> Ecov_lm, int n_NAA_sigma, 
  int do_proj, vector<int> proj_F_opt, array<Type> FAA, matrix<Type> FAA_tot, matrix<Type> MAA, matrix<Type> mature, array<Type> waa, 
  vector<int> waa_pointer_fleets, int waa_pointer_ssb, vector<Type> fracyr_SSB, vector<Type> log_SPR0, vector<int> avg_years_ind, 
  int n_years_model, int n_fleets, vector<int> which_F_age, Type percentSPR, vector<Type> proj_Fcatch, Type percentFXSPR, vector<Type> F_proj_init, Type percentFMSY){

  // Population model (get NAA, numbers-at-age, for all years)
  int ny = log_SR_a.size();
  int n_ages = NAAin.cols();
  vector<Type> SSB(ny);
  SSB(0) = SSBin(0);
  matrix<Type> NAA(ny, n_ages), pred_NAA(ny, n_ages);
  matrix<Type> ZAA = MAA + FAA_tot;
  matrix<Type> log_NAA(ny-1, n_ages);
  log_NAA.setZero();
  NAA.row(0) = pred_NAA.row(0) = NAAin.row(0);
  for(int y = 1; y < ny; y++)
  {
    //use NAA.row(y-1)
    pred_NAA.row(y) = get_pred_NAA_y(y, recruit_model, mean_rec_pars, SSB, NAA, log_SR_a, 
      log_SR_b, Ecov_where, Ecov_how, Ecov_lm, ZAA);
    
    // calculate NAA
    if(n_NAA_sigma > 1){
      // all NAA are estimated (random effects)
      for(int a = 0; a < n_ages; a++) 
      {
        log_NAA(y-1,a) = log(pred_NAA(y,a)) + NAA_devs(y-1,a);
        NAA(y,a) = exp(log_NAA(y-1,a));
      }
    } else {
      // only recruitment estimated (either fixed or random effects)
      log_NAA(y-1,0) = log(pred_NAA(y,0)) + NAA_devs(y-1,0);
      NAA(y,0) = exp(log_NAA(y-1,0));
      for(int a = 1; a < n_ages; a++) NAA(y,a) = pred_NAA(y,a); // survival is deterministic     
    }

    // calculate F and Z in projection years, here bc need NAA(y) if using F from catch
    if(do_proj == 1){ // only need FAA_tot for projections, use average FAA_tot over avg.yrs
      // get selectivity using average over avg.yrs
      if(y > n_years_model-1){
        //vector<Type> waacatch = get_waa_y(waa, y, n_ages, waa_pointer_totcatch);
        matrix<Type> waacatch = get_waacatch_y(waa, y, n_ages, waa_pointer_fleets);
        vector<Type> waassb = get_waa_y(waa, y, n_ages, waa_pointer_ssb);
        matrix<Type> FAA_proj = get_F_proj(y, n_fleets, proj_F_opt, FAA, NAA, MAA, mature, waacatch, waassb, fracyr_SSB, log_SPR0, avg_years_ind, n_years_model,
         which_F_age, percentSPR, proj_Fcatch, percentFXSPR, F_proj_init(y-n_years_model), log_SR_a, log_SR_b, recruit_model, percentFMSY);
        for(int f = 0; f < n_fleets; f++) for(int a = 0; a< n_ages; a++) FAA(y,f,a) = FAA_proj(f,a);
        FAA_tot.row(y) = FAA_proj.colwise().sum();
        ZAA.row(y) = FAA_tot.row(y) + MAA.row(y);
      }
    } // end proj F
    SSB(y) = get_SSB(NAA,ZAA,waa, mature,y, waa_pointer_ssb, fracyr_SSB);
  } // end pop model loop
  
  matrix<Type> out(ny, 2 * n_ages + n_fleets * n_ages + 1);
  for(int i = 0; i < n_ages; i++) out.col(i) = NAA.col(i);
  for(int i = n_ages; i < 2 * n_ages; i++) out.col(i) = pred_NAA.col(i-n_ages);
  for(int y = 0; y < ny; y++) for(int i = 2 ; i < (2 + n_fleets); i++) for(int a = 0; a < n_ages; a++) out(y,i*n_ages + a) = FAA(y,i-2,a);
  out.col(out.cols()-1) = SSB;
  return(out);
}

template <class Type>
array<Type> get_SPR(int y, vector<int> spawn_seasons, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, matrix<Type> fracyr_SSB, array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L, 
  array<Type> mature, array<Type> waa, vector<int> waa_pointer_ssb, vector<Type> fracyr_seasons, int age_specific, int small_dim){
  /* 
    calculate equilibrium spawning biomass per recruit (at age) by stock and region. If movement is set up approriately 
    all fish can be made to return to a single spawning region for each stock.
                  y: the model year for which to use SPR inputs
       spawn_season: vector of indicators telling which season spawning occurs for each stock
      fleet_regions: vector of indicators telling which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_years x n_stocks:  size of interval from beginning of season to time of spawning within that season
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions array of movement matrices
                  L: n_years_model x n_regions. "extra" unobserved mortality
             mature: n_stocks x n_years x n_ages proportion mature at age
                waa: (n_?) x n_years x n_ages_model. weight at age
    waa_pointer_ssb: which waa matrix to use for ssb
     fracyr_seasons: n_seasons: length of intervals for each season
       age_specific: 0/1 telling whether to return SSB/R by age or not. dimensions of returned arrays are different. 
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  int n_ages = log_M.dim(3);
  int n_fleets = FAA.dim(0);
  int n_seasons = can_move.dim(1);
  int P_dim = n_regions + n_fleets + 1;
  array<Type> SPRAA(n_stocks,n_ages,n_regions,n_regions); //SSB/R at age in each region column, given recruited in region row
  SPRAA.setZero();

  matrix<Type> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < P_dim; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<Type> S_ya = get_S(I, n_regions);
    for(int a = 0; a < n_ages-1; a++) {
      matrix<Type> P_ya = I, P_spawn = I; //PTM for year and age and up to time of spawning
      for(int t = 0; t < n_seasons; t++) {
        if(t == spawn_seasons(s)-1) {
          //P(0,t_spawn): PTM over entire year up to time of spawning
          P_spawn = P_ya * get_P_t(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(y,s), FAA, log_M, mu, L);
        }
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P_t(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      } 
      // SSB/R at year and age = prob alive to up to age a-1 x prob spawn at age a x waa x mature
      //should be n_regions x n_regions
      matrix<Type> SPR_ya = S_ya * get_S(P_spawn, n_regions) * mature(s,y,a) * waa(waa_pointer_ssb(s)-1,y,a);
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
        //SSB per Recruit in each region (should only have positive values in spawn_regions(s)-1?)
        SPRAA(s,a,i,j) += SPR_ya(i,j); 
      }
      S_ya = S_ya * get_S(P_ya, n_regions); //accumulate for next age
    }
    //now plus group
    matrix<Type> P_ya = I, P_spawn = I; //PTM for year and age and up to time of spawning
    for(int t = 0; t < n_seasons; t++) {
      if(t == spawn_seasons(s)-1) {
        //P(0,t_spawn): PTM over entire year up to time of spawning
        P_spawn = P_ya * get_P_t(n_ages-1, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(y,s), FAA, log_M, mu, L);
      }
      //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
      P_ya = P_ya * get_P_t(n_ages-1, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
    }
    matrix<Type> fundm(n_regions,n_regions);
    fundm.setZero();
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      fundm(i,j) = -P_ya(i,j);
      if(i==j) fundm(i,j) += 1;
    }
    if(small_dim) fundm = fundm.inverse(); else fundm = atomic::matinv(fundm);
    //for plus group S_ya = S_y,a-1 x (I - S_y,+)^-1
    S_ya = S_ya * fundm;
    // SSB/R at year and age = prob alive to up to age a-1 x prob spawn at age a x waa x mature
    //should be n_regions x n_regions
    matrix<Type> SPR_ya = S_ya * get_S(P_spawn, n_regions) * mature(s,y,n_ages-1) * waa(waa_pointer_ssb(s)-1,y,n_ages-1);
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      //SSB per Recruit in each region (should only have positive values in spawn_regions(s)-1?)
      SPRAA(s,n_ages-1,i,j) += SPR_ya(i,j); 
    }
  }
  if(age_specific) {
    return SPRAA; 
  } else {
    array<Type> SPR(n_stocks,n_regions,n_regions);
    SPR.setZero();
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      SPR(s,i,j) += SPRAA(s,a,i,j);
    }
    return SPR;
  }
}
