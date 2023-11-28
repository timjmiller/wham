template <class T>
vector<T> get_F_t(vector<int> fleet_season, int age, int year, array<T> FAA){
  vector<T> F_t(FAA.dim(0));
  for(int f = 0; f < FAA.dim(0); f++) if(fleet_season(f)) F_t(f) = FAA(f,year,age);
  return F_t;
}

template <class T>
vector<T> get_F_t(vector<int> fleet_season, int age, matrix<T> FAA){
  vector<T> F_t(FAA.rows());
  for(int f = 0; f < FAA.rows(); f++) if(fleet_season(f)) F_t(f) = FAA(f,age);
  return F_t;
}


template<class Type>
matrix<Type> get_log_F(matrix<Type>Fpars, int Fconfig, int n_years_pop){
  int n_y = Fpars.rows();
  int n_fleets = Fpars.cols();
  matrix<Type> log_F(n_years_pop, n_fleets);
  log_F.setZero();
  for(int f = 0; f < n_fleets; f++) {
    if(Fconfig == 1){
      log_F(0,f) = Fpars(0,f);
      for(int y = 1; y < n_y; y++) log_F(y,f) = log_F(y-1,f) + Fpars(y,f);
    }
    if(Fconfig == 2) for(int y = 0; y < n_y; y++) log_F(y,f) = Fpars(y,f);
  }
  return(log_F);
}

template<class Type>
array<Type> get_FAA(matrix<Type> log_F, vector<matrix<Type>> selAA, matrix<int> selblock_pointer, int n_ages, int n_years_model){
  int n_fleets = log_F.cols();
  int n_y = log_F.rows();
  array<Type> FAA(n_fleets,n_y,n_ages);
  FAA.setZero();
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_years_model; y++) for(int a = 0; a < n_ages; a++) {
    FAA(f,y,a) = exp(log_F(y,f)) * selAA(selblock_pointer(y,f)-1)(y,a);
  }
  return(FAA);
}

template<class Type>
array<Type> get_log_FAA(array<Type> FAA){
  array<Type> log_FAA(FAA.dim(0),FAA.dim(1),FAA.dim(2));
  log_FAA.setZero();
  for(int f = 0; f < FAA.dim(0); f++) for(int y = 0; y < FAA.dim(1); y++) for(int a = 0; a < FAA.dim(2); a++) {
    log_FAA(f,y,a) = log(FAA(f,y,a));
  }
  return(log_FAA);
}
//done

template<class Type>
array<Type> get_FAA_by_region(array<Type> FAA, vector<int> fleet_regions, int n_regions){
  int n_fleets = FAA.dim(0);
  int n_ages = FAA.dim(2);
  int n_y = FAA.dim(1);
  array<Type> FAA_tot(n_regions, n_y, n_ages);
  FAA_tot.setZero();
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_y; y++) for(int a = 0; a < n_ages; a++) {
    FAA_tot(fleet_regions(f)-1,y,a) += FAA(f,y,a);
  }
  return(FAA_tot);
}

template <class Type>
array<Type> get_pred_stock_CAA(array<Type> NAA, array<Type> annual_Ps){
  int n_stocks = NAA.dim(0);
  int n_regions = NAA.dim(1);
  int n_fleets = annual_Ps.dim(4)-n_regions-1;
  int n_years = annual_Ps.dim(1);
  int n_ages = NAA.dim(3);

  array<Type> pred_stock_CAA(n_fleets,n_stocks,n_years,n_ages);
  pred_stock_CAA.setZero();

  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) {
    for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) {
      pred_stock_CAA(f,s,y,a) +=  NAA(s,r,y,a) * annual_Ps(s,y,a,r,n_regions + f);
    }
  }
  return pred_stock_CAA;
}

template <class Type>
array<Type> get_pred_CAA(array<Type> pred_stock_CAA){
  int n_fleets = pred_stock_CAA.dim(0);
  int n_stocks = pred_stock_CAA.dim(1);
  int n_years = pred_stock_CAA.dim(2);
  int n_ages = pred_stock_CAA.dim(3);

  array<Type> pred_CAA(n_fleets,n_years,n_ages);
  pred_CAA.setZero(); //(n_years,n_stocks,n_ages,n_regions);
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) {
    for(int s = 0; s < n_stocks; s++) pred_CAA(f,y,a) += pred_stock_CAA(f,s,y,a); 
  }
  return pred_CAA;
}

template <class Type>
array<Type> get_pred_stock_catch(array<Type> pred_stock_CAA, array<Type> waa_catch){
  int n_fleets = pred_stock_CAA.dim(0);
  int n_stocks = pred_stock_CAA.dim(1);
  int n_years = pred_stock_CAA.dim(2);
  int n_ages = pred_stock_CAA.dim(3);

  array<Type> pred_stock_catch(n_fleets,n_stocks,n_years);
  pred_stock_catch.setZero();

  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) {
    for(int s = 0; s < n_stocks; s++) {
        pred_stock_catch(f,s,y) +=  pred_stock_CAA(f,s,y,a) *  waa_catch(f,y,a);
    }
  }
  return pred_stock_catch;
}

template <class Type>
matrix<Type> get_pred_catch(array<Type> pred_stock_catch){
  int n_fleets = pred_stock_catch.dim(0);
  int n_stocks = pred_stock_catch.dim(1);
  int n_years = pred_stock_catch.dim(2);

  matrix<Type> pred_catch(n_years,n_fleets);
  pred_catch.setZero();
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_years; y++) for(int s = 0; s < n_stocks; s++){
    pred_catch(y,f) += pred_stock_catch(f,s,y);
  }
  return pred_catch;
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
array<Type> get_pred_catch_paa(array<Type> pred_CAA, int n_years_model){
  int n_fleets = pred_CAA.dim(0);
  int n_y = pred_CAA.dim(1);
  int n_ages = pred_CAA.dim(2);
  array<Type> pred_catch_paa(n_fleets,n_y, n_ages);

  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_years_model; y++){
    Type tsum = 0.0;
    for(int a = 0; a < n_ages; a++){
      tsum += pred_CAA(f,y,a);
    }
    for(int a = 0; a < n_ages; a++){
      pred_catch_paa(f,y,a) = pred_CAA(f,y,a)/tsum;
    }
  }
  return pred_catch_paa;
}

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
