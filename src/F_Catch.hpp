template<class Type>
array<Type> get_FAA(matrix<Type>F, matrix<int> fleet_seasons, vector<matrix<Type>> selAA, matrix<int> selblock_pointer, int n_ages, int n_seasons){
  int n_fleets = F.cols();
  int n_y = F.rows();
  array<Type> FAA(n_fleets,n_y,n_seasons,n_ages);
  FAA.setZero();
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_y; y++) for(int a = 0; a < n_ages; a++) 
  {
    for(int t = 0; t < n_seasons; t++) if(fleet_seasons(f,t) == 1) 
    {
      FAA(f,y,t,a) = F(y,f) * selAA(selblock_pointer(y,f)-1)(y,a);
    }
  }
  return(FAA);
}
//done

template<class Type>
array<Type> get_FAA_tot(array<Type> FAA, vector<int> fleet_regions, matrix<int> fleet_seasons, int n_regions){
  int n_fleets = FAA.dim(0);
  int n_seasons = FAA.dim(2);
  int n_ages = FAA.dim(3);
  int n_y = FAA.dim(1);
  array<Type> FAA_tot(n_regions, n_y, n_seasons, n_ages);
  FAA_tot.setZero();
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_y; y++) for(int a = 0; a < n_ages; a++) 
  {
    for(int t = 0; t < n_seasons; t++) if(fleet_seasons(f,t) == 1) 
    {
      FAA_tot(fleet_regions(f)-1,y,t,a) += FAA(f,y,t,a);
    }
  }
  return(FAA_tot);
}
//done

template <class Type>
array<Type> get_sel_proj(int y, array<Type> FAA, vector<int> avg_years_ind,
  vector<int> which_F_fleet, vector<int> which_F_season, vector<int> which_F_age){
    /* 
     get selectivity to project for next time step
                   y:  year of projection (>n_years_model)
                 FAA:  FAA array from main code.
       avg_years_ind:  integer vector of years to average F for projection
       which_F_fleet:  define which fleet has max F
      which_F_season:  define which season has max F
         which_F_age:  define which age has max F
    */
  //average F by fleet, season, and age is used to find selectivity (fleet,season,age) to project 
  //full F is the FAA for fleet, season and age defined by which_F_fleet,which_F_season, which_F_age
  int n_toavg = avg_years_ind.size();
  int n_fleets = FAA.dim(0);
  int n_seasons = FAA.dim(2);
  int n_ages = FAA.dim(3);

  array<Type> FAA_avg(n_fleets, n_seasons, n_ages);
  FAA_avg.setZero();
  for(int f = 0; f < n_fleets; f++) for(int t = 0; t < n_seasons; t++)
  {
    for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_toavg; i++){
      FAA_avg(f,t,a) += FAA(f,avg_years_ind(i),t,a)/Type(n_toavg);
    }
  }
  //get selectivity using average over avg.yrs
  array<Type> sel_proj(n_fleets,n_seasons,n_ages);
  //fully selected F across regions, seasons, and ages
  Type F_full = FAA_avg(which_F_fleet(y)-1,which_F_season(y)-1,which_F_age(y)-1);
  for(int f = 0; f < n_fleets; f++) for(int t = 0; t < n_seasons; t++){
    for(int a = 0; a < n_ages; a++) {
      sel_proj(f,t,a) = FAA_avg(f,t,a)/F_full;
    }
  }
  return(sel_proj);
}
//done

template <class Type>
array<Type> get_FAA_proj(int y, vector<int> proj_F_opt, array<Type> sel_proj, 
  array<Type> FAA, array<Type> NAA, array<Type> MAA, 
  array<Type> mature, array<Type> waa, vector<int> waa_pointer_totcatch, vector<int> waa_pointer_ssb, 
  matrix<Type> fracyr_SSB, matrix<Type> log_SPR0, vector<int> avg_years_ind, 
  int n_years_model, vector<int> which_F_fleet, vector<int> which_F_season, vector<int> which_F_age, 
  Type percentSPR, vector<Type> proj_Fcatch, Type percentFXSPR, Type F_init, 
  matrix<Type> log_a, matrix<Type> log_b, vector<int> recruit_model, Type percentFMSY){
    /* 
     get FAA to project for next time step
                   y:  year of projection (>n_years_model)
          proj_F_opt:  for each year, how to specify F for projection (1 to 6)
                 FAA:  FAA array from main code.
                 NAA:  NAA array from main code
                 MAA:  MAA array from main code.
              mature:  maturity array from main code.
                 waa:  weight at age array 
waa_pointer_totcatch:  (n_regions) pointer for waa to use for tot catch (use function get_waa_y defined in helper.hpp)
     waa_pointer_ssb:  (n_stocks) pointer for waa to use for SSB (use function get_waa_y defined in helper.hpp)
          fracyr_SSB:  (n_stocks x n_years) vector of yearly fractions of the year when spawning occurs
            log_SPR0:  matrix (n_stocks x n_years) of yearly log(unfished SSB/R) 
       avg_years_ind:  integer vector of years to average F for projection
       n_years_model:  number of years before projection begins
       which_F_fleet:  define which fleet has max F
      which_F_season:  define which season has max F
         which_F_age:  define which age has max F
          percentSPR:  percentage (0-100) of unfished spawning potential to determine F_percentSPR
         proj_Fcatch:  vector (n_years_proj) of user specified Fishing mortality rates to project
        percentFXSPR:  percentage (0-100) of F_percentSPR to use in catch, e.g. GOM cod uses F = 75% F_40%SPR
              F_init:  initial value to use for FXSPR or FMSY newton method
               log_a:  (n_stocks x n_years) annual log(a) for stock-recruit relationship
               log_b:  (n_stocks x n_years) annual log(b) for stock-recruit relationship
       recruit_model:  (n_stocks) integer for which type of recruit model is assumed (= 3 or 4 for using Fmsy)
         percentFMSY:  percentage (0-100) of FMSY to use in catch.
    */
  int n_fleets = FAA.dim(0);
  int n_seasons = FAA.dim(2);
  int n_ages = FAA.dim(3);
  int n_regions = NAA.dim(1);
  matrix<Type> waacatch = get_waa_y(waa, y, n_ages, waa_pointer_totcatch);
  matrix<Type> waassb = get_waa_y(waa, y, n_ages, waa_pointer_ssb);

  int n_toavg = avg_years_ind.size();
  int proj_F_opt_y = proj_F_opt(y-n_years_model);

  //proj_F_opt == 1, last year F (default)
  array<Type> FAA_proj(n_fleets, n_seasons, n_ages);
  array<Type> FAA_tot_proj(n_regions,n_seasons, n_ages);
  FAA_tot_proj.setZero();
  FAA_proj.setZero();
  if(proj_F_opt_y == 1){ // last year F (default)
    for(int f = 0; f < n_fleets; f++) for(int t = 0; t < n_seasons; t++) for(int a = 0; a < n_ages; a++) {
      FAA_proj(f,t,a) = FAA(f,n_years_model-1,t,a); 
    }
  }
  else { //proj_F_opt_y>1
    //option 2: average F is by fleet and Ftot is sum of fleet averages
    //when there is more than 1 fleet, the sum of predicted catch across fleets will not generally equal the total catch using FAA_tot and waa_totcatch.
    Type F_full_proj = 0.0;
    if(proj_F_opt_y == 2) { //F_full is the same as that used to generate selectivity to project
      F_full_proj = FAA(which_F_fleet(y)-1, which_F_season(y)-1,which_F_age(y)-1); 
    }

    if(proj_F_opt_y == 4){ // user-specified F
      /*if(proj_Fcatch(y-n_years_model) < 1e-10){ // if F = 0, sel_proj is NaN
        FAA_proj.setZero();
      } else { */
        F_full_proj = Type(proj_Fcatch(y-n_years_model));
      //}
    }
    
    array<Type> sel_proj = get_sel_proj(y, FAA, avg_years_ind, which_F_fleet, which_F_season, which_F_age);
    /*
    if(proj_F_opt_y == 3){ // F at X% SPR
      F_full_proj = get_FXSPR(MAA_y, sel_proj, waassb, mat_y, percentSPR, fracyr_SSB_y, log_SPR0_y, F_init) * 0.01* percentFXSPR;
    }
    if(proj_F_opt_y == 5){ // calculate F from user-specified catch
      Type thecatch = proj_Fcatch(y-n_years_model);
      if(thecatch < 1e-10){ // if catch = 0, F = 0 and sel_proj is NaN
        //F_full_proj = 0.0; //already done
        //FAA_proj.setZero();
      } else {
        F_full_proj = get_F_from_log_Catch(thecatch, NAA_y, MAA_y, sel_proj, waacatch, F_init);
      }
    }
    /*if(proj_F_opt_y == 6){ //Fmsy
      vector<Type> log_a_y = log_a.row(y);
      vector<Type> log_b_y = log_b.row(y);
      F_full_proj = get_FMSY(log_a_y, log_b_y, MAA_y, sel_proj, waacatch, waassb, mat_y, fracyr_SSB_y, log_SPR0_y, recruit_model, F_init) * 0.01* percentFMSY;
    } */
    FAA_proj = F_full_proj * sel_proj;
  }
  return(FAA_proj);
}

template <class Type>
array<Type> get_pred_stock_CAA(array<Type> NAA, array<Type> annual_Ps){
  int n_stocks = NAA.dim(0);
  int n_regions = NAA.dim(1);
  int n_fleets = annual_Ps.dim(4)-n_regions-1;
  int n_years = NAA.dim(2);
  int n_ages = NAA.dim(3);

  array<Type> pred_stock_CAA(n_fleets,n_stocks,n_years,n_ages);
  pred_stock_CAA.setZero();

  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) {
    for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) {
      pred_stock_CAA(f,s,y,a) +=  NAA(s,y,a,r) * annual_Ps(s,y,a,r,n_regions + f);
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
array<Type> get_pred_stock_catch(array<Type> pred_stock_CAA, array<Type> waa, vector<int> waa_pointer_fleets){
  int n_fleets = pred_stock_CAA.size();
  int n_stocks = pred_stock_CAA.dim(0);
  int n_years = pred_stock_CAA.dim(2);
  int n_ages = pred_stock_CAA.dim(3);

  array<Type> pred_stock_catch(n_fleets,n_stocks,n_years);
  pred_stock_catch.setZero();

  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) {
    for(int s = 0; s < n_stocks; s++) {
        pred_stock_catch(f,s,y) +=  pred_stock_CAA(f,s,y,a) *  waa(waa_pointer_fleets(f)-1,y,a);
    }
  }
  return pred_stock_catch;
}

template <class Type>
matrix<Type> get_pred_catch(array<Type> pred_stock_catch){
  int n_fleets = pred_stock_catch.size();
  int n_stocks = pred_stock_catch.dim(0);
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
array<Type> get_pred_catch_paa(array<Type> pred_CAA){
  int n_fleets = pred_CAA.dim(0);
  int n_y = pred_CAA.dim(1);
  int n_ages = pred_CAA.dim(2);
  array<Type> pred_catch_paa(n_fleets,n_y, n_ages);

  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_y; y++){
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
  int n_y = pred_catch_paa.dim(1);
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
  array<int> keep_Cpaa, matrix<Type> catch_Neff, vector<int> age_comp_model_fleets, matrix<Type> catch_paa_pars){
  int n_fleets = pred_catch_paa.dim(0);
  int n_y = pred_catch_paa.dim(1);
  int n_ages = pred_catch_paa.dim(2);
  vector<Type> obsvec_out = obsvec;
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_y; y++) if(use_catch_paa(y,f)) {
    vector<Type> t_pred_paa(n_ages);
    for(int a = 0; a < n_ages; a++) t_pred_paa(a) = pred_catch_paa(f,y,a);
    vector<int> ages_obs_y = agesvec.segment(keep_Cpaa(f,y,0), keep_Cpaa(f,y,1));
    vector<Type> tf_paa_obs = sim_acomp(t_pred_paa, catch_Neff(y,f), ages_obs_y, age_comp_model_fleets(f), 
      vector<Type>(catch_paa_pars.row(f)));
    obsvec_out.segment(keep_Cpaa(f,y,0),keep_Cpaa(f,y,1)) = tf_paa_obs;
  }
  return obsvec_out;
}

template <class Type>
array<Type> sim_obsvec_in_catch_paa(vector<Type> obsvec, vector<int> agesvec, array<Type> catch_paa, matrix<int> use_catch_paa, array<int> keep_Cpaa, 
  vector<int> age_comp_model_fleets){
  int n_fleets = catch_paa.dim(0);
  int n_y = catch_paa.dim(1);
  int n_ages = catch_paa.dim(2);
  array<Type> catch_paa_out = catch_paa;
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_y; y++) if(use_catch_paa(y,f)) {
    vector<Type> tf_paa_obs = obsvec.segment(keep_Cpaa(f,y,0),keep_Cpaa(f,y,1));
    vector<int> ages_obs_y = agesvec.segment(keep_Cpaa(f,y,0), keep_Cpaa(f,y,1));
    vector<Type> paa_obs_y = make_paa(tf_paa_obs, age_comp_model_fleets(f), ages_obs_y, paa_obs_y);
    for(int a = 0; a < n_ages; a++) catch_paa_out(f,y,a) = paa_obs_y(a);
  }
  return catch_paa_out;
}
