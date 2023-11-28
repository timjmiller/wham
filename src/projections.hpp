/* calculate total catch at F across multiple fleets*/
template<class Type>
struct log_catch_fleets_F_multi {
  /* Data and parameter objects for calculation: */
  array<Type> NAA;
  array<Type> log_M;
  array<Type> mu;
  vector<Type> L;
  array<Type> sel; //n_fleets x n_ages
  vector<Type> fracyr_season;
  vector<int> fleet_regions;
  matrix<int> fleet_seasons;
  array<int> can_move;
  vector<int> mig_type;
  array<Type> waacatch; //n_fleets x n_ages
  int trace;

  /* Constructor */
  log_catch_fleets_F_multi(
  array<Type> NAA_,
  array<Type> log_M_,
  array<Type> mu_,
  vector<Type> L_,
  array<Type> sel_,
  vector<Type> fracyr_season_,
  vector<int> fleet_regions_,
  matrix<int> fleet_seasons_,
  array<int> can_move_,
  vector<int> mig_type_,
  array<Type> waacatch_,
  int trace_) :
    NAA(NAA_), log_M(log_M_), mu(mu_), L(L_), sel(sel_), fracyr_season(fracyr_season_), 
    fleet_regions(fleet_regions_), fleet_seasons(fleet_seasons_), can_move(can_move_), mig_type(mig_type_), 
    waacatch(waacatch_), trace(trace_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it achieves required catch
    int n_stocks = log_M.dim(0);
    int n_regions = log_M.dim(1);
    int n_seasons = fracyr_season.size();
    int n_ages = log_M.dim(2);
    int n_fleets = waacatch.dim(0);
    int Pdim = n_regions + n_fleets + 1;
  
    if(trace) see(sel);

    matrix<T> FAA_T(n_fleets,n_ages);
    for(int f = 0; f< n_fleets; f++) for(int a = 0; a < n_ages; a++) FAA_T(f,a) = exp(log_F(0)) * T(sel(f,a));
    if(trace) see(FAA_T);
    array<T> logM_T(log_M.dim(0),log_M.dim(1),log_M.dim(2));
    for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) for(int a = 0; a < n_ages; a++){
      logM_T(s,r,a) = T(log_M(s,r,a));
    }
    if(trace) see(logM_T);
    array<T> mu_T(mu.dim(0),mu.dim(1),mu.dim(2),mu.dim(3),mu.dim(4));
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int t = 0; t < mu.dim(2); t++) {
      for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) {
        mu_T(s,a,t,r,rr) = T(mu(s,a,t,r,rr));
      }
    }
    if(trace) see(mu_T);
    vector<T> L_T = L.template cast<T>();
    if(trace) see(L_T);
    //vector<T> fracyrseason_T = fracyr_season.template cast<T>();
    //if(trace) see(fracyrseason_T);

    matrix<T> catch_stock_fleet(n_stocks,n_fleets);
    catch_stock_fleet.setZero();
    matrix<T> I(Pdim,Pdim);
    I.setZero();
    for(int i = 0; i < Pdim; i++) I(i,i) = 1.0;

    for(int s = 0; s < n_stocks; s++) {
      for(int a = 0; a < n_ages; a++) {
        if(trace) see(a);
        matrix<T> P_ya = I;
        for(int t = 0; t < n_seasons; t++) {
          if(trace) see(t);
          //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
          P_ya = P_ya * get_P_t(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, T(fracyr_season(t)), FAA_T, logM_T, mu_T, L_T, trace);
          if(trace) see(P_ya);
        }
        for(int f = 0; f < n_fleets; f++) for(int r = 0; r < n_regions; r++) {
          if(trace) {
            see(f);
            see(r);
            see(a);
            see(NAA.dim);
            see(waacatch.dim(0));
            see(waacatch.dim(1));
          }
          catch_stock_fleet(s,f) +=  T(NAA(s,r,a)) * P_ya(r,n_regions + f) * T(waacatch(f,a));
        }
      }
    }
    T Catch = catch_stock_fleet.sum();
    return log(Catch);
  }
};


//multiple fleets, regions,stocks
template <class Type>
Type get_F_from_Catch(Type Catch, array<Type> NAA, array<Type> log_M, array<Type> mu, vector<Type> L, array<Type> sel,
  vector<Type> fracyr_season, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type,
  array<Type> waacatch, int trace, Type F_init)
{
  int n = 10;
  vector<Type> log_F_i(1);
  vector<Type> log_F_iter(n);
  log_F_iter.fill(log(F_init)); //starting value
  if(trace) see(log_F_iter);
  if(trace) see(Catch);
  log_catch_fleets_F_multi<Type> logcatch_at_F(NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, 
    waacatch, trace);
  if(trace) see("past log_catch_fleets_F_multi");
  for (int i=0; i<n-1; i++) {
    log_F_i(0) = log_F_iter(i);
    if(trace) see(i);
    vector<Type> grad_log_catch_at_F = autodiff::gradient(logcatch_at_F,log_F_i);
    if(trace) see(grad_log_catch_at_F);
    if(trace) see("after autodiff::gradient");
    log_F_iter(i+1) = log_F_iter(i) - (logcatch_at_F(log_F_i) - log(Catch))/grad_log_catch_at_F(0); //uses log(catch) for potentially slower changes in derivatives?
    if(trace) see(log_F_iter);
  }
  Type res = exp(log_F_iter(n-1));
  return res;
}


template <class Type>
array<Type> update_FAA_proj(int y, vector<int> proj_F_opt, array<Type> FAA, array<Type> NAA, array<Type> log_M, array<Type> mu,
  matrix<Type> L, array<Type> mature_proj, array<Type> waa_ssb_proj, array<Type> waa_catch_proj, vector<int> fleet_regions, matrix<int> fleet_seasons, 
  vector<Type> fracyr_SSB_proj, vector<int> spawn_regions, array<int> can_move, array<int> must_move, vector<int> mig_type, 
  vector<int> avg_years_ind, int n_years_model, vector<int> which_F_age, vector<Type> fracyr_seasons, int small_dim,
  Type percentSPR, vector<Type> proj_Fcatch, Type percentFXSPR, Type percentFMSY, matrix<Type> R_XSPR, vector<Type> FXSPR_init, 
  vector<Type> FMSY_init, vector<Type> F_proj_init, matrix<Type> log_a, matrix<Type> log_b, vector<int> spawn_seasons, vector<int> recruit_model, 
  vector<Type> SPR_weights, int SPR_weight_type, int trace){
    /* 
     update FAA in projection year y
                   y:  year of projection (>n_years_model)
          proj_F_opt:  for each projection year, how to specify F for projection. 1: use terminal FAA, 2: use average FAA (avg_years_ind), 
                          3: F at X%SPR, 4: user-specified full-F, 5: user-specified catch, 6: use Fmsy (inputs averaged over avg_years_ind))
                 FAA:  FAA array from main code.
                 NAA:  NAA array from main code
               log_M:  array from main code.
       mu:  array from main code.
              mature_proj:  maturity matrix (n_stocks, n_ages). year y
                 waa_ssb_proj:  weight at age (n_stocks,n_ages) year y
                 waa_catch_proj: weigth at age (n_fleets,n_ages); year y
              fracyr_SSB_proj:  (n_stocks) vector of yearly fractions of the year when spawning occurs. year y
          spawn_regions:
          can_move:
          mig_type:

       avg_years_ind:  integer vector of years to average F for projection
       n_years_model:  number of years before projection begins
         which_F_age:  age to define which age has max F
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
  int n_ages = FAA.dim(2);
  //int n_regions = NAA.dim(1);
  //int n_stocks = NAA.dim(0);
  //int n_seasons = fracyr_seasons.size();

  //int n_toavg = avg_years_ind.size();
  int proj_F_opt_y = proj_F_opt(y-n_years_model);
  if(trace) see(proj_F_opt_y);

  //proj_F_opt == 1, last year F (default)
  array<Type> FAA_proj(n_fleets, n_ages);
  //matrix<Type> FAA_tot_proj(n_regions, n_ages);
  // FAA_tot_proj.setZero();
  FAA_proj.setZero();
  if(proj_F_opt_y == 1){ // last year F (default)
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) {
      FAA_proj(f,a) = FAA(f,n_years_model-1,a); 
    }
    if(trace) see(FAA_proj);
  }
  else { //proj_F_opt_y>1
    //option 2: average F over years defined in avg_years_ind
    if(proj_F_opt_y == 2) { 
      FAA_proj = get_avg_FAA_as_array(FAA,avg_years_ind,0);
      if(trace) see(FAA_proj);
    } else { //proj_F_opt_y > 2

      //need selectivity for projections for rest of options
      array<Type> sel_proj = get_avg_fleet_sel_as_array(FAA, avg_years_ind, which_F_age(y));
      if(trace) see(sel_proj);
      Type Fproj = 0;
      
      if(proj_F_opt_y == 4){ // user-specified F
        if(proj_Fcatch(y-n_years_model) < 1e-10){ // if F = 0, sel_proj is NaN
          // FAA_proj.setZero();
        } else {
          Fproj = Type(proj_Fcatch(y-n_years_model));
          // for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++){ 
          //   FAA_proj(f,a) = sel_proj(f,a) * Type(proj_Fcatch(y-n_years_model));
          // }
        }
        if(trace) see(proj_Fcatch(y-n_years_model));
        // if(trace) see(FAA_proj);
      }
       
      if((proj_F_opt_y == 3) | (proj_F_opt_y == 5) | (proj_F_opt_y == 6) ){

        //These all have been calculated in projection years already so just need to extract for this year
        if(trace) see(L.rows());
        if(trace) see(L.cols());
        vector<Type> L_proj = L.row(y);
        if(trace) see(L_proj);
        array<Type> log_M_proj = get_log_M_y(y, log_M);
        if(trace) see(log_M_proj.dim);
        array<Type> mu_proj = get_mu_y(y, mu);
        if(trace) see(mu_proj.dim);
        if(trace) see(waa_ssb_proj);
        if(trace) see(waa_catch_proj);
        if(trace) see(mature_proj);
        if(trace) see(fracyr_SSB_proj);
        if(trace) see(R_XSPR.row(y));

        if(proj_F_opt_y == 3) {//option 3: use F X%SPR
          vector<Type> FXSPR = get_FXSPR(spawn_seasons, spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB_proj, sel_proj, 
            log_M_proj, mu_proj, L_proj, mature_proj,  waa_ssb_proj, fracyr_seasons, vector<Type> (R_XSPR.row(y)), percentSPR, SPR_weights, 
            SPR_weight_type, small_dim, FXSPR_init(y), 10, trace);
          if(trace) see(FXSPR);
          Fproj = FXSPR(0);
          if(trace) see(FXSPR(0));
          // for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++){ 
          //   FAA_proj(f,a) = sel_proj(f,a) * FXSPR(0);
          // }
          // if(trace) see(FAA_proj);
        }
        
        if(proj_F_opt_y == 5){ // calculate F from user-specified catch
          Type thecatch = proj_Fcatch(y-n_years_model);
          if(trace) see(thecatch);
          if(thecatch < 1e-10){ // if catch = 0, F = 0 and sel_proj is NaN
            // FAA_proj.setZero();
          } else {
            array<Type> NAA_y = get_NAA_y(y, NAA);
            Type F_from_Catch = get_F_from_Catch(thecatch, NAA_y, log_M_proj, mu_proj, L_proj, sel_proj, fracyr_seasons, fleet_regions, 
            fleet_seasons, can_move, mig_type, waa_catch_proj, trace, F_proj_init(y- n_years_model));
            // FAA_proj = sel_proj * get_F_from_Catch(thecatch, NAA_y, log_M_proj, mu_proj, L_proj, sel_proj, fracyr_seasons, fleet_regions, 
            // fleet_seasons, can_move, mig_type, waa_catch_proj, trace, F_proj_init(y- n_years_model));
            // for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++){ 
            //   FAA_proj(f,a) = sel_proj(f,a) * F_from_Catch;
            // }
            if(trace) see(F_from_Catch);
            Fproj = F_from_Catch;
          }
          // if(trace) see(FAA_proj);
        }

        //option 6: use FMSY 
        if(proj_F_opt_y == 6){
          vector<Type> a_proj = exp(vector<Type> (log_a.row(y)));
          vector<Type> b_proj = exp(vector<Type> (log_b.row(y)));
          Type FMSY = get_FMSY(a_proj, b_proj, spawn_seasons, spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, 
            fracyr_SSB_proj, sel_proj, log_M_proj, mu_proj, L_proj, mature_proj,  waa_ssb_proj, waa_catch_proj, fracyr_seasons, recruit_model, small_dim, 
            FMSY_init(y), 10, trace);
          if(trace) see(FMSY);
          // FAA_proj = sel_proj * get_FMSY(a_proj, b_proj, spawn_seasons, spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, 
          //   fracyr_SSB_proj, sel_proj, log_M_proj, mu_proj, L_proj, mature_proj,  waa_ssb_proj, waa_catch_proj, fracyr_seasons, recruit_model, small_dim, 
          //   FMSY_init(y), 10, trace);
          Fproj = FMSY;
            // for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++){ 
            //   FAA_proj(f,a) = sel_proj(f,a) * FMSY;
            // }
          // if(trace) see(FAA_proj);
        //F_full is the same as that used to generate selectivity to project
        }
      }
      for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) FAA_proj(f,a) = sel_proj(f,a) * Fproj;
    }
  }
  array<Type> updated_FAA = FAA;
  if(trace) see(FAA.dim);
  if(trace) see(y);
  for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) updated_FAA(f,y,a) = FAA_proj(f,a);
  
  return updated_FAA;
}
