
template <class T>
array<T> get_SPR(vector<int> spawn_seasons, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, vector<T> fracyr_SSB, matrix<T> FAA, array<T> log_M, array<T> mu, vector<T> L, 
  matrix<T> mature, matrix<T> waa_ssb, vector<T> fracyr_seasons, int age_specific, int small_dim, int trace = 0){
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
  int n_ages = log_M.dim(2);
  int n_fleets = FAA.rows();
  int n_seasons = can_move.dim(1);
  int P_dim = n_regions + n_fleets + 1;
  if(trace) see("inside array<T> get_SPR");
  array<T> SPRAA(n_stocks,n_ages,n_regions,n_regions); //SSB/R at age in each region column, given recruited in region row
  SPRAA.setZero();
  if(trace) see("T1");
  matrix<T> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < P_dim; i++) I(i,i) = 1.0;
  if(trace) see(I);
  for(int s = 0; s < n_stocks; s++) {
    matrix<T> S_ya = get_S(I, n_regions);
    for(int a = 0; a < n_ages; a++) {
      if(trace) see(a);
      matrix<T> P_ya = I, P_spawn = I; //PTM for age and up to time of spawning
      vector<T> M_t(n_regions);
      for(int r = 0; r < n_regions; r++) M_t(r) = exp(log_M(s,r,a));
      if(trace) see(M_t);
      for(int t = 0; t < n_seasons; t++) {
        if(trace) see(t);
        if(trace) see(spawn_seasons(s)-1);
        vector<T> F_t(n_fleets);
        F_t.setZero();
        for(int f = 0; f < n_fleets; f++) if(fleet_seasons(f,t)) F_t(f) = FAA(f,a);
        if(trace) see(F_t);
        matrix<T> mu_t(n_regions,n_regions);
        mu_t.setZero();
        matrix<int> can_move_t(n_regions, n_regions);
        can_move_t.setZero();
        for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) {
          can_move_t(r,rr) = can_move(s,t,r,rr);
          mu_t(r,rr) = mu(s,a,t,r,rr);
        }
        if(trace) see(can_move_t);
        if(trace) see(mu_t);
        if(t == spawn_seasons(s)-1) {
          if(trace) see(P_ya);
          //P(0,t_spawn): PTM over entire year up to time of spawning
          P_spawn = P_ya * get_P_t_base(fleet_regions, can_move_t, mig_type(s), fracyr_SSB(s), F_t, M_t, 
            mu_t, L, trace);
          if(trace) see(P_spawn);
        }
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P_t_base(fleet_regions, can_move_t, mig_type(s), fracyr_seasons(t), F_t, M_t, 
          mu_t, L, trace);
        if(trace) see(P_ya);
      }
      if(a == n_ages-1){ //plus group
        matrix<T> fundm(n_regions,n_regions);
        fundm.setZero();
        for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
          fundm(i,j) = -P_ya(i,j);
          if(i==j) fundm(i,j) += 1;
        }
        if(small_dim) fundm = fundm.inverse(); else fundm = atomic::matinv(fundm);
        if(trace) see("T4");
        //for plus group S_ya = S_y,a-1 x (I - S_y,+)^-1
        S_ya = S_ya * fundm;
      }
      // SSB/R at year and age = prob alive to up to age a-1 x prob spawn at age a x waa x mature
      //should be n_regions x n_regions
      matrix<T> SPR_ya = S_ya * get_S(P_spawn, n_regions) * mature(s,a) * waa_ssb(s,a);
      if(trace) see(SPR_ya);
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
        //SSB per Recruit in each region (should only have positive values in spawn_regions(s)-1?)
        SPRAA(s,a,i,j) += SPR_ya(i,j); 
        if(trace) see(SPRAA);
        if(trace) {
          see(n_regions);
          see(a);
          see(s);
          see(i);
          see(j);
          see(SPRAA.dim);
          see(SPR_ya.rows());
          see(SPR_ya.cols());
        }
      }
      if(trace) see("out");
      if(trace) see(S_ya);
      S_ya = S_ya * get_S(P_ya, n_regions); //accumulate for next age
      if(trace) see(S_ya);
    }
    // if(trace) see("T2");
    // //now plus group
    // matrix<T> P_ya = I, P_spawn = I; //PTM for year and age and up to time of spawning
    // for(int t = 0; t < n_seasons; t++) {
    //   if(t == spawn_seasons(s)-1) {
    //     //P(0,t_spawn): PTM over entire year up to time of spawning
    //     P_spawn = P_ya * get_P(n_ages-1, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(s), FAA, log_M, mu, L, trace);
    //   }
    //   //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
    //   P_ya = P_ya * get_P(n_ages-1, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L, trace);
    // }
    // if(trace) see("T3");
    // matrix<T> fundm(n_regions,n_regions);
    // fundm.setZero();
    // for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
    //   fundm(i,j) = -P_ya(i,j);
    //   if(i==j) fundm(i,j) += 1;
    // }
    // if(small_dim) fundm = fundm.inverse(); else fundm = atomic::matinv(fundm);
    // if(trace) see("T4");
    // //for plus group S_ya = S_y,a-1 x (I - S_y,+)^-1
    // S_ya = S_ya * fundm;
    // // SSB/R at year and age = prob alive to up to age a-1 x prob spawn at age a x waa x mature
    // //should be n_regions x n_regions
    // matrix<T> SPR_ya = S_ya * get_S(P_spawn, n_regions) * mature(s,n_ages-1) * waa_ssb(s,n_ages-1);
    // for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
    //   //SSB per Recruit in each region (should only have positive values in spawn_regions(s)-1?)
    //   SPRAA(s,n_ages-1,i,j) += SPR_ya(i,j); 
    // }
  }
  if(trace) see(SPRAA);
  if(age_specific) {
    return SPRAA; 
  } else {//sum over ages
    array<T> SPR(n_stocks,n_regions,n_regions);
    SPR.setZero();
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      SPR(s,i,j) += SPRAA(s,a,i,j);
    }
    return SPR;
  }
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

template <class T>
array<T> get_YPR(vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, matrix<T> FAA, array<T> log_M, array<T> mu, vector<T> L, matrix<T> waacatch, 
  vector<T> fracyr_seasons, int age_specific, int small_dim){
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
  array<T> YPRAA(n_stocks,n_ages,n_regions,n_fleets); //Yield/R at age in each fleet column, given recruited in region row
  YPRAA.setZero();
  matrix<T> W(n_fleets,n_fleets);
  W.setZero();

  matrix<T> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < P_dim; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<T> S_ya = get_S(I, n_regions);
    for(int a = 0; a < n_ages-1; a++) {
      for(int i = 0; i < n_fleets; i++) W(i,i) = waacatch(i,a);
      matrix<T> P_ya = I; //PTM for year and age
      for(int t = 0; t < n_seasons; t++) {
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P_t(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      }
      // YPR at year and age = prob alive to up to age a-1 x prob caught at age a x waa
      matrix<T> YPR_ya = S_ya * get_D(P_ya, n_regions,n_fleets) * W; //should be n_regions x n_fleets
      S_ya = S_ya * get_S(P_ya, n_regions);  //accumulate for next age
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) YPRAA(s,a,i,j) = YPR_ya(i,j);
    }
    //now plus group
    for(int i = 0; i < n_fleets; i++) W(i,i) = waacatch(i,n_ages-1);
    matrix<T> P_ya = I; //PTM for year and age
    for(int t = 0; t < n_seasons; t++) {
      //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
      P_ya = P_ya * get_P_t(n_ages-1, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
    }
    matrix<T> fundm(n_regions,n_regions);
    fundm.setZero();
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      fundm(i,j) = -P_ya(i,j);
      if(i==j) fundm(i,j) += 1; //I - S_y
    }
    if(small_dim) fundm = fundm.inverse(); else fundm = atomic::matinv(fundm);
    //for plus group S_ya = S_y,a-1 x (I - S_y,+)^-1
    S_ya = S_ya * fundm;
    // YPR at year and age = (prob alive to up to age a-1 + prob alive at older ages) x prob caught at these ages x waa
    matrix<T> YPR_ya = S_ya * get_D(P_ya,n_regions,n_fleets) * W; //should be n_regions x n_fleets
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) YPRAA(s,n_ages-1,i,j) = YPR_ya(i,j);
  }
  if(age_specific) {
    return YPRAA; 
  } else { //sum across age
    array<T> YPR(n_stocks,n_regions,n_fleets);
    YPR.setZero();
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) {
      YPR(s,i,j) += YPRAA(s,a,i,j);
    }
    return YPR;
  }
}

template <class Type>
matrix<Type> get_RXSPR(array<Type> all_NAA, vector<int> spawn_regions, int n_years_model, int n_years_proj, 
  int XSPR_R_opt, vector<int> XSPR_R_avg_yrs){
  //R_XSPR is needed for projections and reference points
  array<Type> NAA = extract_NAA(all_NAA);
  array<Type> pred_NAA = extract_pred_NAA(all_NAA);
  int n_stocks = NAA.dim(0);
  matrix<Type> R_XSPR(n_years_model+n_years_proj, n_stocks);
  for(int s = 0; s< n_stocks; s++) for(int y = 0; y < n_years_model; y++) {
    if(XSPR_R_opt == 1) R_XSPR(y,s) = NAA(s,spawn_regions(s)-1,y,0);
    if(XSPR_R_opt == 3) R_XSPR(y,s) = pred_NAA(s,spawn_regions(s)-1, y,0);
  }
  if((XSPR_R_opt == 2) | (XSPR_R_opt == 4)){
    vector<Type> avg_R(n_stocks);
    for(int s = 0; s< n_stocks; s++) {
      for(int y = 0; y < XSPR_R_avg_yrs.size(); y++) {
        if(XSPR_R_opt == 2) avg_R(s) += NAA(s,spawn_regions(s)-1,XSPR_R_avg_yrs(y),0);
        if(XSPR_R_opt == 4) avg_R(s) += pred_NAA(s,spawn_regions(s)-1,XSPR_R_avg_yrs(y),0);
      }
      avg_R(s) /= Type(XSPR_R_avg_yrs.size());
      for(int y = 0; y < n_years_model + n_years_proj; y++) R_XSPR(y,s) = avg_R(s);
    }
  }
  return R_XSPR;
}


/* calculate single SSB/R at F for spatial model across stocks and regions */
template<class Type>
struct spr_F_spatial {
  // Data and parameter objects for calculation: 
  vector<int> spawn_seasons;
  vector<int> spawn_regions;
  vector<int> fleet_regions; 
  matrix<int> fleet_seasons; 
  array<int> can_move; 
  vector<int> mig_type;
  vector<Type> fracyr_SSB;
  matrix<Type> selectivity;
  array<Type> log_M;
  array<Type> mu; 
  vector<Type> L;
  matrix<Type> mature;
  matrix<Type> waa_ssb; 
  vector<Type> fracyr_seasons;
  vector<Type> SPR_weights; //how to weight stock-specific SSB/R for aggregate SSB/R.
  int age_specific; 
  int small_dim;
  int trace=0;

  // Constructor 
  spr_F_spatial(
  vector<int> spawn_seasons_,
  vector<int> spawn_regions_,
  vector<int> fleet_regions_, 
  matrix<int> fleet_seasons_, 
  array<int> can_move_, 
  vector<int> mig_type_,
  vector<Type> fracyr_SSB_,
  matrix<Type> selectivity_,
  array<Type> log_M_,
  array<Type> mu_, 
  vector<Type> L_,
  matrix<Type> mature_,
  matrix<Type> waa_ssb_, 
  vector<Type> fracyr_seasons_,
  vector<Type> SPR_weights_,
  int age_specific_, 
  int small_dim_,
  int trace_) :
    spawn_seasons(spawn_seasons_), 
    spawn_regions(spawn_regions_), 
    fleet_regions(fleet_regions_),
    fleet_seasons(fleet_seasons_),
    can_move(can_move_),
    mig_type(mig_type_),
    fracyr_SSB(fracyr_SSB_),
    selectivity(selectivity_),
    log_M(log_M_),
    mu(mu_),
    L(L_),
    mature(mature_),
    waa_ssb(waa_ssb_),
    fracyr_seasons(fracyr_seasons_),
    SPR_weights(SPR_weights_),
    age_specific(age_specific_),
    small_dim(small_dim_), 
    trace(trace_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_stocks = log_M.dim(0);
    int n_regions = log_M.dim(1);
    int n_ages = log_M.dim(2);
    int n_fleets = selectivity.rows();
    if(trace) see("in spr_F_spatial");
    //int n_seasons = can_move.dim(2);
    matrix<T> FAA(n_fleets,n_ages);
    FAA.setZero();
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) {
      FAA(f,a) = T(selectivity(f,a)) * exp(log_F(0));
    }
    if(trace) see(FAA);
    vector<T> fracyrssbT = fracyr_SSB.template cast<T>();
    array<T> logMbaseT(log_M.dim(0),log_M.dim(1),log_M.dim(2));
    for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) for(int a = 0; a < n_ages; a++){
      logMbaseT(s,r,a) = T(log_M(s,r,a));
    }
    if(trace) see(logMbaseT);
    if(trace) see(mu.dim);
    array<T> muT(mu.dim(0),mu.dim(1),mu.dim(2),mu.dim(3),mu.dim(4));
    if(trace) see(muT.dim);
    if(trace) see(n_stocks);
    if(trace) see(n_ages);
    if(trace) see(n_regions);
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int t = 0; t < mu.dim(2); t++) {
      for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) {
        muT(s,a,t,r,rr) = T(mu(s,a,t,r,rr));
      }
    }
    if(trace) see(muT);
    vector<T> LT = L.template cast<T>();
    matrix<T> matT = mature.template cast<T>();
    matrix<T> waassbT = waa_ssb.template cast<T>();
    vector<T> fracyrseasonT = fracyr_seasons.template cast<T>();
    //stock-specific SSB/R
    array<T> SPR_sr = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, 
      fracyrssbT, 
      FAA, 
      logMbaseT, 
      muT, 
      LT, 
      matT, 
      waassbT, 
      fracyrseasonT, 
      0, small_dim, trace);
    if(trace) see("end get_SPR spr_F_spatial");

    //weighted-average of SSB/R returned
    T SPR = 0;
    for(int s = 0; s < n_stocks; s++) SPR += T(SPR_weights(s)) * SPR_sr(s,spawn_regions(s)-1,spawn_regions(s)-1); 
    if(trace) see(SPR);
    return SPR;
  }
};

//takes a single year of values for inputs (reduce dimensions appropriately)
//returns just the "solved" log_FXSPR value
template <class Type>
Type get_FXSPR(vector<int> spawn_seasons, vector<int> spawn_regions, vector<int> fleet_regions, vector<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, vector<Type> ssbfrac, matrix<Type> sel, array<Type> log_M, array<Type> mu, 
  vector<Type> L, matrix<Type> mat,  matrix<Type> waassb, vector<Type> fracyr_seasons, vector<Type> R_XSPR, 
  Type percentSPR, vector<Type> SPR_weights, int SPR_weight_type, int small_dim, Type F_init, int n_iter, int trace) {
  int n = n_iter;
  int n_stocks = spawn_seasons.size();
  int n_fleets = fleet_regions.size();
  int n_ages = mat.cols();
  //define SPR weighting 
  if(SPR_weight_type == 0){ //use average Recruitment
    SPR_weights = R_XSPR/R_XSPR.sum();
  } else { //use user-specified weights as provided. do nothing
  }
  if(trace) see(SPR_weights);
  matrix<Type> FAA0(n_fleets,n_ages);
  FAA0.setZero();
  array<Type> SPR0_all = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, FAA0, log_M, mu, L, 
    mat, waassb, fracyr_seasons, 0, small_dim);
  if(trace) see(SPR0_all);
  Type SPR0 = 0;
  for(int s = 0; s < n_stocks; s++) SPR0 += SPR_weights(s) * SPR0_all(s,spawn_regions(s)-1,spawn_regions(s)-1); 
  if(trace) see(SPR0);

  vector<Type> log_FXSPR_i(1), log_FXSPR_iter(n);
  log_FXSPR_iter(0) = log(F_init);

  if(trace) see(log_FXSPR_iter(0));
  spr_F_spatial<Type> sprF(spawn_seasons, spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, sel, log_M,
    mu, L, mat, waassb, fracyr_seasons, SPR_weights, 0, small_dim, trace);
  if(trace) see("after spr_F_spatial sprF defined");
  for(int i=0; i<n-1; i++) {
    if(trace) see(i);
    log_FXSPR_i(0) = log_FXSPR_iter(i);
    vector<Type> grad_spr_F = autodiff::gradient(sprF,log_FXSPR_i);
    log_FXSPR_iter(i+1) = log_FXSPR_iter(i) - (sprF(log_FXSPR_i) - 0.01*percentSPR * SPR0)/grad_spr_F(0);
  }
  Type FXSPR = exp(log_FXSPR_iter(n-1));
  return FXSPR;
}

//takes a single year of values for inputs including log_SPR0 (reduce dimensions appropriately)
//returns just the "solved" log_FXSPR value
template <class Type>
Type get_FXSPR(vector<int> spawn_seasons, vector<int> spawn_regions, vector<int> fleet_regions, vector<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, vector<Type> ssbfrac, matrix<Type> sel, array<Type> log_M, array<Type> mu, 
  vector<Type> L, matrix<Type> mat,  matrix<Type> waassb, vector<Type> fracyr_seasons, vector<Type> R_XSPR, vector<Type> log_SPR0,
  Type percentSPR, vector<Type> SPR_weights, int SPR_weight_type, int small_dim, Type F_init, int n_iter, int trace) {
  int n = n_iter;
  int n_stocks = spawn_seasons.size();
  int n_fleets = fleet_regions.size();
  int n_ages = mat.cols();
  //define SPR weighting 
  if(SPR_weight_type == 0){ //use average Recruitment
    SPR_weights = R_XSPR/R_XSPR.sum();
  } else { //use user-specified weights as provided. do nothing
  }
  if(trace) see(SPR_weights);
  matrix<Type> FAA0(n_fleets,n_ages);
  FAA0.setZero();
  Type SPR0 = 0;
  for(int s = 0; s < n_stocks; s++) SPR0 += SPR_weights(s) * exp(log_SPR0(s)); 
  if(trace) see(SPR0);

  vector<Type> log_FXSPR_i(1), log_FXSPR_iter(n);
  log_FXSPR_iter(0) = log(F_init);

  if(trace) see(log_FXSPR_iter(0));
  spr_F_spatial<Type> sprF(spawn_seasons, spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, sel, log_M,
    mu, L, mat, waassb, fracyr_seasons, SPR_weights, 0, small_dim, trace);
  if(trace) see("after spr_F_spatial sprF defined");
  for(int i=0; i<n-1; i++) {
    if(trace) see(i);
    log_FXSPR_i(0) = log_FXSPR_iter(i);
    vector<Type> grad_spr_F = autodiff::gradient(sprF,log_FXSPR_i);
    log_FXSPR_iter(i+1) = log_FXSPR_iter(i) - (sprF(log_FXSPR_i) - 0.01*percentSPR * SPR0)/grad_spr_F(0);
  }
  Type FXSPR = exp(log_FXSPR_iter(n-1));
  return FXSPR;
}

//returns annual values of 
template <class Type>
vector<Type> get_log_FXSPR(Type percentSPR, array<Type> FAA, vector<int> fleet_regions, vector<int> fleet_seasons, 
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
    matrix<Type> waa_ssb_y = get_matrix_y(waa_ssb, y);
    //matrix<Type> waa_catch_y = get_matrix_y(waa_catch, y);
    matrix<Type> mature_y = get_matrix_y(mature, y);
    matrix<Type> sel_y = get_avg_fleet_sel(FAA, yvec, which_F_age(y));
    if(trace) see(sel_y);
    vector<Type> L_y = L.row(y);
    if(trace) see(L_y);
    array<Type> log_M_y = get_log_M_y(y,log_M);
    if(trace) see(log_M_y);
    array<Type> mu_y = get_mu_y(y,mu);
    if(trace) see(mu_y);
    if(trace) see(mu_y.dim);

    Type FXSPR = get_FXSPR(spawn_seasons, spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, vector<Type> (fracyr_SSB.row(y)),
      sel_y, log_M_y, mu_y, L_y, mature_y,  waa_ssb_y, fracyr_seasons, vector<Type> (R_XSPR.row(y)), vector<Type> (log_SPR0.row(y)),
      percentSPR, SPR_weights, SPR_weight_type, small_dim, FXSPR_init(y), 10, trace);
    log_FXSPR(y) = log(FXSPR);
    if(trace) see(log_FXSPR(y));
  }

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
vector< vector <Type> > get_SPR_res(vector<Type> SPR_weights, array<Type> log_M, array<Type> FAA, vector<int> spawn_seasons,  
  vector<int> spawn_regions,
  vector<int> fleet_regions, 
  matrix<int> fleet_seasons,
  vector<Type> fracyr_seasons,
  array<int> can_move,
  array<int> must_move,
  vector<int> mig_type,
  array<Type> trans_mu_base, 
  matrix<Type> L,
  int which_F_age, array<Type> waa, vector<int> waa_pointer_ssb, 
  vector<int> waa_pointer_fleets,
  array<Type> mature, Type percentSPR, array<Type> NAA, matrix<Type> fracyr_SSB, Type F_init, 
  vector<int> years_M, vector<int> years_mu, vector<int> years_L, vector<int> years_mat, vector<int> years_sel, 
  vector<int> years_waa_ssb, vector<int> years_waa_catch, vector<Type> R_XSPR,
  int small_dim, int SPR_weight_type, int trace = 0) {
  if(trace) see("inside get_SPR_res");
  //see(SPR_weights);
  //see(log_M);
  int n = 10;
  int n_stocks = log_M.dim(0);
  int n_fleets = FAA.dim(0);
  int n_regions = NAA.dim(1);
  int n_seasons = fleet_seasons.cols();
  int n_ages = log_M.dim(3);
  //see(n_ages);
  //see(fleet_seasons(10,10));
  matrix<Type> waa_ssb(n_stocks,n_ages);
  matrix<Type> waa_catch(n_fleets,n_ages);
  matrix<Type> sel(n_fleets, n_ages);
  array<Type> M(n_stocks,n_regions,n_ages);
  array<Type> log_M_avg(n_stocks,n_regions,n_ages);
  array<Type> mu_avg(n_stocks,n_seasons,n_ages, n_regions, n_regions);
  vector<Type> L_avg(n_regions);
  matrix<Type> mat(n_stocks, n_ages);
  vector<Type> ssbfrac(n_stocks); //R(n_stocks); R.setZero(); 

  waa_ssb.setZero(); waa_catch.setZero(); sel.setZero(); M.setZero(); mu_avg.setZero(); L_avg.setZero(); 
  mat.setZero(); ssbfrac.setZero(); log_M_avg.setZero();
  //get average inputs over specified years

  ssbfrac = get_avg_ssbfrac(fracyr_SSB,years_waa_ssb);
  if(trace) see(ssbfrac);

  waa_ssb = get_avg_waa(waa, years_waa_ssb, waa_pointer_ssb);
  if(trace) see(waa_ssb);
  waa_catch = get_avg_waa(waa, years_waa_catch, waa_pointer_fleets);
  if(trace) see(waa_catch);
  L_avg = get_avg_L(L, years_L, 0);
  if(trace) see(L_avg);
  mat = get_avg_mat(mature,years_mat);
  if(trace) see(mat);
  log_M_avg = get_avg_M(log_M, years_M, 1);
  if(trace) see(log_M_avg);
  if(n_regions>1) mu_avg = get_avg_mu(trans_mu_base,years_mu,mig_type, can_move, must_move);
  if(trace) see(mu_avg);

  matrix<Type> FAA_avg = get_avg_FAA(FAA,years_sel,0);
  vector<Type> FAA_avg_tot = FAA_avg.colwise().sum();
  if(trace) see(FAA_avg_tot);

  //which_F_age needs to be set appropriately by user.
  //if an arbitrary index is used, then full F would need to be defined on R side. Could pick out from ADREPORTED FAA
  if(trace) see(which_F_age);
  sel = FAA_avg/FAA_avg_tot(which_F_age-1);
  if(trace) see(sel);

  //define SPR weighting 
  if(SPR_weight_type == 0){ //use average Recruitment
    SPR_weights = R_XSPR/R_XSPR.sum();
  } else { //use user-specified weights as provided. do nothing
  }
  if(trace) see(SPR_weights);
  matrix<Type> FAA0(n_fleets,n_ages);
  FAA0.setZero();
  array<Type> SPR0_all = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, FAA0, log_M_avg, mu_avg, L_avg, 
    mat, waa_ssb, fracyr_seasons, 0, small_dim);
  if(trace) see(SPR0_all);
  Type SPR0 = 0;
  for(int s = 0; s < n_stocks; s++) SPR0 += SPR_weights(s) * SPR0_all(s,spawn_regions(s)-1,spawn_regions(s)-1); 
  if(trace) see(SPR0);

  vector<Type> log_FXSPR_i(1), log_FXSPR_iter(n);
  log_FXSPR_iter(0) = log(F_init);

  if(trace) see(log_FXSPR_iter(0));
  spr_F_spatial<Type> sprF(spawn_seasons, spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, sel, log_M_avg,
    mu_avg, L_avg, mat, waa_ssb, fracyr_seasons, SPR_weights, 0, small_dim, 0);
  if(trace) see("after spr_F_spatial sprF defined");
  for(int i=0; i<n-1; i++) {
    if(trace) see(i);
    log_FXSPR_i(0) = log_FXSPR_iter(i);
    vector<Type> grad_spr_F = autodiff::gradient(sprF,log_FXSPR_i);
    log_FXSPR_iter(i+1) = log_FXSPR_iter(i) - (sprF(log_FXSPR_i) - 0.01*percentSPR * SPR0)/grad_spr_F(0);
  }
  matrix<Type> FAA_XSPR = exp(log_FXSPR_iter(n-1)) * sel;
  if(trace) see(FAA_XSPR);
  vector<Type> log_FAA_XSPR((n_fleets+1)*n_ages);
  log_FAA_XSPR.setZero();
  int k = 0;
  for(int f = 0; f < n_fleets; f++) {
    for(int a = 0; a < n_ages; a++){
      log_FAA_XSPR(k) = log(FAA_XSPR(f,a));
      log_FAA_XSPR(n_fleets*n_ages + a) += FAA_XSPR(f,a); //summing, not log yet
      k++;
    }
  }
  for(int a = 0; a < n_ages; a++) log_FAA_XSPR(n_fleets*n_ages + a) = log(log_FAA_XSPR(n_fleets*n_ages + a)); //log it
  if(trace) see(log_FAA_XSPR);
  array<Type> SPR_all = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, FAA_XSPR, log_M_avg, mu_avg, L_avg, 
    mat, waa_ssb, fracyr_seasons, 0, small_dim);
  if(trace) see(SPR_all);
  array<Type> YPR_all = get_YPR(fleet_regions, fleet_seasons, can_move, mig_type, FAA_XSPR, log_M_avg, mu_avg, L_avg, waa_catch, 
    fracyr_seasons, 0, small_dim); //n_stocks x n_regions x n_fleets (should be 0 for regions not fleet_regions(f)-1)
  //for each stock/fleet and also the weighted average/total
  //see(YPR_all);
  vector<Type> log_SPR(n_stocks + 1), log_SPR0(n_stocks+1), log_Y_XSPR(n_fleets+1), log_SSB_XSPR(n_stocks+1); 
  log_SPR.setZero(); log_SPR0.setZero(); log_Y_XSPR.setZero(); log_SSB_XSPR.setZero();
  
  for(int s = 0; s < n_stocks; s++) {
    log_SPR(s) = log(SPR_all(s,spawn_regions(s)-1,spawn_regions(s)-1));
    log_SPR(n_stocks) += SPR_weights(s) * SPR_all(s,spawn_regions(s)-1,spawn_regions(s)-1);
    log_SPR0(s) = log(SPR0_all(s,spawn_regions(s)-1,spawn_regions(s)-1));
    log_SPR0(n_stocks) += SPR_weights(s) * SPR0_all(s,spawn_regions(s)-1,spawn_regions(s)-1);
    log_SSB_XSPR(s) = log(R_XSPR(s)) + log_SPR(s);
    log_SSB_XSPR(n_stocks) += SPR_weights(s) * R_XSPR(s) * SPR_all(s,spawn_regions(s)-1,spawn_regions(s)-1);
    //log_SSB_XSPR(n_stocks) += R_XSPR(s) * SPR_all(s,spawn_regions(s)-1,spawn_regions(s)-1);
    for(int f = 0; f < n_fleets; f++) {
      //NOTE: If stock s cannot ever be in fleet_region(f)-1 then this will give log(0);
      //need to use can_move flag?
      //for(int r = 0; r < n_regions; r++) 
      log_Y_XSPR(f) += R_XSPR(s) * YPR_all(s,fleet_regions(f)-1,f); //not logged yet
      log_Y_XSPR(n_fleets) += log_Y_XSPR(f); //not logged yet
    }
  }
  log_Y_XSPR = log(log_Y_XSPR);
  log_SPR(n_stocks) = log(log_SPR(n_stocks));
  log_SPR0(n_stocks) = log(log_SPR0(n_stocks));
  //see(log_SPR);
  //see(log_SPR0);
  //log_SSB_XSPR(n_stocks) = log(exp(log_SSB_XSPR.head(n_stocks)).sum());
  //log_Y_XSPR(n_fleets) = log(exp(log_Y_XSPR.head(n_fleets)).sum());
  vector< vector<Type> > res(6); 
  res(0) = log_FAA_XSPR; // log_FAA at FXSPR by fleet and across fleets
  res(1) = log_SSB_XSPR; //log_SSB_FXSPR
  res(2) = log_Y_XSPR; //log_Y_FXSPR
  res(3) = log_SPR; //stock specific log SPRs at FXSPR, (only the weighted sum will be X*SPR0/100)
  res(4) = log_SPR0;
  res(5) = log_FXSPR_iter; //last value is max F at X%SPR across ages and fleets
  if(trace) see("end get_SPR_res")
  return res;
}
