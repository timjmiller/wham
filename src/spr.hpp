
template <class T>
array<T> get_SPR(vector<int> spawn_seasons, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, vector<T> fracyr_SSB, array<T> FAA, array<T> log_M, array<T> mu, vector<T> L, 
  array<T> mature, array<T> waa_ssb, vector<T> fracyr_seasons, int age_specific, int small_dim, int trace = 0, int numbers = 0){
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
        if(n_regions>1) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) {
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
      matrix<T> SPR_ya = S_ya; //eq abundance/recruit (Jan 1)
      if(numbers==0) SPR_ya = SPR_ya * get_S(P_spawn, n_regions) * mature(s,a) * waa_ssb(s,a); //SSB/recruit
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

template <class T>
array<T> get_SPR(vector<int> spawn_seasons, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, vector<T> fracyr_SSB, array<T> FAA, array<T> log_M, array<T> mu, vector<T> L, 
  array<T> mature, array<T> waa_ssb, vector<T> fracyr_seasons, int age_specific, int bias_correct, 
  array<T> marg_NAA_sigma, int small_dim, int trace = 0, int numbers = 0){
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
    matrix<T> cum_S_ya = get_S(I, n_regions); //cumulative survival up to age a
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
        if(n_regions>1) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) {
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
      // SSB/R at year and age = prob alive to up to age a x prob spawn at age a x waa x mature
      //should be n_regions x n_regions
      matrix<T> NPR_ya = cum_S_ya; //eq abundance at age a/recruit (at beginning of year y: Jan 1). Identity matrix for age 1 (0 in c++)
      if(trace) see(NPR_ya);
      
      matrix<T> S_ya = get_S(P_ya, n_regions); //survival over year y for age a up to end of year y (or to age a+1 at beginning of year y+1)
      if(bias_correct) for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
        //bias-correct for next age going to region j?
        int abc = a;
        if(a < n_ages-1) abc += 1;
        S_ya(i,j) = S_ya(i,j) * exp(-0.5*pow(marg_NAA_sigma(s,j,abc),2)); 
      }
      if(a < n_ages-1) cum_S_ya = cum_S_ya * S_ya; //accumulate for next age
      if(trace) see(cum_S_ya);
      if(a == n_ages-1){ //plus group
        matrix<T> fundm = get_S(I, n_regions) - S_ya; //S_ya already has any bias correction by column. 
        if(small_dim) fundm = fundm.inverse(); else fundm = atomic::matinv(fundm); //fundm = (I - S_y,+)^-1
        //for plus group cum_S_ya = S_y,a-1 x (I - S_y,+)^-1
        cum_S_ya = cum_S_ya * fundm;
        NPR_ya = cum_S_ya; //change NPR_ya for plus group
        if(trace) see(cum_S_ya);
      }
      //if numbers == 1 return numbers/recruit
      matrix<T> SPR_ya = NPR_ya; 
      //eq. SSB at age a (at time of spawning)/recruit
      if(numbers==0) SPR_ya = SPR_ya * get_S(P_spawn, n_regions) * mature(s,a) * waa_ssb(s,a); 
      if(trace) see(SPR_ya);
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
        //SSB per Recruit in each region (should only have positive values in spawn_regions(s)-1?)
        SPRAA(s,a,i,j) += SPR_ya(i,j); 
      }
    }
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

template <class T>
array<T> get_YPR_srf(vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, array<T> FAA, array<T> log_M, array<T> mu, vector<T> L, array<T> waacatch, 
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
        P_ya = P_ya * get_P_t(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA.matrix(), log_M, mu, L);
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
      P_ya = P_ya * get_P_t(n_ages-1, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA.matrix(), log_M, mu, L);
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

template <class T>
array<T> get_YPR_srf(vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, array<T> FAA, array<T> log_M, array<T> mu, vector<T> L, array<T> waacatch, 
  vector<T> fracyr_seasons, int age_specific, int bias_correct, 
  array<T> marg_NAA_sigma, int small_dim){
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
    matrix<T> cum_S_ya = get_S(I, n_regions); //cumulative survival up to age a
    for(int a = 0; a < n_ages; a++) {
      for(int i = 0; i < n_fleets; i++) W(i,i) = waacatch(i,a);
      matrix<T> P_ya = I; //PTM for year and age
      for(int t = 0; t < n_seasons; t++) {
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P_t(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA.matrix(), log_M, mu, L);
      }
      // N/recruit at beginning of year (Jan 1) y at age a. Identity matrix for age 1 (0 in c++)
      matrix<T> NPR_ya = cum_S_ya; 
      matrix<T> S_ya = get_S(P_ya, n_regions); //survival over year y for age a up to end of year y (or to age a+1 at beginning of year y+1)
      if(bias_correct) for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
        //bias-correct for next age going to region j?
        int abc = a;
        if(a < n_ages-1) abc += 1;
        S_ya(i,j) = S_ya(i,j) * exp(-0.5*pow(marg_NAA_sigma(s,j,abc),2)); 
      }
      if(a < n_ages-1) cum_S_ya = cum_S_ya * S_ya; //accumulate for next age
      if(a == n_ages-1){ //plus group
        matrix<T> fundm = get_S(I, n_regions) - S_ya; //already has any bias correction by column.
        if(small_dim) fundm = fundm.inverse(); else fundm = atomic::matinv(fundm); //fundm = (I - S_y,+)^-1
        //for plus group cum_S_ya = S_y,a-1 x (I - S_y,+)^-1
        cum_S_ya = cum_S_ya * fundm;
        NPR_ya = cum_S_ya; //revise N/recruit on jan 1 for plus group
      }
      matrix<T> YPR_ya = NPR_ya * get_D(P_ya, n_regions,n_fleets) * W; //(n_regions x n_fleets) yield/recruit starting year in region r and being caught in fleet f.
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) YPRAA(s,a,i,j) = YPR_ya(i,j);
    }
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
  int XSPR_R_opt, vector<int> XSPR_R_avg_yrs, 
  array<Type> marg_NAA_sigma){
  // array<Type> log_NAA_sigma){
  //R_XSPR is needed for projections and reference points
  array<Type> NAA = extract_NAA(all_NAA);
  array<Type> pred_NAA = extract_pred_NAA(all_NAA);
  int n_stocks = NAA.dim(0);
  matrix<Type> R_XSPR(n_years_model+n_years_proj, n_stocks);
  R_XSPR.setZero();
  for(int s = 0; s< n_stocks; s++) for(int y = 0; y < n_years_model; y++) {
    if(XSPR_R_opt == 1) R_XSPR(y,s) = NAA(s,spawn_regions(s)-1,y,0);
    if(XSPR_R_opt == 3) R_XSPR(y,s) = pred_NAA(s,spawn_regions(s)-1, y,0);
  }
  if((XSPR_R_opt == 2) | (XSPR_R_opt == 4) | (XSPR_R_opt == 5)){
    vector<Type> avg_R(n_stocks);
    avg_R.setZero();
    for(int s = 0; s< n_stocks; s++) {
      for(int y = 0; y < XSPR_R_avg_yrs.size(); y++) {
        if(XSPR_R_opt == 2) avg_R(s) += NAA(s,spawn_regions(s)-1,XSPR_R_avg_yrs(y),0);
        if(XSPR_R_opt == 4) avg_R(s) += pred_NAA(s,spawn_regions(s)-1,XSPR_R_avg_yrs(y),0);
        if(XSPR_R_opt == 5) avg_R(s) += pred_NAA(s,spawn_regions(s)-1,XSPR_R_avg_yrs(y),0) * exp(-0.5*pow(marg_NAA_sigma(s,spawn_regions(s)-1,0),2));
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
  array<Type> selectivity;
  array<Type> log_M;
  array<Type> mu; 
  vector<Type> L;
  array<Type> mature;
  array<Type> waa_ssb; 
  vector<Type> fracyr_seasons;
  vector<Type> SPR_weights; //how to weight stock-specific SSB/R for aggregate SSB/R.
  int age_specific; 
  int bias_correct;
  array<Type> marg_NAA_sigma;
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
  array<Type> selectivity_,
  array<Type> log_M_,
  array<Type> mu_, 
  vector<Type> L_,
  array<Type> mature_,
  array<Type> waa_ssb_, 
  vector<Type> fracyr_seasons_,
  vector<Type> SPR_weights_,
  int age_specific_, 
  int bias_correct_,
  array<Type> marg_NAA_sigma_,
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
    bias_correct(bias_correct_),
    marg_NAA_sigma(marg_NAA_sigma_),
    small_dim(small_dim_), 
    trace(trace_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  vector<T> operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_stocks = log_M.dim(0);
    int n_regions = log_M.dim(1);
    int n_ages = log_M.dim(2);
    int n_fleets = selectivity.rows();
    if(trace) see("in spr_F_spatial");
    //int n_seasons = can_move.dim(2);
    array<T> FAA(n_fleets,n_ages);
    FAA.setZero();
    if(trace) see(selectivity);
    int n_F = log_F.size(); //MUST be 1 or n_regions!!!!
    if(trace) see(n_F);
    if(n_F == 1){ //total F across regions or there is just one region
      for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) {
        FAA(f,a) = T(selectivity(f,a)) * exp(log_F(0));
      }
    } else{ //total F by region
      for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) {
        FAA(f,a) = T(selectivity(f,a)) * exp(log_F(fleet_regions(f)-1));
      }
    }
    if(trace) see(FAA);
    vector<T> fracyrssbT = fracyr_SSB.template cast<T>();
    array<T> logMbaseT(n_stocks,n_regions,n_ages), marg_NAA_sigmaT(n_stocks, n_regions, n_ages);
    // array<T> logMbaseT(n_stocks,n_regions,n_ages), log_NAA_sigmaT(n_stocks, n_regions, n_ages);
    for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) for(int a = 0; a < n_ages; a++){
      logMbaseT(s,r,a) = T(log_M(s,r,a));
      marg_NAA_sigmaT(s,r,a) = T(marg_NAA_sigma(s,r,a));
    }
    if(trace) see(logMbaseT.dim);
    if(trace) see(mu.dim);
    array<T> muT(mu.dim(0),mu.dim(1),mu.dim(2),mu.dim(3),mu.dim(4));
    if(trace) see(n_stocks);
    if(trace) see(n_ages);
    if(trace) see(n_regions);
    if(n_regions>1) for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int t = 0; t < mu.dim(2); t++) {
      for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) {
        muT(s,a,t,r,rr) = T(mu(s,a,t,r,rr));
      }
    }
    if(trace) see(muT.dim);
    vector<T> LT = L.template cast<T>();
    if(trace) see(LT);
    array<T> matT(mature.dim(0),mature.dim(1));
    for(int i = 0; i < mature.dim(0); i++) for(int j = 0; j < mature.dim(1); j++) matT(i,j) = T(mature(i,j));
    if(trace) see(matT);
    array<T> waassbT(waa_ssb.dim(0), waa_ssb.dim(1));
    for(int i = 0; i < waa_ssb.dim(0); i++) for(int j = 0; j < waa_ssb.dim(1); j++) waassbT(i,j) = T(waa_ssb(i,j));
    if(trace) see(waassbT);
    vector<T> fracyrseasonT = fracyr_seasons.template cast<T>();
    if(trace) see(fracyrseasonT);
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
      0, 
      bias_correct,
      marg_NAA_sigmaT,
      // log_NAA_sigmaT,
      small_dim, 0, 0);
    if(trace) see("end get_SPR spr_F_spatial");

    //weighted-average of SSB/R returned
    vector<T> SPR(n_F);
    SPR.setZero();
    for(int s = 0; s < n_stocks; s++) {
      if(n_F==1) {
        SPR(0) += T(SPR_weights(s)) * SPR_sr(s,spawn_regions(s)-1,spawn_regions(s)-1); 
      } else {
        SPR(spawn_regions(s)-1) += SPR_sr(s,spawn_regions(s)-1,spawn_regions(s)-1);
      }
    }
    if(trace) see(SPR_weights);
    if(trace) see(SPR);
    
    if(trace) {
      vector<T> SPR_s(n_stocks);      
      for(int s = 0; s < n_stocks; s++) SPR_s(s) = SPR_sr(s,spawn_regions(s)-1,spawn_regions(s)-1);
      see(SPR_s);

      array<T> FAA0(n_fleets,n_ages);
      FAA0.setZero();
      array<T> SPR0_sr = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, 
        fracyrssbT, 
        FAA0, 
        logMbaseT, 
        muT, 
        LT, 
        matT, 
        waassbT, 
        fracyrseasonT, 
        0, 
        bias_correct,
        marg_NAA_sigmaT,
        small_dim, 0, 0);
      T SPR0 = 0;
      for(int s = 0; s < n_stocks; s++) SPR0 += T(SPR_weights(s)) * SPR0_sr(s,spawn_regions(s)-1,spawn_regions(s)-1); 
      if(trace) see(SPR0);
    }

    return SPR;
  }
};

//takes a single year of values for inputs (reduce dimensions appropriately)
//returns just the "solved" log_FXSPR value
template <class Type>
vector<Type> get_FXSPR(vector<int> spawn_seasons, vector<int> spawn_regions, vector<int> fleet_regions, matrix<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, vector<Type> ssbfrac, array<Type> sel, array<Type> log_M, array<Type> mu, 
  vector<Type> L, array<Type> mat,  array<Type> waassb, vector<Type> fracyr_seasons, vector<Type> R_XSPR, 
  Type percentSPR, vector<Type> SPR_weights, int SPR_weight_type, int bias_correct, 
  array<Type> marg_NAA_sigma, 
  // array<Type> log_NAA_sigma, 
  int small_dim, Type F_init, int n_iter, int trace, int by_region = 0) {
  int n_stocks = spawn_seasons.size();
  int n_fleets = fleet_regions.size();
  int n_ages = mat.cols();
  int n_regions = log_M.dim(1);
  //define SPR weighting 
  if(SPR_weight_type == 0){ //use average Recruitment
    SPR_weights = R_XSPR/R_XSPR.sum();
  } else { //use user-specified weights as provided. do nothing
  }
  if(trace) see(SPR_weights);
  array<Type> FAA0(n_fleets,n_ages);
  FAA0.setZero();
  array<Type> SPR0_all = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, FAA0, log_M, mu, L, 
    mat, waassb, fracyr_seasons, 0, bias_correct, 
    marg_NAA_sigma, 
    small_dim, 0, 0);
  if(trace) see(SPR0_all);
  int n_F = 1;
  if(by_region) n_F = n_regions;
  vector<Type> SPR0(n_F);
  SPR0.setZero();
  for(int s = 0; s < n_stocks; s++) {
    if(n_F==1) {
      SPR0(0) += SPR_weights(s) * SPR0_all(s,spawn_regions(s)-1,spawn_regions(s)-1);
    } else {
      SPR0(spawn_regions(s)-1) += SPR_weights(s) * SPR0_all(s,spawn_regions(s)-1,spawn_regions(s)-1);
    }
  }

  if(trace) see(SPR0);
  vector<Type> log_FXSPR_i(n_F);
  matrix<Type> log_FXSPR_iter(n_iter,n_F);
  for(int r = 0; r < n_F; r++) log_FXSPR_iter(0,r) = log(F_init);
  if(trace) see(log_FXSPR_iter.row(0));
  spr_F_spatial<Type> sprF(spawn_seasons, spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, sel, log_M,
    mu, L, mat, waassb, fracyr_seasons, SPR_weights, 0, bias_correct, 
    marg_NAA_sigma, 
    // log_NAA_sigma, 
    small_dim, trace);    
  if(trace) see("after spr_F_spatial sprF defined");
  //trace = 0;
  for(int i=0; i<n_iter-1; i++) {
    // if(trace) see(i);
    log_FXSPR_i = log_FXSPR_iter.row(i);
    matrix<Type> grad_spr_F = autodiff::jacobian(sprF,log_FXSPR_i);
    matrix<Type> inv_grad = grad_spr_F.inverse(); 
    vector<Type> SPR_i = sprF(log_FXSPR_i);
    vector<Type> diff = SPR_i - 0.01*percentSPR * SPR0;
    vector<Type> change = inv_grad * diff; 
    log_FXSPR_iter.row(i+1) = vector<Type> (log_FXSPR_iter.row(i)) - change; // vector<Type>(grad_spr_F.inverse() * (SPR_i - 0.01*percentSPR * SPR0));
  }
  if(trace) see(log_FXSPR_iter);

  if(trace) {
    see(sel);
    see(spawn_seasons);
    see(fleet_regions);
    see( fleet_seasons);
    see( can_move);
    see( mig_type); 
    see(ssbfrac); 
    see(log_M)
    see(L);
    see(mat);
    see(waassb);
    see( fracyr_seasons);
    see(small_dim);
    
    // would need to define this before the code can be uncommented.
    // see(log_FXSPR_static(0));
    // for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++){
    //   FAA0(f,a) = exp(log_FXSPR_static(0)) * sel(f,a);
    // }
    // array<Type> SPR_FXSPR_all = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, FAA0, log_M, mu, L, 
    //   mat, waassb, fracyr_seasons, 0, small_dim, 0, 0);
    // Type SPR_FXSPR_alt = 0;
    // for(int s = 0; s < n_stocks; s++) SPR_FXSPR_alt += SPR_weights(s) * SPR_FXSPR_all(s,spawn_regions(s)-1,spawn_regions(s)-1); 
    // if(trace) see(SPR_FXSPR_alt);
  }

  vector<Type> FXSPR = exp(vector<Type> (log_FXSPR_iter.row(n_iter-1)));
  return FXSPR;
}

//takes a single year of values for inputs including log_SPR0 (reduce dimensions appropriately)
//returns just the "solved" log_FXSPR value
template <class Type>
vector<Type> get_FXSPR(vector<int> spawn_seasons, vector<int> spawn_regions, vector<int> fleet_regions, matrix<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, vector<Type> ssbfrac, array<Type> sel, array<Type> log_M, array<Type> mu, 
  vector<Type> L, array<Type> mat,  array<Type> waassb, vector<Type> fracyr_seasons, vector<Type> R_XSPR, vector<Type> log_SPR0,
  Type percentSPR, vector<Type> SPR_weights, int SPR_weight_type, int bias_correct, 
  array<Type> marg_NAA_sigma, 
  int small_dim, Type F_init, int n_iter, int trace, int by_region = 0) {
  int n_stocks = spawn_seasons.size();
  int n_fleets = fleet_regions.size();
  int n_ages = mat.cols();
  int n_regions = log_M.dim(1);
  //define SPR weighting 
  if(SPR_weight_type == 0){ //use average Recruitment
    SPR_weights = R_XSPR/R_XSPR.sum();
  } else { //use user-specified weights as provided. do nothing
  }
  if(trace) see(SPR_weights);
  int n_F = 1;
  if(by_region) n_F = n_regions;
  vector<Type> SPR0(n_F);
  SPR0.setZero();
  for(int s = 0; s < n_stocks; s++) {
    if(n_F==1) {
      SPR0(0) += SPR_weights(s) * exp(log_SPR0(s));
    } else {
      SPR0(spawn_regions(s)-1) += SPR_weights(s) * exp(log_SPR0(s));
    }
  }
  if(trace) see(SPR0);
  
  vector<Type> log_FXSPR_i(n_F);
  matrix<Type> log_FXSPR_iter(n_iter,n_F);
  for(int r = 0; r < n_F; r++) log_FXSPR_iter(0,r) = log(F_init);
  if(trace) see(log_FXSPR_iter.row(0));
  spr_F_spatial<Type> sprF(spawn_seasons, spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, sel, log_M,
    mu, L, mat, waassb, fracyr_seasons, SPR_weights, 0, bias_correct, 
    marg_NAA_sigma, 
    small_dim, trace);    
  if(trace) see("after spr_F_spatial sprF defined");
  //trace = 0;
  for(int i=0; i<n_iter-1; i++) {
    // if(trace) see(i);
    log_FXSPR_i = log_FXSPR_iter.row(i);
    matrix<Type> grad_spr_F = autodiff::jacobian(sprF,log_FXSPR_i);
    matrix<Type> inv_grad = grad_spr_F.inverse(); 
    vector<Type> SPR_i = sprF(log_FXSPR_i);
    vector<Type> diff = SPR_i - 0.01*percentSPR * SPR0;
    vector<Type> change = inv_grad * diff; 
    log_FXSPR_iter.row(i+1) = vector<Type> (log_FXSPR_iter.row(i)) - change; // vector<Type>(grad_spr_F.inverse() * (SPR_i - 0.01*percentSPR * SPR0));
  }
  if(trace) see(log_FXSPR_iter);

  vector<Type> FXSPR = exp(vector<Type> (log_FXSPR_iter.row(n_iter-1)));
  return FXSPR;
}


template <class Type>
vector< array <Type> > get_SPR_res(vector<Type> SPR_weights, array<Type> log_M, array<Type> FAA, vector<int> spawn_seasons,  
  vector<int> spawn_regions,
  vector<int> fleet_regions, 
  matrix<int> fleet_seasons,
  vector<Type> fracyr_seasons,
  array<int> can_move,
  array<int> must_move,
  vector<int> mig_type,
  array<Type> trans_mu_base, 
  matrix<Type> L,
  int which_F_age, array<Type> waa_ssb, array<Type> waa_catch, 
  array<Type> mature, Type percentSPR, array<Type> NAA, matrix<Type> fracyr_SSB, Type F_init, 
  vector<int> years_M, vector<int> years_mu, vector<int> years_L, vector<int> years_mat, vector<int> years_sel, 
  vector<int> years_waa_ssb, vector<int> years_waa_catch, vector<Type> R_XSPR,
  int small_dim, int SPR_weight_type, int bias_correct, 
  array<Type> marg_NAA_sigma, 
  int trace = 0, int n_iter = 10) {
  //gets SPR-based BRP information for a year, or inputs may be averaged over specified years.  
  if(trace) see("inside get_SPR_res");
  //see(SPR_weights);
  //see(log_M);
  int n_stocks = log_M.dim(0);
  int n_fleets = FAA.dim(0);
  int n_regions = NAA.dim(1);
  int n_seasons = fleet_seasons.cols();
  int n_ages = log_M.dim(3);
  //see(n_ages);
  //see(fleet_seasons(10,10));

  //get average inputs over specified years
  vector<Type> ssbfrac = get_avg_ssbfrac(fracyr_SSB,years_waa_ssb);
  if(trace) see(ssbfrac);
  array<Type> waa_ssb_avg = get_avg_mat_as_array(waa_ssb, years_waa_ssb); //matrix
  if(trace) see(waa_ssb);
  array<Type> waa_catch_avg = get_avg_mat_as_array(waa_catch, years_waa_catch);//matrix
  if(trace) see(waa_catch);
  vector<Type> L_avg = get_avg_L(L, years_L, 0); //vector
  if(trace) see(L_avg);
  array<Type> mat = get_avg_mat_as_array(mature,years_mat);
  if(trace) see(mat);
  array<Type> log_M_avg = get_avg_M(log_M, years_M, 1);
  if(trace) see(log_M_avg);
  array<Type> mu_avg(n_stocks, n_ages, n_seasons, n_regions, n_regions);
  mu_avg.setZero(); 
  if(n_regions>1) mu_avg = get_avg_mu(trans_mu_base,years_mu,mig_type, can_move, must_move);
  if(trace) see(mu_avg);
  array<Type> FAA_avg = get_avg_FAA_as_array(FAA,years_sel,0);
  if(trace) see(FAA_avg);
  vector<Type> FAA_avg_tot = FAA_avg.matrix().colwise().sum();
  if(trace) see(FAA_avg_tot);

  //which_F_age needs to be set appropriately by user.
  //if an arbitrary index is used, then full F would need to be defined on R side. Could pick out from ADREPORTED FAA
  if(trace) see(which_F_age);
  array<Type> sel(n_fleets,n_ages);
  for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) sel(f,a) = FAA_avg(f,a)/FAA_avg_tot(which_F_age-1);
  if(trace) see(sel);

  //define SPR weighting 
  if(SPR_weight_type == 0){ //use average Recruitment
    SPR_weights = R_XSPR/R_XSPR.sum();
  } else { //use user-specified weights as provided. do nothing
  }
  if(trace) see(SPR_weights);
  array<Type> FAA0(n_fleets,n_ages);
  FAA0.setZero();
  //equil abundance/R (Jan 1)
  array<Type> NAAPR0_all = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, FAA0, log_M_avg, mu_avg, L_avg, 
    mat, waa_ssb_avg, fracyr_seasons, 1, bias_correct, 
    marg_NAA_sigma, 
    small_dim, 0, 1); 
  array<Type> SPRAA0_all = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, FAA0, log_M_avg, mu_avg, L_avg, 
    mat, waa_ssb_avg, fracyr_seasons, 1, bias_correct, 
    marg_NAA_sigma, 
    small_dim, 0, 0);
  if(trace) see(SPRAA0_all.dim);
  array<Type> SPR0_all(n_stocks, n_regions,n_regions);
  SPR0_all.setZero();
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) {
    for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) SPR0_all(s,r,rr) += SPRAA0_all(s,a,r,rr);
  }
  if(trace) see(SPR0_all);
  Type SPR0 = 0;
  for(int s = 0; s < n_stocks; s++) SPR0 += SPR_weights(s) * SPR0_all(s,spawn_regions(s)-1,spawn_regions(s)-1); 
  if(trace) see(SPR0);

  vector<Type> log_FXSPR_i(1), log_FXSPR_iter(n_iter);
  log_FXSPR_iter(0) = log(F_init);

  if(trace) see(log_FXSPR_iter(0));
  // trace = 0;
  spr_F_spatial<Type> sprF(spawn_seasons, spawn_regions, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, sel, log_M_avg,
    mu_avg, L_avg, mat, waa_ssb_avg, fracyr_seasons, SPR_weights, 0, bias_correct, 
    marg_NAA_sigma, 
    // log_NAA_sigma, 
    small_dim, trace);
  if(trace) see("after spr_F_spatial sprF defined");
  // trace = 1;
  for(int i=0; i<n_iter-1; i++) {
    if(trace) see(i);
    log_FXSPR_i(0) = log_FXSPR_iter(i);
    matrix<Type> grad_spr_F = autodiff::jacobian(sprF,log_FXSPR_i);
    //vector<Type> grad_spr_F = autodiff::gradient(sprF,log_FXSPR_i);
    Type SPR_i = sprF(log_FXSPR_i)(0);
    log_FXSPR_iter(i+1) = log_FXSPR_iter(i) - (SPR_i - 0.01*percentSPR * SPR0)/grad_spr_F(0,0);
  }
  array<Type> FAA_XSPR(n_fleets, n_ages);
  array<Type> log_FAA_XSPR(n_fleets+n_regions+1, n_ages);
  log_FAA_XSPR.setZero();
  for(int f = 0; f < n_fleets; f++) {
    for(int a = 0; a < n_ages; a++){
      FAA_XSPR(f,a) = sel(f,a) * exp(log_FXSPR_iter(n_iter-1));
      log_FAA_XSPR(f,a) = log(FAA_XSPR(f,a));
      log_FAA_XSPR(n_fleets+fleet_regions(f)-1, a) += FAA_XSPR(f,a); //summing, not log yet
      log_FAA_XSPR(n_fleets+n_regions, a) += FAA_XSPR(f,a); //summing, not log yet
    }
  }
  if(trace) see(FAA_XSPR);
  for(int a = 0; a < n_ages; a++) for(int r = 0; r <=n_regions; r++) {
    log_FAA_XSPR(n_fleets+r, a) = log(log_FAA_XSPR(n_fleets+r,a)); //log it
  }
  if(trace) see(log_FAA_XSPR);
  //eq NAA/R at FXSPR
  array<Type> NAAPR_FXSPR_all = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, FAA_XSPR, log_M_avg, mu_avg, L_avg, 
    mat, waa_ssb_avg, fracyr_seasons, 1, bias_correct, 
    marg_NAA_sigma, 
    small_dim, 0, 1);
  array<Type> SPRAA_all = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, FAA_XSPR, log_M_avg, mu_avg, L_avg, 
    mat, waa_ssb_avg, fracyr_seasons, 1, bias_correct, 
    marg_NAA_sigma, 
    small_dim, 0, 0);
  if(trace) see(SPRAA_all.dim);
  array<Type> SPR_all(n_stocks, n_regions,n_regions);
  SPR_all.setZero();
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) {
    for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) SPR_all(s,r,rr) += SPRAA_all(s,a,r,rr);
  }
  if(trace) see(SPR_all);
  array<Type> YPR_srf = get_YPR_srf(fleet_regions, fleet_seasons, can_move, mig_type, FAA_XSPR, log_M_avg, mu_avg, L_avg, waa_catch_avg, 
    fracyr_seasons, 0, bias_correct, 
    marg_NAA_sigma, 
    small_dim); //n_stocks x n_regions x n_fleets (should be 0 for regions not fleet_regions(f)-1)
  //for each stock/fleet and also the weighted average/total
  if(trace) see(YPR_srf);
  array<Type> log_SPR(1,n_stocks + 1), log_SPR0(1,n_stocks+1),log_Y_XSPR(1,n_fleets+n_regions+1), log_SSB_XSPR(1,n_stocks+1), log_YPR_XSPR(n_stocks,n_fleets+1); 
  log_SPR.setZero(); log_SPR0.setZero(); log_Y_XSPR.setZero(); log_SSB_XSPR.setZero(); log_YPR_XSPR.setZero();
  
  for(int s = 0; s < n_stocks; s++) {
    log_SPR(0,s) = log(SPR_all(s,spawn_regions(s)-1,spawn_regions(s)-1));
    log_SPR(0,n_stocks) += SPR_weights(s) * SPR_all(s,spawn_regions(s)-1,spawn_regions(s)-1);
    log_SPR0(0,s) = log(SPR0_all(s,spawn_regions(s)-1,spawn_regions(s)-1));
    log_SPR0(0,n_stocks) += SPR_weights(s) * SPR0_all(s,spawn_regions(s)-1,spawn_regions(s)-1);
    log_SSB_XSPR(0,s) = log(R_XSPR(s)) + log_SPR(s);
    log_SSB_XSPR(0,n_stocks) += R_XSPR(s) * SPR_all(s,spawn_regions(s)-1,spawn_regions(s)-1);
    for(int f = 0; f < n_fleets; f++) {
      log_YPR_XSPR(s,f) = log(YPR_srf(s,spawn_regions(s)-1,f)); //originate in spawning region and end up captured by fleet f
      log_YPR_XSPR(s,n_fleets) += YPR_srf(s,spawn_regions(s)-1,f); //last n_stocks elements are the totals across fleets. not logged yet
      log_Y_XSPR(0,f) += R_XSPR(s) * YPR_srf(s,spawn_regions(s)-1,f); //not logged yet
      log_Y_XSPR(0,n_fleets+fleet_regions(f)-1) += R_XSPR(s) * YPR_srf(s,spawn_regions(s)-1,f); //yield by region, not logged_yet
    }
    log_YPR_XSPR(s,n_fleets) = log(log_YPR_XSPR(s,n_fleets)); //log total YPR (across fleets) for stock s
  }
  //total yield
  for(int f = 0; f < n_fleets; f++) log_Y_XSPR(0,n_fleets+n_regions) += log_Y_XSPR(0,f);
  // now log the vectors
  log_Y_XSPR = log(log_Y_XSPR); 
  log_SPR(0,n_stocks) = log(log_SPR(0,n_stocks));
  log_SPR0(0,n_stocks) = log(log_SPR0(0,n_stocks));
  log_SSB_XSPR(0,n_stocks) = log(log_SSB_XSPR(0,n_stocks));
  //see(log_SPR);
  //see(log_SPR0);
  vector< array<Type> > res(17); 
  res(0) = log_FAA_XSPR; // log_FAA at FXSPR by fleet and across fleets
  res(1) = log_SSB_XSPR; //log_SSB_FXSPR
  res(2) = log_Y_XSPR; //log_Y_FXSPR
  res(3) = log_SPR; //stock specific log SPRs at FXSPR, (only the weighted sum will be X*SPR0/100)
  res(4) = log_SPR0;
  res(5) = log_YPR_XSPR; 
  array<Type> log_FXSPR_iter_a(1,n_iter);
  for(int i = 0; i < n_iter; i++) log_FXSPR_iter_a(0,i) = log_FXSPR_iter(i);
  if(trace) see(log_FXSPR_iter_a);
  res(6) = log_FXSPR_iter_a;
  if(trace) see(NAAPR0_all.dim);
  res(7) = NAAPR0_all;
  if(trace) see(NAAPR_FXSPR_all.dim);
  res(8) = NAAPR_FXSPR_all;
  if(trace) see(YPR_srf.dim);
  res(9) = YPR_srf;
  if(trace) see(waa_ssb_avg);
  // array<Type> temp1 = waa_ssb_avg;
  res(10) = waa_ssb_avg;
  if(trace) see(waa_catch_avg);
  res(11) = waa_catch_avg;
  if(trace) see(mat);
  res(12) = mat;
  if(trace) see(sel);
  res(13) = sel;
  if(trace) see(FAA_avg);
  res(14) = FAA_avg;
  if(trace) see(log_M_avg.dim);
  res(15) = log_M_avg;
  if(trace) see(mu_avg.dim);
  res(16) = mu_avg;
  if(trace) see("end get_SPR_res")
  return res;
}

template <class Type>
vector< array <Type> > get_annual_SPR_res(vector<Type> SPR_weights, array<Type> log_M, array<Type> FAA, vector<int> spawn_seasons,  
  vector<int> spawn_regions,
  vector<int> fleet_regions, 
  matrix<int> fleet_seasons,
  vector<Type> fracyr_seasons,
  array<int> can_move,
  array<int> must_move,
  vector<int> mig_type,
  array<Type> trans_mu_base, 
  matrix<Type> L,
  vector<int> which_F_age, array<Type> waa_ssb, array<Type> waa_catch,
  array<Type> mature, Type percentSPR, array<Type> NAA, matrix<Type> fracyr_SSB, vector<Type> F_init,  
  matrix<Type> R_XSPR,
  int small_dim, int SPR_weight_type, 
  int bias_correct,
  array<Type> marg_NAA_sigma,
  int trace = 0, int n_iter = 10){
  int ny = which_F_age.size();
  int n_fleets = waa_catch.dim(0);
  int n_regions = can_move.dim(2);
  int n_stocks = waa_ssb.dim(0);
  int n_ages = mature.dim(2);
  vector< array <Type>> all_res(7);
  array<Type> log_FAA_XSPR(n_fleets+n_regions+1,ny,n_ages); //log FAA_XSPR, FAA_XSPR_tot
  array<Type> log_SSB_XSPR(ny,n_stocks+1); //log SSB_XSPR, SSB_XSPR_tot
  array<Type> log_Y_XSPR(ny, n_fleets+n_regions+1); //log Y_XSPR, Y_XSPR_r, Y_XSPR_tot
  array<Type> log_SPR_XSPR(ny,n_stocks+1); //log SPR_XSPR
  array<Type> log_SPR0(ny,n_stocks+1); //log SPR0
  array<Type> log_YPR_XSPR(n_stocks,n_fleets+1,ny); //log YPR at FXSPR by stock and fleet and total across fleets by stock
  array<Type> log_FXSPR_iter(ny,n_iter); //log FXSPR_iter
  //get inputs for each years
  vector<int> yvec(1);

  for(int y = 0; y < ny; y++){
    yvec(0) = y;
    vector< array<Type>> SPR_res_y = get_SPR_res(SPR_weights, log_M, FAA, spawn_seasons,  spawn_regions, fleet_regions, 
      fleet_seasons, fracyr_seasons, can_move, must_move, mig_type, trans_mu_base, L, which_F_age(y), 
      waa_ssb, waa_catch, mature, percentSPR, NAA, fracyr_SSB, F_init(y), yvec, yvec, yvec, yvec, yvec, yvec, yvec, 
      vector<Type> (R_XSPR.row(y)), small_dim, SPR_weight_type, bias_correct, 
      marg_NAA_sigma, 
      trace = trace, n_iter = n_iter);
    for(int f = 0; f <= n_fleets+n_regions; f++) for(int a = 0; a < n_ages; a++){
      log_FAA_XSPR(f,y,a) = SPR_res_y(0)(f,a);
    }
    for(int s = 0; s <= n_stocks; s++) {
      log_SSB_XSPR(y,s) = SPR_res_y(1)(0,s);
      log_SPR_XSPR(y,s) = SPR_res_y(3)(0,s);
      log_SPR0(y,s) = SPR_res_y(4)(0,s);
    }
    for(int f = 0; f <= n_fleets+n_regions; f++) log_Y_XSPR(y,f) = SPR_res_y(2)(0,f);
    for(int s = 0; s < n_stocks; s++) {
      for(int f = 0; f <= n_fleets; f++) log_YPR_XSPR(s,f,y) = SPR_res_y(5)(s,f);
    }
    for(int i = 0; i < n_iter; i++) log_FXSPR_iter(y,i) = SPR_res_y(6)(0,i);
  }
  
  all_res(0) = log_FAA_XSPR;
  all_res(1) = log_SSB_XSPR;
  all_res(2) = log_Y_XSPR;
  all_res(3) = log_SPR_XSPR;
  all_res(4) = log_SPR0;
  all_res(5) = log_YPR_XSPR;
  all_res(6) = log_FXSPR_iter; 
  return all_res;
}

template <class Type>
array <Type> get_annual_SPR0_at_age(array<Type> log_M, vector<int> spawn_seasons,  
  vector<Type> fracyr_seasons,
  array<int> can_move,
  array<int> must_move,
  vector<int> mig_type,
  array<Type> trans_mu_base, 
  matrix<Type> L,
  array<Type> waa_ssb, 
  array<Type> mature, matrix<Type> fracyr_SSB,
  int bias_correct,
  array<Type> marg_NAA_sigma,
  int small_dim, int trace = 0){
  
  int ny = log_M.dim(2);
  int n_seasons = fracyr_seasons.size();
  int n_regions = can_move.dim(2);
  int n_stocks = waa_ssb.dim(0);
  int n_ages = mature.dim(2);
  
  array<Type> SPR0AA(ny, n_stocks, n_ages, n_regions, n_regions);
  //get inputs for each years
  vector<int> yvec(1);
  vector<int> fleet_regions(1);
  fleet_regions(0) = 1;
  matrix<int> fleet_seasons(1,n_seasons);
  fleet_seasons.setZero();
  array<Type> FAA0(1,n_ages);
  FAA0.setZero();
  // see(ny);
  // see(n_stocks);
  // see(n_regions);
  // see(n_ages);
  for(int y = 0; y < ny; y++){
    yvec(0) = y;
    //get average inputs over specified years
    // see(yvec(0));
    vector<Type> ssbfrac = get_avg_ssbfrac(fracyr_SSB,yvec);
    // see(ssbfrac);
    array<Type> waa_ssb_avg = get_avg_mat_as_array(waa_ssb, yvec);
    // see(waa_ssb_avg);
    vector<Type> L_avg = get_avg_L(L, yvec, 0);
    // see(L_avg);
    array<Type> mat = get_avg_mat_as_array(mature,yvec);
    // see(mat);
    array<Type> log_M_avg = get_avg_M(log_M, yvec, 1);
    // see(log_M_avg.dim);
    array<Type> mu_avg(n_stocks, n_ages, n_seasons, n_regions, n_regions);
    mu_avg.setZero(); 
    if(n_regions>1) mu_avg = get_avg_mu(trans_mu_base,yvec,mig_type, can_move, must_move);
    // see(mu_avg.dim);

    array<Type> SPR0AA_y = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, FAA0, log_M_avg, mu_avg, L_avg, 
      mat, waa_ssb_avg, fracyr_seasons, 1, bias_correct, 
      marg_NAA_sigma, 
      small_dim, 0);
    // see(SPR0AA_y.dim);
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++)for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++){
      // see(s);
      // see(a);
      // see(r);
      // see(rr);
      SPR0AA(y,s,a,r,rr) = SPR0AA_y(s,a,r,rr);
    }
  }
  
  return SPR0AA;
}
