template <class Type>
array<Type> get_eq_SAA(int y, vector<int> spawn_seasons, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, matrix<Type> fracyr_SSB, array<Type> FAA, array<Type> log_M_base, array<Type> mu, matrix<Type> L, 
  array<Type> mature, array<Type> waa, vector<int> waa_pointer_ssb, vector<Type> fracyr_seasons, int small_dim){
  /* 
    calculate equilibrium spawning biomass per recruit (at age) by stock and region. If movement is set up approriately 
    all fish can be made to return to a single spawning region for each stock.
                  y: the model year for which to use SPR inputs
       spawn_season: vector of indicators telling which season spawning occurs for each stock
      fleet_regions: vector of indicators telling which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_years x n_stocks:  size of interval from beginning of season to time of spawning within that season
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M_base: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions array of movement matrices
                  L: n_years_model x n_regions. "extra" unobserved mortality
             mature: n_stocks x n_years x n_ages proportion mature at age
                waa: (n_?) x n_years x n_ages_model. weight at age
    waa_pointer_ssb: n_stocks; which waa matrix to use for ssb
     fracyr_seasons: n_seasons: length of intervals for each season
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_stocks = log_M_base.dim(0);
  int n_regions = log_M_base.dim(1);
  int n_ages = log_M_base.dim(3);
  int n_fleets = FAA.dim(0);
  int n_seasons = can_move.dim(2);
  int P_dim = n_regions + n_fleets + 1;
  array<Type> SAA(n_stocks,n_ages,n_regions,n_regions); //SSB/R at age in each region column, given recruited in region row
  SAA.setZero();
  matrix<Type> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < n_regions; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<Type> S_ya = get_S(I, n_regions);
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) SAA(s,0,i,j) = S_ya(i,j);
    for(int a = 0; a < n_ages-1; a++) {
      matrix<Type> P_ya = I; //PTM for year and age and up to time of spawning
      for(int t = 0; t < n_seasons; t++) {
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
      } 
      S_ya = S_ya * get_S(P_ya, n_regions); //accumulate for next age
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) SAA(s,a,i,j) = S_ya(i,j);
    }
    //now plus group
    matrix<Type> P_ya = I; //PTM for year and age and up to time of spawning
    for(int t = 0; t < n_seasons; t++) {
      //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
      P_ya = P_ya * get_P(n_ages-1, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
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
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) SAA(s,n_ages-1,i,j) = S_ya(i,j);
  }
  return SAA;
}

template <class Type>
array<Type> get_YPR(vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, matrix<Type> FAA, array<Type> log_M_base, array<Type> mu, matrix<Type> L, array<Type> waacatch, 
  vector<Type> fracyr_seasons, int age_specific, int small_dim){
  /* 
    calculate equilibrium yield per recruit (at age) by stock and region.
        fleet_regions: vector of indicators telling which region each fleet is operating
        fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
             can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
             mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
                  FAA: fishing mortality: n_fleets x n_ages
           log_M_base: log M (density-independent components): n_stocks x n_regions x n_ages
                   mu: n_stocks x n_ages x n_seasons x n_regions x n_regions array of movement matrices
                    L: n_regions. "extra" unobserved mortality
             waacatch: n_fleets x n_ages. weight at age
       fracyr_seasons: n_seasons: length of intervals for each season
         age_specific: 0/1 telling whether to return SSB/R by age or not. dimensions of returned arrays are different. 
            small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_stocks = log_M_base.dim(0);
  int n_regions = log_M_base.dim(1);
  int n_ages = log_M_base.dim(3);
  int n_fleets = FAA.dim(0);
  int n_seasons = can_move.dim(2);
  int P_dim = n_regions + n_fleets + 1;
  array<Type> YPRAA(n_stocks,n_ages,n_regions,n_fleets); //Yield/R at age in each fleet column, given recruited in region row
  YPRAA.setZero();
  matrix<Type> W(n_fleets,n_fleets);
  W.setZero();

  matrix<Type> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < n_regions; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<Type> S_ya = get_S(I, n_regions);
    for(int a = 0; a < n_ages-1; a++) {
      for(int i = 0; i < n_fleets; i++) W(i,i) = waacatch(i,a);
      matrix<Type> P_ya = I; //PTM for year and age
      for(int t = 0; t < n_seasons; t++) {
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
      }
      // YPR at year and age = prob alive to up to age a-1 x prob caught at age a x waa
      matrix<Type> YPR_ya = S_ya * get_D(P_ya) * W; //should be n_regions x n_fleets
      S_ya = S_ya * get_S(P_ya, n_regions);  //accumulate for next age
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) YPRAA(s,a,i,j) = YPR_ya(i,j);
    }
    //now plus group
    for(int i = 0; i < n_fleets; i++) W(i,i) = waacatch(i,n_ages-1);
    matrix<Type> P_ya = I; //PTM for year and age
    for(int t = 0; t < n_seasons; t++) {
      //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
      P_ya = P_ya * get_P(n_ages-1, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
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
    matrix<Type> YPR_ya = S_ya * get_D(P_ya) * W; //should be n_regions x n_fleets
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
array<Type> get_SPR(int y, vector<int> spawn_seasons, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, matrix<Type> fracyr_SSB, array<Type> FAA, array<Type> log_M_base, array<Type> mu, matrix<Type> L, 
  array<Type> mature, array<Type> waa, vector<int> waa_pointer_ssb, vector<Type> fracyr_seasons, int age_specific, int small_dim){
  /* 
    calculate equilibrium spawning biomass per recruit (at age) by stock and region. If movement is set up approriately 
    all fish can be made to return to a single spawning region for each stock.
                  y: the model year for which to use SPR inputs
       spawn_season: vector of indicators telling which season spawning occurs for each stock
      fleet_regions: vector of indicators telling which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_years x n_stocks:  size of interval from beginning of season to time of spawning within that season
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M_base: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions array of movement matrices
                  L: n_years_model x n_regions. "extra" unobserved mortality
             mature: n_stocks x n_years x n_ages proportion mature at age
                waa: (n_?) x n_years x n_ages_model. weight at age
    waa_pointer_ssb: which waa matrix to use for ssb
     fracyr_seasons: n_seasons: length of intervals for each season
       age_specific: 0/1 telling whether to return SSB/R by age or not. dimensions of returned arrays are different. 
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_stocks = log_M_base.dim(0);
  int n_regions = log_M_base.dim(1);
  int n_ages = log_M_base.dim(3);
  int n_fleets = FAA.dim(0);
  int n_seasons = can_move.dim(2);
  int P_dim = n_regions + n_fleets + 1;
  array<Type> SPRAA(n_stocks,n_ages,n_regions,n_regions); //SSB/R at age in each region column, given recruited in region row
  SPRAA.setZero();

  matrix<Type> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < n_regions; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<Type> S_ya = get_S(I, n_regions);
    for(int a = 0; a < n_ages-1; a++) {
      matrix<Type> P_ya = I, P_spawn = I; //PTM for year and age and up to time of spawning
      for(int t = 0; t < n_seasons; t++) {
        if(t == spawn_seasons(s)-1) {
          //P(0,t_spawn): PTM over entire year up to time of spawning
          P_spawn = P_ya * get_P(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(y,s), FAA, log_M_base, mu, L);
        }
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
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
        P_spawn = P_ya * get_P(n_ages-1, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(y,s), FAA, log_M_base, mu, L);
      }
      //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
      P_ya = P_ya * get_P(n_ages-1, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
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

template <class Type>
array<Type> get_SPR(vector<int> spawn_seasons, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, vector<Type> fracyr_SSB, matrix<Type> FAA, array<Type> log_M_base, array<Type> mu, vector<Type> L, 
  matrix<Type> mature, matrix<Type> waa_ssb, vector<Type> fracyr_seasons, int age_specific, int small_dim){
  /* 
    calculate equilibrium spawning biomass per recruit (at age) by stock and region. If movement is set up approriately 
    all fish can be made to return to a single spawning region for each stock.
    NB: This version of get_SPR does not allow year indices. It is assumed that the yearly (or averaged) values are provided.
       spawn_season: vector of indicators telling which season spawning occurs for each stock
      fleet_regions: vector of indicators telling which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_stocks:  size of interval from beginning of season to time of spawning within that season
                FAA: fishing mortality: n_fleets x n_ages
         log_M_base: log M (density-independent components): n_stocks x n_regions x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_regions x n_regions array of movement matrices
                  L: n_regions. "extra" unobserved mortality
             mature: n_stocks x n_ages proportion mature at age
            waa_ssb: n_stocks x n_ages. weight at age
     fracyr_seasons: n_seasons: length of intervals for each season
       age_specific: 0/1 telling whether to return SSB/R by age or not. dimensions of returned arrays are different. 
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_stocks = log_M_base.dim(0);
  int n_regions = log_M_base.dim(1);
  int n_ages = log_M_base.dim(3);
  int n_fleets = FAA.dim(0);
  int n_seasons = can_move.dim(2);
  int P_dim = n_regions + n_fleets + 1;
  array<Type> SPRAA(n_stocks,n_ages,n_regions,n_regions); //SSB/R at age in each region column, given recruited in region row
  SPRAA.setZero();

  matrix<Type> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < n_regions; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<Type> S_ya = get_S(I, n_regions);
    for(int a = 0; a < n_ages-1; a++) {
      matrix<Type> P_ya = I, P_spawn = I; //PTM for age and up to time of spawning
      for(int t = 0; t < n_seasons; t++) {
        if(t == spawn_seasons(s)-1) {
          //P(0,t_spawn): PTM over entire year up to time of spawning
          P_spawn = P_ya * get_P(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(s), FAA, log_M_base, mu, L);
        }
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
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
        P_spawn = P_ya * get_P(n_ages-1, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(s), FAA, log_M_base, mu, L);
      }
      //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
      P_ya = P_ya * get_P(n_ages-1, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
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
array<Type> get_SPR_0(int y, vector<int> spawn_seasons, array<int> can_move, 
  vector<int> mig_type, matrix<Type> fracyr_SSB, array<Type> log_M_base, array<Type> mu, matrix<Type> L, 
  array<Type> mature, array<Type> waa, vector<int> waa_pointer_ssb, vector<Type> fracyr_seasons, int age_specific, int small_dim){
  /* 
    calculate unfished equilibrium spawning biomass per recruit (at age) by stock and region. If movement is set up approriately 
    all fish can be made to return to a single spawning region for each stock.
                  y: the model year for which to use SPR inputs
       spawn_season: vector of indicators telling which season spawning occurs for each stock
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_years x n_stocks:  size of interval from beginning of season to time of spawning within that season
         log_M_base: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions array of movement matrices
                  L: n_years_model x n_regions. "extra" unobserved mortality
             mature: n_stocks x n_years x n_ages proportion mature at age
                waa: (n_?) x n_years x n_ages_model. weight at age
    waa_pointer_ssb: which waa matrix to use for ssb
     fracyr_seasons: n_seasons: length of intervals for each season
       age_specific: 0/1 telling whether to return SSB/R by age or not. dimensions of returned arrays are different. 
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_seasons = can_move.dim(2);
  int n_ages = can_move.dim(1);
  vector<int> fleet_regions(1);
  fleet_regions(0) = 1;
  matrix<int> fleet_seasons(1,n_seasons);
  fleet_seasons.setZero();
  array<Type> FAA(1,y+1,n_seasons,n_ages);
  FAA.setZero();

  array<Type> SPR_0 = get_SPR(y,spawn_seasons,fleet_regions,fleet_seasons, can_move,mig_type,fracyr_SSB,FAA,log_M_base,mu,L,mature,
    waa,waa_pointer_ssb, fracyr_seasons, age_specific, small_dim);
  
  return SPR_0;
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
  matrix<Type> fracyr_SSB;
  matrix<Type> selectivity;
  array<Type> log_M_base;
  array<Type> mu; 
  matrix<Type> L;
  matrix<Type> mature;
  matrix<Type> waa_ssb; 
  vector<Type> fracyr_seasons;
  vector<Type> SPR_weights; //how to weight stock-specific SSB/R for aggregate SSB/R.
  int age_specific; 
  int small_dim;

  // Constructor 
  spr_F_spatial(
  vector<int> spawn_seasons_,
  vector<int> spawn_regions_,
  vector<int> fleet_regions_, 
  matrix<int> fleet_seasons_, 
  array<int> can_move_, 
  vector<int> mig_type_,
  matrix<Type> fracyr_SSB_,
  matrix<Type> selectivity_,
  array<Type> log_M_base_,
  array<Type> mu_, 
  matrix<Type> L_,
  matrix<Type> mature_,
  matrix<Type> waa_ssb_, 
  vector<Type> fracyr_seasons_,
  vector<Type> SPR_weights_,
  int age_specific_, 
  int small_dim_) :
    spawn_seasons(spawn_seasons_), 
    spawn_regions(spawn_regions_), 
    fleet_regions(fleet_regions_),
    fleet_seasons(fleet_seasons_),
    can_move(can_move_),
    mig_type(mig_type_),
    fracyr_SSB(fracyr_SSB_),
    selectivity(selectivity_),
    log_M_base(log_M_base_),
    mu(mu_),
    L(L_),
    mature(mature_),
    waa_ssb(waa_ssb_),
    fracyr_seasons(fracyr_seasons_),
    SPR_weights(SPR_weights_),
    age_specific(age_specific_),
    small_dim(small_dim_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_stocks = log_M_base.dim(0);
    int n_regions = log_M_base.dim(1);
    int n_ages = log_M_base.dim(3);
    int n_fleets = fleet_regions.size();
    int n_seasons = can_move.dim(2);
    matrix<Type> FAA(n_fleets,n_ages);
    FAA.setZero();
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) {
      FAA(f,a) = selectivity(f,a) * exp(Type(log_F));
    }
    //stock-specific SSB/R
    array<Type> SPR_sr = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB, FAA, log_M_base, mu, L, mature, waa_ssb, 
      fracyr_seasons, 0, small_dim);

    //weighted-average of SSB/R returned
    Type SPR = 0;
    for(int s = 0; s < n_stocks; s++) SPR += SPR_weights(s) * SPR_sr(s,spawn_regions(s)-1,spawn_regions(s)-1); 
    
    return SPR.template cast<T>();
  }
};

template <class Type>
vector< vector <Type> > get_SPR_res(vector<Type> SPR_weights, array<Type> log_M_base, array<Type> FAA, vector<int> spawn_seasons,  
  vector<int> spawn_regions,
  vector<int> fleet_regions, 
  matrix<int> fleet_seasons,
  vector<Type> fracyr_seasons,
  array<int> can_move,
  array<int> must_move,
  vector<int> mig_type,
  array<Type> trans_mu_base, 
  matrix<Type> L,
  matrix<int> which_F_age, array<Type> waa, vector<int> waa_pointer_ssb, 
  vector<int> waa_pointer_fleets,
  array<Type> mature, Type percentSPR, array<Type> NAA, matrix<Type> fracyr_SSB, Type F_init, 
  vector<int> years_M, vector<int> years_mu, vector<int> years_L, vector<int> years_mat, vector<int> years_sel, 
  vector<int> years_waa_ssb, vector<int> years_waa_catch, vector<int> years_R,
  int small_dim, int SPR_weight_type) {
  
  int n = 10;
  int n_stocks = log_M_base.dim(0);
  int n_fleets = FAA.dim(0);
  int n_regions = NAA.dim(1);
  int n_seasons = fleet_seasons.cols();
  int n_ages = log_M_base.dim(4);
  matrix<Type> waa_ssb(n_stocks,n_ages);
  matrix<Type> waa_catch(n_fleets,n_ages);
  matrix<Type> sel(n_fleets, n_ages);
  array<Type> M(n_stocks,n_regions,n_ages);
  array<Type> log_M_base_avg(n_stocks,n_regions,n_ages);
  array<Type> mu_avg(n_stocks,n_seasons,n_ages, n_regions, n_regions);
  vector<Type> L_avg(n_regions);
  matrix<Type> mat(n_stocks, n_ages);
  vector<Type> ssbfrac(n_stocks), R(n_stocks);

  waa_ssb.setZero(); waa_catch.setZero(); sel.setZero(); M.setZero(); mu_avg.setZero(); L_avg.setZero(); 
  mat.setZero(); R.setZero(); ssbfrac.setZero(); log_M_base_avg.setZero();
  //get average inputs over specified years
  for(int y = 0; y < years_L.size(); y++) for(int i = 0; i < n_regions; i++) L_avg(i) += L(years_L(y),i)/years_L.size();
  for(int s = 0; s < n_stocks; s++) {
    for(int y = 0; y < years_R.size(); y++) R(s) += NAA(s,spawn_regions(s)-1,years_R(y),0)/years_R.size();
    for(int y = 0; y < years_mat.size(); y++) for(int a = 0; a < n_ages; a++) mat(s,a) += mature(s,years_mat(y),a)/years_mat.size();
    for(int r = 0; r < n_regions; r++) for(int a = 0; a < n_ages; a++) {
      for(int y = 0; y < years_M.size(); y++) M(s,r,a) += exp(log_M_base(s,r,years_M(y),a))/years_M.size();
      log_M_base_avg(s,r,a) = log(M(s,r,a)); //needed for SPR function
    }
    for(int y = 0; y < years_waa_ssb.size(); y++) {
      ssbfrac(s) += fracyr_SSB(s,years_waa_ssb(y))/years_waa_ssb.size();
      for(int a = 0; a < n_ages; a++) waa_ssb(s,a) += waa(waa_pointer_ssb(s)-1, years_waa_ssb(y),a)/years_waa_ssb.size();
    }
    for(int t = 0; t < n_seasons; t++) for(int a = 0; a < n_ages; a++){
      matrix<Type> tmp = get_avg_mu_matrix(s, a, t, years_mu, mig_type, can_move, must_move, trans_mu_base);
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
        mu_avg(s,t,a,i,j) += tmp(i,j);
      }
    }
  }

  for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++){
    for(int y = 0; y < years_waa_catch.size(); y++)  waa_catch(f,a) += waa(waa_pointer_fleets(f)-1, years_waa_catch(y),a)/years_waa_catch.size();
    for(int y = 0; y < years_sel.size(); y++) for(int a = 0; a < n_ages; a++){
      sel(f,a) += FAA(f,years_sel(y),a)/years_sel.size(); //average F at fleet,season,age over years
    }
  }
  //define SPR weighting 
  if(SPR_weight_type == 0){ //use average Recruitment
    SPR_weights = R/R.sum();
  } else { //use user-specified weights as provided. do nothing
  }

  //which_F_age needs to be set appropriately by user.
  //if an arbitrary index is used, then full F would need to be defined on R side. Could pick out from ADREPORTED FAA
  sel = sel/sel(which_F_age(1)-1,which_F_age(0)-1);
  matrix<Type> FAA0(n_fleets,n_ages);
  FAA0.setZero();
  array<Type> SPR0_all = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, FAA0, log_M_base_avg, mu_avg, L_avg, 
    mat, waa_ssb, fracyr_seasons, 0, small_dim);
  Type SPR0 = 0;
  for(int s = 0; s < n_stocks; s++) SPR0 += SPR_weights(s) * SPR0_all(s,spawn_regions(s)-1,spawn_regions(s)-1); 

  vector<Type> log_FXSPR_i(1), log_FXSPR_iter(n);
  log_FXSPR_iter(0) = log(F_init);

  spr_F_spatial<Type> sprF(spawn_seasons, spawn_regions, fleet_regions, can_move, mig_type, ssbfrac, sel, log_M_base_avg, mu_avg, L_avg, mat, 
    waa_ssb, fracyr_seasons, SPR_weights, 0, small_dim);
  for(int i=0; i<n-1; i++) {
    log_FXSPR_i(0) = log_FXSPR_iter(i);
    vector<Type> grad_spr_F = autodiff::gradient(sprF,log_FXSPR_i);
    log_FXSPR_iter(i+1) = log_FXSPR_iter(i) - (sprF(log_FXSPR_i) - 0.01*percentSPR * SPR0)/grad_spr_F(0);
  }
  //fix this: FAA and sel DO NOT NEED TO BE SEASONAL!
  matrix<Type> FAA_XSPR = exp(log_FXSPR_iter(n-1)) * sel;
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

  array<Type> SPR_all = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, ssbfrac, FAA_XSPR, log_M_base_avg, mu_avg, L_avg, 
    mat, waa_ssb, fracyr_seasons, 0, small_dim);
  array<Type> YPR_all = get_YPR(fleet_regions, fleet_seasons, can_move, mig_type, FAA_XSPR, log_M_base_avg, mu_avg, L_avg, waa_catch, 
    fracyr_seasons, 0, small_dim); //n_stocks x n_regions x n_fleets (should be 0 for regions not fleet_regions(f)-1)
  //for each stock/fleet and also the weighted average/total
  vector<Type> log_SPR(n_stocks + 1), log_SPR0(n_stocks+1), log_Y_XSPR(n_fleets+1), log_SSB_XSPR(n_stocks+1); 
  log_SPR.setZero(); log_SPR0.setZero(); log_Y_XSPR.setZero(); log_SSB_XSPR.setZero();
  
  for(int s = 0; s < n_stocks; s++) {
    log_SPR(s) = log(SPR_all(s,spawn_regions(s)-1,spawn_regions(s)-1));
    log_SPR(n_stocks) += SPR_weights(s) * SPR_all(s,spawn_regions(s)-1,spawn_regions(s)-1);
    log_SPR0(s) = log(SPR0_all(s,spawn_regions(s)-1,spawn_regions(s)-1));
    log_SSB_XSPR(s) = log(R(s)) + log_SPR(s);
    log_SSB_XSPR(n_stocks) += R(s) * SPR_all(s,spawn_regions(s)-1,spawn_regions(s)-1);
    for(int f = 0; f < n_fleets; f++) {
      log_Y_XSPR(f) += log(R(s)) + log(YPR_all(s,fleet_regions(f)-1,f));
      log_Y_XSPR(n_fleets) += exp(log_Y_XSPR(f));
    }
  }
  log_SPR(n_stocks) = log(log_SPR(n_stocks));
  log_SSB_XSPR(n_stocks) = log(log_SSB_XSPR(n_stocks));
  log_Y_XSPR(n_fleets) = log(log_Y_XSPR(n_fleets));
  log_SPR0(n_stocks) = log(SPR0);
  vector< vector<Type> > res(6); 
  res(0) = log_FAA_XSPR; // log_FAA at FXSPR by fleet and across fleets
  res(1) = log_SSB_XSPR; //log_SSB_FXSPR
  res(2) = log_Y_XSPR; //log_Y_FXSPR
  res(3) = log_SPR; //stock specific log SPRs at FXSPR, (only the weighted sum will be X*SPR0/100)
  res(4) = log_SPR0;
  res(5) = log_FXSPR_iter; //last value is max F at X%SPR across ages and fleets
  return res;
}


